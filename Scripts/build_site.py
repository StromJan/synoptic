#!/usr/bin/env python
"""Render the synoptic overview maps and build a static website.

    python build_site.py                      # all products, all CPU cores
    python build_site.py --products ReTop     # just one product (quick check)
    python build_site.py --workers 2

Reads  ../Data/pl_*.grib, ../Data/sfc.grib, ../Data/DWD_analysis/*.png,
       ../Data/meteograms/helsinki_meteogram.png
Writes ../site/index.html and ../site/img/...
"""
import argparse
import json
import os
import re
import shutil
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                     # no display on servers / GitHub runners

import glob
import numpy as np
import pandas as pd
import xarray as xr
import metpy
import metpy.calc
from metpy.units import units

import earthkit.data as ekd
import earthkit.plots as ekp
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmcrameri.cm as cm


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "../Data"
SITE_DIR = HERE / "../site"
TEMPLATE = HERE / "site_template.html"

DOMAIN = "North Atlantic"
LEV = "level"

# Cut the global fields to a box around DOMAIN (lon_min, lon_max, lat_min, lat_max)
# before plotting. earthkit's "North Atlantic" is 50W-35E, 32-72N; the margin
# covers the curved edges of the projection. ~7x less data = much faster maps.
# Set to None to plot the global fields.
CROP_BOX = (-110, 95, 10, 90)

for cmap in (cm.tokyo, cm.batlowW):
    if cmap.name not in mpl.colormaps:
        mpl.colormaps.register(cmap)


# ---------------------------------------------------------------------------
# Data loading (same as the notebook)
# ---------------------------------------------------------------------------
def open_grib(path, pressure_levels=False):
    if pressure_levels:
        # squeeze=False keeps "level" even in single-level files like pl_850.grib
        ds = ekd.from_source("file", str(path)).to_xarray(add_valid_time_coord=True, squeeze=False)
        ds = ds.isel(forecast_reference_time=-1, level_type=0).drop_vars("level_type")
    else:
        ds = ekd.from_source("file", str(path)).to_xarray(add_valid_time_coord=True)
        if "forecast_reference_time" in ds.dims:     # several runs in the file: keep the latest
            ds = ds.isel(forecast_reference_time=-1)
    return ds


def fix_units(ds):
    """earthkit-data >= 1.0 writes units like "m.s**-1", which pint/MetPy can't parse."""
    for v in ds.data_vars:
        u = ds[v].attrs.get("units")
        if isinstance(u, str):
            ds[v].attrs["units"] = re.sub(r"(?<=[A-Za-z0-9])\.(?=[A-Za-z])", " ", u)
    return ds


def load_data():
    pl_files = sorted(glob.glob(str(DATA_DIR / "pl_*.grib")))
    if not pl_files:
        raise FileNotFoundError(f"No pl_*.grib files in {DATA_DIR.resolve()}")
    ds = xr.merge([open_grib(f, pressure_levels=True) for f in pl_files]
                  + [open_grib(DATA_DIR / "sfc.grib")],
                  compat="override", join="outer")
    return crop(fix_units(ds))


def crop(ds):
    if CROP_BOX is None:
        return ds
    lon_min, lon_max, lat_min, lat_max = CROP_BOX
    if float(ds.longitude.max()) > 180:          # 0..360 grid -> -180..180
        ds = ds.assign_coords(longitude=((ds.longitude + 180) % 360) - 180).sortby("longitude")
    lat = slice(lat_max, lat_min) if ds.latitude[0] > ds.latitude[-1] else slice(lat_min, lat_max)
    return ds.sel(longitude=slice(lon_min, lon_max), latitude=lat)


# ---------------------------------------------------------------------------
# Map definitions (copied from Synoptic_Overview notebook)
# ---------------------------------------------------------------------------
def render_SurfaceTemperaturePressure(frame, stamp, fc_step): 

    # do any precalculations necessary 

    temperature = frame["2t"]
    msl = frame["msl"]

    # set up the styles
    t2m_style = ekp.styles.Style(
            units="celsius",
            colors="RdBu_r",
            levels=np.arange(-40,41,2),
            )
    
    msl_style = ekp.styles.Contour(
            units="hPa",
            levels=np.arange(920,1081,4),
            colors="black",
            linewidths=[1.4, 0.8, 0.8, 0.8, 0.8],
            linestyles=["solid",],
            labels=True,
            legend_style=None
            )

    # plot all the variables of interest
    chart = ekp.Map(domain=DOMAIN)
    chart.title(f"{{variable_name}}\n{stamp} (+{fc_step} h)")
    
    

    # Surface temperature
    chart.contourf(temperature, style=t2m_style)

    # Mean sea level pressure
    chart.contour(msl, style=msl_style)

    chart.legend()

    return chart 


def render_ReTop(frame, stamp, fc_step): 

    # do any precalculations necessary 

    gphDiff = (frame["gh"].sel({LEV: 500}) - frame["gh"].sel({LEV: 1000})) 
    gph500 = frame["gh"].sel({LEV: 500})
    msl = frame["msl"]

    # set up the styles
    gphDiff_style = ekp.styles.Style(
            colors="batlowW",
            extend="both", 
            levels=np.arange(4400, 6010, 40),
            )
    
    gph500_style = ekp.styles.Contour(
            levels=np.arange(5000,6000,80),
            colors="black",
            linewidths=[1.4],
            linestyles=["solid",],
            labels=True,
            legend_style=None
            )
    
    msl_style = ekp.styles.Contour(
            units="hPa",
            levels=np.arange(920,1081,5),
            colors="white",
            linewidths=[1.2,],
            linestyles=["solid",],
            labels=True,
            legend_style=None
            )

    # plot all the variables of interest
    chart = ekp.Map(domain=DOMAIN)
    chart.title(f"{{variable_name}}\n{stamp} (+{fc_step} h)")
    

    # Surface temperature
    chart.contourf(gphDiff, style=gphDiff_style)

    # Geopotential height at 500 hPa
    chart.contour(gph500, style=gph500_style)

    # Mean sea level pressure
    chart.contour(msl, style=msl_style)

    chart.legend()

    return chart 


def render_Frontogenesis(frame, stamp, fc_step): 

    msl = frame["msl"].load()

    # do any precalculations necessary 
    frame = frame.sel({LEV: 850}).load()
    pt = metpy.calc.potential_temperature(850 * units.hPa, frame["t"])
    fgp = metpy.calc.frontogenesis(pt, frame["u"], frame["v"]) * 3.6E9
    

    # set up the styles
    fgp_style = ekp.styles.Style(
            units="K/(1000kmh)",
            colors="RdBu_r",
            extend="both",
            levels=np.arange(-15, 15.1, 0.5),
            )
    
    msl_style = ekp.styles.Contour(
            units="hPa",
            levels=np.arange(920,1081,5),
            colors="black",
            linewidths=[1.,],
            linestyles=["solid",],
            labels=True,
            legend_style=None
            )

    # plot all the variables of interest
    chart = ekp.Map(domain=DOMAIN)
    chart.title(f"{{variable_name}}\n{stamp} (+{fc_step} h)")
    

    # Frontogenesis parameter
    chart.contourf(fgp, style=fgp_style)

    # Mean sea level pressure
    chart.contour(msl, style=msl_style)

    chart.legend()

    return chart 


def render_EqPotTemperature(frame, stamp, fc_step): 

    msl = frame["msl"].load()

    # do any precalculations necessary 
    frame = frame.sel({LEV: 850}).load()
    td = metpy.calc.dewpoint_from_specific_humidity(850 * units.hPa, frame["t"],  frame["q"])
    ept = metpy.calc.equivalent_potential_temperature(850 * units.hPa, frame["t"], td).metpy.convert_units("degC")

    # set up the styles
    ept_style = ekp.styles.Style(
            units="degC",
            colors="tokyo",
            extend="both",
            levels=np.arange(-10, 81, 2),
            )
    
    msl_style = ekp.styles.Contour(
            units="hPa",
            levels=np.arange(920,1081,5),
            colors="black",
            linewidths=[1.,],
            linestyles=["solid",],
            labels=True,
            legend_style=None
            )

    # plot all the variables of interest
    chart = ekp.Map(domain=DOMAIN)
    chart.title(f"{{variable_name}}\n{stamp} (+{fc_step} h)")
    

    # Frontogenesis parameter
    chart.contourf(ept, style=ept_style)

    # Mean sea level pressure
    chart.contour(msl, style=msl_style)

    chart.legend()

    return chart 


def render_Tcwv(frame, stamp, fc_step): 

    msl = frame["msl"].load()
    tcwv = frame["tcwv"].load()

    # set up the styles
    tcwv_style = ekp.styles.Style(
            units_label="kg/m2",
            colors="GnBu",
            extend="both",
            levels=np.arange(0, 46, 1),
            )
    
    msl_style = ekp.styles.Contour(
            units="hPa",
            levels=np.arange(920,1081,5),
            colors="black",
            linewidths=[1.,],
            linestyles=["solid",],
            labels=True,
            legend_style=None
            )

    # plot all the variables of interest
    chart = ekp.Map(domain=DOMAIN)
    chart.title(f"{{variable_name}}\n{stamp} (+{fc_step} h)")
    

    # Frontogenesis parameter
    chart.contourf(tcwv, style=tcwv_style)

    # Mean sea level pressure
    chart.contour(msl, style=msl_style)

    chart.legend()

    return chart 


def render_JetDiagnose(frame, stamp, fc_step): 

    msl = frame["msl"]

    # do any precalculations necessary 

    frame = frame.sel({LEV: 200})
    speed = metpy.calc.wind_speed(frame["u"], frame["v"])
    

    # set up the styles
    speed_style = ekp.styles.Style(
            units="m/s",
            colors="YlGn",
            extend="max", 
            levels=np.arange(30,91,5),
            )
    
    msl_style = ekp.styles.Contour(
            units="hPa",
            levels=np.arange(920,1081,4),
            colors="black",
            linewidths=[1.4, 0.8, 0.8, 0.8, 0.8],
            linestyles=["solid",],
            labels=True,
            legend_style=None
            )

    # plot all the variables of interest
    chart = ekp.Map(domain=DOMAIN)
    chart.title(f"{{variable_name}}\n{stamp} (+{fc_step} h)")
    
    

    # Surface temperature
    chart.contourf(speed, style=speed_style)

    # Mean sea level pressure
    chart.contour(msl, style=msl_style)

    chart.legend()

    return chart


# ---------------------------------------------------------------------------
# What goes on the page, in order.  desc may contain simple HTML.
# ---------------------------------------------------------------------------
PRODUCTS = [
    dict(key="SurfaceTemperaturePressure", render=render_SurfaceTemperaturePressure,
         title="2 m temperature and MSLP",
         desc="<b>Shading:</b> 2 m temperature &nbsp; <b>Black contour:</b> MSLP"),
    dict(key="ReTop", render=render_ReTop,
         title="Large-scale synoptic situation",
         desc="<b>Shading:</b> 500&ndash;1000 hPa thickness &rarr; bulk temperature &nbsp; "
              "<b>Black contour:</b> 500 hPa geopotential height &nbsp; <b>White contour:</b> MSLP"),
    dict(key="EqPotTemperature", render=render_EqPotTemperature,
         title="Pseudo-equivalent potential temperature: air masses",
         desc="<b>Shading:</b> &theta;<sub>e</sub> at 850 hPa &nbsp; <b>Black contour:</b> MSLP"),
    dict(key="Frontogenesis", render=render_Frontogenesis,
         title="Frontogenesis: where fronts are forming",
         desc="<b>Shading:</b> frontogenesis parameter at 850 hPa &nbsp; <b>Black contour:</b> MSLP"),
    dict(key="JetDiagnose", render=render_JetDiagnose,
         title="Upper-level jet",
         desc="<b>Shading:</b> wind speed at 200 hPa &nbsp; <b>Black contour:</b> MSLP"),
    dict(key="Tcwv", render=render_Tcwv,
         title="Total column water vapour: how wet are we talking",
         desc="<b>Shading:</b> total column water vapour &nbsp; <b>Black contour:</b> MSLP"),
]
RENDERERS = {p["key"]: p["render"] for p in PRODUCTS}


def frame_times(frame):
    stamp = pd.Timestamp(frame["valid_time"].values).strftime("%d %b %Y, %H:%M UTC")
    fc_step = str(int(frame["step"].values / np.timedelta64(1, "h")))
    return stamp, fc_step


# ---------------------------------------------------------------------------
# Parallel rendering: each worker process opens the data once, then renders
# whichever (product, step) jobs it is given.
# ---------------------------------------------------------------------------
_DATA = None


def _init_worker():
    global _DATA
    _DATA = load_data()


def _render_job(job):
    key, i, path = job
    frame = _DATA.isel(step=i)
    stamp, fc_step = frame_times(frame)

    chart = RENDERERS[key](frame, stamp, fc_step)
    chart.coastlines()
    chart.borders()
    chart.gridlines()
    chart.save(str(path))
    plt.close("all")
    return key


def render_all(keys, n_steps, workers):
    jobs = []
    for key in keys:
        (SITE_DIR / "img" / key).mkdir(parents=True, exist_ok=True)
        jobs += [(key, i, SITE_DIR / "img" / key / f"{i:03d}.png") for i in range(n_steps)]

    t0, done = time.time(), 0
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as pool:
        for _ in pool.map(_render_job, jobs, chunksize=2):
            done += 1
            if done % 10 == 0 or done == len(jobs):
                print(f"  rendered {done}/{len(jobs)} maps ({time.time() - t0:.0f} s)", flush=True)


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------
def copy_image(src, dst_rel):
    src = Path(src)
    if not src.exists():
        print("  missing, skipped:", src)
        return None
    dst = SITE_DIR / dst_rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    return dst_rel


def extra_sections():
    """DWD charts and the ECMWF meteogram, if they were downloaded."""
    sections = []
    dwd = DATA_DIR / "DWD_analysis"

    src = copy_image(dwd / "Analysis.png", "img/dwd/analysis.png")
    if src:
        sections.append(dict(id="dwd-analysis", kind="image", title="Analysis by DWD", src=src))

    frames, labels = [], []
    for f in sorted(dwd.glob("Analysis_p*.png")):
        h = int(re.search(r"_p(\d+)", f.name).group(1))
        frames.append(copy_image(f, f"img/dwd/p{h:03d}.png"))
        labels.append(f"+{h} h")
    if frames:
        sections.append(dict(id="dwd-forecast", kind="viewer", title="Surface forecast charts by DWD",
                             frames=frames, labels=labels))

    src = copy_image(DATA_DIR / "meteograms/helsinki_meteogram.png", "img/meteogram_helsinki.png")
    if src:
        sections.append(dict(id="meteogram", kind="image", title="Ensemble forecast for Helsinki", src=src))
    return sections


def write_page(keys, data):
    run = pd.Timestamp(data["forecast_reference_time"].values)
    version = run.strftime("%Y%m%d%H")          # cache-buster: new URLs for every run

    labels = []
    for i in range(data.sizes["step"]):
        stamp, fc_step = frame_times(data.isel(step=i))
        labels.append(f"{stamp}  (+{fc_step} h)")

    sections = []
    for p in PRODUCTS:
        if p["key"] in keys:
            sections.append(dict(
                id=p["key"], kind="viewer", group="ifs", title=p["title"], desc=p["desc"],
                frames=[f"img/{p['key']}/{i:03d}.png" for i in range(len(labels))],
                labels=labels))
    sections += extra_sections()

    for s in sections:                          # append ?v=... to every image URL
        if "frames" in s:
            s["frames"] = [f"{f}?v={version}" for f in s["frames"]]
        if "src" in s:
            s["src"] = f"{s['src']}?v={version}"

    manifest = dict(
        run=run.strftime("%d %b %Y %H UTC"),
        built=datetime.now(timezone.utc).strftime("%d %b %Y %H:%M UTC"),
        domain=DOMAIN,
        sections=sections,
    )
    blob = json.dumps(manifest).replace("</", "<\\/")      # safe inside <script>
    html = TEMPLATE.read_text().replace("__MANIFEST__", blob)
    (SITE_DIR / "index.html").write_text(html)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--products", nargs="+", choices=list(RENDERERS), default=list(RENDERERS))
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()

    data = load_data()
    n_steps = data.sizes["step"]
    run = pd.Timestamp(data["forecast_reference_time"].values)
    print(f"IFS run {run:%Y-%m-%d %H} UTC, {n_steps} steps, "
          f"{len(args.products)} products, {args.workers} workers")

    if SITE_DIR.exists():                       # start clean so no stale frames survive
        shutil.rmtree(SITE_DIR)
    SITE_DIR.mkdir(parents=True)

    render_all(args.products, n_steps, args.workers)
    write_page(args.products, data)
    size = sum(f.stat().st_size for f in SITE_DIR.rglob("*") if f.is_file())
    print(f"Site written to {SITE_DIR.resolve()} ({size / 1e6:.0f} MB)")


if __name__ == "__main__":
    main()
