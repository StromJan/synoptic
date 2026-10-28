"""Download the IFS fields the maps need, plus the ECMWF meteogram.

    python Retrieve_ECMWF_OpenData.py            # newest complete 00/12 UTC run
    IFS_RUN=2026102800 python Retrieve_ECMWF_OpenData.py   # a specific run
    python Retrieve_ECMWF_OpenData.py --latest   # only print the newest run, e.g. 2026102800
"""
from collections import defaultdict
from datetime import datetime
from pathlib import Path
import os
import sys

import requests
from ecmwf.opendata import Client


# ---------------------------------------------------------------------------
# What the Synoptic_Overview notebook actually uses.
# Add a field/level here when a new map needs it; nothing else has to change.
# ---------------------------------------------------------------------------
STEPS = list(range(0, 169, 6))          # 0, 6, ..., 168 h

# Pressure-level fields: {param: [levels in hPa]}
PL_FIELDS = {
    "gh": [500, 1000],   # ReTop: 500/1000 hPa thickness and 500 hPa height
    "t":  [850],         # Frontogenesis, EqPotTemperature
    "q":  [850],         # EqPotTemperature
    "u":  [200, 850],    # JetDiagnose (200), Frontogenesis (850)
    "v":  [200, 850],
}

# Surface fields
SFC_FIELDS = ["2t", "msl", "tcwv"]

DATA_DIR = Path(__file__).resolve().parent / "../Data"
METEOGRAM_DIR = DATA_DIR / "meteograms"


client = Client(
    source="ecmwf",
    model="ifs",
    resol="0p25",
    preserve_request_order=False,
    retry_after=(10, 300, 2),
)


def latest_run():
    """Latest 00/12 UTC run that already has the last step we need.

    Pinning the run once means every file (and the meteogram) comes from
    the same forecast, even if a new run appears while downloading.
    """
    return client.latest(stream="oper", type="fc", step=STEPS[-1], param="msl")


def group_by_levels(fields):
    """{param: levels} -> {levels: [params]}.

    One request with param=[gh, t] and levelist=[500, 850, 1000] would fetch
    every combination (6 fields). Grouping params that share the same levels
    fetches exactly the fields in PL_FIELDS.

    Each group also goes into its own file, because earthkit's xarray engine
    refuses to put variables with different level sets into one dataset.
    """
    groups = defaultdict(list)
    for param, levels in fields.items():
        groups[tuple(sorted(levels))].append(param)
    return groups


def retrieve_file(run, request, target, expected):
    """Download one request to `target` and check the number of fields.

    Writes via a temporary name, so the notebook never reads a half file.
    """
    target = Path(target)
    tmp = target.with_name(target.name + ".tmp")
    client.retrieve(date=run, stream="oper", type="fc", step=STEPS,
                    target=str(tmp), **request)

    import eccodes                      # imported here so --latest works without eccodes

    with open(tmp, "rb") as f:
        n = eccodes.codes_count_in_file(f)
    if n != expected:
        raise RuntimeError(f"{target.name}: expected {expected} fields, got {n}")

    os.replace(tmp, target)
    print(f"Saved {target.name}: {', '.join(request['param'])} "
          f"({n} fields, {target.stat().st_size / 1e6:.0f} MB)")
    return target


def download_forecast(run):
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    written = set()

    for levels, params in group_by_levels(PL_FIELDS).items():
        name = "pl_" + "-".join(str(l) for l in levels) + ".grib"   # e.g. pl_200-850.grib
        req = dict(levtype="pl", levelist=list(levels), param=params)
        written.add(retrieve_file(run, req, DATA_DIR / name,
                                  len(STEPS) * len(levels) * len(params)))

    req = dict(levtype="sfc", param=SFC_FIELDS)
    retrieve_file(run, req, DATA_DIR / "sfc.grib", len(STEPS) * len(SFC_FIELDS))

    # Remove level-group files left over from an older PL_FIELDS setup,
    # otherwise the notebook would load stale data from a previous run.
    for old in DATA_DIR.glob("pl_*.grib"):
        if old not in written:
            old.unlink()
            print("Removed stale", old.name)


# ---------------------------------------------------------------------------
# ECMWF meteogram for Helsinki
# ---------------------------------------------------------------------------
API = "https://charts.ecmwf.int/opencharts-api/v1/products/opencharts_meteogram/"


def fetch_meteogram(base_time):
    METEOGRAM_DIR.mkdir(parents=True, exist_ok=True)
    params = {
        "base_time": base_time.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "lat": 60.1695,
        "lon": 24.9354,
        "station_name": "Helsinki",
    }
    r = requests.get(API, params=params, timeout=30)
    r.raise_for_status()
    img_url = r.json()["data"]["link"]["href"]

    img = requests.get(img_url, timeout=60)
    img.raise_for_status()
    path = METEOGRAM_DIR / "helsinki_meteogram.png"
    path.write_bytes(img.content)
    print("Saved", path.name)
    return path


if __name__ == "__main__":
    if "--latest" in sys.argv:              # used by the GitHub workflow's check job
        print(f"{latest_run():%Y%m%d%H}")
        sys.exit()

    if os.environ.get("IFS_RUN"):           # the run the check job found
        run = datetime.strptime(os.environ["IFS_RUN"], "%Y%m%d%H")
    else:
        run = latest_run()
    print(f"Using IFS run {run:%Y-%m-%d %H} UTC")
    download_forecast(run)
    try:                                    # a missing meteogram shouldn't stop the site build
        fetch_meteogram(run)
    except Exception as e:
        print("WARNING: meteogram download failed:", e)
