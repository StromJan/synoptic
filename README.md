# Synoptic overview

IFS forecast maps for the North Atlantic, rebuilt twice a day by GitHub Actions
and published on GitHub Pages.

```
Scripts/
  Retrieve_ECMWF_OpenData.py   download only the IFS fields the maps need (+ meteogram)
  Retrieve_Forecast.bash       runs the above, then fetches the DWD charts
  build_site.py                renders all maps in parallel and writes ../site/
  site_template.html           page layout and viewer (the build fills in the data)
  Synoptic_Overview.ipynb      notebook for trying out new maps
Data/                          downloads (not committed)
site/                          generated website (not committed)
requirements.txt               packages for the GitHub build
.github/workflows/pages.yml    schedule: download -> build -> publish
```

## When it updates

A quick check runs every 10 min during 05–12 and 17–24 UTC (every 2 h otherwise).
It compares the newest IFS run that is complete to +168 h with `img/run.txt` on
the live site, and only rebuilds when they differ. Pushing changes to `Scripts/`
or pressing **Run workflow** always rebuilds.

## Build locally

```bash
cd Scripts
CONDA_ENV=<your-env> bash Retrieve_Forecast.bash   # or activate the env yourself first
python build_site.py              # or: --products ReTop --workers 4
python -m http.server -d ../site  # open http://localhost:8000
```

While working on the page or one map, skip re-rendering what's already there:

```bash
python build_site.py --local              # reuse existing maps, render only missing ones
python build_site.py --local --redo ReTop # ...and re-render ReTop after changing it
```

Maps are only reused if they come from the same IFS run as `Data/`; after a new
download everything is rendered again.

## Add a new map

1. Write `render_MyMap(frame, stamp, fc_step)` in `build_site.py` (try it in the notebook first).
2. Add an entry to `PRODUCTS` in `build_site.py`.
3. If it needs new fields or levels, add them to `PL_FIELDS` / `SFC_FIELDS`
   in `Retrieve_ECMWF_OpenData.py`.
