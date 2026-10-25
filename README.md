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

## Build locally

```bash
cd Scripts
CONDA_ENV=<your-env> bash Retrieve_Forecast.bash   # or activate the env yourself first
python build_site.py              # or: --products ReTop --workers 4
python -m http.server -d ../site  # open http://localhost:8000
```

## Add a new map

1. Write `render_MyMap(frame, stamp, fc_step)` in `build_site.py` (try it in the notebook first).
2. Add an entry to `PRODUCTS` in `build_site.py`.
3. If it needs new fields or levels, add them to `PL_FIELDS` / `SFC_FIELDS`
   in `Retrieve_ECMWF_OpenData.py`.
