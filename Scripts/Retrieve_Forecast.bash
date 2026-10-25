#!/bin/bash
set -euo pipefail

# Run from the script's own folder, so this also works from cron
cd "$(dirname "$(readlink -f "$0")")"

# Optionally activate a conda environment: CONDA_ENV=myenv bash Retrieve_Forecast.bash
# Without CONDA_ENV the current Python is used (as on GitHub Actions).
if [[ -n "${CONDA_ENV:-}" ]]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi

# Download the IFS forecast (GRIB) and the ECMWF meteogram
python Retrieve_ECMWF_OpenData.py

# Download the DWD analysis and forecast charts
mkdir -p ../Data/DWD_analysis
pushd ../Data/DWD_analysis > /dev/null

# A failed DWD chart only prints a warning; the site is still built without it.
BASE=https://www.dwd.de/DWD/wetter/wv_spez/hobbymet/wetterkarten
get() { curl -fsS "$BASE/$1" -o "$2" || { echo "WARNING: could not download $1"; rm -f "$2"; }; }

get bwk_bodendruck_na_ana.png Analysis.png
get ico_tkboden_na_v36.png    Analysis_p024.png
get ico_tkboden_na_036.png    Analysis_p036.png
get ico_tkboden_na_048.png    Analysis_p048.png
get ico_tkboden_na_060.png    Analysis_p060.png
get ico_tkboden_na_084.png    Analysis_p084.png
#get ico_tkboden_na_108.png   Analysis_p108.png

popd > /dev/null
