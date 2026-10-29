#!/usr/bin/env bash
# ===========================================================================
# FINN-bwi | DWD 1 km monthly drought grids for scripts/16
# scripts/fetch_dwd_drought.sh [START_YEAR] [END_YEAR]   (default 2000 2022)
#
# Both BWI decades (2002->2012, 2012->2022) need the grids for 2000-2022, so
# the drought covariate carries the 2003 heatwave and the 2018-2020 drought.
# Downloads the four grid types scripts/16 uses; files already present and
# valid are skipped, so re-running only fills gaps.
#
# DWD open data (CC-BY-4.0) has TWO directory layouts, both handled here:
#   drought_index/            -> month subfolders  01_Jan/ .. 12_Dec/
#   soil_moist/ evapo_p/ evapo_r/ -> flat directory
#
#   bash scripts/fetch_dwd_drought.sh
# ===========================================================================
set -uo pipefail

BASE="https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly"
DEST="data/raw/dwd_drought"
START="${1:-2000}"; END="${2:-2022}"
TYPES=(soil_moist evapo_p evapo_r drought_index)
MONTHDIR=(01_Jan 02_Feb 03_Mar 04_Apr 05_May 06_Jun \
          07_Jul 08_Aug 09_Sep 10_Oct 11_Nov 12_Dec)
mkdir -p "$DEST"

get() {  # $1 url  $2 dest ; validates gzip, retries twice
  local url="$1" out="$2" try code
  for try in 1 2 3; do
    code=$(curl -sL --max-time 90 -o "$out" -w '%{http_code}' "$url")
    if [[ "$code" == "200" ]] && gzip -t "$out" 2>/dev/null; then return 0; fi
    rm -f "$out"; sleep 1
  done
  return 1
}

total=0; ok=0; miss=0
for t in "${TYPES[@]}"; do
  for y in $(seq "$START" "$END"); do
    for i in $(seq 0 11); do
      mm=$(printf "%02d" $((i+1)))
      f="grids_germany_monthly_${t}_${y}${mm}.asc.gz"
      total=$((total+1))
      [[ -s "$DEST/$f" ]] && gzip -t "$DEST/$f" 2>/dev/null && { ok=$((ok+1)); continue; }
      if [[ "$t" == "drought_index" ]]; then url="$BASE/$t/${MONTHDIR[$i]}/$f"
      else url="$BASE/$t/$f"; fi
      if get "$url" "$DEST/$f"; then ok=$((ok+1)); else miss=$((miss+1)); fi
    done
    printf "\r%-14s %d: ok=%d miss=%d" "$t" "$y" "$ok" "$miss"
  done
  echo ""
done

echo "Done: $ok/$total present, $miss missing."
for t in "${TYPES[@]}"; do
  n=$(ls "$DEST"/grids_germany_monthly_${t}_20*.asc.gz 2>/dev/null | wc -l | tr -d ' ')
  yr=$(ls "$DEST"/grids_germany_monthly_${t}_*.asc.gz 2>/dev/null | \
       sed -E 's/.*_([0-9]{4})[0-9]{2}\.asc\.gz/\1/' | sort -u | \
       { read a; tac 2>/dev/null | head -1 | xargs -I{} echo "$a-{}"; } 2>/dev/null)
  printf "  %-14s %s files\n" "$t" "$n"
done
echo "Next: Rscript scripts/16_drought_index.R, then scripts/70 for the two-decade inputs."
