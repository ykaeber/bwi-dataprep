#!/usr/bin/env bash
# ===========================================================================
# FINN-bwi | download the Thuenen NFI "inventory archive" (incl. BWI4 / 2022)
# scripts/fetch_bwi4.sh
#
# The public ACCESS2003 download has only BWI2 (2002) + BWI3 (2012). BWI4 (2022)
# is served by the Thuenen inventory-archive REST API at ci.thuenen.de, CC-BY-SA
# 4.0. Access uses a PUBLIC anonymous key that Thuenen documents at
#   https://thuenen-forest-ecosystems.github.io/TFM-Documentation/api/general-api.html
# (role "anon").
#
# The key is read from $FINN_BWI_APIKEY or data/raw/bwi4/.apikey (gitignored).
# Two headers are required on every call:
#   apikey: <key>
#   Accept-Profile: inventory_archive     (PostgREST schema selector)
#
# We download each endpoint IN FULL (all intervals: bwi2012, ci2017, bwi2022,
# ci2027) as flat CSV, plus the lookup tables, then filter locally by joining
# tree.plot_id -> plot.id -> plot.interval_name. Row counts (2026-07): plot 784k,
# tree 1.72M, regeneration 382k, deadwood 494k -- a few hundred MB.
#
#   bash scripts/fetch_bwi4.sh
#   bash scripts/fetch_bwi4.sh --endpoints plot,tree        # subset
#
# Output: data/raw/bwi4/<endpoint>.csv  and  data/raw/bwi4/lookup_<name>.csv
# ===========================================================================
set -uo pipefail

API="https://ci.thuenen.de/rest/v1"
OUT="data/raw/bwi4"
ENDPOINTS=(plot tree regeneration deadwood cluster edges)
LOOKUPS=(lookup_tree_species lookup_federal_state lookup_tree_status \
         lookup_interval)
PAGE=50000                 # wide tables (plot ~70 cols) 500 above ~50k/page
mkdir -p "$OUT"

KEY="${FINN_BWI_APIKEY:-}"
[[ -z "$KEY" && -f "$OUT/.apikey" ]] && KEY="$(tr -d ' \n\r' < "$OUT/.apikey")"
if [[ -z "$KEY" ]]; then
  echo "ERROR: no API key. Set FINN_BWI_APIKEY or create $OUT/.apikey."
  echo "Public anon key: see the TFM-Documentation URL in this script's header."
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --endpoints) IFS=',' read -ra ENDPOINTS <<< "$2"; shift 2 ;;
    *) echo "unknown arg: $1"; exit 2 ;;
  esac
done

hdr=(-H "apikey: $KEY" -H "Accept-Profile: inventory_archive")

# probe
code=$(curl -sS --max-time 30 "${hdr[@]}" -o /dev/null -w '%{http_code}' \
       "$API/plot?limit=1")
[[ "$code" == "200" ]] || { echo "ERROR: key rejected (HTTP $code)."; exit 1; }
echo "[ok] key accepted."

# --- paged download of a full table as CSV -----------------------------------
# Pages by offset until the Content-Range total is reached (NOT until a short
# page: the server caps each page just under PAGE, so "short page" is every
# page). Retries transient 500s. `order=intkey` gives a stable paging order.
fetch_full() {
  local ep="$1"
  local profile="${2:-inventory_archive}"
  local dest="$OUT/$ep.csv"
  local off=0 n=0 rows total h tmp try ok
  local ph=(-H "apikey: $KEY" -H "Accept-Profile: $profile" -H "Accept: text/csv")
  : > "$dest"; total=""
  while : ; do
    tmp="$(mktemp)"; ok=0
    for try in 1 2 3; do
      h=$(curl -sS --max-time 600 "${ph[@]}" \
          -H "Range-Unit: items" -H "Range: ${off}-$((off+PAGE-1))" \
          -D - -o "$tmp" "$API/$ep?order=intkey")
      if grep -qiE '^HTTP/.* (200|206)' <<< "$h"; then ok=1; break; fi
      sleep 2
    done
    if (( ! ok )); then
      echo ""; echo "  [$ep] HTTP error at offset $off:"; grep -iE '^HTTP/' <<< "$h" | head -1
      rm -f "$tmp"; return 1
    fi
    [[ -z "$total" ]] && total=$(grep -i '^content-range:' <<< "$h" | \
                                 grep -oE '/[0-9]+' | tr -d '/\r')
    rows=$(($(wc -l < "$tmp") - 1))
    (( rows < 0 )) && rows=0
    if (( off == 0 )); then cat "$tmp" > "$dest"; else tail -n +2 "$tmp" >> "$dest"; fi
    rm -f "$tmp"
    n=$((n + rows)); off=$((off + rows))
    printf "\r  [%s] %d/%s rows" "$ep" "$n" "${total:-?}"
    (( rows < 1 )) && break
    [[ -n "$total" ]] && (( off >= total )) && break
  done
  echo -e "\r  [$ep] done: $n rows -> $dest"
}

echo "=== inventory_archive endpoints ==="
for ep in "${ENDPOINTS[@]}"; do fetch_full "$ep" inventory_archive || echo "[$ep] FAILED"; done

echo "=== lookup tables ==="
for lk in "${LOOKUPS[@]}"; do fetch_full "$lk" lookup || echo "[$lk] skipped"; done

echo
echo "Downloaded to $OUT/. Intervals present in plot.csv:"
if [[ -f "$OUT/plot.csv" ]]; then
  # interval_name column index, then tally (portable awk)
  awk -F, 'NR==1{for(i=1;i<=NF;i++)if($i=="interval_name")c=i} NR>1&&c{print $c}' \
      "$OUT/plot.csv" | sort | uniq -c
fi
echo "Next: scripts/40_import_new_schema.R builds the BWI3->BWI4 (2012->2022) tree panel."
