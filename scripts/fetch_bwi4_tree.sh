#!/usr/bin/env bash
# ===========================================================================
# FINN-bwi | robust download of the (large) /tree table via KEYSET pagination
# scripts/fetch_bwi4_tree.sh
#
# The /tree table (~1.7M rows) 500s on deep offset paging. Keyset pagination
# (order=intkey, intkey=gt.<last>) is stable and resumable: it never uses a large
# OFFSET, and it restarts from the last intkey already on disk. Safe to re-run.
# ===========================================================================
set -uo pipefail
API="https://ci.thuenen.de/rest/v1"; OUT="data/raw/bwi4"; DEST="$OUT/tree.csv"
PAGE=25000
KEY="${FINN_BWI_APIKEY:-}"; [[ -z "$KEY" && -f "$OUT/.apikey" ]] && KEY="$(tr -d ' \n\r' < "$OUT/.apikey")"
[[ -z "$KEY" ]] && { echo "no api key"; exit 1; }
ph=(-H "apikey: $KEY" -H "Accept-Profile: inventory_archive" -H "Accept: text/csv")

# intkey column index + resume point from any existing partial file
last=""; n=0
if [[ -s "$DEST" ]]; then
  ic=$(head -1 "$DEST" | tr ',' '\n' | grep -n '^intkey$' | cut -d: -f1)
  if [[ -n "$ic" ]]; then
    last=$(tail -n +2 "$DEST" | awk -F, -v c="$ic" 'NF{print $c}' | sort -n | tail -1)
    n=$(($(wc -l < "$DEST") - 1))
    echo "resuming from intkey>$last ($n rows already on disk)"
  else last=""; : > "$DEST"; fi
else : > "$DEST"; fi

while : ; do
  tmp="$(mktemp)"; ok=0
  if [[ -z "$last" ]]; then q="$API/tree?order=intkey&limit=$PAGE"
  else q="$API/tree?order=intkey&intkey=gt.$last&limit=$PAGE"; fi
  for try in 1 2 3 4 5; do
    code=$(curl -sS --max-time 600 "${ph[@]}" -o "$tmp" -w '%{http_code}' "$q")
    [[ "$code" == "200" || "$code" == "206" ]] && { ok=1; break; }
    sleep 3
  done
  (( ! ok )) && { echo ""; echo "failed at intkey>$last (HTTP $code)"; rm -f "$tmp"; exit 1; }
  rows=$(($(wc -l < "$tmp") - 1)); (( rows < 0 )) && rows=0
  if [[ ! -s "$DEST" ]]; then cat "$tmp" > "$DEST"; else tail -n +2 "$tmp" >> "$DEST"; fi
  if (( rows > 0 )); then
    ic=$(head -1 "$DEST" | tr ',' '\n' | grep -n '^intkey$' | cut -d: -f1)
    last=$(tail -n +2 "$tmp" | awk -F, -v c="$ic" 'NF{print $c}' | sort -n | tail -1)
    n=$((n + rows)); printf "\r  tree: %d rows (last intkey %s)" "$n" "$last"
  fi
  rm -f "$tmp"
  ## The server returns PAGE-1 rows per full page, so "rows < PAGE" is ALWAYS
  ## true -- break only when a keyset page comes back empty (data exhausted).
  (( rows == 0 )) && break
done
echo -e "\r  tree: done, $n rows -> $DEST"
