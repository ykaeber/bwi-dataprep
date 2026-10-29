# data/

Gitignored. Everything here is public and re-downloadable.

## Where to get it

<https://bwi.info/Download/de/BWI-Basisdaten/ACCESS2003/>

The BWI *erweiterter Grunddatensatz* is released anonymised (without exact
coordinates) and **requires no data-use agreement**. Cite it as given in
`raw/Zitiervorschlag-BWI-Basisdaten-ACCESS2003-(2015).txt`.

```bash
BASE=https://bwi.info/Download/de/BWI-Basisdaten/ACCESS2003
mkdir -p data/raw && cd data/raw
for f in bwi20150320_alle_daten2002.zip \
         bwi20150320_alle_daten2012.zip \
         bwi20150320_alle_daten2002-2012.zip \
         bwi20150320_alle_metadaten.zip \
         ShapeFile_Tnr_INSPIRE_Poly1000.zip; do
  curl -fSLO "$BASE/$f"
done
for z in bwi20150320_alle_daten2002 bwi20150320_alle_daten2012 \
         bwi20150320_alle_daten2002-2012; do unzip -o "$z.zip" -d "$z"; done
unzip -o bwi20150320_alle_metadaten.zip -d metadaten
unzip -o ShapeFile_Tnr_INSPIRE_Poly1000.zip -d inspire
```

Total ~277 MB.

## Expected layout

```
data/
├── raw/
│   ├── bwi20150320_alle_daten2002/bwi20150320_alle_daten2002.mdb
│   ├── bwi20150320_alle_daten2012/bwi20150320_alle_daten2012.mdb
│   ├── bwi20150320_alle_daten2002-2012/bwi20150320_alle_daten2002-2012.mdb
│   ├── metadaten/bwi20150320_alle_metadaten.mdb
│   ├── inspire/Tnr_INSPIRE_Poly1000.shp   (+ .dbf .prj .shx)
│   ├── bwi4/*.csv                          inventory archive (scripts/fetch_bwi4*.sh)
│   └── dwd_drought/*.asc.gz                DWD monthly grids (scripts/fetch_dwd_drought.sh)
└── bwi_csv/                                created by scripts/00_import_bwi.R
```

Override the root with `FINN_BWI_DIR`.

## Notes

- The `.mdb` files are Access 2003 databases; `scripts/00_import_bwi.R` reads
  them with **mdbtools** (`brew install mdbtools`) and caches every table as
  CSV in `bwi_csv/`. Nothing downstream needs mdbtools.
- `bwi20150320_alle_metadaten.mdb` carries the full data dictionary
  (`z_tab_phys`, `z_col_phys`) and all code tables (`x_ba`, `x_pk`, `x_Av`, …).
  It is the authority for anything schema-related — see
  [`../docs/bwi_data_notes.md`](../docs/bwi_data_notes.md).
- The 2002–2012 change database is the one FINN actually needs; the two
  single-inventory databases are used for the end-of-period tree state and for
  cross-checks.

## Latest BWI (BWI4 / 2022) — Thünen inventory-archive API

The ACCESS2003 download above stops at BWI3 (2012). **BWI4 (2022)** — and with it
the 2012→2022 transition that spans the 2018–2020 drought — is served by the
Thünen inventory-archive REST API at `ci.thuenen.de`, using a **public anonymous
key** (role `anon`) that Thünen documents at
<https://thuenen-forest-ecosystems.github.io/TFM-Documentation/api/general-api.html>.
Check Thünen's own terms / legal notice for the data licence before
redistributing — this repo does not assert one.

**Please note Thünen's usage guidance:** their Getting-Started page says the API
"should only be used for filtered queries and not for bulk data extraction", and
points to a self-hosted **SQL dump** for bulk use
(<https://git-dmz.thuenen.de/tfm-seeds/public>, server:
<https://github.com/Thuenen-Forest-Ecosystems/TFM-Server>). The `fetch_bwi4*.sh`
scripts below pull full tables for convenience; for large or repeated use prefer
the SQL dump and mind the API's rate limits.

Put the anon key in the environment or a (gitignored) key file, then download:

```bash
export FINN_BWI_APIKEY=<public anon key>     # or: echo <key> > data/raw/bwi4/.apikey
bash scripts/fetch_bwi4.sh                   # plot, regeneration, deadwood, lookups
bash scripts/fetch_bwi4_tree.sh              # the large /tree table (keyset-paged, resumable)
```

Output lands in `data/raw/bwi4/*.csv` (all intervals: `bwi2012`, `ci2017`,
`bwi2022`, `ci2027`). Then `scripts/40_import_new_schema.R` builds the BWI3→BWI4
tree panel in the same format as the ACCESS2003 pipeline — see
[`../docs/bwi_data_usage_guide.md`](../docs/bwi_data_usage_guide.md) ("Working
with the latest BWI").

## Drought (DWD)

Only for the two-decade FINN inputs (`scripts/16`, `scripts/70`). The DWD
Climate Data Center publishes 1 km monthly grids for Germany (open data,
CC BY 4.0, <https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/>):
soil moisture, potential and actual evapotranspiration, and the de Martonne index.

```bash
bash scripts/fetch_dwd_drought.sh            # 2000-2022 -> data/raw/dwd_drought/ (~1,100 files)
```

## Climate

Not stored here. WorldClim 2.1 30s `wc2.1_30s_bio` (BIO1, BIO12) and
`wc2.1_30s_elev.tif`, expected at `../worldclim` or `FINN_WORLDCLIM_DIR`.
