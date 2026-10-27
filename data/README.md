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
│   └── inspire/Tnr_INSPIRE_Poly1000.shp   (+ .dbf .prj .shx)
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

## Climate

Not stored here. WorldClim 2.1 30s `wc2.1_30s_bio` (BIO1, BIO12) and
`wc2.1_30s_elev.tif`, expected at `../worldclim` or `FINN_WORLDCLIM_DIR`.
