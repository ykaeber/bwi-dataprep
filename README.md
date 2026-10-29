# bwi-dataprep

Reproducible tools and a step-by-step guide for turning the **German National
Forest Inventory** (Bundeswaldinventur, BWI) public release into analysis-ready
tables.

The BWI is an angle-count (Winkelzählprobe) sample inventory. Its raw release is
powerful but awkward: there is no fixed plot area, "new sample tree" does not
mean "recruit", and harvest has to be separated from natural mortality. This
repository packages the corrections once, so you can get straight to three
products:

1. **A tree-level long table** — one row per tallied tree × inventory, with
   edge-corrected per-hectare densities, the lateral-ingrowth split, and clean
   fate codes (survivor / natural death / harvest / recruit).
2. **Plot-level aggregates** — stand variables (basal area, stem density,
   quadratic mean diameter) and species composition per sample point, derived
   from the tree table with the correct density weighting.
3. **Inputs for dynamic forest models** — an initial stand state plus site
   drivers, ready for cohort/gap models (worked through with
   [FINN](https://github.com/FINNverse/FINN); the same ingredients feed ForClim
   and relatives).

## Start here

- **[`docs/bwi_data_usage_guide.md`](docs/bwi_data_usage_guide.md)** — the
  step-by-step how-to for all three products, with runnable snippets.
- **[`docs/bwi_data_notes.md`](docs/bwi_data_notes.md)** — the *why*: sampling
  design, the ingrowth trap, fate codes, selection-bias corrections, coordinates
  and species, each verified against the actual release.

## Quick start

```r
# from the repository root, after downloading the data (see data/README.md)
source("R/config.R"); source("R/bwi_io.R"); source("R/wzp.R"); source("R/bwi_trees.R")

panel <- bwi_tree_panel()$trees        # the tree-level long table
```

Or run the scripted pipeline end to end:

```bash
Rscript scripts/00_import_bwi.R          # raw Access .mdb -> data/bwi_csv/*.csv (needs mdbtools)
Rscript scripts/01_prepare_metadata.R    # corner population + descriptive tables
Rscript scripts/02_select_calibration_set.R
Rscript scripts/03_attach_environment.R  # WorldClim climate + elevation + coordinates
Rscript scripts/04_prepare_finn_inputs.R # model-ready tables
```

## Data releases covered

- **ACCESS2003** — BWI2 (2002) → BWI3 (2012); the public *erweiterter
  Grunddatensatz*, no data-use agreement. Pipeline: `scripts/00`–`04`.
- **Thünen inventory-archive API** — adds **BWI4 (2022)**, i.e. the
  BWI3 → BWI4 (2012 → 2022) transition over the 2018–2020 drought; public
  anonymous key (check Thünen's terms for the data licence). Download with
  `scripts/fetch_bwi4*.sh` — which also note Thünen's "filtered queries, not bulk"
  guidance — then build the tree panel with `scripts/40_import_new_schema.R`:
  same format as the ACCESS2003 pipeline, so the tree-level and plot-level
  products apply unchanged.
- **Both decades together** — `scripts/70_prepare_two_interval.R` stitches
  2002→2012 and 2012→2022 into one FINN input set with a per-year drought
  covariate (DWD grids via `scripts/fetch_dwd_drought.sh` +
  `scripts/16_drought_index.R`) and spatially blocked train/holdout sets.

The data itself is **not** included here — it is public and re-downloadable. See
[`data/README.md`](data/README.md) for the download commands for both releases.
The `data/` directory is gitignored.

## Requirements

- R with `data.table` (plus `sf` / `terra` for the environment step).
- For the two-decade FINN inputs: the `FINN` R package, `FNN`, and `ggplot2` /
  `patchwork` (figures in `scripts/16`).
- `mdbtools` (only for the one-time Access import in `scripts/00`).
- WorldClim 30″ rasters for `scripts/03` (see the usage guide).

## Citing the BWI

Cite the Thünen BWI *erweiterter Grunddatensatz* as given with the download (see
`data/README.md`).

## License

GPL-3.0 — see [`LICENSE`](LICENSE).
