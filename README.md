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

## Data

The data is **not** included here — it is public, re-downloadable, and needs no
data-use agreement. See [`data/README.md`](data/README.md) for the download
commands. The `data/` directory is gitignored.

## Requirements

- R with `data.table` (plus `sf` / `terra` for the environment step).
- `mdbtools` (only for the one-time Access import in `scripts/00`).
- WorldClim 30″ rasters for `scripts/03` (see the usage guide).

## Citing the BWI

Cite the Thünen BWI *erweiterter Grunddatensatz* as given with the download (see
`data/README.md`).

## License

GPL-3.0 — see [`LICENSE`](LICENSE).
