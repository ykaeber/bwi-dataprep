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
  same columns as the ACCESS2003 panel, so the tree-level and plot-level
  products run unchanged. One difference in meaning: the archive has no
  back-modelled diameters, so `fate == "recruit"` in the 2012→2022 panel is
  *every* new sample tree, lateral ingrowth included (see the usage guide).
- **Both decades together** — `scripts/70_prepare_two_interval.R` stitches
  2002→2012 and 2012→2022 into one FINN input set (`two-interval-v3`) with a
  per-year drought covariate (DWD grids via `scripts/fetch_dwd_drought.sh` +
  `scripts/16_drought_index.R`) and spatially blocked train/holdout sets. See
  [Status](#status-relative-to-the-finn-bwi-study) for how this set differs
  from the one the BWI study now uses.

The data itself is **not** included here — it is public and re-downloadable. See
[`data/README.md`](data/README.md) for the download commands for both releases.
The `data/` directory is gitignored.

## Requirements

- R with `data.table` (plus `sf` / `terra` for the environment step).
- For the two-decade FINN inputs: the `FINN` R package, version 0.2.0 or later
  (`remotes::install_github("FINNverse/FINN")` until the CRAN release; it has
  `makeObsData()`, `resolveSiteIDs()` and the per-tree mortality loss), `FNN`,
  and `ggplot2` / `patchwork` (figures in `scripts/16`).
- `mdbtools` (only for the one-time Access import in `scripts/00`).
- WorldClim 30″ rasters for `scripts/03` (see the usage guide).

## Status relative to the FINN-bwi study

The scripts here are copies of the data-preparation steps of the FINN-bwi study
(FINN fitted to the BWI). `scripts/00`–`04`, `16` and `70` and the helpers in
`R/` match it byte for byte as of 2026-11-09; `scripts/40` and
`scripts/fetch_bwi4.sh` differ only in comments and messages.

The two-decade set this repository writes, `two-interval-v3`, is **not** the
set the study now fits. The study has since moved to `two-interval-v8`, which
corrects three defects that `scripts/40` and `scripts/70` here still have:

1. **Freshly dead trees counted as survivors.** The BWI records a tree dead for
   less than about 12 months like a live sample tree, with a "freshly dead"
   flag (`tot` in the ACCESS2003 `b2_wzp` / `b3_wzp` tables, `damage_dead` in
   the archive's `tree` table). v3 ignores the flag in both decades.
2. **2012–2017 losses missing at the carbon-inventory corners.** Where a corner
   was remeasured in the 2017 carbon inventory (CI2017), BWI2022 refers to
   CI2017 as the previous survey, so trees lost between 2012 and 2017 have no
   2022 record. v3 takes every decade-2 fate from the 2022 records, so these
   losses drop out of both the fates and the 2012 start stand.
3. **Lateral ingrowth counted twice in the stand state.** Decade 1 adds the
   back-modelled 2002 state of laterally ingrown trees to the 2002 start stand,
   on top of the 2002 tally that already represents them (initial basal area
   about 1.26 times BWI2's own 2002 tally). The end state carries the start
   stand forward and adds every new tree, so a laterally ingrown tree counts
   twice there in both decades (end basal area 1.34 and 1.20 times the design
   estimate). In decade 2 this goes with the `recruit` meaning above.

Fixes 1 and 2 arrived in the study's `two-interval-v5`, fix 3 in
`two-interval-v8`. They have not been ported here yet. Until they are, use
`two-interval-v3` with these defects in mind. The tree-level table of
`bwi_tree_panel()` (ACCESS2003, 2002→2012) has defect 1 (it does not read the
`tot` flag) but splits lateral ingrowth from recruits correctly; the
regeneration target from the fixed-area circles is not affected by any of the
three.

## Citing the BWI

Cite the Thünen BWI *erweiterter Grunddatensatz* as given with the download (see
`data/README.md`).

## License

GPL-3.0 — see [`LICENSE`](LICENSE).
