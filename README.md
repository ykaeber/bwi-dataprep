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
- **Both decades together** — one script, `scripts/71_prepare_two_interval.R`,
  stitches 2002→2012 and 2012→2022 into one FINN input set,
  `two-interval-v13-soil`: per-year drought (DWD grids via
  `scripts/fetch_dwd_drought.sh` + `scripts/16_drought_index.R`), static soil
  covariates, and spatially blocked train/holdout sets from a fixed split table
  (`inputs/split_tracts.csv`). Its decade-2 tree panel comes from
  `scripts/91_import_new_schema_chained.R` (the `scripts/40` panel with the
  fates fixed). See [Status](#status-relative-to-the-finn-bwi-study).

The data itself is **not** included here — it is public and re-downloadable. See
[`data/README.md`](data/README.md) for the download commands for both releases.
The `data/` directory is gitignored.

## Requirements

- R with `data.table` (plus `sf` / `terra` for the environment step).
- For the two-decade FINN inputs: the `FINN` R package, version 0.2.0 or later
  (`remotes::install_github("FINNverse/FINN")` until the CRAN release; it has
  `makeObsData()`, `resolveSiteIDs()` and the per-tree mortality loss), `FNN`,
  `sf` (the split in `scripts/71`), and `ggplot2` / `patchwork` (figures in
  `scripts/16`).
- `mdbtools` (the one-time Access import in `scripts/00`, and the soil
  database in `scripts/71`).
- WorldClim 30″ rasters for `scripts/03` (see the usage guide).

## Status relative to the FINN-bwi study

The scripts here are copies of the data-preparation steps of the FINN-bwi study
(FINN fitted to the BWI). `scripts/00`–`04`, `16`, `70` and `91` and the
helpers in `R/` match it byte for byte as of 2026-11-09; `scripts/40` and
`scripts/fetch_bwi4.sh` differ only in comments and messages. (Comments in
`scripts/91` mention study scripts that are not part of this repository.)

**`two-interval-v13-soil` (`scripts/71`) is the two-decade set the study's final
model is fitted on.** The study builds it in three steps (v8, then v12 without
the year-0 env row, then v13-soil with soil); here one script does all of it.
The earlier chain of this repository (`scripts/70` v3, `71` v8, `72` v12, `73`
v13-soil) is in the git history. Against a plain inventory panel the set makes
six decisions (numbered as in the script header):

1. **Freshly dead trees are not survivors.** The BWI records a tree dead for
   less than about 12 months like a live sample tree, with a "freshly dead"
   flag (`tot` in the ACCESS2003 `b2_wzp` / `b3_wzp` tables, `damage_dead` in
   the archive's `tree` table). Dead at the start: not in the start stand; dead
   at the end: a natural death; a new tree already dead: not a recruit.
2. **2012–2017 losses at the carbon-inventory corners are kept.** Where a corner
   was remeasured in the 2017 carbon inventory (CI2017), BWI2022 refers to
   CI2017 as the previous survey, so trees lost between 2012 and 2017 have no
   2022 record; `scripts/91` chains their fate through CI2017.
3. **Lateral ingrowth is counted once.** The start stand is the start tally
   only (no back-modelled 2002 state of laterally ingrown trees on top of it),
   and the end state is the living end tally expanded afresh, so a laterally
   ingrown tree is not counted twice.
4. **One drought index**, `drought2`, the mean of the current and the previous
   year's soil-moisture anomaly.
5. **No env row for the start census.** FINN has no initial-state env row (the
   initial state is `init_cohort_dt`) and every env row drives one annual
   transition, so the env tables hold years 1..10 (calendar t0+1..t0+10; simulated
   year k is calendar t0 + k).
6. **Static soil covariates** per corner from the BWI 2012 environmental
   database (available water capacity, perched-water, groundwater and
   waterlogging degree, nutrient class; gaps filled by natural forest community
   x state; download in `data/README.md`, details in `docs/bwi_data_notes.md`
   section 11).

Earlier versions of the study and of this repository lacked items 1–3:

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

Item 1 is in `scripts/91` (decade 2) and `scripts/71` (decade 1), item 2 in
`scripts/91`, items 3–6 in `scripts/71`. Item 3 leaves 617 decade-1
site-periods without any start tree, so the set has 40,787 site-periods (v3:
41,285; items 1 and 2 alone give 41,404). The tree-level table of
`bwi_tree_panel()` (ACCESS2003, 2002→2012) still lacks item 1 (it does not
read the `tot` flag; `scripts/71` applies it on top) but splits lateral
ingrowth from recruits correctly.

**The split is fixed data.** `inputs/split_tracts.csv` (in git) assigns every
tract its block (100 km on EPSG:3035), set (`train` / `holdout_blocked`) and
`train_quarter` / `pilot` membership; the holdout blocks were drawn once with
seed 42. It is read, not recomputed, so the sets stay those the study's models
were evaluated on.

**Reproducibility.** `outputs/prepared/two-interval-v13-soil.sha256` (in git)
holds the SHA-256 of the study's files. Run against the same raw data and
upstream outputs, `scripts/91` + `scripts/71` reproduce every FINN table of
the study's v13-soil byte for byte. `split_blocks.csv` differs by design (this
set's 12,388 tracts instead of the study's 12,443, with the same block, set,
`train_quarter` and `pilot` on every tract), and the study's four comparison
files against its intermediate versions are not written. The usage guide shows
the check.

**The regeneration target has a known error.** `obs_dt$reg` converts the stock
of 5.0–6.9 cm stems on the regeneration circle into a flux at an assumed
diameter increment of 0.30 cm/yr. The script takes the 2002 circle
(r = 1.75 m) and a 1.9 cm band; from BWI 2012 on the circle is r = 2 m and the
band is 2.0 cm wide, so the target is 1.375 times too high. It stays because the
study's final model was fitted on it. Corrected targets (geometry only; the
circle pattern scaled to the recruits of the angle-count sample) are being
refitted, and the chosen one will become the next version of this script.

## Citing the BWI

Cite the Thünen BWI *erweiterter Grunddatensatz* as given with the download (see
`data/README.md`).

## License

GPL-3.0 — see [`LICENSE`](LICENSE).
