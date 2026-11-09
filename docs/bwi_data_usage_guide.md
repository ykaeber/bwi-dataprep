# Preparing and using the BWI data — a step-by-step guide

This guide shows how to turn the raw German National Forest Inventory
(Bundeswaldinventur, BWI) release into three ready-to-use products:

1. **a tree-level long table** — one row per tallied tree × inventory;
2. **plot-level aggregates** — stand variables and species composition per
   sample point, derived from the tree table;
3. **inputs for a dynamic forest model** — an initial stand state plus site
   drivers (worked through with FINN; the same ingredients feed ForClim and
   other cohort/gap models).

It is grounded in the scripts and helpers in this repository, but products 1–2
are general-purpose: you do not need to care about FINN to use them. For *why*
each step is done the way it is — the sampling design, the angle-count
expansion, the "ingrowth trap", the fate codes and the selection-bias
corrections — see the companion **`docs/bwi_data_notes.md`**, which verifies
every decision against the actual Thünen release. This guide is the *how*.

All snippets assume you start from the project root and use `data.table`.

---

## Step 0 — One-time: raw Access databases → CSV

The public release ships as Microsoft Access `.mdb` files. `scripts/00_import_bwi.R`
exports every needed table to `data/bwi_csv/*.csv` using `mdbtools`; everything
downstream reads the CSVs, so `mdbtools` is needed **only** for this step.

```bash
Rscript scripts/00_import_bwi.R      # writes data/bwi_csv/*.csv
```

In this repository the CSVs are already present, so you can skip straight to
Step 1. From then on you read a cached table with the `bwi()` helper:

```r
source("R/config.R"); source("R/bwi_io.R")
trees_raw <- bwi("b3_baeume")          # end-of-period tree table, as a data.table
species   <- bwi_species_table()       # tidy Ba -> binomial lookup (fixes misspellings)
```

Key facts baked into the constants (`R/config.R`): angle-count factor
`WZP_BAF = 4` m²/ha, caliper threshold `WZP_MIN_BHD = 7` cm, remeasurement
period `BWI_PERIOD_YEARS = 10`, and the diameter→basal-area helper
`bhd2ba(dbh_cm) = π·(dbh/200)²` (m²).

---

## Step 1 — The tree-level long table

**One function does the whole assembly:**

```r
source("R/config.R"); source("R/bwi_io.R"); source("R/wzp.R"); source("R/bwi_trees.R")

panel <- bwi_tree_panel()              # list(trees = ..., corner_flags = ...)
trees <- panel$trees                   # the long table, one row per (tree, date)
```

`bwi_tree_panel()` is doing the hard, BWI-specific work for you — the parts that
are easy to get wrong if you join the raw tables by hand:

- **Angle-count density.** Each tallied stem carries `n_ha` = stems/ha it
  represents (`K/g`), taken from the inventory's own `N_ha` column so the
  Grenzkreis (edge) correction is preserved.
- **The ingrowth trap.** A `Pk == 0` ("new sample tree") code does *not* mean a
  recruit: ~91 % of them are established trees that merely grew into the
  angle-count inclusion zone. The panel splits on the back-modelled start
  diameter — `< 7 cm` → genuine `recruit`, `≥ 7 cm` → lateral ingrowth present at
  both dates — instead of taking `Pk == 0` at face value.
- **Fate codes.** `fate` separates `survivor`, natural death (`dead`), `harvest`
  and `recruit`, so harvest is never mistaken for mortality.
- **Field dates.** The realised remeasurement interval per tract (8–12 yr) comes
  from the field dates, not a nominal 10.

### Columns

| column | meaning | unit |
|---|---|---|
| `Tnr, Enr, Bnr` | tract / corner / tree key | — |
| `treeKey` | unique tree id across both dates | — |
| `siteName` | tract id (the "site") | — |
| `patchName` | corner within the tract (the sample point) | — |
| `time` | period index: `0` = start (BWI2), `1` = end (BWI3) | — |
| `year` | calendar year of the measurement | yr |
| `species_name` | binomial (misspellings corrected) | — |
| `dbh` | diameter at breast height | cm |
| `n_ha` | stems/ha this tallied stem represents | ha⁻¹ |
| `status` | `alive` / `new` (recruit at t1) / `dead` | — |
| `living` | `status %in% c("alive","new")` | logical |
| `mort` | became dead this period (was alive, now dead) | logical |
| `fate` | `survivor` / `dead` / `harvest` / `recruit` / ingrowth | — |
| `is_recruit` | genuine recruit (back-modelled dbh < 7 cm) | logical |

### Common filters

```r
library(data.table)
trees[living == TRUE & time == 1L]          # standing trees at BWI3
trees[fate == "dead"]                        # natural mortality events
trees[fate == "harvest"]                     # removals (not mortality)
trees[is_recruit == TRUE]                    # genuine recruits (appear only at t1)
```

### Quality control

`panel$corner_flags` lists corners carrying unclassifiable trees. They are
flagged rather than silently dropped, because removing a live tree would read
downstream as mortality that never happened. Exclude flagged corners when you
need clean demographic rates.

---

## Step 2 — Plot-level aggregates (derived from the tree table)

**The one rule that matters:** the angle-count sample has **no fixed plot
area**, so you aggregate by **weighting each stem with its `n_ha` density**, not
by counting rows. Summing `n_ha` gives stems/ha; summing `n_ha · bhd2ba(dbh)`
gives basal area per hectare.

### 2a — Stand variables per sample point × date

```r
stand <- trees[living == TRUE, .(
  ba_ha     = sum(n_ha * bhd2ba(dbh)),                 # basal area, m²/ha
  stems_ha  = sum(n_ha),                               # stem density, ha⁻¹
  dg_cm     = sqrt(sum(n_ha * dbh^2) / sum(n_ha)),     # quadratic mean diameter
  dbh_mean  = weighted.mean(dbh, n_ha),                # density-weighted mean dbh
  n_species = uniqueN(species_name)
), by = .(siteName, patchName, time, year)]
```

### 2b — Species composition per sample point × date

```r
comp <- trees[living == TRUE, .(
  ba_ha    = sum(n_ha * bhd2ba(dbh)),
  stems_ha = sum(n_ha)
), by = .(siteName, patchName, time, species_name)]
comp[, ba_share := ba_ha / sum(ba_ha), by = .(siteName, patchName, time)]
```

Aggregate to the **tract** (`siteName`) or to **regions/Bundesländer** by joining
the corner metadata (`scripts/01_prepare_metadata.R`, `R/regions.R`) and
averaging the per-corner stand table.

### 2c — Demographic rates (growth, mortality, recruitment)

Rates need two things the stock aggregates above do not, and both have documented
traps (see `docs/bwi_data_notes.md` §5b–§6):

- **Use tallied counts, not expanded densities.** One *measured* tree must not
  tell a likelihood it is forty. Stocks (BA/ha, N/ha) use `n_ha`; rates and
  sample sizes (`n_at_risk`, `n_died`, growth `n`) use the **count of tallied
  stems**.
- **Competing risks + immortal-time bias.** Natural death and harvest are
  separate removals on the same at-risk cohort; harvest is *censored* out of the
  natural-mortality denominator. Lateral ingrowth must be **excluded** from the
  at-risk cohort — it is conditioned on having survived, so leaving it in
  understates mortality (it is 27 % of the naive cohort, concentrated in small
  trees).

This logic is already encoded — do not re-derive it by hand:

```r
# scripts/04 builds a tallied table `tal` with columns n_ha, at_risk, died,
# growth_ok (+ harvested, rel_growth), then:
eff <- wzp_effective_n(tal, by = c("siteName", "patchName", "species_name"))
#   -> n_at_risk, n_died_eff, n_harvest_eff, growth_n  (tallied counts, censored)
```

Read `scripts/04_prepare_finn_inputs.R` (the block around `wzp_effective_n`) for
the reference implementation, and `wzp_effective_n()` in `R/wzp.R` for the
competing-risk / censoring details.

---

## Step 3 — Inputs for a dynamic forest model

Any cohort- or gap-based dynamic forest model (FINN, ForClim, and relatives)
needs the same three ingredients, all derivable from Steps 1–2:

1. **an initial stand state** — the trees present at one date, per simulated
   patch;
2. **site drivers** — climate, elevation, soil, and any disturbance forcing;
3. *(for calibration only)* **demographic observations** — the growth, mortality
   and recruitment of Step 2c, to fit against.

### 3a — Initial stand state from the tree table

Take one date's living trees per sample point. Models that simulate a
**fixed-area patch** want integer tree counts, so expand the `n_ha` densities
onto a patch of known area (default `PATCH_SIZE_HA = 0.1` ha), preserving group
totals with a largest-remainder allocation:

```r
init <- trees[living == TRUE & time == 1L,
              .(treeKey, siteName, patchName, species_name, n_ha, dbh)]
wzp_allocate(init, by = c("patchName", "species_name"),
             patch_size_ha = PATCH_SIZE_HA)   # adds integer n_trees (by reference)
init_cohorts <- wzp_expand(init)              # replicate to one row per virtual tree
```

Models that work in **per-hectare densities** (many matrix/gap models) can skip
the expansion and use the `n_ha`-weighted species × dbh-class table from Step 2
directly.

### 3b — Site drivers

`scripts/03_attach_environment.R` attaches the environment to each corner:
WorldClim climate normals (`mat_c`, `map_mm`), the inventory's own per-tract
elevation (`elev_m`), anonymised coordinates (`x`, `y`), and — via
`scripts/16_drought_index.R` — a per-year growing-season soil-moisture anomaly.
Coordinates are snapped to the 1 km INSPIRE grid (positional error ≲ 700 m; the
coordinate chain is validated to r = 0.99 against inventory elevation).

Note the resolution mismatch between model families: FINN consumes **annual**
drivers with climate normals, whereas **ForClim needs monthly** temperature and
precipitation plus a soil water-holding capacity and latitude. For a
monthly-climate model, derive monthly normals from WorldClim (`wc2.1_30s_tavg`,
`_prec`) at the same coordinates `scripts/03` already resolves, and supply a
bucket size from a soil map rather than from the inventory.

### 3c — FINN (worked example)

FINN consumes four objects, built by running the pipeline end to end:

```bash
Rscript scripts/02_select_calibration_set.R   # apply filters -> calibration corners
Rscript scripts/03_attach_environment.R       # attach climate / elevation / drought
Rscript scripts/04_prepare_finn_inputs.R      # -> outputs/finn/{obs_dt,env_dt,init_cohort_dt,species_dt}.csv
```

- **`obs_dt`** — one row per site × year × species: `ba, dbh, trees, growth,
  mort, n_at_risk, n_died, reg, species`. The initial state (year 0) is *not* in
  `obs_dt`. Units: `ba` m², `dbh` cm, `reg` trees ha⁻¹; mortality is the count
  pair `(n_at_risk, n_died)` the binomial loss consumes.
- **`env_dt`** — `siteID, year, <drivers>` in natural units, covering every
  site × year (FINN z-standardises internally; do not pre-scale).
- **`init_cohort`** — the year-0 stand (Step 3a), via `makeInitCohorts`.
- **`species_dt`** — contiguous `1..N` species ids aligned across tables.

Two BWI-specific choices live in `scripts/04`: the single observation is placed
at year `P` with a `period_length` column (so fitted rates are per year, not per
decade), and `growth`/`mort` are annualised to match FINN's mean-of-annual-rates
scoring. The data-prep contract is summarised in the `C-Data_preparation`
vignette of the FINN package and in the `finn-dataprep` notes of the FINN
repository
([`.claude/skills/finn-dataprep/SKILL.md`](https://github.com/FINNverse/FINN/blob/main/.claude/skills/finn-dataprep/SKILL.md);
not part of the installed package).

### 3d — ForClim and other gap models

ForClim does not need the BWI change/demographic tables at all — only an
**initial stand** and **site climate** to simulate forward from. Build the init
stand from the BWI3 slice (Step 3a) as a per-patch species + dbh list, attach
monthly climate and a soil bucket (Step 3b), and run the model's own projection.
Use the Step 2c rates only if you want to *validate* the model's output against
the observed BWI2→BWI3 change, rather than to drive it.

---

## Working with the latest BWI (2012→2022)

Steps 1–3 above use the ACCESS2003 release (BWI2→BWI3, 2002→2012). The Thünen
**inventory-archive API** adds **BWI4 (2022)**, so you can build the same
products for the **2012→2022** transition — the drought decade.

1. **Download** the inventory-archive (see `data/README.md`): set the public
   anonymous API key, then
   ```bash
   bash scripts/fetch_bwi4.sh          # plot, regeneration, deadwood, lookups
   bash scripts/fetch_bwi4_tree.sh     # the large /tree table (resumable)
   ```
   → `data/raw/bwi4/*.csv` (intervals `bwi2012`, `ci2017`, `bwi2022`, `ci2027`).

2. **Build the tree panel** in the same format as Step 1:
   ```bash
   Rscript scripts/40_import_new_schema.R    # -> outputs/new_schema/panel.csv (+ corners.csv)
   ```
   `panel.csv` has the **same columns as `bwi_tree_panel()$trees`** (`time` 0/1,
   `species_name`, `dbh`, `n_ha`, `fate`, …), so the code of Steps 2–3 runs
   unchanged — just load it in place of the tree table:
   ```r
   trees <- data.table::fread("outputs/new_schema/panel.csv")
   # the Step 2 plot aggregates and Step 3 model inputs then run as-is
   ```
   **But `recruit` does not mean the same thing.** The archive carries no
   back-modelled diameters, so `scripts/40` cannot split lateral ingrowth from
   genuine recruits: every tree new to the 2022 angle-count sample gets
   `fate = "recruit"` (and `is_recruit = TRUE`). Measured on the panel, the
   median 2022 "recruit" is 32 cm and 90 % are above 12 cm, so most are lateral
   ingrowth. Do not read `is_recruit` as recruitment in this release; take
   recruitment from the regeneration circles (as `scripts/70` and `scripts/71`
   do) or apply a diameter ceiling. The `scripts/40` panel also still has two
   fate defects (freshly dead trees, losses at the CI2017 corners);
   `scripts/91_import_new_schema_chained.R` builds the same panel with both
   fixed (`outputs/new_schema_chained/panel.csv`, same columns; see the README,
   *Status relative to the FINN-bwi study*). It runs after `scripts/40`, whose
   panel it uses for a regression check.
   *Prerequisite:* `scripts/40` reuses `outputs/finn/species_dt.csv` for
   consistent species lumping across releases, so run the ACCESS2003 pipeline
   (Step 0 + `scripts/04`) once first, or point it at your own species list.

**What differs in this release:**

- **Growth is measured directly** — survivors are matched across surveys by
  `(cluster, plot, tree_number)`, so growth is `dbh_2022 − dbh_2012`, with none
  of the Sloboda back-modelling that affected the ACCESS2003 growth (notes §5b).
- **The drought decade is inside the window** — the 2018–2020 dieback is
  captured, the period over which German forests turned from a net carbon sink
  to a source in the national greenhouse-gas inventory.
- **Salvage caveat** — drought-killed trees on managed corners were often
  salvage-logged (`tree_status` 12 = removed), which censors them out of the
  *natural*-mortality signal. For a clean biological mortality rate, restrict to
  **low-harvest corners** (each corner's harvest share is in `corners.csv`),
  where dead trees stay recorded as standing dead (`tree_status` 4/5). The
  BWI4 plot table also records *why* trees were removed
  (`plot.harvest_reason`: 0 regular, 1 storm, 2 insects/drought, 3 fire,
  4 calamity of uncertain cause), which `scripts/70` and `scripts/71` use to
  book salvage as natural death.

### Both decades in one model (FINN, `scripts/71`)

With one mortality observation per corner per decade, a single decade cannot
separate "this species dies more" from "this decade was dry". The same corners
observed over **2002→2012 and 2012→2022** can, which is what makes a per-year
drought covariate worth an annual timestep.
`scripts/71_prepare_two_interval_v8.R` builds the current input set,
`two-interval-v8`, which the FINN-bwi study fits:

```bash
bash scripts/fetch_dwd_drought.sh               # DWD 1 km monthly grids 2000-2022 (see data/README.md)
Rscript scripts/16_drought_index.R              # per-site, per-year soil-moisture anomaly
Rscript scripts/70_prepare_two_interval.R       # -> two-interval-v3/ (71 reuses its split_blocks.csv)
Rscript scripts/91_import_new_schema_chained.R  # -> outputs/new_schema_chained/ (decade-2 panel, fates fixed)
Rscript scripts/71_prepare_two_interval_v8.R    # -> outputs/prepared/two-interval-v8/
( cd outputs/prepared/two-interval-v8 && shasum -a 256 -c ../two-interval-v8.sha256 )
```

It needs the ACCESS2003 pipeline (`scripts/00`–`04`) and `scripts/40` to have
run, and FINN 0.2.0 or later (`remotes::install_github("FINNverse/FINN")` until
the CRAN release), which provides `makeObsData()` and the per-tree mortality
loss that `tree_obs_dt.csv` feeds. (Earlier versions of this guide named the
`loss-aggregation` development branch; that work is part of 0.2.0.) FINN 0.2.0
numbers sites, patches and species in radix order; older versions used the
locale's collation, so ids written under a non-`C` locale by an older FINN may
differ. The last line checks the result against the SHA-256 of the study's
files (`outputs/prepared/two-interval-v8.sha256`, in git). Expect every FINN
table and QC file to be `OK`, `split_blocks.csv` to be `FAILED` and
`v5_v8_state.csv` / `v5_v8_identity.csv` to be missing: `scripts/71` lists
v8's own tracts in `split_blocks.csv` (the study's file is that of its
`two-interval-v5`, with 55 more tracts; every v8 tract has the same block and
set in both) and drops the study's comparisons against its intermediate
versions. The script header lists every difference from the study's script.

`two-interval-v8` differs from `two-interval-v3` (`scripts/70`, kept unchanged)
in the three defects that the README lists under *Status relative to the
FINN-bwi study* — freshly dead trees counted as survivors, the 2012–2017 losses
at the CI2017 corners missing, lateral ingrowth counted twice in the start
stand (decade 1) and in the end state (both decades) — and in a third drought
column, `drought2` = mean of `drought` and `drought_lag1`. Removing the lateral
ingrowth from the 2002 start stand leaves 617 decade-1 site-periods without
any start tree, so they drop out (v8: 40,787 site-periods, 20,092 + 20,695).

**The regeneration target is a documented choice.** `obs_dt$reg` is the stock
of 5.0–6.9 cm stems on the r = 1.75 m regeneration circle, turned into a flux
of recruits per hectare and period at an assumed diameter increment of
0.30 cm/yr (`PARAMS$reg_incr_cm_yr` in `scripts/71`). The increment sets the
level, not the composition, and the level is about twice the decade-1 recruit
flux implied by the angle-count sample. Alternative targets are being
evaluated and are not part of this repository; until that is settled, v8
ships this one.

What `scripts/71` does (as `scripts/70`, unless marked v8):

- **One code path for both decades** — the same expansion, sample sizes,
  total-removal mortality and annualised rates as `scripts/04` (in `scripts/70`,
  interval 1 reproduces `scripts/04`'s output exactly).
- **Start and end stand from their own tallies** (v8) — the initial stand is
  the start tally only; the end state (`ba`, `trees`, `dbh` at year 10) is the
  living end tally expanded at its own end-survey `n_ha`. Growth and mortality
  come from the tallied trees, as before.
- **Freshly dead trees and CI2017 losses** (v8) — see the README.
- **Each decade is its own site-period** — the same corner appears twice with
  a different start stand and different environment years; interval-2 site ids
  are offset, and `siteID_dt` keeps `interval` and the physical `corner`.
- **Pooled scaling** — covariates are z-scaled once over both decades (so the
  same climate gets the same z-score in either decade); `env_scales.csv` holds
  the mean/sd and `env_unscaled_dt.csv` the natural units. Fit with FINN's
  automatic rescaling switched off (`env_autoscale = FALSE`).
- **Realised survey dates** for BWI4 (from `plot.acquisition_date`), so the
  drought years line up with the field work.
- **Two drought covariates** in `env_dt` — `drought`, the growing-season
  (Apr–Sep) soil-moisture anomaly of the model year, and `drought_lag1`, the
  previous calendar year's (v3), because bark-beetle mortality of spruce
  follows a drought by one to two years. Lower values are drier. v8 adds
  `drought2`, their mean, as a third column.
- **The tallied trees** in `tree_obs_dt.csv` (v2): one row per tree of the
  start tally with its removal 0/1 over the decade, survivor growth and the
  stems it stands for on the 0.1 ha patch (v8: no shading-only rows of
  laterally ingrown trees), for FINN's per-tree mortality loss
  (`fit(loss_aggregation = c(mortality = "tree"), tree_data = ...)`). The
  annualised site × species rate in `obs_dt$mort` cannot represent a complete
  removal, which biased the mean annual target well above the real ~3 %/yr.
- **Calamity salvage → natural death** in 2012→2022; the management tier
  counts regular harvest only, so the drought pulse stays in mortality.
- **Regeneration** for both decades from the fixed-area regeneration circles
  (dbh 5.0–6.9 cm, r = 1.75 m), as recruits per hectare per period (the
  documented choice above). The archive's 2012 counts are identical to
  ACCESS2003 `b3_jung`.
- **Spatially blocked sets** — `train/`, `holdout_blocked/`, `train_quarter/` (~25 % of train) and a small
  `pilot/`, next to `national/`, all with the same files. Whole tracts with
  both decades move together, in 100 km blocks; `PREPARATION.md` reports the
  distance from holdout to the nearest training tract.

---

## Reference

### Raw tables → products

| product | main raw inputs | helper / script |
|---|---|---|
| tree-level table | `b23_baeume_a_s`, `b23_baeume_m_s`, `b3_baeume`, `x_ba`, `x_pk` | `bwi_tree_panel()` (`R/bwi_trees.R`) |
| plot aggregates | the tree table + `n_ha` | `bhd2ba`, `wzp_*` (`R/wzp.R`) |
| corner population / filters | `b0_ecke`, `b23_ecke_w`, `b3_bestock` | `scripts/01`, `scripts/02` |
| environment | WorldClim, `b0_tab` (elevation), INSPIRE coords | `scripts/03`, `scripts/16` |
| FINN inputs | all of the above | `scripts/04` |
| FINN inputs, both decades | the above + inventory archive + DWD grids | `scripts/16`, `scripts/40`, `scripts/91`, `scripts/71` (`scripts/70`: v3) |

### Constants (`R/config.R`)

`WZP_BAF = 4` m²/ha · `WZP_MIN_BHD = 7` cm · `BWI_PERIOD_YEARS = 10` ·
`PATCH_SIZE_HA = 0.1` ha · `bhd2ba(dbh) = π·(dbh/200)²`.

### Helper index

- `R/bwi_io.R` — `bwi()` (read cached table), `bwi_species_table()`, `mdb_*`.
- `R/bwi_trees.R` — `bwi_tree_panel()` (the tree table), `bwi_periods()`.
- `R/wzp.R` — `wzp_representation`, `wzp_allocate`, `wzp_expand`,
  `wzp_effective_n` (angle-count expansion and rate counts).
- `R/config.R` — constants, `bhd2ba`, paths.
- `R/regions.R` — map corners to regions / Bundesländer.

### See also

- **`docs/bwi_data_notes.md`** — the verified rationale: sampling design, the
  ingrowth trap, fate codes, the two selection-bias corrections, time handling,
  coordinates, species, sample sizes.
- FINN's `C-Data_preparation` vignette and the `finn-dataprep` notes in the
  FINN repository (`.claude/skills/finn-dataprep/SKILL.md`) — the
  `obs_dt`/`env_dt`/`init` contract.
</content>
