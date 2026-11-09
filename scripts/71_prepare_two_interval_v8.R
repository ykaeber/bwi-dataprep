#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## bwi-dataprep | scripts/71_prepare_two_interval_v8.R
##
## Two-interval FINN inputs, VERSION 8 (two-interval-v8): the set the FINN-bwi
## study fits. Same design as two-interval-v3 (scripts/70, which stays as it
## is): BWI2 -> BWI3 (2002 -> 2012) and BWI3 -> BWI4 (2012 -> 2022) on the same
## corners, annual timestep, per-year drought, spatially blocked sets. v8
## removes three defects of v3 and adds one drought column:
##
##  1. FRESHLY DEAD TREES. The BWI records a tree dead for less than about 12
##     months like a live sample tree, with a "freshly dead" flag. Decade 1:
##     ACCESS2003 b2_wzp / b3_wzp `tot` (fixed below); decade 2: the archive's
##     `damage_dead` (fixed in scripts/91). Dead at the start survey -> not in
##     the start stand; a survivor flagged at the end survey -> natural death; a
##     new tree that is already dead -> not a recruit.
##  2. CI2017 LOSSES (decade 2, scripts/91). Where a corner was remeasured in
##     the 2017 carbon inventory, BWI2022 refers to CI2017, so trees lost in
##     2012-2017 have no 2022 record. scripts/91 chains the fate through CI2017.
##  3. LATERAL INGROWTH COUNTED TWICE IN THE STAND STATE.
##     a. INITIAL STAND = THE START TALLY ONLY. In decade 1 the ACCESS2003 panel
##        (R/bwi_trees.R) gives every tree that grew into the 2012 angle-count
##        sample laterally (`ingrowth_lateral`: new sample tree with a
##        back-modelled 2002 dbh >= 7 cm) a back-modelled 2002 state with
##        N_ha > 0. v3 expanded these trees INTO the 2002 stand on top of the
##        2002 tally that already represents them (init 1.26x BWI2's own 2002
##        tally, b2_baeume). v8 drops them from the start allocation
##        (init_cohort_dt) and from tree_obs_dt (v3's shading rows, died = NA).
##        Decade 2's start stand is the plain 2012 tally and is unchanged.
##     b. END STATE = FRESH EXPANSION OF THE LIVING END TALLY. v3 built ba,
##        trees and dbh at year 10 as "start-tally expansion carried forward +
##        every new tree at its end expansion"; a laterally ingrown tree was
##        already stand at t0, so it was counted twice. v8 expands every living
##        tree of the end tally at its own end-survey n_ha (4 / g(t1); decade 1:
##        Thuenen's N_ha_t1 incl. the edge correction), in both decades,
##        largest-remainder allocation per patch x species as at the start. The
##        t0 and t1 virtual trees are therefore no longer chained one-to-one;
##        growth and mortality come from the tallied trees (steps 6 / 6b), not
##        from makeObsData()'s chained virtual trees.
##     c. CONSEQUENCE, NOT A CHOICE: decade-1 site-periods whose start stand
##        expands to no tree once the lateral ingrowth is gone (2002 angle count
##        without a tree of the 15 species, or tallied trees rounding to 0
##        stems) drop out, because makeObsData() needs both censuses (617 of
##        20,709; counted in PREPARATION.md).
##  4. drought2 = (drought + drought_lag1) / 2 on the natural scale, z-scaled
##     with its pooled national mean / sd, appended as the last column of
##     env_dt / env_unscaled_dt and the last row of env_scales.
##
## Two details keep growth and mortality the same data as in the study's
## two-interval-v5 (= fixes 1 and 2 without fix 3):
##  * The tallied-tree sample of step 6 (rows with n_trees > 0) is selected with
##    the start allocation that still contains the lateral ingrowth (variables
##    `alloc_v5`, `n_trees_v5`). Dropping the lateral trees changes the
##    largest-remainder rounding of the other trees of the same patch x
##    species, and the n_trees > 0 rule would then move a few trees whose
##    expansion rounds near 0 stems in or out of the sample. That allocation is
##    used for nothing else.
##  * The end-survey rows of dead and harvested trees keep their start
##    expansion; they are not part of the end stand, and only makeObsData()'s
##    own mort (overwritten in step 6) reads them.
## The corner filters and the dbh screen run on the full panel, lateral
## ingrowth included.
##
## REGENERATION TARGET -- A DOCUMENTED CHOICE, NOT A SETTLED ONE. obs_dt$reg is
## the stock of stems in the 5.0-6.9 cm band on the r = 1.75 m circle (archive
## size classes 5 + 6), turned into a flux of recruits per ha and period at an
## assumed diameter increment of 0.30 cm/yr through the 1.9 cm band
## (PARAMS reg_*). The increment sets the level only, not the composition. This
## is v8's default target; it is about twice the decade-1 recruit flux implied
## by the angle-count sample. Alternative targets (a lower increment, recruits
## from the angle-count sample, a hybrid of both) are being evaluated and are
## not part of this repository; until that is decided, v8 ships this one.
##
## SPLIT. v3's split (scripts/70, two-interval-v3/split_blocks.csv) is reused:
## every tract keeps its v3 set, block and train_quarter / pilot membership. A
## tract new in v8 (it has v8 site-periods but none in v3, since the fixed fates
## move a few corners across the corner filters; 6 tracts) takes the set of its
## v3 block and joins neither train_quarter nor pilot; a block unknown to v3
## stops the script.
##
## RELATION TO THE FINN-bwi STUDY. The study writes v8 with its
## scripts/131_prepare_two_interval_v8.R, a copy of its scripts/92 (v5) with
## fixes 3 and 4 added, reading v5's split_blocks.csv and checking against the
## v5 and v7 folders. This script is that script with three changes, none of
## which touches a FINN table:
##  a. The split rule of the study's scripts/92 is applied here to v8's tracts
##     (block from the tract's EPSG:3035 coordinates; set from v3) instead of
##     copying v5's split_blocks.csv. Every v8 tract gets the block, set,
##     train_quarter and pilot it has in the study's file; the file lists v8's
##     12,388 tracts instead of v5's 12,443 (55 tracts whose only site-periods
##     are the dropped decade-1 ones), so it is not byte-identical to the
##     study's.
##  b. The regression checks against the study's v5 and v7 folders are dropped
##     (v5_v8_state.csv, v5_v8_identity.csv, the siteID / env / env_scales
##     comparisons and the v5-to-v8 siteID map). Every check on the stitched
##     tables themselves is kept. The dropped decade-1 site-periods are counted
##     from this script's own start-empty sites instead of from v5.
##  c. PREPARATION.md describes v8 relative to v3.
## "v5" in comments and variable names refers to the study's two-interval-v5.
## outputs/prepared/two-interval-v8.sha256 holds the SHA-256 of the study's
## files; the 43 FINN tables and QC files of a run of this script match it (see
## README.md). Check with
##   (cd outputs/prepared/two-interval-v8 && shasum -a 256 -c ../two-interval-v8.sha256)
## which reports split_blocks.csv as FAILED (change a) and v5_v8_state.csv /
## v5_v8_identity.csv as missing (change b).
##
## Inputs (run scripts/00-04, 16, 40, 70 and 91 first)
##   data/bwi_csv/*                        ACCESS2003 tables (bwi_tree_panel,
##                                         species_lookup; b2_wzp/b3_wzp `tot` =
##                                         freshly dead; b2_baeume for the check)
##   data/raw/bwi4/{plot,regeneration}.csv inventory archive (scripts/fetch_bwi4.sh)
##   outputs/new_schema_chained/panel.csv  BWI3->BWI4 tree panel, fates fixed (scripts/91)
##   outputs/calibration_set_corners.csv   interval-1 corners, all tiers (scripts/02)
##   outputs/calibration_set_species.csv   the 15 species (scripts/02)
##   outputs/corner_env.csv                climate, elevation, x/y (scripts/03)
##   outputs/site_drought_annual.csv       per-year drought anomaly (scripts/16)
##   outputs/prepared/two-interval-v3/split_blocks.csv   the split (scripts/70)
##
## Outputs  outputs/prepared/two-interval-v8/  (the layout of v3, plus drought2)
##   PREPARATION.md          what changed vs v3, decisions, checks, PARAMS
##   national/obs_dt.csv     observations, both decades (obs at year 10)
##   national/env_dt.csv     covariates per siteID x year 0..10, z-scaled once
##                           over both decades (management kept on 0..3), incl. drought2
##   national/env_unscaled_dt.csv, env_scales.csv
##   national/init_cohort_dt.csv   start-of-decade stands (start tally only)
##   national/species_dt.csv       the 15 species, one code map for both decades
##   national/siteID_dt.csv        siteID <-> interval, physical corner (crosswalk)
##   national/tree_obs_dt.csv      tallied trees of the start tally: siteID, year,
##                                 patchID, species, dbh (start), trees (stems on
##                                 the 0.1 ha patch), growth (/yr), died (0/1)
##   national/reg_target.csv       regeneration target per corner x species
##   national/wzp_qc_i{1,2}.csv    basal-area conservation of the START expansion
##   train/, holdout_blocked/, train_quarter/, pilot/   the same FINN tables for
##                           the spatially blocked split, siteIDs renumbered
##                           1..n, siteID_dt keeps siteID_national
##   split_blocks.csv        tract -> block -> set (v8's tracts)
##
## Helpers sourced (none of them changes what the data mean):
##   R/config.R    paths and constants (WZP_BAF, PATCH_SIZE_HA, bhd2ba, AV_WZP)
##   R/bwi_io.R    bwi(): read a cached ACCESS2003 table
##   R/wzp.R       angle-count expansion and tallied sample sizes
##   R/bwi_trees.R bwi_tree_panel(): the ACCESS2003 BWI2->BWI3 tree panel
## Packages: data.table, FINN (0.2.0), sf, FNN.
##
## Run from the repository root:  Rscript scripts/71_prepare_two_interval_v8.R
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

suppressPackageStartupMessages({
  library(data.table)
  library(FINN)
})
source("R/config.R")
source("R/bwi_io.R")
source("R/wzp.R")
source("R/bwi_trees.R")

#=#=# PARAMS #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
PARAMS <- list(
  version = "two-interval-v8",

  ## decade-2 panel with the fates fixed (scripts/91)
  panel2_dir = file.path(OUT_DIR, "new_schema_chained"),
  ## v3 (scripts/70): its split_blocks.csv is the split that v8 keeps
  v3_dir     = file.path(OUT_DIR, "prepared", "two-interval-v3"),

  ## Nominal patch area for the angle-count expansion (ha); must equal the
  ## patch_size passed to fit(). Same as scripts/04.
  patch_size_ha = PATCH_SIZE_HA,

  ## FINN reads one aggregation window for all sites, so each decade keeps only
  ## corners whose realised period is exactly this many years (scripts/04).
  period = 10L,

  ## scripts/02 corner filters, applied to interval 2 as they were to interval 1.
  min_ba_ha_t0  = 5,              # drop near-empty corners (m2/ha)
  species_X     = 0.99,           # the 15 species must hold this share of BA
  period_range  = c(8L, 12L),     # realised period window before the exact cut

  ## scripts/04 site screen: drop a site if a tree's diameter change over the
  ## period falls outside this window (cm) -- transcription errors.
  dbh_growth_thresh = c(-10, 60),

  ## Static covariates from corner_env.csv (x, y ride along for hybrids/ALE).
  env_vars    = c("mat_c", "map_mm", "elev_m", "x", "y"),
  ## Per-year anomaly from scripts/16 (standardised per site over 2000-2022, so
  ## both decades share one baseline); stored in env_dt as `drought`.
  drought_var = "z_soil_moist",
  ## the single drought index: the mean of the current and the previous
  ## calendar year's anomaly, pooled-scaled like the others
  drought2_parents = c("drought", "drought_lag1"),

  ## Management tiers exactly as scripts/01: M0 = no tallied harvest, then the
  ## basal-area share removed: <= 5 % M1, <= 20 % M2, above M3 -> 0..3.
  tier_cut = c(M1 = 0.05, M2 = 0.20),
  ## plot.harvest_reason codes booked as natural death in interval 2:
  ## 1 storm, 2 insects/drought, 3 other calamity (fire), 4 calamity, cause
  ## uncertain. 0 regular and 5 unknown stay harvest.
  calamity_reasons = 1:4,

  ## Regeneration target (stock -> flux; a documented choice, see the header):
  ## stems in the 5.0-6.9 cm band on the r = 1.75 m circle, crossing 7 cm at
  ## reg_incr_cm_yr.
  reg_classes    = c(5L, 6L),
  reg_circle_r_m = 1.75,
  reg_incr_cm_yr = 0.30,          # sets the MAGNITUDE only; composition is free of it
  reg_band_cm    = 1.9,

  ## Block size of the spatially blocked split (km on EPSG:3035), as scripts/70;
  ## used to place each tract in its v3 block
  block_km = 100,
  ## the four subsets written beside national/ (membership from v3's split_blocks.csv)
  sets = c("train", "holdout_blocked", "train_quarter", "pilot"),

  bwi4_dir = file.path(BWI_RAW_DIR, "bwi4"),
  seed     = 42L                  # as scripts/70 (nothing below draws at random; kept for parity)
)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

set.seed(PARAMS$seed)
P        <- PARAMS$period
OUT_VER  <- file.path(OUT_DIR, "prepared", PARAMS$version)
OUT_NAT  <- file.path(OUT_VER, "national")
stopifnot(file.exists(file.path(PARAMS$v3_dir, "split_blocks.csv")),   # run scripts/70 first
          file.exists(file.path(PARAMS$panel2_dir, "panel.csv")),      # run scripts/91 first
          !dir.exists(OUT_VER))                       # a version is never modified
dir.create(OUT_NAT, recursive = TRUE, showWarnings = FALSE)
add_site_completeness_flags <- FINN:::add_site_completeness_flags
tier_intensity <- c(M0 = 0, M1 = 1, M2 = 2, M3 = 3)

#' Read an inventory-archive CSV, dropping the few lines that the paged download
#' glued together at page boundaries (they carry twice the field count; a plain
#' fread() stops at the first one and silently returns a fraction of the file).
read_archive <- function(name, select = NULL) {
  f  <- file.path(PARAMS$bwi4_dir, paste0(name, ".csv"))
  nf <- count.fields(f, sep = ",", quote = "\"", comment.char = "",
                     blank.lines.skip = FALSE)
  ln <- readLines(f)
  stopifnot(length(nf) == length(ln))
  keep <- !is.na(nf) & nf == nf[1]
  message(sprintf("[archive] %s: %d lines, %d malformed dropped",
                  name, length(ln) - 1L, sum(!keep)))
  fread(text = ln[keep], select = select, showProgress = FALSE)
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## prepare_interval(): scripts/04 steps 1-7b for ONE decade, unscaled
##
## panel   tallied-tree panel in bwi_tree_panel() format (time 0/1, fate, n_ha)
## corners Tnr, Enr, siteName, patchName, period_length, year_t0, mgmt_tier
## Returns obs, env (unscaled), init, siteID_dt -- with siteIDs 1..n -- QC and
## the design estimates of the start and end stands per site (for the checks).
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
prepare_interval <- function(panel, corners, label, V, env_in, drought) {
  message(sprintf("\n=== %s ===", label))

  ## ---- 1) corners x environment; management as an env predictor (0..3) ----
  corners <- copy(corners)[, .(Tnr, Enr, siteName, patchName, period_length,
                               year_t0,
                               management = as.numeric(tier_intensity[mgmt_tier]))]
  corners <- merge(corners, env_in[, c("Tnr", "Enr", PARAMS$env_vars), with = FALSE],
                   by = c("Tnr", "Enr"))
  env_cols <- c(PARAMS$env_vars, "management")
  message(sprintf("[%s] corners with environment: %d", label, nrow(corners)))

  ## ---- 2) tallied-tree panel restricted to the cohort and to V -----------
  panel <- copy(panel)
  drop <- intersect(c("siteName", "patchName", "period_length"), names(panel))
  if (length(drop)) panel[, (drop) := NULL]
  panel <- merge(panel, corners[, .(Tnr, Enr, siteName, patchName, period_length)],
                 by = c("Tnr", "Enr"))
  n0 <- nrow(panel)
  ## share of start-of-decade basal area that v3 carried as back-modelled
  ## lateral ingrowth (interval 1 only) -- what v8 removes from the start stand
  t0 <- panel[time == 0L]
  lat_share <- if ("ingrowth_lateral" %in% t0$fate)
    t0[, sum((fate == "ingrowth_lateral") * n_ha * bhd2ba(dbh)) / sum(n_ha * bhd2ba(dbh))]
  else 0
  panel <- panel[species_name %in% V]
  message(sprintf("[%s] restricted to V: %d -> %d records", label, n0, nrow(panel)))

  ## the dbh screen runs on the full panel (lateral ingrowth included), as in
  ## v3 / v5
  wide <- dcast(panel[, .(treeKey, time, dbh)], treeKey ~ time, value.var = "dbh")
  setnames(wide, c("treeKey", "dbh_t0", "dbh_t1"), skip_absent = TRUE)
  if (all(c("dbh_t0", "dbh_t1") %in% names(wide))) {
    bad <- wide[is.finite(dbh_t0) & is.finite(dbh_t1) &
                (dbh_t1 - dbh_t0 < PARAMS$dbh_growth_thresh[1] |
                 dbh_t1 - dbh_t0 > PARAMS$dbh_growth_thresh[2]), treeKey]
    if (length(bad)) {
      bad_sites <- unique(panel[treeKey %in% bad, siteName])
      message(sprintf("[%s] dbh screen: %d implausible trees -> dropping %d sites",
                      label, length(bad), length(bad_sites)))
      panel <- panel[!siteName %in% bad_sites]
    }
  }

  ## ---- 3) angle-count expansion -> virtual trees on a fixed-area patch ----
  ## v8: the start and the end stand are expanded SEPARATELY, each from its own
  ## tally (fix 3).
  panel[, patchKey := paste(siteName, patchName, sep = "/")]
  ## (3a) start stand = the start tally only: no back-modelled lateral ingrowth
  alloc <- panel[time == 0L & fate != "ingrowth_lateral",
                 .(treeKey, patchKey, species_name, n_ha, dbh)]
  invisible(wzp_allocate(alloc, by = c("patchKey", "species_name"),
                         patch_size_ha = PARAMS$patch_size_ha))
  qc <- wzp_ba_conservation(alloc, by = "patchKey",
                            patch_size_ha = PARAMS$patch_size_ha)
  message(sprintf("[%s] wzp basal-area conservation, start: median |err| %.2f %%",
                  label, 100 * median(abs(qc$rel_err), na.rm = TRUE)))
  ## (3b) end stand = every LIVING tree of the end tally at its end-survey
  ## n_ha (survivors, new trees, decade-1 lateral ingrowth alike)
  alloc_end <- panel[time == 1L & living %in% TRUE,
                     .(treeKey, patchKey, species_name, n_ha, dbh)]
  invisible(wzp_allocate(alloc_end, by = c("patchKey", "species_name"),
                         patch_size_ha = PARAMS$patch_size_ha))
  qc_end <- wzp_ba_conservation(alloc_end, by = "patchKey",
                                patch_size_ha = PARAMS$patch_size_ha)
  message(sprintf("[%s] wzp basal-area conservation, end: median |err| %.2f %%",
                  label, 100 * median(abs(qc_end$rel_err), na.rm = TRUE)))
  ## the start allocation WITH the lateral ingrowth (the study's v5), used ONLY
  ## to select the tallied-tree sample of step 6 (see the header)
  alloc_v5 <- panel[time == 0L, .(treeKey, patchKey, species_name, n_ha, dbh)]
  invisible(wzp_allocate(alloc_v5, by = c("patchKey", "species_name"),
                         patch_size_ha = PARAMS$patch_size_ha))
  panel <- merge(panel, alloc[, .(treeKey, n_start = n_trees)], by = "treeKey", all.x = TRUE)
  panel <- merge(panel, alloc_end[, .(treeKey, n_end = n_trees)], by = "treeKey", all.x = TRUE)
  panel <- merge(panel, alloc_v5[, .(treeKey, n_trees_v5 = n_trees)], by = "treeKey", all.x = TRUE)
  for (v in c("n_start", "n_end", "n_trees_v5")) panel[is.na(get(v)), (v) := 0L]
  ## living end rows carry the end expansion; start rows, and the end rows of
  ## dead / harvested trees (not part of the end stand), the start expansion
  panel[, n_trees := fifelse(time == 1L & living %in% TRUE, n_end, n_start)]
  ## how many at-risk tallied trees the sample rule keeps that the v8
  ## allocation alone would have moved (in or out): rounding near 0 stems
  n_sample_flip <- panel[time == 0L & fate %in% c("survivor", "dead", "harvest"),
                         sum((n_start > 0L) != (n_trees_v5 > 0L))]
  virtual <- wzp_expand(panel,
                        keep = c("treeKey", "siteName", "patchName",
                                 "species_name", "dbh", "time", "status",
                                 "living", "mort", "fate", "period_length"),
                        index_by = c("treeKey", "time"))
  ## sites whose start stand expands to no tree at all (decade 1: the 2002
  ## tally holds no tree of the 15 species, only back-modelled lateral
  ## ingrowth, or its trees round to 0 stems on the patch). makeObsData() needs
  ## both censuses, so these site-periods drop out (see PREPARATION.md).
  start_empty <- panel[!siteName %chin% unique(virtual[time == 0L, siteName]), .(
    n_start_tally  = uniqueN(treeKey[time == 0L & fate != "ingrowth_lateral"]),
    n_lateral_t0   = uniqueN(treeKey[time == 0L & fate == "ingrowth_lateral"]),
    ba_start_tally = sum((time == 0L & fate != "ingrowth_lateral") * n_ha * bhd2ba(dbh)),
    ba_lateral_t0  = sum((time == 0L & fate == "ingrowth_lateral") * n_ha * bhd2ba(dbh)),
    ba_end_tally   = sum((time == 1L & living %in% TRUE) * n_ha * bhd2ba(dbh))), by = siteName]

  ## ---- 4) tree_dt in FINN's schema, then makeObsData() --------------------
  ## (makeObsData()'s state at each census counts the living virtual trees, so
  ## ba / trees / dbh at year 10 are the fresh end expansion; its chained growth
  ## and mortality are replaced by the tallied-tree values in step 6)
  tree_dt <- virtual[, .(siteName, patchName, treeName, year = time,
                         species_name, dbh, status, living)]
  tree_dt <- merge(tree_dt, unique(virtual[, .(siteName, period_length)]),
                   by = "siteName")
  tree_dt[, year := fifelse(year == 0L, 0L, as.integer(period_length))]
  tree_dt[, period_length := NULL]
  env_dt0 <- unique(corners[, c("siteName", env_cols), with = FALSE])
  env_dt0 <- merge(unique(tree_dt[, .(siteName, year)]), env_dt0, by = "siteName")
  add_site_completeness_flags(env_dt = env_dt0, tree_dt = tree_dt)
  obs_list <- makeObsData(
    tree_dt           = tree_dt[complete == TRUE],
    plotsize          = PARAMS$patch_size_ha,
    aggregate_by_site = FALSE,
    minNyears         = 2L,
    fix_period_length = NULL,
    dbh_growth_thresh = PARAMS$dbh_growth_thresh,
    Npatches          = NULL,
    Nspecies          = NULL,
    NspeciesQuantile  = NULL
  )
  obs_raw  <- obs_list$obs_dt
  tree_obs <- obs_list$tree_dt[siteName %in% unique(obs_raw$siteName)]

  ## ---- 5) integer site / patch / species indices --------------------------
  inputs <- resolveSiteIDs(tree_dt = as.data.table(tree_obs),
                           obs_dt = as.data.table(obs_raw),
                           env_dt = as.data.table(env_dt0),
                           createInitCohorts = FALSE)
  obs_dt     <- inputs$obs_dt
  env_dt     <- inputs$env_dt
  siteID_dt  <- inputs$siteID_dt
  init_trees <- inputs$tree_dt[year == 0]
  init_trees[living == FALSE | is.na(living), trees := 0]
  if (!"trees" %in% names(init_trees)) init_trees[, trees := 1L]

  ## ---- 6) tallied sample sizes; total-removal mortality; survivor growth --
  dbh_w <- dcast(panel[, .(treeKey, time, dbh)], treeKey ~ time, value.var = "dbh")
  setnames(dbh_w, c("0", "1"), c("dbh_t0", "dbh_t1"), skip_absent = TRUE)
  ## the sample is selected with the start allocation that includes the
  ## lateral ingrowth (n_trees_v5), see the header
  tal <- merge(panel[time == 0L & n_trees_v5 > 0L,
                     .(treeKey, siteName, patchName, species_name, n_ha, fate)],
               dbh_w, by = "treeKey", all.x = TRUE)
  tal[, `:=`(
    at_risk    = fate %in% c("survivor", "dead", "harvest"),
    died       = fate == "dead",
    harvested  = fate == "harvest",
    growth_ok  = fate == "survivor",
    rel_growth = fifelse(is.finite(dbh_t0) & is.finite(dbh_t1) & dbh_t0 > 0,
                         dbh_t1 / dbh_t0 - 1, NA_real_)
  )]
  eff <- wzp_effective_n(tal, by = c("siteName", "patchName", "species_name"))
  eff_site <- eff[, .(
    n_at_risk  = sum(n_at_risk, na.rm = TRUE),
    n_died     = sum(n_died_eff, na.rm = TRUE),
    n_harvest  = sum(n_harvest_eff, na.rm = TRUE),
    growth_n   = sum(growth_n, na.rm = TRUE),
    growth_obs = if (sum(growth_n, na.rm = TRUE) > 0)
                   sum(growth_w * growth_n, na.rm = TRUE) /
                   sum(growth_n[is.finite(growth_w)], na.rm = TRUE) else NA_real_
  ), by = .(siteName, species_name)]
  obs_dt <- merge(obs_dt, unique(siteID_dt[, .(siteID, siteName)]),
                  by = "siteID", all.x = TRUE)
  obs_dt <- merge(obs_dt, eff_site, by = c("siteName", "species_name"),
                  all.x = TRUE, suffixes = c("_expanded", ""))
  for (v in c("n_at_risk", "n_died", "n_harvest", "growth_n"))
    obs_dt[is.na(get(v)), (v) := 0L]
  obs_dt[year == 0L, `:=`(n_at_risk = 0L, n_died = 0L, n_harvest = 0L,
                          growth_n = 0L, growth_obs = NA_real_)]
  obs_dt[, mort := fifelse(n_at_risk > 0, (n_died + n_harvest) / n_at_risk,
                           NA_real_)]
  obs_dt[, growth := fifelse(growth_n > 0 & is.finite(growth_obs),
                             growth_obs, NA_real_)]
  obs_dt[, growth_obs := NULL]

  ## ---- 6b) the tallied trees, for FINN's per-tree mortality loss ----------
  ## One row per tree of the start tally: removal 0/1 over the decade (natural
  ## death or harvest, the same total removal as `mort`), annual relative growth
  ## of survivors, and `trees` = the stems on the 0.1 ha patch that the tree
  ## stands for, so FINN computes light from the expanded stand. v8: the start
  ## tally only -- no laterally ingrown shading rows (v3's died = NA rows).
  tree_tab <- merge(panel[time == 0L & fate %in% c("survivor", "dead", "harvest"),
                          .(treeKey, siteName, species_name, n_ha, fate)],
                    dbh_w, by = "treeKey", all.x = TRUE)
  tree_tab <- tree_tab[is.finite(dbh_t0) & dbh_t0 > 0]
  tree_tab[, `:=`(trees  = n_ha * PARAMS$patch_size_ha,
                  died   = fcase(fate %in% c("dead", "harvest"), 1, fate == "survivor", 0, default = NA_real_),
                  growth = fifelse(fate == "survivor" & is.finite(dbh_t1), (dbh_t1 / dbh_t0)^(1 / P) - 1, NA_real_))]
  tree_tab <- tree_tab[, .(siteName, species_name, dbh = dbh_t0, trees, growth, died)]

  drop_cols <- grep("_expanded$", names(obs_dt), value = TRUE)
  if (length(drop_cols)) obs_dt[, (drop_cols) := NULL]

  ## ---- 7) annual time grid + per-year drought -----------------------------
  site_period <- unique(corners[, .(siteName, period_length)])
  obs_dt <- merge(obs_dt, site_period, by = "siteName", all.x = TRUE)
  n_before <- uniqueN(obs_dt$siteID)
  obs_dt <- obs_dt[period_length == P]
  message(sprintf("[%s] period exactly %d yr: %d -> %d sites",
                  label, P, n_before, uniqueN(obs_dt$siteID)))
  obs_dt <- obs_dt[year == 1L]
  obs_dt[, year := P]
  obs_dt[is.finite(growth) & growth > -1, growth := (1 + growth)^(1 / P) - 1]
  obs_dt[is.finite(mort), mort := 1 - (1 - pmin(pmax(mort, 0), 1))^(1 / P)]
  obs_dt[is.finite(mort) & n_at_risk > 0,
         n_died := as.integer(round(mort * n_at_risk))]

  env_dt <- env_dt[year == 0L]
  env_dt <- rbindlist(lapply(0:P, function(k) copy(env_dt)[, year := k]))
  key  <- unique(corners[, .(siteName, year_t0)])
  grid <- CJ(siteName = unique(key$siteName), year = 0:P, unique = TRUE)
  grid <- merge(grid, key, by = "siteName")
  grid[, cal_year := year_t0 + year]
  grid <- merge(grid, drought, by = c("siteName", "cal_year"), all.x = TRUE)
  ## the previous calendar year's anomaly (v3)
  grid <- merge(grid, drought[, .(siteName, cal_year = cal_year + 1L, drought_lag1 = drought)],
                by = c("siteName", "cal_year"), all.x = TRUE)
  drought_cover <- mean(is.finite(grid$drought))
  lag_cover <- mean(is.finite(grid$drought_lag1))
  grid[!is.finite(drought), drought := 0]
  grid[!is.finite(drought_lag1), drought_lag1 := 0]
  if (!"siteName" %in% names(env_dt))
    env_dt <- merge(env_dt, unique(siteID_dt[, .(siteID, siteName)]),
                    by = "siteID", all.x = TRUE)
  env_dt <- merge(env_dt, grid[, .(siteName, year, drought, drought_lag1)],
                  by = c("siteName", "year"), all.x = TRUE)
  env_dt[!is.finite(drought), drought := 0]
  env_dt[!is.finite(drought_lag1), drought_lag1 := 0]
  message(sprintf("[%s] drought attached: %.1f %% (lag %.1f %%) of site-years covered, driest model-year mean %.2f",
                  label, 100 * drought_cover, 100 * lag_cover,
                  min(env_dt[siteID %in% obs_dt$siteID, mean(drought), by = year]$V1)))

  ## ---- 7b) align all tables on the surviving sites, siteIDs 1..n ---------
  keep_sites <- sort(unique(obs_dt$siteID))
  env_dt     <- env_dt[siteID %in% keep_sites]
  init_trees <- init_trees[siteID %in% keep_sites]
  siteID_dt  <- siteID_dt[siteID %in% keep_sites]
  site_map   <- data.table(siteID = keep_sites, siteID_new = seq_along(keep_sites))
  for (dt in list(obs_dt, env_dt, init_trees, siteID_dt))
    dt[site_map, on = "siteID", siteID := i.siteID_new]
  stopifnot(setequal(unique(env_dt$siteID), unique(obs_dt$siteID)),
            setequal(unique(init_trees$siteID), unique(obs_dt$siteID)))

  ## design estimates of the start and end stands on the surviving sites (per
  ## ha, tallied trees at their n_ha), for the checks in PREPARATION.md
  design <- panel[siteName %chin% unique(siteID_dt$siteName), .(
    ba_start_design    = sum((time == 0L & fate != "ingrowth_lateral") * n_ha * bhd2ba(dbh)),
    stems_start_design = sum((time == 0L & fate != "ingrowth_lateral") * n_ha),
    ba_end_design      = sum((time == 1L & living %in% TRUE) * n_ha * bhd2ba(dbh)),
    stems_end_design   = sum((time == 1L & living %in% TRUE) * n_ha)), by = siteName]

  list(obs = obs_dt, env = env_dt[, c("siteID", "year", env_cols, "drought", "drought_lag1"), with = FALSE],
       init = init_trees, siteID_dt = siteID_dt, qc = qc, qc_end = qc_end, trees = tree_tab,
       design = design, start_empty = start_empty,
       diag = list(sites = length(keep_sites), lat_share = lat_share,
                   drought_cover = drought_cover, n_sample_flip = n_sample_flip))
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Shared inputs
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
V <- fread(file.path(OUT_DIR, "calibration_set_species.csv"))$species_name
species_dt <- data.table(species = seq_along(V), species_name = sort(V))
env_in  <- fread(file.path(OUT_DIR, "corner_env.csv"))
drought <- fread(file.path(OUT_DIR, "site_drought_annual.csv"))[
  , .(siteName, cal_year = year, drought = get(PARAMS$drought_var))]
splk <- bwi("species_lookup")[, .(Ba, species_name)]    # archive tree_species == Ba

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Interval 1: BWI2 -> BWI3 (ACCESS2003)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
c1 <- fread(file.path(OUT_DIR, "calibration_set_corners.csv"))
p1 <- bwi_tree_panel()$trees

## ---- fix 1, decade 1: freshly dead trees (ACCESS2003 b2_wzp/b3_wzp `tot`,
## "Frisch abgestorben 1=ja") --------------------------------------------------
fd_end   <- bwi("b3_wzp", select = c("Tnr", "Enr", "Bnr", "tot"))[tot == 1L, .(Tnr, Enr, Bnr, fd_end = TRUE)]
fd_start <- bwi("b2_wzp", select = c("Tnr", "Enr", "Bnr", "tot"))[tot == 1L, .(Tnr, Enr, Bnr, fd_start = TRUE)]
p1 <- merge(p1, fd_end, by = c("Tnr", "Enr", "Bnr"), all.x = TRUE)
p1 <- merge(p1, fd_start, by = c("Tnr", "Enr", "Bnr"), all.x = TRUE)
p1[is.na(fd_end), fd_end := FALSE]
p1[is.na(fd_start), fd_start := FALSE]
## dead at the start survey: not part of the 2002 stand, not at risk
p1 <- p1[fd_start == FALSE]
## a recruit that is dead at its first tally never was a living recruit
p1 <- p1[!(fd_end & fate == "recruit")]
## survivors flagged at 2012 died in the decade; flagged lateral ingrowth is
## dead at t1 (it stays outside the at-risk cohort, as all lateral ingrowth)
p1[fd_end & fate == "survivor", fate := "dead"]
p1[fd_end & time == 1L & fate %in% c("dead", "ingrowth_lateral") & status == "alive",
   `:=`(status = "dead", living = FALSE, mort = TRUE)]
p1[, c("fd_end", "fd_start") := NULL]

i1 <- prepare_interval(p1, c1, "interval 1 (2002->2012)", V, env_in, drought)
rm(p1); gc(verbose = FALSE)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Interval 2: BWI3 -> BWI4 (inventory archive) -- corners built here
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
p2 <- fread(file.path(PARAMS$panel2_dir, "panel.csv"))
pl <- read_archive("plot", c("id", "interval_name", "cluster_name", "plot_name",
                             "acquisition_date", "harvest_reason"))
pl[, `:=`(Tnr = as.integer(cluster_name), Enr = as.integer(plot_name),
          yr = as.integer(substr(acquisition_date, 1, 4)))]

## realised survey years: earliest field year per tract at each inventory
yrs <- pl[interval_name %in% c("bwi2012", "bwi2022") & is.finite(yr),
          .(yr = min(yr)), by = .(Tnr, interval_name)]
yrs <- dcast(yrs, Tnr ~ interval_name, value.var = "yr")
setnames(yrs, c("bwi2012", "bwi2022"), c("year_t0", "year_t1"))
yrs[is.na(year_t0), year_t0 := 2012L]
yrs[is.na(year_t1), year_t1 := 2022L]
yrs[, period_length := year_t1 - year_t0]

## calamity salvage on the BWI4 plot -> the removed trees died naturally
calam <- unique(pl[interval_name == "bwi2022" &
                   harvest_reason %in% PARAMS$calamity_reasons, .(Tnr, Enr)])
calam[, calamity := TRUE]
p2 <- merge(p2, calam, by = c("Tnr", "Enr"), all.x = TRUE)
p2[is.na(calamity), calamity := FALSE]
n_salvage <- p2[time == 0L & calamity & fate == "harvest", .N]
n_harvest_all <- p2[time == 0L & fate == "harvest", .N]
p2[calamity & fate == "harvest", fate := "dead"]
p2[time == 1L & fate == "dead" & status == "harvest", status := "dead"]
p2[, mort := status == "dead"]
p2[, calamity := NULL]
message(sprintf("[interval 2] calamity salvage booked as natural death: %d of %d harvested stems (%.1f %%)",
                n_salvage, n_harvest_all, 100 * n_salvage / n_harvest_all))

## corner metrics at the start of the decade (the 2012 tally)
cs <- p2[time == 0L, .(
  ba_ha_t0      = sum(n_ha * bhd2ba(dbh)),
  ba_V          = sum((species_name %in% V) * n_ha * bhd2ba(dbh)),
  ba_harvest    = sum((fate == "harvest") * n_ha * bhd2ba(dbh)),
  n_tal_harvest = sum(fate == "harvest"),
  n_untrace     = sum(fate == "untraceable")
), by = .(Tnr, Enr)]
cs[, `:=`(share_V       = fifelse(ba_ha_t0 > 0, ba_V / ba_ha_t0, 0),
          harvest_share = fifelse(ba_ha_t0 > 0, ba_harvest / ba_ha_t0, 0))]
cs[, mgmt_tier := fcase(n_tal_harvest == 0L,                     "M0",
                        harvest_share <= PARAMS$tier_cut[["M1"]], "M1",
                        harvest_share <= PARAMS$tier_cut[["M2"]], "M2",
                        default = "M3")]
c2 <- merge(cs, yrs, by = "Tnr", all.x = TRUE)
c2[, `:=`(siteName = paste(Tnr, Enr, sep = "_"), patchName = "1")]   # as scripts/02: one patch per corner

funnel <- data.table(step = "corners in the BWI3->BWI4 panel", n = nrow(c2))
step <- function(d, lab) { funnel <<- rbind(funnel, data.table(step = lab, n = nrow(d))); d }
c2 <- step(c2[n_untrace == 0L],                               "no untraceable trees")
c2 <- step(c2[ba_ha_t0 >= PARAMS$min_ba_ha_t0],               "basal area >= 5 m2/ha")
c2 <- step(c2[share_V >= PARAMS$species_X],                   "15 species >= 99 % of BA")
c2 <- step(c2[period_length >= PARAMS$period_range[1] &
              period_length <= PARAMS$period_range[2]],         "realised period 8-12 yr")
c2 <- step(c2[paste(Tnr, Enr) %in% env_in[, paste(Tnr, Enr)]], "environment available")
print(funnel)
message("[interval 2] tiers: ",
        paste(c2[, .N, by = mgmt_tier][order(mgmt_tier), sprintf("%s=%d", mgmt_tier, N)],
              collapse = ", "))

i2 <- prepare_interval(p2, c2, "interval 2 (2012->2022)", V, env_in, drought)
rm(p2); gc(verbose = FALSE)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Regeneration target, both decades, from the archive circles (the documented
## choice of the header: circle stock -> flux at PARAMS$reg_incr_cm_yr)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
rg <- read_archive("regeneration", c("plot_id", "tree_species",
                                     "tree_size_class", "tree_count"))
rg <- merge(rg, pl[, .(plot_id = id, interval_name, Tnr, Enr)], by = "plot_id")
rg <- rg[interval_name %in% c("bwi2012", "bwi2022") &
         tree_size_class %in% PARAMS$reg_classes]
rg <- merge(rg, splk, by.x = "tree_species", by.y = "Ba", all.x = TRUE)
off_V <- rg[, .(share_off_V = sum(tree_count * !(species_name %in% V)) / sum(tree_count)),
            by = interval_name]
area_ha     <- pi * PARAMS$reg_circle_r_m^2 / 1e4
flux_factor <- PARAMS$reg_incr_cm_yr / PARAMS$reg_band_cm * P   # band stock -> recruits per period
reg_tab <- rg[species_name %in% V,
              .(n_circle = sum(tree_count)), by = .(interval_name, Tnr, Enr, species_name)]
reg_tab[, reg := n_circle / area_ha * flux_factor]                # recruits / ha / period

attach_reg <- function(iv, survey) {
  sid <- unique(iv$siteID_dt[, .(siteID, siteName)])
  sid[, c("Tnr", "Enr") := tstrsplit(siteName, "_", fixed = TRUE, type.convert = TRUE)]
  obs <- merge(iv$obs, sid[, .(siteID, Tnr, Enr)], by = "siteID", all.x = TRUE)
  obs <- merge(obs, reg_tab[interval_name == survey, .(Tnr, Enr, species_name, reg_new = reg)],
               by = c("Tnr", "Enr", "species_name"), all.x = TRUE)
  obs[, reg := fifelse(is.na(reg_new), 0, reg_new)]
  obs[, c("reg_new", "Tnr", "Enr") := NULL]
  obs
}
i1$obs <- attach_reg(i1, "bwi2012")
i2$obs <- attach_reg(i2, "bwi2022")

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Stitch: one species map, unique siteIDs, pooled z-scaling
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
recode_species <- function(dt) {
  if ("species" %in% names(dt)) dt[, species := NULL]
  dt <- merge(dt, species_dt, by = "species_name", all.x = TRUE)
  stopifnot(!anyNA(dt$species))
  dt
}
off <- max(i1$obs$siteID)
qualify <- function(iv, k, offset) {
  iv$obs[, siteID := siteID + offset]
  iv$env[, siteID := siteID + offset]
  iv$init[, siteID := siteID + offset]
  iv$siteID_dt[, siteID := siteID + offset]
  iv$siteID_dt[, `:=`(interval = k, corner = siteName,
                      siteName = paste0(siteName, "@i", k),
                      patchName = paste0(patchName, "@i", k))]
  iv$init[, `:=`(siteName = paste0(siteName, "@i", k),
                 patchName = paste0(patchName, "@i", k))]
  iv$trees[, siteName := paste0(siteName, "@i", k)]
  iv$design[, siteName := paste0(siteName, "@i", k)]
  iv
}
i1 <- qualify(i1, 1L, 0L)
i2 <- qualify(i2, 2L, off)

obs  <- recode_species(rbind(i1$obs, i2$obs, fill = TRUE))
init <- recode_species(rbind(i1$init, i2$init, fill = TRUE))
sids <- rbind(i1$siteID_dt, i2$siteID_dt, fill = TRUE)
env_unscaled <- rbind(i1$env, i2$env)
## tallied trees on the surviving site-periods (one patch per corner, obs at year P)
trees <- recode_species(rbind(i1$trees, i2$trees))
trees <- merge(trees, unique(sids[, .(siteID, siteName)]), by = "siteName")
trees[, `:=`(year = P, patchID = 1L)]
tree_keep <- c("siteID", "year", "patchID", "species", "dbh", "trees", "growth", "died")
setorder(trees, siteID, species, -dbh)
design <- merge(rbind(i1$design, i2$design), unique(sids[, .(siteID, siteName, interval, corner)]),
                by = "siteName")

## the single drought index: mean of the current and the previous year's
## anomaly on the natural scale, appended as the last column
env_unscaled[, drought2 := (get(PARAMS$drought2_parents[1]) + get(PARAMS$drought2_parents[2])) / 2]

## pooled z-scaling over both decades; management stays on its 0..3 scale.
## Row order of env_scales: the v3 rows (management last), then drought2.
scale_vars <- c(PARAMS$env_vars, "drought", "drought_lag1", "drought2")
sc_row <- function(v) data.table(variable = v, mean = mean(env_unscaled[[v]]), sd = sd(env_unscaled[[v]]))
env_scales <- rbindlist(c(
  lapply(c(PARAMS$env_vars, "drought", "drought_lag1"), sc_row),
  list(data.table(variable = "management", mean = 0, sd = 1)),
  list(sc_row("drought2"))))
env_dt <- copy(env_unscaled)
for (v in scale_vars) {
  s <- env_scales[variable == v]
  set(env_dt, j = v, value = (env_dt[[v]] - s$mean) / s$sd)
}

obs[, period_length := P]
obs_keep <- c("siteID", "year", "species", "ba", "dbh", "trees", "growth",
              "growth_n", "mort", "n_at_risk", "n_died", "reg", "period_length")
setorder(obs, siteID, year, species)
setorder(env_dt, siteID, year); setorder(env_unscaled, siteID, year)
setorder(sids, siteID)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Validate the stitched tables
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
n_sites <- uniqueN(obs$siteID)
stopifnot(
  identical(sort(unique(obs$siteID)), seq_len(n_sites)),          # contiguous, unique
  !anyDuplicated(sids, by = c("siteID", "year")),                 # one row per site x inventory
  uniqueN(sids[, .(siteID, siteName)]) == n_sites,
  setequal(unique(env_dt$siteID), unique(obs$siteID)),
  setequal(unique(init$siteID), unique(obs$siteID)),
  all(env_dt[, .N, by = siteID]$N == P + 1L),                       # years 0..P each
  !anyNA(env_dt),
  all(obs$year == P),
  all(obs$species %in% species_dt$species),
  max(i1$obs$siteID) < min(i2$obs$siteID),                          # decades never collide
  setequal(unique(trees$siteID), unique(obs$siteID)),               # every site has its tally
  all(trees$died %in% c(0, 1)),                                     # v8: no shading-only rows
  nrow(design) == n_sites
)
## decade-1 site-periods lost to fix 3 (3c): corners with the exact period whose
## start stand expands to no tree. With the lateral ingrowth in the start
## stand (v3 / the study's v5) they were site-periods; here they are not.
dropped <- i1$start_empty[siteName %chin% c1[period_length == P, siteName]]
stopifnot(!any(paste0(dropped$siteName, "@i1") %chin% sids$siteName))
message(sprintf("[sites] v8 %d site-periods: %d decade-1 site-periods dropped (start stand expands to no tree; %d of them without any tallied start tree)",
                n_sites, nrow(dropped), dropped[n_start_tally == 0L, .N]))

## the per-tree removal must reproduce the decade removal implied by obs_dt$mort
## (both are tallied-tree counts; obs weights trees per stem, so compare loosely)
tree_check <- merge(
  trees[!is.na(died), .(tree_removal = mean(died), at_risk_trees = .N), by = .(interval = fifelse(siteID <= off, 1L, 2L))],
  obs[n_at_risk > 0 & is.finite(mort), .(obs_removal = sum(n_at_risk * (1 - (1 - mort)^P)) / sum(n_at_risk)),
      by = .(interval = fifelse(siteID <= off, 1L, 2L))], by = "interval")
print(tree_check)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Write
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
fwrite(obs[, ..obs_keep],  file.path(OUT_NAT, "obs_dt.csv"))
fwrite(env_dt,             file.path(OUT_NAT, "env_dt.csv"))
fwrite(env_unscaled,       file.path(OUT_NAT, "env_unscaled_dt.csv"))
fwrite(env_scales,         file.path(OUT_NAT, "env_scales.csv"))
fwrite(init,               file.path(OUT_NAT, "init_cohort_dt.csv"))
fwrite(species_dt,         file.path(OUT_NAT, "species_dt.csv"))
fwrite(sids,               file.path(OUT_NAT, "siteID_dt.csv"))
fwrite(trees[, ..tree_keep], file.path(OUT_NAT, "tree_obs_dt.csv"))
fwrite(reg_tab,            file.path(OUT_NAT, "reg_target.csv"))
fwrite(i1$qc,              file.path(OUT_NAT, "wzp_qc_i1.csv"))
fwrite(i2$qc,              file.path(OUT_NAT, "wzp_qc_i2.csv"))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Spatially blocked split: v3's, with the rule for new tracts
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## one location per tract (corner_env x/y are the tract's WGS84 lon/lat)
sids[, Tnr := as.integer(sub("_.*$", "", corner))]
tr <- unique(env_in[Tnr %in% sids$Tnr, .(Tnr, lon = x, lat = y)], by = "Tnr")
xy <- sf::st_coordinates(sf::st_transform(
  sf::st_as_sf(tr, coords = c("lon", "lat"), crs = 4326), 3035))
tr[, `:=`(east = xy[, 1], north = xy[, 2])]
bm <- PARAMS$block_km * 1000
tr[, block := sprintf("E%03d_N%03d", floor(east / bm), floor(north / bm))]

## Every tract keeps the set, block and subset membership it has in
## two-interval-v3/split_blocks.csv. A tract that is new in v8 (a decade-2
## corner that passes the corner filters only with the fixed fates, in a tract
## without v3 site-periods) takes the set of its v3 block and joins neither
## train_quarter nor pilot; a block unknown to v3 stops the script (it would
## need a new holdout decision).
v3_split <- fread(file.path(PARAMS$v3_dir, "split_blocks.csv"))
v3_blocks <- unique(v3_split[, .(block, set)])
stopifnot(!anyDuplicated(v3_blocks$block))
tr <- merge(tr, v3_split[, .(Tnr, block_v3 = block, set, train_quarter, pilot)], by = "Tnr", all.x = TRUE)
new_tracts <- tr[is.na(block_v3), Tnr]
stopifnot(identical(tr[!is.na(block_v3), block], tr[!is.na(block_v3), block_v3]),
          all(tr[is.na(block_v3), block] %in% v3_blocks$block))
set_of_block <- setNames(v3_blocks$set, v3_blocks$block)
tr[is.na(block_v3), `:=`(set = unname(set_of_block[block]),
                         train_quarter = FALSE, pilot = FALSE)]
tr[, block_v3 := NULL]
message(sprintf("[split] %d tracts, %d of them new in v8 (set from their v3 block)",
                nrow(tr), length(new_tracts)))
fwrite(tr[, .(Tnr, lon, lat, east, north, block, set, train_quarter, pilot)],
       file.path(OUT_VER, "split_blocks.csv"))

stopifnot(all(sids$Tnr %in% tr$Tnr), !anyDuplicated(tr$Tnr))
site_set <- merge(unique(sids[, .(siteID, interval, Tnr)]), tr[, .(Tnr, block, set, train_quarter, pilot)],
                  by = "Tnr")
stopifnot(nrow(site_set) == n_sites,
          site_set[, uniqueN(set), by = Tnr][, all(V1 == 1L)])   # tracts never split
blocks  <- sort(unique(tr$block))
hold_bl <- sort(unique(tr[set == "holdout_blocked", block]))
stopifnot(!any(tr[set == "train", block] %in% hold_bl))

## distance from every holdout tract to its nearest training tract (km)
hold_xy  <- as.matrix(tr[set == "holdout_blocked", .(east, north)])
train_xy <- as.matrix(tr[set == "train", .(east, north)])
d_km <- FNN::get.knnx(train_xy, hold_xy, k = 1)$nn.dist[, 1] / 1000
hold_dist <- c(median = median(d_km), min = min(d_km), p10 = unname(quantile(d_km, 0.10)))

set_sites <- list(
  train           = site_set[set == "train", siteID],
  holdout_blocked = site_set[set == "holdout_blocked", siteID],
  train_quarter   = site_set[train_quarter == TRUE, siteID],
  pilot           = site_set[pilot == TRUE, siteID])
stopifnot(identical(names(set_sites), PARAMS$sets))

#' Write the FINN tables for a subset of national siteIDs, renumbered 1..n.
write_set <- function(name, site_ids) {
  d  <- file.path(OUT_VER, name)
  dir.create(d, showWarnings = FALSE)
  mp <- data.table(siteID = sort(site_ids))[, siteID_new := .I]
  pick <- function(x) {
    x <- merge(x[siteID %in% site_ids], mp, by = "siteID")
    x[, siteID := siteID_new][, siteID_new := NULL]
    x
  }
  o  <- pick(obs[, ..obs_keep]); setorder(o, siteID, year, species)
  e  <- pick(env_dt);            setorder(e, siteID, year)
  eu <- pick(env_unscaled);      setorder(eu, siteID, year)
  s  <- merge(sids[siteID %in% site_ids][, siteID_national := siteID],
              mp, by = "siteID")[, siteID := siteID_new][, siteID_new := NULL]
  s[, Tnr := NULL]; setorder(s, siteID)
  fwrite(o,              file.path(d, "obs_dt.csv"))
  fwrite(e,              file.path(d, "env_dt.csv"))
  fwrite(eu,             file.path(d, "env_unscaled_dt.csv"))
  fwrite(env_scales,     file.path(d, "env_scales.csv"))
  fwrite(pick(init),     file.path(d, "init_cohort_dt.csv"))
  tt <- pick(trees[, ..tree_keep]); setorder(tt, siteID, species, -dbh)
  fwrite(tt,             file.path(d, "tree_obs_dt.csv"))
  fwrite(species_dt,     file.path(d, "species_dt.csv"))
  fwrite(s,              file.path(d, "siteID_dt.csv"))
  data.table(set = name, sites = nrow(mp),
             interval_1 = s[, uniqueN(siteID[interval == 1L])],
             interval_2 = s[, uniqueN(siteID[interval == 2L])],
             tracts = uniqueN(sub("_.*$", "", s$corner)))
}
set_tab <- rbindlist(lapply(PARAMS$sets, function(s) write_set(s, set_sites[[s]])))
set_tab[, blocks := c(length(blocks) - length(hold_bl), length(hold_bl),
                      uniqueN(tr[train_quarter == TRUE, block]),
                      uniqueN(tr[pilot == TRUE, block]))]
sids[, Tnr := NULL]
print(set_tab)
message(sprintf("[split] holdout -> nearest training tract: median %.1f km, min %.1f km, p10 %.1f km",
                hold_dist[["median"]], hold_dist[["min"]], hold_dist[["p10"]]))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Checks (for PREPARATION.md)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
ALL_SETS <- c("national", PARAMS$sets)
members  <- rbind(data.table(set = "national", siteID = seq_len(n_sites)),
                  rbindlist(lapply(PARAMS$sets, function(s) data.table(set = s, siteID = set_sites[[s]]))))
members  <- merge(members, unique(sids[, .(siteID, interval, corner)]), by = "siteID")
A <- PARAMS$patch_size_ha

## (a) state per site-period: initial stand (init_cohort_dt), tallied trees
## (tree_obs_dt), end state (obs_dt year 10), the tallied design estimates and
## BWI2's own 2002 tally.
init_site <- function(x) x[living %in% TRUE, .(init_ba = sum(trees * bhd2ba(dbh)) / A,
                                               init_stems = sum(trees) / A), by = siteID]
end_site  <- function(o) o[, .(end_ba = sum(ba) / A, end_stems = sum(trees) / A), by = siteID]
## BWI2's own 2002 angle-count tally (b2_baeume; Av 0 = in the 2002 WZP), dbh >=
## 7 cm, the 15 species, per corner
b2 <- bwi("b2_baeume", select = c("Tnr", "Enr", "Av", "Ba", "Bhd", "G", "N_ha"))
b2 <- merge(b2, splk, by = "Ba", all.x = TRUE)
b2 <- b2[Av == AV_WZP & is.finite(Bhd) & Bhd >= WZP_MIN_BHD * 10 & is.finite(N_ha) & N_ha > 0 &
         species_name %chin% V, .(b2_ba = sum(N_ha * G), b2_stems = sum(N_ha)),
         by = .(corner = paste(Tnr, Enr, sep = "_"))]
S <- Reduce(function(a, b) merge(a, b, by = "siteID", all.x = TRUE), list(
  unique(sids[, .(siteID, interval, corner)]),
  init_site(init), end_site(obs),
  trees[, .(tally_ba = sum(trees * bhd2ba(dbh)) / A, tally_stems = sum(trees) / A), by = siteID],
  design[, .(siteID, ba_start_design, stems_start_design, ba_end_design, stems_end_design)]))
S <- merge(S, b2, by = "corner", all.x = TRUE)
S[interval == 1L & is.na(b2_ba), `:=`(b2_ba = 0, b2_stems = 0)]
num <- setdiff(names(S), c("siteID", "interval", "corner"))
for (v in num) S[is.na(get(v)), (v) := 0]
SM <- merge(members[, .(set, siteID)], S, by = "siteID")
STATE <- SM[, .(
  site_periods            = .N,
  init_ba                 = mean(init_ba),
  `init BA/start design`  = sum(init_ba) / sum(ba_start_design),
  `init BA/b2_baeume`     = if (interval[1] == 1L) sum(init_ba) / sum(b2_ba) else NA_real_,
  `tally BA/b2_baeume`    = if (interval[1] == 1L) sum(tally_ba) / sum(b2_ba) else NA_real_,
  end_ba                  = mean(end_ba),
  `end BA/design`         = sum(end_ba) / sum(ba_end_design),
  `end stems/design`      = sum(end_stems) / sum(stems_end_design),
  `median site abs(end BA/design - 1)` = median(abs(end_ba / ba_end_design - 1)[ba_end_design > 0])
), by = .(set, interval)]
STATE[, set := factor(set, ALL_SETS)]
setorder(STATE, set, interval)
print(STATE)
nat1 <- STATE[set == "national" & interval == 1L]
nat2 <- STATE[set == "national" & interval == 2L]

## (b) basal-area conservation of the start and end expansions
QC <- rbindlist(list(
  data.table(interval = 1L, stand = "start (2002 tally)", i1$qc),
  data.table(interval = 1L, stand = "end (2012 tally, fresh)", i1$qc_end),
  data.table(interval = 2L, stand = "start (2012 tally)", i2$qc),
  data.table(interval = 2L, stand = "end (2022 tally, fresh)", i2$qc_end)), use.names = TRUE)
QC <- QC[, .(patches = .N, `BA expanded / design` = sum(ba_ha_expanded) / sum(ba_ha_design),
             `median abs rel err` = median(abs(rel_err), na.rm = TRUE),
             `p95 abs rel err` = as.numeric(quantile(abs(rel_err), 0.95, na.rm = TRUE))), by = .(interval, stand)]
print(QC)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## PREPARATION.md
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
md_tab <- function(d, digits = 4) {
  d <- copy(as.data.table(d))
  for (v in names(d)) {
    if (is.factor(d[[v]])) set(d, j = v, value = as.character(d[[v]]))
    if (is.double(d[[v]]))
      set(d, j = v, value = ifelse(is.na(d[[v]]), "",
                                   formatC(d[[v]], format = "fg", digits = digits, flag = "#")))
  }
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    apply(d, 1, function(r) paste0("| ", paste(trimws(r), collapse = " | "), " |")))
}
commit <- tryCatch(suppressWarnings(system2("git", c("describe", "--always", "--dirty"),
                                            stdout = TRUE, stderr = FALSE)),
                   error = function(e) character(0))
if (length(commit) != 1L || !nzchar(commit)) commit <- "unknown"
per_iv <- function(iv, k) {
  o <- obs[siteID %in% iv$obs$siteID]
  sprintf("| %d | %d | %d | %.2f | %.4f | %.4f | %.2f | %.1f %% |",
          k, iv$diag$sites, nrow(o), mean(o$ba), mean(o$mort, na.rm = TRUE),
          median(o$growth, na.rm = TRUE), mean(o$reg),
          100 * iv$diag$drought_cover)
}
tier_tab <- function(iv) {
  e <- env_unscaled[siteID %in% iv$obs$siteID & year == 0L]
  paste(e[, .N, by = management][order(management),
          sprintf("M%d %d", as.integer(management), N)], collapse = ", ")
}

md <- c(
  sprintf("# Prepared data: %s", PARAMS$version),
  "",
  sprintf("Written by `scripts/71_prepare_two_interval_v8.R` (decade-2 panel: `scripts/91_import_new_schema_chained.R`) at commit `%s` on %s.",
          commit, format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "This is the set the FINN-bwi study fits. The study's own files are listed with their SHA-256",
  "in `outputs/prepared/two-interval-v8.sha256`; the script header says which files match.",
  "",
  "## What changed vs two-interval-v3 (`scripts/70`)",
  "",
  "1. **Freshly dead trees** (flag `tot` in ACCESS2003 `b2_wzp` / `b3_wzp`, `damage_dead` in the",
  "   archive): dead at the start survey -> not in the start stand; a survivor flagged at the end",
  "   survey -> natural death; a new tree already dead -> not a recruit. Decade 1 here, decade 2 in",
  "   `scripts/91`.",
  "2. **CI2017 losses** (decade 2, `scripts/91`): where BWI2022 refers to the 2017 carbon inventory,",
  "   trees lost in 2012-2017 keep their loss as the decade fate instead of vanishing.",
  "3. **Lateral ingrowth counted once.** (a) Initial stand = the start tally only: decade 1 no longer",
  "   expands the back-modelled 2002 state of laterally ingrown trees into the 2002 stand, and",
  "   `tree_obs_dt.csv` has no shading-only rows (`died = NA`). (b) End state (`ba`, `trees`, `dbh` at",
  "   year 10) = every living tree of the end tally expanded at its own end-survey n_ha,",
  "   largest-remainder allocation per patch x species, both decades. The t0 and t1 virtual trees",
  "   are no longer chained one-to-one; growth and mortality come from the tallied trees.",
  sprintf("   (c) Consequence: **%d decade-1 site-periods drop out**, because their start stand expands to no tree at all: %d corners whose 2002 angle count holds no tree of the 15 species (v3's start stand there was back-modelled lateral ingrowth only, mean %.1f m2/ha; end tally mean %.1f m2/ha) and %d whose few tallied trees round to 0 stems on the 0.1 ha patch. `makeObsData()` needs both censuses.",
          nrow(dropped), dropped[n_start_tally == 0L, .N],
          dropped[n_start_tally == 0L, mean(ba_lateral_t0)], dropped[, mean(ba_end_tally)],
          dropped[n_start_tally > 0L, .N]),
  "4. **`drought2`** = (drought + drought_lag1) / 2 on the natural scale, pooled-scaled over the",
  "   national set, last column of the env tables and last row of `env_scales.csv`.",
  "",
  "The tallied-tree sample (growth, mortality, `tree_obs_dt.csv`) is selected with the start",
  "allocation that still contains the lateral ingrowth, so the sample does not move with the",
  sprintf("rounding of near-zero expansions: with the v8 start allocation alone, %d (decade 1) / %d (decade 2) at-risk tallied trees would move in or out of the sample. End-survey rows of dead and harvested trees keep their start expansion (they are not part of the end stand).",
          i1$diag$n_sample_flip, i2$diag$n_sample_flip),
  "",
  "## Regeneration target: a documented choice",
  "",
  sprintf("`obs_dt$reg` = stems of archive size classes %s (5.0-6.9 cm) on the r = %.2f m circle, turned from stock into a flux of recruits per ha and period at an assumed diameter increment of %.2f cm/yr through the %.1f cm band (factor %.3f per stem/ha and period). The increment sets the level only, not the composition. This is v8's default; it is about twice the decade-1 recruit flux implied by the angle-count sample. Alternative targets (a lower increment, recruits from the angle-count sample, a hybrid of both) are being evaluated and are not part of this repository.",
          paste(PARAMS$reg_classes, collapse = "+"), PARAMS$reg_circle_r_m,
          PARAMS$reg_incr_cm_yr, PARAMS$reg_band_cm, flux_factor),
  sprintf("Off-V share of the circle stems: %s.",
          paste(off_V[, sprintf("%s %.1f %%", interval_name, 100 * share_off_V)], collapse = ", ")),
  "",
  "## Checks",
  "",
  "### Initial stand and end state, per set and decade",
  "",
  "Per hectare, from the prepared tables: `init` = living trees of `init_cohort_dt` (expanded),",
  "`tally` = `tree_obs_dt` (tallied trees at their stems), `end` = `obs_dt` at year 10, summed over",
  "species; `design` = the tally expanded per tallied tree at its n_ha (no rounding);",
  "`b2_baeume` = BWI2's own 2002 tally (Av 0, dbh >= 7 cm, the 15 species) on the same corners.",
  "Ratios are ratios of sums over the site-periods; `*_ba` columns are means in m2/ha.",
  "",
  md_tab(STATE), "",
  sprintf("National: decade-1 init / BWI2's 2002 tally %.3f (tallied %.3f); end-state BA / design estimate %.4f (decade 1) / %.4f (decade 2) -- identical by construction up to the integer rounding of the expansion.",
          nat1$`init BA/b2_baeume`, nat1$`tally BA/b2_baeume`,
          nat1$`end BA/design`, nat2$`end BA/design`),
  "",
  sprintf("Open decision, not applied: the decade-1 corner filter of `scripts/02` (start BA >= 5 m2/ha, the 15 species >= 99 %% of it) was evaluated on the start stand with the lateral ingrowth. Evaluated on the start tally, %d of v8's %d decade-1 site-periods hold < 5 m2/ha of the 15 species (decade 2's filter already sees its plain 2012 tally). v8 keeps the corner selection.",
          S[interval == 1L & ba_start_design < PARAMS$min_ba_ha_t0, .N], S[interval == 1L, .N]),
  "",
  "### Basal-area conservation of the expansions (national)",
  "",
  "`wzp_qc_i{1,2}.csv` hold the start expansion; the end expansion is new in v8.", "",
  md_tab(QC), "",
  "### Per-tree removal vs obs_dt (tallied trees)",
  "",
  "| interval | at-risk tallied trees | removal share, trees | removal share, obs_dt (per stem) |",
  "|---|---|---|---|",
  tree_check[, sprintf("| %d | %d | %.3f | %.3f |", interval, at_risk_trees, tree_removal, obs_removal)],
  "",
  "### Pooled scaling",
  "",
  md_tab(env_scales, digits = 6), "",
  "## What this is",
  "",
  "FINN inputs for **two decades on the same corners**: interval 1 = BWI2 -> BWI3",
  "(2002 -> 2012, ACCESS2003 release) and interval 2 = BWI3 -> BWI4 (2012 -> 2022,",
  "Thuenen inventory archive). Annual timestep, observation at year 10 of each",
  "decade, per-year drought anomaly in `env_dt$drought` (+ `drought_lag1`, `drought2`).",
  "Each decade is its own FINN site-period; `siteID_dt.csv` maps `siteID` to `interval`",
  sprintf("and the physical `corner` (Tnr_Enr). Interval-2 siteIDs are offset by %d.", off),
  "",
  "## Contents (`national/`; every set folder has the same FINN tables)",
  "",
  "| file | content |",
  "|---|---|",
  "| obs_dt.csv | siteID, year, species, ba, dbh, trees (end state: fresh expansion of the living end tally), growth, growth_n, mort, n_at_risk, n_died (tallied trees), reg (circles), period_length |",
  "| env_dt.csv | siteID, year 0..10, mat_c, map_mm, elev_m, x, y, management, drought, drought_lag1, drought2 -- z-scaled ONCE over both decades; management left on 0..3 |",
  "| env_unscaled_dt.csv / env_scales.csv | natural scale and the pooled mean/sd |",
  "| init_cohort_dt.csv | start-of-decade stand = the start tally (BWI2 for i1, BWI3 for i2), for makeInitCohorts() |",
  "| species_dt.csv | the 15 species, one code map for both decades |",
  "| siteID_dt.csv | siteID -> interval, corner, qualified siteName (`Tnr_Enr@i1/@i2`) |",
  "| tree_obs_dt.csv | tallied trees of the start tally: siteID, year, patchID, species, dbh (start, cm), trees (stems on the 0.1 ha patch), growth (/yr, survivors), died (removal 0/1 over the decade; no NA rows) |",
  "| reg_target.csv | regeneration target per survey x corner x species (diagnostic) |",
  "| wzp_qc_i{1,2}.csv | basal-area conservation of the start expansion |",
  "",
  "## Per decade",
  "",
  "| interval | sites | obs rows | mean ba (per 0.1 ha patch) | mean mort (/yr) | median growth (/yr) | mean reg (/ha/period) | drought coverage |",
  "|---|---|---|---|---|---|---|---|",
  per_iv(i1, 1L), per_iv(i2, 2L),
  "",
  sprintf("Management tiers (sites): interval 1 %s; interval 2 %s.", tier_tab(i1), tier_tab(i2)),
  "",
  sprintf("Lateral ingrowth that v3 carried in the decade-1 start stand: %.1f %% of its start basal area (all species of the panel, before the restriction to the 15).",
          100 * i1$diag$lat_share),
  "",
  "Interval-2 corner funnel:",
  "",
  "| step | corners |", "|---|---|",
  funnel[, sprintf("| %s | %d |", step, n)],
  "",
  sprintf("Calamity salvage booked as natural death in interval 2 (harvest_reason %s): %d of %d harvested stems (%.1f %%).",
          paste(PARAMS$calamity_reasons, collapse = ","), n_salvage, n_harvest_all,
          100 * n_salvage / n_harvest_all),
  "",
  "## Sets (spatially blocked, v3's split)",
  "",
  sprintf("Every tract keeps its v3 set, block and train_quarter / pilot membership (`two-interval-v3/split_blocks.csv`); %d tract(s) new in v8 take the set of their v3 block and join neither train_quarter nor pilot. Whole tracts, both decades together, %d km blocks on EPSG:3035, %d of %d blocks held out. Every set folder has the national file layout with siteIDs renumbered 1..n (`siteID_dt$siteID_national` maps back) and the national pooled scaling. `split_blocks.csv` lists v8's %d tracts.",
          length(new_tracts), PARAMS$block_km, length(hold_bl), length(blocks), nrow(tr)),
  "",
  "| set | site-periods | interval 1 | interval 2 | tracts | blocks |", "|---|---|---|---|---|---|",
  set_tab[, sprintf("| %s | %d | %d | %d | %d | %d |", set, sites, interval_1,
                    interval_2, tracts, blocks)],
  "",
  sprintf("Distance from each holdout tract to the nearest training tract: median %.1f km, minimum %.1f km, 10th percentile %.1f km.",
          hold_dist[["median"]], hold_dist[["min"]], hold_dist[["p10"]]),
  "",
  "## Caveats",
  "",
  "- **The end state is a design estimate, not FINN's accounting of the start stand.** A survivor",
  "  stands for 4 / g(t0) stems at the start and 4 / g(t1) at the end (it grew, so fewer); the",
  "  expected stand still balances (start stems - deaths + recruits = end stems in expectation),",
  "  but no longer per site. The state targets are therefore noisier per site-period than v3's",
  "  carried construction, and unbiased instead of inflated.",
  "- **The regeneration target is a choice** (section above). With the inflated end state gone, an",
  "  over-prediction of stems by a fitted model points at the increment assumption, not at ingrowth.",
  "- **Regeneration geometry for 2022** is assumed unchanged (r = 1.75 m for classes 5+6).",
  "- **Tier M1 is near-empty in both decades** (BAF 4: one tallied tree carries 4 m2/ha).",
  "- Environment exists only for interval-1 corners.",
  "",
  "## PARAMS",
  "",
  "```",
  capture.output(str(PARAMS, give.attr = FALSE)),
  "```"
)
writeLines(md, file.path(OUT_VER, "PREPARATION.md"))

message(sprintf("\nDone -> %s/", OUT_VER))
message(sprintf("  sites      : %d  (interval 1: %d, interval 2: %d)",
                n_sites, i1$diag$sites, i2$diag$sites))
message(sprintf("  obs rows   : %d", nrow(obs)))
message(sprintf("  env rows   : %d", nrow(env_dt)))
message(sprintf("  init trees : %d", nrow(init)))
