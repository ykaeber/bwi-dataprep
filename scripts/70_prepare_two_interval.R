#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-bwi | scripts/70_prepare_two_interval.R
##
## Two-interval FINN inputs: the BWI2 -> BWI3 decade (2002 -> 2012) AND the
## BWI3 -> BWI4 decade (2012 -> 2022), on an annual timestep with the per-year
## drought covariate. Each decade enters FINN as its own site-period: the same
## physical corner appears twice, with a different initial stand, different
## environment years and a unique siteID.
##
## WHY TWO INTERVALS
## The second decade contains the 2018-2020 drought. With one mortality
## observation per site per decade, a single-decade fit cannot separate "this
## species dies more" from "this decade was dry"; two decades on the same
## corners can, and that is what makes the annual drought covariate earn its
## timestep.
##
## ONE CODE PATH FOR BOTH DECADES
## prepare_interval() below mirrors scripts/04_prepare_finn_inputs.R, steps 1-7b
## (angle-count expansion, makeObsData, tallied sample sizes, total-removal
## mortality, survivor-only growth, annual grid, drought) and is called once per
## decade, so the two decades cannot diverge. It returns UNSCALED tables; the
## stitch at the end z-scales once over both decades together, so the same
## climate maps to the same z-score in either decade.
##
## DECISIONS MADE HERE (values in PARAMS)
##  * Interval 2 corners pass the corner filters scripts/02 applied to interval
##    1: basal area >= 5 m2/ha at the start, the 15 species >= 99 % of basal
##    area, realised period 8-12 yr; corners with untraceable trees are dropped
##    (interval 1 drops them as "not usable"). They must also have environment
##    in corner_env.csv, which makes them a subset of the interval-1 corners.
##  * Interval 2 survey years come from plot.acquisition_date (earliest per
##    tract, as bwi_periods() does from Datum_feld for interval 1), not from the
##    nominal 2012/2022, so the drought years line up with the real field work.
##  * Interval 2 management tier = basal-area share removed by REGULAR harvest,
##    with the thresholds of scripts/01. Trees removed in calamity salvage
##    (plot.harvest_reason 1-4: storm, insects/drought, fire, uncertain
##    calamity) are booked as natural deaths, so the 2018-2020 drought pulse
##    stays a biological signal instead of being absorbed by the management
##    covariate. BWI2012 records no harvest_reason, so interval 1 keeps all of
##    its harvest as harvest (this includes the Kyrill storm of 2007).
##  * Regeneration target for BOTH decades from the inventory-archive
##    regeneration circles (size classes 5+6 = dbh 5.0-6.9 cm, radius 1.75 m),
##    per HECTARE per period: stock in the 5.0-6.9 cm band -> flux through 7 cm
##    via an assumed diameter increment. obs_dt$reg is per hectare
##    (FINN::makeObsData: recruits / plotsize); an earlier version of this
##    target (FINN-bwi scripts/55) also multiplied by the 0.1 ha patch area,
##    leaving it 10x too low. The archive's bwi2012 counts equal ACCESS2003
##    b3_jung exactly, so interval 1 loses nothing by switching source.
##
## Inputs
##   data/bwi_csv/*                        ACCESS2003 tables (bwi_tree_panel,
##                                         species_lookup)
##   data/raw/bwi4/{plot,regeneration}.csv inventory archive (scripts/fetch_bwi4.sh)
##   outputs/new_schema/panel.csv          BWI3->BWI4 tree panel (scripts/40)
##   outputs/calibration_set_corners.csv   interval-1 corners, all tiers (scripts/02)
##   outputs/calibration_set_species.csv   the 15 species (scripts/02)
##   outputs/corner_env.csv                climate, elevation, x/y (scripts/03)
##   outputs/site_drought_annual.csv       per-year drought anomaly (scripts/16)
##
## Outputs  outputs/prepared/two-interval-v1/
##   PREPARATION.md          what was done, decisions, PARAMS, counts, caveats
##   national/obs_dt.csv     observations, both decades (obs at year 10)
##   national/env_dt.csv     covariates per siteID x year 0..10, z-scaled once
##                           over both decades (management kept on 0..3)
##   national/env_unscaled_dt.csv, env_scales.csv
##   national/init_cohort_dt.csv   start-of-decade stands for both decades
##   national/species_dt.csv       the 15 species, one code map for both decades
##   national/siteID_dt.csv        siteID <-> interval, physical corner (crosswalk)
##   national/reg_target.csv       regeneration target per corner x species
##   national/wzp_qc_i{1,2}.csv    basal-area conservation of the expansion
##   train/, holdout_blocked/, train_quarter/, pilot/   the same FINN tables for
##                           the spatially blocked split (train_quarter ~25 % and
##                           pilot ~2,000 site-periods: nested subsets of train);
##                           siteIDs renumbered 1..n, siteID_dt keeps
##                           siteID_national. Scaling is the national one.
##   split_blocks.csv        tract -> block -> set
##
## Helpers sourced (none of them changes what the data mean):
##   R/config.R    paths and constants (WZP_BAF, PATCH_SIZE_HA, bhd2ba)
##   R/bwi_io.R    bwi(): read a cached ACCESS2003 table
##   R/wzp.R       angle-count expansion and tallied sample sizes
##   R/bwi_trees.R bwi_tree_panel(): the ACCESS2003 BWI2->BWI3 tree panel
##
## Run from the project root:  Rscript scripts/70_prepare_two_interval.R
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
  version = "two-interval-v1",

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

  ## Management tiers exactly as scripts/01: M0 = no tallied harvest, then the
  ## basal-area share removed: <= 5 % M1, <= 20 % M2, above M3 -> 0..3.
  tier_cut = c(M1 = 0.05, M2 = 0.20),
  ## plot.harvest_reason codes booked as natural death in interval 2:
  ## 1 storm, 2 insects/drought, 3 other calamity (fire), 4 calamity, cause
  ## uncertain. 0 regular and 5 unknown stay harvest.
  calamity_reasons = 1:4,

  ## Regeneration target (stock -> flux, as scripts/55): stems in the
  ## 5.0-6.9 cm band on the r = 1.75 m circle, crossing 7 cm at reg_incr_cm_yr.
  reg_classes    = c(5L, 6L),
  reg_circle_r_m = 1.75,
  reg_incr_cm_yr = 0.30,          # sets the MAGNITUDE only; composition is free of it
  reg_band_cm    = 1.9,

  ## Spatially blocked holdout (r-pipeline-conventions §5). The unit that moves
  ## between sets is a whole TRACT with both its decades, so no corner trains on
  ## one decade while its other decade is held out. Tracts fall into square
  ## blocks on the ETRS89-LAEA grid (EPSG:3035, metres). On the 4 km tract grid,
  ## holdout -> nearest training tract (median): 50 km blocks 10.8 km, 75 km
  ## 18.4 km, 100 km 22.1 km (53 blocks, 11 held out), 150 km 40 km but the
  ## holdout jumps to 39 % of tracts from 6 blocks. 100 km (~1 degree) balances
  ## separation against a stable holdout share.
  block_km      = 100,
  holdout_share = 0.20,           # share of BLOCKS held out
  pilot_target  = 2000L,          # site-periods in the pilot subset of train
  quarter_share = 0.25,           # train_quarter: share of train site-periods (faster cluster fits)

  bwi4_dir = file.path(BWI_RAW_DIR, "bwi4"),
  seed     = 42L
)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

set.seed(PARAMS$seed)
P        <- PARAMS$period
OUT_VER  <- file.path(OUT_DIR, "prepared", PARAMS$version)
OUT_NAT  <- file.path(OUT_VER, "national")
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
## Returns obs, env (unscaled), init, siteID_dt -- with siteIDs 1..n -- and QC.
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
  ## share of start-of-decade basal area carried by back-modelled lateral
  ## ingrowth (interval 1 only; see the caveat in PREPARATION.md)
  t0 <- panel[time == 0L]
  lat_share <- if ("ingrowth_lateral" %in% t0$fate)
    t0[, sum((fate == "ingrowth_lateral") * n_ha * bhd2ba(dbh)) / sum(n_ha * bhd2ba(dbh))]
  else 0
  panel <- panel[species_name %in% V]
  message(sprintf("[%s] restricted to V: %d -> %d records", label, n0, nrow(panel)))

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
  panel[, patchKey := paste(siteName, patchName, sep = "/")]
  alloc <- panel[time == 0L, .(treeKey, patchKey, species_name, n_ha, dbh)]
  invisible(wzp_allocate(alloc, by = c("patchKey", "species_name"),
                         patch_size_ha = PARAMS$patch_size_ha))
  qc <- wzp_ba_conservation(alloc, by = "patchKey",
                            patch_size_ha = PARAMS$patch_size_ha)
  message(sprintf("[%s] wzp basal-area conservation: median |err| %.2f %%",
                  label, 100 * median(abs(qc$rel_err), na.rm = TRUE)))
  rec <- panel[time == 1L & fate == "recruit",
               .(treeKey, patchKey, species_name, n_ha, dbh)]
  if (nrow(rec)) {
    invisible(wzp_allocate(rec, by = c("patchKey", "species_name"),
                           patch_size_ha = PARAMS$patch_size_ha))
    alloc <- rbind(alloc, rec)
  }
  panel <- merge(panel, alloc[, .(treeKey, n_trees)], by = "treeKey", all.x = TRUE)
  panel[is.na(n_trees), n_trees := 0L]
  virtual <- wzp_expand(panel,
                        keep = c("treeKey", "siteName", "patchName",
                                 "species_name", "dbh", "time", "status",
                                 "living", "mort", "fate", "period_length"),
                        index_by = c("treeKey", "time"))

  ## ---- 4) tree_dt in FINN's schema, then makeObsData() --------------------
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
  tal <- merge(panel[time == 0L & n_trees > 0L,
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
  drought_cover <- mean(is.finite(grid$drought))
  grid[!is.finite(drought), drought := 0]
  if (!"siteName" %in% names(env_dt))
    env_dt <- merge(env_dt, unique(siteID_dt[, .(siteID, siteName)]),
                    by = "siteID", all.x = TRUE)
  env_dt <- merge(env_dt, grid[, .(siteName, year, drought)],
                  by = c("siteName", "year"), all.x = TRUE)
  env_dt[!is.finite(drought), drought := 0]
  message(sprintf("[%s] drought attached: %.1f %% of site-years covered, driest model-year mean %.2f",
                  label, 100 * drought_cover,
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

  list(obs = obs_dt, env = env_dt[, c("siteID", "year", env_cols, "drought"), with = FALSE],
       init = init_trees, siteID_dt = siteID_dt, qc = qc,
       diag = list(sites = length(keep_sites), lat_share = lat_share,
                   drought_cover = drought_cover))
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
i1 <- prepare_interval(p1, c1, "interval 1 (2002->2012)", V, env_in, drought)
rm(p1); gc(verbose = FALSE)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Interval 2: BWI3 -> BWI4 (inventory archive) -- corners built here
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
p2 <- fread(file.path(OUT_DIR, "new_schema", "panel.csv"))
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
## Regeneration target, both decades, from the archive circles
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

## Cross-check: archive bwi2012 counts vs ACCESS2003 b3_jung (must be identical)
j3 <- bwi("b3_jung")[Gr %in% PARAMS$reg_classes, .(acc = sum(Anz)), by = .(Tnr, Enr)]
a3 <- rg[interval_name == "bwi2012", .(arch = sum(tree_count)), by = .(Tnr, Enr)]
x3 <- merge(a3, j3, by = c("Tnr", "Enr"))
reg_check <- sprintf("%d corners in both; identical on %.1f %%; totals %d vs %d",
                     nrow(x3), 100 * mean(x3$arch == x3$acc), sum(x3$arch), sum(x3$acc))
message("[reg] archive bwi2012 vs b3_jung (classes 5+6): ", reg_check)

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
## Regression check: interval 1 vs the validated single-decade inputs
## (outputs/finn, scripts/04 all tiers + management). Everything except reg
## (unit fix) should be identical; drought does not touch observations.
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
reg_cmp <- tryCatch({
  old  <- fread(file.path(OUT_DIR, "finn", "obs_dt.csv"))
  osid <- unique(fread(file.path(OUT_DIR, "finn", "siteID_dt.csv"))[, .(siteID, siteName)])
  old  <- merge(old, osid, by = "siteID")
  osp  <- fread(file.path(OUT_DIR, "finn", "species_dt.csv"))
  old  <- merge(old, osp, by = "species")
  new  <- copy(i1$obs)                                   # carries siteName from step 6
  m <- merge(old, new, by = c("siteName", "species_name"), suffixes = c(".old", ".new"))
  vars <- c("ba", "dbh", "trees", "growth", "mort", "n_at_risk", "n_died")
  d <- sapply(vars, function(v) max(abs(m[[paste0(v, ".old")]] - m[[paste0(v, ".new")]]),
                                    na.rm = TRUE))
  sprintf("sites old %d / new %d / matched %d; max |diff|: %s",
          uniqueN(old$siteName), uniqueN(new$siteName), uniqueN(m$siteName),
          paste(sprintf("%s %.2g", vars, d), collapse = ", "))
}, error = function(e) paste("not run:", conditionMessage(e)))
message("[regression] interval 1 vs outputs/finn: ", reg_cmp)

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
  iv
}
i1 <- qualify(i1, 1L, 0L)
i2 <- qualify(i2, 2L, off)

obs  <- recode_species(rbind(i1$obs, i2$obs, fill = TRUE))
init <- recode_species(rbind(i1$init, i2$init, fill = TRUE))
sids <- rbind(i1$siteID_dt, i2$siteID_dt, fill = TRUE)
env_unscaled <- rbind(i1$env, i2$env)

## pooled z-scaling over both decades; management stays on its 0..3 scale
scale_vars <- c(PARAMS$env_vars, "drought")
env_scales <- rbindlist(c(
  lapply(scale_vars, function(v) data.table(variable = v,
                                            mean = mean(env_unscaled[[v]]),
                                            sd   = sd(env_unscaled[[v]]))),
  list(data.table(variable = "management", mean = 0, sd = 1))))
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
  max(i1$obs$siteID) < min(i2$obs$siteID)                           # decades never collide
)

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
fwrite(reg_tab,            file.path(OUT_NAT, "reg_target.csv"))
fwrite(i1$qc,              file.path(OUT_NAT, "wzp_qc_i1.csv"))
fwrite(i2$qc,              file.path(OUT_NAT, "wzp_qc_i2.csv"))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Spatially blocked split: train / holdout_blocked / pilot
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## one location per tract (corner_env x/y are the tract's WGS84 lon/lat)
sids[, Tnr := as.integer(sub("_.*$", "", corner))]
tr <- unique(env_in[Tnr %in% sids$Tnr, .(Tnr, lon = x, lat = y)], by = "Tnr")
xy <- sf::st_coordinates(sf::st_transform(
  sf::st_as_sf(tr, coords = c("lon", "lat"), crs = 4326), 3035))
tr[, `:=`(east = xy[, 1], north = xy[, 2])]
bm <- PARAMS$block_km * 1000
tr[, block := sprintf("E%03d_N%03d", floor(east / bm), floor(north / bm))]

set.seed(PARAMS$seed)
blocks  <- sort(unique(tr$block))
hold_bl <- sample(blocks, round(PARAMS$holdout_share * length(blocks)))
tr[, set := fifelse(block %in% hold_bl, "holdout_blocked", "train")]

site_set <- merge(unique(sids[, .(siteID, interval, Tnr)]), tr[, .(Tnr, block, set)],
                  by = "Tnr")
stopifnot(nrow(site_set) == n_sites,
          site_set[, uniqueN(set), by = Tnr][, all(V1 == 1L)])   # tracts never split

## distance from every holdout site to its nearest training site (km)
hold_xy  <- as.matrix(tr[set == "holdout_blocked", .(east, north)])
train_xy <- as.matrix(tr[set == "train", .(east, north)])
d_km <- FNN::get.knnx(train_xy, hold_xy, k = 1)$nn.dist[, 1] / 1000
hold_dist <- c(median = median(d_km), min = min(d_km), p10 = unname(quantile(d_km, 0.10)))

## pilot and train_quarter: whole training tracts (both decades kept together),
## taken in ONE random tract order up to each target, so pilot sits inside
## train_quarter and both are fixed by the seed
train_tr <- site_set[set == "train", .N, by = Tnr]
train_tr <- train_tr[sample(.N)]
pilot_tr <- train_tr[cumsum(N) <= PARAMS$pilot_target, Tnr]
site_set_pilot <- site_set[Tnr %in% pilot_tr]
quarter_tr <- train_tr[cumsum(N) <= PARAMS$quarter_share * sum(N), Tnr]
site_set_quarter <- site_set[Tnr %in% quarter_tr]

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
  fwrite(species_dt,     file.path(d, "species_dt.csv"))
  fwrite(s,              file.path(d, "siteID_dt.csv"))
  data.table(set = name, sites = nrow(mp),
             interval_1 = s[, uniqueN(siteID[interval == 1L])],
             interval_2 = s[, uniqueN(siteID[interval == 2L])],
             tracts = uniqueN(sub("_.*$", "", s$corner)))
}
set_tab <- rbind(
  write_set("train",           site_set[set == "train", siteID]),
  write_set("holdout_blocked", site_set[set == "holdout_blocked", siteID]),
  write_set("train_quarter",   site_set_quarter$siteID),
  write_set("pilot",           site_set_pilot$siteID))
set_tab[, blocks := c(length(blocks) - length(hold_bl), length(hold_bl),
                      uniqueN(tr[Tnr %in% quarter_tr, block]),
                      uniqueN(tr[Tnr %in% pilot_tr, block]))]
fwrite(tr[, .(Tnr, lon, lat, east, north, block, set,
              train_quarter = Tnr %in% quarter_tr,
              pilot = Tnr %in% pilot_tr)], file.path(OUT_VER, "split_blocks.csv"))
sids[, Tnr := NULL]
print(set_tab)
message(sprintf("[split] holdout -> nearest training tract: median %.1f km, min %.1f km, p10 %.1f km",
                hold_dist[["median"]], hold_dist[["min"]], hold_dist[["p10"]]))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## PREPARATION.md
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
commit <- tryCatch(system2("git", c("describe", "--always", "--dirty"), stdout = TRUE),
                   error = function(e) "unknown")
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
reg_comp <- obs[, .(reg = sum(reg)), by = .(interval = fifelse(siteID <= off, 1L, 2L), species)]
reg_comp <- merge(reg_comp, species_dt, by = "species")
reg_comp[, share := reg / sum(reg), by = interval]
reg_comp <- dcast(reg_comp, species_name ~ interval, value.var = "share")
setnames(reg_comp, c("1", "2"), c("i1", "i2"))
setorder(reg_comp, -i1)

md <- c(
  sprintf("# Prepared data: %s", PARAMS$version),
  "",
  sprintf("Written by `scripts/70_prepare_two_interval.R` at commit `%s` on %s.",
          commit, format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "## What this is",
  "",
  "FINN inputs for **two decades on the same corners**: interval 1 = BWI2 -> BWI3",
  "(2002 -> 2012, ACCESS2003 release) and interval 2 = BWI3 -> BWI4 (2012 -> 2022,",
  "Thuenen inventory archive). Annual timestep, observation at year 10 of each",
  "decade, per-year drought anomaly in `env_dt$drought`. Each decade is its own",
  "FINN site-period; `siteID_dt.csv` maps `siteID` to `interval` and the physical",
  sprintf("`corner` (Tnr_Enr). Interval-2 siteIDs are offset by %d.", off),
  "",
  "## Contents (`national/`)",
  "",
  "| file | content |",
  "|---|---|",
  "| obs_dt.csv | siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length |",
  "| env_dt.csv | siteID, year 0..10, mat_c, map_mm, elev_m, x, y, management, drought -- z-scaled ONCE over both decades; management left on 0..3 |",
  "| env_unscaled_dt.csv / env_scales.csv | natural scale and the pooled mean/sd |",
  "| init_cohort_dt.csv | start-of-decade stand (BWI2 state for i1, BWI3 for i2), for makeInitCohorts() |",
  "| species_dt.csv | the 15 species, one code map for both decades |",
  "| siteID_dt.csv | siteID -> interval, corner, qualified siteName (`Tnr_Enr@i1/@i2`) |",
  "| reg_target.csv | regeneration target per survey x corner x species (diagnostic) |",
  "",
  "## Per decade",
  "",
  "| interval | sites | obs rows | mean ba | mean mort (/yr) | median growth (/yr) | mean reg (/ha/period) | drought coverage |",
  "|---|---|---|---|---|---|---|---|",
  per_iv(i1, 1L), per_iv(i2, 2L),
  "",
  sprintf("Management tiers (sites): interval 1 %s; interval 2 %s.", tier_tab(i1), tier_tab(i2)),
  "",
  "Interval-2 corner funnel:",
  "",
  "| step | corners |", "|---|---|",
  funnel[, sprintf("| %s | %d |", step, n)],
  "",
  "## Decisions",
  "",
  "- **Same code path.** Both decades run through one function that mirrors",
  "  `scripts/04` steps 1-7b (angle-count expansion, tallied sample sizes,",
  "  total-removal mortality, survivor-only growth, annualised rates).",
  "- **Interval-2 corners** pass the `scripts/02` filters (BA >= 5 m2/ha, 15 species",
  "  >= 99 % of BA, period 8-12 yr, no untraceable trees) and need environment in",
  "  `corner_env.csv`, so they are a subset of the interval-1 corners.",
  "- **Realised survey years** for interval 2 from `plot.acquisition_date` (earliest",
  "  per tract), matching `bwi_periods()` for interval 1; only period == 10 kept.",
  sprintf("- **Calamity salvage booked as natural death** in interval 2 (harvest_reason %s): %d of %d harvested stems (%.1f %%). Interval-2 tiers use regular harvest only. BWI2012 has no harvest_reason, so interval 1 keeps all harvest (incl. Kyrill 2007).",
          paste(PARAMS$calamity_reasons, collapse = ","), n_salvage, n_harvest_all,
          100 * n_salvage / n_harvest_all),
  "- **Mortality target = total removal** (natural + harvest) on the shared at-risk",
  "  denominator, as in `scripts/04`; management enters as an env predictor.",
  sprintf("- **Regeneration** from the archive circles, classes %s (dbh 5.0-6.9 cm), r = %.2f m, stock -> flux at %.2f cm/yr over a %.1f cm band: recruits per **hectare** per period (factor %.3f per stem/ha). Off-V share of individuals: %s.",
          paste(PARAMS$reg_classes, collapse = "+"), PARAMS$reg_circle_r_m,
          PARAMS$reg_incr_cm_yr, PARAMS$reg_band_cm, flux_factor,
          paste(off_V[, sprintf("%s %.1f %%", interval_name, 100 * share_off_V)], collapse = ", ")),
  sprintf("- Archive bwi2012 vs ACCESS2003 b3_jung (classes 5+6): %s.", reg_check),
  "- **Drought** = `z_soil_moist` from `scripts/16`, standardised per site over",
  "  2000-2022 (one baseline for both decades), placed at calendar year",
  "  `year_t0 + k`; uncovered site-years set to 0.",
  "- **Pooled scaling**: mean/sd computed once over both decades (`env_scales.csv`).",
  "",
  "## Regeneration target composition (share of recruits)",
  "",
  "| species | interval 1 | interval 2 |", "|---|---|---|",
  reg_comp[, sprintf("| %s | %.3f | %.3f |", species_name, i1, i2)],
  "",
  "## Sets (spatially blocked)",
  "",
  sprintf("Whole tracts, both decades together, in %d km blocks on EPSG:3035; %d of %d blocks (%.0f %%) held out at seed %d. Training sites come only from retained blocks. Every set folder has the national file layout with siteIDs renumbered 1..n (`siteID_dt$siteID_national` maps back); all use the national pooled scaling (`env_scales.csv`). `split_blocks.csv` lists tract -> block -> set.",
          PARAMS$block_km, length(hold_bl), length(blocks),
          100 * length(hold_bl) / length(blocks), PARAMS$seed),
  "",
  "| set | site-periods | interval 1 | interval 2 | tracts | blocks |", "|---|---|---|---|---|---|",
  set_tab[, sprintf("| %s | %d | %d | %d | %d | %d |", set, sites, interval_1,
                    interval_2, tracts, blocks)],
  "",
  sprintf("Distance from each holdout tract to the nearest training tract: median %.1f km, minimum %.1f km, 10th percentile %.1f km.",
          hold_dist[["median"]], hold_dist[["min"]], hold_dist[["p10"]]),
  "`train_quarter/` (~25 % of train, for faster fits) and `pilot/` (timing, quick iterations) are fixed random subsets of whole training tracts, drawn from one tract order so `pilot` lies inside `train_quarter`; they are not separate populations.",
  "",
  "## Checks",
  "",
  sprintf("- Interval 1 vs the validated single-decade inputs (`outputs/finn`): %s.", reg_cmp),
  "- Stitched tables: siteIDs contiguous and unique, decades never collide, env",
  "  covers years 0..10 for every site, init sites == obs sites, no NA in env.",
  "",
  "## Caveats",
  "",
  sprintf("- **Ingrowth asymmetry.** Interval 1's start stand includes laterally ingrown stems at their back-modelled 2002 diameter (`bwi_tree_panel()`), %.1f %% of its start basal area; interval 2's start stand is the plain 2012 tally. The interval-1 initial state is therefore denser by construction. Kept so interval 1 stays identical to the validated recipe; harmonising it is an open decision.",
          100 * i1$diag$lat_share),
  "- **Regeneration geometry for 2022** is assumed unchanged (r = 1.75 m for",
  "  classes 5+6). Counts on corners surveyed in both years fall by only ~6 %",
  "  (2012 -> 2022), consistent with the same circle; confirm in the BWI4 manual.",
  sprintf("- The regeneration **magnitude** rests on the assumed %.2f cm/yr increment in the 5-7 cm band; composition does not.",
          PARAMS$reg_incr_cm_yr),
  "- **Tier M1 is near-empty in both decades.** With BAF 4 one tallied tree carries",
  "  4 m2/ha, so any harvest removes > 5 % of BA unless the corner holds > 80 m2/ha;",
  "  the 0..3 code is effectively 0 / 2 / 3. Same definition in both decades.",
  "- Environment exists only for interval-1 corners; extending climate to all",
  "  BWI4 corners would enlarge interval 2.",
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
