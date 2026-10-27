#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-bwi | scripts/04_prepare_finn_inputs.R
##
## Turns the calibration cohort into the four objects FINN's fit() consumes:
## an observation table, an environment table, an initial cohort and a species
## lookup.
##
## THREE BWI-SPECIFIC DECISIONS ARE MADE HERE
##
## 1. ANGLE-COUNT -> FIXED-AREA PATCHES (R/wzp.R)
##    Each tallied tree stands for n_ha = K/g stems per hectare. We expand it
##    into round(n_ha * patch_size) identical virtual trees, allocating by
##    largest remainder within (patch x species) so group totals survive the
##    rounding. On the real data this is unbiased (mean basal-area error
##    -0.01 %, median |error| 1.1 % at 0.1 ha) and far below the angle-count's
##    own ~38 % per-corner sampling error.
##
## 2. SAMPLE SIZES ARE NOT EXPANDED
##    Expansion is right for densities, wrong for information content: one
##    measured tree must not tell a binomial likelihood it is forty. So the
##    rates (mort, growth) are computed from the expanded, density-weighted
##    table, while `n_at_risk`, `n_died` and `growth_n` are overwritten with
##    counts of TALLIED trees. Harvested stems are censored out of the
##    at-risk denominator rather than counted as survivors.
##
## 3. AN ANNUAL TIME GRID
##    FINN advances one timestep per row of the time axis -- there is no
##    sub-stepping -- so if the axis were just {BWI2, BWI3} the fitted growth
##    and mortality rates would be per decade, not per year. We therefore lay
##    the observations on an annual grid: env covers years 0..T-1 and each
##    site's observation sits at its realised period length (8-12 years), with
##    every other year NA and masked out of the loss. Set time_mode = "period"
##    for the fast, coarse alternative.
##
## Run from the project root, after 03:
##   Rscript scripts/04_prepare_finn_inputs.R
##
## Outputs (outputs/finn/):
##   obs_dt.csv           site-aggregated observations -- the table fit() takes
##   obs_dt_patches.csv   patch-level, for diagnostics only
##   env_dt.csv           z-scaled covariates per (siteID, year)
##   env_unscaled_dt.csv  same on the natural scale
##   env_scales.csv       mean + sd used for the z-scaling
##   init_cohort_dt.csv   initial state, fed to makeInitCohorts() at fit time
##   species_dt.csv       species integer code <-> name
##   siteID_dt.csv        siteName <-> siteID index map
##   wzp_qc.csv           per-corner basal-area conservation of the expansion
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

suppressPackageStartupMessages({
  library(data.table)
  library(FINN)
})
source("R/config.R")
source("R/bwi_io.R")
source("R/wzp.R")
source("R/bwi_trees.R")

#=#=# PARAMS -- edit freely #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
PARAMS <- list(
  ## Nominal patch area for the angle-count expansion (ha). Must match the
  ## patch_size passed to fit(). See R/config.R for the rounding threshold.
  patch_size_ha   = PATCH_SIZE_HA,

  ## "annual" : env spans years 0..fix_period_length and the single observation
  ##            sits at the last year, carrying a `period_length` column. FINN
  ##            then simulates annually and averages its rates over the period
  ##            before comparing (finn-class.R, the `accumulate_gradients`
  ##            branch), so the fitted rates are PER YEAR.
  ## "period" : one timestep; rates come out per remeasurement period. Cheaper,
  ##            and directly comparable to the FINN-fia panel setup.
  time_mode       = "annual",

  ## FINN reads the aggregation window as `period = unique(period_length)` and
  ## assumes it is constant across sites ("currently we assume that they are
  ## constant over sites", finn-class.R). BWI intervals realised 8-12 years, so
  ## we keep a single value and drop the rest rather than let a mixed window
  ## silently mis-align the rate aggregation. 10 years is the modal interval.
  ## NULL keeps every site and forces time_mode = "period".
  fix_period_length = 10L,

  ## Optional subsample of sites, for a quick pilot. NULL = keep all.
  ## Override from the shell without editing: FINN_BWI_MAX_SITES=2000
  max_sites       = {
    v <- Sys.getenv("FINN_BWI_MAX_SITES", unset = "")
    if (nzchar(v)) as.integer(v) else NULL
  },
  seed            = 42L,

  ## Environmental covariates carried into env_dt.
  env_vars        = c("mat_c", "map_mm", "elev_m"),

  ## Per-YEAR drought forcing from scripts/16_drought_index.R. The calibration
  ## period contains the 2003 European heatwave (national growing-season soil
  ## moisture at -2.2 sd), which a static 1970-2000 climate normal cannot
  ## represent at all. Because scripts/04 lays the observations on an annual
  ## grid, the drought anomaly can vary year by year within the period -- the
  ## covariate is genuinely time-varying rather than another site constant.
  ## Set to NULL to reproduce the climate-normals-only calibration.
  drought_var     = Sys.getenv("FINN_BWI_DROUGHT", unset = "z_soil_moist"),
  ## Spatial coordinates ride along z-scaled so a hybrid process can consume
  ## them and ALE can perturb them; the climate-only fit ignores them.
  include_xy      = TRUE,
  xy_vars         = c("x", "y"),

  ## Drop a site if any tree's diameter change over the period falls outside
  ## this window (cm) -- guards against transcription errors.
  dbh_growth_thresh = c(-10, 60)
)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

stopifnot(PARAMS$time_mode %in% c("annual", "period"))
set.seed(PARAMS$seed)

OUT_FINN <- file.path(OUT_DIR, "finn")
dir.create(OUT_FINN, recursive = TRUE, showWarnings = FALSE)

add_site_completeness_flags <- FINN:::add_site_completeness_flags
ENV_COLS <- c(PARAMS$env_vars,
              if (isTRUE(PARAMS$include_xy)) PARAMS$xy_vars else character(0))
USE_DROUGHT <- nzchar(PARAMS$drought_var %||% "") &&
  file.exists(file.path(OUT_DIR, "site_drought_annual.csv"))
if (nzchar(PARAMS$drought_var %||% "") && !USE_DROUGHT)
  warning("drought_var set but outputs/site_drought_annual.csv is missing -- ",
          "run scripts/16_drought_index.R. Continuing without it.",
          call. = FALSE, immediate. = TRUE)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 1) Cohort, environment, species vector
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

## A pre-built panel + corners (e.g. the new inventory-archive schema, 2012->2022,
## from scripts/40) can be injected via FINN_BWI_CORNERS / FINN_BWI_PANEL, reusing
## everything below unchanged. Default is the ACCESS2003 BWI2->BWI3 pipeline.
CORNERS_SRC <- Sys.getenv("FINN_BWI_CORNERS",
                          unset = file.path(OUT_DIR, "calibration_set_corners.csv"))
PANEL_SRC   <- Sys.getenv("FINN_BWI_PANEL", unset = "")
message("Corners source: ", CORNERS_SRC,
        if (nzchar(PANEL_SRC)) paste0(" | panel: ", PANEL_SRC) else "")

corners <- fread(CORNERS_SRC)
env_in  <- fread(file.path(OUT_DIR, "corner_env.csv"))
V       <- fread(file.path(OUT_DIR, "calibration_set_species.csv"))$species_name

## Management intensity as an ENVIRONMENT predictor (opt-in FINN_BWI_MANAGEMENT=1).
## It is an exogenous per-corner forcing, like drought: the mortality process
## learns how removal responds to management, and succession sets it (0 =
## unmanaged). Derived from the corner's management tier; NOT a separate response.
USE_MGMT <- nzchar(Sys.getenv("FINN_BWI_MANAGEMENT", unset = ""))
keep_cols <- c("Tnr", "Enr", "siteName", "patchName", "period_length")
if (USE_MGMT) {
  if (!"mgmt_tier" %in% names(corners))
    stop("FINN_BWI_MANAGEMENT set but corners have no mgmt_tier column.", call. = FALSE)
  tier_intensity <- c(M0 = 0, M1 = 1, M2 = 2, M3 = 3)
  corners[, management := as.numeric(tier_intensity[mgmt_tier])]
  keep_cols <- c(keep_cols, "management")
  message("Management intensity attached as an env predictor (M0..M3 -> 0..3).")
}

## Keep only corners that have both a cohort row and complete environment.
corners <- merge(corners[, ..keep_cols],
                 env_in[, c("Tnr", "Enr", ENV_COLS), with = FALSE],
                 by = c("Tnr", "Enr"))
if (USE_MGMT) ENV_COLS <- c(ENV_COLS, "management")
message(sprintf("Corners with environment: %d on %d sites",
                nrow(corners), uniqueN(corners$siteName)))

if (!is.null(PARAMS$max_sites)) {
  keep <- sample(unique(corners$siteName),
                 min(PARAMS$max_sites, uniqueN(corners$siteName)))
  corners <- corners[siteName %in% keep]
  message(sprintf("  subsampled to %d sites", uniqueN(corners$siteName)))
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 2) Tallied-tree panel, restricted to the cohort and to V
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

panel <- if (nzchar(PANEL_SRC)) fread(PANEL_SRC) else bwi_tree_panel()$trees
## bwi_tree_panel() names sites by tract and patches by corner; script 02 has
## already resolved siteName/patchName for the chosen patch_mode, so drop the
## panel's defaults rather than letting the merge produce .x/.y duplicates.
panel[, c("siteName", "patchName") := NULL]
panel <- merge(panel, corners[, .(Tnr, Enr, siteName, patchName, period_length)],
               by = c("Tnr", "Enr"))

n0 <- nrow(panel)
panel <- panel[species_name %in% V]
message(sprintf("[species] restricted to V (%d species): %d -> %d records",
                length(V), n0, nrow(panel)))

## Diameter-change screen, applied at the level of the measured tree.
wide <- dcast(panel[, .(treeKey, time, dbh)], treeKey ~ time, value.var = "dbh")
setnames(wide, c("treeKey", "dbh_t0", "dbh_t1"), skip_absent = TRUE)
if (all(c("dbh_t0", "dbh_t1") %in% names(wide))) {
  bad <- wide[is.finite(dbh_t0) & is.finite(dbh_t1) &
              (dbh_t1 - dbh_t0 < PARAMS$dbh_growth_thresh[1] |
               dbh_t1 - dbh_t0 > PARAMS$dbh_growth_thresh[2]), treeKey]
  if (length(bad)) {
    bad_sites <- unique(panel[treeKey %in% bad, siteName])
    message(sprintf("[dbh screen] %d implausible trees -> dropping %d sites",
                    length(bad), length(bad_sites)))
    panel <- panel[!siteName %in% bad_sites]
  }
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 3) Angle-count expansion -> virtual trees on a fixed-area patch
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

## Allocate integer counts on the START-of-period state and carry that count
## through to t1, so a tree does not change how many stems it represents
## halfway through the period. (n_ha does drift as the stem grows -- that drift
## is a property of the sampling design, not of the forest.)
panel[, patchKey := paste(siteName, patchName, sep = "/")]

alloc <- panel[time == 0L, .(treeKey, patchKey, species_name, n_ha, dbh)]
invisible(wzp_allocate(alloc, by = c("patchKey", "species_name"),
                       patch_size_ha = PARAMS$patch_size_ha))

qc <- wzp_ba_conservation(alloc, by = "patchKey",
                          patch_size_ha = PARAMS$patch_size_ha)
fwrite(qc, file.path(OUT_FINN, "wzp_qc.csv"))
message(sprintf(
  "[wzp] basal-area conservation: median |err| %.2f %%, mean err %+.3f %% over %d patches",
  100 * median(abs(qc$rel_err), na.rm = TRUE),
  100 * mean(qc$rel_err, na.rm = TRUE), nrow(qc)))

## Recruits appear only at t1 and so get no allocation above; give them one.
rec <- panel[time == 1L & fate == "recruit",
             .(treeKey, patchKey, species_name, n_ha, dbh)]
if (nrow(rec)) {
  invisible(wzp_allocate(rec, by = c("patchKey", "species_name"),
                         patch_size_ha = PARAMS$patch_size_ha))
  alloc <- rbind(alloc, rec)
}

panel <- merge(panel, alloc[, .(treeKey, n_trees)], by = "treeKey", all.x = TRUE)
panel[is.na(n_trees), n_trees := 0L]

## index_by includes `time`: the panel holds one row per (tree, date), and the
## replicate index has to restart at each date so that virtual tree k at t0 is
## the same tree as virtual tree k at t1.
virtual <- wzp_expand(panel,
                      keep = c("treeKey", "siteName", "patchName",
                               "species_name", "dbh", "time", "status",
                               "living", "mort", "fate", "period_length"),
                      index_by = c("treeKey", "time"))
message(sprintf("[wzp] %d tallied records -> %d virtual trees (%.1fx)",
                nrow(panel), nrow(virtual), nrow(virtual) / nrow(panel)))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 4) tree_dt in FINN's schema, then makeObsData()
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

## Harvested stems leave the stand but are not mortality. They stay in the t0
## state (they were really there) and are dropped from t1, which is what
## `living == FALSE` already does; the at-risk correction in step 6 then
## censors them instead of scoring them as survivors.
tree_dt <- virtual[, .(
  siteName, patchName, treeName, year = time,
  species_name, dbh, status, living
)]

## makeObsData() reads `year` as a real time coordinate for period_length, so
## use the realised interval rather than the 0/1 index.
tree_dt <- merge(tree_dt, unique(virtual[, .(siteName, period_length)]),
                 by = "siteName")
tree_dt[, year := fifelse(year == 0L, 0L, as.integer(period_length))]
tree_dt[, period_length := NULL]

env_dt0 <- unique(corners[, c("siteName", ENV_COLS), with = FALSE])
site_years <- unique(tree_dt[, .(siteName, year)])
env_dt0 <- merge(site_years, env_dt0, by = "siteName")

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
message(sprintf("[makeObsData] %d rows on %d sites",
                nrow(obs_raw), uniqueN(obs_raw$siteName)))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 5) Integer site / patch / species indices
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

inputs <- resolveSiteIDs(
  tree_dt           = as.data.table(tree_obs),
  obs_dt            = as.data.table(obs_raw),
  env_dt            = as.data.table(env_dt0),
  createInitCohorts = FALSE
)

obs_dt     <- inputs$obs_dt                 # site-aggregated
obs_patch  <- inputs$obs_dt_patches
env_dt     <- inputs$env_dt
siteID_dt  <- inputs$siteID_dt
species_dt <- inputs$species_dt
init_trees <- inputs$tree_dt[year == 0]
init_trees[living == FALSE | is.na(living), trees := 0]
if (!"trees" %in% names(init_trees)) init_trees[, trees := 1L]

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 6) Replace expanded sample sizes with tallied ones
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

## One row per MEASURED tree that was standing at the period start.
##
## Only genuine BWI2 sample trees may enter these cohorts (see R/wzp.R):
##   - harvest          censored -- a competing risk FINN does not model
##   - ingrowth_lateral excluded -- classed as ingrowth *because* it was tallied
##                      at t1, so it is conditioned on survival (immortal time
##                      bias), and its t0 diameter is a Sloboda back-model
##                      rather than a measurement, so its increment is not an
##                      observation either
## That leaves survivors and natural deaths: the trees actually measured twice.
dbh_w <- dcast(panel[, .(treeKey, time, dbh)], treeKey ~ time, value.var = "dbh")
setnames(dbh_w, c("0", "1"), c("dbh_t0", "dbh_t1"), skip_absent = TRUE)

tal <- merge(panel[time == 0L & n_trees > 0L,
                   .(treeKey, siteName, patchName, species_name, n_ha, fate)],
             dbh_w, by = "treeKey", all.x = TRUE)
## Harvest is now a modelled competing risk, not a censored one: every stem
## standing at t0 is at risk of BOTH natural death and harvest, and both events
## share that denominator. For M0 corners no stem is harvested, so this is
## identical to the old `survivor|dead` at-risk set; for managed corners it adds
## the harvested stems to the population and records the harvest event.
tal[, `:=`(
  at_risk    = fate %in% c("survivor", "dead", "harvest"),
  died       = fate == "dead",
  harvested  = fate == "harvest",
  growth_ok  = fate == "survivor",
  rel_growth = fifelse(is.finite(dbh_t0) & is.finite(dbh_t1) & dbh_t0 > 0,
                       dbh_t1 / dbh_t0 - 1, NA_real_)
)]
message(sprintf(
  "[cohorts] at risk %d of %d measured stems (%d natural deaths, %d harvested, excl. %d lateral ingrowth); growth on %d survivors",
  sum(tal$at_risk), nrow(tal), sum(tal$died), sum(tal$harvested),
  sum(tal$fate == "ingrowth_lateral"), sum(tal$growth_ok)))

eff <- wzp_effective_n(tal, by = c("siteName", "patchName", "species_name"))

## Aggregate to the site, matching obs_dt's own aggregation: counts sum, the
## rate is the pooled design-weighted rate.
eff_site <- eff[, .(
  n_at_risk = sum(n_at_risk, na.rm = TRUE),
  n_died    = sum(n_died_eff, na.rm = TRUE),
  n_harvest = sum(n_harvest_eff, na.rm = TRUE),
  growth_n  = sum(growth_n,  na.rm = TRUE),
  ## growth pooled over patches the way the loss will weight it
  growth_obs = if (sum(growth_n, na.rm = TRUE) > 0)
                 sum(growth_w * growth_n, na.rm = TRUE) /
                 sum(growth_n[is.finite(growth_w)], na.rm = TRUE) else NA_real_
), by = .(siteName, species_name)]

## resolveSiteIDs()'s aggregated obs_dt already carries `species_name`; what it
## lacks is the original site label, which is the key `eff_site` is built on.
obs_dt <- merge(obs_dt, unique(siteID_dt[, .(siteID, siteName)]),
                by = "siteID", all.x = TRUE)
obs_dt <- merge(obs_dt, eff_site, by = c("siteName", "species_name"),
                all.x = TRUE, suffixes = c("_expanded", ""))

for (v in c("n_at_risk", "n_died", "n_harvest", "growth_n"))
  obs_dt[is.na(get(v)), (v) := 0L]
## Only the observation year carries a transition; year 0 is the initial state.
obs_dt[year == 0L, `:=`(n_at_risk = 0L, n_died = 0L, n_harvest = 0L, growth_n = 0L,
                        growth_obs = NA_real_)]
## Mortality target = TOTAL removal (natural death + harvest) on the shared
## at-risk denominator -- NOT a separate harvest response. For M0 corners harvest
## is 0 so this is the natural-mortality rate unchanged; with managed corners the
## mortality process models the full removal, and MANAGEMENT INTENSITY enters as
## an environment predictor (FINN_BWI_MANAGEMENT) so succession can toggle it
## (management = 0 recovers the unmanaged/biological trajectory).
obs_dt[, mort := fifelse(n_at_risk > 0, (n_died + n_harvest) / n_at_risk, NA_real_)]
obs_dt[, n_removed := n_died + n_harvest]   # kept for diagnostics

## Replace makeObsData()'s growth with the survivor-only, design-weighted one.
## makeObsData() averages over every living virtual tree, which silently folds
## in the back-modelled increments of laterally ingrown stems.
obs_dt[, growth := fifelse(growth_n > 0 & is.finite(growth_obs),
                           growth_obs, NA_real_)]
obs_dt[, growth_obs := NULL]

drop_cols <- grep("_expanded$", names(obs_dt), value = TRUE)
if (length(drop_cols)) obs_dt[, (drop_cols) := NULL]

message(sprintf("[sample sizes] median n_at_risk %.0f (expanded would have been %.0f)",
                median(obs_dt[n_at_risk > 0]$n_at_risk),
                median(eff$n_at_risk[eff$n_at_risk > 0]) *
                  mean(panel$n_trees[panel$n_trees > 0])))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 7) Time grid
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

## resolveSiteIDs() recodes years to a 0-based rank, so at this point the
## observation sits at year 1 whatever the real interval was. Put it back on a
## calendar footing and pad the intervening years with NA so the simulation
## runs annually.
site_period <- unique(corners[, .(siteName, period_length)])
obs_dt <- merge(obs_dt, site_period, by = "siteName", all.x = TRUE)
obs_dt[is.na(period_length), period_length := BWI_PERIOD_YEARS]

if (identical(PARAMS$time_mode, "annual")) {
  if (is.null(PARAMS$fix_period_length))
    stop("time_mode = 'annual' requires a fixed period length; ",
         "set PARAMS$fix_period_length or use time_mode = 'period'.",
         call. = FALSE)
  P <- as.integer(PARAMS$fix_period_length)

  n_before <- uniqueN(obs_dt$siteID)
  obs_dt <- obs_dt[period_length == P]
  if (!nrow(obs_dt))
    stop("No sites have a realised period length of ", P, " years.", call. = FALSE)
  message(sprintf("[time] fixed period length %d yr: %d -> %d sites",
                  P, n_before, uniqueN(obs_dt$siteID)))

  ## The single observation sits at the last simulated year. Everything before
  ## it is simulated but unobserved; FINN averages its rates over the window
  ## (i-P+1):i, so no NA-padded observation years are needed -- and must not be
  ## supplied, because a year whose observations are entirely NA yields a loss
  ## with no gradient and backward() then fails.
  obs_dt <- obs_dt[year == 1L]
  obs_dt[, year := P]

  ## ---- annualise the rate responses ------------------------------------
  ## fit() scores growth and mortality against the MEAN of the simulated
  ## annual rates over the window (i-P+1):i, and regeneration against their
  ## SUM (finn-class.R). makeObsData() however reports growth as the relative
  ## diameter increment over the WHOLE period and mortality as the proportion
  ## dying over the whole period, so on an annual grid both would be compared
  ## against a per-year quantity and the fitted rates would come out inflated
  ## by roughly a factor P. Convert them to per-year equivalents.
  ## Regeneration is a period total and already matches the summed prediction,
  ## so it is deliberately left alone.
  ##
  ## GROWTH SCALE. Annualising growth (^(1/P)) is correct for the default
  ## mean-of-annual-rates scoring, but it collapses growth's variance (~P^2),
  ## which inflates its "auto" loss weight and FLATTENS the fitted growth-env
  ## response (ALE range ~0.007 vs mortality ~0.7). FINN_BWI_GROWTH_SCALE=period
  ## keeps growth on the PERIOD scale (relative increment over the whole window)
  ## and is scored against FINN's compounded period growth (growth_period_scale
  ## in fit()) -- env stays annual so the drought forcing is unaffected. Mortality
  ## does not show the flattening, so it is annualised regardless.
  GROWTH_SCALE <- Sys.getenv("FINN_BWI_GROWTH_SCALE", unset = "annual")
  if (!identical(GROWTH_SCALE, "period"))
    obs_dt[is.finite(growth) & growth > -1,
           growth := (1 + growth)^(1 / P) - 1]
  else
    message("[time] GROWTH_SCALE=period: growth left as the period relative ",
            "increment (fit with growth_period_scale=TRUE)")
  obs_dt[is.finite(mort),
         mort := 1 - (1 - pmin(pmax(mort, 0), 1))^(1 / P)]
  obs_dt[is.finite(mort) & n_at_risk > 0,
         n_died := as.integer(round(mort * n_at_risk))]
  message(sprintf(
    "[time] annualised growth (median %.4f /yr) and mortality (mean %.5f /yr); reg left as a period total",
    median(obs_dt$growth, na.rm = TRUE), mean(obs_dt$mort, na.rm = TRUE)))

  env_dt <- env_dt[year == 0L]
  env_dt <- rbindlist(lapply(0:P, function(k) {
    e <- copy(env_dt); e[, year := k]; e
  }))

  ## Attach the per-year drought anomaly. Model year k is calendar year
  ## year_t0 + k for that site, so the 2003 anomaly lands on the right timestep
  ## for a tract surveyed in 2001 and for one surveyed in 2002 alike.
  if (USE_DROUGHT) {
    dr <- fread(file.path(OUT_DIR, "site_drought_annual.csv"))
    if (!PARAMS$drought_var %in% names(dr))
      stop("Column '", PARAMS$drought_var, "' not in site_drought_annual.csv",
           call. = FALSE)
    dr <- dr[, .(siteName, cal_year = year, drought = get(PARAMS$drought_var))]
    ## `corners` was subset to the columns needed upstream, so the survey year
    ## usually has to come back from the metadata table.
    key <- if ("year_t0" %in% names(corners)) {
      unique(corners[, .(siteName, year_t0)])
    } else {
      merge(unique(corners[, .(Tnr, Enr, siteName)]),
            fread(file.path(OUT_DIR, "corner_metadata.csv"))[
              , .(Tnr, Enr, year_t0)],
            by = c("Tnr", "Enr"))[, .(siteName, year_t0)]
    }
    grid <- CJ(siteName = unique(key$siteName), year = 0:P, unique = TRUE)
    grid <- merge(grid, key, by = "siteName")
    grid[, cal_year := year_t0 + year]
    grid <- merge(grid, dr, by = c("siteName", "cal_year"), all.x = TRUE)
    ## A site whose calendar window falls outside the DWD extraction gets 0,
    ## i.e. "average year" on the standardised scale, rather than NA.
    grid[!is.finite(drought), drought := 0]

    ## resolveSiteIDs() already carries siteName on env_dt; joining it again
    ## would produce siteName.x / siteName.y and break the key below.
    if (!"siteName" %in% names(env_dt))
      env_dt <- merge(env_dt, unique(siteID_dt[, .(siteID, siteName)]),
                      by = "siteID", all.x = TRUE)
    env_dt <- merge(env_dt, grid[, .(siteName, year, drought)],
                    by = c("siteName", "year"), all.x = TRUE)
    env_dt[!is.finite(drought), drought := 0]
    ENV_COLS <- c(ENV_COLS, "drought")
    message(sprintf(
      "[drought] attached '%s' as a per-year covariate: range %.2f to %.2f, driest model-year mean %.2f",
      PARAMS$drought_var, min(env_dt$drought), max(env_dt$drought),
      min(env_dt[, mean(drought), by = year]$V1)))
  }
  message(sprintf("[time] env spans years 0..%d (%d timesteps); rates are per year",
                  P, P + 1L))
} else {
  obs_dt <- obs_dt[year > 0L]
  obs_dt[, period_length := NA_real_]
  env_dt <- env_dt[year %in% c(0L, 1L)]
  message("[time] period grid: 1 timestep; rates are per remeasurement period")
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 7b) Align every table on the surviving sites, with contiguous siteIDs
##
## Dropping sites above leaves gaps in the siteID sequence and leaves env and
## init_cohort carrying sites the observations no longer mention. FINN's
## obsDF2arrays() would silently renumber them itself -- independently per
## table, which would scramble the site alignment between the observations,
## the environment and the initial state.
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

keep_sites <- sort(unique(obs_dt$siteID))
env_dt     <- env_dt[siteID %in% keep_sites]
init_trees <- init_trees[siteID %in% keep_sites]
siteID_dt  <- siteID_dt[siteID %in% keep_sites]

site_map <- data.table(siteID = keep_sites, siteID_new = seq_along(keep_sites))
for (dt in list(obs_dt, env_dt, init_trees, siteID_dt))
  dt[site_map, on = "siteID", siteID := i.siteID_new]

stopifnot(setequal(unique(env_dt$siteID), unique(obs_dt$siteID)),
          setequal(unique(init_trees$siteID), unique(obs_dt$siteID)),
          identical(sort(unique(obs_dt$siteID)), seq_len(length(keep_sites))))
message(sprintf("[align] %d sites; env rows %d, init trees %d",
                length(keep_sites), nrow(env_dt), nrow(init_trees)))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 8) z-scale the covariates
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

env_unscaled <- copy(env_dt)
env_scales   <- data.table()
for (v in ENV_COLS) {
  if (v == "management") {
    ## management is a controllable ORDINAL forcing (0 = unmanaged .. 3 = most
    ## intensively managed), not a measured gradient. Keep it on its natural
    ## 0..3 scale so that setting management = 0 in simulation genuinely means
    ## "unmanaged" (rather than a z-score). Record an identity scaling so ALE and
    ## the unscaling round-trip leave it unchanged.
    env_scales <- rbind(env_scales, data.table(variable = v, mean = 0, sd = 1))
    next
  }
  mu <- mean(env_dt[[v]], na.rm = TRUE); sdv <- sd(env_dt[[v]], na.rm = TRUE)
  env_scales <- rbind(env_scales, data.table(variable = v, mean = mu, sd = sdv))
  set(env_dt, j = v, value = (env_dt[[v]] - mu) / sdv)
}
keep <- c("siteID", "year", ENV_COLS)
env_dt       <- env_dt[, ..keep]
env_unscaled <- env_unscaled[, ..keep]

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 9) Write
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

## `period_length` must travel with the observations: FINN reads it to decide
## whether to compare against a single simulated year or against the mean of
## the preceding `period` years.
obs_keep <- intersect(c("siteID", "year", "species", "ba", "dbh", "trees",
                        "growth", "growth_n", "mort", "n_at_risk", "n_died",
                        "reg", "period_length"), names(obs_dt))
setorder(obs_dt, siteID, year, species)

fwrite(obs_dt[, ..obs_keep], file.path(OUT_FINN, "obs_dt.csv"))
fwrite(obs_patch,            file.path(OUT_FINN, "obs_dt_patches.csv"))
fwrite(env_dt,               file.path(OUT_FINN, "env_dt.csv"))
fwrite(env_unscaled,         file.path(OUT_FINN, "env_unscaled_dt.csv"))
fwrite(env_scales,           file.path(OUT_FINN, "env_scales.csv"))
fwrite(init_trees,           file.path(OUT_FINN, "init_cohort_dt.csv"))
fwrite(species_dt,           file.path(OUT_FINN, "species_dt.csv"))
fwrite(siteID_dt,            file.path(OUT_FINN, "siteID_dt.csv"))

message(sprintf("\nDone -> %s/", OUT_FINN))
message(sprintf("  sites     : %d", uniqueN(obs_dt$siteID)))
message(sprintf("  species   : %d", uniqueN(species_dt$species)))
message(sprintf("  timesteps : %d", uniqueN(obs_dt$year)))
message(sprintf("  obs rows  : %d  (%d with data)",
                nrow(obs_dt), obs_dt[!is.na(ba), .N]))
message(sprintf("  init trees: %d", nrow(init_trees)))
message("\nNext: Rscript scripts/05_fit_finn.R")
