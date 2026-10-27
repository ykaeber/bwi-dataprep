#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-bwi | scripts/02_select_calibration_set.R
##
## Turns the descriptive tables from 01 into a concrete calibration cohort by
## applying the parameters you choose below. Read outputs/management_tiers.csv,
## outputs/strictness_grid.csv and outputs/trakt_summary.csv first -- they show
## what each knob costs in sample size.
##
## Run from the project root, after 01:
##   Rscript scripts/02_select_calibration_set.R
##
## Outputs (outputs/):
##   calibration_set_corners.csv   the cohort, one row per corner
##   calibration_set_species.csv   the target species vector V
##   calibration_set_summary.txt   parameters used + resulting counts
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

#=#=# PARAMS -- edit freely #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
PARAMS <- list(

  ## ---- What counts as a FINN site and a FINN patch ------------------------
  ## "ecke"  : site = one tract corner, 1 patch. The corner is the unit that
  ##           was actually measured, so this keeps every usable corner --
  ##           ~23k sites at tier M0. Best for estimating environmental
  ##           responses, which is what the wide climate gradient is for.
  ## "trakt" : site = the tract, patches = its corners, requiring exactly
  ##           `n_patches` usable corners. Corners of a tract lie 150 m apart
  ##           and share an environment, so this mirrors FIA's 4 subplots and
  ##           averages down the angle-count noise (a single corner tallies
  ##           ~7 trees, i.e. ~38 % standard error on BA/ha). It costs a lot
  ##           of sample: only 733 tracts have 4 usable corners all at M0.
  patch_mode         = "ecke",
  n_patches          = 4L,          # used only when patch_mode == "trakt"

  ## ---- Management ---------------------------------------------------------
  ## FINN has no harvest process. Over 2002-2012 harvest removed basal area
  ## ~6x faster than natural mortality did (1.56 vs 0.26 %/yr), so admitting
  ## harvested corners means asking the mortality process to absorb
  ## silviculture. M0 = corners with no tree cut at all; it still retains
  ## about half of Germany's corners.
  ## FINN_BWI_MGMT_TIERS overrides this (comma-separated), e.g. "M0,M1,M2,M3" to
  ## include managed corners for the harvest-process calibration. Default M0
  ## reproduces the biological-only set.
  mgmt_tiers         = strsplit(Sys.getenv("FINN_BWI_MGMT_TIERS", unset = "M0"),
                                ",", fixed = TRUE)[[1]],

  ## ---- Species ------------------------------------------------------------
  ## V = the top `n_species` by corner occurrence. A corner is kept only if V
  ## accounts for at least `species_X` of its basal area at the period start;
  ## the residual off-V trees are then dropped in script 04. 15 species covers
  ## 95 % of national basal area and keeps ~40k corners at X = 0.99.
  n_species          = 15L,
  species_X          = 0.99,

  ## ---- Optional geographic / structural restrictions ----------------------
  bundesland         = NULL,        # e.g. c(8, 9) for BW + Bayern; NULL = all
  min_ba_ha_t0       = 5,           # drop near-empty corners (m^2/ha)
  max_ba_ha_t0       = Inf,
  min_period_length  = 8L,
  max_period_length  = 12L
)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

suppressPackageStartupMessages(library(data.table))
source("R/config.R")
source("R/bwi_io.R")
source("R/wzp.R")
source("R/bwi_trees.R")

stopifnot(PARAMS$patch_mode %in% c("ecke", "trakt"))

cm      <- fread(file.path(OUT_DIR, "corner_metadata.csv"))
sp_rank <- fread(file.path(OUT_DIR, "species_rank.csv"))

log_step <- function(label, n_before, n_after) {
  message(sprintf("[%-22s] %7d -> %7d corners", label, n_before, n_after))
  n_after
}

n <- nrow(cm)
message(sprintf("Starting from %d corners.", n))

## 1) usability (Thuenen's Use23 flag + measured at both dates + classifiable)
d <- cm[usable == TRUE]; n <- log_step("usable", n, nrow(d))

## 2) management tier
d <- d[mgmt_tier %in% PARAMS$mgmt_tiers]
n <- log_step(paste("tier", paste(PARAMS$mgmt_tiers, collapse = "+")), n, nrow(d))

## 3) region
if (!is.null(PARAMS$bundesland)) {
  d <- d[Bl %in% PARAMS$bundesland]
  n <- log_step("Bundesland", n, nrow(d))
}

## 4) structure + period length
d <- d[is.finite(ba_ha_t0) &
       ba_ha_t0 >= PARAMS$min_ba_ha_t0 & ba_ha_t0 <= PARAMS$max_ba_ha_t0]
n <- log_step("basal area window", n, nrow(d))

d <- d[period_length >= PARAMS$min_period_length &
       period_length <= PARAMS$max_period_length]
n <- log_step("period length", n, nrow(d))

## 5) species contamination filter
V <- sp_rank[order(rank)][seq_len(min(PARAMS$n_species, .N)), species_name]

panel <- bwi_tree_panel()$trees
panel[, ba_ha := n_ha * bhd2ba(dbh)]
ba_corner <- panel[time == 0L & living == TRUE,
                   .(ba_tot = sum(ba_ha, na.rm = TRUE),
                     ba_V   = sum(ba_ha[species_name %in% V], na.rm = TRUE)),
                   by = .(Tnr, Enr)]
ba_corner[, share_V := fifelse(ba_tot > 0, ba_V / ba_tot, 0)]

d <- merge(d, ba_corner[, .(Tnr, Enr, share_V)], by = c("Tnr", "Enr"))
d <- d[share_V >= PARAMS$species_X]
n <- log_step(sprintf("species V>=%.0f%%", 100 * PARAMS$species_X), n, nrow(d))

## 6) patch geometry
if (identical(PARAMS$patch_mode, "trakt")) {
  d[, n_corners_in_set := .N, by = Tnr]
  d <- d[n_corners_in_set == PARAMS$n_patches]
  n <- log_step(sprintf("tracts with %d corners", PARAMS$n_patches), n, nrow(d))
  d[, siteName  := Tnr]
  d[, patchName := Enr]
} else {
  d[, siteName  := paste(Tnr, Enr, sep = "_")]
  d[, patchName := 1L]
}

if (!nrow(d)) stop("No corners survive the current PARAMS.", call. = FALSE)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Write
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

fwrite(d, file.path(OUT_DIR, "calibration_set_corners.csv"))

sp_out <- sp_rank[species_name %in% V][order(rank)]
fwrite(sp_out, file.path(OUT_DIR, "calibration_set_species.csv"))

summary_txt <- c(
  "FINN-bwi calibration set",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "",
  "PARAMS",
  paste0("  ", names(unlist(PARAMS)), " = ", unlist(PARAMS)),
  "",
  "RESULT",
  sprintf("  corners        : %d", nrow(d)),
  sprintf("  sites          : %d", uniqueN(d$siteName)),
  sprintf("  patches / site : %s",
          if (identical(PARAMS$patch_mode, "trakt")) PARAMS$n_patches else 1L),
  sprintf("  tracts         : %d", uniqueN(d$Tnr)),
  sprintf("  Bundeslaender  : %d", uniqueN(d$Bl)),
  sprintf("  species in V   : %d", length(V)),
  sprintf("  BA/ha t0       : median %.1f  (5-95%%: %.1f - %.1f)",
          median(d$ba_ha_t0), quantile(d$ba_ha_t0, .05), quantile(d$ba_ha_t0, .95)),
  sprintf("  period length  : median %.0f yr (range %d - %d)",
          median(d$period_length), min(d$period_length), max(d$period_length)),
  sprintf("  harvest share  : max %.3f", max(d$harvest_share)),
  sprintf("  mort share     : mean %.4f", mean(d$mort_share)),
  "",
  "SPECIES VECTOR V",
  paste0("  ", sp_out$rank, ". ", sp_out$species_name,
         "  (", sp_out$n_corners_present, " corners, ",
         round(100 * sp_out$ba_share, 2), " % of national BA)")
)
writeLines(summary_txt, file.path(OUT_DIR, "calibration_set_summary.txt"))

message("\n", paste(summary_txt[8:20], collapse = "\n"))
message(sprintf("\nDone. Wrote calibration_set_{corners,species}.csv + summary to %s/",
                OUT_DIR))
