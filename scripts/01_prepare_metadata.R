#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-bwi | scripts/01_prepare_metadata.R
##
## Descriptive tables for the BWI2 -> BWI3 corner population: what was
## measured, how much of it was harvested, which species carry the basal area,
## and how many corners survive each combination of filters.
##
## Nothing here is prescriptive. `corner_metadata.csv` is the selection table
## and cohort definitions are queries against it -- the same paradigm as
## FINN-fia. scripts/02 applies the choice you make after reading these.
##
## THE MANAGEMENT PROBLEM
## Germany's forest is almost entirely managed: over the 2002-2012 period
## harvest removes far more stems than natural mortality does. FINN has no
## harvest process, so on a corner where trees were cut the simulated stand
## cannot track the observed one, and any mortality signal is contaminated by
## silviculture. We therefore grade corners by removal intensity (M0-M3) and
## let the calibration cohort be chosen against a stated sample-size cost.
##
## Run from the project root, after 00:
##   Rscript scripts/01_prepare_metadata.R
##
## Outputs (outputs/):
##   corner_metadata.csv    one row per (Tnr, Enr) -- the selection table
##   trakt_summary.csv      one row per Tnr -- usable corners per tract
##   species_rank.csv       species by occurrence and basal-area share
##   management_tiers.csv   tier definitions + corner counts
##   strictness_grid.csv    corners retained under (tier x species coverage)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

suppressPackageStartupMessages(library(data.table))
source("R/config.R")
source("R/bwi_io.R")
source("R/wzp.R")
source("R/bwi_trees.R")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 1) Tree panel
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

message("Assembling the BWI2 -> BWI3 tree panel ...")
panel  <- bwi_tree_panel()
trees  <- panel$trees
flags  <- panel$corner_flags

message(sprintf("  %d tree-records, %d corners, %d tracts",
                nrow(trees), uniqueN(trees[, .(Tnr, Enr)]), uniqueN(trees$Tnr)))

## Represented basal area per hectare, by corner and date. Each tallied tree
## contributes n_ha * g, which for an uncorrected angle-count tree is just the
## Zaehlfaktor -- so this is the design-based BA/ha estimate.
trees[, ba_ha := n_ha * bhd2ba(dbh)]

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 2) Corner attributes from the inventory tables
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

b0 <- bwi("b0_ecke", select = c("Tnr", "Enr", "Bl", "Use23", "Gneig", "Gexp",
                                "natHoe", "potNatWg", "InvE2", "InvE3"))
e3 <- bwi("b3_ecke_w", select = c("Tnr", "Enr", "Wa", "Begehbar", "Eg", "Ne"))
setnames(e3, c("Wa", "Begehbar", "Eg", "Ne"),
         c("Wa_t1", "Begehbar_t1", "Eg_t1", "Ne_t1"))
e2 <- bwi("b2_ecke_w", select = c("Tnr", "Enr", "Wa", "Begehbar"))
setnames(e2, c("Wa", "Begehbar"), c("Wa_t0", "Begehbar_t0"))

## Elevation is measured per tract (b0_tab.HoeheNN), which beats a DEM lookup
## against the 1 km anonymised coordinate.
tab <- bwi("b0_tab", select = c("Tnr", "HoeheNN"))
tab <- tab[, .(elev_m = suppressWarnings(as.numeric(HoeheNN)[1])), by = Tnr]

per <- bwi_periods()

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 3) Per-corner demography
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

## State at each date. `living` excludes harvested and dead stems, so the t1
## figures describe the stand that was actually left standing.
state <- trees[living == TRUE, .(
  ba_ha    = sum(ba_ha,  na.rm = TRUE),
  n_ha     = sum(n_ha,   na.rm = TRUE),
  n_tallied = .N,
  n_species = uniqueN(species_name),
  dbh_mean = weighted.mean(dbh, n_ha, na.rm = TRUE)
), by = .(Tnr, Enr, time)]

state_w <- dcast(state, Tnr + Enr ~ time,
                 value.var = c("ba_ha", "n_ha", "n_tallied", "n_species", "dbh_mean"))
setnames(state_w, sub("_0$", "_t0", names(state_w)))
setnames(state_w, sub("_1$", "_t1", names(state_w)))

## Fluxes over the period, expressed as the represented BA/ha that left the
## stand by each route. Harvest and natural mortality are kept apart: FINN
## models the second and has no process for the first.
flux <- trees[time == 0L, .(
  ba_ha_harvest = sum(ba_ha[fate == "harvest"], na.rm = TRUE),
  ba_ha_dead    = sum(ba_ha[fate == "dead"],    na.rm = TRUE),
  n_ha_harvest  = sum(n_ha[fate == "harvest"],  na.rm = TRUE),
  n_ha_dead     = sum(n_ha[fate == "dead"],     na.rm = TRUE),
  n_tal_harvest = sum(fate == "harvest"),
  n_tal_dead    = sum(fate == "dead"),
  n_tal_atrisk  = .N
), by = .(Tnr, Enr)]

recr <- trees[time == 1L, .(
  n_ha_recruit  = sum(n_ha[fate == "recruit"], na.rm = TRUE),
  n_tal_recruit = sum(fate == "recruit"),
  n_ha_ingr_lat = sum(n_ha[fate == "ingrowth_lateral"], na.rm = TRUE)
), by = .(Tnr, Enr)]

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 4) Assemble corner_metadata
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

cm <- Reduce(function(x, y) merge(x, y, by = c("Tnr", "Enr"), all.x = TRUE),
             list(state_w, flux, recr, flags, b0, e3, e2))
cm <- merge(cm, tab, by = "Tnr", all.x = TRUE)
cm <- merge(cm, per[, .(Tnr, year_t0, year_t1, period_length)],
            by = "Tnr", all.x = TRUE)
cm <- merge(cm, BL_LOOKUP, by = "Bl", all.x = TRUE)

for (v in c("ba_ha_harvest", "ba_ha_dead", "n_ha_recruit", "n_tal_recruit",
            "n_ha_harvest", "n_ha_dead", "n_tal_harvest", "n_tal_dead"))
  cm[is.na(get(v)), (v) := 0]

## Removal intensity as a share of the standing basal area at the period start.
## Every tree present at t0 carries status "alive", including those later cut
## or found dead, so `ba_ha_t0` IS the full standing stock at the period start
## and must not have the fluxes added back on top of it.
cm[, harvest_share := fifelse(ba_ha_t0 > 0, ba_ha_harvest / ba_ha_t0, 0)]
cm[, mort_share    := fifelse(ba_ha_t0 > 0, ba_ha_dead    / ba_ha_t0, 0)]

## Usability. Use23 == 100 is Thuenen's own verdict that the corner may enter a
## BWI2 -> BWI3 change analysis; it already folds in "same place, forest and
## accessible at both dates". We add: measured at both dates, and no tree we
## could not classify.
cm[, usable := Use23 == USE23_OK &
     !is.na(ba_ha_t0) & !is.na(ba_ha_t1) &
     n_tallied_t0 > 0 & n_tallied_t1 > 0 &
     has_invalid == FALSE]

## Management tiers, in the FINN-fia spirit: M0 is the clean subset, each
## further tier trades contamination for sample size.
cm[, mgmt_tier := fcase(
  n_tal_harvest == 0L,   "M0",
  harvest_share <= 0.05, "M1",
  harvest_share <= 0.20, "M2",
  default                = "M3")]

setcolorder(cm, c("Tnr", "Enr", "Bl", "BL_NAME", "usable", "mgmt_tier",
                  "year_t0", "year_t1", "period_length"))
fwrite(cm, file.path(OUT_DIR, "corner_metadata.csv"))
message(sprintf("\n[corner_metadata] %d corners (%d usable)",
                nrow(cm), sum(cm$usable, na.rm = TRUE)))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 5) Tract summary -- how many usable corners does a tract have?
##    This decides whether a FINN site can be a Trakt with 4 patches or has to
##    be a single corner.
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

ts <- cm[usable == TRUE, .(
  n_corners_usable = .N,
  n_corners_M0     = sum(mgmt_tier == "M0"),
  n_corners_M0M1   = sum(mgmt_tier %in% c("M0", "M1")),
  ba_ha_t0         = mean(ba_ha_t0, na.rm = TRUE),
  period_length    = period_length[1]
), by = Tnr]
fwrite(ts, file.path(OUT_DIR, "trakt_summary.csv"))

message("\n[trakt_summary] usable corners per tract:")
print(ts[, .N, by = n_corners_usable][order(n_corners_usable)])
message("  ... restricted to tier M0:")
print(ts[, .N, by = n_corners_M0][order(n_corners_M0)])

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 6) Species ranking
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

usable_corners <- cm[usable == TRUE, .(Tnr, Enr)]
tr_u <- merge(trees, usable_corners, by = c("Tnr", "Enr"))

sp_rank <- tr_u[time == 0L & living == TRUE, .(
  n_corners_present = uniqueN(paste(Tnr, Enr)),
  ba_ha_total       = sum(ba_ha, na.rm = TRUE),
  n_ha_total        = sum(n_ha,  na.rm = TRUE),
  dbh_mean          = weighted.mean(dbh, n_ha, na.rm = TRUE)
), by = species_name][order(-n_corners_present)]

sp_rank[, ba_share      := ba_ha_total / sum(ba_ha_total)]
sp_rank[, ba_share_cum  := cumsum(ba_share)]
sp_rank[, rank          := .I]
fwrite(sp_rank, file.path(OUT_DIR, "species_rank.csv"))

message("\n[species_rank] top 15 by corner occurrence:")
print(head(sp_rank[, .(rank, species_name, n_corners_present,
                       ba_share = round(ba_share, 4),
                       ba_share_cum = round(ba_share_cum, 4))], 15))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 7) Management tiers + strictness grid
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

tier_def <- data.table(
  tier = c("M0", "M1", "M2", "M3"),
  definition = c(
    "no harvested tree over the period",
    "harvested basal area <= 5 % of the standing stock at period start",
    "harvested basal area <= 20 %",
    "any removal intensity"),
  n_corners = c(
    cm[usable == TRUE & mgmt_tier == "M0", .N],
    cm[usable == TRUE & mgmt_tier %in% c("M0", "M1"), .N],
    cm[usable == TRUE & mgmt_tier %in% c("M0", "M1", "M2"), .N],
    cm[usable == TRUE, .N])
)
tier_def[, share := round(n_corners / max(n_corners), 3)]
fwrite(tier_def, file.path(OUT_DIR, "management_tiers.csv"))
message("\n[management_tiers] cumulative corner counts:")
print(tier_def)

## Species contamination filter, as in FINN-fia: with V = the top-N species by
## corner occurrence, keep a corner only when V accounts for at least X of its
## basal area at the period start.
ba_by_corner_sp <- tr_u[time == 0L & living == TRUE,
                        .(ba = sum(ba_ha, na.rm = TRUE)),
                        by = .(Tnr, Enr, species_name)]
ba_by_corner <- ba_by_corner_sp[, .(ba_tot = sum(ba)), by = .(Tnr, Enr)]

grid <- CJ(n_species = c(5L, 10L, 15L, 20L, 25L, 30L),
           X = c(0.90, 0.95, 0.99, 1.00))
grid[, corners_kept := NA_integer_]
for (i in seq_len(nrow(grid))) {
  V   <- sp_rank[seq_len(min(grid$n_species[i], .N)), species_name]
  shr <- ba_by_corner_sp[species_name %in% V, .(ba_V = sum(ba)), by = .(Tnr, Enr)]
  shr <- merge(ba_by_corner, shr, by = c("Tnr", "Enr"), all.x = TRUE)
  shr[is.na(ba_V), ba_V := 0]
  grid$corners_kept[i] <- shr[ba_tot > 0 & ba_V / ba_tot >= grid$X[i], .N]
}
grid[, share_kept := round(corners_kept / nrow(ba_by_corner), 3)]
fwrite(grid, file.path(OUT_DIR, "strictness_grid.csv"))
message("\n[strictness_grid] corners kept by (|V|, contamination threshold X):")
print(dcast(grid, n_species ~ X, value.var = "corners_kept"))

message(sprintf("\nDone. Wrote 5 tables to %s/", OUT_DIR))
