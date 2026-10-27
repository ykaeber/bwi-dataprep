#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-bwi | R/bwi_trees.R
##
## Assembling the BWI2 -> BWI3 tree panel: one row per (tallied tree, date),
## with the status / living / mort semantics FINN's makeObsData() expects.
##
## SOURCE TABLES
##   b23_baeume_a_s  state at the START of the period (2002). Measured where the
##                   tree was a BWI2 sample tree, back-modelled (Sloboda) where
##                   it was not. Key (Tnr, Enr, Bnr).
##   b23_baeume_m_s  same key, carries `Pk` -- the tree's fate over the period.
##   b3_baeume       state at the END of the period (2012), measured.
##
## THE INGROWTH TRAP
## `Pk == 0` ("neuer Probebaum") does NOT mean "recruit". In an angle-count
## sample a tree is tallied when its distance from the corner is under 25*d, so
## the inclusion zone grows with d^2: a tree that was merely too far away in
## 2002 is tallied in 2012 without anything having happened to it biologically.
##
## In the BWI2->BWI3 data 113,237 trees carry Pk == 0, but their back-modelled
## 2002 diameters show only 8.9 % were below the 7 cm caliper threshold. The
## other 91.1 % already existed as measurable trees (median back-modelled 2002
## dbh 23.6 cm) and simply grew into the inclusion zone. The design predicts
## this: ~15 % diameter growth over the period inflates the inclusion area by
## ~32 %, so ~24 % of the final tally should be newly included; 29.5 % is
## observed.
##
## Treating all of Pk == 0 as recruitment would overstate FINN's regeneration
## response by roughly an order of magnitude. So we split on the back-modelled
## diameter, which is exactly what Thuenen produced it for:
##   back-modelled dbh <  7 cm -> genuine recruit; appears only at t1, "new"
##   back-modelled dbh >= 7 cm -> lateral ingrowth; present at BOTH dates,
##                                using the back-modelled state at t0
##
## FATE CODES (x_pk, see R/config.R)
##   0 new sample tree   -> recruit or lateral ingrowth, split as above
##   1 survivor          -> alive at both dates
##   2 selective removal | 3 clearfell  -> HARVEST, not mortality
##   4 fallen dead       | 5 standing dead -> natural mortality
##   6,8,9,10,11         -> out of sample / untraceable; dropped, and the corner
##                          is flagged so callers can exclude it (silently
##                          dropping a live tree would read as mortality)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

suppressPackageStartupMessages(library(data.table))

#' Per-Trakt field survey dates and the realised period length.
#'
#' `Datum_feld` in b2_tab_w / b3_tab_w is the date the field crew worked the
#' tract. Tracts were spread over two field seasons at each inventory, so the
#' realised interval varies by roughly +/- 1.5 years around the nominal 10.
#'
#' @return data.table(Tnr, year_t0, year_t1, period_length)
bwi_periods <- function() {
  parse_year <- function(x) {
    ## mdbtools writes Access datetimes as "MM/DD/YY HH:MM:SS".
    y <- suppressWarnings(as.integer(sub("^\\d+/\\d+/(\\d+).*$", "\\1",
                                         as.character(x))))
    fifelse(is.na(y), NA_integer_, fifelse(y < 50L, 2000L + y, 1900L + y))
  }
  ## min() over an all-NA integer group returns Inf (a double) and would make
  ## the column type inconsistent across groups, so reduce explicitly.
  min_year <- function(x) {
    v <- parse_year(x)
    v <- v[is.finite(v)]
    if (length(v)) min(v) else NA_integer_
  }
  t2 <- bwi("b2_tab_w", select = c("Tnr", "Datum_feld"))
  t3 <- bwi("b3_tab_w", select = c("Tnr", "Datum_feld"))
  t2 <- t2[, .(year_t0 = min_year(Datum_feld)), by = Tnr]
  t3 <- t3[, .(year_t1 = min_year(Datum_feld)), by = Tnr]
  p  <- merge(t2, t3, by = "Tnr", all = TRUE)
  ## Fall back to the nominal calendar where a tract has no usable date.
  p[is.na(year_t0), year_t0 := BWI_YEARS[["b2"]]]
  p[is.na(year_t1), year_t1 := BWI_YEARS[["b3"]]]
  p[, period_length := year_t1 - year_t0]
  p[!is.finite(period_length) | period_length <= 0,
    `:=`(year_t0 = BWI_YEARS[["b2"]], year_t1 = BWI_YEARS[["b3"]],
         period_length = BWI_PERIOD_YEARS)]
  p[]
}

#' Assemble the BWI2 -> BWI3 tallied-tree panel.
#'
#' @param min_dbh_cm    caliper threshold; also the recruit/ingrowth split
#' @param drop_lying    drop trees recorded as lying (Bs = 9)
#' @param species_map   data.table(Ba, species_name); defaults to the x_ba lookup
#' @return list(trees = panel, corner_flags = per-corner QC flags)
#'
#' `trees` columns:
#'   Tnr, Enr, Bnr, treeKey, siteName, patchName, year, Ba, species_name,
#'   dbh (cm), n_ha (stems/ha), status, living, mort, pk, fate, is_recruit
bwi_tree_panel <- function(min_dbh_cm = WZP_MIN_BHD,
                           drop_lying  = TRUE,
                           species_map = NULL) {

  if (is.null(species_map)) {
    species_map <- bwi("species_lookup")[, .(Ba, species_name)]
  }
  min_bhd_mm <- min_dbh_cm * 10

  a <- bwi("b23_baeume_a_s",
           select = c("Tnr", "Enr", "Bnr", "Av", "Ba", "Bhd", "N_ha", "Bs"))
  m <- bwi("b23_baeume_m_s", select = c("Tnr", "Enr", "Bnr", "Pk"))
  e <- bwi("b3_baeume",
           select = c("Tnr", "Enr", "Bnr", "Av", "Ba", "Bhd", "N_ha", "Bs", "Bz"))

  setnames(a, c("Av", "Ba", "Bhd", "N_ha", "Bs"),
           c("Av_t0", "Ba_t0", "Bhd_t0", "N_ha_t0", "Bs_t0"))
  setnames(e, c("Av", "Ba", "Bhd", "N_ha", "Bs", "Bz"),
           c("Av_t1", "Ba_t1", "Bhd_t1", "N_ha_t1", "Bs_t1", "Bz_t1"))

  d <- merge(a, m, by = c("Tnr", "Enr", "Bnr"), all = TRUE)
  d <- merge(d, e, by = c("Tnr", "Enr", "Bnr"), all.x = TRUE)

  ## --- fate classification -------------------------------------------------
  d[, fate := fcase(
    Pk %in% PK_HARVEST,  "harvest",
    Pk %in% PK_MORT_NAT, "dead",
    Pk == PK_SURVIVOR,   "survivor",
    Pk == PK_NEW,        "entered",
    default = "invalid")]

  ## Split the "entered" class on the back-modelled start-of-period diameter.
  d[fate == "entered", fate := fifelse(is.finite(Bhd_t0) & Bhd_t0 >= min_bhd_mm,
                                       "ingrowth_lateral", "recruit")]

  ## --- corner-level QC flags ----------------------------------------------
  ## A corner carrying trees we cannot classify is suspect: dropping a live tree
  ## would show up downstream as mortality that never happened.
  corner_flags <- d[, .(
    n_trees          = .N,
    n_invalid        = sum(fate == "invalid"),
    n_harvest        = sum(fate == "harvest"),
    n_dead           = sum(fate == "dead"),
    n_recruit        = sum(fate == "recruit"),
    n_ingrowth_lat   = sum(fate == "ingrowth_lateral"),
    n_survivor       = sum(fate == "survivor")
  ), by = .(Tnr, Enr)]
  corner_flags[, has_invalid := n_invalid > 0]
  corner_flags[, has_harvest := n_harvest > 0]

  d <- d[fate != "invalid"]

  ## --- long panel ----------------------------------------------------------
  ## t0 row: every tree that existed at the start of the period. Recruits are
  ## excluded by construction -- they were below the caliper threshold.
  t0 <- d[fate != "recruit", .(
    Tnr, Enr, Bnr, pk = Pk, fate,
    time  = 0L,
    Ba    = Ba_t0,
    Bhd   = Bhd_t0,
    n_ha  = N_ha_t0,
    Bs    = Bs_t0,
    ## Everything alive at t0 is "alive"; makeObsData derives growth and the
    ## at-risk cohort from the transition, not from this label.
    status = "alive"
  )]

  ## t1 row: survivors and lateral ingrowth are alive, recruits are new, and
  ## natural deaths are recorded as dead carrying their last known diameter.
  ## Harvested trees get a t1 row too, tagged so callers can decide; the
  ## management tiers in scripts/01 exclude them by default.
  t1 <- d[, .(
    Tnr, Enr, Bnr, pk = Pk, fate,
    time = 1L,
    Ba   = fifelse(is.finite(Ba_t1), Ba_t1, Ba_t0),
    Bhd  = fifelse(is.finite(Bhd_t1), Bhd_t1, Bhd_t0),
    n_ha = fifelse(is.finite(N_ha_t1) & N_ha_t1 > 0, N_ha_t1, N_ha_t0),
    Bs   = fifelse(is.finite(Bs_t1), Bs_t1, Bs_t0),
    status = fcase(fate == "recruit",          "new",
                   fate %in% c("survivor", "ingrowth_lateral"), "alive",
                   fate == "dead",             "dead",
                   fate == "harvest",          "harvest")
  )]

  p <- rbindlist(list(t0, t1), use.names = TRUE)

  ## --- filters -------------------------------------------------------------
  p <- p[is.finite(Bhd) & Bhd >= min_bhd_mm]
  p <- p[is.finite(n_ha) & n_ha > 0]          # drops Bz == 0 (not in stand)
  if (drop_lying) p <- p[is.na(Bs) | Bs != BS_LYING]

  ## --- naming, units, FINN semantics --------------------------------------
  p[, dbh := mm2cm(Bhd)]
  p <- merge(p, species_map, by = "Ba", all.x = TRUE)
  p[is.na(species_name), species_name := paste0("BA_", Ba)]

  p[, siteName  := Tnr]
  p[, patchName := Enr]
  p[, treeKey   := paste(Tnr, Enr, Bnr, sep = "_")]
  p[, living    := status %chin% c("alive", "new")]
  p[, is_recruit := fate == "recruit"]

  ## Calendar years, so period_length is real time rather than an index.
  per <- bwi_periods()
  p <- merge(p, per[, .(Tnr, year_t0, year_t1)], by = "Tnr", all.x = TRUE)
  p[, year := fifelse(time == 0L, year_t0, year_t1)]
  p[, c("year_t0", "year_t1") := NULL]

  ## `mort` marks a death event at this record: dead now, alive before. The t0
  ## rows are always alive, so this is only ever TRUE on a t1 row.
  setorder(p, treeKey, time)
  p[, status_before := shift(status, 1L, type = "lag"), by = treeKey]
  p[, mort := status == "dead" & status_before %chin% c("alive", "new")]
  p[, status_before := NULL]

  list(trees = p[], corner_flags = corner_flags[])
}
