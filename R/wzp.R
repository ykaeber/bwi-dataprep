#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-bwi | R/wzp.R
##
## Turning angle-count (Winkelzaehlprobe) tallies into fixed-area FINN patches.
##
## THE PROBLEM
## FIA measures every tree on a fixed-area subplot, so a subplot IS a patch and
## a tree count is a count. The BWI instead uses a Winkelzaehlprobe with
## Zaehlfaktor 4 (BAF = 4 m^2/ha): a tree is tallied if its horizontal distance
## from the corner is below 25 * d (d in m), so inclusion probability is
## proportional to basal area. There is no fixed plot area, and big trees are
## massively over-sampled relative to small ones.
##
## What the design does give, per tallied tree i, is an unbiased density:
##   n_ha_i = K / g_i        stems per hectare represented by tree i
## with g_i the tree's basal area (m^2) and K the count factor. K = BAF = 4 for
## 98 % of trees; near a forest or stand edge Thuenen applies a Grenzkreis
## correction, so we take their published `N_ha` rather than recomputing 4/g.
## Each tallied tree represents exactly K m^2/ha of basal area, independent of
## its size.
##
## THE CONVERSION
## We declare a nominal patch area A (PATCH_SIZE_HA) and expand each tallied
## tree into round(n_ha_i * A) identical "virtual" trees that share its species,
## dbh at both dates, and fate. This preserves stem density and basal area per
## hectare, which is what FINN's state variables mean, and it lets us reuse
## FINN's own makeObsData() unchanged: it counts rows, and after expansion a
## row is a tree on a patch of area A, exactly as for FIA.
##
## Rounding is done by largest-remainder allocation within (patch x species), so
## the group total is preserved exactly even when individual trees round to
## zero. A tree is only at risk of vanishing when n_ha_i * A < 0.5, i.e.
## g_i > 8A; at A = 0.1 ha that is dbh > 101 cm and affects 0.14 % of BWI3
## trees, all of which are recovered by the remainder step at group level.
##
## THE SAMPLE-SIZE CAVEAT
## Expansion is right for densities and WRONG for sample sizes: one measured
## tree becoming 40 virtual trees would tell a binomial likelihood it has 40
## independent observations. FINN's obs table carries `n_at_risk` and
## `growth_n` precisely so the loss can weight by sample size, so we compute
## the rates from the expanded (density-weighted) table but reset the counts to
## the number of TALLIED trees -- see wzp_effective_n() and its use in
## scripts/04_prepare_finn_inputs.R. Rate: design-weighted. n: what was really
## measured.
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

suppressPackageStartupMessages(library(data.table))

#' Stems per hectare represented by a tallied angle-count tree.
#'
#' @param g_m2 basal area of the stem, m^2
#' @param baf  Zaehlfaktor, m^2/ha
#' @return numeric, stems/ha
#'
#' Prefer the BWI's own `N_ha` column where it exists: it carries the
#' Grenzkreis (edge) correction that this formula does not.
wzp_representation <- function(g_m2, baf = WZP_BAF) {
  ifelse(is.finite(g_m2) & g_m2 > 0, baf / g_m2, NA_real_)
}

#' Integer allocation that preserves the group total (largest remainder).
#'
#' Every element gets floor(x); the round(sum(x)) - sum(floor(x)) elements with
#' the largest fractional parts are bumped by one. Ties break on the supplied
#' order, so the result is deterministic.
#'
#' @param x non-negative numeric
#' @return integer vector, sum(result) == round(sum(x))
largest_remainder <- function(x) {
  x  <- as.numeric(x)
  x[!is.finite(x) | x < 0] <- 0
  fl <- floor(x)
  k  <- as.integer(round(sum(x)) - sum(fl))
  if (k > 0L) {
    ## order() is stable, so equal remainders keep input order.
    ord <- order(x - fl, decreasing = TRUE)
    bump <- ord[seq_len(min(k, length(x)))]
    fl[bump] <- fl[bump] + 1
  }
  as.integer(fl)
}

#' Allocate integer tree counts to tallied trees within groups.
#'
#' @param dt        data.table of tallied trees
#' @param n_ha_col  column name holding stems/ha
#' @param by        grouping columns within which the total is preserved
#' @param patch_size_ha nominal patch area
#' @return `dt` with an integer column `n_trees` added (by reference)
wzp_allocate <- function(dt, n_ha_col = "n_ha", by = c("patchKey", "species_name"),
                         patch_size_ha = PATCH_SIZE_HA) {
  stopifnot(is.data.table(dt), n_ha_col %in% names(dt), all(by %in% names(dt)))
  ## `dt[[n_ha_col]]` rather than get(): inside a data.table j-expression the
  ## argument name would be shadowed by the identically named column.
  set(dt, j = "n_exact", value = dt[[n_ha_col]] * patch_size_ha)
  dt[, n_trees := largest_remainder(n_exact), by = by]
  dt[]
}

#' Expand tallied trees into virtual trees on a fixed-area patch.
#'
#' Each input row is replicated `n_trees` times. The replicates are identical
#' apart from `treeName`, which appends the replicate index so FINN can chain
#' the same virtual tree across the two inventories.
#'
#' @param dt       data.table with an integer `n_trees` column and a `treeKey`
#'                 column uniquely identifying the tallied tree
#' @param keep     columns to carry onto the virtual trees
#' @param index_by columns within which the replicate index restarts. When the
#'                 input holds one row per (tree, date) this MUST include the
#'                 date column: numbering across dates would give replicate 1 at
#'                 t0 and replicate 1 at t1 to different virtual trees, so no
#'                 virtual tree would exist at both dates and every growth rate
#'                 and survival chain would come out NA.
#' @return data.table with one row per virtual tree and a `treeName` column
wzp_expand <- function(dt, keep = names(dt), index_by = "treeKey") {
  stopifnot(is.data.table(dt), "n_trees" %in% names(dt),
            all(index_by %in% names(dt)), "treeKey" %in% names(dt))
  keep <- setdiff(unique(c(keep, index_by, "treeKey")), c("n_trees", "n_exact"))
  out  <- dt[n_trees > 0L, .SD[rep(seq_len(.N), n_trees)], .SDcols = keep]
  ## Replicate index within each tallied tree (and date) -> a virtual tree keeps
  ## its identity across the inventories.
  out[, rep_i := seq_len(.N), by = c(index_by)]
  out[, treeName := paste0(treeKey, "#", rep_i)]
  out[, rep_i := NULL]
  out[]
}

#' Basal-area conservation check for an expansion.
#'
#' Compares the design-based BA/ha (K summed over tallied trees, which is the
#' unbiased angle-count estimator) with the BA/ha implied by the rounded
#' virtual trees. A large discrepancy means the patch area is too small for the
#' size distribution.
#'
#' @param tallied data.table of tallied trees with `n_ha`, `dbh` (cm),
#'                `n_trees` and the grouping key
#' @param by      grouping columns (typically the patch key)
#' @param patch_size_ha nominal patch area
#' @return data.table with `ba_ha_design`, `ba_ha_expanded`, `rel_err`
wzp_ba_conservation <- function(tallied, by = "patchKey",
                                patch_size_ha = PATCH_SIZE_HA) {
  stopifnot(all(c("n_ha", "dbh", "n_trees") %in% names(tallied)))
  tallied[, .(
    ba_ha_design   = sum(n_ha * bhd2ba(dbh), na.rm = TRUE),
    ba_ha_expanded = sum(n_trees * bhd2ba(dbh), na.rm = TRUE) / patch_size_ha
  ), by = by][, rel_err := fifelse(ba_ha_design > 0,
                                   (ba_ha_expanded - ba_ha_design) / ba_ha_design,
                                   NA_real_)][]
}

#' Sample sizes behind each observation, counted in TALLIED trees.
#'
#' FINN weights the mortality (binomial) and growth losses by `n_at_risk` and
#' `growth_n`. After expansion those counts refer to virtual trees and would
#' overstate the information by the expansion factor, so we recompute them from
#' the measured trees and substitute them back into obs_dt.
#'
#' `n_died_eff` is derived from the design-weighted rate rather than counted
#' directly, so that n_died_eff / n_at_risk reproduces the density-weighted
#' mortality while the denominator stays the true sample size.
#'
#' WHICH TREES MAY ENTER EACH COHORT
#' Two exclusions matter, and both bite hard on BWI data:
#'
#' * Harvested stems are CENSORED from the mortality denominator. They left the
#'   stand by a competing risk FINN does not model, so counting them as
#'   survivors would deflate mortality.
#'
#' * Laterally ingrown stems must be excluded from BOTH cohorts.
#'   - Mortality: a stem is classed as lateral ingrowth precisely because it was
#'     tallied at t1, so it is conditioned on having survived and cannot die.
#'     Leaving it in the denominator is immortal time bias. It is 26.6 % of the
#'     BWI at-risk cohort and 42.9 % of the 7-15 cm class, and it deflates
#'     observed mortality by up to 1.8x -- enough to inverse the size gradient
#'     and make small trees look like the safest stems in the forest.
#'   - Growth: its t0 diameter is a Sloboda back-model, not a measurement, so
#'     its apparent increment is Thuenen's growth model rather than an
#'     observation. Matched on size it runs 13 % higher than real survivors and
#'     visibly smoother (sd 0.0053 vs 0.0075).
#'
#' @param tallied data.table of tallied trees, one row per measured tree, with
#'                columns: the `by` keys, `n_ha`, `at_risk` (logical),
#'                `died` (logical), `growth_ok` (logical), and optionally
#'                `rel_growth` (relative diameter increment over the period)
#' @param by      grouping columns matching obs_dt's key
#' @return data.table with n_at_risk, n_died_eff, mort_w, growth_n, growth_w
wzp_effective_n <- function(tallied, by = c("siteName", "patchName", "year",
                                            "species_name")) {
  stopifnot(all(c("n_ha", "at_risk", "died", "growth_ok") %in% names(tallied)))
  has_growth  <- "rel_growth" %in% names(tallied)
  ## harvest is an optional SECOND removal event on the same at-risk population
  ## (competing risk with natural death); present only when managed corners are in
  ## the set. When absent everything below reduces to the old behaviour exactly.
  has_harvest <- "harvested" %in% names(tallied)
  out <- tallied[, {
    w  <- n_ha * at_risk
    sw <- sum(w, na.rm = TRUE)
    gw <- n_ha * growth_ok * is.finite(if (has_growth) rel_growth else NA_real_)
    sg <- sum(gw, na.rm = TRUE)
    .(n_at_risk = sum(at_risk, na.rm = TRUE),
      ## design-weighted mortality: weight each measured tree by the stems/ha
      ## it stands for, so the rate is per stem rather than per basal area
      mort_w    = if (sw > 0) sum(w * died, na.rm = TRUE) / sw else NA_real_,
      ## design-weighted harvest rate on the SAME at-risk population
      harv_w    = if (has_harvest && sw > 0)
                    sum(w * harvested, na.rm = TRUE) / sw else NA_real_,
      growth_n  = sum(growth_ok & (if (has_growth) is.finite(rel_growth) else TRUE),
                      na.rm = TRUE),
      ## design-weighted mean relative growth over the same measured trees
      growth_w  = if (has_growth && sg > 0)
                    sum(gw * rel_growth, na.rm = TRUE) / sg else NA_real_)
  }, by = by]
  out[, n_died_eff    := fifelse(is.finite(mort_w),
                                 as.integer(round(mort_w * n_at_risk)), 0L)]
  out[, n_harvest_eff := fifelse(is.finite(harv_w),
                                 as.integer(round(harv_w * n_at_risk)), 0L)]
  out[]
}
