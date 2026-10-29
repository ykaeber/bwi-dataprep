#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-bwi | scripts/40_import_new_schema.R
##
## REFIT part 1: build a tree panel from the Thuenen inventory-archive schema
## (data/raw/bwi4), for the BWI3->BWI4 (2012->2022) transition -- the drought
## decade, WITH the sink->source turn -- in the exact column format that
## scripts/04 already consumes (so the whole tested WZP -> FINN-input machinery
## is reused unchanged).
##
## Why this is worth doing over the ACCESS2003 pipeline:
##   * DIRECT growth. Survivors are matched across surveys by (Tnr, Enr,
##     tree_number), 96.7 % matched, 99.3 % species-consistent -- so growth is
##     measured dbh_2022 - dbh_2012, with none of the Sloboda back-modelling that
##     contaminated 22 % of the ACCESS2003 growth (bug #2 in docs/bwi_data_notes).
##   * The calibration period now CONTAINS the 2018-2020 drought.
##   * Cleaner fate codes (tree_status) and richer covariates for later use.
##
## The salvage subtlety: drought-killed trees were often salvage-logged (fate 12,
## harvest), which CENSORS them out of the natural-mortality signal. So a clean
## biological mortality calibration wants LOW-HARVEST corners, where dead trees
## stay recorded as standing dead (fate 4/5). This script records each corner's
## harvest fraction; scripts/41 uses it to define the calibration set.
##
## tree_status: 0 new (ingrowth), 1 survivor, 4/5 natural death, 12 removed
## (harvest), 9 untraceable, 6/7/8/10/11 out of sample.
##
##   Rscript scripts/40_import_new_schema.R
##
## Outputs:
##   outputs/new_schema/panel.csv           long panel (time 0/1), bwi_tree_panel format
##   outputs/new_schema/corners.csv         per corner: year_t0/t1, period, harvest_share
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

suppressPackageStartupMessages({ library(data.table) })
source("R/config.R")

RAW <- "data/raw/bwi4"
OUT <- file.path(OUT_DIR, "new_schema"); dir.create(OUT, recursive = TRUE,
                                                    showWarnings = FALSE)
YEAR_T0 <- 2012L; YEAR_T1 <- 2022L; PERIOD <- 10L
MIN_DBH_MM <- WZP_MIN_BHD * 10        # 70 mm = 7 cm WZP inclusion threshold

read_big <- function(name, select = NULL, lookup = FALSE) {
  f <- file.path(RAW, paste0(name, ".csv"))
  if (lookup) return(fread(f, select = select, showProgress = FALSE))
  cl <- tempfile()
  system2("awk", c("-F,", shQuote('NR==1{n=NF;print;next} NF==n{print}'), f),
          stdout = cl)
  d <- fread(cl, select = select, showProgress = FALSE); unlink(cl); d
}

#=#=#= species map: API code -> our 15; others labelled BA_<code> ============
splk <- read_big("lookup_tree_species", c("code", "genus", "species"),
                 lookup = TRUE)
splk[, latin := paste(genus, species)]
our <- fread(file.path(OUT_DIR, "finn", "species_dt.csv"))
splk <- merge(splk[, .(code, latin)], our, by.x = "latin", by.y = "species_name",
              all.x = TRUE)
splk[, species_name := fifelse(is.na(species), paste0("BA_", code), latin)]
sp_of <- splk[, setNames(species_name, code)]

#=#=#= plot + tree ==========================================================
plot <- read_big("plot", c("id", "interval_name", "cluster_name", "plot_name",
                           "federal_state"))
plot[, `:=`(Tnr = cluster_name, Enr = plot_name)]
setnames(plot, "id", "plot_id")
tree <- read_big("tree", c("plot_id", "tree_number", "tree_status",
                           "tree_species", "dbh"))
tree <- merge(tree, plot[, .(plot_id, interval_name, Tnr, Enr)], by = "plot_id")
tree[, species_name := sp_of[as.character(tree_species)]]
tree[, dbh_cm := dbh / 10]

## 2012 WZP population (dbh >= 7 cm)
b12 <- tree[interval_name == "bwi2012" & is.finite(dbh) & dbh >= MIN_DBH_MM,
            .(Tnr, Enr, tree_number, dbh12 = dbh_cm, sp12 = species_name)]
## 2022 records carry the 2012->2022 fate
b22 <- tree[interval_name == "bwi2022",
            .(Tnr, Enr, tree_number, status = tree_status,
              dbh22 = dbh_cm, sp22 = species_name)]

#=#=#= at-risk 2012 trees (present in 2012) + their fate ====================
atrisk <- b22[status %in% c(1L, 4L, 5L, 12L, 9L)]
atrisk <- merge(atrisk, b12, by = c("Tnr", "Enr", "tree_number"), all.x = TRUE)
## keep only those we can place in 2012 (have a matched 2012 dbh >= 7 cm)
atrisk <- atrisk[is.finite(dbh12)]
atrisk[, species_name := fifelse(!is.na(sp12), sp12, sp22)]
atrisk[, fate := fcase(status == 1L, "survivor",
                       status %in% c(4L, 5L), "dead",
                       status == 12L, "harvest",
                       status == 9L, "untraceable")]
## drop implausible growth for survivors (measurement/matching errors), mirroring
## the ACCESS2003 growth-threshold filter
atrisk[fate == "survivor", grow := dbh22 - dbh12]
atrisk <- atrisk[fate != "survivor" | (grow >= -2 & grow <= 40)]

## ingrowth: new WZP trees in 2022 (crossed the 7 cm threshold)
ingrowth <- b22[status == 0L & is.finite(dbh22) & dbh22 >= WZP_MIN_BHD,
                .(Tnr, Enr, tree_number, dbh22, species_name = sp22)]

message(sprintf("2012 WZP trees: %d | at-risk matched: %d | ingrowth: %d",
                nrow(b12), nrow(atrisk), nrow(ingrowth)))
message("fate distribution: ",
        paste(atrisk[, .N, fate][order(-N)][, sprintf("%s=%d", fate, N)],
              collapse = ", "))

#=#=#= assemble the long panel (time 0 and 1), bwi_tree_panel format =========
wzp_nha <- function(dbh_cm) WZP_BAF / bhd2ba(dbh_cm)   # stems/ha per tallied tree
atrisk[, treeKey := paste(Tnr, Enr, tree_number, sep = "_")]
ingrowth[, treeKey := paste(Tnr, Enr, tree_number, sep = "_")]

## time 0: every at-risk tree, alive, at its 2012 dbh
t0 <- atrisk[, .(Tnr, Enr, treeKey, species_name, time = 0L, dbh = dbh12,
                 status = "alive", living = TRUE, n_ha = wzp_nha(dbh12))]
## time 1: fate-dependent
t1 <- atrisk[, .(Tnr, Enr, treeKey, species_name, time = 1L,
                 dbh = fifelse(fate == "survivor", dbh22, dbh12),
                 status = fcase(fate == "survivor", "alive",
                                fate == "dead", "dead",
                                fate == "harvest", "harvest",
                                default = "gone"),
                 living = fate == "survivor",
                 n_ha = wzp_nha(fifelse(fate == "survivor", dbh22, dbh12)))]
## ingrowth appears only at time 1 (recruit)
ing1 <- ingrowth[, .(Tnr, Enr, treeKey, species_name, time = 1L, dbh = dbh22,
                     status = "new", living = TRUE, n_ha = wzp_nha(dbh22))]

panel <- rbindlist(list(t0, t1, ing1), use.names = TRUE)
## fate + mort + recruit flags, and the calendar the annual grid needs
panel <- merge(panel, atrisk[, .(treeKey, fate)], by = "treeKey", all.x = TRUE)
panel[treeKey %in% ingrowth$treeKey, fate := "recruit"]
## year_t0/year_t1 stay on the panel; period_length is supplied by corners in
## scripts/04, so leaving it here too would collide in the merge.
panel[, `:=`(year_t0 = YEAR_T0, year_t1 = YEAR_T1)]
panel[, year := fifelse(time == 0L, YEAR_T0, YEAR_T1)]
panel[, mort := status == "dead"]
panel[, is_recruit := fate == "recruit"]
## Guard against any residual duplicate tree x time rows (e.g. a data glitch):
## FINN's virtual-tree expansion assumes one row per (treeKey, time).
panel <- unique(panel, by = c("treeKey", "time"))
setorder(panel, Tnr, Enr, treeKey, time)
fwrite(panel, file.path(OUT, "panel.csv"))

#=#=#= corner table: harvest fraction + calendar ============================
corners <- atrisk[, .(n_atrisk = .N,
                      n_harvest = sum(fate == "harvest"),
                      n_dead = sum(fate == "dead")),
                  by = .(Tnr, Enr)]
corners[, `:=`(harvest_share = n_harvest / n_atrisk,
               mort_share = n_dead / n_atrisk,
               siteName = paste(Tnr, Enr, sep = "_"),
               patchName = paste(Tnr, Enr, sep = "_"),
               year_t0 = YEAR_T0, year_t1 = YEAR_T1, period_length = PERIOD)]
fwrite(corners, file.path(OUT, "corners.csv"))

message(sprintf("\nCorners: %d | median harvest share %.1f%% | low-harvest (<5%%): %d",
                nrow(corners), 100 * median(corners$harvest_share),
                corners[harvest_share < 0.05, .N]))
cat("\nWrote", file.path(OUT, "panel.csv"), "and corners.csv\n")
cat("Next: aggregate outputs/new_schema/panel.csv to the plot level, or build\n",
    "model inputs with scripts/04 the same way as for the ACCESS2003 panel\n",
    "(see docs/bwi_data_usage_guide.md, 'Working with the latest BWI').\n", sep = "")
