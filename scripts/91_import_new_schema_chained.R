#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-bwi | scripts/91_import_new_schema_chained.R
##
## The BWI3 -> BWI4 (2012 -> 2022) tree panel of scripts/40, rebuilt with two
## fate defects fixed. scripts/40 stays as it is (it feeds two-interval-v3/v4);
## this copy writes to its own folder and feeds two-interval-v5 (scripts/92)
## and the three-period CI2017 set (scripts/93).
##
## DEFECT 1 -- FRESHLY DEAD TREES WERE SURVIVORS
## BWI2022 field manual, sections 5.3.3 and 5.3.14: a tree dead "presumably for
## less than 12 months" is recorded like a live sample tree (BK = tree_status 1)
## and carries the flag "frisch abgestorben" (tree.csv column `damage_dead`).
## scripts/40 maps status 1 to survivor, so these trees counted as living. Here:
##   * 2022 status 1 + damage_dead  -> natural death in the decade;
##   * 2022 status 0 + damage_dead  -> a new tree that is already dead: not
##     ingrowth (it never enters the living 2022 stand);
##   * 2012 damage_dead (dead at the START survey) -> not part of the 2012
##     stand and not at risk in 2012-2022 (it died in the previous decade).
## The freshly dead tree keeps its measured 2022 diameter at t1, so the
## diameter screens (scripts/40 survivor window, scripts/70 site screen) see
## exactly the trees they saw in v3.
##
## DEFECT 2 -- 2012-2017 LOSSES AT THE CI2017 CORNERS WERE MISSING
## A share of the corners was remeasured in the carbon inventory 2017 (archive
## interval_name "ci2017", field work 2016-2018). For those corners BWI2022
## refers to CI2017 as the previous inventory (manual, footnotes 1/2 of
## chapter 3: the previous inventory "can be BWI 2012, LWI Brandenburg 2014,
## LWI NRW 2012-2014, Kohlenstoffinventur 2017"). A tree removed (CI2017 code
## 12) or dead (4, 5) by 2017 has NO bwi2022 row, and CI2017 ingrowth (code 0)
## appears in 2022 as code 1 (or 12/4/5). scripts/40 took every fate from the
## bwi2022 rows, so these losses vanished -- from the fates AND from the 2012
## start stand (its time-0 rows are the at-risk trees). Here the fate is
## chained: a 2012 tree's decade fate is its CI2017 fate if it was lost by 2017
## (CI2017 code 4, 5, 12 or 9), otherwise its BWI2022 fate (with defect 1
## fixed). CI2017 ingrowth alive in 2022 (CI2017 code 0, 2022 code 1, not
## freshly dead, not a 2012 tree) is decade ingrowth like a 2022 code-0 tree.
## CI2017 records no `damage_dead` flag (every value is "f"), so a tree freshly
## dead in 2017 shows up through its 2022 code (4/5/12) instead.
##
## Everything else is scripts/40 unchanged: the 2012 population is every 2012
## record with dbh >= 7 cm whose (chained) fate is one of 1/4/5/12/9, the
## species comes from 2012 (else 2022), survivors with a decade diameter change
## outside -2..40 cm are dropped, ingrowth is every new tree >= 7 cm in 2022,
## n_ha = BAF / basal area, the long panel has time 0 (2012) and time 1 (2022).
## The archive reader is the R version of scripts/40's awk filter (it drops the
## same 59 tree / 16 plot lines that the paged download glued together;
## tree.csv contains no quoted fields, so the field counts agree).
##
## Inputs
##   data/raw/bwi4/{plot,tree}.csv      Thuenen inventory archive (scripts/fetch_bwi4*.sh)
##   outputs/finn/species_dt.csv        ACCESS2003 species list, lumping as scripts/40
##   outputs/new_schema/panel.csv       scripts/40 output, for the regression check only
##
## Outputs  outputs/new_schema_chained/
##   panel.csv        long panel (time 0/1) in scripts/40's format, fates fixed
##   corners.csv      per corner: harvest/mortality shares, calendar (as scripts/40)
##   tree_fates.csv   every 2012 tree (dbh >= 7 cm) on every corner: raw codes of
##                    the three surveys, v3 fate (scripts/40 rule), v5 fate, the
##                    source of the v5 fate and the sub-interval of the loss
##   tree_chain_ci2017.csv  every tree record of the CI2017 corners, the three
##                    surveys side by side (raw codes, dbh, damage_dead); the
##                    three-period set (scripts/93) builds its fates from it
##   IMPORT.md        counts and the regression check against scripts/40
##
## Helpers sourced: R/config.R (paths, WZP_BAF, WZP_MIN_BHD, bhd2ba).
##
## Run from the project root:  Rscript scripts/91_import_new_schema_chained.R
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

suppressPackageStartupMessages({ library(data.table) })
source("R/config.R")

#=#=# PARAMS #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
PARAMS <- list(
  ## raw archive and the output folder (a new folder: outputs/new_schema keeps
  ## the scripts/40 panel that v3/v4 were built from)
  bwi4_dir = file.path(BWI_RAW_DIR, "bwi4"),
  out_dir  = file.path(OUT_DIR, "new_schema_chained"),

  ## nominal calendar of the decade, as scripts/40 (realised years: scripts/92)
  year_t0 = 2012L, year_t1 = 2022L, period = 10L,

  ## WZP inclusion threshold (mm), as scripts/40
  min_dbh_mm = WZP_MIN_BHD * 10,

  ## end-of-decade codes that put a 2012 tree into the at-risk population
  ## (scripts/40): 1 survivor, 4/5 natural death, 12 removed, 9 untraceable
  atrisk_codes = c(1L, 4L, 5L, 12L, 9L),
  ## CI2017 codes that END a tree's decade by 2017 (defect 2): natural death,
  ## removal, untraceable. Codes 6/8/10/11 (no longer a sample tree, no loss)
  ## do not; such a tree keeps its 2022 fate, as it would without CI2017.
  ci_loss_codes = c(4L, 5L, 12L, 9L),

  ## scripts/40 survivor screen: decade diameter change (cm) outside this window
  ## is a matching/measurement error and the tree is dropped
  survivor_growth_cm = c(-2, 40),

  ## the flag value of "frisch abgestorben" in tree.csv (logical stored as t/f)
  fresh_dead_value = "t"
)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

OUT <- PARAMS$out_dir
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
KEY <- c("Tnr", "Enr", "tree_number")

#' Read an inventory-archive CSV, dropping the lines the paged download glued
#' together (they carry twice the field count; scripts/40 drops the same lines
#' with awk, scripts/70 with this function).
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
wzp_nha <- function(dbh_cm) WZP_BAF / bhd2ba(dbh_cm)   # stems/ha per tallied tree

#=#=#= species map: API code -> our 15; others labelled BA_<code> (scripts/40) =
splk <- fread(file.path(PARAMS$bwi4_dir, "lookup_tree_species.csv"),
              select = c("code", "genus", "species"), showProgress = FALSE)
splk[, latin := paste(genus, species)]
our  <- fread(file.path(OUT_DIR, "finn", "species_dt.csv"))
splk <- merge(splk[, .(code, latin)], our, by.x = "latin",
              by.y = "species_name", all.x = TRUE)
splk[, species_name := fifelse(is.na(species), paste0("BA_", code), latin)]
sp_of <- splk[, setNames(species_name, code)]

#=#=#= plot + tree ==========================================================
plot <- read_archive("plot", c("id", "interval_name", "cluster_name", "plot_name"))
plot[, `:=`(Tnr = cluster_name, Enr = plot_name)]
setnames(plot, "id", "plot_id")
tree <- read_archive("tree", c("plot_id", "tree_number", "tree_status",
                               "tree_species", "dbh", "damage_dead", "distance"))
tree <- merge(tree, plot[, .(plot_id, interval_name, Tnr, Enr)], by = "plot_id")
tree[, species_name := sp_of[as.character(tree_species)]]
tree[, dbh_cm := dbh / 10]
tree[, fresh_dead := as.character(damage_dead) %chin% PARAMS$fresh_dead_value]
stopifnot(!anyDuplicated(tree, by = c(KEY, "interval_name")))

b12 <- tree[interval_name == "bwi2012" & is.finite(dbh) & dbh >= PARAMS$min_dbh_mm,
            .(Tnr, Enr, tree_number, s12 = tree_status, dbh12 = dbh_cm,
              sp12 = species_name, fd12 = fresh_dead, dist12 = distance)]
c17 <-tree[interval_name == "ci2017",
            .(Tnr, Enr, tree_number, s17 = tree_status, dbh17 = dbh_cm,
              sp17 = species_name, fd17 = fresh_dead, dist17 = distance)]
b22 <- tree[interval_name == "bwi2022",
            .(Tnr, Enr, tree_number, s22 = tree_status, dbh22 = dbh_cm,
              sp22 = species_name, fd22 = fresh_dead)]
ci_corners <- unique(c17[, .(Tnr, Enr)])[, ci := TRUE]
message(sprintf("corners with CI2017 tree records: %d; CI2017 records flagged freshly dead: %d",
                nrow(ci_corners), c17[fd17 == TRUE, .N]))

#=#=#= the three surveys side by side for the CI2017 corners (scripts/93) ====
## (b12 holds only dbh >= 7 cm; the chain carries every 2012 record)
b12_all <- tree[interval_name == "bwi2012",
                .(Tnr, Enr, tree_number, s12 = tree_status, dbh12 = dbh_cm,
                  sp12 = species_name, fd12 = fresh_dead, dist12 = distance)]
chain_ci <- merge(merge(b12_all, c17, by = KEY, all = TRUE), b22, by = KEY, all = TRUE)
chain_ci <- chain_ci[ci_corners[, .(Tnr, Enr)], on = c("Tnr", "Enr"), nomatch = NULL]
setorderv(chain_ci, KEY)
fwrite(chain_ci, file.path(OUT, "tree_chain_ci2017.csv"))
rm(b12_all)

#=#=#= fates: v3 rule (scripts/40) and the chained v5 rule ====================
ch <- merge(b12, c17, by = KEY, all.x = TRUE)
ch <- merge(ch, b22, by = KEY, all.x = TRUE)
ch <- merge(ch, ci_corners, by = c("Tnr", "Enr"), all.x = TRUE)
ch[is.na(ci), ci := FALSE]

code_fate <- function(s) fcase(s == 1L, "survivor", s %in% c(4L, 5L), "dead",
                               s == 12L, "harvest", s == 9L, "untraceable",
                               default = NA_character_)
## v3: the bwi2022 code alone
ch[, fate_v3 := code_fate(s22)]
ch[, grow := dbh22 - dbh12]
ch[, drop_v3 := fate_v3 %chin% "survivor" &
       !(grow >= PARAMS$survivor_growth_cm[1] & grow <= PARAMS$survivor_growth_cm[2])]
## v5: lost by 2017 at a CI2017 corner -> the CI2017 fate, else the 2022 fate
## with the freshly dead flag; a tree freshly dead in 2012 is not at risk
ch[, lost17 := ci & s17 %in% PARAMS$ci_loss_codes]
ch[, fate := fifelse(lost17, code_fate(s17), code_fate(s22))]
ch[!lost17 & fate == "survivor" & fd22 %in% TRUE, fate := "dead"]
ch[, fate_src := fcase(lost17, "ci2017",
                       !lost17 & s22 == 1L & fd22 %in% TRUE, "bwi2022_fresh_dead",
                       !is.na(fate), "bwi2022", default = NA_character_)]
ch[fd12 %in% TRUE, `:=`(fate = NA_character_, fate_src = "dead_at_2012")]
## the survivor screen applies to the same trees as in v3: a 2022 code 1 that
## was not lost by 2017 (survivor or freshly dead)
ch[, drop_v5 := !lost17 & s22 %in% 1L &
       !(grow >= PARAMS$survivor_growth_cm[1] & grow <= PARAMS$survivor_growth_cm[2])]
ch[, loss_period := fcase(!fate %chin% c("dead", "harvest", "untraceable"), NA_character_,
                          lost17, "2012-2017",
                          ci, "2017-2022",
                          default = "2012-2022")]

## dist12 = 2012 horizontal distance (cm) -- a BAF-4 sample tree lies within
## 25 * dbh (cm); kept for the PREPARATION.md diagnostics (Brandenburg's BWI2012
## also tallied a BAF-1 sample, whose extra trees have no BWI2022 row)
fwrite(ch[, .(Tnr, Enr, tree_number, ci, s12, s17, s22, fd12, fd22, dbh12, dbh17, dbh22,
              sp12, sp22, dist12, fate_v3, drop_v3, fate, fate_src, drop_v5, loss_period)],
       file.path(OUT, "tree_fates.csv"))

#=#=#= at-risk 2012 trees + ingrowth (scripts/40 logic on the v5 fates) =======
atrisk <- ch[!is.na(fate) & !drop_v5]
atrisk[, species_name := fifelse(!is.na(sp12), sp12, sp22)]
## t1 diameter: survivors and freshly dead trees were measured in 2022; every
## other loss keeps its 2012 diameter (scripts/40)
atrisk[, dbh_t1 := fifelse(fate_src %chin% "bwi2022" & fate == "survivor" |
                             fate_src %chin% "bwi2022_fresh_dead", dbh22, dbh12)]

## ingrowth: new 2022 trees (code 0, not freshly dead) and CI2017 ingrowth that
## is alive in 2022 (CI2017 code 0 -> 2022 code 1) and is not a 2012 tree
ing22 <- b22[s22 == 0L & is.finite(dbh22) & dbh22 >= WZP_MIN_BHD & !fd22]
ing17 <- merge(b22[s22 == 1L & is.finite(dbh22) & dbh22 >= WZP_MIN_BHD & !fd22],
               c17[s17 == 0L, .(Tnr, Enr, tree_number)], by = KEY)
ing17 <- ing17[!b12, on = KEY]
ingrowth <- rbind(ing22[, .(Tnr, Enr, tree_number, dbh22, species_name = sp22, src = "bwi2022")],
                  ing17[, .(Tnr, Enr, tree_number, dbh22, species_name = sp22, src = "ci2017")])
## v3 ingrowth, for the comparison
n_ing_v3 <- b22[s22 == 0L & is.finite(dbh22) & dbh22 >= WZP_MIN_BHD, .N]

message(sprintf("2012 WZP trees: %d | at-risk v3: %d -> v5: %d | ingrowth v3: %d -> v5: %d (CI2017: %d)",
                nrow(b12), ch[!is.na(fate_v3) & !drop_v3, .N], nrow(atrisk),
                n_ing_v3, nrow(ingrowth), nrow(ing17)))

#=#=#= assemble the long panel (time 0 and 1), bwi_tree_panel format =========
atrisk[, treeKey := paste(Tnr, Enr, tree_number, sep = "_")]
ingrowth[, treeKey := paste(Tnr, Enr, tree_number, sep = "_")]

t0 <- atrisk[, .(Tnr, Enr, treeKey, species_name, time = 0L, dbh = dbh12,
                 status = "alive", living = TRUE, n_ha = wzp_nha(dbh12))]
t1 <- atrisk[, .(Tnr, Enr, treeKey, species_name, time = 1L,
                 dbh = dbh_t1,
                 status = fcase(fate == "survivor", "alive",
                                fate == "dead", "dead",
                                fate == "harvest", "harvest",
                                default = "gone"),
                 living = fate == "survivor",
                 n_ha = wzp_nha(dbh_t1))]
ing1 <- ingrowth[, .(Tnr, Enr, treeKey, species_name, time = 1L, dbh = dbh22,
                     status = "new", living = TRUE, n_ha = wzp_nha(dbh22))]

panel <- rbindlist(list(t0, t1, ing1), use.names = TRUE)
panel <- merge(panel, atrisk[, .(treeKey, fate)], by = "treeKey", all.x = TRUE)
panel[treeKey %in% ingrowth$treeKey, fate := "recruit"]
panel[, `:=`(year_t0 = PARAMS$year_t0, year_t1 = PARAMS$year_t1)]
panel[, year := fifelse(time == 0L, PARAMS$year_t0, PARAMS$year_t1)]
panel[, mort := status == "dead"]
panel[, is_recruit := fate == "recruit"]
n_dup <- nrow(panel) - uniqueN(panel, by = c("treeKey", "time"))
panel <- unique(panel, by = c("treeKey", "time"))
setorder(panel, Tnr, Enr, treeKey, time)
fwrite(panel, file.path(OUT, "panel.csv"))

#=#=#= corner table (scripts/40) =============================================
corners <- atrisk[, .(n_atrisk = .N,
                      n_harvest = sum(fate == "harvest"),
                      n_dead = sum(fate == "dead")),
                  by = .(Tnr, Enr)]
corners[, `:=`(harvest_share = n_harvest / n_atrisk,
               mort_share = n_dead / n_atrisk,
               siteName = paste(Tnr, Enr, sep = "_"),
               patchName = paste(Tnr, Enr, sep = "_"),
               year_t0 = PARAMS$year_t0, year_t1 = PARAMS$year_t1,
               period_length = PARAMS$period)]
corners <- merge(corners, ci_corners, by = c("Tnr", "Enr"), all.x = TRUE)
corners[is.na(ci), ci := FALSE]
fwrite(corners, file.path(OUT, "corners.csv"))

#=#=#= regression check: the v3 rule reproduces the scripts/40 panel ==========
reg40 <- tryCatch({
  old <- fread(file.path(OUT_DIR, "new_schema", "panel.csv"),
               select = c("treeKey", "time", "fate", "status", "dbh"))
  v3 <- ch[!is.na(fate_v3) & !drop_v3]
  v3[, treeKey := paste(Tnr, Enr, tree_number, sep = "_")]
  v3_t <- rbind(v3[, .(treeKey, time = 0L, fate = fate_v3, dbh = dbh12)],
                v3[, .(treeKey, time = 1L, fate = fate_v3,
                       dbh = fifelse(fate_v3 == "survivor", dbh22, dbh12))])
  ing_v3 <- b22[s22 == 0L & is.finite(dbh22) & dbh22 >= WZP_MIN_BHD,
                .(treeKey = paste(Tnr, Enr, tree_number, sep = "_"), time = 1L,
                  fate = "recruit", dbh = dbh22)]
  v3_t <- unique(rbind(v3_t, ing_v3), by = c("treeKey", "time"))
  m <- merge(old, v3_t, by = c("treeKey", "time"), all = TRUE, suffixes = c(".40", ".v3"))
  sprintf("scripts/40 panel %d rows, v3 rule here %d rows, in both %d; fate differs %d, dbh differs %d",
          nrow(old), nrow(v3_t), m[!is.na(fate.40) & !is.na(fate.v3), .N],
          m[!is.na(fate.40) & !is.na(fate.v3) & fate.40 != fate.v3, .N],
          m[!is.na(fate.40) & !is.na(fate.v3) & abs(dbh.40 - dbh.v3) > 1e-9, .N])
}, error = function(e) paste("not run:", conditionMessage(e)))
message("[regression] ", reg40)

#=#=#= IMPORT.md ===============================================================
tab <- function(d) paste(capture.output(print(d)), collapse = "\n")
fate_x <- dcast(ch[, .N, by = .(fate_v3 = fifelse(is.na(fate_v3), "(not at risk)", fate_v3),
                               fate_v5 = fifelse(is.na(fate), "(not at risk)", fate))],
                fate_v3 ~ fate_v5, value.var = "N", fill = 0)
src_x <- ch[!is.na(fate), .N, by = .(ci, fate_src, fate)][order(ci, fate_src, fate)]
ci_x <- dcast(ch[ci == TRUE, .N, by = .(s17 = fifelse(is.na(s17), "NA", as.character(s17)),
                                        s22 = fifelse(is.na(s22), "NA", as.character(s22)))],
              s17 ~ s22, value.var = "N", fill = 0)
commit <- tryCatch(system2("git", c("describe", "--always", "--dirty"), stdout = TRUE),
                   error = function(e) "unknown")
md <- c(
  "# outputs/new_schema_chained",
  "",
  sprintf("Written by `scripts/91_import_new_schema_chained.R` at commit `%s` on %s.",
          commit, format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "The 2012 -> 2022 tree panel of scripts/40 with (1) freshly dead trees (BWI2022",
  "`damage_dead`) as natural deaths and (2) fates chained through CI2017. All corners of",
  "the archive (the study-corner selection happens in scripts/92).",
  "",
  sprintf("- corners with CI2017 tree records: %d; CI2017 records flagged freshly dead: %d", nrow(ci_corners), c17[fd17 == TRUE, .N]),
  sprintf("- 2012 trees >= 7 cm: %d; at risk v3 rule %d -> v5 %d", nrow(b12), ch[!is.na(fate_v3) & !drop_v3, .N], nrow(atrisk)),
  sprintf("- freshly dead at the 2012 survey (removed from the 2012 stand): %d", ch[fd12 %in% TRUE, .N]),
  sprintf("- 2022 code 1 + freshly dead (now natural deaths): %d", ch[fate_src %chin% "bwi2022_fresh_dead", .N]),
  sprintf("- lost by 2017 at CI2017 corners (fate from CI2017): %d (dead %d, removed %d, untraceable %d)",
          ch[fate_src %chin% "ci2017", .N],
          ch[fate_src %chin% "ci2017" & fate == "dead", .N], ch[fate_src %chin% "ci2017" & fate == "harvest", .N],
          ch[fate_src %chin% "ci2017" & fate == "untraceable", .N]),
  sprintf("- ingrowth v3 %d -> v5 %d (new 2022 trees freshly dead dropped: %d; CI2017 ingrowth alive in 2022 added: %d)",
          n_ing_v3, nrow(ingrowth), b22[s22 == 0L & is.finite(dbh22) & dbh22 >= WZP_MIN_BHD & fd22, .N], nrow(ing17)),
  sprintf("- duplicate tree x time rows removed: %d", n_dup),
  sprintf("- regression check (v3 rule here vs outputs/new_schema/panel.csv): %s", reg40),
  "",
  "## v3 fate (rows) x v5 fate (columns), every 2012 tree >= 7 cm", "", "```", tab(fate_x), "```",
  "",
  "## v5 fate by source", "", "```", tab(src_x), "```",
  "",
  "## CI2017 corners: CI2017 code (rows) x BWI2022 code (columns), 2012 trees >= 7 cm", "", "```", tab(ci_x), "```",
  "",
  "## PARAMS", "", "```", capture.output(str(PARAMS, give.attr = FALSE)), "```")
writeLines(md, file.path(OUT, "IMPORT.md"))
cat("\nWrote", OUT, "\n")
