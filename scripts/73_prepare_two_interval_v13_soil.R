#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## bwi-dataprep | scripts/73_prepare_two_interval_v13_soil.R
##
## two-interval-v13-soil = two-interval-v12 plus static soil covariates per
## corner. Nothing else changes.
##
## RELATION TO THE FINN-bwi STUDY. This is the study's
## scripts/198_prepare_two_interval_v13_soil.R; only this header, the sourced
## R/config.R (OUT_DIR, BWI_RAW_DIR) and the script name in PREPARATION.md
## differ. Run on the study's v12 it reproduces every file of the study's
## v13-soil byte for byte (outputs/prepared/two-interval-v13-soil.sha256); run
## on this repository's v12 (scripts/72) the same files differ as v12 does
## (split_blocks.csv, the two v5 comparison files). Check with
##   (cd outputs/prepared/two-interval-v13-soil && shasum -a 256 -c ../two-interval-v13-soil.sha256)
## The soil database is a separate download: data/README.md, "Soil profiles".
##
## WHY    The env covariates carry no site soil information: the drought index
##        is a DWD soil-moisture anomaly computed for grass on one sandy loam
##        everywhere. The bare-ground equilibria of the FINN-bwi study end
##        beech-dominated even in climates where the BWI's natural forest
##        community is oak, pine or alder; in Central Europe beech loses to oak,
##        hornbeam, alder and ash on waterlogged soils and to pine and oak on
##        dry, poor sands (Ellenberg & Leuschner 2010). The BWI 2012
##        environmental database gives a soil profile for the 4 x 4 km corners.
##
## SOURCE Benning, R., Petzold, R., Danigel, J., Mayer, S., Profft, I., Steinicke, C.,
##        Hafner, S., Ahrends, B., Heinz, F., Puhlmann, H., von Wilpert, K.,
##        Amberger, H., Gauer, J., Pieper, A., Wirner, M., Wallor, E., Riek, W.,
##        Janott, M., Kölling, C., Mette, T., 2019. BWI 2012 Umweltdatenbank
##        Bodenprofile. https://doi.org/10.3220/DATA20190625100522
##        (local copy data/raw/openagrar_derivate_00022249, read with mdbtools;
##        docs/bwi_data_notes.md section 13; download: data/README.md).
##
## VARIABLES (one value per corner, the same in every env row of a site)
##        soil_nfk           available water capacity to 100 cm (mm), table
##                           04_nFk_Dehner_et_al_2015_h0 column nfk_100
##        soil_stagnic       perched-water (stagnic) degree 0-6 from the
##                           harmonised water-regime class WASSER (S0 .. S6; "S"
##                           very deep = 0.25, S0 deep = 0.5); 0 on terrestrial
##                           and groundwater classes
##        soil_ground        groundwater degree 0-6 from WASSER (G1 .. G6, G0 =
##                           0.5, "G2/3" = 2.5); 0 on terrestrial and stagnic classes
##        soil_wet           waterlogging degree 0-6 = soil_stagnic + soil_ground
##                           (a class is either S or G; a range such as "S3-S6"
##                           or "T3-S1" takes the mean of its two ends)
##        soil_nutrient_poor nutrient class NAEHR, 1 eutrophic .. 6 dystrophic
##        A corner with two mapping units gets their area-weighted mean
##        (SOEH_ANT, normalised to 1; equal weights when the share is missing).
##
## GAPS   The database covers the 4 x 4 km grid; densification corners have no
##        profile, and some profiled units lack a class. Every missing value is
##        filled per variable with the mean over profiled corners of the same
##        natural forest community (potNatWg, outputs/corner_metadata.csv) in
##        the same state (Bl) when at least min_group corners hold it, else of
##        the same community nationally, else of the state. soil_dt.csv records
##        the source of every value.
##
## SCALING env_scales.csv keeps v12's rows unchanged; the soil variables get
##        their mean / sd over the env rows of the national set (as v8 scaled
##        its variables over the pooled env rows).
##
## Inputs   outputs/prepared/two-interval-v12/ (every set; scripts/72);
##          data/raw/openagrar_derivate_00022249/bwi2012_umweltdatenbank_bodenprofile.mdb;
##          outputs/corner_metadata.csv (Tnr, Enr, Bl, potNatWg; scripts/01)
## Helpers  R/config.R (OUT_DIR, BWI_RAW_DIR); mdbtools (mdb-export) on the PATH
## Outputs  outputs/prepared/two-interval-v13-soil/
##            PREPARATION.md, soil_corners.csv (every corner of the national set)
##            <set>/env_dt.csv, env_unscaled_dt.csv, env_scales.csv   v12 + soil columns
##            <set>/soil_dt.csv                                       siteID, corner, unscaled soil values, sources
##            every other file                                        byte copy of v12
##          outputs/prepared/two-interval-v13-soil.sha256  manifest of every file
##          (the manifest in git is the study's; a run overwrites it with its own)
## RUN      Rscript scripts/73_prepare_two_interval_v13_soil.R   (from the repo root, ~2 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))
source("R/config.R")

PARAMS <- list(
  version    = "two-interval-v13-soil",
  base_dir   = file.path(OUT_DIR, "prepared", "two-interval-v12"),   # the base (scripts/72); never modified
  soil_mdb   = file.path(BWI_RAW_DIR, "openagrar_derivate_00022249", "bwi2012_umweltdatenbank_bodenprofile.mdb"),
  corner_meta = file.path(OUT_DIR, "corner_metadata.csv"),           # Tnr, Enr, Bl, potNatWg for the gap filling (scripts/01)
  nfk_col    = "nfk_100",                             # rooting zone of a mature stand taken as 100 cm (the 160 cm value adds deep subsoil few roots reach)
  soil_vars  = c("soil_nfk", "soil_stagnic", "soil_ground", "soil_wet", "soil_nutrient_poor"),
  min_group  = 10L,                                   # profiled corners a community x state (or community) group needs to fill a gap
  scale_set  = "national",                            # the set whose env rows define the soil scales
  env_files  = c("env_dt.csv", "env_unscaled_dt.csv", "env_scales.csv")   # the only v12 files that change
)
BASE <- PARAMS$base_dir
OUT  <- file.path(OUT_DIR, "prepared", PARAMS$version)
stopifnot(dir.exists(BASE), !dir.exists(OUT), file.exists(PARAMS$soil_mdb))   # a version is never modified
sets <- list.dirs(BASE, full.names = FALSE, recursive = FALSE)
message("sets: ", paste(sets, collapse = ", "))
mdb <- function(tab) fread(cmd = paste("mdb-export", shQuote(PARAMS$soil_mdb), shQuote(tab)), na.strings = c("NA", "-9999", ""))

## ---- 1. soil per mapping unit -> per corner -------------------------------------------
units <- mdb("02_KARTIEREINHEITEN")[, .(TNR, ENR, SOEH_NR, SOEH_ANT, WASSER, NAEHR)]
nfk   <- mdb("04_nFk_Dehner_et_al_2015_h0")[, c("TNR", "ENR", "SOEH_NR", PARAMS$nfk_col), with = FALSE]
setnames(nfk, PARAMS$nfk_col, "nfk")
units <- merge(units, nfk, by = c("TNR", "ENR", "SOEH_NR"), all.x = TRUE)

## water-regime class -> stagnic and groundwater degree; a range "A-B" is the mean of its ends
wr_end <- function(code) {
  lvl <- suppressWarnings(as.numeric(sub("^[TSG]", "", sub("2/3", "2.5", code))))
  fam <- substr(code, 1, 1)
  lvl[fam == "S" & code == "S"]  <- 0.25                 # very deep perched water
  lvl[fam == "S" & code == "S0"] <- 0.5                  # deep perched water
  lvl[fam == "G" & code == "G0"] <- 0.5                  # groundwater fresh, depth not given
  list(stagnic = ifelse(fam == "S", lvl, 0), ground = ifelse(fam == "G", lvl, 0))
}
wr_parse <- function(w) {
  out <- data.table(stagnic = NA_real_, ground = NA_real_)[rep(1L, length(w))]
  ok <- !is.na(w) & w != "NA"
  ends <- tstrsplit(w[ok], "-", fixed = TRUE, fill = NA_character_)
  a <- wr_end(ends[[1]]); b <- if (length(ends) > 1) wr_end(ifelse(is.na(ends[[2]]), ends[[1]], ends[[2]])) else a
  out[ok, `:=`(stagnic = (a$stagnic + b$stagnic) / 2, ground = (a$ground + b$ground) / 2)]
  out
}
units[, c("stagnic", "ground") := wr_parse(WASSER)]
stopifnot(units[!is.na(WASSER) & WASSER != "NA", all(is.finite(stagnic) & is.finite(ground))],
          units[, all(is.na(stagnic) | (stagnic >= 0 & stagnic <= 6 & ground >= 0 & ground <= 6 & (stagnic == 0 | ground == 0 | grepl("-", WASSER))))])
units[, w := fifelse(is.na(SOEH_ANT), 1, SOEH_ANT)]
wmean <- function(x, w) { ok <- is.finite(x); if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok]) }
prof <- units[, .(soil_nfk = wmean(nfk, w), soil_stagnic = wmean(stagnic, w), soil_ground = wmean(ground, w),
                  soil_nutrient_poor = wmean(NAEHR, w), n_units = .N), by = .(TNR, ENR)]
prof[, soil_wet := soil_stagnic + soil_ground]
prof[, corner := paste0(TNR, "_", ENR)]

## ---- 2. every corner of the national set, gaps filled by community x state ---------------
meta <- fread(PARAMS$corner_meta)[, .(corner = paste0(Tnr, "_", Enr), Bl, potNatWg)]
meta <- unique(meta, by = "corner")
nat_sites <- unique(fread(file.path(BASE, "national", "siteID_dt.csv"))[, .(siteID, corner)])   # siteID_dt has a row per patch and year
stopifnot(!anyDuplicated(nat_sites$siteID))
SC <- merge(data.table(corner = unique(nat_sites$corner)), meta, by = "corner", all.x = TRUE)
stopifnot(!anyNA(SC$Bl))
SC <- merge(SC, prof[, c("corner", PARAMS$soil_vars, "n_units"), with = FALSE], by = "corner", all.x = TRUE)
SC[, has_profile := !is.na(n_units)]
for (v in PARAMS$soil_vars) {
  src <- paste0(v, "_source")
  SC[, (src) := fifelse(is.finite(get(v)), "profile", NA_character_)]
  g_cs <- SC[is.finite(get(v)), .(m = mean(get(v)), n = .N), by = .(potNatWg, Bl)][n >= PARAMS$min_group]
  g_c  <- SC[is.finite(get(v)), .(m = mean(get(v)), n = .N), by = .(potNatWg)][n >= PARAMS$min_group]
  g_s  <- SC[is.finite(get(v)), .(m = mean(get(v))), by = .(Bl)]
  SC[g_cs, on = .(potNatWg, Bl), `:=`(fill_cs = i.m)]
  SC[g_c,  on = .(potNatWg),     `:=`(fill_c = i.m)]
  SC[g_s,  on = .(Bl),           `:=`(fill_s = i.m)]
  SC[!is.finite(get(v)) & is.finite(fill_cs), c(v, src) := list(fill_cs, "community x state")]
  SC[!is.finite(get(v)) & is.finite(fill_c),  c(v, src) := list(fill_c, "community")]
  SC[!is.finite(get(v)) & is.finite(fill_s),  c(v, src) := list(fill_s, "state")]
  SC[, c("fill_cs", "fill_c", "fill_s") := NULL]
}
stopifnot(SC[, all(vapply(.SD, function(x) all(is.finite(x)), logical(1))), .SDcols = PARAMS$soil_vars])
COVER <- rbindlist(lapply(PARAMS$soil_vars, function(v) SC[, .N, by = c(source = paste0(v, "_source"))][, `:=`(variable = v, share = N / sum(N))]))
print(dcast(COVER, variable ~ source, value.var = "share"))

dir.create(OUT, recursive = TRUE)
fwrite(SC, file.path(OUT, "soil_corners.csv"))
top_files <- setdiff(list.files(BASE), c(sets, "PREPARATION.md"))
invisible(file.copy(file.path(BASE, top_files), OUT, copy.date = TRUE))

## ---- 3. scales from the national env rows ------------------------------------------------
envu_nat <- merge(fread(file.path(BASE, PARAMS$scale_set, "env_unscaled_dt.csv"))[, .(siteID)],
                  merge(nat_sites[, .(siteID, corner)], SC[, c("corner", PARAMS$soil_vars), with = FALSE], by = "corner"), by = "siteID")
SOIL_SCALES <- data.table(variable = PARAMS$soil_vars,
                          mean = vapply(PARAMS$soil_vars, function(v) mean(envu_nat[[v]]), numeric(1)),
                          sd   = vapply(PARAMS$soil_vars, function(v) sd(envu_nat[[v]]), numeric(1)))
stopifnot(all(SOIL_SCALES$sd > 0))
print(SOIL_SCALES)

## ---- 4. every set: v12 + soil columns ----------------------------------------------------
CHECKS <- list()
for (s in sets) {
  dir.create(file.path(OUT, s))
  files <- list.files(file.path(BASE, s))
  invisible(file.copy(file.path(BASE, s, setdiff(files, PARAMS$env_files)), file.path(OUT, s), copy.date = TRUE))
  sites <- fread(file.path(BASE, s, "siteID_dt.csv"))
  site_soil <- merge(unique(sites[, .(siteID, corner)]), SC[, c("corner", PARAMS$soil_vars, paste0(PARAMS$soil_vars, "_source")), with = FALSE],
                     by = "corner", all.x = TRUE)
  stopifnot(!anyNA(site_soil[, PARAMS$soil_vars, with = FALSE]), !anyDuplicated(site_soil$siteID))
  fwrite(site_soil[order(siteID)], file.path(OUT, s, "soil_dt.csv"))

  sc_old <- fread(file.path(BASE, s, "env_scales.csv"))
  stopifnot(!any(PARAMS$soil_vars %in% sc_old$variable))
  fwrite(rbind(sc_old, SOIL_SCALES), file.path(OUT, s, "env_scales.csv"))

  for (f in c("env_unscaled_dt.csv", "env_dt.csv")) {
    env <- fread(file.path(BASE, s, f))
    new <- copy(env)   # an update join keeps v12's row order
    new[site_soil, on = "siteID", (PARAMS$soil_vars) := mget(paste0("i.", PARAMS$soil_vars))]
    if (f == "env_dt.csv") for (v in PARAMS$soil_vars) new[, (v) := (get(v) - SOIL_SCALES[variable == v, mean]) / SOIL_SCALES[variable == v, sd]]
    stopifnot(nrow(new) == nrow(env), !anyNA(new[, PARAMS$soil_vars, with = FALSE]),
              new[, all(vapply(.SD, uniqueN, integer(1)) == 1L), by = siteID, .SDcols = PARAMS$soil_vars][, all(V1)])
    fwrite(new, file.path(OUT, s, f))
    back <- fread(file.path(OUT, s, f))
    CHECKS[[length(CHECKS) + 1L]] <- data.table(
      set = s, file = f, rows = nrow(back), sites = uniqueN(back$siteID),
      v12_columns_identical = isTRUE(all.equal(env, back[, names(env), with = FALSE], check.attributes = FALSE)),
      share_profile = mean(site_soil$soil_nfk_source == "profile"))
  }
}
CHECKS <- rbindlist(CHECKS)
print(CHECKS)
stopifnot(all(CHECKS$v12_columns_identical))
fwrite(CHECKS, file.path(OUT, "v12_v13_env_check.csv"))

## ---- 5. manifest and PREPARATION.md -------------------------------------------------------
all_files <- list.files(OUT, recursive = TRUE)
sha <- vapply(all_files, function(f) {
  out <- system2("shasum", c("-a", "256", shQuote(file.path(OUT, f))), stdout = TRUE)
  strsplit(out, " ")[[1]][1]
}, character(1))
writeLines(sprintf("%s  ./%s", sha, all_files), paste0(OUT, ".sha256"))

md_table <- function(x) {
  x <- copy(x); for (cc in names(x)) if (is.numeric(x[[cc]])) x[, (cc) := formatC(get(cc), digits = 3, format = "g")]
  c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
SUMM <- rbindlist(lapply(PARAMS$soil_vars, function(v) SC[, .(variable = v, profile_mean = mean(get(v)[get(paste0(v, "_source")) == "profile"]),
                                                              filled_mean = mean(get(v)[get(paste0(v, "_source")) != "profile"]),
                                                              q05 = quantile(get(v), 0.05), q50 = median(get(v)), q95 = quantile(get(v), 0.95))]))
commit <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
prep <- c(
  "# two-interval-v13-soil", "",
  sprintf("Written by `scripts/73_prepare_two_interval_v13_soil.R` (bwi-dataprep) at commit `%s` on %s.", commit, Sys.Date()), "",
  "v13-soil = v12 plus static soil covariates per corner. Nothing else changes.", "",
  "**Why.** No v12 covariate describes the site's soil (the drought index is a soil-moisture anomaly for grass on one sandy loam everywhere); the bare-ground equilibria of the FINN-bwi study end beech-dominated even where the BWI's natural forest community is oak, pine or alder.", "",
  "**Source.** Benning, R., Petzold, R., Danigel, J., Mayer, S., Profft, I., Steinicke, C., Hafner, S., Ahrends, B., Heinz, F., Puhlmann, H., von Wilpert, K., Amberger, H., Gauer, J., Pieper, A., Wirner, M., Wallor, E., Riek, W., Janott, M., Kölling, C., Mette, T., 2019. BWI 2012 Umweltdatenbank Bodenprofile. BWI 2012 Umweltdatenbank. https://doi.org/10.3220/DATA20190625100522", "",
  "**Variables** (static per corner, the same value in every env row): `soil_nfk` available water capacity to 100 cm (mm); `soil_stagnic` perched-water degree 0-6 and `soil_ground` groundwater degree 0-6 from the harmonised water-regime class (S0..S6 / G1..G6; S = 0.25, S0 = G0 = 0.5, G2/3 = 2.5; ranges take the mean of their ends); `soil_wet` = stagnic + groundwater degree; `soil_nutrient_poor` nutrient class 1 eutrophic .. 6 dystrophic. Two mapping units per corner: area-weighted mean.", "",
  sprintf("**Gaps.** The database covers the 4 x 4 km corners. Missing values are filled per variable with the mean of profiled corners of the same natural forest community (potNatWg) in the same state when >= %d corners hold it, else of the community nationally, else of the state; the source of every value is in `<set>/soil_dt.csv` and `soil_corners.csv`.", PARAMS$min_group), "",
  md_table(dcast(COVER, variable ~ source, value.var = "share", fill = 0)), "",
  "**Values** over the national corners (filled_mean: the mean of the filled values):", "",
  md_table(SUMM), "",
  "**Scaling.** `env_scales.csv` = v12's rows unchanged plus the soil rows (mean / sd over the national env rows):", "",
  md_table(SOIL_SCALES), "",
  "**Files.** `env_dt.csv`, `env_unscaled_dt.csv` = v12 plus the soil columns (v12 columns identical, `v12_v13_env_check.csv`); `env_scales.csv` = v12 plus the soil rows; new `soil_dt.csv` per set; every other file a byte copy of v12. Manifest: `../two-interval-v13-soil.sha256`.", "",
  "## PARAMS", "",
  sprintf("- `%s` = `%s`", names(PARAMS), vapply(PARAMS, function(x) paste(x, collapse = ", "), character(1)))
)
writeLines(prep, file.path(OUT, "PREPARATION.md"))
message("done: ", OUT)
