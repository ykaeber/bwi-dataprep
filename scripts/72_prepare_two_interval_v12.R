#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## bwi-dataprep | scripts/72_prepare_two_interval_v12.R
##
## two-interval-v12 = two-interval-v8 without env row `year 0`. Nothing else
## changes. v12 is the set the final model of the FINN-bwi study is fitted on.
##
## RELATION TO THE FINN-bwi STUDY. This is the study's
## scripts/148_prepare_two_interval_v12_no_year0.R; only this header, the
## sourced R/config.R (OUT_DIR) and the script name in PREPARATION.md differ.
## Run on the study's v8 it reproduces every file of the study's v12 byte for
## byte (outputs/prepared/two-interval-v12.sha256). Run on this repository's
## v8 (scripts/71), the same files differ only where v8 differs by design:
## split_blocks.csv (v8's 12,388 tracts) and the study's two v5 comparison
## files (v5_v8_state.csv, v5_v8_identity.csv), which are not written. Check
## with
##   (cd outputs/prepared/two-interval-v12 && shasum -a 256 -c ../two-interval-v12.sha256)
## The study's versions v9-v11 (alternative regeneration targets) are not part
## of this repository; v12 builds on v8 directly.
##
## WHY    Every two-interval version (v1-v11) lays the annual env grid out as
##        year = 0:10 (calendar t0..t0+10; scripts/71, the env grid) with the census
##        at obs year 10. Year 0 was meant as the initial state, but FINN has no
##        initial-state env row: the initial state is init_cohort_dt, and every
##        env row drives one annual transition. fit() maps obs year 10 to step
##        11 (year_sequence from the factor levels 0..10, FINN
##        R/finn-class.R:1894), and each step applies growth, mortality and
##        regeneration before the stand is aggregated (:946-1197). So the
##        scored dbh / ba / trees were the state after 11 transitions for a
##        10-year interval, the first driven by the t0 climate. The rate
##        windows (growth, mortality, recruitment) were already the correct 10
##        steps, calendar t0+1..t0+10. Verified with a synthetic run (a
##        certain-death pulse in env row 10 kills the stand at step 11, not 10).
##
## FIX    env_dt and env_unscaled_dt keep only years 1..10 (calendar
##        t0+1..t0+10). obs_dt / tree_obs_dt stay at year 10 with
##        period_length 10, so year_sequence = 10: the scored state is 10
##        transitions after the start census, and the rate window is steps
##        1..10 = the same calendar years as in v8. Calendar year of simulated
##        year k is now t0 + k (v8: t0 + k - 1); scripts reading v12 use
##        sim_year = obs_year, not obs_year + 1.
##
## SCALING (decision of this script)
##        env_scales.csv is v8's, unchanged, and the kept env values are v8's
##        values byte for byte. Recomputing the pooled mean / sd over years
##        1..10 would move every scaled value slightly and mix a second change
##        into the version; with v8's scales the only difference between v8 and
##        v12 is the dropped row, which is what the fix is.
##
## Inputs   outputs/prepared/two-interval-v8/ (every set, split_blocks.csv; scripts/71)
## Outputs  outputs/prepared/two-interval-v12/
##            PREPARATION.md
##            <set>/env_dt.csv, env_unscaled_dt.csv   v8 rows with year >= 1
##            every other file                        byte copy of v8
##          outputs/prepared/two-interval-v12.sha256  manifest of every file
##          (the manifest in git is the study's; a run overwrites it with its own)
## Helpers  R/config.R (OUT_DIR)
## RUN      Rscript scripts/72_prepare_two_interval_v12.R   (from the repo root, ~1 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))
source("R/config.R")

PARAMS <- list(
  version   = "two-interval-v12",
  base_dir  = file.path(OUT_DIR, "prepared", "two-interval-v8"),   # the base (scripts/71); never modified
  drop_year = 0L,                                   # the env row FINN would turn into an extra transition
  P         = 10L,                                  # interval length; obs year and period_length stay P
  env_files = c("env_dt.csv", "env_unscaled_dt.csv")  # the only files that change
)
BASE <- PARAMS$base_dir
OUT  <- file.path(OUT_DIR, "prepared", PARAMS$version)
stopifnot(dir.exists(BASE), !dir.exists(OUT))       # a version is never modified
sets <- list.dirs(BASE, full.names = FALSE, recursive = FALSE)
message("sets: ", paste(sets, collapse = ", "))

dir.create(OUT, recursive = TRUE)
top_files <- setdiff(list.files(BASE), c(sets, "PREPARATION.md"))
file.copy(file.path(BASE, top_files), OUT, copy.date = TRUE)

CHECKS <- list()
for (s in sets) {
  dir.create(file.path(OUT, s))
  files <- list.files(file.path(BASE, s))
  file.copy(file.path(BASE, s, setdiff(files, PARAMS$env_files)), file.path(OUT, s), copy.date = TRUE)

  obs  <- fread(file.path(BASE, s, "obs_dt.csv"))
  stopifnot(all(obs$year == PARAMS$P), all(obs$period_length == PARAMS$P))
  if (file.exists(file.path(BASE, s, "tree_obs_dt.csv"))) {
    tob <- fread(file.path(BASE, s, "tree_obs_dt.csv"))
    stopifnot(all(tob$year == PARAMS$P))
  }

  for (f in PARAMS$env_files) {
    env <- fread(file.path(BASE, s, f))
    stopifnot(setequal(unique(env$year), 0:PARAMS$P))
    new <- env[year != PARAMS$drop_year]
    stopifnot(setequal(unique(new$year), 1:PARAMS$P),
              nrow(new) == nrow(env) * PARAMS$P / (PARAMS$P + 1L),
              uniqueN(new$siteID) == uniqueN(env$siteID))
    fwrite(new, file.path(OUT, s, f))
    CHECKS[[length(CHECKS) + 1L]] <- data.table(
      set = s, file = f, rows_v8 = nrow(env), rows_v12 = nrow(new),
      sites = uniqueN(new$siteID),
      kept_rows_identical = isTRUE(all.equal(env[year != PARAMS$drop_year], fread(file.path(OUT, s, f)),
                                             check.attributes = FALSE)))
  }
}
CHECKS <- rbindlist(CHECKS)
print(CHECKS)
stopifnot(all(CHECKS$kept_rows_identical))
fwrite(CHECKS, file.path(OUT, "v8_v12_env_check.csv"))

## manifest, same format as two-interval-v8.sha256 (paths relative to the version folder)
all_files <- list.files(OUT, recursive = TRUE)
sha <- vapply(all_files, function(f) {
  out <- system2("shasum", c("-a", "256", shQuote(file.path(OUT, f))), stdout = TRUE)
  strsplit(out, " ")[[1]][1]
}, character(1))
writeLines(sprintf("%s  ./%s", sha, all_files), paste0(OUT, ".sha256"))

commit <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
prep <- c(
  "# two-interval-v12",
  "",
  sprintf("Written by `scripts/72_prepare_two_interval_v12.R` (bwi-dataprep) at commit `%s` on %s.", commit, Sys.Date()),
  "",
  "v12 = v8 without env row `year 0`. Nothing else changes.",
  "",
  "**Why.** Every earlier two-interval version laid the env grid out as year 0..10 (calendar t0..t0+10)",
  "with the census at obs year 10. Year 0 was meant as the initial state, but FINN has no initial-state",
  "env row: every env row drives one annual transition, and the initial state is `init_cohort_dt`.",
  "`fit()` therefore scored dbh / ba / trees after 11 transitions for a 10-year interval; the rate",
  "windows were already the correct 10 years. In v12 the census is at step 10.",
  "",
  "**Calendar.** env year k = calendar t0 + k, k = 1..10. Simulated year k = calendar t0 + k.",
  "Scripts reading v12 take the scored state at `sim_year = obs_year` (v8: `obs_year + 1`).",
  "",
  "**Scaling.** `env_scales.csv` is v8's; every kept env value is v8's value unchanged.",
  "",
  "**Files.** `env_dt.csv` and `env_unscaled_dt.csv` lose the year-0 rows; every other file is a byte",
  "copy of v8. Checks: `v8_v12_env_check.csv`. Manifest: `../two-interval-v12.sha256`.",
  "",
  "## PARAMS",
  "",
  sprintf("- `%s` = `%s`", names(PARAMS), vapply(PARAMS, function(x) paste(x, collapse = ", "), character(1)))
)
writeLines(prep, file.path(OUT, "PREPARATION.md"))
message("done: ", OUT)
