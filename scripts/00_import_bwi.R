#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-bwi | scripts/00_import_bwi.R
##
## One-time import: Thuenen ACCESS2003 databases -> data/bwi_csv/*.csv.
## Everything downstream reads the CSVs, so mdbtools is only needed here.
##
## Prerequisites (see README):
##   data/raw/bwi20150320_alle_daten{2002,2012,2002-2012}/*.mdb
##   data/raw/metadaten/bwi20150320_alle_metadaten.mdb
##   data/raw/inspire/Tnr_INSPIRE_Poly1000.shp
##   brew install mdbtools
##
## Run from the project root:
##   Rscript scripts/00_import_bwi.R          # incremental (skips fresh CSVs)
##   Rscript scripts/00_import_bwi.R --force  # re-export everything
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

suppressPackageStartupMessages(library(data.table))
source("R/config.R")
source("R/bwi_io.R")

FORCE <- "--force" %in% commandArgs(trailingOnly = TRUE)

## Which tables to pull from which database. Names are resolved
## case-insensitively against the real table list, because the Thuenen export
## is inconsistent about capitalisation (x_Av vs x_ba vs x_Wald).
IMPORTS <- list(
  ## --- BWI2 -> BWI3 change database: the demographic core -----------------
  list(mdb = MDB_CHANGE, tables = c(
    "b0_ecke",          # permanent corner info incl. the Use23 change flag
    "b23_baeume_a_s",   # tree state at the START of the period (BWI2)
    "b23_baeume_m_s",   # mid-period state + Pk fate code + *Diff columns
    "b23_ecke_w"        # forest-cover change over the period
  )),
  ## --- BWI3 (2012) --------------------------------------------------------
  list(mdb = MDB_2012, tables = c(
    "b3_baeume",        # tree state at the END of the period
    "b3_ecke_w",        # corner attributes (Wa, Eg, Ne, Begehbar, …)
    "b3_jung",          # regeneration, fixed-radius r=1 m / r=2 m circles
    "b3_bestock",       # stand-level attributes
    "b3_tab_w",         # per-Trakt field date (Datum_feld) -> period length
    "b0_tab",           # per-Trakt permanent info incl. HoeheNN (elevation)
    "b3_wzp"            # field record incl. Azi [gon] / Hori [cm] -> stem maps
  )),
  ## --- BWI2 (2002), for cross-checks --------------------------------------
  list(mdb = MDB_2002, tables = c(
    "b2_baeume",
    "b2_ecke_w",
    "b2_jung",
    "b2_bestock",
    "b2_tab_w",         # per-Trakt field date at the start of the period
    "b2_wzp"            # field record incl. Azi [gon] / Hori [cm] -> stem maps
  )),
  ## --- Metadata: data dictionary + the code tables we interpret -----------
  list(mdb = MDB_META, tables = c(
    "z_tab_phys", "z_col_phys",   # data dictionary
    "x_ba",                       # tree species  (Gattung + Art -> binomial)
    "x_pk",                       # Probebaum fate over the period
    "x_Av",                       # Aufnahmeverfahren (WZP vs sample circle)
    "x_Bs",                       # Bestandesschicht
    "x_bl",                       # Bundesland
    "x_Wald",                     # Waldkategorie
    "x_eg",                       # Eigentumsart
    "x_ne",                       # Nutzungseinschraenkung
    "x_Use23",                    # BWI2->BWI3 change-usability reasons
    "x_gr",                       # regeneration size class
    "x_begehbar"                  # accessibility
  ))
)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

dir.create(BWI_CSV_DIR, recursive = TRUE, showWarnings = FALSE)
missing_db <- character(0)
n_ok <- 0L; n_skip <- 0L

for (job in IMPORTS) {
  if (!file.exists(job$mdb)) {
    missing_db <- c(missing_db, job$mdb)
    next
  }
  message("\n== ", basename(job$mdb), " ==")
  available <- mdb_tables(job$mdb)
  lookup    <- setNames(available, tolower(available))

  for (want in job$tables) {
    actual <- lookup[[tolower(want)]]
    if (is.null(actual)) {
      warning("table '", want, "' not present in ", basename(job$mdb),
              call. = FALSE, immediate. = TRUE)
      n_skip <- n_skip + 1L
      next
    }
    ## Always cache under the *requested* (canonical) name so downstream code
    ## does not have to know about the source database's capitalisation.
    mdb_cache(job$mdb, actual, name = want, force = FORCE)
    n_ok <- n_ok + 1L
  }
}

if (length(missing_db)) {
  stop("\nMissing Access database(s):\n  ", paste(missing_db, collapse = "\n  "),
       "\n\nDownload them from\n",
       "  https://bwi.info/Download/de/BWI-Basisdaten/ACCESS2003/\n",
       "and unzip into data/raw/ (see README).", call. = FALSE)
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Derived lookups written alongside the raw exports
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

message("\n== derived ==")

## Species binomials from x_ba.
sp <- bwi_species_table()
fwrite(sp, file.path(BWI_CSV_DIR, "species_lookup.csv"))
message(sprintf("  %-22s %8d rows -> species_lookup.csv", "x_ba (binomials)",
                nrow(sp)))

## Anonymised tract coordinates (WGS84 lon/lat from the 1 km INSPIRE cells).
if (file.exists(INSPIRE_SHP)) {
  co <- bwi_trakt_coords()
  fwrite(co, file.path(BWI_CSV_DIR, "trakt_coords.csv"))
  message(sprintf("  %-22s %8d rows -> trakt_coords.csv", "INSPIRE centroids",
                  nrow(co)))
} else {
  warning("INSPIRE shapefile missing at ", INSPIRE_SHP,
          " -- scripts/03_attach_environment.R will not run.",
          call. = FALSE, immediate. = TRUE)
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Sanity-check the sampling design the whole pipeline rests on.
##
## The BWI stores the representation number as N_ha = K / G, where K is a
## per-tree count factor. K equals the nominal Zaehlfaktor (4 m^2/ha) except
## near a forest or stand edge, where Thuenen apply a Grenzkreis correction.
## So we do NOT expect an exact match with 4/G -- we expect the bulk of trees
## to sit exactly on it and a small edge-corrected minority to deviate. If that
## picture ever changes, the WZP_BAF assumption in R/config.R needs revisiting.
##
## Note also: N_ha is the WITHIN-STAND density and is 0 for trees flagged as not
## belonging to the stand (Bz = 0); N_ha_Hb is the across-stand density on
## wooded land. The pipeline uses N_ha (see R/wzp.R).
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

chk <- bwi("b3_baeume", select = c("Av", "G", "N_ha", "N_ha_Hb", "Bz"))
chk <- chk[Av == AV_WZP & is.finite(G) & G > 0 & Bz == 1L & N_ha > 0]
if (nrow(chk)) {
  ratio    <- chk$N_ha_Hb / (WZP_BAF / chk$G)
  frac_baf <- mean(abs(ratio - 1) < 1e-6)
  message(sprintf(
    "\n[check] WZP design on %d stand trees (BWI3):", nrow(chk)))
  message(sprintf(
    "        N_ha_Hb == %.0f/G exactly for %.1f%% of trees (rest: Grenzkreis edge correction) %s",
    WZP_BAF, 100 * frac_baf, if (frac_baf > 0.9) "OK" else "*** UNEXPECTED ***"))
  message(sprintf(
    "        median edge-corrected count factor K = %.3f m2/ha",
    median((WZP_BAF * ratio)[abs(ratio - 1) >= 1e-6])))
  if (frac_baf <= 0.9)
    warning("Fewer than 90 % of trees sit on the nominal Zaehlfaktor -- ",
            "verify WZP_BAF in R/config.R.", call. = FALSE)
}

message(sprintf("\nDone. %d tables exported (%d skipped) to %s/",
                n_ok, n_skip, BWI_CSV_DIR))
