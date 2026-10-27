#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Project: FINN-bwi
## Shared configuration: paths, BWI sampling-design constants, code semantics.
## Sourced by every script. Author: Yannek Kaeber <y.kaeber@posteo.de>
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

## ---- Paths --------------------------------------------------------------
## BWI_DIR holds the downloaded Thuenen "erweiterte Grunddaten" (ACCESS2003).
## Layout after scripts/00_import_bwi.R has run:
##   data/raw/bwi20150320_alle_daten2002/…2002.mdb          (BWI2, 2002)
##   data/raw/bwi20150320_alle_daten2012/…2012.mdb          (BWI3, 2012)
##   data/raw/bwi20150320_alle_daten2002-2012/…2002-2012.mdb(change / period)
##   data/raw/metadaten/…metadaten.mdb                      (dictionary + codes)
##   data/raw/inspire/Tnr_INSPIRE_Poly1000.shp              (anonymised coords)
##   data/bwi_csv/*.csv                                     (import output)
BWI_DIR     <- Sys.getenv("FINN_BWI_DIR",     unset = "data")
BWI_RAW_DIR <- file.path(BWI_DIR, "raw")
BWI_CSV_DIR <- Sys.getenv("FINN_BWI_CSV_DIR", unset = file.path(BWI_DIR, "bwi_csv"))

MDB_2002      <- file.path(BWI_RAW_DIR, "bwi20150320_alle_daten2002",
                           "bwi20150320_alle_daten2002.mdb")
MDB_2012      <- file.path(BWI_RAW_DIR, "bwi20150320_alle_daten2012",
                           "bwi20150320_alle_daten2012.mdb")
MDB_CHANGE    <- file.path(BWI_RAW_DIR, "bwi20150320_alle_daten2002-2012",
                           "bwi20150320_alle_daten2002-2012.mdb")
MDB_META      <- file.path(BWI_RAW_DIR, "metadaten",
                           "bwi20150320_alle_metadaten.mdb")
INSPIRE_SHP   <- file.path(BWI_RAW_DIR, "inspire", "Tnr_INSPIRE_Poly1000.shp")

## WorldClim 2.1 30-arcsec bundle (bio + elev). Same layout FINN-fia expects.
WORLDCLIM_DIR <- Sys.getenv("FINN_WORLDCLIM_DIR",
                            unset = file.path("..", "worldclim"))

OUT_DIR  <- "outputs"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---- Inventory calendar --------------------------------------------------
## BWI1 1986-88 (old Laender only, different design -> not used by default)
## BWI2 2001-02, BWI3 2011-12. Nominal remeasurement interval = 10 years; the
## per-tree truth is `VegJahreBaum` in b23_baeume_m_s (vegetation years).
BWI_YEARS         <- c(b2 = 2002L, b3 = 2012L)
BWI_PERIOD_YEARS  <- 10L

## ---- Sampling design constants ------------------------------------------
## Winkelzaehlprobe (angle-count / Bitterlich) with Zaehlfaktor 4 m^2/ha, for
## all stand trees with BHD >= 7 cm. Each tallied tree represents exactly
## BAF m^2/ha of basal area and N_ha = BAF / G stems/ha.
## Verified against the data: b3_baeume$N_ha == 4 / b3_baeume$G exactly.
WZP_BAF     <- 4.0          # m^2/ha, Zaehlfaktor ZF4
WZP_MIN_BHD <- 7.0          # cm, inclusion threshold for the WZP

## Fixed-radius circles used alongside the WZP (see code table x_Av):
##   r = 1.00 m  regeneration (b*_jung, Bnr 1xxx)
##   r = 2.00 m  regeneration (b*_jung, Bnr 2xxx)
##   r = 1.75 m / 4.00 m  additional sample circles in some inventories
CIRCLE_RADII_M <- c(r100 = 1.00, r175 = 1.75, r200 = 2.00, r400 = 4.00)
circle_area_ha <- function(r_m) pi * r_m^2 / 10000

## Trakt geometry: a Trakt is a square with 4 corners (Ecken) 150 m apart,
## numbered 1=A, 2=B, 3=C, 4=D, laid out on a 4x4 km national grid (denser in
## Laender with a Verdichtung, see column VBl).
TRAKT_SIDE_M  <- 150
TRAKT_N_ECKEN <- 4L

## Nominal FINN patch area (ha). The WZP has no fixed area, so we expand each
## tallied tree's N_ha to a nominal patch of this size (see R/wzp.R). 0.1 ha
## matches FINN's default patch_size and the classic gap-model scale.
## A tree is lost to integer rounding when N_ha * A < 0.5, i.e. when
## G > 8 * A m^2  <=>  BHD > 200 * sqrt(8 * A / pi) cm.
##   A = 0.10 ha -> BHD > 101 cm ; A = 0.25 ha -> BHD > 160 cm
## R/wzp.R uses largest-remainder allocation so species totals survive anyway.
PATCH_SIZE_HA <- 0.1

## ---- Code semantics ------------------------------------------------------
## Pk = Probebaum-Kennung (b23_baeume_m_s), the fate of a BWI2 sample tree over
## the BWI2->BWI3 period. Source: code table x_pk in the metadata database.
##   0  neuer Probebaum          -> ingrowth into the WZP at BWI3
##   1  wiederh. aufgen.         -> survivor
##   2  selektiv entnommen       -> HARVEST (thinning)
##   3  bei Kahlschlag entnommen -> HARVEST (clearfell)
##   4  am Ort verblieben        -> natural mortality (fallen / not standing)
##   5  abgestorben              -> natural mortality (standing dead > 12 mo)
##   9  nicht auffindbar         -> ambiguous (untraceable; often removal)
##   6,7,8,10,11 out of sample / invalid
PK_NEW       <- 0L
PK_SURVIVOR  <- 1L
PK_HARVEST   <- c(2L, 3L)
PK_MORT_NAT  <- c(4L, 5L)
PK_AMBIGUOUS <- 9L
PK_INVALID   <- c(6L, 7L, 8L, 10L, 11L)

## Av = Aufnahmeverfahren (x_Av). 0 = WZP of the current inventory; 1/2/3 and
## 4/5/6 flag modelled values; >= 100 are the fixed-radius sample circles.
AV_WZP           <- 0L
AV_MODELLED      <- c(1L, 2L, 3L, 4L, 5L, 6L)
AV_CIRCLE        <- c(100L, 101L, 102L, 175L, 200L, 400L)

## Bs = Bestandesschicht (x_Bs). 0 Plenterwald, 1 Hauptbestand, 2 Unterstand,
## 3 Oberstand, 4 Verjuengung unter Schirm, 9 liegend.
BS_LYING <- 9L

## Wa = Waldkategorie / Waldentscheid (x_Wald). Accessible forest (begehbarer
## Holzboden) is what the WZP is measured on.
## Use23 (b0_ecke) is Thuenen's own flag for whether an Ecke may be used in a
## BWI2->BWI3 change analysis; 100 = usable. Always prefer it over hand-rolled
## rules.
USE23_OK <- 100L

## ---- Bundesland lookup (x_bl) -------------------------------------------
BL_LOOKUP <- data.frame(
  Bl = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16),
  BL_NAME = c("Schleswig-Holstein", "Hamburg", "Niedersachsen", "Bremen",
              "Nordrhein-Westfalen", "Hessen", "Rheinland-Pfalz",
              "Baden-Wuerttemberg", "Bayern", "Saarland", "Berlin",
              "Brandenburg", "Mecklenburg-Vorpommern", "Sachsen",
              "Sachsen-Anhalt", "Thueringen"),
  stringsAsFactors = FALSE
)

## ---- Small helpers -------------------------------------------------------
## BWI stores BHD in mm and heights in dm; FINN wants cm and m.
mm2cm <- function(x) x / 10
dm2m  <- function(x) x / 10
## Basal area of a single stem, m^2, from BHD in cm.
bhd2ba <- function(bhd_cm) pi * (bhd_cm / 200)^2

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- Which env columns does a fitted model expect? -----------------------
## FINN records the covariates it was scaled on ON THE MODEL, as
## `m$env_scaling$variable` (written by fit() when env_autoscale = TRUE), and
## predict()/simulate()/ALE() all run apply_env_scaling() internally. Handing
## them a different set fails with
##   "env is missing environmental variable(s) the model was scaled on: ..."
## Downstream scripts must therefore NOT hardcode c("mat_c","map_mm","elev_m"):
## the hpc/configs.tsv matrix fits variants with `drought` and with `x`,`y`, so
## a hardcoded list silently pins every diagnostic to the baseline variant and
## then dies on any other. Always ask the model.
fit_env_vars <- function(m, fallback = c("mat_c", "map_mm", "elev_m")) {
  ## 1) env_scaling$variable is authoritative -- but it is only populated when
  ##    FINN did the scaling (env_autoscale != FALSE). When env is PRE-scaled in
  ##    scripts/04 and the model is fitted with env_autoscale = FALSE, this is
  ##    NULL, and the old fallback silently dropped any predictor beyond the
  ##    climate triple (e.g. `management`), breaking every downstream env rebuild.
  v <- tryCatch(m$env_scaling$variable, error = function(e) NULL)
  if (!is.null(v) && length(v)) return(as.character(v))
  ## 2) train_env carries the exact fitted env columns, in order.
  te <- tryCatch(names(m$train_env), error = function(e) NULL)
  if (!is.null(te) && length(te))
    return(setdiff(as.character(te), c("siteID", "year")))
  ## 3) last resort before the hardcoded default: the process formulas name
  ##    every env covariate the model actually uses.
  ff <- c("growth_formula", "mortality_formula",
          "regeneration_formula", "competition_formula")
  fv <- unique(unlist(lapply(ff, function(f)
    tryCatch(all.vars(m[[f]]), error = function(e) character(0)))))
  if (length(fv)) return(as.character(fv))
  fallback
}

## One row per site, holding each covariate CONSTANT for a long forward run.
## Long-horizon simulations extend far past the observed decade, so a per-year
## anomaly like `drought` has no calendar to be read off; averaging it over the
## observed years is the neutral choice and keeps the run at mean drought.
## This also fixes a real trap: the previous idiom
##   unique(env_s[siteID %in% keep, c("siteID", FIT_VARS), with = FALSE])
## returns ONE row per site only while every covariate is time-invariant. Add a
## per-year covariate and it returns up to one row per site-year, which then gets
## replicated across years into duplicated (siteID, year) keys.
site_env_constant <- function(env_dt, vars, years) {
  stopifnot(all(vars %in% names(env_dt)))
  e <- env_dt[, lapply(.SD, function(x) mean(x, na.rm = TRUE)),
              by = siteID, .SDcols = vars]
  out <- data.table::rbindlist(
    lapply(years, function(k) data.table::copy(e)[, year := k]))
  data.table::setcolorder(out, c("siteID", "year", vars))
  out[]
}

## Subset an env table to exactly those columns, erroring early (with a useful
## message) rather than deep inside torch.
env_for_model <- function(env_dt, m, extra = character()) {
  v <- fit_env_vars(m)
  miss <- setdiff(v, names(env_dt))
  if (length(miss))
    stop("env_dt lacks column(s) the model was fitted on: ",
         paste(miss, collapse = ", "),
         ".\nRe-run scripts/04_prepare_finn_inputs.R with a matching ",
         "FINN_BWI_ENV / FINN_BWI_DROUGHT.", call. = FALSE)
  env_dt[, unique(c("siteID", "year", v, extra)), with = FALSE]
}
