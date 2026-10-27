#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-bwi | scripts/03_attach_environment.R
##
## Attaches environmental covariates to the calibration corners.
##
## COORDINATE PRECISION
## Thuenen do not release exact tract positions -- the sample has to be
## protected. What is public is the 1 km INSPIRE (ETRS89-LAEA) grid cell
## containing the tract centre, so we use the cell centroid and the positional
## error is bounded by ~700 m. WorldClim 30s cells are ~1 km at this latitude,
## so for climate normals this is at worst a one-cell displacement; it is not
## adequate for terrain-derived predictors, which is why elevation comes from
## the inventory itself (b0_tab.HoeheNN, measured per tract) rather than from a
## DEM lookup at the anonymised coordinate.
##
## Run from the project root, after 02:
##   Rscript scripts/03_attach_environment.R
##
## Outputs (outputs/):
##   corner_env.csv   one row per calibration corner: mat_c, map_mm, elev_m, x, y
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

#=#=# PARAMS #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
PARAMS <- list(
  ## WorldClim 2.1 bioclim layers to extract, as name = BIO number.
  ##   BIO1  = annual mean temperature (deg C)
  ##   BIO12 = annual precipitation (mm)
  bio_vars = c(mat_c = 1L, map_mm = 12L),

  ## Elevation source: "bwi" uses the inventory's own per-tract HoeheNN and
  ## falls back to WorldClim where it is missing; "worldclim" always uses the
  ## raster. "bwi" is preferred -- see the note on coordinate precision above.
  elev_source = "bwi"
)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

suppressPackageStartupMessages({
  library(data.table)
  library(terra)
})
source("R/config.R")
source("R/bwi_io.R")

corners <- fread(file.path(OUT_DIR, "calibration_set_corners.csv"))
message(sprintf("Calibration corners: %d on %d tracts",
                nrow(corners), uniqueN(corners$Tnr)))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 1) Coordinates
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

coords <- bwi("trakt_coords")            # Tnr, x (lon), y (lat), WGS84
d <- merge(corners[, .(Tnr, Enr, siteName, patchName, elev_m_bwi = elev_m)],
           coords, by = "Tnr", all.x = TRUE)

n_nocoord <- d[is.na(x) | is.na(y), .N]
if (n_nocoord)
  message(sprintf("  %d corners have no anonymised coordinate; dropping them.",
                  n_nocoord))
d <- d[!is.na(x) & !is.na(y)]
if (!nrow(d)) stop("No corners left with coordinates.", call. = FALSE)

message(sprintf("  extent: lon %.2f..%.2f  lat %.2f..%.2f",
                min(d$x), max(d$x), min(d$y), max(d$y)))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 2) WorldClim extraction
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

bio_dir <- file.path(WORLDCLIM_DIR, "wc2.1_30s_bio")
if (!dir.exists(bio_dir))
  stop("WorldClim bio directory not found at ", bio_dir,
       "\nSet FINN_WORLDCLIM_DIR or download the wc2.1_30s_bio bundle.",
       call. = FALSE)

pts <- vect(as.data.frame(d[, .(x, y)]), geom = c("x", "y"), crs = "EPSG:4326")

## Crop each raster to the German bounding box before extracting: the global
## 30s layers are ~1.7 GB each and we only need a small window.
win <- ext(min(d$x) - 0.5, max(d$x) + 0.5, min(d$y) - 0.5, max(d$y) + 0.5)

for (nm in names(PARAMS$bio_vars)) {
  f <- file.path(bio_dir, sprintf("wc2.1_30s_bio_%d.tif", PARAMS$bio_vars[[nm]]))
  if (!file.exists(f)) stop("Missing raster: ", f, call. = FALSE)
  r <- crop(rast(f), win)
  d[, (nm) := terra::extract(r, pts)[, 2]]
  message(sprintf("  %-8s BIO%-2d  %s", nm, PARAMS$bio_vars[[nm]],
                  paste0("median ", round(median(d[[nm]], na.rm = TRUE), 1))))
}

## Elevation
elev_f <- file.path(WORLDCLIM_DIR, "wc2.1_30s_elev.tif")
if (file.exists(elev_f)) {
  d[, elev_m_wc := terra::extract(crop(rast(elev_f), win), pts)[, 2]]
} else {
  d[, elev_m_wc := NA_real_]
  warning("wc2.1_30s_elev.tif not found; relying on the inventory elevation.",
          call. = FALSE)
}

if (identical(PARAMS$elev_source, "bwi")) {
  d[, elev_m := fifelse(is.finite(elev_m_bwi) & elev_m_bwi > 0,
                        as.numeric(elev_m_bwi), elev_m_wc)]
} else {
  d[, elev_m := elev_m_wc]
}

## Cross-check the two elevation sources: a large disagreement would mean the
## anonymised coordinates are not landing where we think they are.
both <- d[is.finite(elev_m_bwi) & elev_m_bwi > 0 & is.finite(elev_m_wc)]
if (nrow(both) > 100) {
  rho <- cor(both$elev_m_bwi, both$elev_m_wc)
  mad <- median(abs(both$elev_m_bwi - both$elev_m_wc))
  message(sprintf(
    "\n[check] inventory vs WorldClim elevation on %d corners: r = %.3f, median |diff| = %.0f m %s",
    nrow(both), rho, mad, if (rho > 0.95) "OK" else "*** CHECK COORDINATES ***"))
  if (rho <= 0.95)
    warning("Elevation sources disagree -- verify the INSPIRE coordinate join.",
            call. = FALSE)
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 3) Write
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

env_cols <- c(names(PARAMS$bio_vars), "elev_m", "x", "y")
out <- d[, c("Tnr", "Enr", "siteName", "patchName", env_cols), with = FALSE]

n_incomplete <- out[!complete.cases(out), .N]
if (n_incomplete) {
  message(sprintf("  dropping %d corners with incomplete covariates", n_incomplete))
  out <- out[complete.cases(out)]
}

fwrite(out, file.path(OUT_DIR, "corner_env.csv"))

message("\n[corner_env] covariate summary:")
print(out[, lapply(.SD, function(v) round(quantile(v, c(0, .5, 1)), 2)),
          .SDcols = env_cols])
message(sprintf("\nDone. %d corners with environment -> %s/corner_env.csv",
                nrow(out), OUT_DIR))
