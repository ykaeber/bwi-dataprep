#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-bwi | scripts/16_drought_index.R
##
## Per-site, PER-YEAR drought exposure for both BWI decades (2002-2022), from the
## DWD 1 km monthly grids.
##
## WHY THIS MATTERS
## The calibration period is 2002-2012, which contains the 2003 European
## heatwave -- one of the most severe droughts on record -- and the 2006 event.
## Until now the model saw only a static WorldClim 1970-2000 normal, so it could
## not know 2003 happened. Mortality is precisely the process that responds to
## drought, and mortality is the process that fits worst.
##
## The annual time grid built in scripts/04 already carries one env row per
## (site, year 0..10), so a genuinely time-varying forcing drops straight in
## with no change to FINN and no rework of the time axis.
##
## WHICH INDEX
## DWD publish three candidates. We compute all of them and default to the
## soil-moisture one:
##   soil_moist  modelled soil moisture, % of usable field capacity (nFK), under
##               grass over sandy loam. Physiologically the closest to what a
##               tree experiences; the grass/soil assumption is a caveat, not a
##               fatal one, since we use ANOMALIES rather than absolute values.
##   aet_pet     actual / potential evapotranspiration -- the classic water
##               stress ratio, 1 = unstressed, 0 = fully stressed.
##   demartonne  P / (T + 10), DWD's published "drought index". Crude: no
##               evapotranspiration at all, so it misses the heat component of a
##               2003-type event. Kept for comparison, not recommended.
##
## Each is reduced to the GROWING SEASON (Apr-Sep), which is when water
## limitation acts on growth and mortality, then expressed as a standardised
## anomaly against the site's own mean over all extracted years (2000-2022 by
## default, one baseline for both BWI decades). Standardising per site removes
## the strong west-east and elevation gradients already carried by the static
## climate covariates, so the drought variable adds *temporal* information
## rather than duplicating the spatial signal.
##
## Run from the project root, after 03:
##   Rscript scripts/16_drought_index.R
##
## Outputs (outputs/):
##   site_drought_annual.csv   siteName x year: all three indices, raw + anomaly
##   drought_summary.csv       per-year national means -- the sanity check
##   fig10_drought.png         the 2003 signal, and where it hit
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

suppressPackageStartupMessages({
  library(data.table); library(terra); library(ggplot2); library(patchwork)
})
source("R/config.R")

PARAMS <- list(
  ## Both BWI decades (surveys 2001-02, 2011-12, 2021-22) plus lead-in, so the
  ## covariate carries the 2003 heatwave and the 2018-2020 drought, and the
  ## per-site anomaly has one baseline for both decades.
  years        = 2000:2022,
  season       = 4:9,          # growing season months
  ## DWD grids are Gauss-Krueger zone 3 (DHDN). Not documented in the .asc
  ## header, which carries no CRS at all, so it is asserted here and then
  ## VERIFIED below by checking that the extracted values are non-NA for
  ## German forest sites -- a wrong CRS puts every site off-grid.
  dwd_crs      = "EPSG:31467",
  drought_dir  = file.path(BWI_RAW_DIR, "dwd_drought")
)

save_png <- function(path, plot, width, height, dpi = 200) {
  if (requireNamespace("ragg", quietly = TRUE))
    ggsave(path, plot, width = width, height = height, dpi = dpi,
           device = ragg::agg_png)
  else ggsave(path, plot, width = width, height = height, dpi = dpi)
}

corners <- fread(file.path(OUT_DIR, "calibration_set_corners.csv"))
env     <- fread(file.path(OUT_DIR, "corner_env.csv"))
sites   <- unique(merge(corners[, .(Tnr, Enr, siteName)],
                        env[, .(Tnr, Enr, x, y)], by = c("Tnr", "Enr")))
message(sprintf("Sites: %d", nrow(sites)))

pts_wgs <- vect(as.data.frame(sites[, .(x, y)]), geom = c("x", "y"),
                crs = "EPSG:4326")
pts_dwd <- project(pts_wgs, PARAMS$dwd_crs)

#' Read one DWD monthly ASCII grid and extract at the sites.
read_dwd <- function(var, year, month) {
  f <- file.path(PARAMS$drought_dir,
                 sprintf("grids_germany_monthly_%s_%d%02d.asc.gz", var, year, month))
  if (!file.exists(f)) return(NULL)
  tmp <- tempfile(fileext = ".asc"); on.exit(unlink(tmp), add = TRUE)
  R.utils_ok <- requireNamespace("R.utils", quietly = TRUE)
  if (R.utils_ok) R.utils::gunzip(f, destname = tmp, remove = FALSE, overwrite = TRUE)
  else system2("gunzip", c("-c", shQuote(f)), stdout = tmp)
  r <- rast(tmp); crs(r) <- PARAMS$dwd_crs
  v <- terra::extract(r, pts_dwd)[, 2]
  v[v <= -900] <- NA_real_          # DWD NODATA is -999
  v
}

VARS <- c(soil_moist = "soil_moist", evapo_p = "evapo_p", evapo_r = "evapo_r",
          demartonne = "drought_index")

message("Extracting ", length(PARAMS$years) * length(PARAMS$season),
        " grids per variable ...")
long <- rbindlist(lapply(names(VARS), function(nm) {
  rbindlist(lapply(PARAMS$years, function(y) {
    rbindlist(lapply(PARAMS$season, function(mo) {
      v <- read_dwd(VARS[[nm]], y, mo)
      if (is.null(v)) return(NULL)
      data.table(siteName = sites$siteName, year = y, month = mo,
                 var = nm, value = v)
    }))
  }))
}))

if (!nrow(long)) stop("No DWD grids could be read.", call. = FALSE)

## --- verify the asserted CRS ------------------------------------------------
cov <- long[var == "soil_moist", .(frac = mean(!is.na(value)))]$frac
message(sprintf("[check] %.1f %% of site-months have data (a wrong CRS gives ~0 %%)",
                100 * cov))
if (cov < 0.8)
  stop("Extraction coverage is too low -- the DWD CRS assumption is probably wrong.",
       call. = FALSE)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Growing-season aggregation, then standardise per site
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

## soil moisture and de Martonne: mean over the season (low = dry)
## evapotranspiration: sum over the season, then the AET/PET ratio
seas <- long[var %in% c("soil_moist", "demartonne"),
             .(value = mean(value, na.rm = TRUE)), by = .(siteName, year, var)]
et <- dcast(long[var %in% c("evapo_p", "evapo_r"),
                 .(value = sum(value, na.rm = TRUE)), by = .(siteName, year, var)],
            siteName + year ~ var, value.var = "value")
et[, value := fifelse(evapo_p > 0, evapo_r / evapo_p, NA_real_)]
seas <- rbind(seas, et[, .(siteName, year, var = "aet_pet", value)])

## Standardise within site: the anomaly is what carries temporal information;
## the site mean duplicates the static climate covariates already in env_dt.
seas[, `:=`(mu = mean(value, na.rm = TRUE), sdv = sd(value, na.rm = TRUE)),
     by = .(siteName, var)]
seas[, z := fifelse(is.finite(sdv) & sdv > 0, (value - mu) / sdv, 0)]

wide <- dcast(seas, siteName + year ~ var, value.var = c("value", "z"))
setnames(wide, sub("^value_", "", names(wide)))
fwrite(wide, file.path(OUT_DIR, "site_drought_annual.csv"))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Sanity check: 2003 must stand out
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

nat <- seas[, .(mean_z = mean(z, na.rm = TRUE),
                mean_raw = mean(value, na.rm = TRUE)), by = .(var, year)]
fwrite(nat, file.path(OUT_DIR, "drought_summary.csv"))

message("\nNational growing-season anomaly by year (negative = dry):")
print(dcast(nat, year ~ var, value.var = "mean_z")[, lapply(.SD, function(v)
  if (is.numeric(v)) round(v, 2) else v)])

worst <- nat[var == "soil_moist"][order(mean_z)][1:3]
message(sprintf("\nDriest years by soil moisture: %s",
                paste(worst$year, sprintf("(z = %.2f)", worst$mean_z),
                      collapse = ", ")))
if (!2003 %in% worst$year)
  warning("2003 is not among the three driest years -- check the extraction.",
          call. = FALSE)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Figure
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

LAB <- c(soil_moist = "soil moisture (% nFK)", aet_pet = "AET / PET",
         demartonne = "de Martonne P/(T+10)")
p1 <- ggplot(nat[var %in% names(LAB)], aes(year, mean_z)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_col(aes(fill = mean_z < 0)) +
  facet_wrap(~factor(var, names(LAB), LAB), ncol = 1) +
  scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"),
                    guide = "none") +
  scale_x_continuous(breaks = PARAMS$years) +
  labs(title = "A  Growing-season drought anomaly, national mean",
       subtitle = sprintf("standardised per site; negative = drier than that site's %d-%d average",
                          min(PARAMS$years), max(PARAMS$years)),
       x = NULL, y = "anomaly (sd)") +
  theme_bw(base_size = 9)

sm <- merge(wide[, .(siteName, year, z_soil_moist)],
            sites[, .(siteName, x, y)], by = "siteName")
p2 <- ggplot(sm[year %in% c(2003, 2006, 2010)], aes(x, y, colour = z_soil_moist)) +
  geom_point(size = 0.35) +
  facet_wrap(~year) +
  scale_colour_gradient2(low = "#b2182b", mid = "grey90", high = "#2166ac",
                         midpoint = 0, name = "anomaly (sd)") +
  coord_quickmap() +
  labs(title = "B  Where the drought hit",
       subtitle = "growing-season soil-moisture anomaly at the calibration corners",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9)

save_png(file.path(OUT_DIR, "figures", "fig10_drought.png"), p1 / p2, 11, 9)
message(sprintf("\nWrote site_drought_annual.csv, drought_summary.csv and fig10_drought.png"))
