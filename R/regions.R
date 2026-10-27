#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-bwi | R/regions.R
##
## Named German forest regions, defined by a bounding box + elevation window +
## (optionally) Bundesland. Deliberately coarse: the point is to pick sites that
## a German forest scientist would recognise as "a Brandenburg pine forest" or
## "a Black Forest spruce-fir stand", not to delineate Wuchsgebiete precisely.
##
## Site coordinates are the anonymised 1 km INSPIRE cell centroids, so a box has
## roughly +/- 700 m of positional slack -- irrelevant at this scale.
##
## Each region carries the forest type it SHOULD show. That expectation is what
## makes the simulation assessable: scripts/14_regions.R first checks the
## observed composition against it (validating the box), then runs the stands
## forward.
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

suppressPackageStartupMessages(library(data.table))

## Bundesland codes: 1 SH, 3 NI, 5 NW, 6 HE, 7 RP, 8 BW, 9 BY, 12 BB,
##                   13 MV, 14 SN, 15 ST, 16 TH
GERMAN_REGIONS <- rbindlist(list(
  data.table(region = "Schwarzwald (Black Forest)",
             lon_min = 7.6, lon_max = 8.6, lat_min = 47.7, lat_max = 48.9,
             elev_min = 600, elev_max = 1500, bl = 8,
             expected = "Picea / Abies / Fagus montane mixed"),

  data.table(region = "Oberrheinebene (Rhine valley)",
             lon_min = 7.5, lon_max = 8.6, lat_min = 47.8, lat_max = 49.8,
             elev_min = -10, elev_max = 200, bl = NA_integer_,
             expected = "Quercus / Fagus / Pinus warm lowland"),

  data.table(region = "Harz",
             lon_min = 10.2, lon_max = 11.1, lat_min = 51.5, lat_max = 51.95,
             elev_min = 350, elev_max = 1200, bl = NA_integer_,
             expected = "Picea-dominated montane"),

  data.table(region = "Brandenburg lowland",
             lon_min = 11.3, lon_max = 14.8, lat_min = 51.3, lat_max = 53.6,
             elev_min = -10, elev_max = 120, bl = 12,
             expected = "Pinus sylvestris on sandy, continental, dry sites"),

  data.table(region = "Bayerischer Wald",
             lon_min = 12.6, lon_max = 13.9, lat_min = 48.6, lat_max = 49.5,
             elev_min = 600, elev_max = 1500, bl = 9,
             expected = "Picea / Fagus / Abies montane"),

  data.table(region = "Alpenvorland / Alps",
             lon_min = 10.0, lon_max = 13.2, lat_min = 47.3, lat_max = 47.9,
             elev_min = 600, elev_max = 1800, bl = 9,
             expected = "Picea / Abies / Fagus, high precipitation"),

  data.table(region = "Schleswig-Holstein / Baltic",
             lon_min = 8.5, lon_max = 14.5, lat_min = 53.5, lat_max = 55.0,
             elev_min = -10, elev_max = 150, bl = NA_integer_,
             expected = "Fagus / Quercus maritime lowland"),

  data.table(region = "Eifel / Hunsrueck",
             lon_min = 6.1, lon_max = 7.6, lat_min = 49.5, lat_max = 50.8,
             elev_min = 350, elev_max = 900, bl = NA_integer_,
             expected = "Fagus / Quercus / Picea submontane")
), fill = TRUE)

#' Assign sites to the named regions.
#'
#' @param d data.table with columns x (lon), y (lat), elev_m and optionally Bl
#' @param regions region definition table
#' @return `d` with a `region` column (NA where a site matches none)
#'
#' A site matching more than one box takes the first match, so order the table
#' from most specific to least.
assign_regions <- function(d, regions = GERMAN_REGIONS) {
  d <- copy(d)
  d[, region := NA_character_]
  for (i in seq_len(nrow(regions))) {
    r <- regions[i]
    hit <- is.na(d$region) &
      d$x >= r$lon_min & d$x <= r$lon_max &
      d$y >= r$lat_min & d$y <= r$lat_max &
      d$elev_m >= r$elev_min & d$elev_m <= r$elev_max &
      (is.na(r$bl) | (!is.null(d$Bl) & d$Bl == r$bl))
    hit[is.na(hit)] <- FALSE
    d[hit, region := r$region]
  }
  d[]
}
