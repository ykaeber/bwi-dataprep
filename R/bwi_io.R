#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-bwi | R/bwi_io.R
##
## Reading the BWI "erweiterte Grunddaten" out of the Thuenen ACCESS2003
## databases. The .mdb files are read once with mdbtools and cached as CSV in
## data/bwi_csv/; everything downstream reads the CSVs, so mdbtools is an
## import-time dependency only.
##
##   brew install mdbtools        (macOS)
##   apt-get install mdbtools     (Debian/Ubuntu)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

suppressPackageStartupMessages(library(data.table))

## mdbtools lives in /opt/homebrew/bin on Apple silicon, which is not always on
## the PATH a non-interactive Rscript inherits. Find the binaries explicitly.
mdb_bin <- function(cmd = "mdb-export") {
  found <- Sys.which(cmd)
  if (nzchar(found)) return(unname(found))
  for (d in c("/opt/homebrew/bin", "/usr/local/bin", "/usr/bin")) {
    p <- file.path(d, cmd)
    if (file.exists(p)) return(p)
  }
  stop("`", cmd, "` not found. Install mdbtools:\n",
       "  macOS : brew install mdbtools\n",
       "  Debian: sudo apt-get install mdbtools", call. = FALSE)
}

#' List the tables in an Access database.
mdb_tables <- function(mdb) {
  stopifnot(file.exists(mdb))
  out <- system2(mdb_bin("mdb-tables"), c("-1", shQuote(mdb)), stdout = TRUE)
  trimws(out[nzchar(trimws(out))])
}

#' Export one table from an Access database to a data.table.
#'
#' @param mdb   path to the .mdb file
#' @param table table name
#' @param ...   passed to data.table::fread
mdb_read <- function(mdb, table, ...) {
  stopifnot(file.exists(mdb))
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  st <- system2(mdb_bin("mdb-export"), c(shQuote(mdb), shQuote(table)),
                stdout = tmp, stderr = FALSE)
  if (!identical(st, 0L) || !file.exists(tmp) || file.size(tmp) == 0)
    stop("mdb-export failed for table '", table, "' in ", basename(mdb),
         call. = FALSE)
  ## BWI encodes missing values as -1 in coded (nominal/ordinal) columns; we
  ## keep them as-is here and let the callers decide, because 0 and -1 carry
  ## different meanings per code table.
  fread(tmp, showProgress = FALSE, ...)
}

#' Export a table to data/bwi_csv/<name>.csv, reusing the cache when fresh.
#'
#' @param mdb    path to the .mdb file
#' @param table  table name inside the database
#' @param name   output basename (defaults to `table`)
#' @param csv_dir output directory
#' @param force  re-export even when the cached CSV is newer than the .mdb
mdb_cache <- function(mdb, table, name = table, csv_dir = BWI_CSV_DIR,
                      force = FALSE) {
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(csv_dir, paste0(name, ".csv"))
  if (!force && file.exists(out) &&
      file.mtime(out) > file.mtime(mdb)) {
    return(invisible(out))
  }
  dt <- mdb_read(mdb, table)
  fwrite(dt, out)
  message(sprintf("  %-22s %8d rows x %2d cols -> %s",
                  table, nrow(dt), ncol(dt), basename(out)))
  invisible(out)
}

#' Read a cached BWI table from data/bwi_csv/.
bwi <- function(name, csv_dir = BWI_CSV_DIR, ...) {
  p <- file.path(csv_dir, paste0(name, ".csv"))
  if (!file.exists(p))
    stop("Cached table '", name, "' not found at ", p,
         "\nRun: Rscript scripts/00_import_bwi.R", call. = FALSE)
  fread(p, showProgress = FALSE, ...)
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Species lookup
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

#' Build a tidy species table from the BWI code table `x_ba`.
#'
#' `x_ba` carries `Gattung` (genus) and `Art` (epithet) alongside the German
#' name, so we can produce proper binomials. Entries whose epithet is a
#' placeholder such as "(sonstige)" or "-" are aggregate classes ("Sonstige
#' Laubbaeume") and get `species_name = "<Genus> spp."` instead.
#'
#' @return data.table(Ba, ba_short, name_de, genus, epithet, species_name,
#'                    ba_group, is_aggregate)
bwi_species_table <- function(csv_dir = BWI_CSV_DIR) {
  x <- bwi("x_ba", csv_dir = csv_dir)
  setnames(x, names(x), trimws(names(x)))
  sp <- data.table(
    Ba       = as.integer(x$ICode),
    ba_short = trimws(as.character(x$KurzD)),
    name_de  = trimws(as.character(x$LangD)),
    genus    = trimws(as.character(x$Gattung)),
    epithet  = trimws(as.character(x$Art)),
    ba_group = trimws(as.character(x$BaGr))
  )
  sp <- sp[!is.na(Ba) & Ba >= 0]

  ## Aggregate classes ("Sonstige Laubbaeume", "UeLN") carry a bracketed or
  ## underscore-prefixed placeholder in Gattung/Art rather than a real name.
  placeholder <- function(v) is.na(v) | v == "" | v == "-" |
    grepl("^[_(]", v) | grepl("sonstige|spec\\.|spp\\.", v, ignore.case = TRUE)

  sp[, is_aggregate := placeholder(epithet)]
  sp[, species_name := fifelse(
    placeholder(genus), paste0("BA_", Ba),
    fifelse(is_aggregate, paste(genus, "spp."), paste(genus, epithet)))]

  ## Spelling fixes for the source table. x_ba is a forestry code list, not a
  ## nomenclatural authority, and carries a handful of misspelt epithets; left
  ## uncorrected these split a species across two names in any downstream join
  ## against an external trait or distribution database.
  fixes <- c(
    "Fraxinus exelsior"  = "Fraxinus excelsior",
    "Quercus robus"      = "Quercus robur",
    "Prunus avum"        = "Prunus avium",
    "Sorbus aucuparius"  = "Sorbus aucuparia"
  )
  hit <- sp$species_name %in% names(fixes)
  if (any(hit)) sp[hit, species_name := unname(fixes[species_name])]

  setkey(sp, Ba)
  sp[]
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Coordinates
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

#' Anonymised Trakt coordinates as WGS84 lon/lat.
#'
#' Thuenen release tract positions only as the 1 km INSPIRE (ETRS89-LAEA,
#' EPSG:3035) grid cell containing the tract centre. We take the cell centroid,
#' which bounds the positional error at ~700 m — well inside the 1 km cell of
#' WorldClim 30s, so climate extraction is unaffected in any practical sense.
#'
#' @return data.table(Tnr, x, y)  with x = longitude, y = latitude (WGS84)
bwi_trakt_coords <- function(shp = INSPIRE_SHP) {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required to read the INSPIRE shapefile.", call. = FALSE)
  if (!file.exists(shp))
    stop("INSPIRE shapefile not found at ", shp,
         "\nDownload ShapeFile_Tnr_INSPIRE_Poly1000.zip from\n",
         "  https://bwi.info/Download/de/BWI-Basisdaten/ACCESS2003/",
         call. = FALSE)

  g <- sf::st_read(shp, quiet = TRUE)
  nm <- names(g)
  tnr_col <- nm[tolower(nm) %in% c("tnr", "traktnr", "trakt")][1]
  if (is.na(tnr_col))
    stop("No Tnr column in ", basename(shp), "; columns are: ",
         paste(nm, collapse = ", "), call. = FALSE)

  ## Cell centroid -> WGS84. suppressWarnings: st_centroid on lon/lat is only
  ## approximate, but the layer is projected (LAEA) so this is exact here.
  ctr <- suppressWarnings(sf::st_centroid(sf::st_geometry(g)))
  ll  <- sf::st_coordinates(sf::st_transform(ctr, 4326))

  data.table(Tnr = as.integer(g[[tnr_col]]),
             x   = ll[, 1],
             y   = ll[, 2])[!is.na(Tnr)][, .(x = mean(x), y = mean(y)), by = Tnr]
}
