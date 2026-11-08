# =====================================================================
# run_all.R  --  HVI ueber ALLE AHN4-Kacheln, direkt aus dem Original-
#                Ordner (kein Kopieren, kein Download noetig).
# ---------------------------------------------------------------------
# Verarbeitet den kompletten LAScatalog gekachelt (Chunks), leitet je
# Chunk Hecken aus dem CHM ab und sammelt ALLE Hecken zusammen. Der HVI
# wird EINMAL global (ueber alle Hecken) normiert, damit die Werte ueber
# den gesamten Datensatz vergleichbar sind.
#
# Ausgabe in testdata/all/:
#   hedges_all.gpkg   alle Hecken + Metriken + HVI + Arteignung (QGIS)
#   results_all.csv   Tabelle
#   map_all.png       Karte (HVI + Eignung je Art)
#
# Aufruf aus dem Projektordner shrub_div:
#   Rscript run_all.R
# Optional erste Zahl = Kantenlaenge der Verarbeitungs-Chunks in m:
#   Rscript run_all.R 500
# =====================================================================

## ---- Konfiguration --------------------------------------------------
LAZ_DIR <- "C:/Users/Admin/R/agroforstAHN/agroforst_pipeline/data/raw/ahn/ahn4/laz"
OUT_DIR <- file.path("testdata", "all")
EPSG    <- 28992
args    <- commandArgs(trailingOnly = TRUE)
CHUNK   <- if (length(args) >= 1) as.numeric(args[1]) else 500   # m
THIN    <- NA        # z.B. 8 -> auf 8 Pkt/m2 ausduennen (schneller); NA = volle Dichte

## ---- Pakete + HVI-Funktionen ---------------------------------------
need <- c("lidR","sf","terra")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Fehlende Pakete: ", paste(miss, collapse=", "),
   "\n  install.packages(c(", paste(sprintf('\"%s\"',miss),collapse=","), "))")
suppressMessages({library(lidR); library(sf); library(terra)})
source("R/hvi_metrics.R"); source("R/hvi_index.R")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(LAZ_DIR))
req <- read.csv("data/species_requirements.csv", stringsAsFactors = FALSE)

## ---- Catalog aufsetzen ---------------------------------------------
ctg <- readLAScatalog(LAZ_DIR)
opt_select(ctg)        <- "xyzc"
opt_filter(ctg)        <- "-drop_z_below -1 -drop_z_above 60"
opt_chunk_size(ctg)    <- CHUNK
opt_chunk_buffer(ctg)  <- 0          # kein Buffer -> keine Doppelzaehlung
opt_progress(ctg)      <- TRUE
message(sprintf("Katalog: %d Kacheln, Chunk %d m%s",
        length(ctg$filename), CHUNK,
        if (is.na(THIN)) "" else sprintf(", ausgeduennt auf %.0f Pkt/m2", THIN)))

## ---- Pro Chunk: Hecken + ROH-Metriken (noch kein globaler HVI) -----
hvi_chunk <- function(chunk, req, thin) {
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  if (!is.na(thin)) las <- decimate_points(las, homogenize(thin, 1))
  if (!any(las$Classification == 2L)) las <- classify_ground(las, csf())
  las <- normalize_height(las, tin())
  las <- filter_poi(las, Z >= -0.5, Z < 40)
  if (npoints(las) < 500) return(NULL)

  chm  <- rasterize_canopy(las, res = 0.5, algorithm = p2r())
  veg  <- terra::classify(chm, matrix(c(-Inf,1,NA, 1,8,1, 8,Inf,NA),
                                      ncol = 3, byrow = TRUE))
  polys <- sf::st_as_sf(terra::as.polygons(veg, dissolve = TRUE))
  polys <- sf::st_cast(polys, "POLYGON", warn = FALSE)
  if (nrow(polys) == 0) return(NULL)
  polys$area <- as.numeric(sf::st_area(polys))
  polys <- polys[polys$area >= 20 & polys$area <= 2000, ]
  if (nrow(polys) == 0) return(NULL)
  elong <- vapply(sf::st_geometry(polys), function(g){
    b <- sf::st_bbox(g); w <- b["xmax"]-b["xmin"]; ht <- b["ymax"]-b["ymin"]
    mx <- max(w,ht); mn <- min(w,ht); if (mn <= 0) 0 else mx/mn
  }, numeric(1))
  polys <- polys[elong >= 2.5, ]
  if (nrow(polys) == 0) return(NULL)

  rows <- lapply(seq_len(nrow(polys)), function(i){
    sub <- clip_roi(las, polys[i, ])
    if (is.empty(sub)) return(NULL)
    cv  <- terra::extract(chm, terra::vect(polys[i, ]))[, 2]
    m   <- hvi_segment_metrics(sub$Z, chm_vals = cv, ndvi = NA)
    sf::st_sf(as.data.frame(t(m)), geometry = sf::st_geometry(polys[i, ]))
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  sf::st_crs(out) <- EPSG
  out
}

## ---- Ausfuehren -----------------------------------------------------
res <- catalog_apply(ctg, hvi_chunk, req = req, thin = THIN,
                     .options = list(automerge = TRUE))
if (is.list(res) && !inherits(res, "sf")) {
  res <- res[!vapply(res, is.null, logical(1))]
  res <- do.call(rbind, res)
}
if (is.null(res) || nrow(res) == 0) stop("Keine Hecken gefunden.")
res$hedge_id <- seq_len(nrow(res))
message(nrow(res), " Hecken ueber den gesamten Datensatz.")

## ---- HVI + Eignung EINMAL global -----------------------------------
M   <- sf::st_drop_geometry(res)
idx <- hvi_compute(M)
hsi <- hvi_species_suitability(M, req)
spcols <- c("Neuntoeter","Dorngrasmuecke","Goldammer","Fledermaus","Wildbienen")
colnames(hsi) <- spcols[seq_len(ncol(hsi))]

final <- cbind(res["hedge_id"], M,
               idx[, c("vertical_complexity","volume_size","heterogeneity","HVI")],
               hsi)
if (nrow(idx) >= 4) final$hedge_type <- hvi_classify(idx, k = 4)

## ---- Export ---------------------------------------------------------
gpkg <- file.path(OUT_DIR, "hedges_all.gpkg")
csv  <- file.path(OUT_DIR, "results_all.csv")
sf::st_write(final, gpkg, delete_dsn = TRUE, quiet = TRUE)
write.csv(sf::st_drop_geometry(final), csv, row.names = FALSE)
message("geschrieben: ", gpkg, " | ", csv)

## ---- Karte ----------------------------------------------------------
ramp <- colorRamp(c("#e8e8e6","#e7c46a","#5fae74","#1f6b41"))
colf <- function(v){ v[is.na(v)] <- 0
  apply(ramp(pmin(pmax(v,0),1)), 1, function(r) rgb(r[1],r[2],r[3],maxColorValue=255)) }
panels <- c(HVI="HVI", Neuntoeter="Neuntöter", Dorngrasmuecke="Dorngrasmücke",
            Goldammer="Goldammer", Fledermaus="Fledermaus", Wildbienen="Wildbienen")
panels <- panels[names(panels) %in% names(final)]
png(file.path(OUT_DIR, "map_all.png"), width = 1600, height = 1100, res = 130)
op <- par(mfrow = c(2, 3), mar = c(1,1,2.4,1))
for (nm in names(panels)) {
  plot(sf::st_geometry(final), col = colf(final[[nm]]), border = NA,
       main = panels[[nm]]); }
par(op); dev.off()
message("geschrieben: ", file.path(OUT_DIR, "map_all.png"))

cat("\nTop-10 Hecken nach HVI:\n")
print(head(final[order(-final$HVI),
   c("hedge_id","HVI","h_p95","fhd","n_layers", intersect(spcols, names(final)))], 10),
   row.names = FALSE)
cat("\nFertig – alle Kacheln verarbeitet.\n")
