# =====================================================================
# hvi_ahn.R
# Adapter: HVI auf AHN4-Kacheln (Niederlande) anwenden.
# ---------------------------------------------------------------------
# AHN4-LAZ sind bereits klassifiziert (Boden = Klasse 2), daher braucht
# es KEIN separates DTM - lidR::normalize_height() nutzt die Bodenpunkte
# direkt. CRS = RD New / Amersfoort (EPSG:28992), Einheiten Meter.
#
# Da (noch) keine Heckenpolygone vorliegen, leitet dieses Skript
# Hecken-KANDIDATEN aus dem Canopy Height Model (CHM) ab: Vegetation in
# Heckenhoehe (1-8 m), dann Filter auf kleine, LINEARE Objekte (schliesst
# Waldbloecke und Gebaeude weitgehend aus). Das ersetzt keine echte
# Kartierung, macht die Kachel aber sofort testbar.
#
# Benoetigt: lidR, sf, terra   (install.packages(c("lidR","sf","terra")))
# =====================================================================

source(file.path("R", "hvi_metrics.R"))
source(file.path("R", "hvi_index.R"))

suppressMessages({library(lidR); library(sf); library(terra)})

# ---------------------------------------------------------------------
# Hecken-Kandidaten aus einer normalisierten LAS ableiten
# ---------------------------------------------------------------------
# chm_res      : CHM-Aufloesung (m)
# h_min,h_max  : Heckenhoehenfenster
# max_area     : max. Flaeche eines Objekts (m2) - grosse Bloecke = Wald
# min_area     : min. Flaeche (Rauschen entfernen)
# min_elong    : min. Laenge/Breite-Verhaeltnis (Linearitaet) einer Hecke
# ---------------------------------------------------------------------
hvi_ahn_hedges <- function(las_norm, chm_res = 0.5,
                           h_min = 1, h_max = 8,
                           max_area = 1500, min_area = 20,
                           min_elong = 2.5) {
  chm <- rasterize_canopy(las_norm, res = chm_res, algorithm = p2r())
  veg <- terra::classify(chm,
            rbind(c(-Inf, h_min, NA), c(h_min, h_max, 1),
                  c(h_max, Inf, NA)))
  polys <- sf::st_as_sf(terra::as.polygons(veg, dissolve = TRUE))
  polys <- sf::st_cast(polys, "POLYGON", warn = FALSE)
  polys$area <- as.numeric(sf::st_area(polys))
  polys <- polys[polys$area >= min_area & polys$area <= max_area, ]
  # Linearitaet ueber das umschliessende Rechteck (OBB-Naeherung)
  elong <- vapply(sf::st_geometry(polys), function(g) {
    bb <- sf::st_bbox(g); w <- bb["xmax"]-bb["xmin"]; h <- bb["ymax"]-bb["ymin"]
    mx <- max(w, h); mn <- min(w, h); if (mn <= 0) return(0); mx/mn
  }, numeric(1))
  polys <- polys[elong >= min_elong, ]
  polys$hedge_id <- seq_len(nrow(polys))
  polys[, "hedge_id"]
}

# ---------------------------------------------------------------------
# Voll-Lauf auf einer AHN-Kachel
# ---------------------------------------------------------------------
# Fuer grosse Kacheln: erst auf eine Test-AOI zuschneiden (clip_rectangle)
# ODER den lidR-Catalog nutzen. Beispiel unten schneidet einen 1x1 km
# Ausschnitt heraus, damit ein erster Lauf schnell ist.
# ---------------------------------------------------------------------
hvi_ahn_run <- function(laz_file, aoi = NULL, requirements = NULL,
                        k_types = 4, thin = NULL) {
  las <- readLAS(laz_file, select = "xyzc",
                 filter = "-drop_z_below -1 -drop_z_above 60")
  if (is.empty(las)) stop("leere LAS: ", laz_file)
  if (!is.null(aoi)) las <- clip_roi(las, aoi)
  if (!is.null(thin)) las <- decimate_points(las, homogenize(thin, 1))
  las <- normalize_height(las, tin())
  las <- filter_poi(las, Z >= -0.5, Z < 40)

  hedges <- hvi_ahn_hedges(las)
  message(nrow(hedges), " Hecken-Kandidaten gefunden.")
  chm <- rasterize_canopy(las, res = 0.5, algorithm = p2r())
  M   <- hvi_segment_table(las, hedges, chm = chm)
  idx <- hvi_compute(M)
  res <- cbind(M["hedge_id"], idx)
  if (nrow(idx) >= k_types)
    res$hedge_type <- hvi_classify(idx, k = k_types)
  if (!is.null(requirements))
    res <- cbind(res, hvi_species_suitability(M, requirements))
  list(hedges = hedges, metrics = M, result = res)
}

# ---------------------------------------------------------------------
# Beispiel (auskommentiert):
# ---------------------------------------------------------------------
# source("R/hvi_pipeline.R")   # fuer hvi_segment_table()
# req <- read.csv("data/species_requirements.csv")
# # 1x1 km Test-AOI in der Kachelmitte (RD-Koordinaten anpassen):
# aoi <- sf::st_as_sfc(sf::st_bbox(c(xmin=256000, ymin=494000,
#          xmax=257000, ymax=495000), crs = 28992))
# out <- hvi_ahn_run("testdata/C_28FN2.LAZ", aoi = aoi, requirements = req)
# print(head(out$result))
# sf::st_write(out$hedges, "testdata/hedges_28FN2.gpkg", append = FALSE)
