# =====================================================================
# hvi_pipeline.R
# End-to-End-Workflow von roher ALS-Punktwolke zum HVI je Heckensegment
# ---------------------------------------------------------------------
# Benoetigt: lidR, sf, terra  (fuer echte Daten)
#   install.packages(c("lidR","sf","terra"))
# Frei verfuegbare ALS-Quellen (DE, CC-BY / OpenData):
#   NRW  : https://www.opengeodata.nrw.de  (3dm_l_las, ~ >10 Pkt/m2)
#   Bayern: https://geodaten.bayern.de       (OpenData seit 2023)
#   Niedersachsen, Hamburg, Berlin, Brandenburg, Thueringen ...
#   bundesweit: BKG DigiZ-DE
# Zusatz: Sentinel-2 / DOP-CIR fuer NDVI (Bluehsaum-Proxy).
# =====================================================================

source(file.path("R", "hvi_metrics.R"))
source(file.path("R", "hvi_index.R"))

# ---------------------------------------------------------------------
# 1) Einlesen + Hoehennormalisierung
# ---------------------------------------------------------------------
# las_file: .las/.laz. Gibt eine hoehennormalisierte LAS zurueck,
# in der Z = Hoehe ueber Grund.
# ---------------------------------------------------------------------
hvi_prepare_las <- function(las_file, dtm_res = 0.5) {
  stopifnot(requireNamespace("lidR", quietly = TRUE))
  las <- lidR::readLAS(las_file, select = "xyzc")
  if (lidR::is.empty(las)) stop("leere/ungueltige LAS: ", las_file)
  # Bodenklassifikation nur falls noch nicht klassifiziert
  if (!any(las$Classification == 2L))
    las <- lidR::classify_ground(las, lidR::csf())
  las <- lidR::normalize_height(las, lidR::tin())
  las <- lidR::filter_poi(las, Z >= -0.5, Z < 40)   # Ausreisser weg
  las
}

# ---------------------------------------------------------------------
# 2) Heckenmaske / Segmentierung
# ---------------------------------------------------------------------
# Praxis: Hecken aus einem nDOM (CHM) ableiten -> Vegetation zwischen
# ~0.5 und ~8 m, linear/schmal. Hier zwei Wege:
#   (a) uebergebene Polygone (sf) der Hecken (z.B. aus Biotopkataster)
#   (b) einfache CHM-Schwelle + Segmentierung (Fallback)
# Rueckgabe: sf-Polygone mit Spalte hedge_id.
# ---------------------------------------------------------------------
hvi_chm <- function(las, res = 0.5) {
  lidR::rasterize_canopy(las, res = res, algorithm = lidR::p2r())
}

# ---------------------------------------------------------------------
# 3) Metriken je Heckensegment
# ---------------------------------------------------------------------
# hedges: sf-Polygone mit hedge_id. ndvi_rast: optional terra-Raster.
# Rueckgabe: data.frame (eine Zeile je hedge_id) mit allen Kennwerten.
# ---------------------------------------------------------------------
hvi_segment_table <- function(las, hedges, chm = NULL,
                              ndvi_rast = NULL, dz = 0.5) {
  stopifnot(requireNamespace("sf", quietly = TRUE))
  if (is.null(chm)) chm <- hvi_chm(las, res = dz)
  ids <- hedges$hedge_id
  rows <- lapply(seq_along(ids), function(i) {
    poly <- hedges[i, ]
    sub  <- lidR::clip_roi(las, poly)
    if (lidR::is.empty(sub)) return(NULL)
    z <- sub$Z
    chm_vals <- terra::extract(chm, terra::vect(poly))[, 2]
    ndvi <- if (!is.null(ndvi_rast))
      mean(terra::extract(ndvi_rast, terra::vect(poly))[, 2],
           na.rm = TRUE) else NA
    m <- hvi_segment_metrics(z, chm_vals = chm_vals,
                             ndvi = ndvi, dz = dz)
    data.frame(hedge_id = ids[i], t(m))
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------
# 4) Voll-Workflow
# ---------------------------------------------------------------------
hvi_run <- function(las_file, hedges, ndvi_rast = NULL,
                    requirements = NULL, k_types = 4) {
  las <- hvi_prepare_las(las_file)
  chm <- hvi_chm(las)
  M   <- hvi_segment_table(las, hedges, chm = chm,
                           ndvi_rast = ndvi_rast)
  idx <- hvi_compute(M)
  res <- cbind(M["hedge_id"], idx)
  res$hedge_type <- hvi_classify(idx, k = k_types)
  if (!is.null(requirements)) {
    hsi <- hvi_species_suitability(M, requirements)
    res <- cbind(res, hsi)
  }
  list(metrics = M, result = res)
}
