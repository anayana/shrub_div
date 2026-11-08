# =====================================================================
# add_ndvi.R  --  NDVI fuer eine AOI aus PDOK Color-Infrarot-Luftbild
# ---------------------------------------------------------------------
# NL-weit kostenlos, kein Login. Das CIR-Bild ist Falschfarben:
#   Kanal 1 = Nahinfrarot (NIR), Kanal 2 = Rot, Kanal 3 = Gruen.
#   NDVI = (NIR - Rot) / (NIR + Rot)
# Ergebnis: testdata/ndvi_<tile>.tif (EPSG:28992). run_real() nutzt es
# automatisch, sobald die Datei existiert (Zustands-Subindex aktiv).
#
# In RStudio:
#   source("add_ndvi.R")
#   fetch_ndvi("C_28FN2.LAZ", 255000, 494000, 2000)   # gleiche AOI wie run_real!
#   # danach:  run_real("C_28FN2.LAZ", 255000, 494000, 2000)
# =====================================================================

suppressMessages({library(terra)})

# PDOK WMS (Actueel, 25 cm). Falls PDOK die URL/Layer aendert, hier anpassen.
PDOK_CIR   <- "https://service.pdok.nl/hwh/luchtfotocir/wms/v1_0"
PDOK_LAYER <- "Actueel_ortho25"

fetch_ndvi <- function(TILE = "C_28FN2.LAZ",
                       AOI_XMIN, AOI_YMIN, AOI_SIZE = 2000,
                       px_res = 1) {                 # m/Pixel
  STAMP <- sub("\\.LAZ$", "", TILE, ignore.case = TRUE)
  xmax <- AOI_XMIN + AOI_SIZE; ymax <- AOI_YMIN + AOI_SIZE
  px   <- min(2500, round(AOI_SIZE / px_res))        # WMS-Grenze beachten
  url  <- sprintf(paste0("%s?service=WMS&request=GetMap&version=1.1.1",
                         "&layers=%s&styles=&srs=EPSG:28992",
                         "&bbox=%f,%f,%f,%f&width=%d&height=%d&format=image/png"),
                  PDOK_CIR, PDOK_LAYER, AOI_XMIN, AOI_YMIN, xmax, ymax, px, px)

  tmp <- tempfile(fileext = ".png")
  message("lade CIR-Luftbild von PDOK ...")
  utils::download.file(url, tmp, mode = "wb", quiet = TRUE)

  img <- terra::rast(tmp)
  if (terra::nlyr(img) < 3) stop("CIR-Bild hat <3 Kanaele – Layer/URL pruefen. ",
                                 "Antwort evtl. XML-Fehler; URL im Browser testen:\n", url)
  terra::ext(img) <- terra::ext(AOI_XMIN, xmax, AOI_YMIN, ymax)
  terra::crs(img) <- "EPSG:28992"

  nir <- img[[1]]; red <- img[[2]]
  ndvi <- (nir - red) / (nir + red + 1e-6)           # [-1, 1]
  names(ndvi) <- "ndvi"

  out <- file.path("testdata", paste0("ndvi_", STAMP, ".tif"))
  dir.create("testdata", showWarnings = FALSE)
  terra::writeRaster(ndvi, out, overwrite = TRUE)
  message("geschrieben: ", out,
          sprintf("  (NDVI min %.2f / mean %.2f / max %.2f)",
                  terra::global(ndvi,"min",na.rm=TRUE)[1,1],
                  terra::global(ndvi,"mean",na.rm=TRUE)[1,1],
                  terra::global(ndvi,"max",na.rm=TRUE)[1,1]))
  invisible(out)
}
