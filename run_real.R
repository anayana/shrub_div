# =====================================================================
# run_real.R  --  HVI auf ECHTEN AHN4-Kacheln  (RStudio-Nutzung)
# ---------------------------------------------------------------------
# In RStudio EINMAL sourcen, dann die Funktion aufrufen:
#
#   source("run_real.R")                       # laedt run_real()
#   run_real()                                 # C_28FN2, 500-m-Fenster (Mitte)
#   run_real("C_28FN2.LAZ", 256000, 494000, 1000)   # eigene AOI (1000 m)
#   run_real("C_28FN1.LAZ")                    # andere Kachel, Mitte
#
# Schreibt nach testdata/:
#   hedges_<tile>.gpkg   Geometrie + alle Werte (in QGIS oeffnen)
#   results_<tile>.csv   Tabelle
#   map_<tile>.png       Karte (HVI + Eignung je Art)
# =====================================================================

## ---- Pakete + HVI-Funktionen (einmalig beim Sourcen) ----------------
local({
  need <- c("lidR","sf","terra")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Fehlende Pakete: ", paste(miss, collapse=", "),
     "\n  install.packages(c(", paste(sprintf('\"%s\"',miss),collapse=","), "))")
})
suppressMessages({library(lidR); library(sf); library(terra)})
source("R/hvi_metrics.R")
source("R/hvi_index.R")

## ---- Hauptfunktion --------------------------------------------------
run_real <- function(TILE = "C_28FN2.LAZ",
                     AOI_XMIN = NA, AOI_YMIN = NA, AOI_SIZE = 500) {
  EPSG  <- 28992                     # AHN: RD New / Amersfoort
  LAZ   <- file.path("testdata", TILE)
  STAMP <- sub("\\.LAZ$", "", TILE, ignore.case = TRUE)
  stopifnot(file.exists(LAZ))

  ## AOI bestimmen (Kachelmitte, falls nicht angegeben)
  h  <- readLASheader(LAZ)
  if (is.na(AOI_XMIN) || is.na(AOI_YMIN)) {
    cx <- mean(c(h@PHB[["Min X"]], h@PHB[["Max X"]]))
    cy <- mean(c(h@PHB[["Min Y"]], h@PHB[["Max Y"]]))
    AOI_XMIN <- cx - AOI_SIZE/2; AOI_YMIN <- cy - AOI_SIZE/2
  }
  AOI_XMAX <- AOI_XMIN + AOI_SIZE; AOI_YMAX <- AOI_YMIN + AOI_SIZE
  message(sprintf("Kachel %s | AOI %.0f,%.0f – %.0f,%.0f (%.0f m)",
                  TILE, AOI_XMIN, AOI_YMIN, AOI_XMAX, AOI_YMAX, AOI_SIZE))

  ## Nur die AOI einlesen (Filter beim Lesen -> kein 3-GB-Load)
  las <- readLAS(LAZ, select = "xyzc",
                 filter = sprintf("-keep_xy %f %f %f %f -drop_z_below -1 -drop_z_above 60",
                                  AOI_XMIN, AOI_YMIN, AOI_XMAX, AOI_YMAX))
  if (is.empty(las)) stop("Keine Punkte in der AOI – andere Koordinaten waehlen.")
  if (is.na(crs(las)) || crs(las) == "") crs(las) <- paste0("EPSG:", EPSG)
  message(sprintf("%s Punkte gelesen (%.1f Pkt/m2)",
                  format(npoints(las), big.mark=" "), npoints(las)/(AOI_SIZE^2)))

  ## Hoehennormalisierung (AHN4-Boden = Klasse 2)
  if (!any(las$Classification == 2L)) las <- classify_ground(las, csf())
  las <- normalize_height(las, tin())
  las <- filter_poi(las, Z >= -0.5, Z < 40)

  ## optional: NDVI-Raster (aus add_ndvi.R), aktiviert den Zustands-Subindex
  ndvi_path <- file.path("testdata", paste0("ndvi_", STAMP, ".tif"))
  ndvi_rast <- if (file.exists(ndvi_path)) terra::rast(ndvi_path) else NULL
  message(if (is.null(ndvi_rast)) "kein NDVI (Zustand aus Gewichtung genommen)"
          else paste("NDVI genutzt:", ndvi_path))

  ## CHM + Hecken-Kandidaten (Vegetation 1–8 m, klein & linear)
  chm <- rasterize_canopy(las, res = 0.5, algorithm = p2r())
  veg <- classify(chm, matrix(c(-Inf,1,NA, 1,8,1, 8,Inf,NA), ncol = 3, byrow = TRUE))
  polys <- st_as_sf(as.polygons(veg, dissolve = TRUE))
  polys <- st_cast(polys, "POLYGON", warn = FALSE)
  polys$area <- as.numeric(st_area(polys))
  polys <- polys[polys$area >= 20 & polys$area <= 2000, ]
  elong <- vapply(st_geometry(polys), function(g){
    b <- st_bbox(g); w <- b["xmax"]-b["xmin"]; ht <- b["ymax"]-b["ymin"]
    mx <- max(w,ht); mn <- min(w,ht); if (mn <= 0) 0 else mx/mn
  }, numeric(1))
  polys <- polys[elong >= 2.5, ]
  if (nrow(polys) == 0) stop("Keine Hecken-Kandidaten – Parameter lockern.")
  polys$hedge_id <- seq_len(nrow(polys)); polys <- polys[, "hedge_id"]
  st_crs(polys) <- EPSG
  message(nrow(polys), " Hecken-Kandidaten gefunden.")

  ## Metriken je Hecke (echte Punkte je Polygon)
  M <- do.call(rbind, lapply(seq_len(nrow(polys)), function(i){
    sub <- clip_roi(las, polys[i, ]); if (is.empty(sub)) return(NULL)
    cv  <- terra::extract(chm, vect(polys[i, ]))[, 2]
    nd  <- if (!is.null(ndvi_rast))
             mean(terra::extract(ndvi_rast, vect(polys[i, ]))[, 2], na.rm = TRUE) else NA
    m   <- hvi_segment_metrics(sub$Z, chm_vals = cv, ndvi = nd)
    data.frame(hedge_id = polys$hedge_id[i], t(m))
  }))
  polys <- polys[polys$hedge_id %in% M$hedge_id, ]

  ## HVI + artspezifische Eignung
  idx <- hvi_compute(M)
  req <- read.csv("data/species_requirements.csv", stringsAsFactors = FALSE)
  hsi <- hvi_species_suitability(M, req)
  spcols <- c("Neuntoeter","Dorngrasmuecke","Goldammer","Fledermaus","Wildbienen")
  colnames(hsi) <- spcols[seq_len(ncol(hsi))]
  res <- cbind(M, idx[, c("vertical_complexity","volume_size","heterogeneity","HVI")], hsi)
  if (nrow(idx) >= 4) res$hedge_type <- hvi_classify(idx, k = 4)

  ## Export GPKG + CSV
  out  <- merge(polys, res, by = "hedge_id")
  gpkg <- file.path("testdata", paste0("hedges_", STAMP, ".gpkg"))
  csv  <- file.path("testdata", paste0("results_", STAMP, ".csv"))
  st_write(out, gpkg, delete_dsn = TRUE, quiet = TRUE)
  write.csv(res, csv, row.names = FALSE)
  message("geschrieben: ", gpkg, " | ", csv)

  ## Karte: HVI + Eignung je Art
  ramp <- colorRamp(c("#f4e79e","#c9de84","#7cc08a","#2f8f5b","#14532d"))
  col_for <- function(v){ v[is.na(v)] <- 0
    apply(ramp(pmin(pmax(v,0),1)), 1, function(r) rgb(r[1],r[2],r[3],maxColorValue=255)) }
  legbar <- function(){ xs <- seq(0,1,length.out=100); usr <- par("usr")
    x0<-usr[1]+0.04*(usr[2]-usr[1]); x1<-usr[1]+0.30*(usr[2]-usr[1])
    y0<-usr[3]+0.04*(usr[4]-usr[3]); y1<-y0+0.03*(usr[4]-usr[3])
    rect(x0+(x1-x0)*xs[-100], y0, x0+(x1-x0)*xs[-1], y1, col=col_for(xs[-1]), border=NA)
    text(c(x0,x1), y1, c("0","1"), pos=3, cex=0.7, col="grey30", offset=0.15) }
  # auf die Hecken zoomen (+ Rand), sonst verschwinden sie in leerer Flaeche
  e  <- st_bbox(out); mrg <- 0.08*max(e[3]-e[1], e[4]-e[2]) + 10
  xl <- c(e[1]-mrg, e[3]+mrg); yl <- c(e[2]-mrg, e[4]+mrg)
  panels <- c(HVI="HVI", Neuntoeter="Neuntöter", Dorngrasmuecke="Dorngrasmücke",
              Goldammer="Goldammer", Fledermaus="Fledermaus", Wildbienen="Wildbienen")
  panels <- panels[names(panels) %in% names(res)]
  png(file.path("testdata", paste0("map_", STAMP, ".png")),
      width = 1500, height = 1000, res = 130)
  op <- par(mfrow = c(2,3), mar = c(1,1,2.4,1))
  for (nm in names(panels)) {
    plot(st_geometry(out), col = col_for(out[[nm]]), border="grey40", lwd=0.6,
         xlim=xl, ylim=yl, main = panels[[nm]], reset = FALSE)
    legbar()
  }
  par(op); dev.off()
  message("geschrieben: testdata/map_", STAMP, ".png")

  cat("\nTop-5 Hecken nach HVI:\n")
  print(head(res[order(-res$HVI),
      c("hedge_id","HVI","h_p95","fhd","n_layers", intersect(spcols, names(res)))], 5),
      row.names = FALSE)
  invisible(res)
}

## ---- Auto-Start nur im Terminal (Rscript); in RStudio ignoriert -----
if (!interactive()) {
  a <- commandArgs(trailingOnly = TRUE)
  run_real(if (length(a) >= 1) a[1] else "C_28FN2.LAZ",
           if (length(a) >= 3) as.numeric(a[2]) else NA,
           if (length(a) >= 3) as.numeric(a[3]) else NA,
           if (length(a) >= 4) as.numeric(a[4]) else 500)
}
