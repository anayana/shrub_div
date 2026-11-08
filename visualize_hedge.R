# =====================================================================
# visualize_hedge.R  --  ECHTE Punktwolke einer Hecke + Parameter-Profil
# ---------------------------------------------------------------------
# Clippt die realen LiDAR-Punkte EINER Hecke aus der Kachel, zeigt das
# Seitenprofil (Punktwolke, nach Hoehenschicht gefaerbt) und daneben,
# wie hoch die Hecke bei jedem Indexparameter liegt (0..1-Balken + Wert).
# Voraussetzung: run_real() wurde gelaufen (es gibt testdata/hedges_<tile>.gpkg).
#
# In RStudio:
#   source("visualize_hedge.R")
#   show_hedge("C_28FN2.LAZ")            # Hecke mit hoechstem HVI
#   show_hedge("C_28FN2.LAZ", id = 18)   # bestimmte Hecke
# Schreibt testdata/hedge_<tile>_<id>.png und plottet ins RStudio-Fenster.
# =====================================================================

suppressMessages({library(lidR); library(sf); library(terra)})
source("R/hvi_metrics.R")

STRATA_COL <- c("#b79b6a","#9bbf6a","#6fae5f","#3f8f52","#245f39")
STRATA_LAB <- c("0–0,5 m","0,5–1,5 m","1,5–3 m","3–6 m",">6 m")
band_of <- function(z) findInterval(z, c(0.5,1.5,3,6)) + 1

show_hedge <- function(TILE = "C_28FN2.LAZ", id = NULL, buffer = 4) {
  STAMP <- sub("\\.LAZ$", "", TILE, ignore.case = TRUE)
  LAZ   <- file.path("testdata", TILE)
  gpkg  <- file.path("testdata", paste0("hedges_", STAMP, ".gpkg"))
  stopifnot(file.exists(LAZ), file.exists(gpkg))
  hedges <- sf::st_read(gpkg, quiet = TRUE)
  if (is.null(id)) id <- hedges$hedge_id[which.max(hedges$HVI)]
  h <- hedges[hedges$hedge_id == id, ]
  if (nrow(h) == 0) stop("hedge_id ", id, " nicht gefunden. Vorhanden: ",
                         paste(sort(hedges$hedge_id), collapse = ", "))

  ## echte Punkte der Hecke holen (nur bbox lesen -> schnell)
  b <- as.numeric(sf::st_bbox(h))
  las <- readLAS(LAZ, select = "xyzc",
     filter = sprintf("-keep_xy %f %f %f %f", b[1]-buffer, b[2]-buffer,
                      b[3]+buffer, b[4]+buffer))
  if (is.empty(las)) stop("keine Punkte fuer Hecke ", id)
  if (is.na(crs(las)) || crs(las) == "") crs(las) <- "EPSG:28992"
  las <- normalize_height(las, tin())
  sub <- clip_roi(las, h)
  sub <- filter_poi(sub, Z >= 0, Z < 40)
  if (is.empty(sub) || npoints(sub) < 20) stop("zu wenige Punkte in Hecke ", id)

  X <- sub$X; Y <- sub$Y; Z <- sub$Z
  ## Laengsachse der Hecke via PCA -> Punkte darauf projizieren
  pc  <- prcomp(cbind(X, Y)); along <- pc$x[, 1]; along <- along - min(along)

  ## Kennwerte aus dem GPKG (real)
  M <- sf::st_drop_geometry(h)
  getn <- function(k) suppressWarnings(as.numeric(M[[k]]))
  params <- list(
    c("FHD (Schichtungsvielfalt)",        getn("fhd"),                 2.5, "%.2f"),
    c("VCI (norm. Schichtung)",           getn("vci"),                 1.0, "%.2f"),
    c("Schichten (Anzahl)",               getn("n_layers"),            5.0, "%.0f"),
    c("Evenness (Gleichverteilung)",      getn("strata_evenness"),     1.0, "%.2f"),
    c("Unterwuchs-Offenheit",             getn("understory_openness"), 1.0, "%.2f"),
    c("Vertikale Lueckigkeit",            getn("vertical_gap"),        1.0, "%.2f"),
    c("P95-Hoehe (m)",                    getn("h_p95"),              10.0, "%.1f"),
    c("Kronen-Rauigkeit (SD)",            getn("canopy_sd"),           2.5, "%.2f"),
    c("Deckungsgrad",                     getn("cover_frac"),          1.0, "%.2f")
  )

  ## ---- Plot -------------------------------------------------------
  png(file.path("testdata", sprintf("hedge_%s_%s.png", STAMP, id)),
      width = 1600, height = 850, res = 135)
  op <- par(mfrow = c(1,2), mar = c(4,4,3,1))

  ## (1) Seitenprofil der echten Punktwolke
  zmax <- max(4, ceiling(max(Z)))
  plot(along, Z, pch = 16, cex = 0.35, col = STRATA_COL[band_of(Z)],
       xlab = "Laenge entlang der Hecke (m)", ylab = "Hoehe ueber Grund (m)",
       ylim = c(0, zmax),
       main = sprintf("Hecke %s (Typ %s) – %s echte LiDAR-Punkte",
                      id, M$hedge_type, format(npoints(sub), big.mark=" ")))
  abline(h = c(0.5,1.5,3,6), col = "grey85", lty = 3)
  legend("topright", legend = STRATA_LAB, col = STRATA_COL, pch = 16,
         cex = 0.7, bg = "white", title = "Hoehenschicht")

  ## (2) Parameterprofil: wie hoch liegt die Hecke je Index-Groesse
  vals  <- vapply(params, function(p) as.numeric(p[[2]]), numeric(1))
  refs  <- vapply(params, function(p) as.numeric(p[[3]]), numeric(1))
  lvl   <- pmin(pmax(vals/refs, 0), 1)
  labs  <- vapply(params, function(p) p[[1]], character(1))
  fmt   <- vapply(params, function(p) p[[4]], character(1))
  txt   <- mapply(function(f,v) sprintf(f, v), fmt, vals)
  ramp  <- colorRamp(c("#f4e79e","#7cc08a","#14532d"))
  cols  <- apply(ramp(lvl), 1, function(r) rgb(r[1],r[2],r[3],maxColorValue=255))
  par(mar = c(4,11,3,2))
  bp <- barplot(rev(lvl), horiz = TRUE, col = rev(cols), border = NA,
                xlim = c(0,1.15), names.arg = rev(labs), las = 1,
                cex.names = 0.75, xlab = "Auspraegung (0 = niedrig, 1 = hoch)",
                main = "Strukturprofil je Indexparameter")
  text(rev(lvl), bp, labels = rev(txt), pos = 4, cex = 0.72, xpd = TRUE)
  mtext(sprintf("HVI gesamt: %.2f", getn("HVI")), side = 3, line = -1.2,
        adj = 0.98, cex = 0.85, font = 2, col = "#14532d")
  par(op); dev.off()

  ## auch ins RStudio-Fenster
  message(sprintf("Hecke %s: HVI %.2f, P95 %.1f m, %d Schichten, %s Punkte",
                  id, getn("HVI"), getn("h_p95"), getn("n_layers"),
                  format(npoints(sub), big.mark=" ")))
  message("geschrieben: testdata/hedge_", STAMP, "_", id, ".png")
  invisible(sub)
}
