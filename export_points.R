# =====================================================================
# export_points.R  --  DICHTE 3D-Punktwolke je Hecke fuers Dashboard
# ---------------------------------------------------------------------
# Behaelt viele Punkte (bis Nmax) in echtem 3D: a = Laengsachse der
# Hecke, b = quer, z = Hoehe. Das Dashboard zeichnet daraus eine
# drehbare Punktwolke. Gemeinsame Normalisierung ueber den ganzen
# Heckenbereich (wie run_real) -> Hoehen korrekt.
# Voraussetzung: run_real() lief (testdata/hedges_<tile>.gpkg).
#
# In RStudio:
#   source("export_points.R")
#   export_points("C_28FN2.LAZ")            # alle Hecken im GPKG
#   export_points("C_28FN2.LAZ", Nmax=8000) # mehr Punkte je Hecke
# =====================================================================

suppressMessages({library(lidR); library(sf)})

export_points <- function(TILE = "C_28FN2.LAZ", Nmax = 100000, buffer = 3) {
  STAMP <- sub("\\.LAZ$", "", TILE, ignore.case = TRUE)
  LAZ   <- file.path("testdata", TILE)
  gpkg  <- file.path("testdata", paste0("hedges_", STAMP, ".gpkg"))
  stopifnot(file.exists(LAZ), file.exists(gpkg))
  hedges <- sf::st_read(gpkg, quiet = TRUE)

  ## einmal gemeinsam einlesen + normalisieren (korrekte Hoehen)
  bb  <- as.numeric(sf::st_bbox(hedges))
  las <- readLAS(LAZ, select = "xyzc",
      filter = sprintf("-keep_xy %f %f %f %f -drop_z_below -1 -drop_z_above 60",
                       bb[1]-buffer, bb[2]-buffer, bb[3]+buffer, bb[4]+buffer))
  if (is.empty(las)) stop("keine Punkte im Heckenbereich")
  if (is.na(crs(las)) || crs(las) == "") crs(las) <- "EPSG:28992"
  las <- normalize_height(las, tin())
  las <- filter_poi(las, Z >= 0, Z < 40)

  parts <- character(0); tot <- 0
  for (i in seq_len(nrow(hedges))) {
    h   <- hedges[i, ]
    bh  <- as.numeric(sf::st_bbox(h))          # ganze Hecke, nicht nur CHM-Kontur
    sub <- clip_rectangle(las, bh[1]-1, bh[2]-1, bh[3]+1, bh[4]+1)
    if (is.empty(sub) || npoints(sub) < 10) next
    X <- sub$X; Y <- sub$Y; Z <- sub$Z
    if (length(X) > Nmax) { k <- sample(length(X), Nmax); X<-X[k]; Y<-Y[k]; Z<-Z[k] }
    pc <- prcomp(cbind(X, Y))
    a  <- pc$x[, 1] - min(pc$x[, 1])     # Laengsachse (0-basiert)
    b  <- pc$x[, 2]                       # Querachse (zentriert)
    tot <- tot + length(Z)
    coords <- paste(sprintf("[%.2f,%.2f,%.1f]", a, b, round(Z,1)), collapse = ",")
    parts <- c(parts, sprintf('"%s":[%s]', h$hedge_id, coords))
  }
  js  <- paste0("window.HEDGE_POINTS={", paste(parts, collapse = ","), "};")
  out <- file.path("testdata", paste0("points_", STAMP, ".js"))
  writeLines(js, out)
  message("geschrieben: ", out, "  (", length(parts), " Hecken, ",
          format(tot, big.mark=" "), " Punkte gesamt)")
  invisible(out)
}
