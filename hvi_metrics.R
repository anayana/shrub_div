# =====================================================================
# hvi_metrics.R
# Hedge Vertical Index (HVI) - core structural metrics from ALS LiDAR
# ---------------------------------------------------------------------
# Alle Metriken werden pro Hecken-Segment (bzw. pro Rasterzelle/
# Voxelsaeule) aus einer HOEHENNORMALISIERTEN Punktwolke berechnet.
# Abhaengigkeiten: nur base R fuer die Metriken selbst; lidR wird in
# hvi_pipeline.R fuer Einlesen/Normalisierung/Segmentierung benutzt.
# =====================================================================

# --- Standard-Hoehenschichten (Grenzen in m ueber Grund) -------------
# Angelehnt an die typische Heckenschichtung Mitteleuropas:
# Kraut/Streu, niedriger Strauch, mittlerer Strauch, hoher Strauch/
# Ueberhaelter, Baumschicht.
HVI_STRATA <- c(0, 0.5, 1.5, 3.0, 6.0, Inf)
HVI_STRATA_NAMES <- c("krautschicht", "niedrigstrauch",
                      "mittelstrauch", "hochstrauch", "baumschicht")

# --- Grundhoehe fuer Vegetation (Bodenpunkte/Rauschen ausblenden) ----
HVI_GROUND_CUTOFF <- 0.15  # m; Returns darunter = Boden

# ---------------------------------------------------------------------
# 1) Vertikales LAD-Profil nach MacArthur-Horn
# ---------------------------------------------------------------------
# Schaetzt die relative Blattflaechendichte (Leaf Area Density) pro
# Hoehenband aus dem vertikalen Extinktionsverlauf der Returns.
# LAD_i = ln( N(>=untere Grenze) / N(>=obere Grenze) )
# Rueckgabe: data.frame mit z (Bandmitte) und lad (rel. LAD, normiert).
# ---------------------------------------------------------------------
hvi_lad_profile <- function(z, dz = 0.5, z_max = NULL,
                            k = 1, min_returns = 5) {
  z <- z[is.finite(z) & z >= HVI_GROUND_CUTOFF]
  if (length(z) < min_returns) return(NULL)
  if (is.null(z_max)) z_max <- ceiling(max(z) / dz) * dz
  breaks <- seq(0, z_max, by = dz)
  mids   <- breaks[-length(breaks)] + dz / 2
  # kumulative Returns oberhalb jeder unteren Bandgrenze
  n_above <- sapply(breaks, function(b) sum(z >= b))
  # MacArthur-Horn: Verhaeltnis benachbarter kumulativer Zaehlungen
  # +1 gegen log(0); k = Extinktionskoeffizient
  lad <- (1 / k) * log((n_above[-length(n_above)] + 1) /
                       (n_above[-1] + 1))
  lad[!is.finite(lad) | lad < 0] <- 0
  data.frame(z = mids, lad = lad)
}

# ---------------------------------------------------------------------
# 2) Foliage Height Diversity (Shannon-Entropie des LAD-Profils)
# ---------------------------------------------------------------------
# FHD = -sum(p_i * ln p_i), p_i = Anteil des LAD in Band i.
# Klassiker: MacArthur & MacArthur (1961). Hoch = Blattmasse
# gleichmaessig ueber viele Hoehen verteilt.
# ---------------------------------------------------------------------
hvi_fhd <- function(lad) {
  p <- lad[lad > 0]
  if (length(p) == 0) return(0)
  p <- p / sum(p)
  -sum(p * log(p))
}

# ---------------------------------------------------------------------
# 3) Vertical Complexity Index (auf [0,1] normierte FHD)
# ---------------------------------------------------------------------
# VCI = FHD / ln(n_bins). Macht FHD zwischen Segmenten mit
# unterschiedlicher Hoehe/Bandzahl vergleichbar (van Ewijk et al. 2011).
# ---------------------------------------------------------------------
hvi_vci <- function(lad) {
  n <- length(lad[lad >= 0])
  if (n < 2) return(0)
  fhd <- hvi_fhd(lad)
  fhd / log(n)
}

# ---------------------------------------------------------------------
# 4) Schichtbelegung & Schichtungs-Evenness
# ---------------------------------------------------------------------
# Anteil der Returns je Standardschicht -> Zahl belegter Schichten
# und Pielou-Evenness der Schichtbelegung. Trennt "geschichtet
# strukturreich" von "monolithisch dicht".
# ---------------------------------------------------------------------
hvi_strata <- function(z, strata = HVI_STRATA,
                       names = HVI_STRATA_NAMES, occ_thr = 0.02) {
  z <- z[is.finite(z) & z >= HVI_GROUND_CUTOFF]
  if (length(z) == 0)
    return(list(occ_frac = setNames(rep(0, length(names)), names),
                n_layers = 0, evenness = 0))
  cnt <- table(cut(z, breaks = strata, right = FALSE, labels = names))
  frac <- as.numeric(cnt) / sum(cnt)
  names(frac) <- names
  n_layers <- sum(frac >= occ_thr)
  p <- frac[frac > 0]
  H <- if (length(p) > 1) -sum(p * log(p)) else 0
  evenness <- if (length(p) > 1) H / log(length(p)) else 0
  list(occ_frac = frac, n_layers = n_layers, evenness = evenness)
}

# ---------------------------------------------------------------------
# 5) Vertikale Lueckigkeit / Durchlaessigkeit (Porosity)
# ---------------------------------------------------------------------
# understory_openness = Anteil Returns < 1 m an allen Vegetations-
#   returns (offener Unterwuchs -> Nahrungssuche am Boden).
# vertical_gap = 1 - Anteil belegter dz-Baender bis P95 (vertikale
#   Lueckigkeit; hoch = zerrissenes Profil).
# ---------------------------------------------------------------------
hvi_porosity <- function(z, dz = 0.5) {
  z <- z[is.finite(z) & z >= HVI_GROUND_CUTOFF]
  if (length(z) < 3)
    return(list(understory_openness = NA, vertical_gap = NA))
  p95 <- as.numeric(quantile(z, 0.95))
  understory_openness <- mean(z < 1.0)
  breaks <- seq(0, ceiling(p95 / dz) * dz, by = dz)
  if (length(breaks) < 2)
    return(list(understory_openness = understory_openness,
                vertical_gap = 0))
  occ <- table(cut(z, breaks = breaks, right = FALSE))
  vertical_gap <- mean(occ == 0)
  list(understory_openness = understory_openness,
       vertical_gap = vertical_gap)
}

# ---------------------------------------------------------------------
# 6) Hoehen-Kennwerte
# ---------------------------------------------------------------------
hvi_height_stats <- function(z) {
  z <- z[is.finite(z) & z >= HVI_GROUND_CUTOFF]
  if (length(z) == 0)
    return(list(h_max = NA, h_mean = NA, h_p95 = NA, h_sd = NA))
  list(h_max  = max(z),
       h_mean = mean(z),
       h_p95  = as.numeric(quantile(z, 0.95)),
       h_sd   = sd(z))
}

# ---------------------------------------------------------------------
# 7) Horizontale Heterogenitaet (Rugositaet des Kronenprofils)
# ---------------------------------------------------------------------
# Aus einem CHM-Ausschnitt (Vektor der Zellhoehen des Segments):
# SD und Variationskoeffizient der Oberhoehe + Lueckenanteil.
# Bester Einzel-Praediktor fuer Gesamtdiversitaet in Heckenstudien.
# ---------------------------------------------------------------------
hvi_horizontal <- function(chm_vals, gap_height = 0.5) {
  v <- chm_vals[is.finite(chm_vals)]
  if (length(v) < 4)
    return(list(canopy_sd = NA, canopy_cv = NA, cover_frac = NA))
  cover <- v >= gap_height
  list(canopy_sd   = sd(v),
       canopy_cv   = sd(v) / (mean(v) + 1e-6),
       cover_frac  = mean(cover))
}
