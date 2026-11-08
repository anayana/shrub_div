# =====================================================================
# hvi_index.R
# Zusammengesetzter Hedge Vertical Index + artspezifische Eignung
# ---------------------------------------------------------------------
# Baut auf den Roh-Metriken aus hvi_metrics.R auf. Zwei Ebenen:
#   (A) HVI  - datengetriebener Struktur-Komplexitaetsindex (0..1)
#   (B) HSI  - artspezifische Habitateignung ueber Fuzzy-Antwortkurven
# =====================================================================

# ---------------------------------------------------------------------
# Metrik-Tabelle je Hecken-Segment erzeugen
# ---------------------------------------------------------------------
# Erwartet: z (Vegetationshoehen der Punkte im Segment) und optional
# chm_vals (CHM-Zellhoehen des Segments) sowie ndvi (mittlere Greenness).
# Rueckgabe: benannter Vektor aller Kennwerte fuer EIN Segment.
# ---------------------------------------------------------------------
hvi_segment_metrics <- function(z, chm_vals = NULL, ndvi = NA,
                                dz = 0.5) {
  lad_df <- hvi_lad_profile(z, dz = dz)
  lad    <- if (is.null(lad_df)) numeric(0) else lad_df$lad
  st     <- hvi_strata(z)
  po     <- hvi_porosity(z, dz = dz)
  hs     <- hvi_height_stats(z)
  hz     <- if (!is.null(chm_vals)) hvi_horizontal(chm_vals)
            else list(canopy_sd = NA, canopy_cv = NA, cover_frac = NA)

  c(fhd                 = hvi_fhd(lad),
    vci                 = hvi_vci(lad),
    n_layers            = st$n_layers,
    strata_evenness     = st$evenness,
    understory_openness = po$understory_openness,
    vertical_gap        = po$vertical_gap,
    h_max               = hs$h_max,
    h_mean              = hs$h_mean,
    h_p95               = hs$h_p95,
    h_sd                = hs$h_sd,
    canopy_sd           = hz$canopy_sd,
    canopy_cv           = hz$canopy_cv,
    cover_frac          = hz$cover_frac,
    ndvi                = ndvi)
}

# ---------------------------------------------------------------------
# Robuste Min-Max-Normierung (Perzentil-getrimmt) auf [0,1]
# ---------------------------------------------------------------------
hvi_rescale <- function(x, lo = 0.02, hi = 0.98) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  qs <- quantile(x, c(lo, hi), na.rm = TRUE)
  if (!is.finite(qs[1]) || !is.finite(qs[2]) || diff(qs) == 0)
    return(rep(0.5, length(x)))
  y <- (x - qs[1]) / (qs[2] - qs[1])
  pmin(pmax(y, 0), 1)
}

# ---------------------------------------------------------------------
# (A) Zusammengesetzter HVI
# ---------------------------------------------------------------------
# Vier Subindizes, jeweils aus normierten Metriken gemittelt:
#   vertical_complexity : fhd, vci, n_layers, strata_evenness
#   volume_size         : h_p95, cover_frac
#   heterogeneity       : canopy_sd, canopy_cv, vertical_gap
#   condition           : ndvi  (optional; NA -> ignoriert)
# Gesamt-HVI = gewichtetes Mittel der Subindizes.
# Default-Gewichte sind literaturmotiviert, aber ueberschreibbar bzw.
# durch PCA (siehe hvi_pca_weights) ersetzbar.
# ---------------------------------------------------------------------
# Default: vertikale Komplexitaet dominiert (es ist ein *Vertikal*-Index).
# Im Test kippte 0.40/0.20/0.30/0.10 die Rangfolge zugunsten hoher, aber
# einschichtiger Baumreihen; 0.50/0.15/0.25/0.10 stellt die mehrschichtige
# Hecke wieder an die Spitze. Gewichte sind gegen Felddaten kalibrierbar.
hvi_default_weights <- c(vertical_complexity = 0.50,
                         volume_size         = 0.15,
                         heterogeneity       = 0.25,
                         condition           = 0.10)

hvi_compute <- function(M, weights = hvi_default_weights) {
  M <- as.data.frame(M)
  R <- data.frame(
    fhd             = hvi_rescale(M$fhd),
    vci             = M$vci,                       # bereits [0,1]
    n_layers        = hvi_rescale(M$n_layers),
    strata_evenness = M$strata_evenness,           # bereits [0,1]
    h_p95           = hvi_rescale(M$h_p95),
    cover_frac      = ifelse(is.na(M$cover_frac), NA, M$cover_frac),
    canopy_sd       = hvi_rescale(M$canopy_sd),
    canopy_cv       = hvi_rescale(M$canopy_cv),
    vertical_gap    = hvi_rescale(M$vertical_gap),
    ndvi            = hvi_rescale(M$ndvi)
  )
  rmean <- function(...) rowMeans(cbind(...), na.rm = TRUE)
  sub <- data.frame(
    vertical_complexity = rmean(R$fhd, R$vci, R$n_layers,
                                R$strata_evenness),
    volume_size         = rmean(R$h_p95, R$cover_frac),
    heterogeneity       = rmean(R$canopy_sd, R$canopy_cv, R$vertical_gap),
    condition           = R$ndvi
  )
  w <- weights[colnames(sub)]
  w[is.na(w)] <- 0
  # Subindizes, die komplett NA sind (z.B. keine Spektraldaten),
  # aus Gewichtung nehmen
  valid <- colSums(!is.na(sub)) > 0
  w[!valid] <- 0
  w <- w / sum(w)
  hvi <- as.numeric(as.matrix(
    replace(sub, is.na(sub), 0)) %*% w)
  cbind(sub, HVI = hvi)
}

# PCA-basierte Gewichte (erste Achse) als datengetriebene Alternative
hvi_pca_weights <- function(M) {
  cols <- c("fhd", "vci", "n_layers", "strata_evenness",
            "h_p95", "canopy_sd", "vertical_gap")
  X <- as.matrix(M[, cols])
  X <- X[complete.cases(X), , drop = FALSE]
  pc <- prcomp(X, scale. = TRUE)
  load <- abs(pc$rotation[, 1])
  load / sum(load)
}

# ---------------------------------------------------------------------
# (B) Artspezifische Habitateignung (Fuzzy-Antwortkurven)
# ---------------------------------------------------------------------
# Jede Art wird durch eine Tabelle von Antwortkurven auf Kennwerte
# beschrieben. Kurventypen:
#   trapz(a,b,c,d)   Trapez: 0 bis a, linear bis 1 in [b,c], 0 ab d
#   gauss(mu,sigma)  Glockenkurve, Optimum bei mu
#   mono_inc(a,b)    monoton steigend, 0 bei a, 1 ab b
#   mono_dec(a,b)    monoton fallend, 1 bis a, 0 ab b
# Gesamteignung = geometrisches Mittel der Teil-Zugehoerigkeiten
# (Liebig-Prinzip: ein fehlender Schluesselfaktor zieht Eignung stark).
# ---------------------------------------------------------------------
hvi_membership <- function(x, type, p) {
  if (is.na(x)) return(NA)
  switch(type,
    trapz = {
      a<-p[1]; b<-p[2]; c<-p[3]; d<-p[4]
      if (x <= a || x >= d) 0
      else if (x < b) (x - a) / (b - a)
      else if (x <= c) 1
      else (d - x) / (d - c)
    },
    gauss = exp(-((x - p[1])^2) / (2 * p[2]^2)),
    mono_inc = { a<-p[1]; b<-p[2]
      if (x <= a) 0 else if (x >= b) 1 else (x - a) / (b - a) },
    mono_dec = { a<-p[1]; b<-p[2]
      if (x <= a) 1 else if (x >= b) 0 else (b - x) / (b - a) },
    stop("unbekannter Kurventyp: ", type)
  )
}

# requirements: data.frame mit Spalten
#   species, metric, type, p1, p2, p3, p4, weight
hvi_species_suitability <- function(M, requirements) {
  M <- as.data.frame(M)
  species <- unique(requirements$species)
  out <- matrix(NA_real_, nrow = nrow(M), ncol = length(species),
                dimnames = list(NULL, species))
  for (sp in species) {
    req <- requirements[requirements$species == sp, ]
    memb <- matrix(NA_real_, nrow(M), nrow(req))
    for (j in seq_len(nrow(req))) {
      pars <- as.numeric(req[j, c("p1","p2","p3","p4")])
      pars <- pars[!is.na(pars)]
      memb[, j] <- vapply(M[[req$metric[j]]],
                          hvi_membership, numeric(1),
                          type = req$type[j], p = pars)
    }
    w <- req$weight; w <- w / sum(w)
    # gewichtetes geometrisches Mittel (Liebig)
    lg <- log(pmax(memb, 1e-6))
    out[, sp] <- exp(as.numeric(
      replace(lg, is.na(lg), 0) %*%
      ifelse(is.na(memb[1, ]), 0, w)))
  }
  as.data.frame(out)
}

# ---------------------------------------------------------------------
# Habitat-Klassifikation: Segmente in Struktur-Typen clustern
# ---------------------------------------------------------------------
hvi_classify <- function(sub, k = 4, seed = 1) {
  n <- nrow(sub)
  k <- max(1, min(k, n - 1))            # k muss < Zahl der Hecken sein
  if (n < 2) return(factor(rep(1L, n)))
  X <- scale(sub[, c("vertical_complexity","volume_size","heterogeneity")])
  X[is.na(X)] <- 0
  set.seed(seed)
  km <- kmeans(X, centers = k, nstart = 25)
  factor(km$cluster)
}
