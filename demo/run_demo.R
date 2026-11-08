# =====================================================================
# run_demo.R
# HVI-Demo OHNE echte LiDAR-Daten: erzeugt synthetische Hecken-
# Punktwolken verschiedener Struktur-Typen und rechnet den kompletten
# Index + artspezifische Eignung durch.
#   Ausfuehren:  Rscript demo/run_demo.R   (aus dem Projektordner)
# Benoetigt nur base R.
# =====================================================================

source(file.path("R", "hvi_metrics.R"))
source(file.path("R", "hvi_index.R"))

set.seed(42)

# --- Synthetische Hecke: Mischung von Hoehenschichten ----------------
# layers: Liste aus (mean_height, sd, n_returns) je Schicht.
# chm_rough: Rauigkeit der Oberhoehe (horizontale Heterogenitaet).
make_hedge <- function(layers, chm_top, chm_rough, n_chm = 200, ...) {
  z <- unlist(lapply(layers, function(L)
    pmax(0, rnorm(L[3], L[1], L[2]))))
  chm <- pmax(0, rnorm(n_chm, chm_top, chm_rough))
  list(z = z, chm = chm)
}

# --- Vier idealtypische Heckenstrukturen -----------------------------
hedges <- list(
  # 1: niedrige, dicht geschorene Hecke (monolithisch, eine Schicht)
  "geschoren_niedrig" = make_hedge(
    layers = list(c(0.8, 0.25, 1500)),
    chm_top = 1.0, chm_rough = 0.15, ndvi = 0.5),
  # 2: mehrschichtige strukturreiche Hecke mit Krautsaum
  "strukturreich" = make_hedge(
    layers = list(c(0.3,0.15,400), c(1.2,0.4,600),
                  c(2.5,0.6,600), c(4.5,0.9,400)),
    chm_top = 4.8, chm_rough = 1.4, ndvi = 0.6),
  # 3: hohe, lueckige Baumreihe (leerer Unterwuchs)
  "baumreihe_hoch" = make_hedge(
    layers = list(c(7.5, 1.5, 1200)),
    chm_top = 8.0, chm_rough = 2.2, ndvi = 0.65),
  # 4: mittlere dornige Hecke mit offenem Saum (Neuntoeter-Typ)
  "dornig_mittel" = make_hedge(
    layers = list(c(0.3,0.15,300), c(1.8,0.5,900), c(2.8,0.5,500)),
    chm_top = 2.6, chm_rough = 0.7, ndvi = 0.55)
)
ndvi_vals <- c(0.5, 0.6, 0.65, 0.55)

# --- Metriktabelle aufbauen ------------------------------------------
M <- do.call(rbind, Map(function(h, nd) {
  as.data.frame(t(hvi_segment_metrics(h$z, chm_vals = h$chm,
                                      ndvi = nd)))
}, hedges, ndvi_vals))
M <- cbind(hedge = names(hedges), M)

cat("\n== Roh-Metriken je Hecke ==\n")
Mp <- M[, c("hedge","fhd","vci","n_layers","strata_evenness",
            "understory_openness","vertical_gap","h_p95",
            "canopy_sd","cover_frac")]
num <- sapply(Mp, is.numeric)
Mp[num] <- round(Mp[num], 3)
print(Mp, row.names = FALSE)

# --- Zusammengesetzter HVI -------------------------------------------
idx <- hvi_compute(M)
res <- cbind(hedge = M$hedge, round(idx, 3))
cat("\n== Subindizes + HVI ==\n")
print(res, row.names = FALSE)

# --- Habitat-Typ (Cluster) -------------------------------------------
# bei nur 4 Segmenten k=2 zur Demonstration
res$type <- hvi_classify(idx, k = 2)

# --- Artspezifische Eignung ------------------------------------------
req <- read.csv(file.path("data", "species_requirements.csv"),
                stringsAsFactors = FALSE)
hsi <- hvi_species_suitability(M, req)
hsi <- cbind(hedge = M$hedge, round(hsi, 2))
cat("\n== Habitateignung je Art (0..1) ==\n")
print(hsi, row.names = FALSE)

cat("\nInterpretation: 'strukturreich' sollte den hoechsten HVI und",
    "\nhohe Vogel-Eignung zeigen; 'baumreihe_hoch' punktet bei der",
    "\nFledermaus-Leitfunktion; 'dornig_mittel' beim Neuntoeter;",
    "\n'geschoren_niedrig' bei Bestaeubern (offener Saum).\n")
