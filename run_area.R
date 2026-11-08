# =====================================================================
# run_area.R  --  NUR einen kleinen Bereich neu rechnen (schnell)
# ---------------------------------------------------------------------
# Rechnet HVI + Eignung + exportiert die Punktwolken fuer eine kleine
# AOI in einem Rutsch. 300 m = wenige Hecken, Sekunden statt Stunden.
#
# In RStudio:
#   source("run_area.R")
#   run_area("C_28FN2.LAZ", 255050, 494259, 300)   # repraesentativer Ausschnitt
#
# Danach: mir Bescheid geben -> ich baue das Dashboard aus dem neuen,
# kleinen (und hoehen-korrekten) Ergebnis neu zusammen.
# =====================================================================

source("run_real.R")
source("export_points.R")

run_area <- function(TILE = "C_28FN2.LAZ", xmin, ymin, size = 300) {
  run_real(TILE, xmin, ymin, size)     # -> hedges_/results_/map_ (nur diese AOI)
  export_points(TILE)                  # -> points_ (nur diese Hecken, schnell)
  message("\nFertig – nur der Bereich ", xmin, ",", ymin, " (", size, " m) neu gerechnet.")
}
