# Hedge Vertical Index (HVI)

Zusammengesetzter LiDAR-Strukturindex zur Beschreibung,
Habitat-Klassifikation und artspezifischen Bewertung von Hecken.

## Inhalt

    HVI_Konzept.docx          Methodenkonzept (Motivation, Novelty, Formeln, Validierung)
    R/hvi_metrics.R           Roh-Metriken je Segment (base R): LAD-Profil, FHD, VCI,
                              Schichtbelegung, Lueckigkeit, Hoehen, horizontale Heterogenitaet
    R/hvi_index.R             Composite-HVI, PCA-Gewichte, Fuzzy-Suitability, Habitat-Cluster
    R/hvi_pipeline.R          End-to-End mit echten ALS-Daten (lidR/sf/terra)
    data/species_requirements.csv   Artanforderungen als Antwortkurven (kalibrierbar)
    demo/run_demo.R           Demo ohne echte Daten (synthetische Hecken), nur base R

## Schnellstart (Demo, ohne LiDAR-Daten)

    cd shrub_div
    Rscript demo/run_demo.R

Erzeugt vier idealtypische Heckenstrukturen und rechnet Metriken, HVI,
Habitattyp und artspezifische Eignung durch.

## Echte Daten

1. Freie ALS-Punktwolke laden (z. B. opengeodata.nrw.de, Bayern OpenData).
2. Heckenpolygone bereitstellen (`sf`, Spalte `hedge_id`) - z. B. aus dem
   Biotopkataster oder per CHM-Schwelle.
3. Optional NDVI-Raster (Sentinel-2 / DOP-CIR) als `terra`-Raster.

```r
source("R/hvi_pipeline.R")
req <- read.csv("data/species_requirements.csv")
out <- hvi_run("kachel.laz", hedges, ndvi_rast = ndvi, requirements = req)
out$result   # HVI, Subindizes, Habitattyp und Arteignungen je Segment
```

## Index in Kurzform

- **Vertikale Komplexitaet** (Gewicht 0,50): FHD, VCI, Schichtzahl, Schicht-Evenness
- **Volumen/Groesse** (0,15): P95-Hoehe, Deckungsgrad
- **Heterogenitaet** (0,25): Kronen-SD, Kronen-CV, vertikale Luecke
- **Zustand** (0,10, optional): NDVI

HVI = gewichtetes Mittel der auf [0,1] normierten Subindizes.
Gewichte sind gegen Felddaten kalibrierbar (siehe `hvi_pca_weights`).

## Abhaengigkeiten

- Demo/Metriken: nur base R
- Echte Daten: `lidR`, `sf`, `terra`

## Status

Prototyp. Die Index-Mathematik wurde per Python-Portierung mit synthetischem
Szenario verifiziert (plausible Rangfolge; zwei Fehler korrigiert - Details in
Abschnitt 8 des Konzeptpapiers). Antwortkurven sind Startwerte und muessen
gegen reale Vorkommensdaten kalibriert werden.
