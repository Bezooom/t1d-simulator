# Benchmark: VEGF/Angiogenesis (Berney 2016)

## Paper
- Authors: Berney, Holdsworth, Kulkarni, et al.
- Year: 2016
- Journal: Nature Biomedical Engineering, 1:0026
- DOI: 10.1038/s41551-016-0026-6

## Prediction
- Method: solve_coupled_neovascularization(species="Mouse"/"Human")
- Input params:
  - VEGF threshold = 50.0 pg/mL
  - Angiogenesis rate = 0.055 1/day
  - K_vegf = 0.1
  - Organoid radius = 150 μm
  - Cell density = 1e+08 cells/cm³
- Predicted time to vascularization (mouse): 5.6 days
- Predicted time to vascularization (human): 8.8 days

## Reported
- Reported time to vascularization (mouse): 7.0 days (range 7-14)
- Reported time to vascularization (human): 31.5 days (range 21-42)
- Source: "Vascularization occurred within 10 days in mice and ~30 days in humans."

## Agreement
- Mouse: ratio = 0.80 (✅ PASS)
- Human: ratio = 0.28 (❌ FAIL)
- Notes: Human scaling factor of 2.5x applied to angiogenesis rate (per Berney 2016).

## Figure
- Path: /home/bezoom/storage/Projects/Диабет/reports/benchmarks/fig_vegf.png
