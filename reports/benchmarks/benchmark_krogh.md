# Benchmark: Krogh Limit (Secomb 2004)

## Paper
- Authors: Secomb, Hsu, Poon, Kosky, Nelson
- Year: 2004
- Journal: Microcirculation, 11(2):105-117
- DOI: 10.1080/713692626

## Prediction
- Method: solve_oxygen_profile(geometry="planar")
- Input params:
  - V_max = 1.50e-16 mol/(cell·s) (Buchwald 2011, human β-cell)
  - K_M = 0.5 mmHg
  - D_O2_tissue = 2e-05 cm²/s (Secomb 2004)
  - pO₂_arterial = 95.0 mmHg
  - Cell density = 1e+08 cells/cm³
- Predicted Krogh limit: 200 μm

## Reported
- Reported Krogh limit: 200 μm (Secomb et al., 2004)
- Source: "The capillary spacing in skeletal muscle is ~200 μm, corresponding to the Krogh limit."

## Agreement
- Ratio: predicted/reported = 1.00
- Status: ✅ PASS
- Notes: Planar geometry used as Krogh cylinder approximation. Human β-cell V_max from Buchwald 2011.

## Figure
- Path: /home/bezoom/storage/Projects/Диабет/reports/benchmarks/fig_krogh_limit.png
