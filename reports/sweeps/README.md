# Design-space sweeps

Default planar / subcutaneous / no-heparin grid from:

```bash
python3 t1d_simulator/screen_design.py --sweep \
  --geometry planar --site subcutaneous --no-heparin \
  --densities 20,50,80,100 --radii 100,150,200,250 \
  --export-csv reports/sweeps/planar_sq.csv \
  --export-json reports/sweeps/planar_sq.json \
  --heatmap reports/sweeps/heatmap_viability_planar_sq.png \
  --heatmap-po2 reports/sweeps/heatmap_po2_planar_sq.png
```

On this grid (16 designs): **0 PASS / 8 WARNING / 8 FAIL**. Thick, dense SQ slabs hit Krogh anoxia. Uncoated IBMIR retention is ~81% on every row (heparin is the IBMIR lever, not density).
