# T1D Digital Twin & Hypoimmune Organoid Simulator

A multiphysics **in silico** engine for designing and screening cell replacement therapies for Type 1 Diabetes (T1D). Models macroencapsulation failure modes — hypoxia, fibrosis, immune activation — and hypoimmune organoid dynamics with IBMIR kinetics, site selection, and dose calculation.

## What the Model Does

- **O₂ PDE solver** — 1D reaction-diffusion with Michaelis–Menten kinetics across planar, cylindrical, and spherical geometries; couples FBR, ROS, pH, and gel swelling.
- **PINN surrogate** — mesh-free physics-informed neural network solving the same boundary-value problems at ~100× speed.
- **GNN biocompatibility** — PyTorch Geometric graph network mapping polymer SMILES to anti-fibrotic scores ($L_{fib}$).
- **Organoid simulator** — 0–48 h IBMIR kinetics, gene-edit sets ($B2M^{-/-}$, $CIITA^{-/-}$, $CD47^{KI}$, $PD-L1^{KI}$, $CD55/CD59^{KI}$), iCasp9/apoptosis, and site selection matrix (portal vs. omental vs. subcutaneous).
- **Interactive UI** — Streamlit app for parameter exploration, CGM metrics (TIR/TBR/TAR/HbA1c), and CAD scaffold export.

## Repository Structure

```
t1d_simulator/
├── simulator.py              # Core multiphysics PDE + neovascularization
├── organoid_simulator.py     # Organoid dynamics, IBMIR, CGM, dosing
├── pinn_solver.py            # PINN surrogate (torch)
├── gnn_pipeline.py           # SMILES → biocompatibility score
├── mesh_generator.py         # TPMS Gyroid/Schwarz meshes
├── organoid_cad_exporter.py  # STL export + clinical passport
├── param_loader.py           # YAML config loader
├── verify_model.py           # 35 integration tests (standalone + pytest)
├── app.py                    # Streamlit web interface
├── aid_controller.py         # Closed-loop insulin pump (PID)
├── parameters.yaml           # Production runtime parameters
├── requirements.txt          # Runtime dependencies
├── pyproject.toml            # Project metadata + tool configs
├── biocompatibility_gnn.pt   # Trained GNN weights (N=46)
└── config/                   # Configuration subpackage
```

## Installation

```bash
# Option 1: pip install from repo
pip install -e .

# Option 2: via requirements.txt
pip install -r t1d_simulator/requirements.txt
```

Optional dev dependencies: `pip install -e ".[dev]"` for pytest + ruff.

## Quickstart

```bash
# 1. Verify the model (35 tests, ~30 s)
python3 t1d_simulator/verify_model.py

# 2. Launch the interactive UI
streamlit run t1d_simulator/app.py

# 3. Run with ruff (lint check)
ruff check t1d_simulator/
```

## Configuration

Parameters are loaded via `param_loader.py` from `parameters.yaml` (rooted at `t1d_simulator/parameters.yaml`) and cross-referenced with `data/literature_params.yaml` in the repo root.

Key calibrated values:

| Parameter | Value | Source |
|-----------|-------|--------|
| O₂ solubility ($S_{O2}$) | $1.34 \times 10^{-9}$ mol/(cm³·mmHg) | — |
| V$_{max}$ (mouse β-cell) | $1.2 \times 10^{-16}$ mol/(cell·s) | Dionne et al., 1993 |
| V$_{max}$ (human β-cell) | $1.5 \times 10^{-16}$ mol/(cell·s) | Buchwald et al., 2011 |
| K$_M$ (viability) | 0.5 mmHg | — |
| K$_M$ (insulin secretion) | 5.0 mmHg | — |

## Architecture

```
Parameters (YAML)
       │
       ▼
┌─────────────┐    O₂ / cytokine / VEGF    ┌──────────────┐
│  simulator   │◄─────────────────────────►│  organoid     │
│  (PDE solver)│                             │  simulator    │
└──────┬──────┘                             └──────┬───────┘
       │                                            │
       ▼                                            ▼
  pinn_solver.py                              gnn_pipeline.py
  (fast surrogate)                            (SMILES → score)
                                                    │
                                                    ▼
                                              mesh_generator.py
                                              (3D scaffold export)
```

Data flows: parameters → `simulator.py` (PDE) ↔ `organoid_simulator.py` (IBMIR/CGM) ↔ `gnn_pipeline.py` (materials) → `mesh_generator.py` (geometry) → STL/CAD.

## Scientific Validation

| Benchmark | Status | Key Result |
|-----------|--------|------------|
| Krogh limit (R=200 µm, 100M cells/ml) | ✅ PASS | Sphere > Cylinder > Planar viability |
| IBMIR 0–48 h kinetics | ✅ PASS | Heparinization increases 48 h retention |
| iCasp9/apoptosis (AP1903) | ✅ PASS | >95% elimination in 4 h |
| CGM 30-day (transplanted) | ✅ PASS | TIR >90%, TBR=0%, HbA1c <6.0% |

## Limitations

- **In silico only** — failure mode screening and design-space exploration; does not replace in vivo studies.
- **Simplified fluid dynamics** — IBMIR uses lumped ODE kinetics, not full 3D CFD.
- **GNN dataset** — trained on N≈46 curated polymers; novel chemistries need experimental validation.
- **1D geometry** — PDE solvers are 1D; 3D effects captured via spherical/cylindrical corrections.
- **No wet-lab calibration yet** — parameters are literature-calibrated but not jointly fitted to a single experimental dataset.

## Roadmap

Full roadmap: [ROADMAP.md](../ROADMAP.md)

| Milestone | Status | Target |
|-----------|--------|--------|
| M0 Engineering hygiene | ✅ Done | — |
| M1 Literature calibration | ✅ Done | — |
| M2 Open science release | 🔄 In progress | Preprint + Zenodo |
| M3 First external signal | — | Q1 2027 |
| M4 Joint science | — | Q2 2027 |

## Citation

If you use this model or code in your research:

```bibtex
@software{t1d_digital_twin_2026,
  title = {T1D Digital Twin & Hypoimmune Organoid Simulator},
  author = {Author Name},
  year = {2026},
  url = {https://github.com/your-username/Диабет},
  version = {1.0.0},
  license = {MIT}
}
```

## License

MIT — see [LICENSE](../LICENSE) for details. Academic co-development welcome.
