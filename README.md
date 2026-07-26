# T1D Digital Twin & Hypoimmune Organoid Simulator

[![Tests](https://github.com/Bezooom/t1d-simulator/actions/workflows/ci.yml/badge.svg)](https://github.com/Bezooom/t1d-simulator/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-Preprint-red.svg)](docs/manuscript_biorxiv_en.md)
[![Zenodo](https://img.shields.io/badge/Zenodo-Archive-blue.svg)](zenodo.json)

An open-source, literature-calibrated multiphysics digital twin engine (`t1d_simulator`) developed by **Pavel V. Naumov** to simulate, evaluate, and de-risk cell replacement therapies for Type 1 Diabetes (T1D). The platform maps the physical, biochemical, and immunological design space of macroencapsulated stem cell-derived $\beta$-cell grafts and hypoimmune organoid implants.

![Graphical Abstract](docs/figures/graphical_abstract.png)

---

## 🔮 Core Features

*   **Oxygen Transport PDE Solver:** Solves the non-linear reaction-diffusion equations under spherical, cylindrical, and planar geometries to identify diffusion barriers and necrotic core formation ($pO_2 < 0.5\ \mathrm{mmHg}$).
*   **Physics-Informed Neural Network (PINN):** A mesh-free neural solver surrogate that computes steady-state oxygen concentration profiles in under 10 milliseconds, enabling real-time parameter sweeps.
*   **IBMIR Coagulation ODE Cascade:** A compartmental model of the 0–48h Instant Blood-Mediated Inflammatory Reaction (*Tissue Factor → Thrombin → Fibrin Clot*) which dynamically reduces oxygen permeability at the boundary.
*   **GNN Anti-Fibrotic Screening:** A Graph Neural Network that maps chemically modified alginate monomers (SMILES) to foreign body reaction (FBR) fibrotic capsule thickness ($L_{\mathrm{fib}}$) predictions.
*   **Clinical Dose Calculator:** Translates virtual patient mass and daily insulin requirements (TDI) into target islet equivalents (IEQ), organoid counts, and omental scaffold volume.
*   **Bergman OGTT Glycemic Dynamics:** Simulates glucose-stimulated insulin secretion (GSIS) and glucose clearance curves during oral glucose tolerance tests (OGTT).
*   **CGM 30-Day Monitoring:** Simulates long-term continuous glucose metrics, returning GMI (HbA1c), Time-in-Range (TIR), Time-Below-Range (TBR), and Time-Above-Range (TAR).
*   **CAD/STL Scaffold Exporter:** Procedurally generates and exports 3D porous TPMS scaffolds (STL format) matching patient-specific omental pouch area requirements.

---

## 📂 Repository Structure

```
.
├── .github/workflows/          # CI/CD pipelines
├── data/
│   ├── literature_params.yaml  # Curated parameters with literature citations
│   └── benchmarks.py           # Verification script for validation datasets
├── docs/
│   ├── figures/                # Manuscripts schematics and graphical abstract
│   ├── manuscript_biorxiv_en.md# Primary English bioRxiv preprint (v1.2)
│   ├── manuscript_biorxiv_ru.md# Russian parallel preprint draft (v1.2)
│   ├── site_comparison_matrix.md# Comparative matrix of transplant locations
│   └── *.md                    # Clinical passports, MoU templates, protocols
├── reports/
│   └── benchmarks/             # Target benchmark figures and reports
├── scripts/                    # Helper scripts for screening and email templates
├── t1d_simulator/              # Python simulator source package
│   ├── app.py                  # Main Streamlit dashboard (monolith)
│   ├── ibmir_module.py         # IBMIR kinetics solver
│   ├── verify_model.py         # Integration test suite (40 automated tests)
│   └── parameters.yaml         # Default runtime parameters
└── pyproject.toml              # Build configurations and dependency lists
```

---

## 🚀 Quickstart

### 1. Installation
The package runs on Python 3.9–3.12. Clone the repository and install dependencies:

```bash
git clone https://github.com/Bezooom/t1d-simulator.git
cd t1d-simulator
pip install -r t1d_simulator/requirements.txt
```

### 2. Run Interactive Dashboard
Launch the Streamlit graphical user interface:

```bash
streamlit run t1d_simulator/app.py
```

### 3. Verify Model and Run Tests
Run the 40-test integration test suite on CPU:

```bash
export CUDA_VISIBLE_DEVICES=""
python3 t1d_simulator/verify_model.py
```

### 4. Execute Construct Screening
Run a CLI construct design screen:

```bash
python3 t1d_simulator/screen_design.py
```

---

## 📊 Scientific Validation

The digital twin's transport and cell survival solvers have been benchmarked and validated against independent published datasets:

| Benchmark | Reference Paper | Target Metric | Agreement Status |
|-----------|-----------------|---------------|------------------|
| **Krogh Limit** | Secomb et al. (2004) | Viability boundary limit (~200 $\mu\text{m}$) | **PASS** (1.00 ratio) |
| **Spherical Alginate** | Papas et al. (2007) | Viability sweep ($R = 400\ \mu\text{m}$) | **PASS** (RMSE = 10.42%) |
| **Planar Macrodevice** | Papabathini et al. (2023) | Hypoxic core profile ($L = 250\ \mu\text{m}$) | **PASS** (RMSE = 11.42%) |
| **IBMIR Time-Course** | Hackett et al. (2013) | 48h acute cell retention | **PASS** (1.03 ratio) |
| **Angiogenesis** | Berney et al. (2016) | Capillary perfusion time | **PASS** (0.80 ratio) |

---

## 📑 Citation & Academic Use

If you use this simulator or reference the preprint findings in your research, please cite the work as follows:

```bibtex
@article{naumov2026multiphysics,
  title={Multiphysics Failure Modes of \beta-Cell Encapsulation: From the Krogh Diffusion Limit to Hypoimmune Organoid Digital Twins},
  author={Naumov, Pavel V.},
  journal={bioRxiv},
  year={2026},
  doi={To be assigned},
  url={https://github.com/Bezooom/t1d-simulator}
}
```

---

## 📄 License
This project is licensed under the **MIT License** — see the `LICENSE` file for details. Open for academic and commercial co-development. For research collaborations and custom in silico design sweeps, please contact `naumov122@gmail.com`.
