# Release notes — `v1.0.0`

**Tag:** `v1.0.0`  
**Package:** `t1d_simulator` 1.0.0  
**Date:** 2026-07-26  
**License:** MIT  

Use this text as the GitHub **Release description** (paste below the horizontal rule into the release form). After Zenodo minting and bioRxiv posting, fill the placeholders in **Citation**.

---

## t1d_simulator v1.0.0 — Multiphysics Digital Twin for T1D Cell Therapy

Open-source, literature-calibrated multiphysics digital twin for screening **macroencapsulation failure modes** and **hypoimmune organoid** design choices *before* wet-lab experiments.

**Positioning:** independent multiphysics failure analysis + design-space screening; model co-development with academic and industry labs.  
**Not:** a clinical product, wet-lab protocol, or claim of a T1D cure.

### Highlights

- **Oxygen PDE engine** — planar / cylindrical / spherical reaction–diffusion with Michaelis–Menten kinetics; FBR, ROS, pH, OGM, cytokine coupling
- **PINN surrogate** — mesh-free steady-state solver for fast design sweeps
- **IBMIR 0–48 h** — TF → thrombin → clot → O₂ drop → viability (heparin / site effects)
- **Site decision matrix** — portal vs omental pouch vs subcutaneous
- **TPMS CAD** — Gyroid / Schwarz meshes + STL export for bioprinting workflows
- **GNN biocompatibility** — exploratory SMILES → \(L_{fib}\) ranking (N≈46; see limitations)
- **Design-space tools** — CLI `screen_design`, tornado uncertainty analysis, partner batch screening
- **Interactive UI** — Streamlit app for scenario exploration
- **Open science pack** — literature parameters, benchmark suite, preprint drafts (EN/RU), Zenodo metadata

### Validation (v1)

| Benchmark | Metric | Result |
|-----------|--------|--------|
| Papas et al. (spherical capsules) | viability RMSE | **10.42%** |
| Papabathini et al. (planar slabs) | viability RMSE | **11.42%** |
| Integration suite | `verify_model.py` | **40 / 40 green** |

Reproduce:

```bash
python3 reports/benchmarks/reproduce_benchmarks.py
CUDA_VISIBLE_DEVICES="" python3 t1d_simulator/verify_model.py
```

### Milestone coverage (in-repo)

| Milestone | Status |
|-----------|--------|
| M0 Engineering (`parameters.yaml`, `literature_params.yaml`, loaders, CI) | Done |
| M1 Calibration (Krogh/O₂, IBMIR, sites, lit benchmarks) | Done |
| M2 Open Science pack (preprint drafts, `screen_design`, uncertainty, `zenodo.json`) | Done |
| M3 Outreach Ready (UI screening, partner script, email drafts) | Done |

Public **bioRxiv** posting and **Zenodo DOI** are external steps after this tag (see `ROADMAP.md` v3.0).

### Install

```bash
git clone <REPO_URL>
cd <repo>
pip install -e t1d_simulator
# or
pip install -r t1d_simulator/requirements.txt
```

Requires Python ≥ 3.9. Optional: `pip install -e "t1d_simulator[dev]"` for pytest/ruff.

### Quickstart

```bash
# Verify (prefer CPU if GPU memory is tight)
CUDA_VISIBLE_DEVICES="" python3 t1d_simulator/verify_model.py

# Interactive UI
streamlit run t1d_simulator/app.py

# Design-space screen
python3 -m t1d_simulator.screen_design --help
# or, from package dir depending on install:
python3 t1d_simulator/screen_design.py --help

# Partner batch screening
python3 scripts/run_partner_screening.py

# Uncertainty (tornado)
python3 t1d_simulator/uncertainty_analysis.py
```

### Repository map

| Path | Role |
|------|------|
| `t1d_simulator/simulator.py` | Core O₂ multiphysics PDE |
| `t1d_simulator/organoid_simulator.py` / `ibmir_module.py` | Organoid + IBMIR |
| `t1d_simulator/pinn_solver.py` | PINN surrogate |
| `t1d_simulator/gnn_pipeline.py` | Polymer GNN (exploratory) |
| `t1d_simulator/screen_design.py` | Failure-mode screening CLI/API |
| `t1d_simulator/uncertainty_analysis.py` | Parameter sensitivity |
| `data/literature_params.yaml` | Literature-calibrated parameters |
| `reports/benchmarks/` | Predicted vs reported reports + figures |
| `docs/manuscript_biorxiv_en.md` | Preprint manuscript **v1.1** (EN, primary for bioRxiv) |
| `docs/outreach_*` / `reports/outreach_emails/` | Collab pack |
| `zenodo.json` | Citeable software metadata |
| `ROADMAP.md` | Project roadmap (M0–M5, external actions) |

### Limitations (please read)

- **In silico only** — not a substitute for peer review, animal studies, or clinical data
- **IBMIR ODE** — semi-mechanistic compartment model, not a full coagulation cascade
- **GNN** — small dataset (N≈46); ranking aid, not a validated materials predictor
- **Organoid CGM / TIR** — scenario tooling, not clinical trial outcomes
- **Gene-edit / iCasp9 narratives** — design checklists grounded in literature, not wet CRISPR SOPs

Full critical analysis: `Critical_Analysis_*.md` and manuscript Methods/Limitations.

### Citation

Until Zenodo and bioRxiv IDs are live, cite the software as:

```text
Naumov, P. V. (2026). t1d_simulator: Multiphysics Digital Twin Engine for
Encapsulated Beta-Cell & Hypoimmune Organoid Transplants (v1.0.0) [Computer software].
https://github.com/<ORG_OR_USER>/<REPO>
```

After deposit, prefer:

```text
Naumov, P. V. (2026). t1d_simulator (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.XXXX
```

Preprint (when posted): `docs/manuscript_biorxiv_en.md` → bioRxiv URL TBD.

### What's next

1. Link this release to Zenodo → record DOI in `ROADMAP.md` / citation block  
2. Submit preprint to bioRxiv (Synthetic Biology / Bioengineering)  
3. Outreach with live code + preprint URLs (`reports/outreach_emails/`)  
4. M4 joint science after partner signal; M5 wet-lab only with a lab partner  

### Acknowledgments & contact

Author: **Pavel V. Naumov** (Independent Bioengineering Researcher)  
Collaboration offers: free design-space screen on partner geometry / density / site — see `docs/outreach_1pager_en.md`.

---

**Checksum note for maintainers:** run `CUDA_VISIBLE_DEVICES="" python3 t1d_simulator/verify_model.py` on a clean clone before announcing the release.
