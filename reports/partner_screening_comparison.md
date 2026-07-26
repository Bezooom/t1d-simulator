# Comparative Partner Construct Failure Screening Report

**Generated:** July 2026  
**Engine:** `t1d_simulator` Multiphysics Screening API (`screen_design`)  

---  

## Executive Summary Table

| Partner / Group | Construct Topology | Site | Viability | Core pO₂ | 48h IBMIR Retention | Status | Recommendation |
|---|---|---|---|---|---|---|---|
| **Harvard SEAS (Mooney Lab)** | `SPHERICAL (400.0 µm, 50.0 M/ml)` | Подкожная клетчатка (SQ) | **70.6%** | **0.00 mmHg** | **81.1%** | **FAIL 🚨** | Critical failure predicted! Reduce outer radius / half-thickness < 150 µm or switch to Omental Pouch site. |
| **MIT Koch Institute (Anderson Lab)** | `SPHERICAL (250.0 µm, 60.0 M/ml)` | Подкожная клетчатка (SQ) | **97.7%** | **0.02 mmHg** | **81.1%** | **FAIL 🚨** | Critical failure predicted! Reduce outer radius / half-thickness < 150 µm or switch to Omental Pouch site. |
| **Sana Biotechnology** | `PLANAR (150.0 µm, 80.0 M/ml)` | Сальник (Omental Pouch, реваскуляризованный) | **100.0%** | **4.36 mmHg** | **94.6%** | **WARNING ⚠️** | Sub-optimal parameters. Consider adding VEGF neovascularization or reducing cell density. |
| **Seraxis / Sernova** | `PLANAR (250.0 µm, 100.0 M/ml)` | Сальник (Omental Pouch, реваскуляризованный) | **55.0%** | **0.00 mmHg** | **94.6%** | **FAIL 🚨** | Critical failure predicted! Reduce outer radius / half-thickness < 150 µm or switch to Omental Pouch site. |

---

## Detailed Construct Diagnostics

### 🔬 Harvard SEAS (Mooney Lab)
- **Notes:** Standard macroporous alginate in subcutaneous space (uncoated).
- **Detected Failures:**
  - 🚨 Krogh Limit Core Anoxia (pO2 < 0.5 mmHg)
  - 🚨 Acute IBMIR Thrombin Clot Lysis (Retention < 85%)
- **Actionable Insight:** Critical failure predicted! Reduce outer radius / half-thickness < 150 µm or switch to Omental Pouch site.

### 🔬 MIT Koch Institute (Anderson Lab)
- **Notes:** Zwitterionic SBMA modified alginate microcapsules (low FBR).
- **Detected Failures:**
  - 🚨 Krogh Limit Core Anoxia (pO2 < 0.5 mmHg)
  - 🚨 Acute IBMIR Thrombin Clot Lysis (Retention < 85%)
- **Actionable Insight:** Critical failure predicted! Reduce outer radius / half-thickness < 150 µm or switch to Omental Pouch site.

### 🔬 Sana Biotechnology
- **Notes:** Hypoimmune iPSC organoid sheet with heparin surface in omentum.
- **Detected Failures:**
  - 🚨 Severe Core Hypoxia (pO2 < 5.0 mmHg)
- **Actionable Insight:** Sub-optimal parameters. Consider adding VEGF neovascularization or reducing cell density.

### 🔬 Seraxis / Sernova
- **Notes:** Thick pre-vascularized pouch macro-device with heparin.
- **Detected Failures:**
  - 🚨 Krogh Limit Core Anoxia (pO2 < 0.5 mmHg)
  - 🚨 Overall Cell Viability Failure (< 70%)
- **Actionable Insight:** Critical failure predicted! Reduce outer radius / half-thickness < 150 µm or switch to Omental Pouch site.

