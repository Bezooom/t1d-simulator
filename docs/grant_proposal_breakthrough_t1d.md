# Research Grant Proposal Blueprint: Breakthrough T1D / NIH R01 Application

**Project Title:** Calibrated Multiphysics Digital Twin Platform for Rational Engineering and Failure-Mode Screening of Hypoimmune $\beta$-Cell Therapies  
**Target Funding Agency:** Breakthrough T1D (formerly JDRF) Strategic Research Grant / NIH NIDDK (R01)  
**Proposed Project Period:** 3 Years (2027 – 2030)  
**Total Direct Requested Budget:** $1,250,000 USD  
**Date:** July 2026  
**Document File:** [grant_proposal_breakthrough_t1d.md](file:///home/bezoom/storage/Projects/%D0%94%D0%B8%D0%B0%D0%B1%D0%B5%D1%82/docs/grant_proposal_breakthrough_t1d.md)

---

## 1. Project Summary / Abstract

Stem cell-derived insulin-producing $\beta$-cell replacement therapies represent a transformational curative approach for Type 1 Diabetes (T1D). However, clinical translation is constrained by high attrition rates during preclinical testing caused by unresolved physical and biological bottlenecks: severe core hypoxia exceeding the Krogh limit, foreign body reaction (FBR) encapsulation, and acute Instant Blood-Mediated Inflammatory Reaction (IBMIR). 

This proposal establishes a multi-institutional computational-experimental framework to de-risk T1D cell therapy development. Building on our open-source, literature-calibrated digital twin engine (`t1d_simulator`), we propose to:
1. Advance 3D multiphysics transport models coupling oxygen diffusion, ROS toxicity, cytokine kinetics, and neovascularization.
2. Calibrate machine learning pipelines (GNN surface chemistry screening) using high-throughput microfluidic hydrogel data.
3. Perform closed-loop *in vitro* and *in vivo* preclinical validation in partnership with leading cell therapy laboratories.

---

## 2. Specific Aims

### Specific Aim 1: Multiphysics Model Expansion & High-Resolution 3D Transport Dynamics (Months 1–12)
- Extend 1D reaction-diffusion solvers into full 3D boundary-element CFD models incorporating patient-specific microvascular geometries.
- Integrate multi-scale models of hypoimmune gene editing ($B2M^{-/-}$, $CIITA^{-/-}$, $CD47^{KI}$, $PD-L1^{KI}$) and $iCasp9$ suicide switch elimination kinetics under inflammatory cytokines (IL-1$\beta$, TNF-$\alpha$, IFN-$\gamma$).

### Specific Aim 2: Machine Learning Guided Bio-Interface & Anti-Fibrotic Screening (Months 13–24)
- Expand the Graph Neural Network (GNN) biocompatibility predictor from 52 to $>500$ zwitterionic and PEGylated hydrogel formulations.
- Predict fibrotic capsule deposition ($L_{fib}$) and mass transfer degradation as a function of polymer SMILES structures.

### Specific Aim 3: Preclinical In Vitro and In Vivo Validation in Omental Pouch Models (Months 25–36)
- Perform 3D bioprinting of TPMS Gyroid scaffolds ($200\ \mu\text{m}$ pores) seeded with human iPSC-derived islets.
- Validate predicted core oxygenation ($pO_2 > 40\ \text{mmHg}$), anti-IBMIR retention ($>90\%$), and 30-day CGM Time in Range ($TIR > 95\%$) in mouse and porcine omental pouch models.

---

## 3. Significance & Innovation

- **Pioneering Digital Twin in T1D Cell Therapy**: Provides the first open-source, literature-calibrated multiphysics simulation engine dedicated to failure-mode screening.
- **De-Risking Preclinical Pipelines**: Reduces expensive trial-and-error animal studies by predicting hypoxic thresholds and IBMIR clot kinetics prior to implantation.
- **Open Science & Community Impact**: Delivers citeable parameter databases (`literature_params.yaml`), automated screening tools, and web dashboards for global academic and industry researchers.

---

## 4. Budget & Work Package Overview

| Year | Work Package Focus | Key Deliverables | Direct Budget |
|---|---|---|---|
| **Year 1** | 3D Multiphysics PDE & Gene Edit Kinetics | 3D CFD solver, expanded parameter database | $400,000 |
| **Year 2** | GNN Polymer Screening & High-Throughput Microfluidics | Trained GNN model, microfluidic hydrogel assays | $425,000 |
| **Year 3** | Preclinical *In Vivo* Validation (Omental Pouch) | Mouse/Porcine efficacy data, joint publications | $425,000 |

---

## 5. Participating Institutions & Co-Investigators

- **Lead Computational Group**: In Silico Bioengineering Consortium (`t1d_simulator` maintainers)
- **Experimental Wet-Lab Co-PI**: Biomaterials & Cell Therapy Laboratory
- **Clinical Endocrine Advisor**: Academic Diabetes Center
