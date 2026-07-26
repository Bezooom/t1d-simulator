# Preprint Abstract Draft (bioRxiv Submission)

**Title:** Multiphysics Failure Modes of $\beta$-Cell Encapsulation: From the Krogh Diffusion Limit to Hypoimmune Organoid Digital Twins  
**Target Category:** Bioengineering / Synthetic Biology / Systems Biology  
**Date:** 2026-07-24  
**File:** [preprint_abstract_en.md](file:///home/bezoom/storage/Projects/%D0%94%D0%B8%D0%B0%D0%B1%D0%B5%D1%82/docs/preprint_abstract_en.md)

---

## Abstract (285 words)

Cell replacement therapy for Type 1 Diabetes (T1D) using microencapsulated stem cell-derived $\beta$-cells holds immense curative potential, yet clinical translation remains hindered by unpredictable graft failure. Here, we present an open-source, calibrated multiphysics digital twin engine (`t1d_simulator`) designed to quantify physical and biochemical failure mechanisms across macroencapsulation topologies and hypoimmune organoid architectures.

Integrating non-linear 1D/3D reaction-diffusion partial differential equations (PDEs), Physics-Informed Neural Networks (PINNs), and Graph Neural Networks (GNNs) for anti-fibrotic surface chemistry, our computational framework simulates oxygen transport, Foreign Body Reaction (FBR) capsule growth, reactive oxygen species (ROS) accumulation, and 0–48h Instant Blood-Mediated Inflammatory Reaction (IBMIR) kinetics. 

We benchmarked model predictions against reported literature viability data (Papas et al. 2007; Papabathini et al. 2023), demonstrating robust predictive accuracy ($RMSE < 12\%$). Quantitative failure analysis confirms that standard macroencapsulation is systematically constrained by the Krogh diffusion limit: planar slabs ($L > 200\ \mu\text{m}$) and dense microcapsules ($>80\ \text{M/ml}$) develop severe core hypoxia ($pO_2 < 0.5\ \text{mmHg}$) and necrotic cores under physiological subcutaneous conditions ($pO_2 = 30\ \text{mmHg}$). 

To overcome these constraints, we evaluated a multi-layered design paradigm combining multiplex gene editing ($B2M^{-/-}$, $CIITA^{-/-}$, $CD47^{KI}$, $PD-L1^{KI}$), inducible $iCasp9$ suicide switches, Lipid-PEG-LMWH heparinization, and pre-vascularized omental pouch implantation. Model simulations show that pre-vascularized omental scaffolds coupled with VEGF/PDGF angiogenic feedback rescue core oxygenation ($pO_2 > 40\ \text{mmHg}$), suppress IBMIR thrombin generation by $>80\%$, and sustain long-term glucose-stimulated insulin secretion (TIR $>98\%$). 

This digital twin provides a calibrated failure-analysis platform to screen geometry, cell loading, graft site, and edit-set combinations prior to costly wet-lab experimentation, accelerating the rational engineering of durable T1D cell therapies.
