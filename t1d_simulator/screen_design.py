# -*- coding: utf-8 -*-
"""
Construct Failure Screening API & CLI (WP-Screen)
Calculates core pO2, Krogh limit safety, IBMIR 48h cell retention, and generates structured CSV/MD reports.
"""

import sys
import os
import argparse
import json
import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from simulator import solve_oxygen_profile, HYDROGELS, IMPLANTATION_SITES
from organoid_simulator import solve_ibmir_kinetics

def screen_design(
    geometry="spherical",
    radius_microns=250.0,
    density_million_per_ml=50.0,
    site_key="omental_pouch",
    heparin_coated=True,
    l_fibrosis_microns=0.0,
    hydrogel_key="alginate_2%"
):
    """
    Executes multiphysics failure screening for a given construct design.
    Returns a dict with viability, core pO2, IBMIR retention, status badge, and recommendations.
    """
    hydrogel = HYDROGELS.get(hydrogel_key, HYDROGELS["alginate_2%"])
    site = IMPLANTATION_SITES.get(site_key, IMPLANTATION_SITES["omental_pouch"])
    
    # 1. Oxygen Diffusion PDE
    sim_res = solve_oxygen_profile(
        R_outer_microns=radius_microns,
        rho_million_per_ml=density_million_per_ml,
        p_boundary=site["pO2"],
        D_coefficient=hydrogel["D"],
        geometry=geometry,
        L_fibrosis_microns=l_fibrosis_microns
    )
    
    # 2. IBMIR Coagulation ODE
    peg_density = 1.0 if heparin_coated else 0.0
    ibmir_res = solve_ibmir_kinetics(t_hours=[0, 12, 24, 48], peg_lmwh_density=peg_density)
    
    vf = float(sim_res["viable_fraction"])
    min_pO2 = float(sim_res["min_pO2"])
    ins_cap = float(sim_res["insulin_capacity"])
    retention_48h = float(ibmir_res["retention_48h_percent"])
    
    # Failure Mode Assessment
    failures = []
    if min_pO2 < 0.5:
        failures.append("Krogh Limit Core Anoxia (pO2 < 0.5 mmHg)")
    elif min_pO2 < 5.0:
        failures.append("Severe Core Hypoxia (pO2 < 5.0 mmHg)")
        
    if retention_48h < 85.0:
        failures.append("Acute IBMIR Thrombin Clot Lysis (Retention < 85%)")
        
    if vf < 70.0:
        failures.append("Overall Cell Viability Failure (< 70%)")
        
    if not failures:
        status_badge = "PASS ✅"
        recommendation = "Optimal construct topology. Suitable for preclinical trial design."
    elif vf >= 70.0 and min_pO2 >= 0.5:
        status_badge = "WARNING ⚠️"
        recommendation = "Sub-optimal parameters. Consider adding VEGF neovascularization or reducing cell density."
    else:
        status_badge = "FAIL 🚨"
        recommendation = "Critical failure predicted! Reduce outer radius / half-thickness < 150 µm or switch to Omental Pouch site."
        
    report = {
        "status_badge": status_badge,
        "geometry": geometry,
        "radius_microns": radius_microns,
        "density_million_per_ml": density_million_per_ml,
        "implantation_site": site["name"],
        "boundary_pO2_mmHg": site["pO2"],
        "hydrogel_material": hydrogel_key,
        "heparin_coated": heparin_coated,
        "fibrotic_thickness_microns": l_fibrosis_microns,
        "viable_fraction_percent": vf,
        "min_core_pO2_mmHg": min_pO2,
        "insulin_secretion_percent": ins_cap,
        "ibmir_48h_retention_percent": retention_48h,
        "detected_failure_modes": failures,
        "recommendation": recommendation
    }
    
    return report

def export_markdown_report(report_dict, output_path=None):
    heparin_str = "Yes (Lipid-PEG-LMWH)" if report_dict['heparin_coated'] else "No"
    vf_pass = "PASS" if report_dict['viable_fraction_percent'] >= 70.0 else "FAIL"
    po2_pass = "PASS" if report_dict['min_core_pO2_mmHg'] >= 0.5 else "FAIL"
    ins_pass = "PASS" if report_dict['insulin_secretion_percent'] >= 80.0 else "WARN"
    ibmir_pass = "PASS" if report_dict['ibmir_48h_retention_percent'] >= 90.0 else "FAIL"
    
    md = f"""# Construct Failure Screening Report (`screen_design`)

**Status:** {report_dict['status_badge']}  
**Timestamp:** July 2026  

---

## Design Parameters
- **Topology:** `{report_dict['geometry'].upper()}`
- **Radius / Half-Thickness:** `{report_dict['radius_microns']} µm`
- **Cell Loading Density:** `{report_dict['density_million_per_ml']} M cells/mL`
- **Hydrogel Matrix:** `{report_dict['hydrogel_material']}`
- **Implantation Site:** `{report_dict['implantation_site']}` (pO2 = {report_dict['boundary_pO2_mmHg']} mmHg)
- **Heparinized Surface:** `{heparin_str}`
- **Fibrotic Capsule (L_fib):** `{report_dict['fibrotic_thickness_microns']} µm`

---

## Multiphysics Performance Metrics
| Metric | Value | Reference Threshold | Assessment |
|---|---|---|---|
| **Cell Viability (f_viab)** | **{report_dict['viable_fraction_percent']:.1f}%** | > 70.0% | {vf_pass} |
| **Minimum Core pO2** | **{report_dict['min_core_pO2_mmHg']:.2f} mmHg** | > 0.5 mmHg | {po2_pass} |
| **Insulin Capacity** | **{report_dict['insulin_secretion_percent']:.1f}%** | > 80.0% | {ins_pass} |
| **48h IBMIR Cell Retention** | **{report_dict['ibmir_48h_retention_percent']:.1f}%** | > 90.0% | {ibmir_pass} |

---

## Detected Failure Modes
"""
    if report_dict['detected_failure_modes']:
        for f in report_dict['detected_failure_modes']:
            md += f"- 🚨 **{f}**\n"
    else:
        md += "- ✅ **No physical or biological failure modes detected.**\n"
        
    md += f"\n**Recommendation:** {report_dict['recommendation']}\n"
    
    if output_path:
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(md)
            
    return md

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct Failure Screening Utility (WP-Screen)")
    parser.add_argument("--geometry", type=str, default="spherical", choices=["planar", "cylindrical", "spherical"])
    parser.add_argument("--radius", type=float, default=250.0, help="Outer radius or half-thickness in microns")
    parser.add_argument("--density", type=float, default=50.0, help="Cell seeding density in M cells/ml")
    parser.add_argument("--site", type=str, default="omental_pouch", help="Target implantation site key")
    parser.add_argument("--no-heparin", action="store_true", help="Disable surface heparinization")
    parser.add_argument("--fibrosis", type=float, default=0.0, help="Fibrotic capsule thickness in microns")
    parser.add_argument("--export-md", type=str, default=None, help="Output markdown report filepath")
    
    args = parser.parse_args()
    
    res = screen_design(
        geometry=args.geometry,
        radius_microns=args.radius,
        density_million_per_ml=args.density,
        site_key=args.site,
        heparin_coated=not args.no_heparin,
        l_fibrosis_microns=args.fibrosis
    )
    
    md_output = export_markdown_report(res, args.export_md)
    print(md_output)
