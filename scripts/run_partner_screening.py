# -*- coding: utf-8 -*-
"""
Partner Construct Batch Screening Utility
Evaluates 4 target partner constructs in batch mode using screen_design()
and generates a comparative evaluation markdown report.
"""

import os
import sys
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "t1d_simulator"))

from screen_design import screen_design, export_markdown_report

PARTNER_CONSTRUCTS = [
    {
        "partner": "Harvard SEAS (Mooney Lab)",
        "geometry": "spherical",
        "radius": 400.0,
        "density": 50.0,
        "site": "subcutaneous",
        "heparin": False,
        "fibrosis": 25.0,
        "notes": "Standard macroporous alginate in subcutaneous space (uncoated)."
    },
    {
        "partner": "MIT Koch Institute (Anderson Lab)",
        "geometry": "spherical",
        "radius": 250.0,
        "density": 60.0,
        "site": "subcutaneous",
        "heparin": False,
        "fibrosis": 5.0,
        "notes": "Zwitterionic SBMA modified alginate microcapsules (low FBR)."
    },
    {
        "partner": "Sana Biotechnology",
        "geometry": "planar",
        "radius": 150.0,
        "density": 80.0,
        "site": "omental_pouch",
        "heparin": True,
        "fibrosis": 0.0,
        "notes": "Hypoimmune iPSC organoid sheet with heparin surface in omentum."
    },
    {
        "partner": "Seraxis / Sernova",
        "geometry": "planar",
        "radius": 250.0,
        "density": 100.0,
        "site": "omental_pouch",
        "heparin": True,
        "fibrosis": 0.0,
        "notes": "Thick pre-vascularized pouch macro-device with heparin."
    }
]

def run_batch_partner_screening():
    output_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reports")
    os.makedirs(output_dir, exist_ok=True)
    report_path = os.path.join(output_dir, "partner_screening_comparison.md")
    
    results = []
    print("=== Running Batch Partner Construct Multiphysics Screening ===")
    
    for c in PARTNER_CONSTRUCTS:
        res = screen_design(
            geometry=c["geometry"],
            radius_microns=c["radius"],
            density_million_per_ml=c["density"],
            site_key=c["site"],
            heparin_coated=c["heparin"],
            l_fibrosis_microns=c["fibrosis"]
        )
        res["partner"] = c["partner"]
        res["notes"] = c["notes"]
        results.append(res)
        print(f"  [{res['status_badge']}] {c['partner']} -> Viability: {res['viable_fraction_percent']:.1f}%, Core pO2: {res['min_core_pO2_mmHg']:.2f} mmHg, IBMIR: {res['ibmir_48h_retention_percent']:.1f}%")
        
    # Generate Markdown Summary Table
    md = "# Comparative Partner Construct Failure Screening Report\n\n"
    md += "**Generated:** July 2026  \n"
    md += "**Engine:** `t1d_simulator` Multiphysics Screening API (`screen_design`)  \n\n"
    md += "---  \n\n"
    md += "## Executive Summary Table\n\n"
    md += "| Partner / Group | Construct Topology | Site | Viability | Core pO₂ | 48h IBMIR Retention | Status | Recommendation |\n"
    md += "|---|---|---|---|---|---|---|---|\n"
    
    for r in results:
        top_str = f"{r['geometry'].upper()} ({r['radius_microns']} µm, {r['density_million_per_ml']} M/ml)"
        md += f"| **{r['partner']}** | `{top_str}` | {r['implantation_site']} | **{r['viable_fraction_percent']:.1f}%** | **{r['min_core_pO2_mmHg']:.2f} mmHg** | **{r['ibmir_48h_retention_percent']:.1f}%** | **{r['status_badge']}** | {r['recommendation']} |\n"
        
    md += "\n---\n\n"
    md += "## Detailed Construct Diagnostics\n\n"
    for r in results:
        md += f"### 🔬 {r['partner']}\n"
        md += f"- **Notes:** {r['notes']}\n"
        md += f"- **Detected Failures:**\n"
        if r['detected_failure_modes']:
            for f in r['detected_failure_modes']:
                md += f"  - 🚨 {f}\n"
        else:
            md += "  - ✅ No multiphysics failure modes detected.\n"
        md += f"- **Actionable Insight:** {r['recommendation']}\n\n"
        
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(md)
        
    print(f"\n[SUCCESS] Comparative partner screening report saved to: {report_path}")
    return report_path

if __name__ == "__main__":
    run_batch_partner_screening()
