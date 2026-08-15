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
        "status": status_badge.split()[0],
        "geometry": geometry,
        "radius_microns": radius_microns,
        "density_million_per_ml": density_million_per_ml,
        "site_key": site_key if site_key in IMPLANTATION_SITES else "omental_pouch",
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


DEFAULT_DENSITIES = (20.0, 50.0, 80.0, 100.0)
DEFAULT_RADII = (100.0, 150.0, 200.0, 250.0)
SWEEP_CSV_COLUMNS = [
    "geometry",
    "site_key",
    "radius_microns",
    "density_million_per_ml",
    "heparin_coated",
    "fibrotic_thickness_microns",
    "boundary_pO2_mmHg",
    "viable_fraction_percent",
    "min_core_pO2_mmHg",
    "insulin_secretion_percent",
    "ibmir_48h_retention_percent",
    "status",
    "detected_failure_modes",
]


def _parse_float_list(raw, fallback):
    if raw is None or str(raw).strip() == "":
        return list(fallback)
    values = [float(part.strip()) for part in str(raw).split(",") if part.strip()]
    if not values:
        raise ValueError("Expected a comma-separated list of numbers")
    return values


def sweep_designs(
    geometry="planar",
    densities=None,
    radii=None,
    site_key="subcutaneous",
    heparin_coated=False,
    l_fibrosis_microns=0.0,
    hydrogel_key="alginate_2%",
):
    """Cartesian sweep of density × radius for one geometry/site."""
    densities = list(densities) if densities is not None else list(DEFAULT_DENSITIES)
    radii = list(radii) if radii is not None else list(DEFAULT_RADII)
    rows = []
    for density in densities:
        for radius in radii:
            report = screen_design(
                geometry=geometry,
                radius_microns=radius,
                density_million_per_ml=density,
                site_key=site_key,
                heparin_coated=heparin_coated,
                l_fibrosis_microns=l_fibrosis_microns,
                hydrogel_key=hydrogel_key,
            )
            rows.append({
                "geometry": report["geometry"],
                "site_key": report["site_key"],
                "radius_microns": report["radius_microns"],
                "density_million_per_ml": report["density_million_per_ml"],
                "heparin_coated": report["heparin_coated"],
                "fibrotic_thickness_microns": report["fibrotic_thickness_microns"],
                "boundary_pO2_mmHg": report["boundary_pO2_mmHg"],
                "viable_fraction_percent": report["viable_fraction_percent"],
                "min_core_pO2_mmHg": report["min_core_pO2_mmHg"],
                "insulin_secretion_percent": report["insulin_secretion_percent"],
                "ibmir_48h_retention_percent": report["ibmir_48h_retention_percent"],
                "status": report["status"],
                "detected_failure_modes": "; ".join(report["detected_failure_modes"]),
            })
    return pd.DataFrame(rows, columns=SWEEP_CSV_COLUMNS)


def export_sweep_csv(df, output_path):
    parent = os.path.dirname(os.path.abspath(output_path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    df.to_csv(output_path, index=False)
    return output_path


def export_sweep_json(df, output_path):
    parent = os.path.dirname(os.path.abspath(output_path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    payload = {
        "n_designs": int(len(df)),
        "n_fail": int((df["status"] == "FAIL").sum()),
        "n_warning": int((df["status"] == "WARNING").sum()),
        "n_pass": int((df["status"] == "PASS").sum()),
        "designs": df.to_dict(orient="records"),
    }
    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    return output_path


def plot_sweep_heatmap(df, output_path, value_col="viable_fraction_percent"):
    """Save a density × radius heatmap PNG. Requires matplotlib."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    if df.empty:
        raise ValueError("Sweep dataframe is empty")
    if value_col not in df.columns:
        raise KeyError(value_col)

    pivot = df.pivot_table(
        index="radius_microns",
        columns="density_million_per_ml",
        values=value_col,
        aggfunc="mean",
    ).sort_index(ascending=False)
    parent = os.path.dirname(os.path.abspath(output_path))
    if parent:
        os.makedirs(parent, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    data = pivot.to_numpy(dtype=float)
    if value_col == "min_core_pO2_mmHg":
        vmax = max(5.0, float(np.nanmax(data)))
        cmap = "magma"
        label = "min core pO2 (mmHg)"
        fmt = "{:.1f}"
    else:
        vmax = 100.0
        cmap = "RdYlGn"
        label = "viable fraction (%)"
        fmt = "{:.0f}"

    image = ax.imshow(data, aspect="auto", cmap=cmap, vmin=0.0, vmax=vmax)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([f"{c:g}" for c in pivot.columns])
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels([f"{r:g}" for r in pivot.index])
    ax.set_xlabel("Cell density (M/mL)")
    ax.set_ylabel("Radius / half-thickness (µm)")
    site = df["site_key"].iloc[0]
    geom = df["geometry"].iloc[0]
    ax.set_title(f"{geom} · {site} · {label}")
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label=label)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j, i, fmt.format(data[i, j]), ha="center", va="center", fontsize=8, color="black")
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct Failure Screening Utility (WP-Screen)")
    parser.add_argument("--geometry", type=str, default="spherical", choices=["planar", "cylindrical", "spherical"])
    parser.add_argument("--radius", type=float, default=250.0, help="Outer radius or half-thickness in microns")
    parser.add_argument("--density", type=float, default=50.0, help="Cell seeding density in M cells/ml")
    parser.add_argument("--site", type=str, default="omental_pouch", help="Target implantation site key")
    parser.add_argument("--no-heparin", action="store_true", help="Disable surface heparinization")
    parser.add_argument("--fibrosis", type=float, default=0.0, help="Fibrotic capsule thickness in microns")
    parser.add_argument("--export-md", type=str, default=None, help="Output markdown report filepath")
    parser.add_argument("--sweep", action="store_true", help="Run density × radius sweep instead of a single design")
    parser.add_argument("--densities", type=str, default=None, help="Comma-separated densities (M/mL)")
    parser.add_argument("--radii", type=str, default=None, help="Comma-separated radii / half-thicknesses (µm)")
    parser.add_argument("--export-csv", type=str, default=None, help="Sweep CSV output path")
    parser.add_argument("--export-json", type=str, default=None, help="Sweep JSON output path")
    parser.add_argument("--heatmap", type=str, default=None, help="Viability heatmap PNG path")
    parser.add_argument("--heatmap-po2", type=str, default=None, help="Core pO2 heatmap PNG path")

    args = parser.parse_args()

    if args.sweep:
        df = sweep_designs(
            geometry=args.geometry,
            densities=_parse_float_list(args.densities, DEFAULT_DENSITIES),
            radii=_parse_float_list(args.radii, DEFAULT_RADII),
            site_key=args.site,
            heparin_coated=not args.no_heparin,
            l_fibrosis_microns=args.fibrosis,
        )
        print(df.to_string(index=False))
        print(
            f"\nSummary: {len(df)} designs | "
            f"PASS={int((df['status']=='PASS').sum())} "
            f"WARNING={int((df['status']=='WARNING').sum())} "
            f"FAIL={int((df['status']=='FAIL').sum())}"
        )
        if args.export_csv:
            export_sweep_csv(df, args.export_csv)
            print(f"CSV: {args.export_csv}")
        if args.export_json:
            export_sweep_json(df, args.export_json)
            print(f"JSON: {args.export_json}")
        if args.heatmap:
            plot_sweep_heatmap(df, args.heatmap, value_col="viable_fraction_percent")
            print(f"Heatmap: {args.heatmap}")
        if args.heatmap_po2:
            plot_sweep_heatmap(df, args.heatmap_po2, value_col="min_core_pO2_mmHg")
            print(f"Heatmap pO2: {args.heatmap_po2}")
    else:
        res = screen_design(
            geometry=args.geometry,
            radius_microns=args.radius,
            density_million_per_ml=args.density,
            site_key=args.site,
            heparin_coated=not args.no_heparin,
            l_fibrosis_microns=args.fibrosis,
        )
        print(export_markdown_report(res, args.export_md))
