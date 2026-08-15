# -*- coding: utf-8 -*-
"""
Dose–transport conflict for T1D cell replacement.

Clinical mass (IEQ) wants volume. Oxygen L* forbids thickness.
Therefore the implant must grow in *area*, not in height.
"""

from __future__ import annotations

import json
import os
import sys

import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from design_bounds import V_MAX_HUMAN, analytical_l_crit_um
from simulator import HYDROGELS, IMPLANTATION_SITES

CELLS_PER_IEQ = 1560.0
IEQ_PER_KG = 11000.0
P_SURVIVAL = 0.5
P_FUNCTION = 5.0  # K_M for GSIS in this codebase
SAFETY = 0.90


def therapeutic_cell_load(weight_kg=70.0, ieq_per_kg=IEQ_PER_KG):
    ieq = weight_kg * ieq_per_kg
    cells = ieq * CELLS_PER_IEQ
    return {
        "weight_kg": weight_kg,
        "ieq_per_kg": ieq_per_kg,
        "total_ieq": ieq,
        "total_cells": cells,
        "total_cells_millions": cells / 1e6,
    }


def l_star_pair(density, site_key, geometry, hydrogel_key="alginate_2%", v_max=V_MAX_HUMAN):
    p_b = IMPLANTATION_SITES[site_key]["pO2"]
    d = HYDROGELS[hydrogel_key]["D"]
    l_surv = analytical_l_crit_um(density, p_b, d, v_max, geometry, p_crit=P_SURVIVAL)
    l_func = analytical_l_crit_um(density, p_b, d, v_max, geometry, p_crit=P_FUNCTION)
    return p_b, l_surv, l_func


def sheet_area_cm2(volume_ml, half_thickness_um, two_sided=True):
    """Planar sheet area. Two-sided oxygenation ⇒ full thickness = 2L."""
    full_um = 2.0 * half_thickness_um if two_sided else half_thickness_um
    full_cm = full_um * 1e-4
    if full_cm <= 0:
        return float("inf")
    return volume_ml / full_cm


def dose_transport_table(
    weight_kg=70.0,
    densities=(20.0, 50.0, 80.0, 100.0),
    sites=("subcutaneous", "venous", "omental_pouch"),
    geometries=("planar", "spherical"),
):
    load = therapeutic_cell_load(weight_kg)
    rows = []
    for site in sites:
        for geom in geometries:
            for dens in densities:
                p_b, l_surv, l_func = l_star_pair(dens, site, geom)
                l_use = SAFETY * l_func  # functional certificate
                vol_ml = load["total_cells_millions"] / dens
                area_2s = sheet_area_cm2(vol_ml, l_use, two_sided=True) if geom == "planar" else None
                area_1s = sheet_area_cm2(vol_ml, l_use, two_sided=False) if geom == "planar" else None
                n_spheres = None
                if geom == "spherical":
                    r_cm = l_use * 1e-4
                    v_one_ml = (4.0 / 3.0) * 3.141592653589793 * r_cm**3 * 1000.0
                    n_spheres = int(vol_ml / v_one_ml) if v_one_ml > 0 else None
                rows.append({
                    "weight_kg": weight_kg,
                    "site_key": site,
                    "p_boundary_mmHg": p_b,
                    "geometry": geom,
                    "density_M_per_ml": dens,
                    "L_star_survival_um": round(l_surv, 1),
                    "L_star_function_um": round(l_func, 1),
                    "certified_L_func_um": round(l_use, 1),
                    "graft_volume_ml": round(vol_ml, 2),
                    "planar_area_two_sided_cm2": None if area_2s is None else round(area_2s, 1),
                    "planar_area_one_sided_cm2": None if area_1s is None else round(area_1s, 1),
                    "n_spheres_at_Rstar": n_spheres,
                })
    return load, pd.DataFrame(rows)


def score_literature_devices():
    """Public-typical construct sizes vs L* (functional, human V_max)."""
    specs = [
        {
            "name": "Alginate microcapsule (Papas-like)",
            "geometry": "spherical",
            "site_key": "subcutaneous",
            "density_M_per_ml": 50.0,
            "length_um": 400.0,
            "heparin": False,
        },
        {
            "name": "Planar SQ macrodevice, 500 µm full (2-sided → L=250)",
            "geometry": "planar",
            "site_key": "subcutaneous",
            "density_M_per_ml": 50.0,
            "length_um": 250.0,
            "heparin": False,
        },
        {
            "name": "Planar SQ thin, 200 µm full (L=100)",
            "geometry": "planar",
            "site_key": "subcutaneous",
            "density_M_per_ml": 50.0,
            "length_um": 100.0,
            "heparin": False,
        },
        {
            "name": "Omental sphere (this study, recommended)",
            "geometry": "spherical",
            "site_key": "omental_pouch",
            "density_M_per_ml": 80.0,
            "length_um": 200.0,
            "heparin": True,
        },
        {
            "name": "Portal naked / uncoated sphere",
            "geometry": "spherical",
            "site_key": "venous",
            "density_M_per_ml": 80.0,
            "length_um": 150.0,
            "heparin": False,
        },
    ]
    rows = []
    for spec in specs:
        p_b, l_surv, l_func = l_star_pair(spec["density_M_per_ml"], spec["site_key"], spec["geometry"])
        o2_surv = spec["length_um"] <= SAFETY * l_surv
        o2_func = spec["length_um"] <= SAFETY * l_func
        ibmir_ok = spec["heparin"] or spec["site_key"] in ("subcutaneous", "omental_pouch")
        if spec["site_key"] == "omental_pouch" and not spec["heparin"]:
            ibmir_ok = True  # low blood contact; still not zero
        if spec["site_key"] == "venous" and not spec["heparin"]:
            ibmir_ok = False
        if o2_func and ibmir_ok:
            verdict = "PASS"
        elif o2_surv and not o2_func:
            verdict = "SURVIVES / POOR GSIS"
        else:
            verdict = "FAIL"
        rows.append({
            "device": spec["name"],
            "L_or_R_um": spec["length_um"],
            "L_star_surv_um": round(l_surv, 1),
            "L_star_func_um": round(l_func, 1),
            "oxygen_survival": o2_surv,
            "oxygen_function": o2_func,
            "ibmir_ok": ibmir_ok,
            "verdict": verdict,
        })
    return pd.DataFrame(rows)


def run_study(output_dir=None, weight_kg=70.0):
    if output_dir is None:
        output_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "reports",
            "dose_transport",
        )
    os.makedirs(output_dir, exist_ok=True)
    load, table = dose_transport_table(weight_kg=weight_kg)
    devices = score_literature_devices()
    table.to_csv(os.path.join(output_dir, "dose_vs_area.csv"), index=False)
    devices.to_csv(os.path.join(output_dir, "literature_device_scores.csv"), index=False)
    path = _write_report(load, table, devices, output_dir)
    return {"load": load, "report": path, "n_rows": int(len(table))}


def _write_report(load, table, devices, output_dir):
    path = os.path.join(output_dir, "dose_vs_krogh.md")

    def md(df):
        cols = list(df.columns)
        lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
        for _, row in df.iterrows():
            lines.append("| " + " | ".join("" if pd.isna(row[c]) else str(row[c]) for c in cols) + " |")
        return "\n".join(lines)

    om_plan = table[(table.site_key == "omental_pouch") & (table.geometry == "planar") & (table.density_M_per_ml == 80.0)].iloc[0]
    sq_plan = table[(table.site_key == "subcutaneous") & (table.geometry == "planar") & (table.density_M_per_ml == 50.0)].iloc[0]
    om_sph = table[(table.site_key == "omental_pouch") & (table.geometry == "spherical") & (table.density_M_per_ml == 80.0)].iloc[0]

    highlight = table[
        ((table.geometry == "planar") & (table.density_M_per_ml.isin([50.0, 80.0])))
        | ((table.geometry == "spherical") & (table.density_M_per_ml == 80.0))
    ][
        [
            "site_key",
            "geometry",
            "density_M_per_ml",
            "L_star_survival_um",
            "L_star_function_um",
            "certified_L_func_um",
            "graft_volume_ml",
            "planar_area_two_sided_cm2",
            "n_spheres_at_Rstar",
        ]
    ]

    text = f"""# Dose vs oxygen: why encapsulation trials run out of area

**Date:** 2026-08-15  
**Patient:** {load['weight_kg']:.0f} kg · Edmonton-like **{load['ieq_per_kg']:.0f} IEQ/kg**  
**Mass:** {load['total_ieq']:.0f} IEQ = **{load['total_cells_millions']:.0f} million cells** ({CELLS_PER_IEQ:.0f} cells/IEQ)  
**Oxygen:** human $V_{{\\max}}=1.5\\times 10^{{-16}}$ mol/(cell·s), 2% alginate, certificate $0.9\\,L^*$  
**Two $p_{{\\mathrm{{crit}}}}$:** survival $0.5$ mmHg vs **function / GSIS $5$ mmHg** ($K_{{M,\\mathrm{{ins}}}}$)

## The biotech conflict

To replace insulin you need a **fixed cell mass**.  
Krogh forbids packing that mass into a thick lump.  
So the only remaining degree of freedom is **surgical area** (omental sheet, pouch, or many microcapsules).

$$
V = \\frac{{N_{{\\mathrm{{cells}}}}}}{{\\rho}},\\qquad
A_{{\\mathrm{{planar,\\ 2-sided}}}} = \\frac{{V}}{{2\\,L^{{\\ast}}_{{\\mathrm{{func}}}}}}
$$

$L^{{\\ast}}_{{\\mathrm{{func}}}}$ is the same zero-order bound as before, but with $p_c=5$ mmHg. It is **strictly smaller** than the survival bound ($p_c=0.5$). A graft can be “alive” on paper and still not secrete.

## Numbers that matter

| Scenario | $L^*_{{\\mathrm{{surv}}}}$ | $L^*_{{\\mathrm{{func}}}}$ | Graft volume | Surgical implication |
|----------|---------------------------|----------------------------|--------------|----------------------|
| Omentum planar 80 M/ml | {om_plan['L_star_survival_um']} µm | {om_plan['L_star_function_um']} µm | {om_plan['graft_volume_ml']} ml | **{om_plan['planar_area_two_sided_cm2']} cm²** two-sided sheet |
| SQ planar 50 M/ml | {sq_plan['L_star_survival_um']} µm | {sq_plan['L_star_function_um']} µm | {sq_plan['graft_volume_ml']} ml | **{sq_plan['planar_area_two_sided_cm2']} cm²** two-sided (worse if one-sided: {sq_plan['planar_area_one_sided_cm2']} cm²) |
| Omentum sphere 80 M/ml | {om_sph['L_star_survival_um']} µm | {om_sph['L_star_function_um']} µm | {om_sph['graft_volume_ml']} ml | **{om_sph['n_spheres_at_Rstar']}** spheres at $R=0.9 L^*_{{\\mathrm{{func}}}}$ |

A credit-card is ~46 cm². A **~{om_plan['planar_area_two_sided_cm2']:.0f} cm²** omental sheet is large but surgically thinkable.  
A **~{sq_plan['planar_area_one_sided_cm2']:.0f} cm²** one-sided SQ plaque is why thick subcutaneous macrodevices fail clinically: people thicken the slab instead of growing area, and $L$ crosses $L^*$.

## Literature-typical devices vs $L^*$

{md(devices)}

Reading: Papas-like **400 µm** capsules at 50 M/ml SQ sit **above** functional $L^*$ (~190 µm class) — core hypoxia is not a surprise, it is the bound. The recommended omental 200 µm / 80 M/ml sphere is the first row that is oxygen-functional *and* not portal-naked.

## Mouse vs human (why papers look better than trials)

| | Mouse | Human |
|--|-------|-------|
| $V_{{\\max}}$ used here | $1.2\\times 10^{{-16}}$ | $1.5\\times 10^{{-16}}$ (we certified on human) |
| Time to vessels | ~7–10 d | ~14–21 d |
| Device scale | mm | must cover Edmonton mass → tens of cm² |

$L^*\\propto 1/\\sqrt{{Q}}$. Human OCR + longer avascular window + 1000× more cells is a **quantitative** reason encapsulation “works in mice” and dies in SQ human devices. Not a new slogan — the same formula.

## FBR (order-of-magnitude)

A 50 µm fibrotic shell does not change $Q$, it drops boundary $pO_2$. Cutting $p_b$ from 30 → 20 mmHg shrinks planar $L^*_{{\\mathrm{{func}}}}$ by $\\sqrt{{(20-5)/(30-5)}}=\\sqrt{{0.6}}\\approx 0.77$. Certified thickness falls ~23%; required area rises ~30%. FBR is an **area tax**, not a separate mystery.

## What a lab can use tomorrow

1. Do not design planar SQ thicker than **~$L^*_{{\\mathrm{{func}}}}$** at the intended density (50 M/ml SQ ≈ **{sq_plan['L_star_function_um']} µm** half-thickness before safety factor).  
2. If you need Edmonton mass, budget **area first** (table), then print/pouch that area.  
3. Portal oxygen can pass $L^*$ and still lose the graft to IBMIR — heparin or leave the blood.  
4. Report **functional** core $pO_2$ (GSIS), not only live/dead.

## Reproduce

```bash
python3 t1d_simulator/dose_transport.py
```

Oxygen certificate details: `reports/design_bounds/verified_design_space.md`.
"""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return path


if __name__ == "__main__":
    out = run_study()
    print(json.dumps({"cells_millions": out["load"]["total_cells_millions"], "n_rows": out["n_rows"], "report": out["report"]}, indent=2))
