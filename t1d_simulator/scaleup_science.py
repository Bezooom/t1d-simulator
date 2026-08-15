# -*- coding: utf-8 -*-
"""
Remaining biotech bounds that a continuum capsule model hides:

1. Single-islet / pseudoislet scale (discrete cluster, not smear density).
2. Irreversible death in the avascular window (mouse 8 d vs human 18 d).
3. Published device sizes scored against L* (virtual review).
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from design_bounds import V_MAX_HUMAN, analytical_l_crit_um
from dose_transport import CELLS_PER_IEQ, SAFETY
from simulator import HYDROGELS, IMPLANTATION_SITES

P_SURVIVAL = 0.5
P_FUNCTION = 5.0
K_DEATH_PER_DAY = 0.35  # irreversible loss of anoxic volume, ~2-day half-life


def sphere_volume_cm3(radius_um):
    r_cm = radius_um * 1e-4
    return (4.0 / 3.0) * np.pi * r_cm**3


def islet_local_density_m_per_ml(diameter_um, cells=CELLS_PER_IEQ):
    """Mean intra-islet density if `cells` occupy a sphere of given diameter."""
    vol = sphere_volume_cm3(0.5 * diameter_um)
    return (cells / vol) / 1e6 if vol > 0 else float("inf")


def anoxic_volume_fraction(length_um, l_star_um, geometry):
    """Zero-order anoxic core fraction when length > L*."""
    if length_um <= 0:
        return 1.0
    ratio = l_star_um / length_um
    if ratio >= 1.0:
        return 0.0
    if geometry == "planar":
        return float(np.sqrt(max(0.0, 1.0 - ratio**2)))
    if geometry == "cylindrical":
        return float(max(0.0, 1.0 - ratio**2))
    # spherical: r_c^2 = R^2 - L*^2  →  V_anox/V = (r_c/R)^3
    return float(max(0.0, 1.0 - ratio**2) ** 1.5)


def remaining_after_window(f_anox, days, k=K_DEATH_PER_DAY):
    """Anoxic volume is lost irreversibly; oxygenated shell stays."""
    lost = f_anox * (1.0 - np.exp(-k * days))
    return float(max(0.0, 1.0 - lost))


def single_islet_table(
    diameters_um=(100.0, 150.0, 200.0, 250.0, 300.0),
    sites=("subcutaneous", "omental_pouch"),
    v_max=V_MAX_HUMAN,
):
    d_coeff = HYDROGELS["alginate_2%"]["D"]
    rows = []
    for d_um in diameters_um:
        r_um = 0.5 * d_um
        rho = islet_local_density_m_per_ml(d_um)
        for site in sites:
            p_b = IMPLANTATION_SITES[site]["pO2"]
            l_s = analytical_l_crit_um(rho, p_b, d_coeff, v_max, "spherical", p_crit=P_SURVIVAL)
            l_f = analytical_l_crit_um(rho, p_b, d_coeff, v_max, "spherical", p_crit=P_FUNCTION)
            f_an = anoxic_volume_fraction(r_um, l_s, "spherical")
            rows.append({
                "diameter_um": d_um,
                "radius_um": r_um,
                "local_density_M_per_ml": round(rho, 1),
                "site_key": site,
                "L_star_surv_um": round(l_s, 1),
                "L_star_func_um": round(l_f, 1),
                "single_islet_anoxic_frac": round(f_an, 3),
                "single_islet_ok_surv": r_um <= SAFETY * l_s,
                "single_islet_ok_func": r_um <= SAFETY * l_f,
            })
    return pd.DataFrame(rows)


def packing_table(
    capsule_radius_um=200.0,
    packing_densities=(20.0, 50.0, 80.0, 100.0),
    islet_diameter_um=150.0,
    site_key="subcutaneous",
    geometry="spherical",
    v_max=V_MAX_HUMAN,
):
    """Continuum ρ is a smear. True local ρ is islet density; packing fraction φ = ρ/ρ_islet."""
    d_coeff = HYDROGELS["alginate_2%"]["D"]
    p_b = IMPLANTATION_SITES[site_key]["pO2"]
    rho_islet = islet_local_density_m_per_ml(islet_diameter_um)
    rows = []
    for rho in packing_densities:
        phi = rho / rho_islet
        l_smear = analytical_l_crit_um(rho, p_b, d_coeff, v_max, geometry, p_crit=P_SURVIVAL)
        # capsule-scale uses smear Q; cluster-scale uses full islet Q at islet R
        r_islet = 0.5 * islet_diameter_um
        l_islet = analytical_l_crit_um(rho_islet, p_b, d_coeff, v_max, "spherical", p_crit=P_SURVIVAL)
        rows.append({
            "site_key": site_key,
            "geometry": geometry,
            "capsule_R_um": capsule_radius_um,
            "packing_M_per_ml": rho,
            "islet_diameter_um": islet_diameter_um,
            "islet_local_M_per_ml": round(rho_islet, 1),
            "packing_fraction": round(phi, 4),
            "smear_L_star_um": round(l_smear, 1),
            "smear_capsule_ok": capsule_radius_um <= SAFETY * l_smear,
            "islet_L_star_um": round(l_islet, 1),
            "islet_ok": r_islet <= SAFETY * l_islet,
            "both_scales_ok": (capsule_radius_um <= SAFETY * l_smear) and (r_islet <= SAFETY * l_islet),
        })
    return pd.DataFrame(rows)


def death_window_table(
    geometries=("planar", "spherical"),
    sites=("subcutaneous", "omental_pouch"),
    densities=(50.0, 80.0),
    lengths=(150.0, 200.0, 250.0),
    mouse_days=8.0,
    human_days=18.0,
    v_max=V_MAX_HUMAN,
):
    d_coeff = HYDROGELS["alginate_2%"]["D"]
    rows = []
    for site in sites:
        p_b = IMPLANTATION_SITES[site]["pO2"]
        for geom in geometries:
            for dens in densities:
                for length in lengths:
                    l_s = analytical_l_crit_um(dens, p_b, d_coeff, v_max, geom, p_crit=P_SURVIVAL)
                    f_an = anoxic_volume_fraction(length, l_s, geom)
                    rows.append({
                        "site_key": site,
                        "geometry": geom,
                        "density_M_per_ml": dens,
                        "length_um": length,
                        "L_star_um": round(l_s, 1),
                        "anoxic_frac": round(f_an, 3),
                        "remain_mouse_8d": round(remaining_after_window(f_an, mouse_days), 3),
                        "remain_human_18d": round(remaining_after_window(f_an, human_days), 3),
                        "remain_if_vegf_day7": round(remaining_after_window(f_an, 7.0), 3),
                    })
    return pd.DataFrame(rows)


def virtual_review_table():
    """Score published / typical constructs against L* (human V_max, 2% alginate)."""
    d_coeff = HYDROGELS["alginate_2%"]["D"]
    cases = [
        {
            "source": "Papas et al. 2007 (benchmark)",
            "geometry": "spherical",
            "site_note": "in vitro unstirred ~28 mmHg",
            "p_boundary": 28.0,
            "length_um": 400.0,
            "density_M_per_ml": 50.0,
            "reported_viability": 68.0,
        },
        {
            "source": "Papas et al. 2007 high density",
            "geometry": "spherical",
            "site_note": "in vitro unstirred ~28 mmHg",
            "p_boundary": 28.0,
            "length_um": 400.0,
            "density_M_per_ml": 80.0,
            "reported_viability": 45.0,
        },
        {
            "source": "Papabathini et al. 2023 (benchmark)",
            "geometry": "planar",
            "site_note": "in vitro / SQ-like 28 mmHg",
            "p_boundary": 28.0,
            "length_um": 250.0,
            "density_M_per_ml": 50.0,
            "reported_viability": 32.0,
        },
        {
            "source": "Papabathini et al. 2023 thin-load",
            "geometry": "planar",
            "site_note": "in vitro / SQ-like 28 mmHg",
            "p_boundary": 28.0,
            "length_um": 250.0,
            "density_M_per_ml": 10.0,
            "reported_viability": 85.0,
        },
        {
            "source": "Typical SQ macrodevice 500 um full (L=250)",
            "geometry": "planar",
            "site_note": "subcutaneous 30 mmHg",
            "p_boundary": 30.0,
            "length_um": 250.0,
            "density_M_per_ml": 50.0,
            "reported_viability": None,
        },
        {
            "source": "This study recommended omental sphere",
            "geometry": "spherical",
            "site_note": "omental 55 mmHg",
            "p_boundary": 55.0,
            "length_um": 200.0,
            "density_M_per_ml": 80.0,
            "reported_viability": None,
        },
    ]
    rows = []
    for c in cases:
        l_s = analytical_l_crit_um(
            c["density_M_per_ml"], c["p_boundary"], d_coeff, V_MAX_HUMAN, c["geometry"], p_crit=P_SURVIVAL
        )
        l_f = analytical_l_crit_um(
            c["density_M_per_ml"], c["p_boundary"], d_coeff, V_MAX_HUMAN, c["geometry"], p_crit=P_FUNCTION
        )
        f_an = anoxic_volume_fraction(c["length_um"], l_s, c["geometry"])
        above = c["length_um"] > l_s
        rows.append({
            "source": c["source"],
            "geometry": c["geometry"],
            "p_boundary_mmHg": c["p_boundary"],
            "length_um": c["length_um"],
            "density_M_per_ml": c["density_M_per_ml"],
            "L_star_surv_um": round(l_s, 1),
            "L_star_func_um": round(l_f, 1),
            "anoxic_frac": round(f_an, 3),
            "bound_says": "OVER L*" if above else "UNDER L*",
            "reported_viability": c["reported_viability"],
            "remain_human_18d": round(remaining_after_window(f_an, 18.0), 3),
        })
    return pd.DataFrame(rows)


def run_study(output_dir=None):
    if output_dir is None:
        output_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "reports",
            "scaleup",
        )
    os.makedirs(output_dir, exist_ok=True)

    islets = single_islet_table()
    pack = packing_table()
    death = death_window_table()
    review = virtual_review_table()

    islets.to_csv(os.path.join(output_dir, "single_islet_scale.csv"), index=False)
    pack.to_csv(os.path.join(output_dir, "packing_two_scale.csv"), index=False)
    death.to_csv(os.path.join(output_dir, "death_window.csv"), index=False)
    review.to_csv(os.path.join(output_dir, "virtual_review.csv"), index=False)
    _plot_islet_and_window(islets, death, output_dir)
    report = _write_report(islets, pack, death, review, output_dir)
    return {"report": report, "n_islets": int(len(islets)), "n_review": int(len(review))}


def _plot_islet_and_window(islets, death, output_dir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))

    sq = islets[islets.site_key == "subcutaneous"].sort_values("diameter_um")
    om = islets[islets.site_key == "omental_pouch"].sort_values("diameter_um")
    axes[0].plot(sq["diameter_um"], sq["radius_um"], "k--", label="islet radius R")
    axes[0].plot(sq["diameter_um"], sq["L_star_surv_um"], "o-", label="SQ L* surv")
    axes[0].plot(om["diameter_um"], om["L_star_surv_um"], "s-", label="omentum L* surv")
    axes[0].axvline(150, color="grey", lw=0.8, ls=":")
    axes[0].set_xlabel("Islet / organoid diameter (µm)")
    axes[0].set_ylabel("µm")
    axes[0].set_title("Single-cluster bound: R vs L*")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.3)

    sub = death[(death.site_key == "subcutaneous") & (death.geometry == "planar") & (death.density_M_per_ml == 50.0)]
    sub = sub.sort_values("length_um")
    axes[1].plot(sub["length_um"], 100 * sub["remain_mouse_8d"], "o-", label="mouse 8 d")
    axes[1].plot(sub["length_um"], 100 * sub["remain_human_18d"], "s-", label="human 18 d")
    axes[1].plot(sub["length_um"], 100 * sub["remain_if_vegf_day7"], "^--", label="vessels by day 7")
    axes[1].set_xlabel("Planar half-thickness (µm), SQ 50 M/ml")
    axes[1].set_ylabel("Remaining mass (%)")
    axes[1].set_title("Avascular window (k=0.35/d)")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.3)
    axes[1].set_ylim(0, 105)

    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "fig_islet_and_death_window.png"), dpi=160)
    plt.close(fig)


def _md(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def _write_report(islets, pack, death, review, output_dir):
    path = os.path.join(output_dir, "remaining_biotech_bounds.md")
    native = islets[(islets.diameter_um == 150.0)]
    native_sq = native[native.site_key == "subcutaneous"].iloc[0]
    native_om = native[native.site_key == "omental_pouch"].iloc[0]
    small = islets[(islets.diameter_um == 100.0) & (islets.site_key == "subcutaneous")].iloc[0]
    big = islets[(islets.diameter_um == 300.0) & (islets.site_key == "subcutaneous")].iloc[0]

    dw = death[
        (death.site_key == "subcutaneous")
        & (death.geometry == "planar")
        & (death.density_M_per_ml == 50.0)
        & (death.length_um == 250.0)
    ].iloc[0]
    dw_rec = death[
        (death.site_key == "omental_pouch")
        & (death.geometry == "spherical")
        & (death.density_M_per_ml == 80.0)
        & (death.length_um == 200.0)
    ].iloc[0]

    pack_show = pack[[
        "packing_M_per_ml", "packing_fraction", "smear_L_star_um",
        "smear_capsule_ok", "islet_L_star_um", "islet_ok", "both_scales_ok",
    ]]

    text = f"""# Three remaining biotech bounds

**Date:** 2026-08-15  
**Stack:** human $V_{{\\max}}$, 2% alginate, zero-order $L^*$, irreversible anoxic loss $k={K_DEATH_PER_DAY}$/day  
**Together with:** `design_bounds.py` (capsule $L^*$), `dose_transport.py` (Edmonton area)

Continuum smear density hides three things a cell-therapy group actually hits.

---

## 1. Discrete islets, not 50 M/ml jelly

One IEQ ($\\approx {CELLS_PER_IEQ:.0f}$ cells) in a 150 µm sphere is **not** 50 M/ml.  
Local density $\\approx$ **{native_sq['local_density_M_per_ml']} M/ml**.

| Diameter | Local ρ | SQ $L^*$ | Omentum $L^*$ | SQ anoxic core |
|----------|---------|----------|---------------|----------------|
| 100 µm | {small['local_density_M_per_ml']} | {small['L_star_surv_um']} µm | — | {small['single_islet_anoxic_frac']} |
| **150 µm (native IEQ)** | {native_sq['local_density_M_per_ml']} | **{native_sq['L_star_surv_um']} µm** | **{native_om['L_star_surv_um']} µm** | **{native_sq['single_islet_anoxic_frac']}** |
| 300 µm | {big['local_density_M_per_ml']} | {big['L_star_surv_um']} µm | — | {big['single_islet_anoxic_frac']} |

Native islet **radius = 75 µm**. SQ $L^*$ at intra-islet ρ is **{native_sq['L_star_surv_um']} µm** — a single unvascularized islet in SQ $pO_2$ already has a core. Omentum $L^*$ = **{native_om['L_star_surv_um']} µm > 75** — the same islet is oxygen-legal.

**Design rule:** make **pseudoislets ≤ 100–120 µm**, not 250–300 µm aggregates. Capsule $L^*$ on smear density can look PASS while every cluster is FAIL.

Packing fraction $\\phi = \\rho_{{\\mathrm{{smear}}}}/\\rho_{{\\mathrm{{islet}}}}$ at 150 µm IEQ (SQ, capsule $R=200$ µm):

{_md(pack_show)}

Both scales must pass. Smear-only certificates are necessary, not sufficient.

---

## 2. Death window before vessels (mouse vs human)

Anoxic volume is not a snapshot. It is lost while you wait for capillaries.

$$
f_{{\\mathrm{{anox}}}}=
\\begin{{cases}}
\\sqrt{{1-(L^*/L)^2}} & \\text{{planar}}\\\\
(1-(L^*/R)^2)^{{3/2}} & \\text{{sphere}}
\\end{{cases}}
\\qquad
m(T)=1-f_{{\\mathrm{{anox}}}}(1-e^{{-kT}})
$$

| | Mouse | Human | VEGF / prevascular pouch |
|--|-------|-------|--------------------------|
| Typical wait | **8 d** | **18 d** | vessels by ~**7 d** |

SQ planar 50 M/ml, $L=250$ µm (Papabathini-class): $f_{{\\mathrm{{anox}}}}={dw['anoxic_frac']}$  
remaining **mouse {100*dw['remain_mouse_8d']:.0f}%** vs **human {100*dw['remain_human_18d']:.0f}%** vs vessels@7d **{100*dw['remain_if_vegf_day7']:.0f}%**.

Recommended omental sphere 200 µm / 80 M/ml: $f_{{\\mathrm{{anox}}}}={dw_rec['anoxic_frac']}$ → remaining **human {100*dw_rec['remain_human_18d']:.0f}%** (bound already says no core).

That is the quantitative “works in mice, dies in SQ human devices” without a new slogan: same $L^*$, longer $T$, more mass already shown in `dose_transport.py`.

![Islet L* and death window](fig_islet_and_death_window.png)

---

## 3. Virtual review of published / typical sizes

Same $L^*$ against numbers we already calibrated (Papas RMSE 10.42%, Papabathini RMSE 11.42%) plus the recommended point.

{_md(review)}

**Reading.** Every OVER $L^*$ row has a large anoxic fraction and a low reported viability where the paper measured one. The 400 µm capsule and 250 µm planar slab are not near-misses — they sit **~2×** above $L^*$. The omental 200 µm / 80 M/ml sphere is the first listed design that is under both survival and function $L^*$.

We are not claiming we reproduced ViaCyte/Vertex internals (those thicknesses are not fully public). We are claiming: **any** planar SQ half-thickness ≥ 250 µm at ≥ 50 M/ml is forbidden by the same inequality that matches Papas/Papabathini.

---

## What to put in a protocol / preprint

1. Organoid **diameter ≤ 120 µm** (single-cluster $R < L^*$ even in SQ).  
2. Capsule / sheet **$L \\le 0.9 L^*_{{\\mathrm{{func}}}}$** at the *smear* density (previous reports).  
3. Edmonton mass → **budget area**, not thickness (`dose_transport.py`: ~650 cm² omental sheet at 80 M/ml).  
4. Do not wait 18 days avascular at SQ $pO_2$ with $L>L^*$ — prevascularize or go omentum.  
5. Portal still needs anti-IBMIR even if oxygen-certified.

## Reproduce

```bash
python3 t1d_simulator/scaleup_science.py
```
"""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return path


if __name__ == "__main__":
    out = run_study()
    print(json.dumps(out, indent=2))
