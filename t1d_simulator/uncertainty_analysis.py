# -*- coding: utf-8 -*-
"""
Uncertainty & Sensitivity Analysis Module (WP-Uncertainty)
Calculates local parameter sensitivities and Tornado analysis for top-5 physical constants.
"""

import sys
import os
import pandas as pd
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from simulator import solve_oxygen_profile, HYDROGELS, IMPLANTATION_SITES

def run_tornado_sensitivity(
    base_radius_microns=200.0,
    base_density_million_per_ml=80.0,
    geometry="planar",
    variation_pct=20.0
):
    """
    Performs +/- variation_pct parameter sweeps on top-5 physical constants:
    1. V_max (Oxygen consumption rate)
    2. K_M (Michaelis-Menten constant)
    3. D_eff (Hydrogel O2 diffusion coefficient)
    4. pO2_boundary (Implantation site O2 partial pressure)
    5. L_fibrosis (Fibrotic capsule thickness)
    """
    base_D = HYDROGELS["alginate_2%"]["D"]
    base_pO2 = IMPLANTATION_SITES["omental_pouch"]["pO2"]
    
    # Baseline simulation
    base_res = solve_oxygen_profile(
        R_outer_microns=base_radius_microns,
        rho_million_per_ml=base_density_million_per_ml,
        p_boundary=base_pO2,
        D_coefficient=base_D,
        geometry=geometry
    )
    base_viability = base_res["viable_fraction"]

    results = []
    factor_low = 1.0 - (variation_pct / 100.0)
    factor_high = 1.0 + (variation_pct / 100.0)

    base_Vmax = 1.5e-16
    parameters = [
        {"name": "V_max (Metabolic Consumption)", "key": "V_max", "base": base_Vmax},
        {"name": "K_M (Michaelis Constant)", "key": "K_M", "base": 0.5},
        {"name": "D_eff (Diffusion Coefficient)", "key": "D", "base": base_D},
        {"name": "pO2_boundary (Tissue Oxygen)", "key": "pO2", "base": base_pO2},
        {"name": "L_fibrosis (Capsule Thickness)", "key": "fibrosis", "base": 20.0},
    ]

    for p in parameters:
        v_low = p["base"] * factor_low
        v_high = p["base"] * factor_high

        def _eval(val):
            D_val = val if p["key"] == "D" else base_D
            pO2_val = val if p["key"] == "pO2" else base_pO2
            fib_val = val if p["key"] == "fibrosis" else 0.0
            vmax_val = val if p["key"] == "V_max" else base_Vmax
            # K_M is module-level in simulator; vary via temporary monkeypatch when requested
            kwargs = dict(
                R_outer_microns=base_radius_microns,
                rho_million_per_ml=base_density_million_per_ml,
                p_boundary=pO2_val,
                D_coefficient=D_val,
                geometry=geometry,
                L_fibrosis_microns=fib_val,
                V_max_cell=vmax_val,
            )
            if p["key"] == "K_M":
                import simulator as sim_mod
                old_km = sim_mod.K_M
                sim_mod.K_M = val
                try:
                    return solve_oxygen_profile(**kwargs)
                finally:
                    sim_mod.K_M = old_km
            return solve_oxygen_profile(**kwargs)

        res_low = _eval(v_low)
        res_high = _eval(v_high)
        
        viab_low = res_low["viable_fraction"]
        viab_high = res_high["viable_fraction"]
        delta_viab = abs(viab_high - viab_low)
        
        results.append({
            "parameter": p["name"],
            "base_value": p["base"],
            "viability_low": viab_low,
            "viability_high": viab_high,
            "delta_viability": delta_viab,
            "sensitivity_index": (delta_viab / base_viability) if base_viability > 0 else 0.0,
        })
        
    df_results = pd.DataFrame(results).sort_values(by="delta_viability", ascending=False)
    
    print("\n=================================================================")
    print(f"      TORNADO PARAMETER SENSITIVITY ANALYSIS (±{variation_pct}%)   ")
    print("=================================================================")
    print(f"Baseline Viability: {base_viability:.2f}% | Geometry: {geometry.upper()}")
    print("-----------------------------------------------------------------")
    for idx, row in df_results.iterrows():
        print(f" {row['parameter']:<35} | Delta Viab: {row['delta_viability']:5.2f}% | Sensitivity: {row['sensitivity_index']:5.3f}")
    print("=================================================================\n")
    
    return {
        "base_viability": base_viability,
        "dataframe": df_results,
        "top_sensitive_parameter": df_results.iloc[0]["parameter"]
    }

if __name__ == "__main__":
    run_tornado_sensitivity()
