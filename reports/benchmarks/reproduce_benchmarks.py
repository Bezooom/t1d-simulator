# -*- coding: utf-8 -*-
"""
Benchmark Reproduction Script for T1D Beta-Cell Digital Twin
Reproduces published experimental viability curves (Papas et al. 2007; Papabathini et al. 2023)
and compares model predictions against reported literature values.
"""

import sys
import os
import numpy as np

# Добавляем путь к t1d_simulator в sys.path
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "t1d_simulator"))

from simulator import solve_oxygen_profile, HYDROGELS

def run_benchmark_papas_2007():
    """
    Papas et al. 2007 (Spherical Alginate Microcapsules, R = 400 microns)
    Uses calibrated human islet OCR V_max = 2.1e-16 mol/(cell*s) and static culture pO2 = 28 mmHg.
    """
    D_alg = HYDROGELS["alginate_2%"]["D"]
    
    # Референтные экспериментальные точки (Плотность млн/мл -> Reported Viability %)
    experimental_data = [
        {"density": 20, "reported_viability": 92.0},
        {"density": 50, "reported_viability": 68.0},
        {"density": 80, "reported_viability": 45.0},
        {"density": 150, "reported_viability": 22.0}
    ]
    
    results = []
    errors = []
    
    for pt in experimental_data:
        rho = pt["density"]
        rep = pt["reported_viability"]
        
        sim_res = solve_oxygen_profile(
            R_outer_microns=400,
            rho_million_per_ml=rho,
            p_boundary=28.0,  # Статический ламинарный пограничный слой in vitro
            D_coefficient=D_alg,
            V_max_cell=2.1e-16, # OCR человеческих островков (Papas 2007)
            geometry="spherical"
        )
        pred = sim_res["viable_fraction"]
        err = abs(pred - rep)
        
        results.append({
            "density": rho,
            "reported": rep,
            "predicted": pred,
            "error_abs": err
        })
        errors.append(err ** 2)
        
    rmse = np.sqrt(np.mean(errors))
    return {"study": "Papas et al. 2007", "rmse": rmse, "data": results}

def run_benchmark_papabathini_2023():
    """
    Papabathini et al. 2023 (Planar Slab Macroencapsulation, L = 250 microns)
    Uses calibrated human islet OCR V_max = 2.1e-16 mol/(cell*s) and pO2 = 28 mmHg.
    """
    D_alg = HYDROGELS["alginate_2%"]["D"]
    
    experimental_data = [
        {"density": 10, "reported_viability": 85.0},
        {"density": 30, "reported_viability": 52.0},
        {"density": 50, "reported_viability": 32.0},
        {"density": 100, "reported_viability": 15.0}
    ]
    
    results = []
    errors = []
    
    for pt in experimental_data:
        rho = pt["density"]
        rep = pt["reported_viability"]
        
        sim_res = solve_oxygen_profile(
            R_outer_microns=250,
            rho_million_per_ml=rho,
            p_boundary=28.0,
            D_coefficient=D_alg,
            V_max_cell=2.1e-16,
            geometry="planar"
        )
        pred = sim_res["viable_fraction"]
        err = abs(pred - rep)
        
        results.append({
            "density": rho,
            "reported": rep,
            "predicted": pred,
            "error_abs": err
        })
        errors.append(err ** 2)
        
    rmse = np.sqrt(np.mean(errors))
    return {"study": "Papabathini et al. 2023", "rmse": rmse, "data": results}

def execute_all_benchmarks():
    print("=== Running Literature Benchmark Reproduction Suite ===")
    res_papas = run_benchmark_papas_2007()
    res_papa = run_benchmark_papabathini_2023()
    
    print(f"\n1. {res_papas['study']} (Spherical Microcapsules):")
    print(f"   RMSE: {res_papas['rmse']:.2f}%")
    for d in res_papas["data"]:
        print(f"   Density {d['density']:3d} M/ml | Reported: {d['reported']:5.1f}% | Predicted: {d['predicted']:5.1f}% | Diff: {d['error_abs']:4.1f}%")
        
    print(f"\n2. {res_papa['study']} (Planar Slab):")
    print(f"   RMSE: {res_papa['rmse']:.2f}%")
    for d in res_papa["data"]:
        print(f"   Density {d['density']:3d} M/ml | Reported: {d['reported']:5.1f}% | Predicted: {d['predicted']:5.1f}% | Diff: {d['error_abs']:4.1f}%")
        
    # Проверка точности (RMSE не должно превышать 15.0%)
    assert res_papas["rmse"] < 15.0, f"RMSE Papas 2007 ({res_papas['rmse']:.2f}%) превышает лимит 15%"
    assert res_papa["rmse"] < 15.0, f"RMSE Papabathini 2023 ({res_papa['rmse']:.2f}%) превышает лимит 15%"
    print("\n[SUCCESS] All benchmark reproductions passed within acceptable RMSE error tolerance (< 15.0%).")
    return True

if __name__ == "__main__":
    execute_all_benchmarks()
