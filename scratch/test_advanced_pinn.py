import sys
import os
sys.path.append("t1d_simulator")
sys.path.append("scratch")

import numpy as np
from pinn_solver import solve_advanced_oxygen_profile_pinn, solve_cytokine_profile_pinn
from simulate_advanced_biodesign import solve_advanced_oxygen_profile_fd, solve_cytokine_profile

def run_pinn_advanced_test():
    print("=== Тестирование усовершенствованного PINN-решателя ===")
    
    # 1. Сравнение кислорода (BVP-FD vs PINN)
    print("Запуск сопоставления профилей кислорода (активное воспаление + OGM)...")
    res_fd = solve_advanced_oxygen_profile_fd(
        R_outer_microns=150,
        rho_million_per_ml=80,
        p_boundary_base=30.0,
        L_fibrosis_microns=50.0,
        rho_mac_million_per_ml=80.0,
        q_ogm_mmHg_per_sec=0.05
    )
    
    res_pinn = solve_advanced_oxygen_profile_pinn(
        R_outer_microns=150,
        rho_million_per_ml=80,
        p_boundary_base=30.0,
        L_fibrosis_microns=50.0,
        rho_mac_million_per_ml=80.0,
        q_ogm_mmHg_per_sec=0.05
    )
    
    mae_oxygen = np.mean(np.abs(np.interp(res_fd["r"], res_pinn["z"], res_pinn["pO2"]) - res_fd["p"]))
    print(f"  Выживаемость клеток: FD = {res_fd['viable_fraction']:.1f}% | PINN = {res_pinn['viable_fraction']:.1f}%")
    print(f"  Мин pO2: FD = {res_fd['min_pO2']:.2f} mmHg | PINN = {res_pinn['min_pO2']:.2f} mmHg")
    print(f"  Средняя абсолютная ошибка (MAE) профиля: {mae_oxygen:.4f} mmHg")
    
    # 2. Сравнение цитокинов (FD vs PINN)
    print("\nЗапуск сопоставления профилей цитокинов (диффузия + ловушки)...")
    res_cyt_fd = solve_cytokine_profile(
        R_outer_microns=150,
        C_ext=10.0,
        D_cyt=1.0e-6,
        k_bind_scav=0.5,
        k_deg=0.01
    )
    
    res_cyt_pinn = solve_cytokine_profile_pinn(
        R_outer_microns=150,
        C_ext=10.0,
        D_cyt=1.0e-6,
        k_bind_scav=0.5,
        k_deg=0.01
    )
    
    mae_cyt = np.mean(np.abs(np.interp(res_cyt_fd["r"], res_cyt_pinn["z"], res_cyt_pinn["C"]) - res_cyt_fd["C"]))
    print(f"  Защищенная фракция клеток: FD = {res_cyt_fd['protected_fraction']:.1f}% | PINN = {res_cyt_pinn['protected_fraction']:.1f}%")
    print(f"  Средняя абсолютная ошибка (MAE) профиля цитокинов: {mae_cyt:.4f} ng/ml")
    
    assert mae_oxygen < 2.0, "Ошибка по кислороду слишком велика!"
    assert mae_cyt < 1.0, "Ошибка по цитокинам слишком велика!"
    print("\n=== Все тесты усовершенствованного PINN успешно пройдены! ===")

if __name__ == "__main__":
    run_pinn_advanced_test()
