import sys
import numpy as np
from simulator import solve_oxygen_profile
from pinn_solver import solve_oxygen_profile_pinn

def main():
    print("=== Тестирование точности PINN по токсикологии (pH и H2O2) ===")
    
    # Задаем сложные параметры для OGM, макрофагов, буфера, PLGA и т.д.
    params = {
        "R_outer_microns": 150.0,
        "rho_million_per_ml": 80.0,
        "p_boundary": 30.0,
        "D_coefficient": 1.5e-5,
        "geometry": "spherical",
        "L_fibrosis_microns": 50.0,
        "D_fibrosis": 0.3 * 3e-5,
        "rho_mac_million_per_ml": 50.0,
        "q_ogm_mmHg_per_sec": 0.02, # OGM активен
        "catalase_activity_relative": 1.0,
        "catalase_half_life_days": 1.5,
        "buffer_capacity_mM": 10.0,
        "swelling_ratio": 1.2,
        "plga_acidification_factor": 0.3,
        "t_days": 3.0
    }
    
    print("Запуск BVP-решателя (1D BVP)...")
    res_bvp = solve_oxygen_profile(**params)
    
    print("Запуск PINN-решателя (PyTorch)...")
    res_pinn = solve_oxygen_profile_pinn(**params)
    
    # Проверка совпадения сеток
    # Поскольку сетки могут слегка отличаться по количеству точек или координатам,
    # интерполируем PINN-решение на сетку BVP
    z_bvp = res_bvp["z"]
    pO2_bvp = res_bvp["pO2"]
    pH_bvp = res_bvp["pH"]
    H2O2_bvp = res_bvp["H2O2"]
    
    pO2_pinn_interp = np.interp(z_bvp, res_pinn["z"], res_pinn["pO2"])
    pH_pinn_interp = np.interp(z_bvp, res_pinn["z"], res_pinn["pH"])
    H2O2_pinn_interp = np.interp(z_bvp, res_pinn["z"], res_pinn["H2O2"])
    
    # Расчет MAE
    mae_pO2 = np.mean(np.abs(pO2_bvp - pO2_pinn_interp))
    mae_pH = np.mean(np.abs(pH_bvp - pH_pinn_interp))
    mae_H2O2 = np.mean(np.abs(H2O2_bvp - H2O2_pinn_interp))
    
    print(f"MAE по pO2: {mae_pO2:.4f} mmHg")
    print(f"MAE по pH: {mae_pH:.4f}")
    print(f"MAE по H2O2: {mae_H2O2:.4f} uM")
    
    # Проверка порогов точности
    assert mae_pO2 < 2.5, f"Ошибка pO2 превышает лимит 2.5 mmHg: {mae_pO2:.4f}"
    assert mae_pH < 0.1, f"Ошибка pH превышает лимит 0.1: {mae_pH:.4f}"
    assert mae_H2O2 < 1.0, f"Ошибка H2O2 превышает лимит 1.0 uM: {mae_H2O2:.4f}"
    
    print("  [OK] Все точностные характеристики соответствуют физическому канону!")
    print("=== Тест успешно завершен! ===")

if __name__ == "__main__":
    try:
        main()
        sys.exit(0)
    except AssertionError as e:
        print(f"\n  [ERROR] Ошибка валидации точности: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\n  [ERROR] Внутреннее исключение: {e}", file=sys.stderr)
        sys.exit(1)
