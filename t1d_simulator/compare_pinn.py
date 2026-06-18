import time
import numpy as np
from simulator import solve_oxygen_profile, HYDROGELS
from pinn_solver import solve_oxygen_profile_pinn

def run_comparison_test():
    print("================================================================================")
    print("🔬 СРАВНИТЕЛЬНОЕ ТЕСТИРОВАНИЕ: РЕШАТЕЛЬ BVP (SciPy) VS NEURAL SOLVER (PINN)")
    print("================================================================================")
    
    D_alginate = HYDROGELS["alginate_2%"]["D"]
    
    # Сценарии сравнения:
    # 1. Плоский лист, L=150 мкм, плотность=80 млн/мл, pO2=30 mmHg, без фиброза
    # 2. Цилиндрическое волокно, R=200 мкм, плотность=100 млн/мл, pO2=30 mmHg, фиброз 20 мкм
    # 3. Сферическая микрокапсула, R=250 мкм, плотность=120 млн/мл, pO2=30 mmHg, фиброз 50 мкм
    scenarios = [
        {
            "name": "1. Плоская мембрана (Slab)",
            "geometry": "planar",
            "R": 150,
            "rho": 80,
            "pO2": 30.0,
            "L_fib": 0.0,
            "D_fib": 1.0e-5
        },
        {
            "name": "2. Цилиндрическая нить (Fiber)",
            "geometry": "cylindrical",
            "R": 200,
            "rho": 100,
            "pO2": 30.0,
            "L_fib": 20.0,
            "D_fib": 1.0e-5
        },
        {
            "name": "3. Сферическая гранула (Microsphere)",
            "geometry": "spherical",
            "R": 250,
            "rho": 120,
            "pO2": 30.0,
            "L_fib": 50.0,
            "D_fib": 1.0e-5
        }
    ]
    
    print(f"\nПараметры среды: 2% Альгинат (D = {D_alginate:.1e} см²/с)\n")
    
    for sc in scenarios:
        print(f"🔹 Тестируем {sc['name']}...")
        print(f"   Параметры: R={sc['R']} мкм, плотность={sc['rho']} млн/мл, pO₂={sc['pO2']} mmHg, фиброз={sc['L_fib']} мкм")
        
        # Замер времени классического решателя
        t0 = time.time()
        res_bvp = solve_oxygen_profile(
            R_outer_microns=sc["R"],
            rho_million_per_ml=sc["rho"],
            p_boundary=sc["pO2"],
            D_coefficient=D_alginate,
            geometry=sc["geometry"],
            L_fibrosis_microns=sc["L_fib"],
            D_fibrosis=sc["D_fib"]
        )
        t_bvp = (time.time() - t0) * 1000.0 # в миллисекунды
        
        # Замер времени нейросетевого решателя
        t0 = time.time()
        res_pinn = solve_oxygen_profile_pinn(
            R_outer_microns=sc["R"],
            rho_million_per_ml=sc["rho"],
            p_boundary=sc["pO2"],
            D_coefficient=D_alginate,
            geometry=sc["geometry"],
            L_fibrosis_microns=sc["L_fib"],
            D_fibrosis=sc["D_fib"]
        )
        t_pinn = (time.time() - t0) * 1000.0 # в миллисекунды
        
        # Расчет ошибок по pO2 профилю
        mae = np.mean(np.abs(res_bvp["pO2"] - res_pinn["pO2"]))
        max_err = np.max(np.abs(res_bvp["pO2"] - res_pinn["pO2"]))
        
        # Вывод результатов сравнения
        print(f"   📊 Метрика точности:")
        print(f"      - Средняя абсолютная ошибка (MAE): {mae:.4f} mmHg (относительная: {mae / sc['pO2'] * 100.0:.2f}%)")
        print(f"      - Максимальная ошибка профиля: {max_err:.4f} mmHg")
        print(f"   📊 Метрики выживаемости и инсулина:")
        print(f"      - Выживаемость клеток:   BVP = {res_bvp['viable_fraction']:.1f}% | PINN = {res_pinn['viable_fraction']:.1f}%")
        print(f"      - Секреция инсулина:      BVP = {res_bvp['insulin_capacity']:.1f}% | PINN = {res_pinn['insulin_capacity']:.1f}%")
        print(f"      - Мин. pO₂ в центре:     BVP = {res_bvp['min_pO2']:.4f} mmHg | PINN = {res_pinn['min_pO2']:.4f} mmHg")
        print(f"   ⏱️ Время выполнения:")
        print(f"      - Классический SciPy (solve_bvp): {t_bvp:.2f} мс")
        print(f"      - Нейросеть PyTorch (PINN):       {t_pinn:.2f} мс (включая обучение)")
        print("-" * 80)
        
        # Небольшая проверка сходимости: MAE по pO2 должно быть меньше 1.5 mmHg
        assert mae < 1.5, f"Точность PINN вышла за допустимый предел! MAE = {mae:.4f} mmHg"
        
    print("\n✅ Сравнительное тестирование успешно завершено. PINN показывает высокую сходимость с классическими моделями.")

if __name__ == "__main__":
    run_comparison_test()
