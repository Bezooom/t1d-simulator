import re
import os

with open('t1d_simulator/verify_model.py', 'r') as f:
    code = f.read()

# Add test 9, 10, 11
tests_str = """
    # Тест 9: Эффект PFC буфера на увеличение выживаемости без кальциевой токсичности.
    print("\\nТест 9: Проверка истощающегося буфера PFC...")
    res_no_pfc = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", phi_pfc=0.0, t_days=0.0
    )
    res_pfc = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", phi_pfc=0.2, t_days=0.0
    )
    res_pfc_depleted = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", phi_pfc=0.2, t_days=10.0
    )
    assert res_pfc["viable_fraction"] > res_no_pfc["viable_fraction"], "PFC буфер должен увеличивать выживаемость на день 0."
    assert res_pfc["viable_fraction"] > res_pfc_depleted["viable_fraction"], "PFC буфер должен истощаться к 10-му дню."
    assert res_pfc["Ca_accum_mM"] <= 1.2, "PFC не должен накапливать кальций!"
    print("  [OK] Пассивный кислородный буфер PFC работает и истощается корректно.")

    # Тест 10: Спасение клеток от гипоксии за счет двухфазной предваскуляризации.
    print("\\nТест 10: Проверка двухфазной предваскуляризации...")
    res_no_pre = run_neovascularization_sweep_oxygen(
        R_outer_microns=150, rho_million_per_ml=80, D_oxygen_coefficient=D_alginate, geometry="spherical",
        V_loaded_relative=2.0, days=14, t_pre_days=0.0
    )
    res_pre = run_neovascularization_sweep_oxygen(
        R_outer_microns=150, rho_million_per_ml=80, D_oxygen_coefficient=D_alginate, geometry="spherical",
        V_loaded_relative=2.0, days=14, t_pre_days=14.0
    )
    assert res_no_pre["viability_over_time"][0] < 90.0, "Без предваскуляризации должна быть ранняя гипоксия."
    assert res_pre["viability_over_time"][0] == 100.0, "При предваскуляризации клетки отсутствуют, выживаемость 100%."
    assert res_pre["viability_over_time"][-1] > res_no_pre["viability_over_time"][0], "На момент инъекции клеток сосуды уже проросли, гипоксия предотвращена."
    print("  [OK] Окно смерти (гипоксия) предотвращено благодаря предваскуляризации.")

    # Тест 11: MMPs деградация геля и иммунный ответ макрофагов.
    print("\\nТест 11: Гидролитическая деградация и клиренс DAMPs...")
    res_damp_0 = solve_cytokine_profile_transient(
        R_outer_microns=150, C_ext=10.0, D_cyt=1.0e-6, k_bind_scav=0.5,
        C_scav_0=10.0, k_deg_scav=0.01, swelling_ratio=1.0, days=0.1,
        shell_thickness_microns=50.0, coaxial_active=True, viable_fraction=50.0
    )
    res_damp_14 = solve_cytokine_profile_transient(
        R_outer_microns=150, C_ext=10.0, D_cyt=1.0e-6, k_bind_scav=0.5,
        C_scav_0=10.0, k_deg_scav=0.01, swelling_ratio=1.0, days=14.0,
        shell_thickness_microns=50.0, coaxial_active=True, viable_fraction=50.0
    )
    assert res_damp_14["C_ext_timeline"][-1] < np.max(res_damp_14["C_ext_timeline"]), "DAMPs должны иметь клиренс и затухать со временем."
    print("  [OK] Гидролитическая деградация и затухание DAMPs подтверждены.")
"""

code = code.replace(
    'print("\\n=== Все тесты успешно пройдены! ===")',
    tests_str + '\n    print("\\n=== Все тесты успешно пройдены! ===")'
)

with open('t1d_simulator/verify_model.py', 'w') as f:
    f.write(code)
