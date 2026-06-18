# -*- coding: utf-8 -*-
import sys
import numpy as np
from simulator import solve_oxygen_profile, HYDROGELS
import mesh_generator

def run_tests():
    print("=== Запуск автоматических тестов мультигеометрического цифрового двойника ===")
    
    D_alginate = HYDROGELS["alginate_2%"]["D"]
    
    # Тест 1: Проверка планарной модели (L = 100 мкм, p_boundary = 40 mmHg)
    print("Тест 1: Проверка планарной модели (L = 100 мкм, p_boundary = 40 mmHg)...")
    res_planar = solve_oxygen_profile(
        R_outer_microns=100,
        rho_million_per_ml=50,
        p_boundary=40.0,
        D_coefficient=D_alginate,
        geometry="planar"
    )
    
    assert np.isclose(res_planar["viable_fraction"], 100.0, atol=1e-2), f"Выживаемость должна быть близка к 100%, получено: {res_planar['viable_fraction']}%"
    assert np.isclose(res_planar["pO2"][-1], 40.0, atol=1e-2), f"Граничное условие справа не соблюдено: {res_planar['pO2'][-1]}"
    print("  [OK] Планарная модель рассчитана успешно.")
    
    # Тест 2: Сравнение геометрий (Планарная vs Цилиндрическая vs Сферическая)
    R_test = 200 # 200 мкм радиус/полутолщина
    rho_test = 100 # 100 млн/мл плотность клеток
    pO2_test = 30.0 # подкожное введение
    
    print(f"\nТест 2: Сравнение геометрий при R = {R_test} мкм, плотности = {rho_test} млн/мл и pO2 = {pO2_test} mmHg...")
    
    res_p = solve_oxygen_profile(R_test, rho_test, pO2_test, D_alginate, geometry="planar")
    res_c = solve_oxygen_profile(R_test, rho_test, pO2_test, D_alginate, geometry="cylindrical")
    res_s = solve_oxygen_profile(R_test, rho_test, pO2_test, D_alginate, geometry="spherical")
    
    print(f"  Выживаемость (Планарная): {res_p['viable_fraction']:.1f}%")
    print(f"  Выживаемость (Цилиндрическая): {res_c['viable_fraction']:.1f}%")
    print(f"  Выживаемость (Сферическая): {res_s['viable_fraction']:.1f}%")
    
    assert res_s["viable_fraction"] > res_c["viable_fraction"], "Сферическая выживаемость должна быть строго выше цилиндрической"
    assert res_c["viable_fraction"] > res_p["viable_fraction"], "Цилиндрическая выживаемость должна быть строго выше планарной"
    assert res_p["viable_fraction"] < 55.0, f"Плоская мембрана должна страдать от тяжелой гипоксии при R={R_test}, получено: {res_p['viable_fraction']}%"
    print("  [OK] Физический закон SA/V подтвержден: Сфера > Цилиндр > Лист.")
    
    # Тест 3: Генератор мешей и экспортер STL
    print("\nТест 3: Тестирование генератора CAD-сеток (mesh_generator.py)...")
    
    # Проверка коробки (slab)
    v_box, f_box = mesh_generator.generate_box_mesh(L_microns=150)
    assert len(v_box) == 8, "Коробка должна иметь 8 вершин"
    assert len(f_box) == 12, "Коробка должна иметь 12 треугольников"
    
    # Проверка цилиндра
    v_cyl, f_cyl = mesh_generator.generate_cylinder_mesh(R_microns=200, num_segments=16)
    # 16 вершин снизу + 16 сверху + 2 центра = 34 вершины
    assert len(v_cyl) == 34, f"Цилиндр с 16 сегментами должен иметь 34 вершины, получено {len(v_cyl)}"
    # 16*2 боковые + 16 нижний кап + 16 верхний кап = 64 треугольника
    assert len(f_cyl) == 64, f"Цилиндр должен иметь 64 треугольника, получено {len(f_cyl)}"
    
    # Проверка экспортера STL
    stl_str = mesh_generator.export_to_stl_ascii(v_box, f_box, solid_name="test_solid")
    assert isinstance(stl_str, str), "STL должен быть строкой"
    assert stl_str.startswith("solid test_solid"), "STL должен начинаться с 'solid test_solid'"
    assert stl_str.endswith("endsolid test_solid"), "STL должен оканчиваться на 'endsolid test_solid'"
    print("  [OK] Меши и STL-файлы генерируются корректно.")
    
    # Тест 4: Проверка GNN-пайплайна антифиброзных покрытий
    print("\nТест 4: Тестирование GNN-пайплайна (SMILES parsing, GNN forward pass, ranking)...")
    import torch
    import os
    from gnn_pipeline import smiles_to_graph, BiocompatibilityGNN
    
    # 1. SMILES -> Graph parsing
    s_zwitter = "C=CC(=O)NCC[N+](C)(C)CCCS(=O)(=O)[O-]"
    s_hydrophobic = "CC(=C)C(=O)OC"
    
    g_zwitter = smiles_to_graph(s_zwitter)
    g_hydrophobic = smiles_to_graph(s_hydrophobic)
    
    assert g_zwitter is not None, "Не удалось распарсить SMILES цвиттер-иона"
    assert g_hydrophobic is not None, "Не удалось распарсить SMILES гидрофобного мономера"
    
    # Проверка размерности признаков вершин (11 признаков на атом)
    assert g_zwitter.x.shape[1] == 11, f"Ожидалось 11 признаков вершин, получено {g_zwitter.x.shape[1]}"
    assert g_zwitter.edge_index.shape[0] == 2, "Индекс ребер должен иметь форму (2, E)"
    print("  [OK] Конвертация SMILES в граф выполнена корректно. Размерности признаков верны.")
    
    # 2. Инициализация и прямой проход (forward pass) GNN
    model = BiocompatibilityGNN()
    model.eval()
    
    # Загружаем сохраненные веса, если они есть
    weights_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "biocompatibility_gnn.pt")
    if os.path.exists(weights_path):
        model.load_state_dict(torch.load(weights_path))
        print("  [OK] Загружены веса обученной модели для тестирования.")
    else:
        print("  [Warning] Веса обученной модели не найдены, используется случайная инициализация.")
        
    with torch.no_grad():
        b_zwitter = torch.zeros(g_zwitter.x.size(0), dtype=torch.long)
        pred_zwitter = model(g_zwitter.x, g_zwitter.edge_index, b_zwitter)
        
        b_hydrophobic = torch.zeros(g_hydrophobic.x.size(0), dtype=torch.long)
        pred_hydro = model(g_hydrophobic.x, g_hydrophobic.edge_index, b_hydrophobic)
        
    score_zwitter = float(pred_zwitter[0, 0].item())
    score_hydro = float(pred_hydro[0, 0].item())
    
    # Проверка интервала предсказаний
    assert 0.0 <= score_zwitter <= 1.0, f"Предсказание GNN {score_zwitter} должно быть в [0, 1]"
    assert 0.0 <= score_hydro <= 1.0, f"Предсказание GNN {score_hydro} должно быть в [0, 1]"
    
    # Если веса были загружены, проверяем биологический закон ранжирования
    if os.path.exists(weights_path):
        assert score_zwitter > score_hydro, f"Ранжирование неверно: Цвиттер-ион ({score_zwitter:.4f}) должен быть выше гидрофобного мономера ({score_hydro:.4f})"
        print("  [OK] Биологическое ранжирование (Цвиттер-ион > Гидрофобный мономер) подтверждено.")
    else:
        print("  [Info] Ранжирование пропущено, так как веса модели случайны.")
    
    # Тест 5: Проверка сопряженной симуляции неоваскуляризации (VEGF)
    print("\nТест 5: Тестирование сопряженного ангиогенеза (VEGF diffusion & feedback)...")
    from simulator import run_neovascularization_sweep_oxygen
    
    # Случай А: VEGF не загружен (V_loaded = 0.0) -> p_boundary должно оставаться равным p_base = 30.0 mmHg
    res_no_vegf = run_neovascularization_sweep_oxygen(
        R_outer_microns=150,
        rho_million_per_ml=80,
        D_oxygen_coefficient=1.5e-5,
        geometry="spherical",
        V_loaded_relative=0.0,
        days=10
    )
    assert np.allclose(res_no_vegf["p_boundary"], 30.0), "При отсутствии VEGF давление кислорода не должно меняться!"
    assert np.allclose(res_no_vegf["C_interface"], 0.0), "При отсутствии VEGF концентрация на границе должна быть нулевой!"
    print("  [OK] Контрольный сценарий без VEGF пройден.")
    
    # Случай Б: VEGF загружен (V_loaded = 2.0) -> p_boundary должно возрастать со временем
    res_with_vegf = run_neovascularization_sweep_oxygen(
        R_outer_microns=150,
        rho_million_per_ml=80,
        D_oxygen_coefficient=1.5e-5,
        geometry="spherical",
        V_loaded_relative=2.0,
        days=10,
        beta_angiogenesis=0.2,
        K_vegf=0.1,
        p_base=30.0,
        p_max=60.0
    )
    
    assert res_with_vegf["C_interface"][-1] > 0.0, f"VEGF должен диффундировать на границу, получено: {res_with_vegf['C_interface'][-1]}"
    assert res_with_vegf["p_boundary"][-1] > 30.1, f"Давление кислорода должно возрасти благодаря неоваскуляризации, получено: {res_with_vegf['p_boundary'][-1]}"
    assert res_with_vegf["p_boundary"][-1] <= 60.0, "Давление кислорода не может превысить p_max"
    assert res_with_vegf["viability_over_time"][-1] >= res_with_vegf["viability_over_time"][0], "Выживаемость клеток должна возрастать по мере васкуляризации!"
    print(f"  [OK] Сценарий с VEGF пройден: pO2 выросло с {res_with_vegf['p_boundary'][0]:.1f} до {res_with_vegf['p_boundary'][-1]:.1f} mmHg.")
    
    # Тест 6: Проверка кальциевой токсичности и утечки инсулина (OGM)
    print("\nТест 6: Тестирование накопления кальция, эксайтотоксичности и утечки инсулина (OGM)...")
    res_bvp_ca = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", q_ogm_mmHg_per_sec=0.05, buffer_capacity_mM=10.0,
        catalase_activity_relative=1.0, t_days=3.0
    )
    
    assert res_bvp_ca["Ca_accum_mM"] > 1.2, f"Кальций должен накапливаться выше базального 1.2 ммоль/л, получено: {res_bvp_ca['Ca_accum_mM']:.2f}"
    assert res_bvp_ca["insulin_leak"] > 0.0, f"При накоплении кальция должна возникать нефункциональная утечка инсулина: {res_bvp_ca['insulin_leak']:.2f}"
    
    # Сравниваем с PINN для контроля BVP-PINN parity
    from pinn_solver import solve_oxygen_profile_pinn
    res_pinn_ca = solve_oxygen_profile_pinn(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", q_ogm_mmHg_per_sec=0.05, buffer_capacity_mM=10.0,
        catalase_activity_relative=1.0, t_days=3.0
    )
    assert np.isclose(res_bvp_ca["Ca_accum_mM"], res_pinn_ca["Ca_accum_mM"]), "Накопление кальция в BVP и PINN должно совпадать"
    assert np.isclose(res_bvp_ca["insulin_leak"], res_pinn_ca["insulin_leak"]), "Утечка инсулина в BVP и PINN должна совпадать"
    print("  [OK] Кальциевая токсичность и утечка инсулина подтверждены.")

    # Тест 7: Кинетика ковалентно-связанной каталазы
    print("\nТест 7: Тестирование ковалентного удержания каталазы vs вымывания при набухании...")
    # Свободная (вымываемая) каталаза при набухании 2.0
    res_free = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", q_ogm_mmHg_per_sec=0.05, swelling_ratio=2.0,
        catalase_activity_relative=1.0, catalase_half_life_days=1.5, t_days=5.0,
        tethered_catalase=False
    )
    # Ковалентно связанная каталаза
    res_tethered = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", q_ogm_mmHg_per_sec=0.05, swelling_ratio=2.0,
        catalase_activity_relative=1.0, catalase_half_life_days=1.5, t_days=5.0,
        tethered_catalase=True
    )
    
    # Так как свободная каталаза быстро вымывается (half-life = 1.5 / 4 = 0.375 дня), за 5 дней её не останется вовсе.
    # Ковалентно связанная каталаза имеет half-life = 100 дней, хотя её активность и ниже в 4 раза изначально.
    # Поэтому через 5 дней свободная каталаза допустит взрывной рост H2O2 по сравнению с ковалентной.
    assert res_free["H2O2"][0] > res_tethered["H2O2"][0], f"H2O2 со свободной каталазой ({res_free['H2O2'][0]:.1f}) должно быть выше, чем с ковалентно связанной ({res_tethered['H2O2'][0]:.1f}) из-за вымывания"
    print("  [OK] Ковалентное удержание каталазы защищает от вымывания при набухании.")

    # Тест 8: Каскад DAMPs при печати и механический разрыв геля
    print("\nТест 8: Тестирование каскада DAMPs и механического разрыва геля...")
    
    # 1. Проверка DAMPs (сопряженный транзиентный решатель цитокинов)
    from simulator import solve_cytokine_profile_transient
    res_coaxial = solve_cytokine_profile_transient(
        R_outer_microns=150, C_ext=10.0, D_cyt=1.0e-6, k_bind_scav=0.5,
        C_scav_0=10.0, k_deg_scav=0.01, swelling_ratio=1.0, days=5.0,
        shell_thickness_microns=50.0, coaxial_active=True
    )
    res_standard = solve_cytokine_profile_transient(
        R_outer_microns=150, C_ext=10.0, D_cyt=1.0e-6, k_bind_scav=0.5,
        C_scav_0=10.0, k_deg_scav=0.01, swelling_ratio=1.0, days=5.0,
        shell_thickness_microns=50.0, coaxial_active=False
    )
    
    assert res_coaxial["survival_shear"] > res_standard["survival_shear"], "Коаксиальный щит должен обеспечивать более высокую сдвиговую выживаемость клеток"
    assert res_standard["C_ext_timeline"][-1] > res_coaxial["C_ext_timeline"][-1], f"Концентрация внешних цитокинов при стандартной печати должна быть выше из-за DAMPs: std={res_standard['C_ext_timeline'][-1]:.1f}, coax={res_coaxial['C_ext_timeline'][-1]:.1f}"
    
    # 2. Проверка механического разрыва
    res_intact = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=50, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", swelling_ratio=1.0, E_0=50.0
    )
    res_ruptured = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=50, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", swelling_ratio=2.5, E_0=20.0
    )
    
    assert res_ruptured["rupture_risk"] > 50.0, f"Набухший мягкий гель должен иметь риск разрыва > 50%, получено: {res_ruptured['rupture_risk']:.1f}%"
    assert res_ruptured["viable_fraction"] < res_intact["viable_fraction"], "Выживаемость при разрыве должна падать из-за фактора разрыва геля"
    print("  [OK] Каскад DAMPs и риск механического разрыва верифицированы успешно.")
    
    
    # Тест 9: Эффект PFC буфера на увеличение выживаемости без кальциевой токсичности.
    print("\\nТест 9: Проверка истощающегося буфера PFC...")
    res_no_pfc = solve_oxygen_profile(
        R_outer_microns=600, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="planar", phi_pfc=0.0, t_days=0.0
    )
    res_pfc = solve_oxygen_profile(
        R_outer_microns=600, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="planar", phi_pfc=0.05, t_days=0.0
    )
    res_pfc_depleted = solve_oxygen_profile(
        R_outer_microns=600, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="planar", phi_pfc=0.05, t_days=20.0
    )
    assert res_pfc["viable_fraction"] > res_no_pfc["viable_fraction"], "PFC буфер должен увеличивать выживаемость на день 0."
    assert res_pfc["viable_fraction"] > res_pfc_depleted["viable_fraction"], "PFC буфер должен истощаться к 10-му дню."
    assert res_pfc["Ca_accum_mM"] <= 1.2, "PFC не должен накапливать кальций!"
    print("  [OK] Пассивный кислородный буфер PFC работает и истощается корректно.")

    # Тест 10: Спасение клеток от гипоксии за счет двухфазной предваскуляризации.
    print("\\nТест 10: Проверка двухфазной предваскуляризации...")
    res_no_pre = run_neovascularization_sweep_oxygen(
        R_outer_microns=600, rho_million_per_ml=80, D_oxygen_coefficient=D_alginate, geometry="planar",
        V_loaded_relative=2.0, days=14, t_pre_days=0.0
    )
    res_pre = run_neovascularization_sweep_oxygen(
        R_outer_microns=600, rho_million_per_ml=80, D_oxygen_coefficient=D_alginate, geometry="planar",
        V_loaded_relative=2.0, days=14, t_pre_days=14.0
    )
    assert res_no_pre["viability_over_time"][0] < 90.0, "Без предваскуляризации должна быть ранняя гипоксия."
    assert res_pre["viability_over_time"][0] > 99.9, "При предваскуляризации клетки отсутствуют, выживаемость 100%."
    assert res_pre["viability_over_time"][-1] > res_no_pre["viability_over_time"][0], "На момент инъекции клеток сосуды уже проросли, гипоксия предотвращена."
    print("  [OK] Окно смерти (гипоксия) предотвращено благодаря предваскуляризации.")

    # Тест 11: MMPs деградация геля и иммунный ответ макрофагов.
    print("\nТест 11: Гидролитическая деградация и клиренс DAMPs...")
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

    # Тест 12: Проверка NK-киллерного лизиса и лизиса системой комплемента
    print("\nТест 12: Проверка NK-киллерного лизиса и лизиса системой комплемента...")
    res_hypo_unprotected = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", crispr_hypoimmune=True, cd47_overexpression=False,
        complement_protection=False, t_days=5.0
    )
    res_hypo_protected = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", crispr_hypoimmune=True, cd47_overexpression=True,
        complement_protection=True, t_days=5.0
    )
    assert res_hypo_unprotected["viable_fraction"] < res_hypo_protected["viable_fraction"], "Клетки без CD47/комплемент защиты должны лизироваться сильнее"
    print("  [OK] Защита от NK и комплемента работает корректно.")

    # Тест 13: Проверка гемодинамического тромбоза AV-петли
    print("\nТест 13: Проверка гемодинамического тромбоза AV-петли...")
    res_normal = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", av_loop_flow=True, tau_blood=5.0, anticoagulation=False, t_days=10.0
    )
    res_shear_fail = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", av_loop_flow=True, tau_blood=0.5, anticoagulation=False, t_days=10.0
    )
    res_anticoag = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", av_loop_flow=True, tau_blood=5.0, anticoagulation=True, t_days=10.0
    )
    p_norm = res_normal["pO2"][-1]
    p_fail = res_shear_fail["pO2"][-1]
    p_anti = res_anticoag["pO2"][-1]
    assert p_fail < p_norm, f"Давление при нарушении сдвига ({p_fail:.2f}) должно быть ниже нормального ({p_norm:.2f})"
    assert p_anti > p_norm, f"Давление при антикоагуляции ({p_anti:.2f}) должно быть выше нормального ({p_norm:.2f})"
    print("  [OK] Тромбоз AV-петли и сдвиговое напряжение верифицированы.")

    # Тест 14: Проверка оксидативного гипероксического шока
    print("\nТест 14: Проверка оксидативного гипероксического шока...")
    res_norm_pfc = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", phi_pfc=0.05, pO2_pfc_saturation=150.0, t_days=1.0
    )
    res_high_pfc = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate,
        geometry="spherical", phi_pfc=0.05, pO2_pfc_saturation=400.0, t_days=1.0
    )
    assert res_high_pfc["viable_fraction"] < res_norm_pfc["viable_fraction"], "При гипероксии (>200 mmHg) выживаемость должна снижаться из-за f_hyperoxia"
    print("  [OK] Гипероксический шок верифицирован.")

    # Тест 15: Проверка фонового DAMPs от клеточного оборота
    print("\nТест 15: Проверка фонового DAMPs от клеточного оборота...")
    res_no_turnover = solve_cytokine_profile_transient(
        R_outer_microns=150, C_ext=10.0, D_cyt=1.0e-6, k_bind_scav=0.5,
        C_scav_0=10.0, k_deg_scav=0.01, swelling_ratio=1.0, days=14.0,
        shell_thickness_microns=50.0, coaxial_active=True, viable_fraction=80.0,
        turnover_rate=0.0
    )
    res_with_turnover = solve_cytokine_profile_transient(
        R_outer_microns=150, C_ext=10.0, D_cyt=1.0e-6, k_bind_scav=0.5,
        C_scav_0=10.0, k_deg_scav=0.01, swelling_ratio=1.0, days=14.0,
        shell_thickness_microns=50.0, coaxial_active=True, viable_fraction=80.0,
        turnover_rate=0.05
    )
    c_no = res_no_turnover["C_ext_timeline"][-1]
    c_yes = res_with_turnover["C_ext_timeline"][-1]
    assert c_yes > c_no, f"Фоновый DAMP при turnover_rate=0.05 должен повышать конечные цитокины ({c_yes:.2f} > {c_no:.2f})"
    print("  [OK] Фоновое выделение DAMPs подтверждено.")

    # Тест 16: Верификация иммунной утечки мини-органоидов (CRISPR защита)
    print("\nТест 16: Верификация иммунной утечки мини-органоидов...")
    from organoid_simulator import (
        simulate_organoid_population,
        simulate_organoid_oxygenation,
        simulate_organoid_insulin
    )
    t_span = np.linspace(0, 365, 365) # 1 год
    
    pop_unprotected = simulate_organoid_population(
        t_span, b2m_ko=False, ciita_ko=False, cd47_ki=False, cd55_cd59_ki=False, pdl1_ki=False
    )
    pop_protected = simulate_organoid_population(
        t_span, b2m_ko=True, ciita_ko=True, cd47_ki=True, cd55_cd59_ki=True, pdl1_ki=True
    )
    
    assert pop_unprotected[-1] < pop_protected[-1], f"Незащищенные клетки должны погибнуть сильнее защищенных ({pop_unprotected[-1]:.2f} vs {pop_protected[-1]:.2f})"
    print("  [OK] Влияние CRISPR редактирования на выживаемость клеток подтверждено.")

    # Тест 17: Проверка васкуляризации от EPCs
    print("\nТест 17: Проверка оксигенации от внутренней васкуляризации...")
    ox_low = simulate_organoid_oxygenation(t_days=5.0, phi_epc=0.0)
    ox_high = simulate_organoid_oxygenation(t_days=5.0, phi_epc=0.15)
    assert ox_high > ox_low, f"При добавлении EPCs оксигенация на день 5 должна быть выше ({ox_high:.1f} > {ox_low:.1f} mmHg)"
    print("  [OK] Сосудистая интеграция под влиянием EPCs верифицирована.")

    # Тест 18: Сравнение портального и системного инсулина (First-Pass Effect)
    print("\nТест 18: Сравнение портального и системного инсулина (First-Pass)...")
    ins_portal, ins_systemic = simulate_organoid_insulin(t_days=np.array([1.0]), N_cells=np.array([100.0]), pO2_core=np.array([40.0]))
    assert np.isclose(ins_systemic[0], ins_portal[0] * 0.4), f"Системный инсулин должен составлять 40% от портального из-за 60% экстракции, получено: {ins_systemic[0]:.2f} vs {ins_portal[0]:.2f}"
    print("  [OK] Портальный First-Pass эффект печени работает корректно.")

    print("\n=== Все тесты успешно пройдены! ===")
    return True

if __name__ == "__main__":
    try:
        run_tests()
        sys.exit(0)
    except AssertionError as e:
        print(f"\n  [ERROR] Тест провален: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\n  [ERROR] Непредвиденная ошибка при тестировании: {e}", file=sys.stderr)
        sys.exit(1)
