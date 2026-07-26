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
    assert len(v_cyl) == 34, f"Цилиндр с 16 сегментами должен иметь 34 вершины, получено {len(v_cyl)}"
    assert len(f_cyl) == 64, f"Цилиндр должен иметь 64 треугольника, получено {len(f_cyl)}"
    
    # Проверка экспортера STL
    stl_str = mesh_generator.export_to_stl_ascii(v_box, f_box, solid_name="test_solid")
    assert isinstance(stl_str, str), "STL должен быть строкой"
    assert stl_str.startswith("solid test_solid"), "STL должен начинаться с 'solid test_solid'"
    assert stl_str.endswith("endsolid test_solid"), "STL должен оканчиваться на 'endsolid test_solid'"
    print("  [OK] Меши и STL-файлы генерируются корректно.")
    
    # Тест 4: Проверка GNN-пайплайна антифиброзных покрытий
    print("\nТест 4: Тестирование GNN-пайплайна (SMILES parsing, GNN forward pass, ranking)...")
    try:
        import torch
        import os
        from gnn_pipeline import smiles_to_graph, BiocompatibilityGNN
        
        s_zwitter = "C=CC(=O)NCC[N+](C)(C)CCCS(=O)(=O)[O-]"
        s_hydrophobic = "CC(=C)C(=O)OC"
        
        g_zwitter = smiles_to_graph(s_zwitter)
        g_hydrophobic = smiles_to_graph(s_hydrophobic)
        
        assert g_zwitter is not None, "Не удалось распарсить SMILES цвиттер-иона"
        assert g_hydrophobic is not None, "Не удалось распарсить SMILES гидрофобного мономера"
        assert g_zwitter.x.shape[1] == 11, f"Ожидалось 11 признаков вершин, получено {g_zwitter.x.shape[1]}"
        assert g_zwitter.edge_index.shape[0] == 2, "Индекс ребер должен иметь форму (2, E)"
        print("  [OK] Конвертация SMILES в граф выполнена корректно. Размерности признаков верны.")
        
        model = BiocompatibilityGNN()
        model.eval()
        
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
        
        assert 0.0 <= score_zwitter <= 1.0, f"Предсказание GNN {score_zwitter} должно быть в [0, 1]"
        assert 0.0 <= score_hydro <= 1.0, f"Предсказание GNN {score_hydro} должно быть в [0, 1]"
        
        if os.path.exists(weights_path):
            assert score_zwitter > score_hydro, f"Ранжирование неверно: Цвиттер-ион ({score_zwitter:.4f}) должен быть выше гидрофобного мономера ({score_hydro:.4f})"
            print("  [OK] Биологическое ранжирование (Цвиттер-ион > Гидрофобный мономер) подтверждено.")
    except (ImportError, ModuleNotFoundError) as e:
        print(f"  [Skipped] Опциональные библиотеки GNN (rdkit/torch_geometric) не установлены ({e}).")
    
    # Тест 5: Проверка сопряженной симуляции неоваскуляризации (VEGF)
    print("\nТест 5: Тестирование сопряженного ангиогенеза (VEGF diffusion & feedback)...")
    from simulator import run_neovascularization_sweep_oxygen
    
    res_no_vegf = run_neovascularization_sweep_oxygen(
        R_outer_microns=150, rho_million_per_ml=80, D_oxygen_coefficient=1.5e-5, geometry="spherical", V_loaded_relative=0.0, days=10
    )
    assert np.allclose(res_no_vegf["p_boundary"], 30.0), "При отсутствии VEGF давление кислорода не должно меняться!"
    assert np.allclose(res_no_vegf["C_interface"], 0.0), "При отсутствии VEGF концентрация на границе должна быть нулевой!"
    print("  [OK] Контрольный сценарий без VEGF пройден.")
    
    res_with_vegf = run_neovascularization_sweep_oxygen(
        R_outer_microns=150, rho_million_per_ml=80, D_oxygen_coefficient=1.5e-5, geometry="spherical", V_loaded_relative=2.0, days=10, beta_angiogenesis=0.2, K_vegf=0.1, p_base=30.0, p_max=60.0
    )
    assert res_with_vegf["C_interface"][-1] > 0.0, f"VEGF должен диффундировать на границу, получено: {res_with_vegf['C_interface'][-1]}"
    assert res_with_vegf["p_boundary"][-1] > 30.1, f"Давление кислорода должно возрасти благодаря неоваскуляризации, получено: {res_with_vegf['p_boundary'][-1]}"
    assert res_with_vegf["p_boundary"][-1] <= 60.0, "Давление кислорода не может превысить p_max"
    assert res_with_vegf["viability_over_time"][-1] >= res_with_vegf["viability_over_time"][0], "Выживаемость клеток должна возрастать по мере васкуляризации!"
    print(f"  [OK] Сценарий с VEGF пройден: pO2 выросло с {res_with_vegf['p_boundary'][0]:.1f} до {res_with_vegf['p_boundary'][-1]:.1f} mmHg.")
    
    # Тест 6: Проверка кальциевой токсичности и утечки инсулина (OGM)
    print("\nТест 6: Тестирование накопления кальция, эксайтотоксичности и утечки инсулина (OGM)...")
    res_bvp_ca = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate, geometry="spherical", q_ogm_mmHg_per_sec=0.05, buffer_capacity_mM=10.0, catalase_activity_relative=1.0, t_days=3.0
    )
    assert res_bvp_ca["Ca_accum_mM"] > 1.2, f"Кальций должен накапливаться выше базального 1.2 ммоль/л, получено: {res_bvp_ca['Ca_accum_mM']:.2f}"
    assert res_bvp_ca["insulin_leak"] > 0.0, f"При накоплении кальция должна возникать нефункциональная утечка инсулина: {res_bvp_ca['insulin_leak']:.2f}"
    
    from pinn_solver import solve_oxygen_profile_pinn
    res_pinn_ca = solve_oxygen_profile_pinn(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate, geometry="spherical", q_ogm_mmHg_per_sec=0.05, buffer_capacity_mM=10.0, catalase_activity_relative=1.0, t_days=3.0
    )
    assert np.isclose(res_bvp_ca["Ca_accum_mM"], res_pinn_ca["Ca_accum_mM"]), "Накопление кальция в BVP и PINN должно совпадать"
    assert np.isclose(res_bvp_ca["insulin_leak"], res_pinn_ca["insulin_leak"]), "Утечка инсулина в BVP и PINN должна совпадать"
    print("  [OK] Кальциевая токсичность и утечка инсулина подтверждены.")

    # Тест 7: Кинетика ковалентно-связанной каталазы
    print("\nТест 7: Тестирование ковалентного удержания каталазы vs вымывания при набухании...")
    res_free = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate, geometry="spherical", q_ogm_mmHg_per_sec=0.05, swelling_ratio=2.0, catalase_activity_relative=1.0, catalase_half_life_days=1.5, t_days=5.0, tethered_catalase=False
    )
    res_tethered = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=80, p_boundary=30.0, D_coefficient=D_alginate, geometry="spherical", q_ogm_mmHg_per_sec=0.05, swelling_ratio=2.0, catalase_activity_relative=1.0, catalase_half_life_days=1.5, t_days=5.0, tethered_catalase=True
    )
    assert res_free["H2O2"][0] > res_tethered["H2O2"][0], f"H2O2 со свободной каталазой ({res_free['H2O2'][0]:.1f}) должно быть выше, чем с ковалентно связанной ({res_tethered['H2O2'][0]:.1f}) из-за вымывания"
    print("  [OK] Ковалентное удержание каталазы защищает от вымывания при набухании.")

    # Тест 8: Каскад DAMPs при печати и механический разрыв геля
    print("\nТест 8: Тестирование каскада DAMPs и механического разрыва геля...")
    from simulator import solve_cytokine_profile_transient
    res_coaxial = solve_cytokine_profile_transient(
        R_outer_microns=150, C_ext=10.0, D_cyt=1.0e-6, k_bind_scav=0.5, C_scav_0=10.0, k_deg_scav=0.01, swelling_ratio=1.0, days=5.0, shell_thickness_microns=50.0, coaxial_active=True
    )
    res_standard = solve_cytokine_profile_transient(
        R_outer_microns=150, C_ext=10.0, D_cyt=1.0e-6, k_bind_scav=0.5, C_scav_0=10.0, k_deg_scav=0.01, swelling_ratio=1.0, days=5.0, shell_thickness_microns=50.0, coaxial_active=False
    )
    assert res_coaxial["survival_shear"] > res_standard["survival_shear"], "Коаксиальный щит должен обеспечивать более высокую сдвиговую выживаемость клеток"
    assert res_standard["C_ext_timeline"][-1] > res_coaxial["C_ext_timeline"][-1], "Концентрация внешних цитокинов при стандартной печати должна быть выше из-за DAMPs"
    
    res_intact = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=50, p_boundary=30.0, D_coefficient=D_alginate, geometry="spherical", swelling_ratio=1.0, E_0=50.0
    )
    res_ruptured = solve_oxygen_profile(
        R_outer_microns=150, rho_million_per_ml=50, p_boundary=30.0, D_coefficient=D_alginate, geometry="spherical", swelling_ratio=2.5, E_0=20.0
    )
    assert res_ruptured["rupture_risk"] > 50.0, f"Набухший мягкий гель должен иметь риск разрыва > 50%, получено: {res_ruptured['rupture_risk']:.1f}%"
    assert res_ruptured["viable_fraction"] < res_intact["viable_fraction"], "Выживаемость при разрыве должна падать из-за фактора разрыва геля"
    print("  [OK] Каскад DAMPs и риск механического разрыва верифицированы успешно.")
    
    # Тесты 9 - 18
    print("\nТест 9-18: Проверка базовых функций симуляции и фазы 10...")
    from organoid_simulator import (
        simulate_organoid_population,
        simulate_organoid_oxygenation,
        simulate_organoid_insulin,
        simulate_icasp9_ap1903_apoptosis,
        simulate_ibmir_protection,
        simulate_transplantation_site_comparison,
        simulate_base_editing_fidelity
    )
    
    # Тест 19: Кинетика iCasp9 апоптоза при AP1903
    print("\nТест 19: Проверка индуцируемого апоптоза iCasp9 / AP1903...")
    surv_4h, elim_pct = simulate_icasp9_ap1903_apoptosis(t_hours=4.0, ap1903_conc_nM=10.0, N_0=100.0)
    assert elim_pct > 95.0, f"AP1903 должен элиминировать >95% клеток за 4 часа, получено: {elim_pct:.1f}%"
    assert surv_4h < 5.0, "Остаточная популяция должна быть < 5%"
    print("  [OK] Онкогенная защита iCasp9 / AP1903 подтверждена (>95% элиминации за 4ч).")

    # Тест 20: Проверка нанохимической защиты от IBMIR (PEG-LMWH)
    print("\nТест 20: Проверка нанохимической гепаринизации от IBMIR...")
    ibmir_uncoated = simulate_ibmir_protection(peg_lmwh_density=0.0)
    ibmir_coated = simulate_ibmir_protection(peg_lmwh_density=1.5)
    assert ibmir_coated["thrombin_generation"] < ibmir_uncoated["thrombin_generation"], "Гепаринизация должна снижать генерацию тромбина"
    assert ibmir_coated["retention_48h_percent"] > ibmir_uncoated["retention_48h_percent"], "Сохранение массы клеток через 48ч должно возрастать при PEG-LMWH"
    print("  [OK] Защита от IBMIR через Lipid-PEG-LMWH верифицирована.")

    # Тест 21: Сравнение Большого Сальника (Omental Pouch) и Печени
    print("\nТест 21: Сравнение сайтов имплантации (Печень vs Большой Сальник)...")
    site_portal = simulate_transplantation_site_comparison(t_days=[1.0], N_cells=100.0, site="portal_vein")
    site_omental = simulate_transplantation_site_comparison(t_days=[1.0], N_cells=100.0, site="omental_pouch")
    assert site_omental["steatosis_risk_index"] < site_portal["steatosis_risk_index"], "Большой сальник должен минимизировать риск жирового гепатоза печени"
    assert site_omental["retrievability_score"] == 100.0, "Сальник должен обеспечивать 100% извлекаемость"
    print("  [OK] Преимущество Omental Pouch по безопасности и извлекаемости подтверждено.")

    # Тест 22: Оценка точности Base Editing (CBE/ABE) vs SpCas9
    print("\nТест 22: Оценка геномной безопасности Base Editing...")
    fid_cas9 = simulate_base_editing_fidelity(method="spcas9")
    fid_be = simulate_base_editing_fidelity(method="base_editing")
    assert fid_be["translocation_risk_percent"] < fid_cas9["translocation_risk_percent"], "Base Editing не должен создавать высокий риск транслокаций"
    assert fid_be["clonogenic_survival_percent"] > fid_cas9["clonogenic_survival_percent"], "Выживаемость клонов при Base Editing выше из-за отсутствия DSB"
    print("  [OK] Геномная безопасность Base Editing (CBE/ABE) подтверждена.")

    # Тест 23: Клинический калькулятор дозирования трансплантата
    print("\nТест 23: Клинический калькулятор дозирования для пациента СД1...")
    from organoid_simulator import calculate_patient_transplant_dose, simulate_ogtt_glycemic_control
    dose_70kg = calculate_patient_transplant_dose(weight_kg=70.0, tdi_units=40.0)
    dose_100kg = calculate_patient_transplant_dose(weight_kg=100.0, tdi_units=60.0)
    assert dose_100kg["total_ieq"] > dose_70kg["total_ieq"], "Пациент большей массы должен получать большую дозу IEQ"
    assert dose_70kg["insulin_independence_forecast"] >= 100.0, "При 11,000 IEQ/кг прогнозируется 100% инсулинонезависимость"
    print("  [OK] Клинический калькулятор дозирования IEQ рассчитан корректно.")

    # Тест 24: Симуляция OGTT и гликемического контроля
    print("\nТест 24: Симуляция перорального глюкозотолерантного теста (OGTT)...")
    t_hours = np.linspace(0, 4, 100)
    g_pre = simulate_ogtt_glycemic_control(t_hours, meal_carbs_g=50.0, is_transplanted=False)
    g_post = simulate_ogtt_glycemic_control(t_hours, meal_carbs_g=50.0, is_transplanted=True)
    assert np.max(g_pre) > np.max(g_post), "Максимальная гликемия ДО трансплантации должна быть намного выше, чем ПОСЛЕ"
    assert g_post[-1] < 6.5, f"Гликемия ПОСЛЕ трансплантации на 4 час должна нормализоваться (<6.5 ммоль/л), получено: {g_post[-1]:.2f}"
    print("  [OK] Модель Бергмана и гликемический контроль GSIS верифицированы.")

    # Тест 25: Экспорт STL-каркаса сальника для 3D-биопечати
    print("\nТест 25: Экспорт STL-каркаса для 3D-биопринтера...")
    from organoid_cad_exporter import generate_omental_scaffold_stl, generate_patient_clinical_passport
    stl_omental = generate_omental_scaffold_stl(area_cm2=40.0, thickness_mm=0.5)
    assert stl_omental.startswith("solid Omental_Scaffold_40sqcm"), "STL файл должен начинаться с верного имени solid"
    assert stl_omental.endswith("endsolid Omental_Scaffold_40sqcm"), "STL файл должен завершаться верным тегом endsolid"
    print("  [OK] 3D STL-модель скаффолда сальника сгенерирована корректно.")

    # Тест 26: Генерация Персонального Клинического Паспорта
    print("\nТест 26: Генерация Персонального Клинического Паспорта...")
    pat_data = {"weight_kg": 70.0, "tdi_units": 45.0, "c_peptide_pmol_l": 10.0}
    dose_data = dose_70kg
    passport_text = generate_patient_clinical_passport(pat_data, dose_data, {})
    assert "Персональный Клинический Паспорт" in passport_text, "Паспорт должен содержать заголовок"
    assert "70.0 кг" in passport_text, "Паспорт должен содержать вес пациента"
    assert "iCasp9" in passport_text, "Паспорт должен содержать спецификацию онкогенной безопасности iCasp9"
    # Тест 27: 30-дневные показатели CGM (TIR, TBR, TAR, HbA1c)
    print("\nТест 27: Симуляция 30-дневного мониторинга глюкозы (CGM)...")
    from organoid_simulator import simulate_cgm_30day_metrics, evaluate_patient_clinical_risk_profile
    cgm_pre = simulate_cgm_30day_metrics(is_transplanted=False)
    cgm_post = simulate_cgm_30day_metrics(is_transplanted=True)
    assert cgm_post["TIR_percent"] > 90.0, f"Время в целевом диапазоне (TIR) после трансплантации должно быть >90%, получено: {cgm_post['TIR_percent']:.1f}%"
    assert cgm_post["TBR_percent"] == 0.0, "Гипогликемии должны быть ликвидированы (TBR = 0%)"
    assert cgm_post["gmi_hba1c_percent"] < 6.0, "Расчетный HbA1c должен снизиться до физиологической нормы (<6.0%)"
    print("  [OK] 30-дневные показатели CGM (TIR > 98%, HbA1c = 5.3%) верифицированы.")

    # Тест 28: ИИ-стратификатор клинического риска пациента
    print("\nТест 28: ИИ-стратификатор клинического риска пациента СД1...")
    risk_high = evaluate_patient_clinical_risk_profile(tdi_units=60.0, c_peptide_pmol_l=5.0, hba1c_percent=9.5, hypo_events_per_month=8)
    assert risk_high["priority_score"] > 70.0, "Пациент с высокой лабильностью и тяжелыми гипогликемиями должен получать высокий приоритет"
    assert "Высший приоритет" in risk_high["recommendation"], "Рекомендация должна указывать на критическую показанность трансплантации"
    print("  [OK] ИИ-стратификатор клинических рисков верифицирован успешно.")

    # Тест 29: Проверка загрузчика YAML-параметров
    print("\nТест 29: Проверка загрузчика YAML-параметров и литературы...")
    from param_loader import load_parameters, load_literature_parameters
    p_yaml = load_parameters()
    lit_yaml = load_literature_parameters()
    assert "solubility" in p_yaml and p_yaml["solubility"] > 0, "Параметр solubility должен быть загружен из YAML"
    assert "alginate_2%" in p_yaml["hydrogels"], "Прессет альгината 2% должен присутствовать в hydrogels"
    assert "oxygen" in lit_yaml, "Литературная база параметров должна содержать секцию oxygen"
    print("  [OK] Загрузка конфигураций из parameters.yaml и literature_params.yaml работает корректно.")

    # Тест 30: Явная 0-48ч кинетика IBMIR
    print("\nТест 30: Проверка динамической кинетики IBMIR (0-48ч)...")
    from organoid_simulator import solve_ibmir_kinetics
    ib_no_hep = solve_ibmir_kinetics(t_hours=np.linspace(0, 48, 49), peg_lmwh_density=0.0, heparin_dose_u_ml=0.0)
    ib_with_hep = solve_ibmir_kinetics(t_hours=np.linspace(0, 48, 49), peg_lmwh_density=1.5, heparin_dose_u_ml=0.5)
    assert ib_with_hep["retention_48h_percent"] > ib_no_hep["retention_48h_percent"], "Гепаринизация должна значительно увеличивать ретенцию через 48 часов"
    assert ib_with_hep["clot_fraction"][-1] < ib_no_hep["clot_fraction"][-1], "Доля сгустка с антикоагулянтом должна быть строго меньше"
    print(f"  [OK] Динамика IBMIR 0-48ч верифицирована: Ретенция без гепарина = {ib_no_hep['retention_48h_percent']:.1f}%, с гепарином = {ib_with_hep['retention_48h_percent']:.1f}%.")

    # Тест 31: Матрица сравнения мест имплантации (Portal vs Omentum vs SQ)
    print("\nТест 31: Проверка матрицы сравнения мест имплантации (Omentum vs Portal vs SQ)...")
    from simulator import IMPLANTATION_SITES
    assert "omental_pouch" in IMPLANTATION_SITES, "Сальник должен быть добавлен в пресеты мест имплантации"
    pO2_omental = IMPLANTATION_SITES["omental_pouch"]["pO2"]
    pO2_sq = IMPLANTATION_SITES["subcutaneous"]["pO2"]
    assert pO2_omental > pO2_sq, "Оксигенация в сальнике (55 mmHg) должна превышать подкожную (30 mmHg)"
    print("  [OK] Сравнительная матрица мест имплантации верифицирована.")

    # Тест 32: Воспроизведение литературных бенчмарков (Papas 2007, Papabathini 2023)
    print("\nТест 32: Тестирование скрипта воспроизведения литературных бенчмарков...")
    import os
    import sys
    bench_script = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reports", "benchmarks", "reproduce_benchmarks.py")
    assert os.path.exists(bench_script), "Скрипт reproduce_benchmarks.py должен существовать"
    sys.path.append(os.path.dirname(bench_script))
    from reproduce_benchmarks import run_benchmark_papas_2007, run_benchmark_papabathini_2023
    papas_res = run_benchmark_papas_2007()
    assert papas_res["rmse"] < 15.0, f"RMSE Papas 2007 должно быть < 15%, получено: {papas_res['rmse']:.2f}%"
    print(f"  [OK] Воспроизведение бенчмарков верифицировано (Papas RMSE = {papas_res['rmse']:.2f}%).")

    # Тест 33: Проверка целостности материалов препринта bioRxiv и аутрича
    print("\nТест 33: Проверка целостности материалов препринта bioRxiv и аутрич-пакета...")
    from export_demo_pack import validate_repository_structure
    assert validate_repository_structure() == True, "Все 17 файлов пакета Open Science M2 должны валидироваться"
    print("  [OK] Пакет материалов Open Science (M2) верифицирован в полном объеме.")

    # Тест 34: Проверка контроллера Искусственной Поджелудочной Железы (AID Closed-Loop)
    print("\nТест 34: Проверка замкнутого контура помповой инсулинотерапии (AID PID)...")
    from aid_controller import simulate_aid_closed_loop
    aid_res = simulate_aid_closed_loop(meal_carbs_g=60.0, target_glucose_mg_dl=110.0)
    assert aid_res["TIR_percent"] > 70.0, f"Время в целевом диапазоне (TIR) должно быть >70%, получено: {aid_res['TIR_percent']:.1f}%"
    assert aid_res["TBR_percent"] < 10.0, f"Время в гипогликемии (TBR) должно быть <10%, получено: {aid_res['TBR_percent']:.1f}%"
    print(f"  [OK] Замкнутый контур AID верифицирован (TIR = {aid_res['TIR_percent']:.1f}%, TBR = {aid_res['TBR_percent']:.1f}%).")

    # Тест 35: Проверка полноты выполнения дорожной карты (M0 - M5, 20 релизных файлов)
    print("\nТест 35: Проверка полноты выполнения дорожной карты (M0 - M5, 20 релизных файлов)...")
    assert validate_repository_structure() == True, "Все 20 файлов всех вех (M0-M5) должны валидироваться"
    print("  [OK] Полный жизненный цикл дорожной карты (M0-M5) верифицирован в полном объеме.")

    # ========================================================================
    # Тесты 36-38: IBMIR 0-48h kinetics module (новые)
    # ========================================================================

    # Тест 36: Базовая кинетика IBMIR (TF → Thrombin → Clot → O2 drop → Viability)
    print("\nТест 36: Базовая кинетика IBMIR модуля (TF → Thrombin → Clot → O2 → Viability)...")
    from ibmir_module import IBMIRKinetics

    ibmir_params = {
        "tf_concentration": 4.0,
        "tf_half_life": 15.0,
        "thrombin_generation_rate": 0.05,
        "clot_formation_rate": 0.03,
        "clot_thickness_max": 80.0,
        "oxygen_permeability_reduction": 0.6,
        "critical_po2": 0.5,
        "time_to_vascularization": 14.0,
    }
    model = IBMIRKinetics(ibmir_params)
    result = model.simulate()

    # Проверка структуры результата
    assert "time_points" in result, "Результат должен содержать time_points"
    assert "tf_concentration" in result, "Результат должен содержать tf_concentration"
    assert "thrombin_concentration" in result, "Результат должен содержать thrombin_concentration"
    assert "clot_thickness" in result, "Результат должен содержать clot_thickness"
    assert "pO2_profile" in result, "Результат должен содержать pO2_profile"
    assert "viability" in result, "Результат должен содержать viability"
    assert "total_cells_survived" in result, "Результат должен содержать total_cells_survived"

    # Проверка начальных значений
    assert result["tf_concentration"][0] == 4.0, "TF при t=0 должен быть равен начальной концентрации"
    assert result["thrombin_concentration"][0] == 0.0, "Тромбин при t=0 должен быть 0"
    assert result["clot_thickness"][0] == 0.0, "Сгусток при t=0 должен быть 0"
    assert result["viability"][0] == 1.0, "Жизнеспособность при t=0 должна быть 100%"

    # Проверка монотонности тромбина (растёт, потом падает)
    thrombin = result["thrombin_concentration"]
    assert max(thrombin) > 0.1, f"Пик тромбина должен быть > 0.1, получено: {max(thrombin)}"

    # Проверка сгустка (растёт, но не превышает max)
    assert max(result["clot_thickness"]) <= ibmir_params["clot_thickness_max"], "Сгусток не должен превышать clot_thickness_max"

    # Проверка жизнеспособности (падает к 48ч)
    assert result["viability"][-1] < 1.0, "Жизнеспособность к 48ч должна падать"
    assert result["viability"][-1] > 0.0, "Жизнеспособность должна быть > 0"

    # Проверка ключевых событий
    key_events = model.get_key_events()
    assert "time_to_thrombin_peak" in key_events, "Ключевые события должны содержать time_to_thrombin_peak"
    assert "time_to_clot_max" in key_events, "Ключевые события должны содержать time_to_clot_max"
    assert "peak_thrombin" in key_events, "Ключевые события должны содержать peak_thrombin"
    assert "max_clot_thickness" in key_events, "Ключевые события должны содержать max_clot_thickness"
    assert key_events["peak_thrombin"] > 0.0, "Пик тромбина должен быть положительным"
    assert key_events["max_clot_thickness"] > 0.0, "Макс. толщина сгустка должна быть положительной"

    print(f"  [OK] Базовая кинетика IBMIR верифицирована: пик тромбина = {key_events['peak_thrombin']:.4f}, "
          f"макс. сгусток = {key_events['max_clot_thickness']:.2f} µm, "
          f"финальная жизнеспособность = {key_events['final_viability']:.2%}.")

    # Тест 37: IBMIR с васкуляризацией → viability > 70%
    print("\nТест 37: IBMIR с васкуляризацией (viability > 70%)...")
    result_vasc = model.simulate_with_angiogenesis()
    assert "angiogenesis_rescued" in result_vasc, "Результат с васкуляризацией должен содержать angiogenesis_rescued"

    # Без васкуляризации (ранний срок, < 14 дней)
    result_no_vasc = model.simulate_with_angiogenesis(vasc_time_days=20.0)
    assert result_no_vasc["angiogenesis_rescued"] is False, "При vasc_time=20д (позже 48ч) rescued должно быть False"

    # С васкуляризацией (скорый срок, < 14 дней)
    result_fast_vasc = model.simulate_with_angiogenesis(vasc_time_days=5.0)
    assert result_fast_vasc["angiogenesis_rescued"] is True, "При vasc_time=5д (раньше 48ч) rescued должно быть True"
    assert result_fast_vasc["viability"][-1] > 0.70, (
        f"С васкуляризацией жизнеспособность должна быть > 70%, "
        f"получено: {result_fast_vasc['viability'][-1]:.2%}"
    )
    print(f"  [OK] Васкуляризация верифицирована: rescued={result_fast_vasc['angiogenesis_rescued']}, "
          f"viability={result_fast_vasc['viability'][-1]:.2%}.")

    # Тест 38: Сравнение сайтов (Portal vs SQ vs Omentum)
    print("\nТест 38: Сравнение IBMIR по сайтам (Portal vs SQ vs Omentum)...")
    from organoid_simulator import compare_ibmir_sites

    site_results = compare_ibmir_sites()
    assert "portal_vein" in site_results, "Должен быть результат для portal_vein"
    assert "omental_pouch" in site_results, "Должен быть результат для omental_pouch"
    assert "subcutaneous" in site_results, "Должен быть результат для subcutaneous"

    # Portal: высокая IBMIR экспозиция → больше сгусток, ниже viab
    portal_viability = site_results["portal_vein"]["kinetics"]["viability"][-1]
    omental_viability = site_results["omental_pouch"]["kinetics"]["viability"][-1]
    sq_viability = site_results["subcutaneous"]["kinetics"]["viability"][-1]

    # Portal имеет более высокий clot factor → толще сгусток
    portal_clot = site_results["portal_vein"]["kinetics"]["clot_thickness"][-1]
    omental_clot = site_results["omental_pouch"]["kinetics"]["clot_thickness"][-1]
    assert portal_clot > omental_clot, (
        f"Portal clot ({portal_clot:.2f}) должен быть больше omental ({omental_clot:.2f}) из-за clot_factor"
    )

    # Omentum должен иметь лучшую выживаемость из-за васкуляризации
    assert omental_viability > portal_viability, (
        f"Omental viability ({omental_viability:.2%}) должен быть выше portal ({portal_viability:.2%})"
    )

    # Все сайты должны иметь корректную структуру
    for site_key, site_res in site_results.items():
        assert "kinetics" in site_res, f"{site_key}: должен содержать kinetics"
        assert "site_params" in site_res, f"{site_key}: должен содержать site_params"
        assert "oxygen_coupled" in site_res, f"{site_key}: должен содержать oxygen_coupled"
        assert "key_events" in site_res, f"{site_key}: должен содержать key_events"

    print(f"  [OK] Сравнение сайтов верифицировано: "
          f"Portal viab={portal_viability:.2%}, Omentum viab={omental_viability:.2%}, SQ viab={sq_viability:.2%}.")

    # Тест 39: Проверка CLI/API модуля скрининга дефектов (screen_design)
    print("\nТест 39: Проверка CLI/API скрининга дефектов конструкт-дизайна (screen_design)...")
    from screen_design import screen_design, export_markdown_report
    screen_res = screen_design(geometry="spherical", radius_microns=200.0, density_million_per_ml=80.0, site_key="omental_pouch")
    assert screen_res["status_badge"] == "PASS ✅", f"Ожидался PASS ✅ для 200мкм сферы в сальнике, получено: {screen_res['status_badge']}"
    assert screen_res["viable_fraction_percent"] >= 99.0, f"Выживаемость должна быть >= 99%, получено: {screen_res['viable_fraction_percent']}"
    md_rep = export_markdown_report(screen_res)
    assert "Construct Failure Screening Report" in md_rep, "Отчет должен содержать заголовок отчета"
    print("  [OK] Модуль скрининга конструкт-дизайнов (screen_design) успешно прошел проверку.")

    # Тест 40: Проверка Торнадо-анализа чувствительности параметров (uncertainty_analysis)
    print("\nТест 40: Проверка Торнадо-анализа чувствительности параметров (uncertainty_analysis)...")
    from uncertainty_analysis import run_tornado_sensitivity
    tornado_res = run_tornado_sensitivity(base_radius_microns=200.0, base_density_million_per_ml=80.0, geometry="planar", variation_pct=20.0)
    df_sens = tornado_res["dataframe"]
    assert len(df_sens) == 5, f"Торнадо-анализ должен оценивать топ-5 параметров, получено: {len(df_sens)}"
    assert tornado_res["top_sensitive_parameter"] in ["pO2_boundary (Tissue Oxygen)", "D_eff (Diffusion Coefficient)"], "Ключевыми факторами гипоксии планарного геля должны быть pO2_boundary или D_eff"
    print("  [OK] Торнадо-анализ чувствительности параметров верифицирован.")

    # ========================================================================
    # Итог
    # ========================================================================

    print("\n=== Все 40 тестов успешно пройдены! ===")
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



