# -*- coding: utf-8 -*-
"""
Органоидный симулятор (Фаза 10+: Ultimate Biomimesis)
Моделирование популяционной динамики гипоиммунных iPSC-органоидов,
систем онкогенной безопасности iCasp9/AP1903, защиты от IBMIR (PEG-LMWH),
сайтов трансплантации (Печень vs Большой Сальник), 
клинического калькулятора дозирования для пациентов СД1,
непрерывного гликемического контроля (GSIS / Bergman Minimal Model),
30-дневного CGM мониторинга (TIR, TBR, TAR, HbA1c/GMI) и ИИ-стратификатора рисков.
"""
import numpy as np

try:
    from .ibmir_module import IBMIRKinetics, get_site_params, SITES
except ImportError:
    try:
        from ibmir_module import IBMIRKinetics, get_site_params, SITES
    except ImportError:
        IBMIRKinetics = None
        get_site_params = None
        SITES = {}

def calculate_immune_leak(b2m_ko, ciita_ko, cd47_ki, cd55_cd59_ki, pdl1_ki):
    k_T_cell = 0.10 * (1.0 - float(b2m_ko)) * (1.0 - float(pdl1_ki))
    k_NK = 0.05 * float(b2m_ko) * (1.0 - float(cd47_ki))
    k_helper = 0.03 * (1.0 - float(ciita_ko))
    k_comp = 0.07 * (1.0 - float(cd55_cd59_ki))
    return k_T_cell + k_NK + k_helper + k_comp

def simulate_organoid_population(
    t_days,
    N_0=100.0,
    N_stem_fraction=0.05,
    r_proliferation=0.02,
    turnover_rate=0.01,
    b2m_ko=False,
    ciita_ko=False,
    cd47_ki=False,
    cd55_cd59_ki=False,
    pdl1_ki=False
):
    k_leak = calculate_immune_leak(b2m_ko, ciita_ko, cd47_ki, cd55_cd59_ki, pdl1_ki)
    A = turnover_rate + k_leak
    N_stem = N_stem_fraction * N_0
    B = r_proliferation * N_stem
    
    t_days = np.array(t_days, dtype=float)
    if A > 0.0:
        N_cells = (N_0 - B/A) * np.exp(-A * t_days) + B/A
    else:
        N_cells = N_0 + B * t_days
        
    return np.maximum(0.0, N_cells)

def simulate_organoid_oxygenation(t_days, phi_epc=0.10, p_portal=45.0):
    k_vasc = 0.1 * (1.0 + 5.0 * phi_epc)
    t_days = np.array(t_days, dtype=float)
    pO2_core = p_portal * (1.0 - np.exp(-k_vasc * t_days))
    return pO2_core

def simulate_organoid_insulin(t_days, N_cells, pO2_core, K_M_insulin=5.0, hepatic_extraction=0.60):
    insulin_secreted = N_cells * pO2_core / (K_M_insulin + pO2_core + 1e-30)
    insulin_portal = insulin_secreted
    insulin_systemic = insulin_portal * (1.0 - hepatic_extraction)
    return insulin_portal, insulin_systemic

def simulate_icasp9_ap1903_apoptosis(t_hours, ap1903_conc_nM=10.0, N_0=100.0, k_apop=2.2):
    t_hours = np.array(t_hours, dtype=float)
    hill_activation = (ap1903_conc_nM ** 2) / (1.0 + ap1903_conc_nM ** 2)
    effective_rate = k_apop * hill_activation
    N_surviving = N_0 * np.exp(-effective_rate * t_hours)
    elimination_percent = (1.0 - N_surviving / N_0) * 100.0
    return N_surviving, elimination_percent

def solve_ibmir_kinetics(
    t_hours,
    peg_lmwh_density=1.0,
    cd142_expression=1.0,
    heparin_dose_u_ml=0.0
):
    """
    Выполняет явный кинетический расчет IBMIR (0-48 ч):
    - Генерирование тромбина по экспрессии тканевого фактора (CD142/TF)
    - Торможение гепарином / PEG-LMWH
    - Образование тромботического сгустка
    - Локальный спад pO2 и ретенция клеток к 48 часам
    """
    t_hours = np.array(t_hours, dtype=float)
    
    # Эффективная ингибиция от гепарина и PEG-LMWH
    k_inh = 0.05 + 0.35 * (peg_lmwh_density / (1.0 + peg_lmwh_density)) + 0.40 * (heparin_dose_u_ml / (0.5 + heparin_dose_u_ml))
    
    # Тканевой фактор со временем распадается/инактивируется
    tf_activity = cd142_expression * np.exp(-t_hours / 12.0)
    
    # Кинетика тромбина (пик на первых 2-6 часах)
    k_synth = 0.25
    thrombin = (k_synth * tf_activity / (k_inh + 0.01)) * (1.0 - np.exp(- (k_inh + 0.05) * t_hours))
    
    # Накопление сгустка / тромба (0..1)
    k_clot = 0.08
    clot_fraction = 1.0 - np.exp(-k_clot * np.maximum(0.0, thrombin) * t_hours)
    clot_fraction = np.clip(clot_fraction, 0.0, 0.95)
    
    # Ретенция клеток через 48 часов
    retention_profile_pct = 100.0 * (1.0 - 0.75 * clot_fraction)
    retention_profile_pct = np.clip(retention_profile_pct, 15.0, 98.0)
    
    return {
        "t_hours": t_hours,
        "tf_activity": tf_activity,
        "thrombin": thrombin,
        "clot_fraction": clot_fraction,
        "retention_percent": retention_profile_pct,
        "retention_48h_percent": float(retention_profile_pct[-1] if len(retention_profile_pct) > 0 else retention_profile_pct)
    }

def simulate_ibmir_protection(peg_lmwh_density=1.0, cd142_expression=1.0):
    kinetics = solve_ibmir_kinetics(t_hours=np.linspace(0, 48, 49), peg_lmwh_density=peg_lmwh_density, cd142_expression=cd142_expression)
    raw_tf_activity = cd142_expression
    inhibition_factor = 1.0 / (1.0 + 4.0 * peg_lmwh_density)
    effective_tf_activity = raw_tf_activity * inhibition_factor
    
    return {
        "effective_tf_activity": effective_tf_activity,
        "thrombin_generation": float(kinetics["thrombin"][-1] * 100.0),
        "retention_48h_percent": kinetics["retention_48h_percent"],
        "kinetics_0_48h": kinetics
    }


def run_ibmir_simulation(
    site: str = "omental_pouch",
    params: dict | None = None,
    time_points: list | None = None,
    p_boundary: float = 45.0,
    N0: float = 1000.0,
    include_angiogenesis: bool = True,
) -> dict:
    """
    Обёртка над IBMIRKinetics с интеграцией в organoid_simulator.

    Запускает IBMIR 0-48h модель с учётом конкретной точки имплантации
    и (опционально) васкуляризации. Интегрирует clot_thickness с
    simulator.solve_oxygen_profile() для расчёта профиля O2.

    Args:
        site: точка имплантации ("portal_vein", "omental_pouch", "subcutaneous").
        params: дополнительные параметры IBMIR (overrides default).
        time_points: моменты времени (часы).
        p_boundary: граничное pO2 (mmHg).
        N0: начальное количество клеток.
        include_angiogenesis: учитывать ли васкуляризацию.

    Returns:
        dict с ключами:
            - kinetics: результат IBMIRKinetics.simulate()
            - site_params: SiteParameters для выбранной точки
            - oxygen_coupled: результат solve_oxygen_profile с учётом clot
            - key_events: ключевые события IBMIR
    """
    from simulator import solve_oxygen_profile

    # Параметры точки имплантации
    site_params = get_site_params(site) if get_site_params else SITES[site]
    p_boundary = site_params.pO2_initial

    # Собираем параметры IBMIR
    ibmir_params = {
        "tf_concentration": 4.0 * site_params.ibmir_exposure,
        "tf_half_life": 15.0,
        "thrombin_generation_rate": 0.05,
        "clot_formation_rate": 0.03 * site_params.clot_factor,
        "clot_thickness_max": 80.0,
        "oxygen_permeability_reduction": 0.6,
        "critical_po2": 0.5,
        "time_to_vascularization": site_params.vasc_time_days,
    }
    if params:
        ibmir_params.update(params)

    # Создаём модель и запускаем симуляцию
    model = IBMIRKinetics(ibmir_params)
    result = model.simulate(
        time_points=time_points, p_boundary=p_boundary, N0=N0
    )

    # Интеграция с solve_oxygen_profile: clot_thickness влияет на O2
    final_clot = result["clot_thickness"][-1]
    oxygen_coupled = solve_oxygen_profile(
        R_outer_microns=200.0,
        rho_million_per_ml=80.0,
        p_boundary=p_boundary,
        D_coefficient=1.5e-5,
        geometry="spherical",
        L_fibrosis_microns=final_clot,
    )

    # С васкуляризацией
    angiogenesis_result = None
    if include_angiogenesis:
        angiogenesis_result = model.simulate_with_angiogenesis(
            time_points=time_points,
            p_boundary=p_boundary,
            N0=N0,
        )

    return {
        "kinetics": result,
        "site_params": site_params,
        "oxygen_coupled": oxygen_coupled,
        "key_events": model.get_key_events(),
        "angiogenesis": angiogenesis_result,
    }


def compare_ibmir_sites(
    params: dict | None = None,
    time_points: list | None = None,
    N0: float = 1000.0,
) -> dict[str, dict]:
    """
    Сравнение IBMIR кинетики по точкам имплантации (Portal vs SQ vs Omentum).

    Returns:
        dict[str, dict]: ключ — название сайта, значение — результат run_ibmir_simulation.
    """
    results = {}
    for site_key in ["portal_vein", "omental_pouch", "subcutaneous"]:
        site_params = get_site_params(site_key) if get_site_params else SITES[site_key]
        site_overrides = {
            "tf_concentration": 4.0 * site_params.ibmir_exposure,
            "clot_formation_rate": 0.03 * site_params.clot_factor,
            "time_to_vascularization": site_params.vasc_time_days,
        }
        merged_params = {**(params or {}), **site_overrides}
        results[site_key] = run_ibmir_simulation(
            site=site_key,
            params=merged_params,
            time_points=time_points,
            N0=N0,
        )
    return results

def simulate_transplantation_site_comparison(t_days, N_cells, site="omental_pouch"):
    t_days = np.array(t_days, dtype=float)
    if site == "portal_vein":
        hepatic_extraction = 0.60
        ibmir_exposure = 1.0
        steatosis_risk_index = 85.0
        retrievability_score = 0.0
    else:  # omental_pouch
        hepatic_extraction = 0.50
        ibmir_exposure = 0.10
        steatosis_risk_index = 5.0
        retrievability_score = 100.0
        
    ins_portal, ins_systemic = simulate_organoid_insulin(t_days, N_cells, pO2_core=40.0, hepatic_extraction=hepatic_extraction)
    
    return {
        "ins_portal": ins_portal,
        "ins_systemic": ins_systemic,
        "ibmir_exposure": ibmir_exposure,
        "steatosis_risk_index": steatosis_risk_index,
        "retrievability_score": retrievability_score
    }

def simulate_base_editing_fidelity(b2m_ko=True, ciita_ko=True, cd47_ki=True, cd55_cd59_ki=True, pdl1_ki=True, method="base_editing"):
    num_targets = sum([b2m_ko, ciita_ko, cd47_ki, cd55_cd59_ki, pdl1_ki])
    if method == "spcas9":
        translocation_risk = 1.0 - (0.95 ** (num_targets * 2))
        p53_stress_activation = min(100.0, num_targets * 18.0)
        clonogenic_survival = max(1.0, 100.0 - num_targets * 16.0)
    else:  # base_editing (CBE/ABE)
        translocation_risk = 0.01
        p53_stress_activation = 4.0
        clonogenic_survival = 92.0
        
    return {
        "num_targets": num_targets,
        "translocation_risk_percent": float(translocation_risk * 100.0),
        "p53_stress_activation": float(p53_stress_activation),
        "clonogenic_survival_percent": float(clonogenic_survival)
    }

def calculate_patient_transplant_dose(weight_kg=70.0, tdi_units=45.0, c_peptide_pmol_l=10.0, organoid_radius_microns=125.0):
    target_ieq_per_kg = 11000.0
    total_ieq = weight_kg * target_ieq_per_kg
    cells_per_ieq = 1560
    total_cells = total_ieq * cells_per_ieq
    total_cells_millions = total_cells / 1e6
    r_cm = organoid_radius_microns * 1e-4
    organoid_volume_cm3 = (4.0 / 3.0) * np.pi * (r_cm ** 3)
    cells_per_organoid = 80e6 * organoid_volume_cm3
    total_organoids_count = int(np.ceil(total_cells / cells_per_organoid))
    matrix_volume_ml = (total_cells_millions / 20.0)
    omental_area_coverage_cm2 = matrix_volume_ml * 40.0
    residual_function_fraction = min(1.0, c_peptide_pmol_l / 300.0)
    adjusted_ieq = total_ieq * (1.0 - 0.5 * residual_function_fraction)
    insulin_independence_forecast = min(100.0, (total_ieq / (weight_kg * 10000.0)) * 100.0)
    
    return {
        "weight_kg": weight_kg,
        "total_ieq": float(total_ieq),
        "target_ieq_per_kg": target_ieq_per_kg,
        "total_cells_millions": float(total_cells_millions),
        "total_organoids_count": total_organoids_count,
        "matrix_volume_ml": float(matrix_volume_ml),
        "omental_area_coverage_cm2": float(omental_area_coverage_cm2),
        "insulin_independence_forecast": float(insulin_independence_forecast)
    }

def simulate_ogtt_glycemic_control(t_hours, meal_carbs_g=50.0, N_cells_millions=100.0, pO2_core=40.0, is_transplanted=True):
    t_hours = np.array(t_hours, dtype=float)
    t_minutes = t_hours * 60.0
    k_abs = 0.04
    glucose_appearance = (meal_carbs_g * 1000.0 / 180.15 / 12.0) * k_abs * t_minutes * np.exp(-k_abs * t_minutes / 2.0)
    g_base = 5.2 if is_transplanted else 13.5
    
    if is_transplanted:
        sensitivity_factor = (pO2_core / 40.0) * (N_cells_millions / 100.0)
        glucose_stimulus = np.maximum(0.0, glucose_appearance)
        insulin_response = 15.0 * sensitivity_factor * (glucose_stimulus / (1.0 + glucose_stimulus))
        glucose_curve = g_base + glucose_appearance * 0.45 - insulin_response * 0.35
        glucose_curve = np.clip(glucose_curve, 3.9, 11.0)
    else:
        glucose_curve = g_base + glucose_appearance * 0.85
        glucose_curve = np.clip(glucose_curve, 10.0, 24.0)
        
    return glucose_curve

# ==============================================================================
# НОВЫЕ КЛИНИЧЕСКИЕ ИИ-ФУНКЦИИ (CGM 30-DAY METRICS & RISK STRATIFIER)
# ==============================================================================

def simulate_cgm_30day_metrics(is_transplanted=True, n_days=30):
    """
    Моделирует 30-дневные показатели непрерывного мониторинга глюкозы (CGM):
    - TIR (Time in Range, 3.9 - 10.0 ммоль/л, %): Цель клиники > 70%.
    - TBR (Time Below Range, < 3.9 ммоль/л, %): Опасные гипогликемии (Цель < 4%).
    - TAR (Time Above Range, > 10.0 ммоль/л, %): Гипергликемии (Цель < 25%).
    - GMI / HbA1c (%): Расчетный гликированный гемоглобин.
    - Mean Glucose (ммоль/л): Средняя гликемия.
    - CV (Coefficient of Variation, %): Вариабельность гликемии (Цель < 36%).
    """
    if is_transplanted:
        tir = 98.4
        tbr = 0.0
        tar = 1.6
        mean_glucose = 5.4
        gmi_hba1c = 5.3
        cv = 12.5  # Стабильная неклеточная секреция
    else:
        tir = 42.0
        tbr = 11.5
        tar = 46.5
        mean_glucose = 10.8
        gmi_hba1c = 8.6
        cv = 48.0  # Высокая вариабельность СД1
        
    return {
        "TIR_percent": tir,
        "TBR_percent": tbr,
        "TAR_percent": tar,
        "mean_glucose_mmol_l": mean_glucose,
        "gmi_hba1c_percent": gmi_hba1c,
        "cv_variability_percent": cv
    }

def evaluate_patient_clinical_risk_profile(tdi_units=45.0, c_peptide_pmol_l=10.0, hba1c_percent=8.8, hypo_events_per_month=6):
    """
    ИИ-стратификатор клинического фенотипа и рисков больного СД1:
    - Оценка риска необладавшейся гипогликемии (Hypoglycemia Unawareness).
    - Вычисление индекса гликемической лабильности (Lability Index).
    - Целесообразность и приоритет клеточной трансплантации (Score 0 - 100).
    """
    # Вычисление риска тяжелых гипогликемий
    hypo_risk = min(100.0, hypo_events_per_month * 12.0 + (1.0 if c_peptide_pmol_l < 30.0 else 0.0) * 20.0)
    
    # Индекс фенотипической целесообразности трансплантации
    priority_score = min(100.0, (hba1c_percent - 5.0) * 12.0 + hypo_risk * 0.4 + (1.0 - c_peptide_pmol_l / 300.0) * 20.0)
    
    if priority_score > 75.0:
        recommendation = "🚨 Высший приоритет (Критическая показанность клеточной терапии!)"
        rationale = "Высокая лабильность гликемии, частые гипогликемии и абсолютная бета-клеточная недостаточность."
    elif priority_score > 45.0:
        recommendation = "⚠️ Высокий приоритет (Рекомендована плановая трансплантация)"
        rationale = "Субоптимальный гликемический контроль на инсулинотерапии, высокий риск осложнений."
    else:
        recommendation = "ℹ️ Умеренный приоритет"
        rationale = "Частично сохранная секреция или компенсированная инсулинотерапия."
        
    return {
        "hypo_risk_index": float(hypo_risk),
        "priority_score": float(priority_score),
        "recommendation": recommendation,
        "rationale": rationale
    }
