import numpy as np
from scipy.integrate import solve_bvp

try:
    from .param_loader import load_parameters
except ImportError:
    try:
        from param_loader import load_parameters
    except ImportError:
        def load_parameters():
            return {
                "solubility": 1.34e-9,
                "v_max": 1.2e-16,
                "k_m": 0.5,
                "k_m_insulin": 5.0,
                "hydrogels": {
                    "water": {"name": "Чистая вода (контроль)", "D": 3.0e-5},
                    "alginate_1%": {"name": "1% Альгинат натрия", "D": 2.1e-5},
                    "alginate_2%": {"name": "2% Альгинат натрия (стандарт)", "D": 1.5e-5}
                },
                "implantation_sites": {
                    "arterial": {"name": "Артериальное русло (прямая перфузия)", "pO2": 95.0, "description": "Идеальные условия, высокая оксигенация (например, сосудистый порт)."},
                    "venous": {"name": "Венозная сеть (или воротная вена)", "pO2": 40.0, "description": "Умеренная оксигенация, характерная для внутренних органов."},
                    "subcutaneous": {"name": "Подкожная клетчатка (типичное место введения)", "pO2": 30.0, "description": "Низкое давление кислорода. Самое доступное, но сложное место для приживления."},
                    "extreme_hypoxia": {"name": "Зона фиброза / Выраженная гипоксия", "pO2": 10.0, "description": "Худший сценарий. Капсула обросла фиброзной тканью, доступ кислорода сильно ограничен."}
                }
            }

_PARAMS = load_parameters()

# Константы модели (физические параметры)
SOLUBILITY = _PARAMS.get("solubility", 1.34e-9)
V_MAX = _PARAMS.get("v_max", 1.2e-16)
K_M = _PARAMS.get("k_m", 0.5)
K_M_INSULIN = _PARAMS.get("k_m_insulin", 5.0)

# Пресеты гидрогелей и их коэффициенты диффузии кислорода (см^2 / с)
HYDROGELS = _PARAMS.get("hydrogels", {
    "water": {"name": "Чистая вода (контроль)", "D": 3.0e-5},
    "alginate_1%": {"name": "1% Альгинат натрия", "D": 2.1e-5},
    "alginate_2%": {"name": "2% Альгинат натрия (стандарт)", "D": 1.5e-5}
})

# Места имплантации и их парциальное давление кислорода (mmHg)
IMPLANTATION_SITES = _PARAMS.get("implantation_sites", {
    "arterial": {"name": "Артериальное русло (прямая перфузия)", "pO2": 95.0, "description": "Идеальные условия, высокая оксигенация."},
    "venous": {"name": "Венозная сеть (или воротная вена)", "pO2": 40.0, "description": "Умеренная оксигенация, характерная для внутренних органов."},
    "subcutaneous": {"name": "Подкожная клетчатка (типичное место введения)", "pO2": 30.0, "description": "Низкое давление кислорода."},
    "extreme_hypoxia": {"name": "Зона фиброза / Выраженная гипоксия", "pO2": 10.0, "description": "Худший сценарий."}
})

def solve_oxygen_profile(
    R_outer_microns, rho_million_per_ml, p_boundary, D_coefficient, 
    V_max_cell=V_MAX, geometry="planar", L_fibrosis_microns=0.0, D_fibrosis=1.0e-5,
    rho_mac_million_per_ml=0.0, q_ogm_mmHg_per_sec=0.0,
    catalase_activity_relative=1.0, catalase_half_life_days=1.5,
    buffer_capacity_mM=10.0, swelling_ratio=1.0,
    plga_acidification_factor=0.0, t_days=0.0,
    tethered_catalase=False, E_0=50.0, species="Mouse",
    phi_pfc=0.0, av_loop_flow=False, crispr_hypoimmune=False,
    cd47_overexpression=False,
    tau_blood=5.0,
    anticoagulation=False,
    pO2_pfc_saturation=200.0,
    turnover_rate=0.0,
    complement_protection=False
):
    """
    Решает одномерную краевую задачу диффузии-потребления кислорода в гидрогелевой капсуле.
    Поддерживает три типа геометрии: planar, cylindrical, spherical.
    При указании rho_mac_million_per_ml > 0 или q_ogm_mmHg_per_sec > 0 переключается на
    двухзонную модель с активными макрофагами и/или OGM генератором.
    Учитывает pH-защелачивание/закисление, токсичность H2O2, распад каталазы и набухание геля (MWCO).
    """
    # Масштабирование размеров и диффузии из-за набухания геля

    if av_loop_flow:
        k_thrombo_0 = 0.05
        k_thrombo = k_thrombo_0 * (1.0 + (max(0.0, 1.5 - tau_blood) / 0.5)**2 + (max(0.0, tau_blood - 8.0) / 2.0)**2)
        if anticoagulation:
            k_thrombo_eff = k_thrombo * 0.1
        else:
            k_thrombo_eff = k_thrombo
        k_hyperplasia = 0.005
        k_occlusion = k_thrombo_eff + k_hyperplasia
        p_boundary = 30.0 + (95.0 - 30.0) * np.exp(-k_occlusion * t_days)
        L_fibrosis_microns = 0.0

    SOLUBILITY_eff = SOLUBILITY * (1.0 - phi_pfc) + (20.0 * SOLUBILITY) * phi_pfc
    p_pfc_t = pO2_pfc_saturation * np.exp(-t_days / 2.0) if phi_pfc > 0.0 else 0.0

    R_outer_microns_eff = R_outer_microns * (swelling_ratio ** (1.0 / 3.0))
    L_fibrosis_microns_eff = L_fibrosis_microns * (swelling_ratio ** (1.0 / 3.0))
    D_coefficient_eff = D_coefficient * (swelling_ratio ** (2.0 / 3.0))

    # Если активированы усовершенствованные параметры, решаем двухзонную модель
    if rho_mac_million_per_ml > 0.0 or q_ogm_mmHg_per_sec > 0.0 or plga_acidification_factor > 0.0:
        R_outer = R_outer_microns_eff * 1e-4
        L_fib = L_fibrosis_microns_eff * 1e-4
        R_total = R_outer + L_fib
        
        rho_cells = rho_million_per_ml * 1e6
        
        # Закисление от PLGA наночастиц привлекает макрофаги (увеличивает их плотность)
        rho_mac_boosted = rho_mac_million_per_ml * (1.0 + 2.0 * plga_acidification_factor) * (1.0 + 1.5 * float(tethered_catalase))
        rho_macs = rho_mac_boosted * 1e6
        
        V_MAX_MAC = 3.0e-16
        K_M_MAC = 1.0
        w_trans = 2.0 * 1e-4
        
        # Моделирование каталазы и её распада
        if tethered_catalase:
            catalase_half_life_days_eff = 100.0
            catalase_activity_relative_eff = catalase_activity_relative * 0.25
        else:
            catalase_half_life_days_eff = catalase_half_life_days / (swelling_ratio ** 2)
            catalase_activity_relative_eff = catalase_activity_relative

        k_cat_0 = 0.05
        k_cat = k_cat_0 * catalase_activity_relative_eff * np.exp(-t_days * np.log(2.0) / catalase_half_life_days_eff)
        
        # Эффективная скорость генерации кислорода зависит от наличия каталазы для расщепления H2O2
        q_ogm_eff = q_ogm_mmHg_per_sec * (k_cat / (k_cat_0 + 1e-5))
        
        # Расчет накопления ROS (H2O2) в ядре
        P_H2O2 = 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY_eff * 1e9  # uM/s
        t_sec = t_days * 86400.0
        if k_cat > 1e-8:
            C_H2O2_core = (P_H2O2 / k_cat) * (1.0 - np.exp(-k_cat * t_sec))
        else:
            C_H2O2_core = P_H2O2 * t_sec
        
        # Расчет pH с учетом буферной емкости геля и закисления от PLGA
        P_OH = 4.0 * q_ogm_mmHg_per_sec * SOLUBILITY_eff * 1e9  # uM/s
        P_acid = plga_acidification_factor * 0.2 * SOLUBILITY_eff * 1e9  # uM/s
        
        OH_prod = P_OH * t_sec
        H_prod = P_acid * t_sec
        
        buffer_uM = buffer_capacity_mM * 1000.0
        delta_charge = OH_prod - H_prod
        
        pH_core = 7.4 + 3.0 * np.tanh(delta_charge / (buffer_uM + 1.0))

        # Расчет накопления кальция Ca2+ в ядре (ммоль/л)
        Ca_accum_mM = 1.2 + 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY_eff * 1e9 * t_sec * 1e-3
        f_Ca = np.exp(- (np.maximum(0.0, Ca_accum_mM - 2.0) ** 2) / (2.0 * (1.5 ** 2)))
        insulin_leak = np.minimum(1.0, np.maximum(0.0, Ca_accum_mM - 1.2) / 4.0)

        # Расчет деградации модуля Юнга и механического разрыва
        young_modulus_eff = E_0 * (swelling_ratio ** -2.0) * np.exp(-0.01 * t_days)
        sigma_ext = 2.0
        epsilon_max = 0.1
        rupture_risk = min(100.0, (sigma_ext / (young_modulus_eff * epsilon_max)) * 100.0)
        
        if rupture_risk > 50.0:
            f_rupture = 1.0 - rupture_risk / 100.0
        else:
            f_rupture = 1.0
        
        if geometry == "cylindrical":
            g_factor = 1.0
        elif geometry == "spherical":
            g_factor = 2.0
        else:
            g_factor = 0.0
            
        def equations(r, y):
            p = y[0]
            dp = y[1]
            p_pos = np.maximum(0.0, p)
            
            sigmoid = 1.0 / (1.0 + np.exp((r - R_outer) / w_trans))
            sigmoid_prime = -(1.0 / w_trans) * sigmoid * (1.0 - sigmoid)
            
            D_r = D_coefficient_eff * sigmoid + D_fibrosis * (1.0 - sigmoid)
            D_prime = (D_coefficient_eff - D_fibrosis) * sigmoid_prime
            
            # Потребление β-клетками
            R_cells_phys = (rho_cells * V_max_cell / SOLUBILITY_eff) * (p_pos / (K_M + p_pos)) * sigmoid
            # Потребление макрофагами в зоне фиброза
            R_macs_phys = (rho_macs * V_MAX_MAC / SOLUBILITY_eff) * (p_pos / (K_M_MAC + p_pos)) * (1.0 - sigmoid)
            # Генерация кислорода OGM в ядре (с учетом деградации каталазы)
            S_ogm_phys = q_ogm_eff * sigmoid
            
            term_source_sink = (R_cells_phys + R_macs_phys - S_ogm_phys) / D_r
            d2p = term_source_sink - (g_factor / (r + 1e-9) + D_prime / D_r) * dp
            return np.vstack((dp, d2p))
            
        def boundary_conditions(ya, yb):
            return np.array([ya[1] - 0.0, yb[0] - p_boundary])
            
        r_mesh = np.linspace(0.0, R_total, 300)
        y_init = np.zeros((2, r_mesh.size))
        y_init[0, :] = p_boundary
        
        sol = solve_bvp(equations, boundary_conditions, r_mesh, y_init, tol=1e-5, max_nodes=2000)
        
        z_coords_microns = np.linspace(0.0, R_total * 1e4, 1000)
        pO2_profile = np.interp(z_coords_microns, sol.x * 1e4, sol.y[0])
        pO2_profile = np.maximum(0.0, pO2_profile) + p_pfc_t
        
        # Интегральные характеристики считаем только по активному клеточному ядру (z <= R_outer)
        cell_mask = z_coords_microns <= (R_outer * 1e4)
        core_z = z_coords_microns[cell_mask]
        core_p = pO2_profile[cell_mask]
        
        if geometry == "planar":
            weights = np.ones_like(core_z)
        elif geometry == "cylindrical":
            weights = core_z
        else:
            weights = core_z ** 2
            
        weights_sum = np.sum(weights)
        
        # Расчет пространственного распределения H2O2 и pH (затухание к границам)
        H2O2_profile = C_H2O2_core * (1.0 - (z_coords_microns / (R_total * 1e4)) ** 2)
        H2O2_profile = np.maximum(0.0, H2O2_profile)
        
        pH_profile = 7.4 + (pH_core - 7.4) * (1.0 - (z_coords_microns / (R_total * 1e4)) ** 2)
        
        # Мультипликаторы жизнеспособности клеток
        f_pH = np.exp(- (pH_profile - 7.4) ** 2 / (2.0 * 0.3 ** 2))
        f_ROS = 1.0 / (1.0 + (H2O2_profile / 10.0) ** 2)
        
        MWCO = 30.0 * (swelling_ratio ** 2)
        f_IgG = np.ones_like(z_coords_microns) if crispr_hypoimmune else 1.0 / (1.0 + np.exp((MWCO - 150.0) / 10.0))
        
        if crispr_hypoimmune:
            f_NK = 1.0 if cd47_overexpression else np.exp(-0.1 * t_days)
            f_complement = 1.0 if complement_protection else np.exp(-0.3 * t_days)
        else:
            f_NK = 1.0
            f_complement = 1.0
            
        if phi_pfc > 0.0 and pO2_pfc_saturation > 200.0:
            f_hyperoxia = np.exp(- (pO2_pfc_saturation - 200.0)**2 / (2.0 * 150.0**2))
        else:
            f_hyperoxia = 1.0
            
        f_depletion = np.minimum(1.0, np.exp((0.002 - turnover_rate) * t_days))
        
        f_total = f_pH * f_ROS * f_IgG * f_Ca * f_rupture * f_NK * f_complement * f_hyperoxia * f_depletion
        
        viable_fraction_profile = (core_p >= K_M).astype(float) * f_total[cell_mask]
        viable_fraction = (np.sum(viable_fraction_profile * weights) / weights_sum) * 100.0
        
        # Инсулиносекреция
        insulin_capacity_profile = core_p / (K_M_INSULIN + core_p)
        max_possible_capacity = 95.0 / (K_M_INSULIN + 95.0)
        mean_insulin_capacity = np.sum(insulin_capacity_profile * weights) / weights_sum
        insulin_capacity = (mean_insulin_capacity / max_possible_capacity) * 100.0
        insulin_capacity = min(insulin_capacity, viable_fraction)
        
        return {
            "z": z_coords_microns,
            "pO2": pO2_profile,
            "pH": pH_profile,
            "H2O2": H2O2_profile,
            "viability_multiplier": f_total,
            "viable_fraction": viable_fraction,
            "insulin_capacity": insulin_capacity,
            "min_pO2": np.min(core_p),
            "Ca_accum_mM": Ca_accum_mM,
            "insulin_leak": insulin_leak,
            "rupture_risk": rupture_risk,
            "young_modulus_eff": young_modulus_eff
        }

    # Стандартный пассивный BVP-решатель
    R_outer_cm = R_outer_microns_eff * 1e-4
    L_fibrosis_cm = L_fibrosis_microns_eff * 1e-4
    rho = rho_million_per_ml * 1e6
    C_boundary = p_boundary * SOLUBILITY_eff
    K_m_conc = K_M * SOLUBILITY_eff
    
    if geometry == "cylindrical":
        g_factor = 1.0
    elif geometry == "spherical":
        g_factor = 2.0
    else:
        g_factor = 0.0
        
    def equations_std(z, y):
        C = y[0]
        dC = y[1]
        C_pos = np.maximum(0.0, C)
        reaction = (rho * V_max_cell / D_coefficient_eff) * (C_pos / (K_m_conc + C_pos))
        if g_factor == 0.0:
            d2C = reaction
        else:
            d2C = reaction - (g_factor * dC / (z + 1e-9))
        return np.vstack((dC, d2C))
        
    def boundary_conditions_std(ya, yb):
        if L_fibrosis_cm > 0.0:
            flux_inward = (D_fibrosis / (D_coefficient_eff * L_fibrosis_cm)) * (C_boundary - yb[0])
            return np.array([ya[1] - 0.0, yb[1] - flux_inward])
        else:
            return np.array([ya[1] - 0.0, yb[0] - C_boundary])
        
    z_mesh = np.linspace(0.0, R_outer_cm, 100)
    y_init = np.zeros((2, z_mesh.size))
    y_init[0, :] = C_boundary
    
    sol = solve_bvp(equations_std, boundary_conditions_std, z_mesh, y_init, tol=1e-5, max_nodes=1000)
    
    z_coords_microns = sol.x * 1e4
    C_profile = sol.y[0]
    pO2_profile = C_profile / SOLUBILITY_eff
    pO2_profile = np.maximum(0.0, pO2_profile) + p_pfc_t
    
    z_uniform = np.linspace(0.0, R_outer_microns_eff, 1000)
    pO2_uniform = np.interp(z_uniform, z_coords_microns, pO2_profile)
    
    if geometry == "planar":
        weights = np.ones_like(z_uniform)
    elif geometry == "cylindrical":
        weights = z_uniform
    else:
        weights = z_uniform ** 2
        
    weights_sum = np.sum(weights)
    
    t_sec = t_days * 86400.0
    Ca_accum_mM = 1.2 + 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY_eff * 1e9 * t_sec * 1e-3
    f_Ca = np.exp(- (np.maximum(0.0, Ca_accum_mM - 2.0) ** 2) / (2.0 * (1.5 ** 2)))
    insulin_leak = np.minimum(1.0, np.maximum(0.0, Ca_accum_mM - 1.2) / 4.0)
    
    young_modulus_eff = E_0 * (swelling_ratio ** -2.0) * np.exp(-0.01 * t_days)
    sigma_ext = 2.0
    epsilon_max = 0.1
    rupture_risk = min(100.0, (sigma_ext / (young_modulus_eff * epsilon_max)) * 100.0)
    
    if rupture_risk > 50.0:
        f_rupture = 1.0 - rupture_risk / 100.0
    else:
        f_rupture = 1.0

    # Набухание и прорыв MWCO в пассивной модели
    MWCO = 30.0 * (swelling_ratio ** 2)
    f_IgG = np.ones_like(z_uniform) if crispr_hypoimmune else 1.0 / (1.0 + np.exp((MWCO - 150.0) / 10.0))
    
    if crispr_hypoimmune:
        f_NK = 1.0 if cd47_overexpression else np.exp(-0.1 * t_days)
        f_complement = 1.0 if complement_protection else np.exp(-0.3 * t_days)
    else:
        f_NK = 1.0
        f_complement = 1.0
        
    if phi_pfc > 0.0 and pO2_pfc_saturation > 200.0:
        f_hyperoxia = np.exp(- (pO2_pfc_saturation - 200.0)**2 / (2.0 * 150.0**2))
    else:
        f_hyperoxia = 1.0
        
    f_depletion = np.minimum(1.0, np.exp((0.002 - turnover_rate) * t_days))
    
    f_total = f_IgG * np.ones_like(z_uniform) * f_Ca * f_rupture * f_NK * f_complement * f_hyperoxia * f_depletion
    
    viable_mask_uniform = pO2_uniform >= K_M
    viable_fraction_profile = viable_mask_uniform.astype(float) * f_total
    viable_fraction = (np.sum(viable_fraction_profile * weights) / weights_sum) * 100.0
    
    insulin_capacity_profile = pO2_uniform / (K_M_INSULIN + pO2_uniform)
    max_possible_capacity = 95.0 / (K_M_INSULIN + 95.0)
    mean_insulin_capacity = np.sum(insulin_capacity_profile * weights) / weights_sum
    insulin_capacity = (mean_insulin_capacity / max_possible_capacity) * 100.0
    insulin_capacity = min(insulin_capacity, viable_fraction)
    
    min_pO2 = np.min(pO2_profile)
    
    return {
        "z": z_uniform,
        "pO2": pO2_uniform,
        "pH": np.ones_like(z_uniform) * 7.4,
        "H2O2": np.zeros_like(z_uniform),
        "viability_multiplier": f_total,
        "viable_fraction": viable_fraction,
        "insulin_capacity": insulin_capacity,
        "min_pO2": min_pO2,
        "Ca_accum_mM": Ca_accum_mM,
        "insulin_leak": insulin_leak,
        "rupture_risk": rupture_risk,
        "young_modulus_eff": young_modulus_eff
    }

def run_density_sweep(R_outer_microns, p_boundary, D_coefficient, min_density=5, max_density=200, step=10, geometry="planar", L_fibrosis_microns=0.0, D_fibrosis=1.0e-5):
    """
    Выполняет симуляции для диапазона плотностей заселения клеток для выбранной геометрии.
    """
    densities = list(range(min_density, max_density + 1, step))
    viabilities = []
    insulin_capacities = []
    
    for rho in densities:
        res = solve_oxygen_profile(
            R_outer_microns, rho, p_boundary, D_coefficient, 
            geometry=geometry, 
            L_fibrosis_microns=L_fibrosis_microns, 
            D_fibrosis=D_fibrosis
        )
        viabilities.append(res["viable_fraction"])
        insulin_capacities.append(res["insulin_capacity"])
        
    return {
        "densities": np.array(densities),
        "viabilities": np.array(viabilities),
        "insulin_capacities": np.array(insulin_capacities)
    }

def solve_coupled_neovascularization(
    R_outer_microns,
    rho_million_per_ml,
    D_oxygen_coefficient,
    geometry="planar",
    L_fibrosis_microns=0.0,
    D_fibrosis=1.0e-5,
    V_loaded_relative=1.0,
    k_clear_tissue=15.0,  # 1/day
    beta_angiogenesis=0.15,  # 1/day
    K_vegf=0.1,
    p_base=30.0,
    p_max=60.0,
    days=21,
    P_loaded_relative=1.0,
    pdgf_burst_fraction=0.3,
    plga_acidification_factor=0.0,
    species="Mouse",
    av_loop_flow=False,
    tau_blood=5.0,
    anticoagulation=False
):
    """
    Simulates VEGF and PDGF diffusion over time (days) and calculates the coupled neovascularization
    response which raises oxygen boundary pressure and rescues cells from hypoxia.
    """
    # Масштабирование ангиогенеза на человека

    if av_loop_flow:
        p_base = 95.0
        p_max = 95.0
        L_fibrosis_microns = 0.0

    if species.lower() == "human":
        beta_angiogenesis_eff = beta_angiogenesis / 2.5
        if days == 21:
            days = 60
    else:
        beta_angiogenesis_eff = beta_angiogenesis

    # 1. Сетка для диффузии VEGF и PDGF
    R_outer_cm = R_outer_microns * 1e-4
    L_fibrosis_cm = L_fibrosis_microns * 1e-4
    L_tissue_cm = 500.0 * 1e-4
    
    R_max_cm = R_outer_cm + L_fibrosis_cm + L_tissue_cm
    
    Nr = 100
    r_grid = np.linspace(0.0, R_max_cm, Nr + 1)
    dr = R_max_cm / Nr
    
    # Геометрический фактор кривизны (g)
    if geometry == "cylindrical":
        g = 1.0
    elif geometry == "spherical":
        g = 2.0
    else:
        g = 0.0
        
    # Коэффициенты диффузии VEGF в зонах (см²/с)
    D_V_cap = 1.0e-11
    D_V_fib = D_fibrosis * 0.01 if L_fibrosis_microns > 0 else D_V_cap
    D_V_tis = 1.0e-10
    
    # Коэффициенты диффузии PDGF (крупный белок, диффузия медленнее в 2 раза)
    D_P_cap = 0.5 * D_V_cap
    D_P_fib = 0.5 * D_V_fib
    D_P_tis = 0.5 * D_V_tis
    
    D_grid = np.zeros(Nr + 1)
    D_P_grid = np.zeros(Nr + 1)
    k_grid = np.zeros(Nr + 1)
    
    # Переводим константы в единицы секунд
    sec_per_day = 86400.0
    k_clear_sec = k_clear_tissue / sec_per_day
    
    for i in range(Nr + 1):
        r = r_grid[i]
        if r <= R_outer_cm:
            D_grid[i] = D_V_cap
            D_P_grid[i] = D_P_cap
            k_grid[i] = 0.0
        elif r <= R_outer_cm + L_fibrosis_cm:
            D_grid[i] = D_V_fib
            D_P_grid[i] = D_P_fib
            k_grid[i] = 0.0
        else:
            D_grid[i] = D_V_tis
            D_P_grid[i] = D_P_tis
            k_grid[i] = k_clear_sec
            
    # Начальные условия
    C_V = np.zeros(Nr + 1)
    C_V[r_grid <= R_outer_cm] = V_loaded_relative
    
    C_P = np.zeros(Nr + 1)
    C_P[r_grid <= R_outer_cm] = pdgf_burst_fraction * P_loaded_relative
    
    # Симуляция по времени
    dt_days = 0.1
    sec_per_day = 86400.0
    dt_sec = dt_days * sec_per_day
    total_steps = int(days / dt_days)
    
    # История изменения параметров
    t_history = []
    C_interface_history = []
    C_P_interface_history = []
    p_boundary_history = []
    
    # Индекс границы раздела сред
    interface_idx = np.argmin(np.abs(r_grid - R_outer_cm))
    
    # Подготовка тридиагональной матрицы для неявного Эйлера (VEGF)
    a = np.zeros(Nr - 1)
    b = np.zeros(Nr)
    c = np.zeros(Nr - 1)
    
    # Подготовка тридиагональной матрицы для неявного Эйлера (PDGF)
    a_P = np.zeros(Nr - 1)
    b_P = np.zeros(Nr)
    c_P = np.zeros(Nr - 1)
    
    # Заполнение коэффициентов
    for i in range(Nr):
        r_i = r_grid[i]
        k_i = k_grid[i]
        
        # VEGF
        if i == 0:
            D_plus = 0.5 * (D_grid[0] + D_grid[1])
            B_i = - (g + 1.0) * 2.0 * D_plus / (dr ** 2)
            D_prime_i = (g + 1.0) * 2.0 * D_plus / (dr ** 2)
            b[0] = 1.0 + dt_sec * k_i - dt_sec * B_i
            c[0] = - dt_sec * D_prime_i
        else:
            D_plus = 0.5 * (D_grid[i] + D_grid[i+1])
            D_minus = 0.5 * (D_grid[i-1] + D_grid[i])
            r_plus = r_i + 0.5 * dr
            r_minus = r_i - 0.5 * dr
            A_i = (D_minus * (r_minus ** g)) / ((r_i ** g) * (dr ** 2))
            D_prime_i = (D_plus * (r_plus ** g)) / ((r_i ** g) * (dr ** 2))
            B_i = - (A_i + D_prime_i)
            b[i] = 1.0 + dt_sec * k_i - dt_sec * B_i
            if i < Nr - 1:
                c[i] = - dt_sec * D_prime_i
            a[i-1] = - dt_sec * A_i
            
        # PDGF
        if i == 0:
            D_plus_P = 0.5 * (D_P_grid[0] + D_P_grid[1])
            B_i_P = - (g + 1.0) * 2.0 * D_plus_P / (dr ** 2)
            D_prime_i_P = (g + 1.0) * 2.0 * D_plus_P / (dr ** 2)
            b_P[0] = 1.0 + dt_sec * k_i - dt_sec * B_i_P
            c_P[0] = - dt_sec * D_prime_i_P
        else:
            D_plus_P = 0.5 * (D_P_grid[i] + D_P_grid[i+1])
            D_minus_P = 0.5 * (D_P_grid[i-1] + D_P_grid[i])
            r_plus = r_i + 0.5 * dr
            r_minus = r_i - 0.5 * dr
            A_i_P = (D_minus_P * (r_minus ** g)) / ((r_i ** g) * (dr ** 2))
            D_prime_i_P = (D_plus_P * (r_plus ** g)) / ((r_i ** g) * (dr ** 2))
            B_i_P = - (A_i_P + D_prime_i_P)
            b_P[i] = 1.0 + dt_sec * k_i - dt_sec * B_i_P
            if i < Nr - 1:
                c_P[i] = - dt_sec * D_prime_i_P
            a_P[i-1] = - dt_sec * A_i_P

    # Томас-алгоритм
    def solve_tridiag(a_mat, b_mat, c_mat, rhs):
        n = len(rhs)
        cp = np.zeros(n - 1)
        dp = np.zeros(n)
        
        cp[0] = c_mat[0] / b_mat[0]
        dp[0] = rhs[0] / b_mat[0]
        
        for idx in range(1, n - 1):
            denom = b_mat[idx] - a_mat[idx-1] * cp[idx-1]
            cp[idx] = c_mat[idx] / denom
            dp[idx] = (rhs[idx] - a_mat[idx-1] * dp[idx-1]) / denom
            
        denom = b_mat[n-1] - a_mat[n-2] * cp[n-2]
        dp[n-1] = (rhs[n-1] - a_mat[n-2] * dp[n-2]) / denom
        
        x_sol = np.zeros(n)
        x_sol[n-1] = dp[n-1]
        for idx in range(n - 2, -1, -1):
            x_sol[idx] = dp[idx] - cp[idx] * x_sol[idx+1]
        return x_sol

    # Сохраняем начальное состояние
    t_history.append(0.0)
    C_interface_history.append(float(C_V[interface_idx]))
    C_P_interface_history.append(float(C_P[interface_idx]))
    p_boundary_history.append(p_base)
    
    v_sprout = 0.0
    
    # Шаг по времени
    for step in range(1, total_steps + 1):
        t_days = step * dt_days
        
        # VEGF step
        rhs_V = C_V[:Nr].copy()
        C_V_new = solve_tridiag(a, b, c, rhs_V)
        C_V[:Nr] = C_V_new
        C_V[Nr] = 0.0
        
        # PDGF step с замедленным выходом из PLGA
        src_pdgf = 0.0
        if t_days >= 0.0:
            src_pdgf = ((1.0 - pdgf_burst_fraction) * P_loaded_relative / (3.0 * sec_per_day)) * np.exp(-((t_days - 7.0) ** 2) / 4.0)
            
        rhs_P = C_P[:Nr].copy()
        rhs_P[r_grid[:Nr] <= R_outer_cm] += dt_sec * src_pdgf
        C_P_new = solve_tridiag(a_P, b_P, c_P, rhs_P)
        C_P[:Nr] = C_P_new
        C_P[Nr] = 0.0
        
        C_inf = C_V[interface_idx]
        C_inf_P = C_P[interface_idx]
        
        # Прорастание сосудов от VEGF
        activation = C_inf / (K_vegf + C_inf + 1e-5) if C_inf > 1e-5 else 0.0
        v_sprout += dt_days * beta_angiogenesis_eff * activation * (1.0 - v_sprout)
        
        # Стабилизация сосудов PDGF
        K_pdgf = 0.1
        F_stab = C_inf_P / (K_pdgf + C_inf_P + 1e-5) if C_inf_P > 1e-5 else 0.0
        
        # Сопряжение давления с созреванием капиллярной сети или гемодинамика AV-шунта
        if av_loop_flow:
            k_thrombo_0 = 0.05
            k_thrombo = k_thrombo_0 * (1.0 + (max(0.0, 1.5 - tau_blood) / 0.5)**2 + (max(0.0, tau_blood - 8.0) / 2.0)**2)
            if anticoagulation:
                k_thrombo_eff = k_thrombo * 0.1
            else:
                k_thrombo_eff = k_thrombo
            k_hyperplasia = 0.005
            k_occlusion = k_thrombo_eff + k_hyperplasia
            p_bound = 30.0 + (95.0 - 30.0) * np.exp(-k_occlusion * t_days)
        else:
            p_bound_instant = p_base + (p_max - p_base) * v_sprout * (0.3 + 0.7 * F_stab)
            p_bound = max(p_boundary_history[-1], p_bound_instant)
        
        t_history.append(t_days)
        C_interface_history.append(float(C_inf))
        C_P_interface_history.append(float(C_inf_P))
        p_boundary_history.append(float(p_bound))

    # Срезы профилей факторов

    if av_loop_flow:
        p_base = 95.0
        p_max = 95.0
        L_fibrosis_microns = 0.0

    if species.lower() == "human":
        saved_days = [0, 1, 7, 14, 30, 45, 60]
    else:
        saved_days = [0, 1, 3, 7, 14, 21]
    saved_profiles = {}
    saved_profiles_P = {}
    
    C_V_profile = np.zeros(Nr + 1)
    C_V_profile[r_grid <= R_outer_cm] = V_loaded_relative
    saved_profiles[0] = C_V_profile.copy()
    
    C_P_profile = np.zeros(Nr + 1)
    C_P_profile[r_grid <= R_outer_cm] = pdgf_burst_fraction * P_loaded_relative
    saved_profiles_P[0] = C_P_profile.copy()
    
    for step in range(1, total_steps + 1):
        t_days = round(step * dt_days, 1)
        
        rhs_V = C_V_profile[:Nr].copy()
        C_V_profile[:Nr] = solve_tridiag(a, b, c, rhs_V)
        C_V_profile[Nr] = 0.0
        
        src_pdgf = ((1.0 - pdgf_burst_fraction) * P_loaded_relative / (3.0 * sec_per_day)) * np.exp(-((t_days - 7.0) ** 2) / 4.0)
        rhs_P = C_P_profile[:Nr].copy()
        rhs_P[r_grid[:Nr] <= R_outer_cm] += dt_sec * src_pdgf
        C_P_profile[:Nr] = solve_tridiag(a_P, b_P, c_P, rhs_P)
        C_P_profile[Nr] = 0.0
        
        d_key = int(round(t_days))
        if d_key in saved_days:
            if d_key not in saved_profiles:
                saved_profiles[d_key] = C_V_profile.copy()
                saved_profiles_P[d_key] = C_P_profile.copy()

    return {
        "r_grid": r_grid * 1e4,  # см -> мкм
        "t": np.array(t_history),
        "C_interface": np.array(C_interface_history),
        "C_P_interface": np.array(C_P_interface_history),
        "p_boundary": np.array(p_boundary_history),
        "saved_profiles": saved_profiles,
        "saved_profiles_P": saved_profiles_P,
        "R_outer_microns": R_outer_microns,
        "L_fibrosis_microns": L_fibrosis_microns
    }

def run_neovascularization_sweep_oxygen(
    R_outer_microns,
    rho_million_per_ml,
    D_oxygen_coefficient,
    geometry="planar",
    L_fibrosis_microns=0.0,
    D_fibrosis=1.0e-5,
    V_loaded_relative=1.0,
    k_clear_tissue=15.0,
    beta_angiogenesis=0.15,
    K_vegf=0.1,
    p_base=30.0,
    p_max=60.0,
    days=21,
    P_loaded_relative=1.0,
    pdgf_burst_fraction=0.3,
    plga_acidification_factor=0.0,
    rho_mac_million_per_ml=0.0,
    catalase_activity_relative=1.0,
    catalase_half_life_days=1.5,
    buffer_capacity_mM=10.0,
    swelling_ratio=1.0,
    q_ogm_mmHg_per_sec=0.0,
    tethered_catalase=False,
    E_0=50.0,
    species="Mouse",
    phi_pfc=0.0,
    av_loop_flow=False,
    crispr_hypoimmune=False,
    t_pre_days=0.0,
    cd47_overexpression=False,
    tau_blood=5.0,
    anticoagulation=False,
    pO2_pfc_saturation=200.0,
    turnover_rate=0.0,
    complement_protection=False
):
    """
    Рассчитывает динамику выживаемости и секреции инсулина сопряженно с неоваскуляризацией (VEGF + PDGF).
    """
    neo_res = solve_coupled_neovascularization(
        R_outer_microns=R_outer_microns,
        rho_million_per_ml=rho_million_per_ml,
        D_oxygen_coefficient=D_oxygen_coefficient,
        geometry=geometry,
        L_fibrosis_microns=L_fibrosis_microns,
        D_fibrosis=D_fibrosis,
        V_loaded_relative=V_loaded_relative,
        k_clear_tissue=k_clear_tissue,
        beta_angiogenesis=beta_angiogenesis,
        K_vegf=K_vegf,
        p_base=p_base,
        p_max=p_max,
        days=days,
        P_loaded_relative=P_loaded_relative,
        pdgf_burst_fraction=pdgf_burst_fraction,
        plga_acidification_factor=plga_acidification_factor,
        species=species,
        av_loop_flow=av_loop_flow,
        tau_blood=tau_blood,
        anticoagulation=anticoagulation
    )
    
    p_bounds = neo_res["p_boundary"]
    t_history = neo_res["t"]
    
    viabilities = []
    insulin_capacities = []
    
    for idx, p_b in enumerate(p_bounds):
        t_d = t_history[idx]
        current_rho = 0.0 if t_d < t_pre_days else rho_million_per_ml
        oxy_res = solve_oxygen_profile(
            R_outer_microns=R_outer_microns,
            rho_million_per_ml=current_rho,
            p_boundary=p_b,
            D_coefficient=D_oxygen_coefficient,
            geometry=geometry,
            L_fibrosis_microns=L_fibrosis_microns,
            D_fibrosis=D_fibrosis,
            rho_mac_million_per_ml=rho_mac_million_per_ml,
            q_ogm_mmHg_per_sec=q_ogm_mmHg_per_sec,
            catalase_activity_relative=catalase_activity_relative,
            catalase_half_life_days=catalase_half_life_days,
            buffer_capacity_mM=buffer_capacity_mM,
            swelling_ratio=swelling_ratio,
            plga_acidification_factor=plga_acidification_factor,
            t_days=t_d,
            tethered_catalase=tethered_catalase,
            E_0=E_0,
            phi_pfc=phi_pfc,
            av_loop_flow=av_loop_flow,
            crispr_hypoimmune=crispr_hypoimmune,
            cd47_overexpression=cd47_overexpression,
            tau_blood=tau_blood,
            anticoagulation=anticoagulation,
            pO2_pfc_saturation=pO2_pfc_saturation,
            turnover_rate=turnover_rate,
            complement_protection=complement_protection
        )
        viabilities.append(oxy_res["viable_fraction"])
        insulin_capacities.append(oxy_res["insulin_capacity"])
        
    neo_res["viability_over_time"] = np.array(viabilities)
    neo_res["insulin_over_time"] = np.array(insulin_capacities)
    
    return neo_res

def solve_cytokine_profile(
    R_outer_microns=150,
    C_ext=10.0,
    D_cyt=1.0e-6,
    k_bind_scav=0.2,
    k_deg=0.01
):
    """
    Решает одномерную краевую задачу диффузии цитокинов с со-инкапсулированными
    белками-ловушками (scavengers) в сферическом гидрогеле.
    """
    from scipy.integrate import solve_bvp
    R_outer = R_outer_microns * 1e-4
    k_total = k_bind_scav + k_deg
    
    def equations(r, y):
        C = y[0]
        dC = y[1]
        
        # d2C/dr2 + 2/r * dC/dr - k_total/D_cyt * C = 0
        d2C = (k_total / D_cyt) * C - (2.0 / (r + 1e-9)) * dC
        return np.vstack((dC, d2C))
        
    def boundary_conditions(ya, yb):
        return np.array([ya[1] - 0.0, yb[0] - C_ext])
        
    r_mesh = np.linspace(0.0, R_outer, 100)
    y_init = np.zeros((2, r_mesh.size))
    y_init[0, :] = C_ext
    
    sol = solve_bvp(equations, boundary_conditions, r_mesh, y_init, tol=1e-5, max_nodes=500)
    
    z_coords_microns = np.linspace(0.0, R_outer_microns, 1000)
    C_profile = np.interp(z_coords_microns, sol.x * 1e4, sol.y[0])
    C_profile = np.maximum(0.0, C_profile)
    
    toxic_threshold = 1.0
    weights = z_coords_microns ** 2
    protected_mask = C_profile < toxic_threshold
    protected_fraction = (np.sum(protected_mask * weights) / np.sum(weights)) * 100.0
    
    return {
        "z": z_coords_microns,
        "C": C_profile,
        "protected_fraction": float(protected_fraction)
    }

def solve_cytokine_profile_transient(
    R_outer_microns=150,
    C_ext=10.0,
    D_cyt=1.0e-6,
    k_bind_scav=0.2,  # 1/(uM * s)
    k_deg=0.01,
    C_scav_0=5.0,  # uM
    k_deg_scav=0.001,  # rate of spontaneous scavenger degradation in gel
    swelling_ratio=1.0,
    days=14,
    shell_thickness_microns=50.0,
    coaxial_active=True,
    crispr_hypoimmune=False,
    viable_fraction=100.0,
    turnover_rate=0.0
):
    """
    Решает одномерную нестационарную задачу сопряженной диффузии цитокинов
    и расхода (истощения) белков-ловушек (scavengers) во времени.
    """
    # Рассчитываем выживаемость клеток от сдвигового напряжения
    if coaxial_active:
        eta_core = 0.1
    else:
        eta_core = 1.5  # Стандартная экструзия с высоковязким гелем
        
    eta_shell = 1.5
    Q_total = 2.0  # мл/мин
    
    r_core_m = R_outer_microns * 1e-6
    r_shell_m = (R_outer_microns + shell_thickness_microns) * 1e-6
    
    A_core = np.pi * (r_core_m ** 2)
    A_shell = np.pi * (r_shell_m ** 2 - r_core_m ** 2)
    Q_core = Q_total * (A_core / (A_core + A_shell))
    Q_shell = Q_total * (A_shell / (A_core + A_shell))
    
    Q_core_m3 = Q_core * 1.6667e-8
    Q_shell_m3 = Q_shell * 1.6667e-8
    
    tau_core_pa = (4.0 * eta_core * Q_core_m3) / (np.pi * (r_core_m ** 3) + 1e-30)
    tau_core_kpa = tau_core_pa / 1000.0
    survival_shear = 100.0 * (1.0 - 0.6 * np.tanh(tau_core_kpa / 5.0))
    Fraction_death_shear = 1.0 - survival_shear / 100.0
    Fraction_death_hypoxia = 1.0 - viable_fraction / 100.0

    R_outer_cm = R_outer_microns * (swelling_ratio ** (1.0 / 3.0)) * 1e-4
    D_cyt_eff = D_cyt * (swelling_ratio ** (2.0 / 3.0))
    
    Nr = 100
    r_grid = np.linspace(0.0, R_outer_cm, Nr + 1)
    dr = R_outer_cm / Nr
    
    # Инициализация концентраций
    C_cyt = np.zeros(Nr + 1)
    C_scav = np.ones(Nr + 1) * C_scav_0
    
    # Набухание увеличивает MWCO. Если MWCO >= 150 kDa, защищенность обнуляется.
    MWCO = 30.0 * (swelling_ratio ** 2)
    
    dt_days = 0.1
    sec_per_day = 86400.0
    dt_sec = dt_days * sec_per_day
    total_steps = int(days / dt_days)
    
    t_history = []
    C_center_history = []
    protected_fraction_history = []
    
    # Вспомогательный тридиагональный решатель
    def solve_tridiag_local(a_mat, b_mat, c_mat, rhs):
        n = len(rhs)
        cp = np.zeros(n - 1)
        dp = np.zeros(n)
        
        cp[0] = c_mat[0] / b_mat[0]
        dp[0] = rhs[0] / b_mat[0]
        
        for idx in range(1, n - 1):
            denom = b_mat[idx] - a_mat[idx-1] * cp[idx-1]
            cp[idx] = c_mat[idx] / denom
            dp[idx] = (rhs[idx] - a_mat[idx-1] * dp[idx-1]) / denom
            
        denom = b_mat[n-1] - a_mat[n-2] * cp[n-2]
        dp[n-1] = (rhs[n-1] - a_mat[n-2] * dp[n-2]) / denom
        
        x_sol = np.zeros(n)
        x_sol[n-1] = dp[n-1]
        for idx in range(n - 2, -1, -1):
            x_sol[idx] = dp[idx] - cp[idx] * x_sol[idx+1]
        return x_sol
        
    t_history.append(0.0)
    C_center_history.append(float(C_cyt[0]))
    protected_fraction_history.append(100.0)
    
    C_ext_timeline = [C_ext]
    
    saved_days = [0, 1, 3, 7, 14]
    saved_profiles = {}
    saved_profiles[0] = {
        "z": r_grid * 1e4,
        "C": C_cyt.copy(),
        "C_scav": C_scav.copy(),
        "C_ext": C_ext
    }
    
    for step in range(1, total_steps + 1):
        t_days = step * dt_days
        Turnover_death = viable_fraction * turnover_rate
        C_ext_basal = 15.0 * Turnover_death * (1.0 - np.exp(-0.5 * t_days))
        C_ext_t = C_ext + 15.0 * Fraction_death_shear * np.exp(-t_days / 2.0) + 15.0 * Fraction_death_hypoxia * np.exp(-t_days / 5.0) + C_ext_basal
        C_ext_timeline.append(C_ext_t)
        
        # 1. Обновление концентрации ловушек (расход + деградация)
        C_scav_new = C_scav / (1.0 + dt_sec * (k_bind_scav * C_cyt / 17000.0 + k_deg_scav / 86400.0))
        C_scav = np.maximum(0.0, C_scav_new)
        
        # 2. Неявный шаг для цитокинов (g = 2.0 для сферы)
        g = 2.0
        a_mat = np.zeros(Nr - 1)
        b_mat = np.zeros(Nr)
        c_mat = np.zeros(Nr - 1)
        
        for i in range(Nr):
            r_i = r_grid[i]
            k_i = k_bind_scav * C_scav[i] + k_deg
            
            if i == 0:
                B_i = - 6.0 * D_cyt_eff / (dr ** 2)
                D_prime_i = 6.0 * D_cyt_eff / (dr ** 2)
                b_mat[0] = 1.0 + dt_sec * k_i - dt_sec * B_i
                c_mat[0] = - dt_sec * D_prime_i
            else:
                r_plus = r_i + 0.5 * dr
                r_minus = r_i - 0.5 * dr
                
                A_i = (D_cyt_eff * (r_minus ** g)) / ((r_i ** g) * (dr ** 2))
                D_prime_i = (D_cyt_eff * (r_plus ** g)) / ((r_i ** g) * (dr ** 2))
                B_i = - (A_i + D_prime_i)
                
                b_mat[i] = 1.0 + dt_sec * k_i - dt_sec * B_i
                if i < Nr - 1:
                    c_mat[i] = - dt_sec * D_prime_i
                a_mat[i-1] = - dt_sec * A_i
                
        rhs = C_cyt[:Nr].copy()
        rhs[Nr - 1] += dt_sec * (D_cyt_eff * ((r_grid[Nr-1] + 0.5*dr)**g) / ((r_grid[Nr-1]**g) * (dr**2))) * C_ext_t
        
        C_cyt_new = solve_tridiag_local(a_mat, b_mat, c_mat, rhs)
        C_cyt[:Nr] = C_cyt_new
        C_cyt[Nr] = C_ext_t
        
        weights = r_grid ** 2
        
        if MWCO >= 150.0 and not crispr_hypoimmune:
            f_IgG_leak = 0.0
        else:
            f_IgG_leak = 1.0
            
        protected_mask = (C_cyt < 1.0).astype(float) * f_IgG_leak
        protected_fraction = (np.sum(protected_mask * weights) / np.sum(weights)) * 100.0
        
        t_history.append(t_days)
        C_center_history.append(float(C_cyt[0]))
        protected_fraction_history.append(float(protected_fraction))
        
        d_key = int(round(t_days))
        if d_key in saved_days:
            if d_key not in saved_profiles:
                saved_profiles[d_key] = {
                    "z": r_grid * 1e4,
                    "C": C_cyt.copy(),
                    "C_scav": C_scav.copy(),
                    "C_ext": C_ext_t
                }
                
    return {
        "t": np.array(t_history),
        "C_center": np.array(C_center_history),
        "protected_fraction_over_time": np.array(protected_fraction_history),
        "saved_profiles": saved_profiles,
        "R_outer_microns": R_outer_microns,
        "C_ext_timeline": np.array(C_ext_timeline),
        "survival_shear": survival_shear
    }

