import numpy as np
import torch
import torch.nn as nn
from simulator import SOLUBILITY, V_MAX, K_M, K_M_INSULIN

# --- 1. Архитектура нейросети PINN ---
class PINN(nn.Module):
    """
    Полносвязная нейросеть для аппроксимации безразмерного профиля кислорода u(x).
    Принимает на вход координату x в [0, 1] и возвращает u(x) в [0, 1].
    Использует Tanh для гладкости вторых производных и Sigmoid для удержания диапазона.
    """
    def __init__(self, hidden_dim=32):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(1, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, 1),
            nn.Sigmoid()
        )
        
    def forward(self, x):
        return self.net(x)

# --- 2. Функция обучения PINN ---
def train_pinn_model(geometry, R_outer_microns, rho_million_per_ml, p_boundary, D_coefficient, 
                     V_max_cell=V_MAX, L_fibrosis_microns=0.0, D_fibrosis=1.0e-5, 
                     epochs_adam=1200, lr_adam=0.005, max_iter_lbfgs=300, phi_pfc=0.0, av_loop_flow=False,
                     tau_blood=5.0, anticoagulation=False, t_days=0.0):
    """
    Обучает модель PINN для решения уравнения диффузии кислорода в безразмерных координатах.
    """
    # Установка семян для детерминизма обучения
    torch.manual_seed(42)
    np.random.seed(42)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(42)

    # Перевод в безразмерные переменные
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
    R_outer_cm = R_outer_microns * 1e-4
    L_fibrosis_cm = L_fibrosis_microns * 1e-4
    rho = rho_million_per_ml * 1e6
    SOLUBILITY_eff = SOLUBILITY * (1.0 - phi_pfc) + (20.0 * SOLUBILITY) * phi_pfc
    SOLUBILITY_eff = SOLUBILITY * (1.0 - phi_pfc) + (20.0 * SOLUBILITY) * phi_pfc
    C_boundary = p_boundary * SOLUBILITY_eff
    K_m_conc = K_M * SOLUBILITY_eff
    
    # 1. Модуль Тиле (квадрат)
    phi2 = (rho * V_max_cell * (R_outer_cm ** 2)) / (D_coefficient * C_boundary)
    # 2. Безразмерная константа Михаэлиса
    kappa = K_m_conc / C_boundary
    # 3. Число Био (диффузионное сопротивление)
    if L_fibrosis_microns > 0:
        Biot = (D_fibrosis * R_outer_cm) / (D_coefficient * L_fibrosis_cm)
    else:
        Biot = 0.0
        
    # Фактор геометрии кривизны (g): 0 - slab, 1 - cylinder, 2 - sphere
    if geometry == "cylindrical":
        g_factor = 1.0
    elif geometry == "spherical":
        g_factor = 2.0
    else:
        g_factor = 0.0
        
    # Инициализация модели
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = PINN().to(device)
    
    # Сэмплирование точек коллокации в интервале (0, 1]
    # Начинаем с 1e-4 во избежание деления на ноль при g_factor / x
    x_col = torch.linspace(1e-4, 1.0, 150, requires_grad=True).view(-1, 1).to(device)
    
    # Оптимизатор 1: Adam для глобального поиска минимума
    optimizer_adam = torch.optim.Adam(model.parameters(), lr=lr_adam)
    
    # Точки для граничных условий
    x_0 = torch.tensor([[0.0]], requires_grad=True, device=device)
    x_1 = torch.tensor([[1.0]], requires_grad=True, device=device)
    
    for epoch in range(epochs_adam):
        optimizer_adam.zero_grad()
        
        # Расчет невязки ОДУ
        u_col = model(x_col)
        du_dx = torch.autograd.grad(u_col, x_col, torch.ones_like(u_col), create_graph=True)[0]
        d2u_dx2 = torch.autograd.grad(du_dx, x_col, torch.ones_like(du_dx), create_graph=True)[0]
        
        # Регуляризованное деление (x_col + 1e-4) для исключения сингулярности в центре
        ode_residual = d2u_dx2 + (g_factor / (x_col + 1e-4)) * du_dx - phi2 * (u_col / (kappa + u_col))
        loss_ode = torch.mean(ode_residual ** 2)
        
        # Граничное условие симметрии в центре (x = 0): du/dx = 0
        u_0 = model(x_0)
        du_dx_0 = torch.autograd.grad(u_0, x_0, torch.ones_like(u_0), create_graph=True)[0]
        loss_sym = du_dx_0 ** 2
        
        # Граничное условие на внешней границе (x = 1)
        u_1 = model(x_1)
        if Biot > 0.0:
            # Условие Робина: du/dx = Bi * (1 - u)
            du_dx_1 = torch.autograd.grad(u_1, x_1, torch.ones_like(u_1), create_graph=True)[0]
            loss_bc = (du_dx_1 - Biot * (1.0 - u_1)) ** 2
        else:
            # Условие Дирихле: u(1) = 1
            loss_bc = (u_1 - 1.0) ** 2
            
        # Общая функция потерь
        loss = loss_ode + loss_sym + 15.0 * loss_bc
        loss.backward()
        optimizer_adam.step()
        
    # Оптимизатор 2: L-BFGS для точной сходимости
    optimizer_lbfgs = torch.optim.LBFGS(model.parameters(), max_iter=max_iter_lbfgs, lr=0.1, 
                                        tolerance_grad=1e-7, tolerance_change=1e-9, 
                                        line_search_fn="strong_wolfe")
    
    def closure():
        optimizer_lbfgs.zero_grad()
        
        # Точки коллокации
        u_col_lb = model(x_col)
        du_dx_lb = torch.autograd.grad(u_col_lb, x_col, torch.ones_like(u_col_lb), create_graph=True)[0]
        d2u_dx2_lb = torch.autograd.grad(du_dx_lb, x_col, torch.ones_like(du_dx_lb), create_graph=True)[0]
        
        ode_res_lb = d2u_dx2_lb + (g_factor / (x_col + 1e-4)) * du_dx_lb - phi2 * (u_col_lb / (kappa + u_col_lb))
        loss_ode_lb = torch.mean(ode_res_lb ** 2)
        
        # Граничные условия
        u_0_lb = model(x_0)
        du_dx_0_lb = torch.autograd.grad(u_0_lb, x_0, torch.ones_like(u_0_lb), create_graph=True)[0]
        loss_sym_lb = du_dx_0_lb ** 2
        
        u_1_lb = model(x_1)
        if Biot > 0.0:
            du_dx_1_lb = torch.autograd.grad(u_1_lb, x_1, torch.ones_like(u_1_lb), create_graph=True)[0]
            loss_bc_lb = (du_dx_1_lb - Biot * (1.0 - u_1_lb)) ** 2
        else:
            loss_bc_lb = (u_1_lb - 1.0) ** 2
            
        loss_val = loss_ode_lb + loss_sym_lb + 15.0 * loss_bc_lb
        loss_val.backward()
        return loss_val
        
    optimizer_lbfgs.step(closure)
    
    return model

# --- 3. Обертка-решатель (Drop-in Replacement) ---
def solve_oxygen_profile_pinn(
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
    Интерфейсный обертчик над обученной нейросетью PINN.
    При указании rho_mac_million_per_ml > 0, q_ogm_mmHg_per_sec > 0 или plga_acidification_factor > 0
    делегирует вычисления усовершенствованному multi-physics решателю.
    """
    if rho_mac_million_per_ml > 0.0 or q_ogm_mmHg_per_sec > 0.0 or plga_acidification_factor > 0.0:
        return solve_advanced_oxygen_profile_pinn(
            R_outer_microns=R_outer_microns,
            rho_million_per_ml=rho_million_per_ml,
            p_boundary_base=p_boundary,
            D_gel=D_coefficient,
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
            t_days=t_days,
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
        
    # 1. Обучаем стандартную модель PINN
    model = train_pinn_model(
        geometry=geometry,
        R_outer_microns=R_outer_microns_eff,
        rho_million_per_ml=rho_million_per_ml,
        p_boundary=p_boundary,
        D_coefficient=D_coefficient_eff,
        V_max_cell=V_max_cell,
        L_fibrosis_microns=L_fibrosis_microns_eff,
        D_fibrosis=D_fibrosis,
        phi_pfc=phi_pfc,
        av_loop_flow=av_loop_flow,
        tau_blood=tau_blood,
        anticoagulation=anticoagulation,
        t_days=t_days
    )
    
    # 2. Генерируем равномерную сетку (1000 точек) для предсказания и интегрирования
    z_coords_microns = np.linspace(0.0, R_outer_microns_eff, 1000)
    x_test = torch.tensor(z_coords_microns / R_outer_microns_eff, dtype=torch.float32).view(-1, 1)
    
    device = next(model.parameters()).device
    x_test = x_test.to(device)
    
    with torch.no_grad():
        u_pred = model(x_test).cpu().numpy().flatten()
        
    # Перевод безразмерной концентрации обратно в физическое pO2 (mmHg)
    # pO2 = (C_boundary * u_pred) / SOLUBILITY = p_boundary * u_pred
    pO2_profile = p_boundary * u_pred
    pO2_profile = np.maximum(0.0, pO2_profile) + p_pfc_t
    
    # 3. Вычисление пространственно-взвешенных интегральных характеристик
    if geometry == "planar":
        weights = np.ones_like(z_coords_microns)
    elif geometry == "cylindrical":
        weights = z_coords_microns
    elif geometry == "spherical":
        weights = z_coords_microns ** 2
    else:
        weights = np.ones_like(z_coords_microns)
        
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
    
    f_total = f_IgG * np.ones_like(z_coords_microns) * f_Ca * f_rupture * f_NK * f_complement * f_hyperoxia * f_depletion
    
    # Доля выживших клеток (pO2 >= K_M = 0.5 mmHg)
    viable_mask = pO2_profile >= K_M
    viable_fraction_profile = viable_mask.astype(float) * f_total
    viable_fraction = (np.sum(viable_fraction_profile * weights) / weights_sum) * 100.0
    
    # Функциональная емкость инсулиносекреции
    insulin_capacity_profile = pO2_profile / (K_M_INSULIN + pO2_profile)
    max_possible_capacity = 95.0 / (K_M_INSULIN + 95.0)
    mean_insulin_capacity = np.sum(insulin_capacity_profile * weights) / weights_sum
    insulin_capacity = (mean_insulin_capacity / max_possible_capacity) * 100.0
    insulin_capacity = min(insulin_capacity, viable_fraction)
    
    min_pO2 = np.min(pO2_profile)
    
    return {
        "z": z_coords_microns,
        "pO2": pO2_profile,
        "pH": np.ones_like(z_coords_microns) * 7.4,
        "H2O2": np.zeros_like(z_coords_microns),
        "viability_multiplier": f_total,
        "viable_fraction": viable_fraction,
        "insulin_capacity": insulin_capacity,
        "min_pO2": min_pO2,
        "Ca_accum_mM": Ca_accum_mM,
        "insulin_leak": insulin_leak,
        "rupture_risk": rupture_risk,
        "young_modulus_eff": young_modulus_eff
    }

# --- 4. Усовершенствованный multi-physics PINN-решатель для OGM и активного воспаления ---
def train_pinn_model_advanced(
    geometry, R_outer_microns, rho_million_per_ml, p_boundary, D_gel, 
    L_fibrosis_microns, D_fibrosis, rho_mac_million_per_ml, q_ogm_mmHg_per_sec,
    catalase_activity_relative=1.0, catalase_half_life_days=1.5,
    buffer_capacity_mM=10.0, swelling_ratio=1.0,
    plga_acidification_factor=0.0, t_days=0.0,
    tethered_catalase=False, E_0=50.0,
    epochs_adam=1500, lr_adam=0.005, max_iter_lbfgs=300, phi_pfc=0.0, av_loop_flow=False,
    tau_blood=5.0, anticoagulation=False, pO2_pfc_saturation=200.0
):
    torch.manual_seed(42)
    np.random.seed(42)
    
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
        p_boundary_base = 30.0 + (95.0 - 30.0) * np.exp(-k_occlusion * t_days)
        p_boundary = p_boundary_base
        L_fibrosis_microns = 0.0
    SOLUBILITY_eff = SOLUBILITY * (1.0 - phi_pfc) + (20.0 * SOLUBILITY) * phi_pfc
    p_pfc_t = pO2_pfc_saturation * np.exp(-t_days / 2.0) if phi_pfc > 0.0 else 0.0
    R_outer_microns_eff = R_outer_microns * (swelling_ratio ** (1.0 / 3.0))
    L_fibrosis_microns_eff = L_fibrosis_microns * (swelling_ratio ** (1.0 / 3.0))
    D_gel_eff = D_gel * (swelling_ratio ** (2.0 / 3.0))
    
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
        L_fibrosis_microns_eff = 0.0
    R_outer = R_outer_microns_eff * 1e-4
    L_fib = L_fibrosis_microns_eff * 1e-4
    R_total = R_outer + L_fib
    
    rho_cells = rho_million_per_ml * 1e6
    
    # Закисление PLGA усиливает плотность макрофагов
    rho_mac_boosted = rho_mac_million_per_ml * (1.0 + 2.0 * plga_acidification_factor) * (1.0 + 1.5 * float(tethered_catalase))
    rho_macs = rho_mac_boosted * 1e6
    
    SOLUBILITY_eff = SOLUBILITY * (1.0 - phi_pfc) + (20.0 * SOLUBILITY) * phi_pfc
    C_boundary = p_boundary * SOLUBILITY_eff
    V_MAX_MAC = 3.0e-16
    K_M_MAC = 1.0
    
    # Каталаза и ее полураспад
    if tethered_catalase:
        catalase_half_life_days_eff = 100.0
        catalase_activity_relative_eff = catalase_activity_relative * 0.25
    else:
        catalase_half_life_days_eff = catalase_half_life_days / (swelling_ratio ** 2)
        catalase_activity_relative_eff = catalase_activity_relative

    k_cat_0 = 0.05
    k_cat = k_cat_0 * catalase_activity_relative_eff * np.exp(-t_days * np.log(2.0) / catalase_half_life_days_eff)
    q_ogm_eff = q_ogm_mmHg_per_sec * (k_cat / (k_cat_0 + 1e-5))
    Q_ogm_mol = q_ogm_eff * SOLUBILITY_eff
    
    # Pre-calculate dimensionless scale factors (O(1) values)
    phi2_cells = (rho_cells * V_MAX * (R_total ** 2)) / (D_gel_eff * C_boundary)
    phi2_macs = (rho_macs * V_MAX_MAC * (R_total ** 2)) / (D_fibrosis * C_boundary)
    phi2_ogm = (Q_ogm_mol * (R_total ** 2)) / (D_gel_eff * C_boundary)
    
    kappa = K_M / p_boundary
    kappa_mac = K_M_MAC / p_boundary
    
    if geometry == "cylindrical":
        g_factor = 1.0
    elif geometry == "spherical":
        g_factor = 2.0
    else:
        g_factor = 0.0
        
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = PINN().to(device)
    
    x_col = torch.linspace(1e-4, 1.0, 200, requires_grad=True).view(-1, 1).to(device)
    optimizer_adam = torch.optim.Adam(model.parameters(), lr=lr_adam)
    
    x_0 = torch.tensor([[0.0]], requires_grad=True, device=device)
    x_1 = torch.tensor([[1.0]], requires_grad=True, device=device)
    
    w_trans = 2.0 * 1e-4
    bc_weight = 150.0
    
    for epoch in range(epochs_adam):
        optimizer_adam.zero_grad()
        
        u = model(x_col)
        du_dx = torch.autograd.grad(u, x_col, torch.ones_like(u), create_graph=True)[0]
        d2u_dx2 = torch.autograd.grad(du_dx, x_col, torch.ones_like(du_dx), create_graph=True)[0]
        
        r = x_col * R_total
        sigmoid = 1.0 / (1.0 + torch.exp((r - R_outer) / w_trans))
        D_r = D_gel_eff * sigmoid + D_fibrosis * (1.0 - sigmoid)
        
        R_cells_dim = phi2_cells * (D_gel_eff / D_r) * (u / (kappa + u)) * sigmoid
        R_macs_dim = phi2_macs * (D_fibrosis / D_r) * (u / (kappa_mac + u)) * (1.0 - sigmoid)
        S_ogm_dim = phi2_ogm * (D_gel_eff / D_r) * sigmoid
        
        ode_residual = d2u_dx2 + (g_factor / (x_col + 1e-4)) * du_dx - (R_cells_dim + R_macs_dim - S_ogm_dim)
        loss_ode = torch.mean(ode_residual ** 2)
        
        u_0 = model(x_0)
        du_dx_0 = torch.autograd.grad(u_0, x_0, torch.ones_like(u_0), create_graph=True)[0]
        loss_sym = du_dx_0 ** 2
        
        u_1 = model(x_1)
        loss_bc = (u_1 - 1.0) ** 2
        
        loss = loss_ode + loss_sym + bc_weight * loss_bc
        loss.backward()
        optimizer_adam.step()
        
    optimizer_lbfgs = torch.optim.LBFGS(
        model.parameters(), max_iter=max_iter_lbfgs, lr=0.1, 
        tolerance_grad=1e-7, tolerance_change=1e-9, line_search_fn="strong_wolfe"
    )
    
    def closure():
        optimizer_lbfgs.zero_grad()
        u = model(x_col)
        du_dx = torch.autograd.grad(u, x_col, torch.ones_like(u), create_graph=True)[0]
        d2u_dx2 = torch.autograd.grad(du_dx, x_col, torch.ones_like(du_dx), create_graph=True)[0]
        
        r = x_col * R_total
        sigmoid = 1.0 / (1.0 + torch.exp((r - R_outer) / w_trans))
        D_r = D_gel_eff * sigmoid + D_fibrosis * (1.0 - sigmoid)
        
        R_cells_dim = phi2_cells * (D_gel_eff / D_r) * (u / (kappa + u)) * sigmoid
        R_macs_dim = phi2_macs * (D_fibrosis / D_r) * (u / (kappa_mac + u)) * (1.0 - sigmoid)
        S_ogm_dim = phi2_ogm * (D_gel_eff / D_r) * sigmoid
        
        ode_res = d2u_dx2 + (g_factor / (x_col + 1e-4)) * du_dx - (R_cells_dim + R_macs_dim - S_ogm_dim)
        loss_ode_lb = torch.mean(ode_res ** 2)
        
        u_0_lb = model(x_0)
        du_dx_0_lb = torch.autograd.grad(u_0_lb, x_0, torch.ones_like(u_0_lb), create_graph=True)[0]
        loss_sym_lb = du_dx_0_lb ** 2
        
        u_1_lb = model(x_1)
        loss_bc_lb = (u_1_lb - 1.0) ** 2
        
        loss_val = loss_ode_lb + loss_sym_lb + bc_weight * loss_bc_lb
        loss_val.backward()
        return loss_val
        
    optimizer_lbfgs.step(closure)
    return model

def solve_advanced_oxygen_profile_pinn(
    R_outer_microns=150,
    rho_million_per_ml=80,
    p_boundary_base=30.0,
    D_gel=1.5e-5,
    geometry="spherical",
    L_fibrosis_microns=50.0,
    D_fibrosis=0.3*3e-5,
    rho_mac_million_per_ml=50.0,
    q_ogm_mmHg_per_sec=0.0,
    catalase_activity_relative=1.0,
    catalase_half_life_days=1.5,
    buffer_capacity_mM=10.0,
    swelling_ratio=1.0,
    plga_acidification_factor=0.0,
    t_days=0.0,
    tethered_catalase=False,
    E_0=50.0,
    phi_pfc=0.0,
    av_loop_flow=False,
    crispr_hypoimmune=False,
    cd47_overexpression=False,
    tau_blood=5.0,
    anticoagulation=False,
    pO2_pfc_saturation=200.0,
    turnover_rate=0.0,
    complement_protection=False
):
    model = train_pinn_model_advanced(
        geometry=geometry,
        R_outer_microns=R_outer_microns,
        rho_million_per_ml=rho_million_per_ml,
        p_boundary=p_boundary_base,
        D_gel=D_gel,
        L_fibrosis_microns=L_fibrosis_microns,
        D_fibrosis=D_fibrosis,
        rho_mac_million_per_ml=rho_mac_million_per_ml,
        q_ogm_mmHg_per_sec=q_ogm_mmHg_per_sec,
        catalase_activity_relative=catalase_activity_relative,
        catalase_half_life_days=catalase_half_life_days,
        buffer_capacity_mM=buffer_capacity_mM,
        swelling_ratio=swelling_ratio,
        plga_acidification_factor=plga_acidification_factor,
        t_days=t_days,
        tethered_catalase=tethered_catalase,
        E_0=E_0,
        phi_pfc=phi_pfc,
        av_loop_flow=av_loop_flow,
        tau_blood=tau_blood,
        anticoagulation=anticoagulation,
        pO2_pfc_saturation=pO2_pfc_saturation
    )
    
    if av_loop_flow:
        k_thrombo_0 = 0.05
        k_thrombo = k_thrombo_0 * (1.0 + (max(0.0, 1.5 - tau_blood) / 0.5)**2 + (max(0.0, tau_blood - 8.0) / 2.0)**2)
        if anticoagulation:
            k_thrombo_eff = k_thrombo * 0.1
        else:
            k_thrombo_eff = k_thrombo
        k_hyperplasia = 0.005
        k_occlusion = k_thrombo_eff + k_hyperplasia
        p_boundary_base = 30.0 + (95.0 - 30.0) * np.exp(-k_occlusion * t_days)
        L_fibrosis_microns = 0.0
    SOLUBILITY_eff = SOLUBILITY * (1.0 - phi_pfc) + (20.0 * SOLUBILITY) * phi_pfc
    p_pfc_t = pO2_pfc_saturation * np.exp(-t_days / 2.0) if phi_pfc > 0.0 else 0.0
    R_outer_microns_eff = R_outer_microns * (swelling_ratio ** (1.0 / 3.0))
    L_fibrosis_microns_eff = L_fibrosis_microns * (swelling_ratio ** (1.0 / 3.0))
    R_total_microns = R_outer_microns_eff + L_fibrosis_microns_eff
    
    z_coords_microns = np.linspace(0.0, R_total_microns, 1000)
    x_test = torch.tensor(z_coords_microns / R_total_microns, dtype=torch.float32).view(-1, 1)
    
    device = next(model.parameters()).device
    x_test = x_test.to(device)
    
    with torch.no_grad():
        u_pred = model(x_test).cpu().numpy().flatten()
        
    pO2_profile = p_boundary_base * u_pred
    pO2_profile = np.maximum(0.0, pO2_profile) + p_pfc_t
    
    # Calculate viable fraction (only in the cell core r <= R_outer_microns_eff)
    cell_mask = z_coords_microns <= R_outer_microns_eff
    core_z = z_coords_microns[cell_mask]
    core_p = pO2_profile[cell_mask]
    
    if geometry == "planar":
        weights = np.ones_like(core_z)
    elif geometry == "cylindrical":
        weights = core_z
    else:
        weights = core_z ** 2
        
    weights_sum = np.sum(weights)
    
    # Моделирование каталазы и её распада
    if tethered_catalase:
        catalase_half_life_days_eff = 100.0
        catalase_activity_relative_eff = catalase_activity_relative * 0.25
    else:
        catalase_half_life_days_eff = catalase_half_life_days / (swelling_ratio ** 2)
        catalase_activity_relative_eff = catalase_activity_relative

    k_cat_0 = 0.05
    k_cat = k_cat_0 * catalase_activity_relative_eff * np.exp(-t_days * np.log(2.0) / catalase_half_life_days_eff)
    
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
    
    # Расчет пространственного распределения H2O2 и pH (затухание к границам)
    H2O2_profile = C_H2O2_core * (1.0 - (z_coords_microns / R_total_microns) ** 2)
    H2O2_profile = np.maximum(0.0, H2O2_profile)
    
    pH_profile = 7.4 + (pH_core - 7.4) * (1.0 - (z_coords_microns / R_total_microns) ** 2)
    
    f_pH = np.exp(- (pH_profile - 7.4) ** 2 / (2.0 * 0.3 ** 2))
    f_ROS = 1.0 / (1.0 + (H2O2_profile / 10.0) ** 2)
    
    # Расчет накопления кальция Ca2+ в ядре (ммоль/л)
    Ca_accum_mM = 1.2 + 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY_eff * 1e9 * t_sec * 1e-3
    f_Ca = np.exp(- (np.maximum(0.0, Ca_accum_mM - 2.0) ** 2) / (2.0 * (1.5 ** 2)))
    insulin_leak = np.minimum(1.0, np.maximum(0.0, Ca_accum_mM - 1.2) / 4.0)
    
    # Механические расчеты
    young_modulus_eff = E_0 * (swelling_ratio ** -2.0) * np.exp(-0.01 * t_days) * np.exp(-0.01 * t_days)
    sigma_ext = 2.0
    epsilon_max = 0.1
    rupture_risk = min(100.0, (sigma_ext / (young_modulus_eff * epsilon_max)) * 100.0)
    
    if rupture_risk > 50.0:
        f_rupture = 1.0 - rupture_risk / 100.0
    else:
        f_rupture = 1.0

    MWCO = 30.0 * (swelling_ratio ** 2)
    f_IgG = np.ones_like(z_coords_microns) if crispr_hypoimmune else 1.0 / (1.0 + np.exp((MWCO - 150.0) / 10.0))  # прорыв IgG при MWCO > 150 kDa
    
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

# --- 5. PINN-решатель для диффузии и связывания цитокинов ---
def train_pinn_cytokines(R_outer_microns, D_cyt, k_bind_scav, k_deg, epochs_adam=1200, max_iter_lbfgs=200):
    torch.manual_seed(42)
    np.random.seed(42)
    
    R_outer = R_outer_microns * 1e-4
    k_total = k_bind_scav + k_deg
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = PINN().to(device)
    
    # Clustered grid near boundary (1.0) to resolve boundary layer
    t = torch.linspace(0.0, 1.0, 150, device=device)
    x_col = (1.0 - (1.0 - t)**2).view(-1, 1).requires_grad_(True)
    
    optimizer_adam = torch.optim.Adam(model.parameters(), lr=0.005)
    
    x_0 = torch.tensor([[0.0]], requires_grad=True, device=device)
    x_1 = torch.tensor([[1.0]], requires_grad=True, device=device)
    
    coeff = (R_outer**2) * k_total / D_cyt
    bc_weight = 500.0
    
    for epoch in range(epochs_adam):
        optimizer_adam.zero_grad()
        
        u = model(x_col)
        du_dx = torch.autograd.grad(u, x_col, torch.ones_like(u), create_graph=True)[0]
        d2u_dx2 = torch.autograd.grad(du_dx, x_col, torch.ones_like(du_dx), create_graph=True)[0]
        
        ode_residual = d2u_dx2 + (2.0 / (x_col + 1e-4)) * du_dx - coeff * u
        loss_ode = torch.mean(ode_residual ** 2)
        
        u_0 = model(x_0)
        du_dx_0 = torch.autograd.grad(u_0, x_0, torch.ones_like(u_0), create_graph=True)[0]
        loss_sym = du_dx_0 ** 2
        
        u_1 = model(x_1)
        loss_bc = (u_1 - 1.0) ** 2
        
        loss = loss_ode + loss_sym + bc_weight * loss_bc
        loss.backward()
        optimizer_adam.step()
        
    optimizer_lbfgs = torch.optim.LBFGS(
        model.parameters(), max_iter=max_iter_lbfgs, lr=0.1,
        tolerance_grad=1e-7, tolerance_change=1e-9, line_search_fn="strong_wolfe"
    )
    
    def closure():
        optimizer_lbfgs.zero_grad()
        u = model(x_col)
        du_dx = torch.autograd.grad(u, x_col, torch.ones_like(u), create_graph=True)[0]
        d2u_dx2 = torch.autograd.grad(du_dx, x_col, torch.ones_like(du_dx), create_graph=True)[0]
        
        ode_res = d2u_dx2 + (2.0 / (x_col + 1e-4)) * du_dx - coeff * u
        loss_ode_lb = torch.mean(ode_res ** 2)
        
        u_0_lb = model(x_0)
        du_dx_0_lb = torch.autograd.grad(u_0_lb, x_0, torch.ones_like(u_0_lb), create_graph=True)[0]
        loss_sym_lb = du_dx_0_lb ** 2
        
        u_1_lb = model(x_1)
        loss_bc_lb = (u_1_lb - 1.0) ** 2
        
        loss_val = loss_ode_lb + loss_sym_lb + bc_weight * loss_bc_lb
        loss_val.backward()
        return loss_val
        
    optimizer_lbfgs.step(closure)
    return model

def solve_cytokine_profile_pinn(R_outer_microns=150, C_ext=10.0, D_cyt=1.0e-6, k_bind_scav=0.2, k_deg=0.01):
    model = train_pinn_cytokines(
        R_outer_microns=R_outer_microns,
        D_cyt=D_cyt,
        k_bind_scav=k_bind_scav,
        k_deg=k_deg
    )
    
    z_coords_microns = np.linspace(0.0, R_outer_microns, 1000)
    x_test = torch.tensor(z_coords_microns / R_outer_microns, dtype=torch.float32).view(-1, 1)
    
    device = next(model.parameters()).device
    x_test = x_test.to(device)
    
    with torch.no_grad():
        u_pred = model(x_test).cpu().numpy().flatten()
        
    C_profile = C_ext * u_pred
    C_profile = np.maximum(0.0, C_profile)
    
    # Calculate protected fraction (where cytokines < 1.0 ng/ml)
    toxic_threshold = 1.0
    weights = z_coords_microns ** 2
    protected_mask = C_profile < toxic_threshold
    protected_fraction = (np.sum(protected_mask * weights) / np.sum(weights)) * 100.0
    
    return {
        "z": z_coords_microns,
        "C": C_profile,
        "protected_fraction": protected_fraction
    }
