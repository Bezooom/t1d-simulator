import numpy as np
import os
import json

# Physical and biological constants
SOLUBILITY = 1.34e-9  # mol/(cm^3 * mmHg)
V_MAX = 1.2e-16        # mol/(cell * s) for beta-cells
K_M = 0.5             # mmHg
K_M_INSULIN = 5.0     # mmHg

# Macrophage parameters (inflamed fibrosis)
V_MAX_MAC = 3.0e-16    # mol/(cell * s) - highly active macrophages
K_M_MAC = 1.0         # mmHg

def solve_advanced_oxygen_profile_fd(
    R_outer_microns=150,
    rho_million_per_ml=80,
    p_boundary_base=30.0,
    D_gel=1.5e-5,
    geometry="spherical",
    L_fibrosis_microns=50.0,
    D_fibrosis=0.3*3e-5,
    rho_mac_million_per_ml=50.0, # Active macrophages in fibrosis
    q_ogm_mmHg_per_sec=0.0       # Active oxygen generation from CaO2
):
    """
    Solves advanced BVP using SciPy's solve_bvp for reference.
    """
    from scipy.integrate import solve_bvp
    R_outer = R_outer_microns * 1e-4
    L_fib = L_fibrosis_microns * 1e-4
    R_total = R_outer + L_fib
    
    rho_cells = rho_million_per_ml * 1e6
    rho_macs = rho_mac_million_per_ml * 1e6
    
    w_trans = 2.0 * 1e-4
    
    if geometry == "cylindrical":
        g_factor = 1.0
    elif geometry == "spherical":
        g_factor = 2.0
    else:
        g_factor = 0.0
        
    def equations(r, y):
        # y[0] = p (mmHg)
        # y[1] = dp/dr (mmHg/cm)
        p = y[0]
        dp = y[1]
        
        p_pos = np.maximum(0.0, p)
        sigmoid = 1.0 / (1.0 + np.exp((r - R_outer) / w_trans))
        sigmoid_prime = -(1.0 / w_trans) * sigmoid * (1.0 - sigmoid)
        
        D_r = D_gel * sigmoid + D_fibrosis * (1.0 - sigmoid)
        D_prime = (D_gel - D_fibrosis) * sigmoid_prime
        
        # Consumption in cells (mmHg/s)
        R_cells_phys = (rho_cells * V_MAX / SOLUBILITY) * (p_pos / (K_M + p_pos)) * sigmoid
        # Consumption in macrophages (mmHg/s)
        R_macs_phys = (rho_macs * V_MAX_MAC / SOLUBILITY) * (p_pos / (K_M_MAC + p_pos)) * (1.0 - sigmoid)
        # Source from OGM (mmHg/s)
        S_ogm_phys = q_ogm_mmHg_per_sec * sigmoid
        
        # ODE: d2p_dr2 + (g/r + D_prime/D) * dp - (R_cells + R_macs - S_ogm) / D = 0
        term_source_sink = (R_cells_phys + R_macs_phys - S_ogm_phys) / D_r
        d2p = term_source_sink - (g_factor / (r + 1e-9) + D_prime / D_r) * dp
        
        return np.vstack((dp, d2p))
        
    def boundary_conditions(ya, yb):
        return np.array([ya[1] - 0.0, yb[0] - p_boundary_base])
        
    r_mesh = np.linspace(0.0, R_total, 300)
    y_init = np.zeros((2, r_mesh.size))
    y_init[0, :] = p_boundary_base
    
    sol = solve_bvp(equations, boundary_conditions, r_mesh, y_init, tol=1e-5, max_nodes=2000)
    
    # Interpolate to 1000 points
    z_coords_microns = np.linspace(0.0, R_total * 1e4, 1000)
    p_profile = np.interp(z_coords_microns, sol.x * 1e4, sol.y[0])
    p_profile = np.maximum(0.0, p_profile)
    
    # Calculate cell viability in the core (r <= R_outer)
    cell_mask = z_coords_microns <= R_outer_microns
    core_z = z_coords_microns[cell_mask]
    core_p = p_profile[cell_mask]
    
    if geometry == "planar":
        weights = np.ones_like(core_z)
    elif geometry == "cylindrical":
        weights = core_z
    else:
        weights = core_z ** 2
        
    viable_mask = core_p >= K_M
    viable_fraction = (np.sum(viable_mask * weights) / np.sum(weights)) * 100.0
    
    return {
        "r": z_coords_microns.tolist(),
        "p": p_profile.tolist(),
        "viable_fraction": float(viable_fraction),
        "min_pO2": float(np.min(core_p))
    }

def solve_cytokine_profile(
    R_outer_microns=150,
    C_ext=10.0,            # ng/ml
    D_cyt=1.0e-6,          # cm^2/s
    k_bind_scav=0.2,       # 1/s
    k_deg=0.01             # 1/s
):
    """
    Solves BVP for cytokine diffusion and scavenging inside the hydrogel core
    """
    R_outer = R_outer_microns * 1e-4
    N = 100
    r = np.linspace(1e-9, R_outer, N)
    dr = R_outer / (N - 1)
    
    C = np.ones(N) * C_ext
    k_total = k_bind_scav + k_deg
    
    for it in range(5000):
        C_old = C.copy()
        C[0] = C[1]
        
        for i in range(1, N - 1):
            # d2C/dr2 + 2/r * dC/dr = k_total/D_cyt * C
            term = (k_total / D_cyt) * C_old[i]
            C[i] = 0.5 * (C_old[i+1] + C_old[i-1] + (dr / r[i]) * (C_old[i+1] - C_old[i-1]) - (dr**2) * term)
            C[i] = np.maximum(0.0, C[i])
            
        C[N-1] = C_ext
        C = 0.9 * C + 0.1 * C_old
        if np.max(np.abs(C - C_old)) < 1e-6:
            break
            
    toxic_threshold = 1.0
    weights = r ** 2
    protected_mask = C < toxic_threshold
    protected_fraction = (np.sum(protected_mask * weights) / np.sum(weights)) * 100.0
    
    return {
        "r": (r * 1e4).tolist(),
        "C": C.tolist(),
        "protected_fraction": float(protected_fraction)
    }

def run_simulations():
    print("=== Запуск симуляции OGM и активного воспаления (FBR) ===")
    
    # Сценарий 1: Стандартный пассивный фиброз (мало макрофагов)
    res_control = solve_advanced_oxygen_profile_fd(
        R_outer_microns=150,
        rho_million_per_ml=80,
        p_boundary_base=30.0,
        L_fibrosis_microns=50.0,
        rho_mac_million_per_ml=5.0,   # Низкая воспаленность
        q_ogm_mmHg_per_sec=0.0
    )
    
    # Сценарий 2: Активное воспаление (макрофаги активно дышат)
    res_inflamed = solve_advanced_oxygen_profile_fd(
        R_outer_microns=150,
        rho_million_per_ml=80,
        p_boundary_base=30.0,
        L_fibrosis_microns=50.0,
        rho_mac_million_per_ml=120.0, # Огромная плотность активированных макрофагов
        q_ogm_mmHg_per_sec=0.0
    )
    
    # Сценарий 3: Спасение с помощью OGM (генерация O2)
    # CaO2 выделяет 0.05 mmHg/s кислорода
    res_saved = solve_advanced_oxygen_profile_fd(
        R_outer_microns=150,
        rho_million_per_ml=80,
        p_boundary_base=30.0,
        L_fibrosis_microns=50.0,
        rho_mac_million_per_ml=120.0,
        q_ogm_mmHg_per_sec=0.045       # Выделение кислорода OGM
    )
    
    print(f"Выживаемость клеток в ядре:")
    print(f"  - Пассивный фиброз: {res_control['viable_fraction']:.1f}% (мин pO2: {res_control['min_pO2']:.2f} mmHg)")
    print(f"  - Активное воспаление (макрофаги): {res_inflamed['viable_fraction']:.1f}% (мин pO2: {res_inflamed['min_pO2']:.2f} mmHg)")
    print(f"  - Воспаление + OGM (CaO2): {res_saved['viable_fraction']:.1f}% (мин pO2: {res_saved['min_pO2']:.2f} mmHg)")
    
    print("\n=== Запуск симуляции проникновения цитокинов ===")
    res_cyt_no_scav = solve_cytokine_profile(k_bind_scav=0.0)
    res_cyt_scav = solve_cytokine_profile(k_bind_scav=0.5)
    
    print(f"Доля защищенных клеток от цитокинов:")
    print(f"  - Без нейтрализации: {res_cyt_no_scav['protected_fraction']:.1f}%")
    print(f"  - С нейтрализацией (IL-1Ra / k=0.5): {res_cyt_scav['protected_fraction']:.1f}%")
    
    # Сохраняем в JSON
    data_to_save = {
        "oxygen": {
            "control": res_control,
            "inflamed": res_inflamed,
            "saved": res_saved
        },
        "cytokines": {
            "no_scavenger": res_cyt_no_scav,
            "scavenger": res_cyt_scav
        }
    }
    
    os.makedirs("scratch", exist_ok=True)
    with open("scratch/advanced_simulations.json", "w") as f:
        json.dump(data_to_save, f, indent=2)
    print("\nРезультаты симуляции успешно сохранены в scratch/advanced_simulations.json")

if __name__ == "__main__":
    run_simulations()
