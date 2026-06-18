import sys
sys.path.append("t1d_simulator")
import torch
import numpy as np
from scipy.integrate import solve_bvp
from simulator import SOLUBILITY, V_MAX
from diagnose_pinn import train_pinn_model_advanced_new

# Exact BVP solver using solve_bvp with sigmoid transition
def solve_advanced_oxygen_profile_bvp(
    R_outer_microns=150,
    rho_million_per_ml=80,
    p_boundary_base=30.0,
    D_gel=1.5e-5,
    geometry="spherical",
    L_fibrosis_microns=50.0,
    D_fibrosis=0.3*3e-5,
    rho_mac_million_per_ml=80.0,
    q_ogm_mmHg_per_sec=0.05
):
    R_outer = R_outer_microns * 1e-4
    L_fib = L_fibrosis_microns * 1e-4
    R_total = R_outer + L_fib
    
    rho_cells = rho_million_per_ml * 1e6
    rho_macs = rho_mac_million_per_ml * 1e6
    
    V_MAX_MAC = 3.0e-16
    K_M_MAC = 1.0
    
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
        R_cells_phys = (rho_cells * V_MAX / SOLUBILITY) * (p_pos / (0.5 + p_pos)) * sigmoid
        # Consumption in macrophages (mmHg/s)
        R_macs_phys = (rho_macs * V_MAX_MAC / SOLUBILITY) * (p_pos / (K_M_MAC + p_pos)) * (1.0 - sigmoid)
        # Source from OGM (mmHg/s)
        S_ogm_phys = q_ogm_mmHg_per_sec * sigmoid
        
        # ODE: d2p_dr2 + (g/r + D_prime/D) * dp - (R_cells + R_macs - S_ogm) / D = 0
        term_source_sink = (R_cells_phys + R_macs_phys - S_ogm_phys) / D_r
        d2p = term_source_sink - (g_factor / (r + 1e-9) + D_prime / D_r) * dp
        
        return np.vstack((dp, d2p))
        
    def boundary_conditions(ya, yb):
        # Symmetry at center: dp/dr(0) = 0
        # Dirichlet at boundary: p(R_total) = p_boundary
        return np.array([ya[1] - 0.0, yb[0] - p_boundary_base])
        
    r_mesh = np.linspace(0.0, R_total, 300)
    y_init = np.zeros((2, r_mesh.size))
    y_init[0, :] = p_boundary_base
    
    sol = solve_bvp(equations, boundary_conditions, r_mesh, y_init, tol=1e-5, max_nodes=2000)
    
    # Interpolate to 1000 points
    z_coords_microns = np.linspace(0.0, R_total * 1e4, 1000)
    p_profile = np.interp(z_coords_microns, sol.x * 1e4, sol.y[0])
    p_profile = np.maximum(0.0, p_profile)
    
    # Calculate cell viability in the core
    cell_mask = z_coords_microns <= R_outer_microns
    core_z = z_coords_microns[cell_mask]
    core_p = p_profile[cell_mask]
    
    if geometry == "planar":
        weights = np.ones_like(core_z)
    elif geometry == "cylindrical":
        weights = core_z
    else:
        weights = core_z ** 2
        
    viable_mask = core_p >= 0.5
    viable_fraction = (np.sum(viable_mask * weights) / np.sum(weights)) * 100.0
    
    return {
        "z": z_coords_microns,
        "pO2": p_profile,
        "viable_fraction": viable_fraction,
        "min_pO2": np.min(core_p)
    }

# Run BVP solver
res_bvp = solve_advanced_oxygen_profile_bvp()

# Run PINN solver
model = train_pinn_model_advanced_new(
    geometry="spherical",
    R_outer_microns=150,
    rho_million_per_ml=80,
    p_boundary=30.0,
    D_gel=1.5e-5,
    L_fibrosis_microns=50.0,
    D_fibrosis=0.3*3e-5,
    rho_mac_million_per_ml=80.0,
    q_ogm_mmHg_per_sec=0.05,
    epochs_adam=1500,
    lr_adam=0.005,
    max_iter_lbfgs=300
)

# Predict PINN on BVP grid points
device = next(model.parameters()).device
x_test = torch.tensor(res_bvp["z"] / 200.0, dtype=torch.float32).view(-1, 1).to(device)
with torch.no_grad():
    u_pred = model(x_test).cpu().numpy().flatten()
p_pinn = u_pred * 30.0

# Calculate MAE
mae = np.mean(np.abs(res_bvp["pO2"] - p_pinn))
print(f"BVP min pO2: {res_bvp['min_pO2']:.4f} | PINN min pO2: {np.min(p_pinn[res_bvp['z'] <= 150.0]):.4f}")
print(f"BVP viable: {res_bvp['viable_fraction']:.2f}% | PINN viable: {np.sum(p_pinn[res_bvp['z'] <= 150.0] >= 0.5) / np.sum(res_bvp['z'] <= 150.0) * 100:.2f}%")
print(f"Overall MAE: {mae:.4f} mmHg")

# Let's print some points
indices = [0, 200, 400, 600, 750, 900, 999]
for idx in indices:
    z_val = res_bvp["z"][idx]
    p_bvp_val = res_bvp["pO2"][idx]
    p_pinn_val = p_pinn[idx]
    print(f"z: {z_val:.1f} um | BVP pO2: {p_bvp_val:.4f} | PINN pO2: {p_pinn_val:.4f} | Diff: {abs(p_bvp_val - p_pinn_val):.4f}")
