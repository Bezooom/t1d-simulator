import sys
sys.path.append("t1d_simulator")
sys.path.append("scratch")
import numpy as np
from pinn_solver import solve_cytokine_profile_pinn
from simulate_advanced_biodesign import solve_cytokine_profile

res_fd = solve_cytokine_profile(
    R_outer_microns=150,
    C_ext=10.0,
    D_cyt=1.0e-6,
    k_bind_scav=0.5,
    k_deg=0.01
)

res_pinn = solve_cytokine_profile_pinn(
    R_outer_microns=150,
    C_ext=10.0,
    D_cyt=1.0e-6,
    k_bind_scav=0.5,
    k_deg=0.01
)

mae = np.mean(np.abs(np.interp(res_fd["r"], res_pinn["z"], res_pinn["C"]) - res_fd["C"]))
print(f"FD boundary: {res_fd['C'][-1]:.4f} | PINN boundary: {res_pinn['C'][-1]:.4f}")
print(f"FD center: {res_fd['C'][0]:.4f} | PINN center: {res_pinn['C'][0]:.4f}")
print(f"MAE: {mae:.4f}")
