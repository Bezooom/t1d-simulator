import sys
sys.path.append("t1d_simulator")
import torch
import torch.nn as nn
import numpy as np
from pinn_solver import PINN

def train_pinn_cytokines_perfect(R_outer_microns, D_cyt, k_bind_scav, k_deg, epochs_adam=1200):
    torch.manual_seed(42)
    np.random.seed(42)
    
    R_outer = R_outer_microns * 1e-4
    k_total = k_bind_scav + k_deg
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = PINN().to(device)
    
    # Grid clustered near 1.0
    t = torch.linspace(0.0, 1.0, 150, device=device)
    x_col = (1.0 - (1.0 - t)**2).view(-1, 1).requires_grad_(True)
    
    optimizer_adam = torch.optim.Adam(model.parameters(), lr=0.005)
    
    x_0 = torch.tensor([[0.0]], requires_grad=True, device=device)
    x_1 = torch.tensor([[1.0]], requires_grad=True, device=device)
    
    coeff = (R_outer**2) * k_total / D_cyt
    
    for epoch in range(epochs_adam):
        optimizer_adam.zero_grad()
        
        u_raw = model(x_col)
        u_1_raw = model(x_1)
        u = u_raw + (1.0 - u_1_raw) * (x_col ** 2)
        
        du_dx = torch.autograd.grad(u, x_col, torch.ones_like(u), create_graph=True)[0]
        d2u_dx2 = torch.autograd.grad(du_dx, x_col, torch.ones_like(du_dx), create_graph=True)[0]
        
        ode_residual = d2u_dx2 + (2.0 / (x_col + 1e-4)) * du_dx - coeff * u
        loss_ode = torch.mean(ode_residual ** 2)
        
        u_0_raw = model(x_0)
        u_0 = u_0_raw + (1.0 - u_1_raw) * (x_0 ** 2)
        du_dx_0 = torch.autograd.grad(u_0, x_0, torch.ones_like(u_0), create_graph=True)[0]
        loss_sym = du_dx_0 ** 2
        
        loss = loss_ode + loss_sym
        loss.backward()
        optimizer_adam.step()
        
    optimizer_lbfgs = torch.optim.LBFGS(
        model.parameters(), max_iter=250, lr=0.1,
        tolerance_grad=1e-7, tolerance_change=1e-9, line_search_fn="strong_wolfe"
    )
    
    def closure():
        optimizer_lbfgs.zero_grad()
        u_raw = model(x_col)
        u_1_raw = model(x_1)
        u = u_raw + (1.0 - u_1_raw) * (x_col ** 2)
        
        du_dx = torch.autograd.grad(u, x_col, torch.ones_like(u), create_graph=True)[0]
        d2u_dx2 = torch.autograd.grad(du_dx, x_col, torch.ones_like(du_dx), create_graph=True)[0]
        
        ode_res = d2u_dx2 + (2.0 / (x_col + 1e-4)) * du_dx - coeff * u
        loss_ode_lb = torch.mean(ode_res ** 2)
        
        u_0_raw = model(x_0)
        u_0 = u_0_raw + (1.0 - u_1_raw) * (x_0 ** 2)
        du_dx_0 = torch.autograd.grad(u_0, x_0, torch.ones_like(u_0), create_graph=True)[0]
        loss_sym_lb = du_dx_0 ** 2
        
        loss_val = loss_ode_lb + loss_sym_lb
        loss_val.backward()
        return loss_val
        
    optimizer_lbfgs.step(closure)
    return model

model = train_pinn_cytokines_perfect(R_outer_microns=150, D_cyt=1e-6, k_bind_scav=0.5, k_deg=0.01)
x_test = torch.tensor([[0.0], [0.5], [0.8], [0.9], [1.0]], dtype=torch.float32).to(next(model.parameters()).device)
with torch.no_grad():
    u_pred_raw = model(x_test)
    u_1_raw = model(torch.tensor([[1.0]], device=x_test.device))
    u_pred = (u_pred_raw + (1.0 - u_1_raw) * (x_test ** 2)).cpu().numpy().flatten()
    
# Exact solution
lmbda = np.sqrt((0.5+0.01)/1e-6) * 150e-4
exact = lambda x: (1.0 / (x + 1e-9)) * np.sinh(lmbda * x) / np.sinh(lmbda) if x > 0 else lmbda / np.sinh(lmbda)

print(f"x=0.0: PINN={u_pred[0]:.6f} | Exact={exact(0.0):.6f} | Diff={abs(u_pred[0]-exact(0.0)):.6f}")
print(f"x=0.5: PINN={u_pred[1]:.6f} | Exact={exact(0.5):.6f} | Diff={abs(u_pred[1]-exact(0.5)):.6f}")
print(f"x=0.8: PINN={u_pred[2]:.6f} | Exact={exact(0.8):.6f} | Diff={abs(u_pred[2]-exact(0.8)):.6f}")
print(f"x=0.9: PINN={u_pred[3]:.6f} | Exact={exact(0.9):.6f} | Diff={abs(u_pred[3]-exact(0.9)):.6f}")
print(f"x=1.0: PINN={u_pred[4]:.6f} | Exact={exact(1.0):.6f} | Diff={abs(u_pred[4]-exact(1.0)):.6f}")
