import sys
sys.path.append("t1d_simulator")
import torch
import torch.nn as nn
import numpy as np
from simulator import SOLUBILITY, V_MAX

class PINN(nn.Module):
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

def train_pinn_model_advanced_new(
    geometry, R_outer_microns, rho_million_per_ml, p_boundary, D_gel, 
    L_fibrosis_microns, D_fibrosis, rho_mac_million_per_ml, q_ogm_mmHg_per_sec,
    epochs_adam=1500, lr_adam=0.005, max_iter_lbfgs=300, bc_weight=150.0
):
    torch.manual_seed(42)
    np.random.seed(42)
    
    R_outer = R_outer_microns * 1e-4
    L_fib = L_fibrosis_microns * 1e-4
    R_total = R_outer + L_fib
    
    rho_cells = rho_million_per_ml * 1e6
    rho_macs = rho_mac_million_per_ml * 1e6
    
    C_boundary = p_boundary * SOLUBILITY
    V_MAX_MAC = 3.0e-16
    K_M_MAC = 1.0
    Q_ogm_mol = q_ogm_mmHg_per_sec * SOLUBILITY
    
    # Pre-calculate scaled phi2 (O(1) values)
    phi2_cells = (rho_cells * V_MAX * (R_total ** 2)) / (D_gel * C_boundary)
    phi2_macs = (rho_macs * V_MAX_MAC * (R_total ** 2)) / (D_fibrosis * C_boundary)
    phi2_ogm = (Q_ogm_mol * (R_total ** 2)) / (D_gel * C_boundary)
    
    kappa = 0.5 / p_boundary
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
    
    for epoch in range(epochs_adam):
        optimizer_adam.zero_grad()
        
        u = model(x_col)
        du_dx = torch.autograd.grad(u, x_col, torch.ones_like(u), create_graph=True)[0]
        d2u_dx2 = torch.autograd.grad(du_dx, x_col, torch.ones_like(du_dx), create_graph=True)[0]
        
        r = x_col * R_total
        sigmoid = 1.0 / (1.0 + torch.exp((r - R_outer) / w_trans))
        D_r = D_gel * sigmoid + D_fibrosis * (1.0 - sigmoid)
        
        R_cells_dim = phi2_cells * (D_gel / D_r) * (u / (kappa + u)) * sigmoid
        R_macs_dim = phi2_macs * (D_fibrosis / D_r) * (u / (kappa_mac + u)) * (1.0 - sigmoid)
        S_ogm_dim = phi2_ogm * (D_gel / D_r) * sigmoid
        
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
        
        if epoch % 300 == 0:
            print(f"Epoch {epoch} | Loss ODE: {loss_ode.item():.4e} | Loss BC: {loss_bc.item():.4e} | Total: {loss.item():.4e}")
        
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
        D_r = D_gel * sigmoid + D_fibrosis * (1.0 - sigmoid)
        
        R_cells_dim = phi2_cells * (D_gel / D_r) * (u / (kappa + u)) * sigmoid
        R_macs_dim = phi2_macs * (D_fibrosis / D_r) * (u / (kappa_mac + u)) * (1.0 - sigmoid)
        S_ogm_dim = phi2_ogm * (D_gel / D_r) * sigmoid
        
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

# Test outputs
z_coords = np.linspace(0.0, 200.0, 10)
x_test = torch.tensor(z_coords / 200.0, dtype=torch.float32).view(-1, 1).to(next(model.parameters()).device)
with torch.no_grad():
    u_pred = model(x_test).cpu().numpy().flatten()
print("x values (scaled):", z_coords / 200.0)
print("u_pred values (scaled):", u_pred)
print("pO2 values:", u_pred * 30.0)
