import re

with open('t1d_simulator/pinn_solver.py', 'r') as f:
    code = f.read()

# 1. train_pinn_model
code = code.replace(
    'epochs_adam=1200, lr_adam=0.005, max_iter_lbfgs=300):',
    'epochs_adam=1200, lr_adam=0.005, max_iter_lbfgs=300, phi_pfc=0.0, av_loop_flow=False):'
)
code = code.replace(
    '    R_outer_cm = R_outer_microns * 1e-4\n    L_fibrosis_cm = L_fibrosis_microns * 1e-4',
    '    if av_loop_flow:\n        p_boundary = 95.0\n        L_fibrosis_microns = 0.0\n    R_outer_cm = R_outer_microns * 1e-4\n    L_fibrosis_cm = L_fibrosis_microns * 1e-4'
)
code = code.replace(
    'C_boundary = p_boundary * SOLUBILITY\n    K_m_conc = K_M * SOLUBILITY',
    'SOLUBILITY_eff = SOLUBILITY * (1.0 - phi_pfc) + (20.0 * SOLUBILITY) * phi_pfc\n    C_boundary = p_boundary * SOLUBILITY_eff\n    K_m_conc = K_M * SOLUBILITY_eff'
)

# 2. solve_oxygen_profile_pinn
code = code.replace(
    'tethered_catalase=False, E_0=50.0, species="Mouse"',
    'tethered_catalase=False, E_0=50.0, species="Mouse",\n    phi_pfc=0.0, av_loop_flow=False, crispr_hypoimmune=False'
)
code = code.replace(
    'tethered_catalase=tethered_catalase,\n            E_0=E_0',
    'tethered_catalase=tethered_catalase,\n            E_0=E_0,\n            phi_pfc=phi_pfc,\n            av_loop_flow=av_loop_flow,\n            crispr_hypoimmune=crispr_hypoimmune'
)
code = code.replace(
    '    # Масштабирование размеров и диффузии из-за набухания геля',
    '    if av_loop_flow:\n        p_boundary = 95.0\n        L_fibrosis_microns = 0.0\n    SOLUBILITY_eff = SOLUBILITY * (1.0 - phi_pfc) + (20.0 * SOLUBILITY) * phi_pfc\n    p_pfc_t = 760.0 * np.exp(-t_days / 2.0) if phi_pfc > 0.0 else 0.0\n    # Масштабирование размеров и диффузии из-за набухания геля'
)
code = code.replace(
    'D_fibrosis=D_fibrosis\n    )',
    'D_fibrosis=D_fibrosis,\n        phi_pfc=phi_pfc,\n        av_loop_flow=av_loop_flow\n    )'
)
code = code.replace(
    'pO2_profile = p_boundary * u_pred',
    'pO2_profile = p_boundary * u_pred'
)
code = code.replace(
    'pO2_profile = np.maximum(0.0, pO2_profile)',
    'pO2_profile = np.maximum(0.0, pO2_profile) + p_pfc_t'
)
code = code.replace(
    'young_modulus_eff = E_0 * (swelling_ratio ** -2.0)',
    'young_modulus_eff = E_0 * (swelling_ratio ** -2.0) * np.exp(-0.01 * t_days)'
)
code = code.replace(
    'f_IgG = 1.0 / (1.0 + np.exp((MWCO - 150.0) / 10.0))',
    'f_IgG = np.ones_like(z_coords_microns) if crispr_hypoimmune else 1.0 / (1.0 + np.exp((MWCO - 150.0) / 10.0))'
)
code = code.replace(
    'Ca_accum_mM = 1.2 + 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY',
    'Ca_accum_mM = 1.2 + 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY_eff'
)

# 3. train_pinn_model_advanced
code = code.replace(
    'epochs_adam=1500, lr_adam=0.005, max_iter_lbfgs=300):',
    'epochs_adam=1500, lr_adam=0.005, max_iter_lbfgs=300, phi_pfc=0.0, av_loop_flow=False):'
)
code = code.replace(
    '    R_outer = R_outer_microns_eff * 1e-4',
    '    if av_loop_flow:\n        p_boundary = 95.0\n        L_fibrosis_microns_eff = 0.0\n    R_outer = R_outer_microns_eff * 1e-4'
)
code = code.replace(
    'rho_mac_boosted = rho_mac_million_per_ml * (1.0 + 2.0 * plga_acidification_factor)',
    'rho_mac_boosted = rho_mac_million_per_ml * (1.0 + 2.0 * plga_acidification_factor) * (1.0 + 1.5 * float(tethered_catalase))'
)
code = code.replace(
    'C_boundary = p_boundary * SOLUBILITY',
    'SOLUBILITY_eff = SOLUBILITY * (1.0 - phi_pfc) + (20.0 * SOLUBILITY) * phi_pfc\n    C_boundary = p_boundary * SOLUBILITY_eff'
)
code = code.replace(
    'Q_ogm_mol = q_ogm_eff * SOLUBILITY',
    'Q_ogm_mol = q_ogm_eff * SOLUBILITY_eff'
)

# 4. solve_advanced_oxygen_profile_pinn
code = code.replace(
    'tethered_catalase=False,\n    E_0=50.0\n):',
    'tethered_catalase=False,\n    E_0=50.0,\n    phi_pfc=0.0,\n    av_loop_flow=False,\n    crispr_hypoimmune=False\n):'
)
code = code.replace(
    'E_0=E_0\n    )',
    'E_0=E_0,\n        phi_pfc=phi_pfc,\n        av_loop_flow=av_loop_flow\n    )'
)
code = code.replace(
    '    R_outer_microns_eff = R_outer_microns * (swelling_ratio ** (1.0 / 3.0))',
    '    if av_loop_flow:\n        p_boundary_base = 95.0\n        L_fibrosis_microns = 0.0\n    SOLUBILITY_eff = SOLUBILITY * (1.0 - phi_pfc) + (20.0 * SOLUBILITY) * phi_pfc\n    p_pfc_t = 760.0 * np.exp(-t_days / 2.0) if phi_pfc > 0.0 else 0.0\n    R_outer_microns_eff = R_outer_microns * (swelling_ratio ** (1.0 / 3.0))'
)
code = code.replace(
    'pO2_profile = np.maximum(0.0, pO2_profile)',
    'pO2_profile = np.maximum(0.0, pO2_profile) + p_pfc_t'
)
code = code.replace(
    'young_modulus_eff = E_0 * (swelling_ratio ** -2.0)',
    'young_modulus_eff = E_0 * (swelling_ratio ** -2.0) * np.exp(-0.01 * t_days)'
)
code = code.replace(
    'f_IgG = 1.0 / (1.0 + np.exp((MWCO - 150.0) / 10.0))',
    'f_IgG = np.ones_like(z_coords_microns) if crispr_hypoimmune else 1.0 / (1.0 + np.exp((MWCO - 150.0) / 10.0))'
)
code = code.replace('P_H2O2 = 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY', 'P_H2O2 = 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY_eff')
code = code.replace('P_OH = 4.0 * q_ogm_mmHg_per_sec * SOLUBILITY', 'P_OH = 4.0 * q_ogm_mmHg_per_sec * SOLUBILITY_eff')
code = code.replace('P_acid = plga_acidification_factor * 0.2 * SOLUBILITY', 'P_acid = plga_acidification_factor * 0.2 * SOLUBILITY_eff')
code = code.replace('Ca_accum_mM = 1.2 + 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY', 'Ca_accum_mM = 1.2 + 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY_eff')

with open('t1d_simulator/pinn_solver.py', 'w') as f:
    f.write(code)

