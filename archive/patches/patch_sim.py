import re

with open('t1d_simulator/simulator.py', 'r') as f:
    code = f.read()

# 1. Add parameters to solve_oxygen_profile
code = code.replace(
    'tethered_catalase=False, E_0=50.0, species="Mouse"',
    'tethered_catalase=False, E_0=50.0, species="Mouse",\n    phi_pfc=0.0, av_loop_flow=False, crispr_hypoimmune=False'
)

# 2. Add av_loop_flow logic and SOLUBILITY_eff
av_loop_str = """
    if av_loop_flow:
        p_boundary = 95.0
        L_fibrosis_microns = 0.0

    SOLUBILITY_eff = SOLUBILITY * (1.0 - phi_pfc) + (20.0 * SOLUBILITY) * phi_pfc
    p_pfc_t = 760.0 * np.exp(-t_days / 2.0) if phi_pfc > 0.0 else 0.0
"""
code = code.replace(
    '    R_outer_microns_eff = R_outer_microns * (swelling_ratio ** (1.0 / 3.0))',
    av_loop_str + '\n    R_outer_microns_eff = R_outer_microns * (swelling_ratio ** (1.0 / 3.0))'
)

# Replace SOLUBILITY with SOLUBILITY_eff in solve_oxygen_profile ONLY (we do it carefully)
code = code.replace('P_H2O2 = 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY * 1e9', 'P_H2O2 = 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY_eff * 1e9')
code = code.replace('P_OH = 4.0 * q_ogm_mmHg_per_sec * SOLUBILITY * 1e9', 'P_OH = 4.0 * q_ogm_mmHg_per_sec * SOLUBILITY_eff * 1e9')
code = code.replace('P_acid = plga_acidification_factor * 0.2 * SOLUBILITY * 1e9', 'P_acid = plga_acidification_factor * 0.2 * SOLUBILITY_eff * 1e9')
code = code.replace('Ca_accum_mM = 1.2 + 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY * 1e9 * t_sec * 1e-3', 'Ca_accum_mM = 1.2 + 2.0 * q_ogm_mmHg_per_sec * SOLUBILITY_eff * 1e9 * t_sec * 1e-3')
code = code.replace('R_cells_phys = (rho_cells * V_max_cell / SOLUBILITY)', 'R_cells_phys = (rho_cells * V_max_cell / SOLUBILITY_eff)')
code = code.replace('R_macs_phys = (rho_macs * V_MAX_MAC / SOLUBILITY)', 'R_macs_phys = (rho_macs * V_MAX_MAC / SOLUBILITY_eff)')

# Update young modulus degradation
code = code.replace(
    'young_modulus_eff = E_0 * (swelling_ratio ** -2.0)',
    'young_modulus_eff = E_0 * (swelling_ratio ** -2.0) * np.exp(-0.01 * t_days)'
)

# Update tethered_catalase effect on rho_mac_boosted
code = code.replace(
    'rho_mac_boosted = rho_mac_million_per_ml * (1.0 + 2.0 * plga_acidification_factor)',
    'rho_mac_boosted = rho_mac_million_per_ml * (1.0 + 2.0 * plga_acidification_factor) * (1.0 + 1.5 * float(tethered_catalase))'
)

# Update pO2_profile with p_pfc_t
code = code.replace(
    'pO2_profile = np.maximum(0.0, pO2_profile)',
    'pO2_profile = np.maximum(0.0, pO2_profile) + p_pfc_t'
)
code = code.replace(
    'pO2_profile = C_profile / SOLUBILITY',
    'pO2_profile = C_profile / SOLUBILITY_eff'
)

# Update crispr_hypoimmune
code = code.replace(
    'f_IgG = 1.0 / (1.0 + np.exp((MWCO - 150.0) / 10.0))  # прорыв IgG при MWCO > 150 kDa',
    'f_IgG = np.ones_like(z_coords_microns) if crispr_hypoimmune else 1.0 / (1.0 + np.exp((MWCO - 150.0) / 10.0))'
)
code = code.replace(
    'f_IgG = 1.0 / (1.0 + np.exp((MWCO - 150.0) / 10.0))',
    'f_IgG = np.ones_like(z_uniform) if crispr_hypoimmune else 1.0 / (1.0 + np.exp((MWCO - 150.0) / 10.0))'
)

# Update solve_coupled_neovascularization
code = code.replace(
    'plga_acidification_factor=0.0,\n    species="Mouse"',
    'plga_acidification_factor=0.0,\n    species="Mouse",\n    av_loop_flow=False'
)
av_loop_neo_str = """
    if av_loop_flow:
        p_base = 95.0
        p_max = 95.0
        L_fibrosis_microns = 0.0

"""
code = code.replace(
    '    if species.lower() == "human":',
    av_loop_neo_str + '    if species.lower() == "human":'
)

# Update run_neovascularization_sweep_oxygen
code = code.replace(
    'E_0=50.0,\n    species="Mouse"',
    'E_0=50.0,\n    species="Mouse",\n    phi_pfc=0.0,\n    av_loop_flow=False,\n    crispr_hypoimmune=False,\n    t_pre_days=0.0'
)
code = code.replace(
    '        oxy_res = solve_oxygen_profile(\n            R_outer_microns=R_outer_microns,\n            rho_million_per_ml=rho_million_per_ml,',
    '        current_rho = 0.0 if t_d < t_pre_days else rho_million_per_ml\n        oxy_res = solve_oxygen_profile(\n            R_outer_microns=R_outer_microns,\n            rho_million_per_ml=current_rho,'
)
code = code.replace(
    'species=species\n    )',
    'species=species,\n        av_loop_flow=av_loop_flow\n    )'
)
code = code.replace(
    'tethered_catalase=tethered_catalase,\n            E_0=E_0\n        )',
    'tethered_catalase=tethered_catalase,\n            E_0=E_0,\n            phi_pfc=phi_pfc,\n            av_loop_flow=av_loop_flow,\n            crispr_hypoimmune=crispr_hypoimmune\n        )'
)

# Update solve_cytokine_profile_transient
code = code.replace(
    'shell_thickness_microns=50.0,\n    coaxial_active=True',
    'shell_thickness_microns=50.0,\n    coaxial_active=True,\n    crispr_hypoimmune=False,\n    viable_fraction=100.0'
)
code = code.replace(
    'Fraction_death = 1.0 - survival_shear / 100.0',
    'Fraction_death_shear = 1.0 - survival_shear / 100.0\n    Fraction_death_hypoxia = 1.0 - viable_fraction / 100.0'
)
code = code.replace(
    'C_ext_t = C_ext + 15.0 * Fraction_death * (np.exp(0.5 * t_days) - 1.0)',
    'C_ext_t = C_ext + 15.0 * Fraction_death_shear * np.exp(-t_days / 2.0) + 15.0 * Fraction_death_hypoxia * np.exp(-t_days / 5.0)'
)
code = code.replace(
    'if MWCO >= 150.0:\n            f_IgG_leak = 0.0',
    'if MWCO >= 150.0 and not crispr_hypoimmune:\n            f_IgG_leak = 0.0'
)

# Standard BVP solver fix
code = code.replace(
    'C_boundary = p_boundary * SOLUBILITY',
    'C_boundary = p_boundary * SOLUBILITY_eff'
)
code = code.replace(
    'K_m_conc = K_M * SOLUBILITY',
    'K_m_conc = K_M * SOLUBILITY_eff'
)

with open('t1d_simulator/simulator.py', 'w') as f:
    f.write(code)

