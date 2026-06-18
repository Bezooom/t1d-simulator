# -*- coding: utf-8 -*-
import numpy as np

def calculate_immune_leak(b2m_ko, ciita_ko, cd47_ki, cd55_cd59_ki, pdl1_ki):
    """
    Рассчитывает остаточную скорость гибели клеток от иммунной системы (k_immune_leak, 1/сутки)
    на основе конфигурации генетических модификаций CRISPR.
    """
    # Т-клеточный киллинг (блокируется при B2M-KO или защите PD-L1)
    k_T_cell = 0.10 * (1.0 - float(b2m_ko)) * (1.0 - float(pdl1_ki))
    
    # NK-киллерный лизис ("Missing Self" при MHC-I KO без CD47)
    k_NK = 0.05 * float(b2m_ko) * (1.0 - float(cd47_ki))
    
    # Активация хелперов и макрофагов (блокируется при CIITA-KO)
    k_helper = 0.03 * (1.0 - float(ciita_ko))
    
    # Разрушение системой комплемента (блокируется мембранными ингибиторами CD55/CD59)
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
    """
    Аналитически рассчитывает популяционную динамику N_cells(t) во времени.
    """
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
    """
    Рассчитывает парциальное давление кислорода pO2 в ядре мини-органоида во времени
    под влиянием прорастания внутренней сосудистой сети (зависит от EPCs).
    """
    k_vasc = 0.1 * (1.0 + 5.0 * phi_epc)
    t_days = np.array(t_days, dtype=float)
    pO2_core = p_portal * (1.0 - np.exp(-k_vasc * t_days))
    return pO2_core

def simulate_organoid_insulin(t_days, N_cells, pO2_core, K_M_insulin=5.0, hepatic_extraction=0.60):
    """
    Рассчитывает секрецию инсулина и уровни в воротной (portal) и системной (systemic) венах
    с учетом 60% First-Pass эффекта печени.
    """
    insulin_secreted = N_cells * pO2_core / (K_M_insulin + pO2_core + 1e-30)
    insulin_portal = insulin_secreted
    insulin_systemic = insulin_portal * (1.0 - hepatic_extraction)
    return insulin_portal, insulin_systemic
