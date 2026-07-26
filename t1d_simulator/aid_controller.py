# -*- coding: utf-8 -*-
"""
Automated Insulin Delivery (AID) / Closed-Loop Artificial Pancreas Simulator
Bergman Minimal Model + PID Controller for CGM Time in Range (TIR) optimization.
"""

import numpy as np

def simulate_aid_closed_loop(
    simulation_hours=24.0,
    meal_carbs_g=60.0,
    meal_time_hours=2.0,
    kp=0.03,
    ki=0.0005,
    kd=0.05,
    target_glucose_mg_dl=110.0,
    basal_rate_u_h=0.72
):
    """
    Simulates a 24-hour closed-loop AID artificial pancreas system.
    Converts glucose units (mg/dL <-> mmol/L) and calculates TIR, TBR, TAR metrics.
    """
    dt_min = 1.0  # 1-minute time step
    total_steps = int(simulation_hours * 60)
    time_mins = np.arange(total_steps) * dt_min
    time_hours = time_mins / 60.0
    
    # Bergman Minimal Model Constants
    p1 = 0.028   # 1/min
    p2 = 0.025   # 1/min
    p3 = 0.000013 # (uU/mL)^-1 min^-2
    n = 0.09     # 1/min
    Gb = 110.0   # mg/dL
    Ib = 15.0    # uU/mL
    
    # State vectors
    G = np.zeros(total_steps) # Glucose (mg/dL)
    X = np.zeros(total_steps) # Remote insulin action
    I = np.zeros(total_steps) # Plasma insulin (uU/mL)
    infusion_rate = np.zeros(total_steps) # U/hr
    
    G[0] = Gb
    X[0] = 0.0
    I[0] = Ib
    
    integral_err = 0.0
    prev_err = 0.0
    
    meal_step = int(meal_time_hours * 60)
    
    for step in range(total_steps - 1):
        # 1. Meal appearance (mg/dL / min)
        if step >= meal_step and step < (meal_step + 180):
            t_m = (step - meal_step)
            D_meal = (meal_carbs_g * 1000.0 / 12000.0) * (t_m / 40.0) * np.exp(-t_m / 40.0) * 18.0
        else:
            D_meal = 0.0
            
        # 2. PID Closed-Loop Controller
        current_g = G[step]
        err = current_g - target_glucose_mg_dl
        integral_err += err * dt_min
        integral_err = np.clip(integral_err, -5000.0, 5000.0)
        deriv_err = (err - prev_err) / dt_min
        prev_err = err
        
        u_pid = basal_rate_u_h + kp * err + ki * integral_err + kd * deriv_err
        u_pid = max(0.0, min(12.0, u_pid)) # Clamp pump infusion 0-12 U/hr
        infusion_rate[step] = u_pid
        
        # Convert U/hr to uU/mL/min plasma input
        u_plasma_input = (u_pid * 1e6 / 60.0) / 12000.0
        
        # 3. ODE Step (Euler integration)
        dG = - (p1 + X[step]) * G[step] + p1 * Gb + D_meal
        dX = - p2 * X[step] + p3 * (I[step] - Ib)
        dI = - n * (I[step] - Ib) + u_plasma_input
        
        G[step + 1] = max(30.0, G[step] + dG * dt_min)
        X[step + 1] = max(0.0, X[step] + dX * dt_min)
        I[step + 1] = max(0.0, I[step] + dI * dt_min)
        
    infusion_rate[-1] = infusion_rate[-2]
    
    # Convert Glucose mg/dL to mmol/L (divide by 18.015)
    glucose_mmol_l = G / 18.015
    
    # CGM Time in Range metrics
    tir_pct = float(np.mean((glucose_mmol_l >= 3.9) & (glucose_mmol_l <= 10.0)) * 100.0)
    tbr_pct = float(np.mean(glucose_mmol_l < 3.9) * 100.0)
    tar_pct = float(np.mean(glucose_mmol_l > 10.0) * 100.0)
    mean_g = float(np.mean(glucose_mmol_l))
    gmi_hba1c = float(12.71 + 4.705 * mean_g) / 3.303 # Approximate HbA1c %
    
    return {
        "time_hours": time_hours,
        "glucose_mmol_l": glucose_mmol_l,
        "glucose_mg_dl": G,
        "infusion_rate_u_h": infusion_rate,
        "TIR_percent": tir_pct,
        "TBR_percent": tbr_pct,
        "TAR_percent": tar_pct,
        "mean_glucose_mmol_l": mean_g
    }
