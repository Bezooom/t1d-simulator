import re
import os

with open('t1d_simulator/app.py', 'r') as f:
    code = f.read()

# 1. Add phi_pfc in 1D Simulation mode
code = code.replace(
    '        q_ogm = st.sidebar.slider(',
    '        phi_pfc = st.sidebar.slider(\n            "Фторуглеродная эмульсия PFC (Доля)",\n            0.0, 0.3, 0.0,\n            step=0.05,\n            help="Пассивный кислородный буфер (истощается за 2-4 дня)."\n        )\n        q_ogm = st.sidebar.slider('
)

# 2. Add crispr_hypoimmune in 1D Simulation mode
code = code.replace(
    '    st.sidebar.header("🛡️ Защита от цитокинов (Ловушки)")\n    C_ext = st.sidebar.slider(',
    '    st.sidebar.header("🛡️ Защита от цитокинов (Ловушки)")\n    crispr_hypoimmune = st.sidebar.checkbox(\n        "CRISPR Hypoimmune Клетки",\n        value=False,\n        help="Клетки с нокаутом MHC/HLA. Невидимы для IgG и защищены от антител."\n    )\n    C_ext = st.sidebar.slider('
)

# 3. Add parameters to solve_oxygen_profile and solve_cytokine_profile_transient
code = code.replace(
    'tethered_catalase=tethered_catalase,\n                E_0=E_0,\n                species="Mouse"\n            )',
    'tethered_catalase=tethered_catalase,\n                E_0=E_0,\n                species="Mouse",\n                phi_pfc=phi_pfc,\n                crispr_hypoimmune=crispr_hypoimmune\n            )'
)
code = code.replace(
    'coaxial_active=coaxial_fiber\n            )',
    'coaxial_active=coaxial_fiber,\n                crispr_hypoimmune=crispr_hypoimmune,\n                viable_fraction=res["viable_fraction"]\n            )'
)

# Replace the non-PINN branch
code = code.replace(
    'tethered_catalase=tethered_catalase,\n            E_0=E_0,\n            species="Mouse"\n        )',
    'tethered_catalase=tethered_catalase,\n            E_0=E_0,\n            species="Mouse",\n            phi_pfc=phi_pfc,\n            crispr_hypoimmune=crispr_hypoimmune\n        )'
)
code = code.replace(
    'coaxial_active=coaxial_fiber\n        )',
    'coaxial_active=coaxial_fiber,\n            crispr_hypoimmune=crispr_hypoimmune,\n            viable_fraction=res["viable_fraction"]\n        )'
)

# Hide IgG warning if crispr_hypoimmune
code = code.replace(
    '        if MWCO >= 150.0:\n            st.error(\n                "⚠️ **КРИТИЧЕСКИЙ ПРОРЫВ IgG!** Гель сильно набух ($S_f \\\\ge 2.24$). "\n                "Поры расширились до размера антител IgG (150 кДа), что полностью уничтожает иммуноизоляцию "\n                "и вызывает атаку иммуноглобулинов на β-клетки."\n            )',
    '        if MWCO >= 150.0:\n            if crispr_hypoimmune:\n                st.success("✅ **CRISPR Hypoimmune защита!** Гель сильно набух, но клетки генетически невидимы для антител (IgG).")\n            else:\n                st.error(\n                    "⚠️ **КРИТИЧЕСКИЙ ПРОРЫВ IgG!** Гель сильно набух ($S_f \\\\ge 2.24$). "\n                    "Поры расширились до размера антител IgG (150 кДа), что полностью уничтожает иммуноизоляцию "\n                    "и вызывает атаку иммуноглобулинов на β-клетки."\n                )'
)

# 4. Modify Angiogenesis mode (else branch)
angio_sidebar_code = """
    st.sidebar.header("🩸 Продвинутая имплантация")
    
    prevascularization = st.sidebar.checkbox(
        "Двухэтапная Предваскуляризация",
        value=False,
        help="Сначала имплантируется пустое устройство, сосуды прорастают, затем вводятся клетки."
    )
    
    t_pre_days = 0
    if prevascularization:
        t_pre_days = st.sidebar.slider(
            "Время предваскуляризации (дней)",
            1, 60, 14,
            step=1,
            help="Количество дней до инъекции клеток."
        )

    av_loop_flow = st.sidebar.checkbox(
        "Подключить прямой кровоток (AV loop)",
        value=False,
        help="Хирургическое подключение артерии и вены. Гарантирует pO2 = 95.0 mmHg без лаг-фазы."
    )
"""

code = code.replace(
    '    st.sidebar.header("⏱️ Параметры Ангиогенеза")',
    angio_sidebar_code + '\n    st.sidebar.header("⏱️ Параметры Ангиогенеза")'
)

# Update run_neovascularization_sweep_oxygen in Angiogenesis mode
code = code.replace(
    '        res_vegf = run_neovascularization_sweep_oxygen(\n            R_outer_microns=r_vegf,\n            rho_million_per_ml=rho_vegf,\n            D_oxygen_coefficient=d_vegf,\n            geometry=geom_vegf,\n            L_fibrosis_microns=l_fib,\n            D_fibrosis=1.0e-5 * d_fib_mult,\n            V_loaded_relative=v_loaded,\n            k_clear_tissue=k_clear,\n            beta_angiogenesis=beta_angio,\n            K_vegf=k_vegf,\n            p_base=p_base,\n            p_max=p_max_angio,\n            days=sim_days,\n            P_loaded_relative=p_loaded,\n            pdgf_burst_fraction=pdgf_burst,\n            plga_acidification_factor=plga_acid,\n            species=species_selected\n        )',
    '        res_vegf = run_neovascularization_sweep_oxygen(\n            R_outer_microns=r_vegf,\n            rho_million_per_ml=rho_vegf,\n            D_oxygen_coefficient=d_vegf,\n            geometry=geom_vegf,\n            L_fibrosis_microns=l_fib,\n            D_fibrosis=1.0e-5 * d_fib_mult,\n            V_loaded_relative=v_loaded,\n            k_clear_tissue=k_clear,\n            beta_angiogenesis=beta_angio,\n            K_vegf=k_vegf,\n            p_base=p_base,\n            p_max=p_max_angio,\n            days=sim_days,\n            P_loaded_relative=p_loaded,\n            pdgf_burst_fraction=pdgf_burst,\n            plga_acidification_factor=plga_acid,\n            species=species_selected,\n            av_loop_flow=av_loop_flow,\n            t_pre_days=t_pre_days\n        )'
)

with open('t1d_simulator/app.py', 'w') as f:
    f.write(code)
