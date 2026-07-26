"""
Panel 1: 1D Oxygen Diffusion Simulation (PDE + PINN modes).
Handles sidebar parameters, solver dispatch, metrics, and all plots.
"""
import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from t1d_simulator.simulator import (
    solve_oxygen_profile,
    solve_cytokine_profile_transient,
    run_density_sweep,
)
from t1d_simulator.simulator import (
    HYDROGELS,
    IMPLANTATION_SITES,
    K_M,
    K_M_INSULIN,
    V_MAX,
)
from t1d_simulator.ui.layout import GEOMETRY_FORMAT


def render_oxygen_panel():
    """Render the full 1D Oxygen Simulation panel. Returns the result dict."""
    # --- Solver type ---
    st.sidebar.write("---")
    st.sidebar.header("🧠 Вычислительное ядро")
    solver_type = st.sidebar.radio(
        "Решатель ОДУ:",
        ["SciPy (solve_bvp)", "PyTorch (PINN)"],
        help="PINN обучается в реальном времени под физические параметры системы.",
    )

    # --- Geometry ---
    st.sidebar.header("📐 Выбор топологии (Геометрия)")
    geometry_key = st.sidebar.radio(
        "Форма имплантата:",
        ["planar", "cylindrical", "spherical"],
        format_func=lambda x: GEOMETRY_FORMAT[x],
    )

    # --- Implantation site ---
    st.sidebar.header("📍 Место имплантации (Диффузия O₂)")
    site_key = st.sidebar.selectbox(
        "Выберите зону пересадки:",
        options=list(IMPLANTATION_SITES.keys()),
        format_func=lambda x: IMPLANTATION_SITES[x]["name"],
    )
    site = IMPLANTATION_SITES[site_key]
    st.sidebar.info(site["description"])

    default_pO2 = float(site["pO2"])

    # --- AV-loop ---
    av_loop_flow = st.sidebar.checkbox(
        "Подключить прямой кровоток (AV-loop)",
        value=False,
        help="Хирургическое подключение артерии и вены к капсуле. Обеспечивает pO2 = 95.0 mmHg без лаг-фазы.",
    )
    if av_loop_flow:
        tau_blood = st.sidebar.slider(
            "Сдвиговое напряжение крови (tau_blood, Па)",
            0.1, 12.0, 5.0, step=0.1,
            help="Нормальное сдвиговое напряжение составляет 1.5 - 8.0 Па. Выход за рамки ускоряет тромбоз.",
        )
        anticoagulation = st.sidebar.checkbox(
            "Системная антикоагуляция",
            value=False,
            help="Использование антикоагулянтов замедляет тромбообразование в 10 раз.",
        )
        p_boundary = 95.0
    else:
        tau_blood = 5.0
        anticoagulation = False
        p_boundary = st.sidebar.slider(
            "Граничное давление O₂ (pO₂, mmHg)",
            2.0, 100.0, default_pO2, step=1.0,
        )

    # --- Fibrosis layer ---
    st.sidebar.header("🧱 Слой фиброза / Оболочка капсулы")

    coaxial_fiber = False
    if geometry_key == "cylindrical":
        coaxial_fiber = st.sidebar.checkbox(
            "Коаксиальная Core-Shell нить",
            value=True,
            help="Моделирует двухслойное волокно (клетки в ядре, защитный гель в оболочке).",
        )

    if coaxial_fiber:
        fib_label = "Толщина оболочки Shell (мкм)"
        fib_help = "Толщина защитного внешнего слоя гидрогеля без клеток."
    else:
        fib_label = "Толщина фиброзного слоя (мкм)"
        fib_help = "Толщина рубцовой ткани. 0 означает чистую капсулу."

    L_fibrosis_microns = st.sidebar.slider(
        fib_label,
        0.0, 300.0,
        value=float(st.session_state.get("l_fibrosis_microns", 0.0)),
        step=1.0,
        help=fib_help,
    )
    st.session_state["l_fibrosis_microns"] = L_fibrosis_microns

    D_fibrosis = 1.0e-5
    rho_mac = 0.0
    if L_fibrosis_microns > 0 or coaxial_fiber:
        D_fibrosis_multiplier = st.sidebar.slider(
            "Проницаемость внешней зоны (отн. воды)",
            0.1, 1.0, 0.3, step=0.1,
            help="0.3 означает, что диффузия во внешнем слое/рубце в 3+ раза хуже, чем в воде.",
        )
        D_fibrosis = D_fibrosis_multiplier * 3.0e-5
        rho_mac = st.sidebar.slider(
            "Плотность макрофагов в фиброзе (млн/мл)",
            0.0, 200.0, 50.0, step=10.0,
            help="Моделирует активное воспаление и поглощение кислорода макрофагами в рубце.",
        )

    # --- OGM ---
    st.sidebar.header("⚡ Активная оксигенация (OGM)")
    q_ogm = st.sidebar.slider(
        "Генерация O₂ OGM-генератором (mmHg/s)",
        0.0, 0.1, 0.0, step=0.01,
        help="Выделение кислорода частицами перекиси кальция внутри ядра капсулы.",
    )
    tethered_catalase = st.sidebar.checkbox(
        "Ковалентно связанная каталаза",
        value=False,
        help="Ковалентное сшивание каталазы снижает ее активность на 75%, но предотвращает ее вымывание под действием набухания (полураспад возрастает до 100 дней).",
    )

    # --- Cytokine traps ---
    st.sidebar.header("🛡️ Защита от цитокинов (Ловушки)")
    crispr_hypoimmune = st.sidebar.checkbox(
        "CRISPR Hypoimmune Клетки",
        value=False,
        help="Клетки с нокаутом MHC/HLA. Невидимы для IgG и защищены от антител.",
    )
    if crispr_hypoimmune:
        cd47_overexpression = st.sidebar.checkbox(
            "Гиперэкспрессия CD47 ('Don't eat me')",
            value=False,
            help="Защищает гипоиммунные клетки от лизиса NK-клетками.",
        )
        complement_protection = st.sidebar.checkbox(
            "CD55/CD59 защита комплемента",
            value=False,
            help="Защищает гипоиммунные клетки от активации системы комплемента и лизиса.",
        )
    else:
        cd47_overexpression = False
        complement_protection = False

    C_ext = st.sidebar.slider(
        "Внешний уровень цитокинов (ng/ml)",
        1.0, 50.0, 10.0, step=1.0,
        help="Концентрация провоспалительных цитокинов в ткани.",
    )
    k_bind_scav = st.sidebar.slider(
        "Константа ловушек (IL-1Ra, 1/(uM*s))",
        0.0, 2.0, 0.5, step=0.1,
        help="Скорость связывания цитокинов со-инкапсулированными антагонистами.",
    )

    # --- Physiology & Toxicology ---
    with st.sidebar.expander("🛡️ Физиология & Токсикология in vivo"):
        swelling_ratio = st.sidebar.slider(
            "Коэффициент набухания геля",
            1.0, 3.0, 1.0, step=0.1,
            help="Набухание геля увеличивает MWCO пор и разрушает иммунный барьер IgG.",
        )
        buffer_capacity_mM = st.sidebar.slider(
            "Буферная емкость геля (ммоль/л)",
            1.0, 50.0, 10.0, step=1.0,
            help="Нейтрализует щелочной pH-сдвиг от гидролиза пероксида кальция Ca(OH)2.",
        )
        catalase_activity_relative = st.sidebar.slider(
            "Активность каталазы (отн. нормы)",
            0.0, 2.0, 1.0, step=0.1,
            help="Каталаза расщепляет токсичную перекись H2O2 в кислород O2.",
        )
        catalase_half_life_days = st.sidebar.slider(
            "Полураспад каталазы (дней)",
            0.5, 10.0, 1.5, step=0.5,
            help="Скорость естественной деградации фермента каталазы in vivo.",
        )
        plga_acidification_factor = st.sidebar.slider(
            "Закисление среды PLGA",
            0.0, 1.0, 0.0, step=0.1,
            help="Моделирует кислоту от деградации PLGA. Снижает pH и усиливает FBR (макрофаги).",
        )
        t_days = st.sidebar.slider(
            "Срез времени для 1D анализа (дней)",
            0.0, 14.0, 0.0, step=0.5,
            help="Сутки симуляции для отображения распределения pH и ROS (перекиси H2O2).",
        )
        C_scav_0 = st.sidebar.slider(
            "Начальный пул ловушек C_scav (мкМ)",
            0.0, 20.0, 5.0, step=1.0,
            help="Начальный запас инкапсулированных ловушек цитокинов.",
        )
        k_deg_scav = st.sidebar.slider(
            "Деградация ловушек (1/день)",
            0.0, 0.1, 0.01, step=0.01,
            help="Скорость естественного вымывания/распада ловушек цитокинов.",
        )

    # --- Hydrogel ---
    st.sidebar.header("🧪 Свойства гидрогеля")
    hydrogel_key = st.sidebar.selectbox(
        "Материал капсулы:",
        options=list(HYDROGELS.keys()),
        format_func=lambda x: HYDROGELS[x]["name"],
    )
    D_coeff = HYDROGELS[hydrogel_key]["D"]
    st.sidebar.caption(f"Коэффициент диффузии O₂: {D_coeff:.1e} см²/с")

    phi_pfc = st.sidebar.slider(
        "Фторуглеродная эмульсия PFC (Доля)",
        0.0, 0.3, 0.0, step=0.05,
        help="Пассивный кислородный буфер (истощается за 2-4 дня).",
    )
    if phi_pfc > 0.0:
        pO2_pfc_saturation = st.sidebar.slider(
            "Насыщение PFC кислородом (pO₂, mmHg)",
            150.0, 760.0, 200.0, step=10.0,
            help="Давление насыщения PFC. Свыше 200 mmHg вызывает гипероксический оксидативный шок.",
        )
    else:
        pO2_pfc_saturation = 200.0

    E_0 = st.sidebar.slider(
        "Начальный модуль Юнга E₀ (кПа)",
        10.0, 200.0, 50.0, step=5.0,
        help="Начальная жесткость гидрогелевого каркаса. Снижается при набухании геля.",
    )

    # --- Cell properties ---
    st.sidebar.header("📏 Размеры и плотность клеток")
    R_outer_microns = st.sidebar.slider(
        "Внешний радиус / Полутолщина (R, мкм)",
        50, 600, 250, step=10,
        help="Для плоского листа это половина толщины. Для нити и сферы это их радиус.",
    )

    rho_million = st.sidebar.slider(
        "Плотность заселения клеток (млн/мл)",
        5, 250, 80, step=5,
        help="Количество клеток на миллилитр гидрогеля.",
    )

    turnover_rate = st.sidebar.slider(
        "Скорость апоптоза (Turnover, %/сутки)",
        0.0, 5.0, 1.5, step=0.1,
        help="Естественная скорость гибели бета-клеток в сутки. Вызывает постоянный DAMPs-воспалительный фон.",
    ) / 100.0

    with st.sidebar.expander("⚙️ Физиологические константы"):
        v_max_multiplier = st.slider("Множитель потребления O₂ (OCR)", 0.2, 3.0, 1.0, step=0.1)
        custom_V_max = V_MAX * v_max_multiplier

    # =====================================================================
    # SOLVER DISPATCH
    # =====================================================================
    if solver_type == "PyTorch (PINN)":
        from t1d_simulator.pinn_solver import solve_oxygen_profile_pinn, solve_cytokine_profile_pinn
        with st.spinner("Обучение физико-информированной нейросети (PINN) в реальном времени..."):
            res = solve_oxygen_profile_pinn(
                R_outer_microns=R_outer_microns,
                rho_million_per_ml=rho_million,
                p_boundary=p_boundary,
                D_coefficient=D_coeff,
                V_max_cell=custom_V_max,
                geometry=geometry_key,
                L_fibrosis_microns=L_fibrosis_microns,
                D_fibrosis=D_fibrosis,
                rho_mac_million_per_ml=rho_mac,
                q_ogm_mmHg_per_sec=q_ogm,
                catalase_activity_relative=catalase_activity_relative,
                catalase_half_life_days=catalase_half_life_days,
                buffer_capacity_mM=buffer_capacity_mM,
                swelling_ratio=swelling_ratio,
                plga_acidification_factor=plga_acidification_factor,
                t_days=t_days,
                tethered_catalase=tethered_catalase,
                E_0=E_0,
                species="Mouse",
                phi_pfc=phi_pfc,
                crispr_hypoimmune=crispr_hypoimmune,
                cd47_overexpression=cd47_overexpression,
                tau_blood=tau_blood,
                anticoagulation=anticoagulation,
                pO2_pfc_saturation=pO2_pfc_saturation,
                turnover_rate=turnover_rate,
                complement_protection=complement_protection,
            )
            res_cyt = solve_cytokine_profile_transient(
                R_outer_microns=R_outer_microns,
                C_ext=C_ext,
                D_cyt=1.0e-6,
                k_bind_scav=k_bind_scav,
                k_deg=0.01,
                C_scav_0=C_scav_0,
                k_deg_scav=k_deg_scav,
                swelling_ratio=swelling_ratio,
                days=14,
                shell_thickness_microns=L_fibrosis_microns,
                coaxial_active=coaxial_fiber,
                crispr_hypoimmune=crispr_hypoimmune,
                viable_fraction=res["viable_fraction"],
                turnover_rate=turnover_rate,
            )
    else:
        res = solve_oxygen_profile(
            R_outer_microns=R_outer_microns,
            rho_million_per_ml=rho_million,
            p_boundary=p_boundary,
            D_coefficient=D_coeff,
            V_max_cell=custom_V_max,
            geometry=geometry_key,
            L_fibrosis_microns=L_fibrosis_microns,
            D_fibrosis=D_fibrosis,
            rho_mac_million_per_ml=rho_mac,
            q_ogm_mmHg_per_sec=q_ogm,
            catalase_activity_relative=catalase_activity_relative,
            catalase_half_life_days=catalase_half_life_days,
            buffer_capacity_mM=buffer_capacity_mM,
            swelling_ratio=swelling_ratio,
            plga_acidification_factor=plga_acidification_factor,
            t_days=t_days,
            tethered_catalase=tethered_catalase,
            E_0=E_0,
            species="Mouse",
            phi_pfc=phi_pfc,
            crispr_hypoimmune=crispr_hypoimmune,
            cd47_overexpression=cd47_overexpression,
            tau_blood=tau_blood,
            anticoagulation=anticoagulation,
            pO2_pfc_saturation=pO2_pfc_saturation,
            turnover_rate=turnover_rate,
            complement_protection=complement_protection,
        )
        res_cyt = solve_cytokine_profile_transient(
            R_outer_microns=R_outer_microns,
            C_ext=C_ext,
            D_cyt=1.0e-6,
            k_bind_scav=k_bind_scav,
            k_deg=0.01,
            C_scav_0=C_scav_0,
            k_deg_scav=k_deg_scav,
            swelling_ratio=swelling_ratio,
            days=14,
            shell_thickness_microns=L_fibrosis_microns,
            coaxial_active=coaxial_fiber,
            crispr_hypoimmune=crispr_hypoimmune,
            viable_fraction=res["viable_fraction"],
            turnover_rate=turnover_rate,
        )

    # --- Cytokine compatibility wrapper ---
    avail_days_cyt = sorted(list(res_cyt["saved_profiles"].keys()))
    nearest_d_cyt = avail_days_cyt[np.argmin(np.abs(np.array(avail_days_cyt) - t_days))]
    cyt_snap = res_cyt["saved_profiles"][nearest_d_cyt]

    res_cyt_compat = {
        "z": cyt_snap["z"],
        "C": cyt_snap["C"],
        "protected_fraction": res_cyt["protected_fraction_over_time"][np.argmin(np.abs(res_cyt["t"] - t_days))],
    }
    res_cyt_original = res_cyt
    res_cyt = res_cyt_compat

    # --- Specific surface area ---
    R_cm = R_outer_microns * 1e-4
    if geometry_key == "planar":
        sav = 1.0 / R_cm
        sav_name = "1 / L"
    elif geometry_key == "cylindrical":
        sav = 2.0 / R_cm
        sav_name = "2 / R"
    else:
        sav = 3.0 / R_cm
        sav_name = "3 / R"

    # =====================================================================
    # METRICS
    # =====================================================================
    m1, m2, m3, m4, m5 = st.columns(5)
    vf = res["viable_fraction"]
    if vf > 95.0:
        vf_status = "Отлично"
    elif vf > 70.0:
        vf_status = "Гипоксический стресс"
    else:
        vf_status = "Массовый некроз"

    m1.metric(
        label="Выживаемость клеток",
        value=f"{vf:.1f}%",
        delta=vf_status,
        delta_color="normal" if vf > 70.0 else "inverse",
    )
    m2.metric(
        label="Секреция инсулина",
        value=f"{res['insulin_capacity']:.1f}%",
        delta="От нормы",
    )
    m3.metric(
        label="Мин. O₂ в центре",
        value=f"{res['min_pO2']:.2f} mmHg",
        delta="Порог нормы: 5.0",
    )
    m4.metric(
        label="Защита от цитокинов",
        value=f"{res_cyt['protected_fraction']:.1f}%",
        delta="Порог: < 1.0 ng/ml",
    )
    m5.metric(
        label="Удельная площадь SA/V",
        value=f"{sav:.1f} см⁻¹",
        delta=f"Формула: {sav_name}",
    )

    # =====================================================================
    # O2 PROFILE PLOT
    # =====================================================================
    st.subheader("📊 Распределение кислорода по сечению капсулы")
    st.markdown("График показывает парциальное давление кислорода ($pO_2$) по всему поперечному сечению капсулы (от $-R$ до $+R$).")

    z_coords = res["z"]
    pO2_profile = res["pO2"]
    x_full = np.concatenate((-z_coords[::-1], z_coords))
    pO2_full = np.concatenate((pO2_profile[::-1], pO2_profile))

    fig1 = go.Figure()
    fig1.add_trace(
        go.Scatter(
            x=x_full, y=pO2_full,
            name="Давление кислорода (pO₂)",
            line=dict(color="#60A5FA", width=4),
            fill="tozeroy",
            fillcolor="rgba(96, 165, 250, 0.1)",
        )
    )

    max_y_display = max(p_boundary * 1.1, 15.0)

    fig1.add_hrect(
        y0=0.0, y1=K_M,
        fillcolor="#EF4444", opacity=0.15, line_width=0,
        annotation_text="Зона некроза (гибель клеток) pO₂ < 0.5",
        annotation_position="top left",
        annotation_font=dict(color="#EF4444", size=11),
    )
    fig1.add_hrect(
        y0=K_M, y1=K_M_INSULIN,
        fillcolor="#F59E0B", opacity=0.12, line_width=0,
        annotation_text="Зона гипоксического стресса (нет секреции инсулина)",
        annotation_position="top left",
        annotation_font=dict(color="#F59E0B", size=11),
    )
    fig1.add_hrect(
        y0=K_M_INSULIN, y1=max_y_display,
        fillcolor="#10B981", opacity=0.08, line_width=0,
        annotation_text="Зона нормы",
        annotation_position="top left",
        annotation_font=dict(color="#10B981", size=11),
    )
    fig1.add_vline(x=0, line_width=1, line_dash="dash", line_color="rgba(255, 255, 255, 0.3)")

    fig1.update_layout(
        xaxis_title="Координата поперечного сечения (мкм)",
        yaxis_title="Парциальное давление O₂ (pO₂, mmHg)",
        template="plotly_dark",
        height=450,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=20, r=20, t=20, b=20),
        yaxis=dict(range=[0, max_y_display]),
    )
    st.plotly_chart(fig1, use_container_width=True)

    # =====================================================================
    # CYTOKINE PROFILE PLOT
    # =====================================================================
    st.subheader("🛡️ Распределение воспалительных цитокинов в капсуле")
    st.markdown("График показывает профиль концентрации цитокинов ($C$) по радиальному сечению гидрогеля от центра (0) до внешней границы ($R$).")

    z_cyt = res_cyt["z"]
    C_profile = res_cyt["C"]
    x_cyt_full = np.concatenate((-z_cyt[::-1], z_cyt))
    C_full = np.concatenate((C_profile[::-1], C_profile))

    fig_cyt = go.Figure()
    fig_cyt.add_trace(
        go.Scatter(
            x=x_cyt_full, y=C_full,
            name="Концентрация цитокинов (C)",
            line=dict(color="#F87171", width=4),
            fill="tozeroy",
            fillcolor="rgba(248, 113, 113, 0.1)",
        )
    )

    max_c_display = max(C_ext * 1.1, 5.0)

    fig_cyt.add_hrect(
        y0=1.0, y1=max_c_display,
        fillcolor="#EF4444", opacity=0.15, line_width=0,
        annotation_text="Зона цитотоксичности C >= 1.0 ng/ml",
        annotation_position="top left",
        annotation_font=dict(color="#EF4444", size=11),
    )
    fig_cyt.add_hrect(
        y0=0.0, y1=1.0,
        fillcolor="#10B981", opacity=0.1, line_width=0,
        annotation_text="Безопасная зона C < 1.0 ng/ml",
        annotation_position="top left",
        annotation_font=dict(color="#10B981", size=11),
    )
    fig_cyt.add_vline(x=0, line_width=1, line_dash="dash", line_color="rgba(255, 255, 255, 0.3)")

    fig_cyt.update_layout(
        xaxis_title="Координата поперечного сечения (мкм)",
        yaxis_title="Концентрация цитокинов (ng/ml)",
        template="plotly_dark",
        height=380,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=20, r=20, t=20, b=20),
        yaxis=dict(range=[0, max_c_display]),
    )
    st.plotly_chart(fig_cyt, use_container_width=True)

    # =====================================================================
    # pH & ROS PLOTS
    # =====================================================================
    if "pH" in res and "H2O2" in res:
        st.subheader("🧪 Распределение pH и окислительного стресса (H₂O₂)")
        st.markdown(f"Профили pH и концентрации H₂O₂ по сечению капсулы на **{t_days:.1f}** сутки.")

        col_ph, col_ros = st.columns(2)

        with col_ph:
            fig_ph = go.Figure()
            x_ph_full = np.concatenate((-res["z"][::-1], res["z"]))
            pH_full = np.concatenate((res["pH"][::-1], res["pH"]))
            fig_ph.add_trace(
                go.Scatter(
                    x=x_ph_full, y=pH_full,
                    name="pH",
                    line=dict(color="#10B981", width=3),
                    fill="tozeroy",
                    fillcolor="rgba(16, 185, 129, 0.05)",
                )
            )
            fig_ph.add_hrect(
                y0=6.8, y1=7.8,
                fillcolor="#10B981", opacity=0.1, line_width=0,
                annotation_text="Физиологический диапазон",
                annotation_position="top left",
                annotation_font=dict(color="#10B981", size=10),
            )
            fig_ph.update_layout(
                xaxis_title="Координата (мкм)",
                yaxis_title="pH",
                template="plotly_dark",
                height=350,
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                margin=dict(l=20, r=20, t=20, b=20),
                yaxis=dict(range=[5.0, 11.0]),
            )
            st.plotly_chart(fig_ph, use_container_width=True)

        with col_ros:
            fig_ros = go.Figure()
            x_ros_full = np.concatenate((-res["z"][::-1], res["z"]))
            ros_full = np.concatenate((res["H2O2"][::-1], res["H2O2"]))
            fig_ros.add_trace(
                go.Scatter(
                    x=x_ros_full, y=ros_full,
                    name="H₂O₂ (ROS)",
                    line=dict(color="#F59E0B", width=3),
                    fill="tozeroy",
                    fillcolor="rgba(245, 158, 11, 0.05)",
                )
            )
            max_ros_val = max(np.max(ros_full) * 1.1, 15.0)
            fig_ros.add_hrect(
                y0=10.0, y1=max_ros_val,
                fillcolor="#EF4444", opacity=0.12, line_width=0,
                annotation_text="Цитотоксический порог (10 мкМ)",
                annotation_position="top left",
                annotation_font=dict(color="#EF4444", size=10),
            )
            fig_ros.update_layout(
                xaxis_title="Координата (мкм)",
                yaxis_title="Концентрация H₂O₂ (мкМ)",
                template="plotly_dark",
                height=350,
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                margin=dict(l=20, r=20, t=20, b=20),
                yaxis=dict(range=[0, max_ros_val]),
            )
            st.plotly_chart(fig_ros, use_container_width=True)

        # Calcium & insulin leak
        st.subheader("🥛 Кальциевый гомеостаз и утечка инсулина (14-дневный таймлайн)")
        st.markdown("Моделирует накопление свободных ионов Ca²⁺ при растворении OGM-генератора на основе CaO₂ и последующий неконтролируемый выброс инсулина.")

        t_ca_days = np.linspace(0, 14, 100)
        t_ca_sec = t_ca_days * 24 * 3600
        ca_timeline = 1.2 + 2.0 * q_ogm * 1.34e-9 * 1e9 * t_ca_sec * 1e-3
        leak_timeline = np.minimum(1.0, np.maximum(0.0, ca_timeline - 1.2) / 4.0) * 100.0

        fig_ca = make_subplots(specs=[[{"secondary_y": True}]])
        fig_ca.add_trace(
            go.Scatter(x=t_ca_days, y=ca_timeline,
                       name="Концентрация Ca²⁺ (ммоль/л)",
                       line=dict(color="#3B82F6", width=3)),
            secondary_y=False,
        )
        fig_ca.add_trace(
            go.Scatter(x=t_ca_days, y=leak_timeline,
                       name="Утечка инсулина (%)",
                       line=dict(color="#EF4444", width=3, dash="dash")),
            secondary_y=True,
        )
        fig_ca.add_vline(
            x=t_days, line_width=2, line_dash="dot", line_color="#F59E0B",
            annotation_text=f"Выбранный срез: день {t_days}",
            annotation_position="top right",
        )
        fig_ca.update_layout(
            xaxis_title="Время (сутки)",
            template="plotly_dark",
            height=350,
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            margin=dict(l=20, r=20, t=20, b=20),
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        )
        fig_ca.update_yaxes(title_text="Ca²⁺ (ммоль/л)", secondary_y=False)
        fig_ca.update_yaxes(title_text="Утечка инсулина (%)", range=[0, 105], secondary_y=True)
        st.plotly_chart(fig_ca, use_container_width=True)

    # =====================================================================
    # SCVENGER DEPLETION TIMELINE
    # =====================================================================
    st.subheader("⏳ Динамика истощения ловушек цитокинов (14-дневный таймлайн)")
    st.markdown("Показывает изменение защищенной доли капсулы по мере связывания и вымывания ловушек цитокинов.")

    fig_timeline_scav = make_subplots(specs=[[{"secondary_y": True}]])
    fig_timeline_scav.add_trace(
        go.Scatter(
            x=res_cyt_original["t"],
            y=res_cyt_original["protected_fraction_over_time"],
            name="Защищенная фракция (%)",
            line=dict(color="#60A5FA", width=4),
            fill="tozeroy",
            fillcolor="rgba(96, 165, 250, 0.1)",
        ),
        secondary_y=False,
    )
    fig_timeline_scav.add_trace(
        go.Scatter(
            x=res_cyt_original["t"],
            y=res_cyt_original["C_ext_timeline"],
            name="Внешний уровень цитокинов C_ext (нг/мл)",
            line=dict(color="#F43F5E", width=3, dash="dash"),
        ),
        secondary_y=True,
    )
    fig_timeline_scav.add_vline(
        x=t_days, line_width=2, line_dash="dot", line_color="#F59E0B",
        annotation_text=f"Выбранный срез: день {t_days}",
        annotation_position="top right",
    )
    fig_timeline_scav.update_layout(
        xaxis_title="Время (сутки)",
        template="plotly_dark",
        height=370,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=20, r=20, t=20, b=20),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    fig_timeline_scav.update_yaxes(title_text="Защищенная доля ядра (%)", range=[0, 105], secondary_y=False)
    fig_timeline_scav.update_yaxes(title_text="C_ext (нг/мл)", secondary_y=True)
    st.plotly_chart(fig_timeline_scav, use_container_width=True)

    # =====================================================================
    # THROMBOSIS PLOT (AV-loop)
    # =====================================================================
    if av_loop_flow:
        st.subheader("🩸 Гемодинамическая окклюзия и тромбоз AV-шунта")
        st.markdown("Показывает прогрессирующий спад пропускной способности шунта и граничного давления кислорода из-за тромбоза и неоинтимальной гиперплазии гладких мыщц.")

        t_sweep = np.linspace(0, 30, 100)
        k_th_0 = 0.05
        k_th = k_th_0 * (1.0 + (max(0.0, 1.5 - tau_blood) / 0.5) ** 2 + (max(0.0, tau_blood - 8.0) / 2.0) ** 2)
        k_th_eff = k_th * 0.1 if anticoagulation else k_th
        k_hyp = 0.005
        k_occ = k_th_eff + k_hyp

        occlusion_timeline = (1.0 - np.exp(-k_occ * t_sweep)) * 100.0
        p_bound_timeline = 30.0 + (95.0 - 30.0) * np.exp(-k_occ * t_sweep)

        fig_thrombo = make_subplots(specs=[[{"secondary_y": True}]])
        fig_thrombo.add_trace(
            go.Scatter(x=t_sweep, y=occlusion_timeline,
                       name="Степень окклюзии шунта (%)",
                       line=dict(color="#EF4444", width=3)),
            secondary_y=False,
        )
        fig_thrombo.add_trace(
            go.Scatter(x=t_sweep, y=p_bound_timeline,
                       name="Граничное давление pO₂ (mmHg)",
                       line=dict(color="#3B82F6", width=3, dash="dash")),
            secondary_y=True,
        )
        fig_thrombo.add_vline(
            x=t_days, line_width=2, line_dash="dot", line_color="#F59E0B",
            annotation_text=f"День {t_days}",
            annotation_position="top right",
        )
        fig_thrombo.update_layout(
            xaxis_title="Время (сутки)",
            template="plotly_dark",
            height=350,
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            margin=dict(l=20, r=20, t=20, b=20),
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        )
        fig_thrombo.update_yaxes(title_text="Окклюзия шунта (%)", range=[0, 105], secondary_y=False)
        fig_thrombo.update_yaxes(title_text="Граничное pO₂ (mmHg)", range=[0, 105], secondary_y=True)
        st.plotly_chart(fig_thrombo, use_container_width=True)

    # =====================================================================
    # CELL DEPLETION
    # =====================================================================
    if turnover_rate > 0.0:
        st.subheader("💀 Истощение популяции клеток (Долгосрочный таймлайн)")
        st.markdown(f"Показывает неизбежное выгорание клеточного резервуара (жизнь β-клетки ~30-60 суток) из-за несоответствия скорости апоптоза ({turnover_rate*100.0:.1f}%/день) и митоза (~0.2%/день).")

        t_turn = np.linspace(0, 180, 200)
        depletion_timeline = np.exp((0.002 - turnover_rate) * t_turn) * 100.0

        fig_depletion = go.Figure()
        fig_depletion.add_trace(
            go.Scatter(
                x=t_turn, y=depletion_timeline,
                name="Выжившие клетки (%)",
                line=dict(color="#F59E0B", width=4),
                fill="tozeroy",
                fillcolor="rgba(245, 158, 11, 0.05)",
            )
        )
        fig_depletion.update_layout(
            xaxis_title="Время (сутки)",
            yaxis_title="Оставшаяся популяция клеток (%)",
            template="plotly_dark",
            height=350,
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            margin=dict(l=20, r=20, t=20, b=20),
            yaxis=dict(range=[0, 105]),
        )
        st.plotly_chart(fig_depletion, use_container_width=True)

    # =====================================================================
    # WARNINGS
    # =====================================================================
    survival_shear = res_cyt_original.get("survival_shear", 100.0)
    Fraction_death = 1.0 - survival_shear / 100.0
    if survival_shear < 90.0 or Fraction_death > 0.1:
        st.error(
            f"⚠️ **КРИТИЧЕСКИЙ КАСКАД DAMPs!** Выживаемость клеток при печати составила всего {survival_shear:.1f}% "
            f"(гибель клеток: {Fraction_death*100.0:.1f}%). Высвобождение DAMPs вызывает лавинообразное "
            f"воспаление ткани in vivo, увеличивая внешнюю концентрацию цитокинов C_ext экспоненциально."
        )

    # MWCO
    MWCO = 30.0 * (swelling_ratio ** 2)
    st.subheader("🛡️ Оценка иммуноизоляции и молекулярного сита (MWCO)")
    col_mwco_l, col_mwco_r = st.columns([1, 2])
    with col_mwco_l:
        st.metric(
            label="Эффективный предел MWCO",
            value=f"{MWCO:.1f} кДа",
            delta=f"Коэф. набухания: {swelling_ratio:.1f}",
        )
    with col_mwco_r:
        if MWCO >= 150.0:
            if crispr_hypoimmune:
                st.success("✅ **CRISPR Hypoimmune защита!** Гель сильно набух, но клетки генетически невидимы для антител (IgG).")
            else:
                st.error(
                    "⚠️ **КРИТИЧЕСКИЙ ПРОРЫВ IgG!** Гель сильно набух ($S_f \\ge 2.24$). "
                    "Поры расширились до размера антител IgG (150 кДа), что полностью уничтожает иммуноизоляцию "
                    "и вызывает атаку иммуноглобулинов на β-клетки."
                )
        elif MWCO > 50.0:
            st.warning(
                "⚠️ **Умеренное набухание.** Поры превысили стандартные 30 кДа. "
                "Диффузия цитокинов и белков средней массы ускорена, иммунный барьер ослаблен."
            )
        else:
            st.success(
                "✅ **Иммуноизоляция стабильна.** Предел MWCO удерживается в безопасном диапазоне (< 50 кДа). "
                "Крупные иммуноглобулины IgG надежно блокируются внешней оболочкой."
            )

        if res.get("rupture_risk", 0.0) > 50.0:
            st.error(
                f"⚠️ **РИСК МЕХАНИЧЕСКОГО РАЗРЫВА: {res['rupture_risk']:.1f}%!** "
                f"Вследствие сильного набухания геля эффективный модуль Юнга снизился до {E_0 * (swelling_ratio**-2):.2f} кПа. "
                f"Напряжение со стороны окружающих тканей превышает критический порог прочности, "
                f"что может привести к фрагментации капсулы in vivo и гибели β-клеток."
            )

        # CD47 & complement warnings
        if crispr_hypoimmune:
            if not cd47_overexpression:
                st.error("⚠️ **АТАКА NK-КИЛЛЕРОВ ('Missing Self')!** Полное отсутствие MHC-I на поверхности клеток активирует NK-клетки. Трансплантат будет уничтожен без гиперэкспрессии CD47.")
            if not complement_protection:
                st.error("⚠️ **ЛИЗИС СИСТЕМОЙ КОМПЛЕМЕНТА!** Без мембранных ингибиторов комплемента (CD55/CD59) белки крови пробивают поры в клетках за минуты.")
            if cd47_overexpression and complement_protection:
                st.success("✅ **Генетический щит активен!** Клетки защищены гиперэкспрессией CD47 от NK-клеток и CD55/CD59 от системы комплемента.")

        # Hyperoxia ROS
        if phi_pfc > 0.0 and pO2_pfc_saturation > 200.0:
            st.error(f"⚠️ **ГИПЕРОКСИЧЕСКИЙ ROS-ШОК!** Начальное насыщение PFC кислородом ({pO2_pfc_saturation:.1f} mmHg) превышает безопасные 200 mmHg. Бета-клетки страдают от мгновенного коллапса митохондрий под действием ROS.")

        # AV-loop warnings
        if av_loop_flow:
            if tau_blood < 1.5 or tau_blood > 8.0:
                st.error(f"⚠️ **ГЕМОДИНАМИЧЕСКИЙ КРИЗ (Тромбоз AV-петли)!** Сдвиговое напряжение ({tau_blood:.1f} Па) вне физиологической нормы [1.5, 8.0] Па. Поток крови замедляется тромбом.")
            st.warning(f"🩸 **Интимальная гиперплазия:** В месте анастомоза разрастается гладкая мускулатура (k = 0.005/день). Пропускная способность шунта снизится в долгосрочной перспективе.")

    # =====================================================================
    # COAXIAL 3D PRINTING PANEL
    # =====================================================================
    if geometry_key == "cylindrical" and coaxial_fiber:
        st.write("---")
        st.subheader("🖨️ Аналитическая панель коаксиальной 3D-биопечати")
        st.markdown("Моделирование напряжений сдвига в сопле коаксиального экструдера и оценка выживаемости клеток.")

        eta_core = 0.1
        eta_shell = 1.5
        Q_total = 2.0  # ml/min

        r_core_m = R_outer_microns * 1e-6
        r_shell_m = (R_outer_microns + L_fibrosis_microns) * 1e-6

        A_core = np.pi * (r_core_m ** 2)
        A_shell = np.pi * (r_shell_m ** 2 - r_core_m ** 2)
        Q_core = Q_total * (A_core / (A_core + A_shell))
        Q_shell = Q_total * (A_shell / (A_core + A_shell))

        Q_core_m3 = Q_core * 1.6667e-8
        Q_shell_m3 = Q_shell * 1.6667e-8

        tau_core_pa = (4.0 * eta_core * Q_core_m3) / (np.pi * (r_core_m ** 3) + 1e-30)
        tau_shell_pa = (4.0 * eta_shell * Q_shell_m3) / (np.pi * ((r_shell_m - r_core_m) ** 3) + 1e-30)

        tau_core_kpa = tau_core_pa / 1000.0
        tau_shell_kpa = tau_shell_pa / 1000.0

        Pr = 0.8 + 0.3 * np.tanh(tau_shell_kpa / 10.0)
        survival_shear = 100.0 * (1.0 - 0.6 * np.tanh(tau_core_kpa / 5.0))

        col_sh1, col_sh2, col_sh3, col_sh4 = st.columns(4)

        col_sh1.metric(
            label="Сдвиг в ядре (Core Shear)",
            value=f"{tau_core_kpa:.2f} кПа",
            delta="Безопасно < 5.0 кПа" if tau_core_kpa < 5.0 else "Критический сдвиг!",
            delta_color="normal" if tau_core_kpa < 5.0 else "inverse",
        )
        col_sh2.metric(
            label="Сдвиг в оболочке (Shell)",
            value=f"{tau_shell_kpa:.1f} кПа",
            delta="Сдвиговый щит активен",
        )
        col_sh3.metric(
            label="Выживаемость при печати",
            value=f"{survival_shear:.1f}%",
            delta="Механическая защита" if survival_shear > 90.0 else "Высокая смертность!",
        )
        col_sh4.metric(
            label="Индекс Pr (Printability)",
            value=f"{Pr:.2f}",
            delta="Норма: 0.9 - 1.1",
        )

        if tau_core_kpa < 5.0:
            st.success(
                "✅ **Shear Shielding работает!** Клетки находятся в низковязком ядре (Core) "
                "и защищены от высокого сдвигового напряжения внешней вязкой оболочкой (Shell)."
            )
        else:
            st.error(
                "❌ **Внимание! Сдвиговое напряжение в ядре превышает лимит!** "
                "Слишком тонкое внутреннее сопло или высокая плотность клеток вызывают гибель клеток при экструзии."
            )

    # =====================================================================
    # DENSITY SWEEP
    # =====================================================================
    st.write("---")
    st.subheader("📈 Зависимость выживаемости клеток от плотности упаковки")
    st.markdown("Этот график показывает процент выживаемости клеток в текущей капсуле при разных плотностях заселения.")

    sweep_res = run_density_sweep(
        R_outer_microns, p_boundary, D_coeff,
        geometry=geometry_key,
        L_fibrosis_microns=L_fibrosis_microns,
        D_fibrosis=D_fibrosis,
    )

    fig2 = go.Figure()
    fig2.add_trace(
        go.Scatter(
            x=sweep_res["densities"], y=sweep_res["viabilities"],
            name=f"Выбранная форма ({sav_name})",
            line=dict(color="#10B981", width=4),
        )
    )
    fig2.add_trace(
        go.Scatter(
            x=sweep_res["densities"], y=sweep_res["insulin_capacities"],
            name="Функциональность секреции инсулина (%)",
            line=dict(color="#3B82F6", width=2, dash="dash"),
        )
    )

    if geometry_key != "planar":
        sweep_planar = run_density_sweep(
            R_outer_microns, p_boundary, D_coeff,
            geometry="planar",
            L_fibrosis_microns=L_fibrosis_microns,
            D_fibrosis=D_fibrosis,
        )
        fig2.add_trace(
            go.Scatter(
                x=sweep_planar["densities"], y=sweep_planar["viabilities"],
                name="Контроль: Плоский лист (1/L)",
                line=dict(color="#EF4444", width=2, dash="dot"),
                opacity=0.6,
            )
        )

    fig2.add_vline(
        x=rho_million, line_width=2, line_dash="dot", line_color="#F59E0B",
        annotation_text=f"Выбранная плотность: {rho_million} млн/мл",
        annotation_position="top right",
    )
    fig2.update_layout(
        xaxis_title="Плотность заселения клеток (млн/мл)",
        yaxis_title="Процент (%)",
        template="plotly_dark",
        height=400,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=20, r=20, t=20, b=20),
        yaxis=dict(range=[0, 105]),
    )
    st.plotly_chart(fig2, use_container_width=True)

    # =====================================================================
    # SCIENTIFIC BLOCK
    # =====================================================================
    st.write("---")
    st.subheader("💡 Физико-математическое преимущество топологической оптимизации")
    st.markdown(f"""
    Удельная площадь поверхности капсулы к её объему ($SA/V$) определяет скорость притока кислорода на единицу клеточной массы.

    При выбранном радиусе $R = {R_outer_microns} \\ \\mu\\text{{m}}$:
    *   **Плоский лист (Планарная геометрия):** $SA/V = {1.0 / R_cm:.1f} \\ \\text{{см}}^{{-1}}$. Кислород поступает только с двух плоских сторон.
    *   **Цилиндр (Волокно / Микронить):** $SA/V = {2.0 / R_cm:.1f} \\ \\text{{см}}^{{-1}}$ (в **2 раза** выше). Кислород поступает радиально по всей окружности.
    *   **Сфера (Микрогранула):** $SA/V = {3.0 / R_cm:.1f} \\ \\text{{см}}^{{-1}}$ (в **3 раза** выше). Максимальная скорость диффузии.

    **Вывод для инженеров:**
    Переход от плоских пластин к 3D-биопечати в форме нитей или микросфер позволяет повысить оксигенацию клеток и предотвратить гипоксию.
    """)

    # =====================================================================
    # CAD EXPORT
    # =====================================================================
    st.write("---")
    st.subheader("💾 Экспорт CAD-модели для 3D-биопечати")
    st.markdown("Вы можете экспортировать геометрическую модель текущей капсулы в формат STL.")

    import mesh_generator

    col_cad_left, col_cad_right = st.columns([1, 2])
    with col_cad_left:
        if geometry_key == "planar":
            vertices, faces = mesh_generator.generate_box_mesh(L_microns=R_outer_microns)
            file_name = f"slab_L{R_outer_microns}um.stl"
            graft_name = "planar_slab"
        elif geometry_key == "cylindrical":
            vertices, faces = mesh_generator.generate_cylinder_mesh(R_microns=R_outer_microns)
            file_name = f"fiber_R{R_outer_microns}um.stl"
            graft_name = "cylindrical_fiber"
        else:
            vertices, faces = mesh_generator.generate_sphere_mesh(R_microns=R_outer_microns)
            file_name = f"microsphere_R{R_outer_microns}um.stl"
            graft_name = "spherical_microsphere"

        stl_data = mesh_generator.export_to_stl_ascii(vertices, faces, solid_name=graft_name)

        st.download_button(
            label="Скачать модель STL",
            data=stl_data,
            file_name=file_name,
            mime="application/sla",
        )

    with col_cad_right:
        st.info(f"**Экспортируемые параметры:** {file_name} | Вершин: {vertices.shape[0]} | Треугольников: {faces.shape[0]}. Единицы измерения: **микрометры**.")

    # Return result for external use
    return res
