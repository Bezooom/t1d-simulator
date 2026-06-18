import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as gr_obj
import mesh_generator

from simulator import (
    solve_oxygen_profile,
    solve_cytokine_profile,
    solve_cytokine_profile_transient,
    run_density_sweep,
    HYDROGELS,
    IMPLANTATION_SITES,
    K_M,
    K_M_INSULIN,
    V_MAX
)
from organoid_simulator import (
    simulate_organoid_population,
    simulate_organoid_oxygenation,
    simulate_organoid_insulin,
    calculate_immune_leak
)

# Инициализация состояния сессии для фиброза
if "l_fibrosis_microns" not in st.session_state:
    st.session_state["l_fibrosis_microns"] = 0.0

# Настройка страницы
st.set_page_config(
    page_title="In Silico Beta-Cell Encapsulation Twin",
    page_icon="🔮",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Кастомный CSS для дизайна
st.markdown("""
<style>
    .stMetric {
        background: rgba(255, 255, 255, 0.05);
        padding: 15px;
        border-radius: 12px;
        border: 1px solid rgba(255, 255, 255, 0.1);
    }
</style>
""", unsafe_allow_html=True)

st.title("🔮 Цифровой двойник капсулы β-клеток (Моделирование гипоксии)")
st.markdown("""
Симулятор решает нелинейное уравнение реакции-диффузии кислорода в различных типах геометрии (плоская пластина, цилиндрическая нить, сферическая гранула) и позволяет процедурно генерировать 3D-каркасы TPMS для оптимизации биопечати.
""")

# --- САЙДБАР: РЕЖИМ РАБОТЫ ---
st.sidebar.title("🛠️ Выберите режим")
app_mode = st.sidebar.radio(
    "Режим работы:",
    [
        "1D Симуляция диффузии O₂",
        "🔮 Генеративный 3D-дизайн (TPMS)",
        "🧪 ML-подбор антифиброзных покрытий (GNN)",
        "🩸 Неоваскуляризация (VEGF / Ангиогенез)",
        "🧫 Мини-органоиды (Фаза 10: Biomimesis)"
    ]
)

# ==============================================================================
# РЕЖИМ 1: 1D СИМУЛЯЦИЯ ДИФФУЗИИ O2
# ==============================================================================
if app_mode == "1D Симуляция диффузии O₂":
    st.sidebar.write("---")
    st.sidebar.header("🧠 Вычислительное ядро")
    solver_type = st.sidebar.radio(
        "Решатель ОДУ:",
        ["SciPy (solve_bvp)", "PyTorch (PINN)"],
        help="PINN обучается в реальном времени под физические параметры системы."
    )

    st.sidebar.header("📐 Выбор топологии (Геометрия)")
    geometry_key = st.sidebar.radio(
        "Форма имплантата:",
        ["planar", "cylindrical", "spherical"],
        format_func=lambda x: {
            "planar": "Плоская пластина (Slab / Лист)",
            "cylindrical": "Цилиндрическая нить (Fiber / Волокно)",
            "spherical": "Сферическая микрокапсула (Microsphere)"
        }[x]
    )

    st.sidebar.header("📍 Место имплантации (Диффузия O₂)")
    site_key = st.sidebar.selectbox(
        "Выберите зону пересадки:",
        options=list(IMPLANTATION_SITES.keys()),
        format_func=lambda x: IMPLANTATION_SITES[x]["name"]
    )
    site = IMPLANTATION_SITES[site_key]
    st.sidebar.info(site["description"])

    default_pO2 = float(site["pO2"])
    
    av_loop_flow = st.sidebar.checkbox(
        "Подключить прямой кровоток (AV-loop)",
        value=False,
        help="Хирургическое подключение артерии и вены к капсуле. Обеспечивает pO2 = 95.0 mmHg без лаг-фазы."
    )
    if av_loop_flow:
        tau_blood = st.sidebar.slider(
            "Сдвиговое напряжение крови (tau_blood, Па)",
            0.1, 12.0, 5.0,
            step=0.1,
            help="Нормальное сдвиговое напряжение составляет 1.5 - 8.0 Па. Выход за рамки ускоряет тромбоз."
        )
        anticoagulation = st.sidebar.checkbox(
            "Системная антикоагуляция",
            value=False,
            help="Использование антикоагулянтов замедляет тромбообразование в 10 раз."
        )
        p_boundary = 95.0
    else:
        tau_blood = 5.0
        anticoagulation = False
        p_boundary = st.sidebar.slider(
            "Граничное давление O₂ (pO₂, mmHg)",
            2.0, 100.0, default_pO2,
            step=1.0
        )

    st.sidebar.header("🧱 Слой фиброза / Оболочка капсулы")
    
    coaxial_fiber = False
    if geometry_key == "cylindrical":
        coaxial_fiber = st.sidebar.checkbox(
            "Коаксиальная Core-Shell нить",
            value=True,
            help="Моделирует двухслойное волокно (клетки в ядре, защитный гель в оболочке)."
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
        value=float(st.session_state["l_fibrosis_microns"]),
        step=1.0,
        help=fib_help
    )
    st.session_state["l_fibrosis_microns"] = L_fibrosis_microns

    D_fibrosis = 1.0e-5
    rho_mac = 0.0
    if L_fibrosis_microns > 0 or coaxial_fiber:
        D_fibrosis_multiplier = st.sidebar.slider(
            "Проницаемость внешней зоны (отн. воды)",
            0.1, 1.0, 0.3,
            step=0.1,
            help="0.3 означает, что диффузия во внешнем слое/рубце в 3+ раза хуже, чем в воде."
        )
        D_fibrosis = D_fibrosis_multiplier * 3.0e-5
        
        rho_mac = st.sidebar.slider(
            "Плотность макрофагов в фиброзе (млн/мл)",
            0.0, 200.0, 50.0,
            step=10.0,
            help="Моделирует активное воспаление и поглощение кислорода макрофагами в рубце."
        )

    st.sidebar.header("⚡ Активная оксигенация (OGM)")
    q_ogm = st.sidebar.slider(
        "Генерация O₂ OGM-генератором (mmHg/s)",
        0.0, 0.1, 0.0,
        step=0.01,
        help="Выделение кислорода частицами перекиси кальция внутри ядра капсулы."
    )
    tethered_catalase = st.sidebar.checkbox(
        "Ковалентно связанная каталаза",
        value=False,
        help="Ковалентное сшивание каталазы снижает ее активность на 75%, но предотвращает ее вымывание под действием набухания (полураспад возрастает до 100 дней)."
    )

    st.sidebar.header("🛡️ Защита от цитокинов (Ловушки)")
    crispr_hypoimmune = st.sidebar.checkbox(
        "CRISPR Hypoimmune Клетки",
        value=False,
        help="Клетки с нокаутом MHC/HLA. Невидимы для IgG и защищены от антител."
    )
    if crispr_hypoimmune:
        cd47_overexpression = st.sidebar.checkbox(
            "Гиперэкспрессия CD47 ('Don't eat me')",
            value=False,
            help="Защищает гипоиммунные клетки от лизиса NK-клетками."
        )
        complement_protection = st.sidebar.checkbox(
            "CD55/CD59 защита комплемента",
            value=False,
            help="Защищает гипоиммунные клетки от активации системы комплемента и лизиса."
        )
    else:
        cd47_overexpression = False
        complement_protection = False
        
    C_ext = st.sidebar.slider(
        "Внешний уровень цитокинов (ng/ml)",
        1.0, 50.0, 10.0,
        step=1.0,
        help="Концентрация провоспалительных цитокинов в ткани."
    )
    k_bind_scav = st.sidebar.slider(
        "Константа ловушек (IL-1Ra, 1/(uM*s))",
        0.0, 2.0, 0.5,
        step=0.1,
        help="Скорость связывания цитокинов со-инкапсулированными антагонистами."
    )

    with st.sidebar.expander("🛡️ Физиология & Токсикология in vivo"):
        swelling_ratio = st.slider(
            "Коэффициент набухания геля",
            1.0, 3.0, 1.0,
            step=0.1,
            help="Набухание геля увеличивает MWCO пор и разрушает иммунный барьер IgG."
        )
        buffer_capacity_mM = st.slider(
            "Буферная емкость геля (ммоль/л)",
            1.0, 50.0, 10.0,
            step=1.0,
            help="Нейтрализует щелочной pH-сдвиг от гидролиза пероксида кальция Ca(OH)2."
        )
        catalase_activity_relative = st.slider(
            "Активность каталазы (отн. нормы)",
            0.0, 2.0, 1.0,
            step=0.1,
            help="Каталаза расщепляет токсичную перекись H2O2 в кислород O2."
        )
        catalase_half_life_days = st.slider(
            "Полураспад каталазы (дней)",
            0.5, 10.0, 1.5,
            step=0.5,
            help="Скорость естественной деградации фермента каталазы in vivo."
        )
        plga_acidification_factor = st.slider(
            "Закисление среды PLGA",
            0.0, 1.0, 0.0,
            step=0.1,
            help="Моделирует кислоту от деградации PLGA. Снижает pH и усиливает FBR (макрофаги)."
        )
        t_days = st.slider(
            "Срез времени для 1D анализа (дней)",
            0.0, 14.0, 0.0,
            step=0.5,
            help="Сутки симуляции для отображения распределения pH и ROS (перекиси H2O2)."
        )
        C_scav_0 = st.slider(
            "Начальный пул ловушек C_scav (мкМ)",
            0.0, 20.0, 5.0,
            step=1.0,
            help="Начальный запас инкапсулированных ловушек цитокинов."
        )
        k_deg_scav = st.slider(
            "Деградация ловушек (1/день)",
            0.0, 0.1, 0.01,
            step=0.01,
            help="Скорость естественного вымывания/распада ловушек цитокинов."
        )

    st.sidebar.header("🧪 Свойства гидрогеля")
    hydrogel_key = st.sidebar.selectbox(
        "Материал капсулы:",
        options=list(HYDROGELS.keys()),
        format_func=lambda x: HYDROGELS[x]["name"]
    )
    D_coeff = HYDROGELS[hydrogel_key]["D"]
    st.sidebar.caption(f"Коэффициент диффузии O₂: {D_coeff:.1e} см²/с")
    
    phi_pfc = st.sidebar.slider(
        "Фторуглеродная эмульсия PFC (Доля)",
        0.0, 0.3, 0.0,
        step=0.05,
        help="Пассивный кислородный буфер (истощается за 2-4 дня)."
    )
    if phi_pfc > 0.0:
        pO2_pfc_saturation = st.sidebar.slider(
            "Насыщение PFC кислородом (pO₂, mmHg)",
            150.0, 760.0, 200.0,
            step=10.0,
            help="Давление насыщения PFC. Свыше 200 mmHg вызывает гипероксический оксидативный шок."
        )
    else:
        pO2_pfc_saturation = 200.0
        
    E_0 = st.sidebar.slider(
        "Начальный модуль Юнга E₀ (кПа)",
        10.0, 200.0, 50.0,
        step=5.0,
        help="Начальная жесткость гидрогелевого каркаса. Снижается при набухании геля."
    )

    st.sidebar.header("📏 Размеры и плотность клеток")
    R_outer_microns = st.sidebar.slider(
        "Внешний радиус / Полутолщина (R, мкм)",
        50, 600, 250,
        step=10,
        help="Для плоского листа это половина толщины. Для нити и сферы это их радиус."
    )

    rho_million = st.sidebar.slider(
        "Плотность заселения клеток (млн/мл)",
        5, 250, 80,
        step=5,
        help="Количество клеток на миллилитр гидрогеля."
    )
    
    turnover_rate = st.sidebar.slider(
        "Скорость апоптоза (Turnover, %/сутки)",
        0.0, 5.0, 1.5,
        step=0.1,
        help="Естественная скорость гибели бета-клеток в сутки. Вызывает постоянный DAMPs-воспалительный фон."
    ) / 100.0

    with st.sidebar.expander("⚙️ Физиологические константы"):
        v_max_multiplier = st.slider("Множитель потребления O₂ (OCR)", 0.2, 3.0, 1.0, step=0.1)
        custom_V_max = V_MAX * v_max_multiplier

    # --- РАСЧЕТ ---
    if solver_type == "PyTorch (PINN)":
        from pinn_solver import solve_oxygen_profile_pinn, solve_cytokine_profile_pinn
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
                complement_protection=complement_protection
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
                turnover_rate=turnover_rate
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
            complement_protection=complement_protection
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
            turnover_rate=turnover_rate
        )

    # Извлекаем срез цитокинов для выбранного t_days
    avail_days_cyt = sorted(list(res_cyt["saved_profiles"].keys()))
    nearest_d_cyt = avail_days_cyt[np.argmin(np.abs(np.array(avail_days_cyt) - t_days))]
    cyt_snap = res_cyt["saved_profiles"][nearest_d_cyt]
    
    # Модифицируем res_cyt для совместимости с кодом отрисовки
    res_cyt_compat = {
        "z": cyt_snap["z"],
        "C": cyt_snap["C"],
        "protected_fraction": res_cyt["protected_fraction_over_time"][np.argmin(np.abs(res_cyt["t"] - t_days))]
    }
    res_cyt_original = res_cyt
    res_cyt = res_cyt_compat

    # Удельная площадь
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

    # Вывод метрик
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
        delta_color="normal" if vf > 70.0 else "inverse"
    )

    m2.metric(
        label="Секреция инсулина",
        value=f"{res['insulin_capacity']:.1f}%",
        delta="От нормы"
    )

    m3.metric(
        label="Мин. O₂ в центре",
        value=f"{res['min_pO2']:.2f} mmHg",
        delta="Порог нормы: 5.0"
    )

    m4.metric(
        label="Защита от цитокинов",
        value=f"{res_cyt['protected_fraction']:.1f}%",
        delta="Порог: < 1.0 ng/ml"
    )

    m5.metric(
        label="Удельная площадь SA/V",
        value=f"{sav:.1f} см⁻¹",
        delta=f"Формула: {sav_name}"
    )

    # График распределения
    st.subheader("📊 Распределение кислорода по сечению капсулы")
    st.markdown("График показывает парциальное давление кислорода ($pO_2$) по всему поперечному сечению капсулы (от $-R$ до $+R$).")

    z_coords = res["z"]
    pO2_profile = res["pO2"]
    x_full = np.concatenate((-z_coords[::-1], z_coords))
    pO2_full = np.concatenate((pO2_profile[::-1], pO2_profile))

    fig1 = gr_obj.Figure()
    fig1.add_trace(
        gr_obj.Scatter(
            x=x_full, 
            y=pO2_full, 
            name="Давление кислорода (pO₂)", 
            line=dict(color="#60A5FA", width=4),
            fill="tozeroy",
            fillcolor="rgba(96, 165, 250, 0.1)"
        )
    )

    max_y_display = max(p_boundary * 1.1, 15.0)

    fig1.add_hrect(
        y0=0.0, y1=K_M, 
        fillcolor="#EF4444", opacity=0.15, line_width=0,
        annotation_text="Зона некроза (гибель клеток) pO₂ < 0.5", 
        annotation_position="top left",
        annotation_font=dict(color="#EF4444", size=11)
    )

    fig1.add_hrect(
        y0=K_M, y1=K_M_INSULIN, 
        fillcolor="#F59E0B", opacity=0.12, line_width=0,
        annotation_text="Зона гипоксического стресса (нет секреции инсулина)", 
        annotation_position="top left",
        annotation_font=dict(color="#F59E0B", size=11)
    )

    fig1.add_hrect(
        y0=K_M_INSULIN, y1=max_y_display, 
        fillcolor="#10B981", opacity=0.08, line_width=0,
        annotation_text="Зона нормы", 
        annotation_position="top left",
        annotation_font=dict(color="#10B981", size=11)
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
        yaxis=dict(range=[0, max_y_display])
    )
    st.plotly_chart(fig1, use_container_width=True)

    # График распределения цитокинов
    st.subheader("🛡️ Распределение воспалительных цитокинов в капсуле")
    st.markdown("График показывает профиль концентрации цитокинов ($C$) по радиальному сечению гидрогеля от центра (0) до внешней границы ($R$).")
    
    z_cyt = res_cyt["z"]
    C_profile = res_cyt["C"]
    
    x_cyt_full = np.concatenate((-z_cyt[::-1], z_cyt))
    C_full = np.concatenate((C_profile[::-1], C_profile))
    
    fig_cyt = gr_obj.Figure()
    fig_cyt.add_trace(
        gr_obj.Scatter(
            x=x_cyt_full,
            y=C_full,
            name="Концентрация цитокинов (C)",
            line=dict(color="#F87171", width=4),
            fill="tozeroy",
            fillcolor="rgba(248, 113, 113, 0.1)"
        )
    )
    
    max_c_display = max(C_ext * 1.1, 5.0)
    
    fig_cyt.add_hrect(
        y0=1.0, y1=max_c_display,
        fillcolor="#EF4444", opacity=0.15, line_width=0,
        annotation_text="Зона цитотоксичности C >= 1.0 ng/ml",
        annotation_position="top left",
        annotation_font=dict(color="#EF4444", size=11)
    )
    
    fig_cyt.add_hrect(
        y0=0.0, y1=1.0,
        fillcolor="#10B981", opacity=0.1, line_width=0,
        annotation_text="Безопасная зона C < 1.0 ng/ml",
        annotation_position="top left",
        annotation_font=dict(color="#10B981", size=11)
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
        yaxis=dict(range=[0, max_c_display])
    )
    st.plotly_chart(fig_cyt, use_container_width=True)

    # График распределения pH и ROS (H2O2)
    if "pH" in res and "H2O2" in res:
        st.subheader("🧪 Распределение pH и окислительного стресса (H₂O₂)")
        st.markdown(f"Профили pH и концентрации H₂O₂ по сечению капсулы на **{t_days:.1f}** сутки.")
        
        col_ph, col_ros = st.columns(2)
        
        with col_ph:
            fig_ph = gr_obj.Figure()
            x_ph_full = np.concatenate((-res["z"][::-1], res["z"]))
            pH_full = np.concatenate((res["pH"][::-1], res["pH"]))
            fig_ph.add_trace(
                gr_obj.Scatter(
                    x=x_ph_full,
                    y=pH_full,
                    name="pH",
                    line=dict(color="#10B981", width=3),
                    fill="tozeroy",
                    fillcolor="rgba(16, 185, 129, 0.05)"
                )
            )
            fig_ph.add_hrect(
                y0=6.8, y1=7.8,
                fillcolor="#10B981", opacity=0.1, line_width=0,
                annotation_text="Физиологический диапазон",
                annotation_position="top left",
                annotation_font=dict(color="#10B981", size=10)
            )
            fig_ph.update_layout(
                xaxis_title="Координата (мкм)",
                yaxis_title="pH",
                template="plotly_dark",
                height=350,
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                margin=dict(l=20, r=20, t=20, b=20),
                yaxis=dict(range=[5.0, 11.0])
            )
            st.plotly_chart(fig_ph, use_container_width=True)
            
        with col_ros:
            fig_ros = gr_obj.Figure()
            x_ros_full = np.concatenate((-res["z"][::-1], res["z"]))
            ros_full = np.concatenate((res["H2O2"][::-1], res["H2O2"]))
            fig_ros.add_trace(
                gr_obj.Scatter(
                    x=x_ros_full,
                    y=ros_full,
                    name="H₂O₂ (ROS)",
                    line=dict(color="#F59E0B", width=3),
                    fill="tozeroy",
                    fillcolor="rgba(245, 158, 11, 0.05)"
                )
            )
            max_ros_val = max(np.max(ros_full) * 1.1, 15.0)
            fig_ros.add_hrect(
                y0=10.0, y1=max_ros_val,
                fillcolor="#EF4444", opacity=0.12, line_width=0,
                annotation_text="Цитотоксический порог (10 мкМ)",
                annotation_position="top left",
                annotation_font=dict(color="#EF4444", size=10)
            )
            fig_ros.update_layout(
                xaxis_title="Координата (мкм)",
                yaxis_title="Концентрация H₂O₂ (мкМ)",
                template="plotly_dark",
                height=350,
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                margin=dict(l=20, r=20, t=20, b=20),
                yaxis=dict(range=[0, max_ros_val])
            )
            st.plotly_chart(fig_ros, use_container_width=True)

        # Накопление кальция и утечка инсулина
        st.subheader("🥛 Кальциевый гомеостаз и утечка инсулина (14-дневный таймлайн)")
        st.markdown("Моделирует накопление свободных ионов Ca²⁺ при растворении OGM-генератора на основе CaO₂ и последующий неконтролируемый выброс инсулина.")
        
        t_ca_days = np.linspace(0, 14, 100)
        t_ca_sec = t_ca_days * 24 * 3600
        ca_timeline = 1.2 + 2.0 * q_ogm * 1.34e-9 * 1e9 * t_ca_sec * 1e-3
        leak_timeline = np.minimum(1.0, np.maximum(0.0, ca_timeline - 1.2) / 4.0) * 100.0 # в процентах
        
        from plotly.subplots import make_subplots
        fig_ca = make_subplots(specs=[[{"secondary_y": True}]])
        
        fig_ca.add_trace(
            gr_obj.Scatter(
                x=t_ca_days,
                y=ca_timeline,
                name="Концентрация Ca²⁺ (ммоль/л)",
                line=dict(color="#3B82F6", width=3)
            ),
            secondary_y=False
        )
        
        fig_ca.add_trace(
            gr_obj.Scatter(
                x=t_ca_days,
                y=leak_timeline,
                name="Утечка инсулина (%)",
                line=dict(color="#EF4444", width=3, dash="dash")
            ),
            secondary_y=True
        )
        
        fig_ca.add_vline(
            x=t_days,
            line_width=2,
            line_dash="dot",
            line_color="#F59E0B",
            annotation_text=f"Выбранный срез: день {t_days}",
            annotation_position="top right"
        )
        
        fig_ca.update_layout(
            xaxis_title="Время (сутки)",
            template="plotly_dark",
            height=350,
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            margin=dict(l=20, r=20, t=20, b=20),
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
        )
        
        fig_ca.update_yaxes(title_text="Ca²⁺ (ммоль/л)", secondary_y=False)
        fig_ca.update_yaxes(title_text="Утечка инсулина (%)", range=[0, 105], secondary_y=True)
        
        st.plotly_chart(fig_ca, use_container_width=True)

    # График динамики истощения ловушек за 14 дней
    st.subheader("⏳ Динамика истощения ловушек цитокинов (14-дневный таймлайн)")
    st.markdown("Показывает изменение защищенной доли капсулы по мере связывания и вымывания ловушек цитокинов.")
    
    from plotly.subplots import make_subplots
    fig_timeline_scav = make_subplots(specs=[[{"secondary_y": True}]])
    
    fig_timeline_scav.add_trace(
        gr_obj.Scatter(
            x=res_cyt_original["t"],
            y=res_cyt_original["protected_fraction_over_time"],
            name="Защищенная фракция (%)",
            line=dict(color="#60A5FA", width=4),
            fill="tozeroy",
            fillcolor="rgba(96, 165, 250, 0.1)"
        ),
        secondary_y=False
    )
    
    fig_timeline_scav.add_trace(
        gr_obj.Scatter(
            x=res_cyt_original["t"],
            y=res_cyt_original["C_ext_timeline"],
            name="Внешний уровень цитокинов C_ext (нг/мл)",
            line=dict(color="#F43F5E", width=3, dash="dash")
        ),
        secondary_y=True
    )
    
    fig_timeline_scav.add_vline(
        x=t_days,
        line_width=2,
        line_dash="dot",
        line_color="#F59E0B",
        annotation_text=f"Выбранный срез: день {t_days}",
        annotation_position="top right"
    )
    
    fig_timeline_scav.update_layout(
        xaxis_title="Время (сутки)",
        template="plotly_dark",
        height=370,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=20, r=20, t=20, b=20),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )
    
    fig_timeline_scav.update_yaxes(title_text="Защищенная доля ядра (%)", range=[0, 105], secondary_y=False)
    fig_timeline_scav.update_yaxes(title_text="C_ext (нг/мл)", secondary_y=True)
    
    st.plotly_chart(fig_timeline_scav, use_container_width=True)

    # Phase 9: AV-loop Thrombosis dynamic plot
    if av_loop_flow:
        st.subheader("🩸 Гемодинамическая окклюзия и тромбоз AV-шунта")
        st.markdown("Показывает прогрессирующий спад пропускной способности шунта и граничного давления кислорода из-за тромбоза и неоинтимальной гиперплазии гладких мыщц.")
        
        t_sweep = np.linspace(0, 30, 100)
        k_th_0 = 0.05
        k_th = k_th_0 * (1.0 + (max(0.0, 1.5 - tau_blood) / 0.5)**2 + (max(0.0, tau_blood - 8.0) / 2.0)**2)
        k_th_eff = k_th * 0.1 if anticoagulation else k_th
        k_hyp = 0.005
        k_occ = k_th_eff + k_hyp
        
        occlusion_timeline = (1.0 - np.exp(-k_occ * t_sweep)) * 100.0
        p_bound_timeline = 30.0 + (95.0 - 30.0) * np.exp(-k_occ * t_sweep)
        
        fig_thrombo = make_subplots(specs=[[{"secondary_y": True}]])
        fig_thrombo.add_trace(
            gr_obj.Scatter(
                x=t_sweep,
                y=occlusion_timeline,
                name="Степень окклюзии шунта (%)",
                line=dict(color="#EF4444", width=3)
            ),
            secondary_y=False
        )
        fig_thrombo.add_trace(
            gr_obj.Scatter(
                x=t_sweep,
                y=p_bound_timeline,
                name="Граничное давление pO₂ (mmHg)",
                line=dict(color="#3B82F6", width=3, dash="dash")
            ),
            secondary_y=True
        )
        fig_thrombo.add_vline(
            x=t_days,
            line_width=2,
            line_dash="dot",
            line_color="#F59E0B",
            annotation_text=f"День {t_days}",
            annotation_position="top right"
        )
        fig_thrombo.update_layout(
            xaxis_title="Время (сутки)",
            template="plotly_dark",
            height=350,
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            margin=dict(l=20, r=20, t=20, b=20),
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
        )
        fig_thrombo.update_yaxes(title_text="Окклюзия шунта (%)", range=[0, 105], secondary_y=False)
        fig_thrombo.update_yaxes(title_text="Граничное pO₂ (mmHg)", range=[0, 105], secondary_y=True)
        st.plotly_chart(fig_thrombo, use_container_width=True)

    # Phase 9: Cell depletion curve due to turnover
    if turnover_rate > 0.0:
        st.subheader("💀 Истощение популяции клеток (Долгосрочный таймлайн)")
        st.markdown(f"Показывает неизбежное выгорание клеточного резервуара (жизнь β-клетки ~30-60 суток) из-за несоответствия скорости апоптоза ({turnover_rate*100.0:.1f}%/день) и митоза (~0.2%/день).")
        
        t_turn = np.linspace(0, 180, 200)
        depletion_timeline = np.exp((0.002 - turnover_rate) * t_turn) * 100.0
        
        fig_depletion = gr_obj.Figure()
        fig_depletion.add_trace(
            gr_obj.Scatter(
                x=t_turn,
                y=depletion_timeline,
                name="Выжившие клетки (%)",
                line=dict(color="#F59E0B", width=4),
                fill="tozeroy",
                fillcolor="rgba(245, 158, 11, 0.05)"
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
            yaxis=dict(range=[0, 105])
        )
        st.plotly_chart(fig_depletion, use_container_width=True)
    
    survival_shear = res_cyt_original.get("survival_shear", 100.0)
    Fraction_death = 1.0 - survival_shear / 100.0
    if survival_shear < 90.0 or Fraction_death > 0.1:
        st.error(
            f"⚠️ **КРИТИЧЕСКИЙ КАСКАД DAMPs!** Выживаемость клеток при печати составила всего {survival_shear:.1f}% "
            f"(гибель клеток: {Fraction_death*100.0:.1f}%). Высвобождение DAMPs вызывает лавинообразное "
            f"воспаление ткани in vivo, увеличивая внешнюю концентрацию цитокинов C_ext экспоненциально."
        )

    # Индикатор риска иммунного прорыва IgG
    MWCO = 30.0 * (swelling_ratio ** 2)
    st.subheader("🛡️ Оценка иммуноизоляции и молекулярного сита (MWCO)")
    col_mwco_l, col_mwco_r = st.columns([1, 2])
    with col_mwco_l:
        st.metric(
            label="Эффективный предел MWCO",
            value=f"{MWCO:.1f} кДа",
            delta=f"Коэф. набухания: {swelling_ratio:.1f}"
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
            
        # Phase 9: CD47 & complement warnings
        if crispr_hypoimmune:
            if not cd47_overexpression:
                st.error("⚠️ **АТАКА NK-КИЛЛЕРОВ ('Missing Self')!** Полное отсутствие MHC-I на поверхности клеток активирует NK-клетки. Трансплантат будет уничтожен без гиперэкспрессии CD47.")
            if not complement_protection:
                st.error("⚠️ **ЛИЗИС СИСТЕМОЙ КОМПЛЕМЕНТА!** Без мембранных ингибиторов комплемента (CD55/CD59) белки крови пробивают поры в клетках за минуты.")
            if cd47_overexpression and complement_protection:
                st.success("✅ **Генетический щит активен!** Клетки защищены гиперэкспрессией CD47 от NK-клеток и CD55/CD59 от системы комплемента.")

        # Phase 9: Hyperoxia ROS stress warning
        if phi_pfc > 0.0 and pO2_pfc_saturation > 200.0:
            st.error(f"⚠️ **ГИПЕРОКСИЧЕСКИЙ ROS-ШОК!** Начальное насыщение PFC кислородом ({pO2_pfc_saturation:.1f} mmHg) превышает безопасные 200 mmHg. Бета-клетки страдают от мгновенного коллапса митохондрий под действием ROS.")

        # Phase 9: AV-loop thrombosis & hyperplasia warnings
        if av_loop_flow:
            if tau_blood < 1.5 or tau_blood > 8.0:
                st.error(f"⚠️ **ГЕМОДИНАМИЧЕСКИЙ КРИЗ (Тромбоз AV-петли)!** Сдвиговое напряжение ({tau_blood:.1f} Па) вне физиологической нормы [1.5, 8.0] Па. Поток крови замедляется тромбом.")
            st.warning(f"🩸 **Интимальная гиперплазия:** В месте анастомоза разрастается гладкая мускулатура (k = 0.005/день). Пропускная способность шунта снизится в долгосрочной перспективе.")

    # Панель коаксиальной 3D-печати
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
            delta_color="normal" if tau_core_kpa < 5.0 else "inverse"
        )
        
        col_sh2.metric(
            label="Сдвиг в оболочке (Shell)",
            value=f"{tau_shell_kpa:.1f} кПа",
            delta="Сдвиговый щит активен"
        )
        
        col_sh3.metric(
            label="Выживаемость при печати",
            value=f"{survival_shear:.1f}%",
            delta="Механическая защита" if survival_shear > 90.0 else "Высокая смертность!"
        )
        
        col_sh4.metric(
            label="Индекс Pr (Printability)",
            value=f"{Pr:.2f}",
            delta="Норма: 0.9 - 1.1"
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

    # График зависимости от плотности
    st.write("---")
    st.subheader("📈 Зависимость выживаемости клеток от плотности упаковки")
    st.markdown("Этот график показывает процент выживаемости клеток в текущей капсуле при разных плотностях заселения.")

    sweep_res = run_density_sweep(
        R_outer_microns, p_boundary, D_coeff, 
        geometry=geometry_key,
        L_fibrosis_microns=L_fibrosis_microns,
        D_fibrosis=D_fibrosis
    )

    fig2 = gr_obj.Figure()
    fig2.add_trace(
        gr_obj.Scatter(
            x=sweep_res["densities"], 
            y=sweep_res["viabilities"], 
            name=f"Выбранная форма ({sav_name})", 
            line=dict(color="#10B981", width=4)
        )
    )

    fig2.add_trace(
        gr_obj.Scatter(
            x=sweep_res["densities"], 
            y=sweep_res["insulin_capacities"], 
            name="Функциональность секреции инсулина (%)", 
            line=dict(color="#3B82F6", width=2, dash="dash")
        )
    )

    if geometry_key != "planar":
        sweep_planar = run_density_sweep(
            R_outer_microns, p_boundary, D_coeff, 
            geometry="planar",
            L_fibrosis_microns=L_fibrosis_microns,
            D_fibrosis=D_fibrosis
        )
        fig2.add_trace(
            gr_obj.Scatter(
                x=sweep_planar["densities"], 
                y=sweep_planar["viabilities"], 
                name="Контроль: Плоский лист (1/L)", 
                line=dict(color="#EF4444", width=2, dash="dot"),
                opacity=0.6
            )
        )

    fig2.add_vline(
        x=rho_million, 
        line_width=2, 
        line_dash="dot", 
        line_color="#F59E0B",
        annotation_text=f"Выбранная плотность: {rho_million} млн/мл",
        annotation_position="top right"
    )

    fig2.update_layout(
        xaxis_title="Плотность заселения клеток (млн/мл)",
        yaxis_title="Процент (%)",
        template="plotly_dark",
        height=400,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=20, r=20, t=20, b=20),
        yaxis=dict(range=[0, 105])
    )
    st.plotly_chart(fig2, use_container_width=True)

    # Научный блок
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

    # Экспорт CAD
    st.write("---")
    st.subheader("💾 Экспорт CAD-модели для 3D-биопечати")
    st.markdown("Вы можете экспортировать геометрическую модель текущей капсулы в формат STL.")

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
            mime="application/sla"
        )
        
    with col_cad_right:
        st.info(f"**Экспортируемые параметры:** {file_name} | Вершин: {vertices.shape[0]} | Треугольников: {faces.shape[0]}. Единицы измерения: **микрометры**.")

# ==============================================================================
# РЕЖИМ 2: ГЕНЕРАТИВНЫЙ 3D-ДИЗАЙН (TPMS)
# ==============================================================================
elif app_mode == "🔮 Генеративный 3D-дизайн (TPMS)":
    st.sidebar.write("---")
    st.sidebar.header("📐 Геометрия TPMS каркаса")
    tpms_type = st.sidebar.radio(
        "Тип минимальной поверхности:",
        ["gyroid", "schwarz_p"],
        format_func=lambda x: "Гироид (Gyroid)" if x == "gyroid" else "Поверхность Шварца P (Schwarz P)"
    )

    box_size = st.sidebar.slider(
        "Размер контейнера куба (мкм)",
        200, 1500, 800,
        step=50,
        help="Размер ребра 3D-каркаса в микрометрах."
    )

    unit_cell = st.sidebar.slider(
        "Элементарная ячейка / Период (мкм)",
        100, 600, 250,
        step=10,
        help="Размер одного повторяющегося паттерна структуры."
    )

    thickness = st.sidebar.slider(
        "Толщина стенок (смещение t)",
        -1.0, 1.0, 0.0,
        step=0.1,
        help="Меньше нуля - толстые гидрогелевые стенки, больше нуля - тонкие пористые стенки."
    )

    resolution = st.sidebar.slider(
        "Разрешение voxel-сетки",
        30, 80, 50,
        step=5,
        help="30-50 подходит для мгновенного рендеринга на сайте. 60-80 - для качественного экспорта на 3D-биопринтер."
    )

    # --- ГЕНЕРАЦИЯ ---
    from generator_3d import TPMSGenerator
    with st.spinner("Генерация водонепроницаемой 3D-структуры TPMS..."):
        generator = TPMSGenerator(
            tpms_type=tpms_type,
            box_size=box_size,
            unit_cell=unit_cell,
            thickness=thickness,
            resolution=resolution
        )
        mesh = generator.build_mesh()

    if mesh is not None:
        # Вывод биофизических параметров
        st.subheader("📊 Аналитические характеристики 3D-каркаса")
        
        m1, m2, m3, m4 = st.columns(4)
        m1.metric(
            label="Пористость капсулы",
            value=f"{generator.porosity:.1f}%",
            delta="Свободный объем для клеток"
        )
        
        # Эквивалентный плоский лист 300 мкм имеет SA/V = 66.7 см⁻¹.
        # Сравним наш гироид с листом
        improvement = generator.sav_ratio / 66.7
        m2.metric(
            label="Удельная площадь SA/V",
            value=f"{generator.sav_ratio:.1f} см⁻¹",
            delta=f"В {improvement:.1f} раз лучше плоского листа"
        )
        
        m3.metric(
            label="Объем гидрогеля (Solid Volume)",
            value=f"{generator.solid_volume/1e6:.2f}M мкм³",
            delta=f"Доля пластика/геля: {generator.volume_fraction*100:.1f}%"
        )
        
        m4.metric(
            label="Общая площадь поверхности",
            value=f"{generator.surface_area/1e6:.2f}M мкм²",
            delta="Watertight manifold"
        )

        # 3D Рендеринг Plotly
        st.write("---")
        st.subheader("🔮 Интерактивный 3D-просмотр каркаса")
        st.markdown("Вы можете вращать, приближать и исследовать внутреннюю пористость сгенерированного каркаса.")

        verts = mesh.vertices
        faces = mesh.faces

        # Для ускорения рендеринга на веб-страницах, если граней слишком много, trimesh оставляет исходное
        fig3d = gr_obj.Figure(data=[
            gr_obj.Mesh3d(
                x=verts[:, 0],
                y=verts[:, 1],
                z=verts[:, 2],
                i=faces[:, 0],
                j=faces[:, 1],
                k=faces[:, 2],
                color='#3B82F6',
                opacity=0.85,
                flatshading=True
            )
        ])
        
        fig3d.update_layout(
            template="plotly_dark",
            scene=dict(
                xaxis=dict(title="X (мкм)", backgroundcolor="rgba(0,0,0,0)", gridcolor="rgba(255,255,255,0.1)"),
                yaxis=dict(title="Y (мкм)", backgroundcolor="rgba(0,0,0,0)", gridcolor="rgba(255,255,255,0.1)"),
                zaxis=dict(title="Z (мкм)", backgroundcolor="rgba(0,0,0,0)", gridcolor="rgba(255,255,255,0.1)"),
                aspectmode="data"
            ),
            margin=dict(l=10, r=10, t=10, b=10),
            height=600,
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )
        st.plotly_chart(fig3d, use_container_width=True)

        # Скачивание файла
        st.write("---")
        st.subheader("💾 Скачать CAD-модель (STL)")
        st.markdown("Данный файл полностью замкнут и готов к импорту в слайсер (CURA, PrusaSlicer) или САПР (SolidWorks, Autodesk Fusion).")
        
        col_down_l, col_down_r = st.columns([1, 2])
        with col_down_l:
            # Получаем ASCII строку
            stl_string = generator.get_stl_string()
            st.download_button(
                label=f"Скачать файл {tpms_type.upper()} STL",
                data=stl_string,
                file_name=f"tpms_{tpms_type}_box{box_size}um_cell{unit_cell}um.stl",
                mime="application/sla"
            )
        with col_down_r:
            st.info(f"Вершин: {verts.shape[0]} | Треугольников: {faces.shape[0]}. Единицы измерения в STL: **микрометры**.")
    else:
        st.error("Ошибка при генерации 3D-модели. Пожалуйста, измените граничные условия или уменьшите сдвиг толщины.")

# ==============================================================================
# РЕЖИМ 3: ML-ПОДБОР АНТИФИБРОЗНЫХ ПОКРЫТИЙ (GNN)
# ==============================================================================
elif app_mode == "🧪 ML-подбор антифиброзных покрытий (GNN)":
    st.header("🧪 Машинное обучение для борьбы с фиброзом (FBR)")
    st.markdown("""
    Когда имплантат попадает в организм, иммунная система реагирует реакцией на чужеродное тело (Foreign Body Response - FBR), 
    окутывая капсулу плотным слоем фиброзной ткани. Этот рубец блокирует диффузию кислорода и питательных веществ.
    
    Математическое моделирование показывает, что **даже идеальная сферическая микрокапсула теряет выживаемость клеток до 60.7% при толщине фиброза 50 мкм, а при 100 мкм наступает полный некроз**.
    
    Цель данного ML-модуля — **направленный скрининг молекулярных покрытий**, которые подавляют адгезию макрофагов и белков, 
    что позволяет снизить толщину фиброза до безопасного уровня ($L_{fib} \\le 20$ мкм), гарантирующего высокую жизнеспособность клеток.
    """)
    
    # Загружаем обученную модель
    @st.cache_resource
    def load_or_train_gnn():
        from gnn_pipeline import BiocompatibilityGNN, build_training_dataset, train_gnn_model
        model = BiocompatibilityGNN()
        weights_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "biocompatibility_gnn.pt")
        if os.path.exists(weights_path):
            model.load_state_dict(torch.load(weights_path))
        else:
            dataset = build_training_dataset()
            model, _ = train_gnn_model(dataset, epochs=150, lr=0.01, batch_size=8)
            torch.save(model.state_dict(), weights_path)
        model.eval()
        return model
        
    try:
        import os
        import torch
        from rdkit import Chem
        from rdkit.Chem import Draw
        from gnn_pipeline import smiles_to_graph
        
        model = load_or_train_gnn()
        st.success("🤖 Графовая нейросеть (GNN) успешно инициализирована и готова к работе.")
        
        # Предопределенный список кандидатов
        candidates = [
            {"name": "Zwitterion (Sulfobetaine SBAA)", "smiles": "C=CC(=O)NCC[N+](C)(C)CCCS(=O)(=O)[O-]", "class": "Zwitterionic"},
            {"name": "Zwitterion (Carboxybetaine CBAA)", "smiles": "C=CC(=O)NCC[N+](C)(C)CC(=O)[O-]", "class": "Zwitterionic"},
            {"name": "PEG8-acrylate", "smiles": "C=CC(=O)OCCOCCOCCOCCOCCOCCOCCOCCO C", "class": "PEGylated"},
            {"name": "PEG4-methacrylate", "smiles": "CC(=C)C(=O)OCCOCCOCCOCCO C", "class": "PEGylated"},
            {"name": "N,N-dimethylacrylamide (DMAA)", "smiles": "C=CC(=O)N(C)C", "class": "Neutral Hydrophilic"},
            {"name": "Hydroxyethyl methacrylate (HEMA)", "smiles": "CC(=C)C(=O)OCCO", "class": "Neutral Hydrophilic"},
            {"name": "Methyl methacrylate (MMA)", "smiles": "CC(=C)C(=O)OC", "class": "Hydrophobic"},
            {"name": "Styrene", "smiles": "C=CC1=CC=CC=C1", "class": "Hydrophobic"},
            {"name": "METAC (quaternary cationic)", "smiles": "CC(=C)C(=O)OCC[N+](C)(C)C", "class": "Cationic"},
            {"name": "Pentafluoropropyl methacrylate (PFPMA)", "smiles": "CC(=C)C(=O)OCC(F)(F)C(F)(F)F", "class": "Fluorinated Hydrophobic"}
        ]
        
        col_input, col_vis = st.columns([2, 1])
        
        with col_input:
            st.subheader("🧪 Выберите или спроектируйте покрытие")
            
            selection_mode = st.radio(
                "Способ ввода молекулы:",
                ["Выбрать из готовой базы кандидатов", "Ввести свой SMILES (Конструктор молекул)"]
            )
            
            selected_smiles = ""
            selected_name = ""
            selected_class = "Пользовательская молекула"
            
            if selection_mode == "Выбрать из готовой базы кандидатов":
                selected_idx = st.selectbox(
                    "База покрытий:",
                    range(len(candidates)),
                    format_func=lambda i: f"{candidates[i]['name']} ({candidates[i]['class']})"
                )
                selected_cand = candidates[selected_idx]
                selected_smiles = selected_cand["smiles"]
                selected_name = selected_cand["name"]
                selected_class = selected_cand["class"]
            else:
                custom_smiles = st.text_input(
                    "SMILES строка молекулы (мономера):",
                    value="C=CC(=O)NCC[N+](C)(C)CCCS(=O)(=O)[O-]",
                    help="Введите валидный SMILES. Например, C=CC(=O)N для акриламида."
                )
                # Проверяем валидность
                mol = Chem.MolFromSmiles(custom_smiles)
                if mol is None:
                    st.error("❌ Невалидный формат SMILES. Пожалуйста, проверьте структуру.")
                else:
                    selected_smiles = custom_smiles
                    selected_name = "Спроектированное покрытие"
                    selected_class = "Пользовательская молекула"
                    
            if selected_smiles and model is not None:
                # Оцениваем через GNN
                data = smiles_to_graph(selected_smiles)
                if data is not None:
                    with torch.no_grad():
                        batch = torch.zeros(data.x.size(0), dtype=torch.long)
                        pred = model(data.x, data.edge_index, batch)
                        biocompatibility = float(pred[0, 0].item())
                        
                    # Отображение на толщину фиброза
                    max_fib = 150.0
                    L_fib = max(0.0, max_fib * (1.0 - biocompatibility))
                    
                    st.write("---")
                    st.subheader("🎯 ML-Анализ покрытия:")
                    
                    # Метрики
                    m_biocompatibility, m_fibrosis = st.columns(2)
                    
                    status_color = "normal"
                    if biocompatibility > 0.90:
                        status_text = "Идеально (Цвиттер-ион)"
                    elif biocompatibility > 0.70:
                        status_text = "Отлично (Гидрофильный ПЭГ)"
                    elif biocompatibility > 0.40:
                        status_text = "Умеренно (HEMA/Гидрофильный)"
                    else:
                        status_text = "Критически (Вызывает FBR)"
                        status_color = "inverse"
                        
                    m_biocompatibility.metric(
                        label="Индекс биосовместимости (GNN)",
                        value=f"{biocompatibility:.4f}",
                        delta=status_text,
                        delta_color=status_color
                    )
                    
                    # Оценка выживаемости на основе L_fib
                    if L_fib <= 20.0:
                        fib_status = "Безопасно (L <= 20 мкм)"
                    elif L_fib <= 50.0:
                        fib_status = "Умеренный фиброз"
                    else:
                        fib_status = "Тяжелый фиброз (Критический)"
                        
                    m_fibrosis.metric(
                        label="Ожидаемая толщина фиброза (L_fib)",
                        value=f"{L_fib:.1f} мкм",
                        delta=fib_status,
                        delta_color="normal" if L_fib <= 20.0 else "inverse"
                    )
                    
                    # Кнопка для импорта в 1D симуляцию
                    st.write("")
                    if st.button("🚀 Применить это покрытие (экспортировать L_fib в симулятор)"):
                        st.session_state["l_fibrosis_microns"] = float(np.round(L_fib, 1))
                        st.success(f"Толщина фиброза {L_fib:.1f} мкм успешно перенесена в настройки симулятора! Перейдите на вкладку '1D Симуляция диффузии O₂', чтобы рассчитать профиль.")
                else:
                    st.error("Не удалось построить граф из SMILES.")
                    
        with col_vis:
            st.subheader("🖼️ Структурная формула")
            if selected_smiles:
                mol = Chem.MolFromSmiles(selected_smiles)
                if mol is not None:
                    img = Draw.MolToImage(mol, size=(300, 300))
                    st.image(img, use_column_width=True, caption=selected_name)
                else:
                    st.warning("Молекулярная визуализация недоступна.")
                    
        # Общая таблица скрининга
        st.write("---")
        st.subheader("📊 Рейтинг кандидатов (Сводная таблица скрининга)")
        st.markdown("Ниже приведена таблица сравнения всех 10 эталонных кандидатов нашей базы, предсказанных графовой нейросетью.")
        
        if model is not None:
            table_data = []
            with torch.no_grad():
                for c in candidates:
                    d = smiles_to_graph(c["smiles"])
                    if d is not None:
                        batch = torch.zeros(d.x.size(0), dtype=torch.long)
                        pred = model(d.x, d.edge_index, batch)
                        biocompatibility = float(pred[0, 0].item())
                        L_fib = max(0.0, 150.0 * (1.0 - biocompatibility))
                        table_data.append({
                            "Название покрытия": c["name"],
                            "Класс химического соединения": c["class"],
                            "SMILES": c["smiles"],
                            "Индекс биосовместимости": round(biocompatibility, 4),
                            "Ожидаемый фиброз L_fib (мкм)": round(L_fib, 1)
                        })
            df = pd.DataFrame(table_data)
            df = df.sort_values(by="Индекс биосовместимости", ascending=False).reset_index(drop=True)
            st.dataframe(df, use_container_width=True)
            
    except Exception as ex:
        st.error(f"Произошла ошибка при запуске GNN скрининга: {ex}")

# ==============================================================================
# РЕЖИМ 4: СИМУЛЯЦИЯ НЕОВАСКУЛЯРИЗАЦИИ (VEGF)
# ==============================================================================
elif app_mode == "🩸 Неоваскуляризация (VEGF / Ангиогенез)":
    st.header("🩸 Моделирование неоваскуляризации (Angiogenesis & VEGF Release)")
    st.markdown(r"""
    Даже при отсутствии плотного фиброзного рубца ($L_{fib} = 0$), подкожная жировая ткань характеризуется низким уровнем оксигенации ($pO_2 \\approx 30$ mmHg). 
    Для обеспечения долгосрочной выживаемости терапевтических доз клеток ($\ge 80$ млн/мл) необходимо стимулировать прорастание капилляров непосредственно к границе капсулы.
    
    В данном режиме симулируется выделение фактора роста эндотелия сосудов (VEGF), загруженного в гидрогелевую мембрану, его диффузия в ткань 
    и динамический рост давления кислорода на внешней границе по мере приживления трансплантата (angiogenesis feedback loop).
    """)
    
    # Сайдбар настройки для VEGF
    st.sidebar.write("---")
    st.sidebar.header("🩸 Параметры VEGF & Ангиогенеза")
    
    species_selected = st.sidebar.selectbox(
        "Биологический вид (Species):",
        ["Mouse", "Human"],
        help="Для человека таймлайн ангиогенеза дольше (60 суток), а скорость прорастания сосудов снижена in vivo."
    )
    
    v_loaded = st.sidebar.slider(
        "Начальная загрузка VEGF в гель",
        0.0, 5.0, 1.0,
        step=0.5,
        help="Относительная начальная концентрация фактора роста внутри капсулы."
    )
    
    k_clear = st.sidebar.slider(
        "Скорость тканевого клиренса VEGF (1/день)",
        5.0, 30.0, 15.0,
        step=1.0,
        help="Описывает скорость вымывания и связывания VEGF рецепторами."
    )
    
    beta_angio = st.sidebar.slider(
        "Скорость роста сосудов beta (1/день)",
        0.05, 0.5, 0.15,
        step=0.01,
        help="Определяет биологическую скорость прорастания капилляров (лаг приживления)."
    )
    
    p_max_angio = st.sidebar.slider(
        "Предельное pO2 сосудов (mmHg)",
        40.0, 95.0, 60.0,
        step=5.0,
        help="Максимальное давление кислорода при полной васкуляризации границы капсулы."
    )

    with st.sidebar.expander("🛡️ Параметры Фазы 9 (Гемодинамика и Иммунитет)"):
        av_loop_flow = st.checkbox(
            "Подключить прямой кровоток (AV-loop)",
            value=False,
            help="Артериовенозный шунт."
        )
        if av_loop_flow:
            tau_blood = st.slider(
                "Сдвиговое напряжение AV-петли (Па)",
                0.1, 12.0, 5.0,
                step=0.1
            )
            anticoagulation = st.checkbox(
                "Антикоагулянты",
                value=False
            )
        else:
            tau_blood = 5.0
            anticoagulation = False
            
        crispr_hypoimmune = st.checkbox(
            "CRISPR Hypoimmune Клетки",
            value=False
        )
        if crispr_hypoimmune:
            cd47_overexpression = st.checkbox(
                "Гиперэкспрессия CD47",
                value=False
            )
            complement_protection = st.checkbox(
                "CD55/CD59 ингибиторы",
                value=False
            )
        else:
            cd47_overexpression = False
            complement_protection = False
            
        turnover_rate = st.slider(
            "Скорость апоптоза (Turnover, %/сутки)",
            0.0, 5.0, 1.5,
            step=0.1
        ) / 100.0
        
        phi_pfc = st.slider(
            "Доля PFC в геле",
            0.0, 0.3, 0.0,
            step=0.05
        )
        if phi_pfc > 0.0:
            pO2_pfc_saturation = st.slider(
                "Насыщение PFC (pO₂, mmHg)",
                150.0, 760.0, 200.0,
                step=10.0
            )
        else:
            pO2_pfc_saturation = 200.0
    
    st.sidebar.header("📐 Геометрические параметры")
    geo_vegf = st.sidebar.radio(
        "Форма капсулы (VEGF):",
        ["planar", "cylindrical", "spherical"],
        index=2,
        format_func=lambda x: {
            "planar": "Плоская пластина (Slab)",
            "cylindrical": "Цилиндрическая нить (Fiber)",
            "spherical": "Сферическая микрокапсула"
        }[x]
    )
    
    r_vegf = st.sidebar.slider(
        "Радиус / Полутолщина R (мкм):",
        50, 400, 150,
        step=10
    )
    
    rho_vegf = st.sidebar.slider(
        "Плотность клеток (млн/мл):",
        5, 200, 80,
        step=5
    )
    
    l_fib = st.sidebar.slider(
        "Толщина фиброза L_fib (мкм):",
        0.0, 150.0, float(st.session_state["l_fibrosis_microns"]),
        step=1.0,
        help="Фиброзный рубец замедляет диффузию VEGF в ткань."
    )
    st.session_state["l_fibrosis_microns"] = l_fib
    
    # Расчет
    from simulator import run_neovascularization_sweep_oxygen
    with st.spinner("Численное интегрирование уравнений сопряженного ангиогенеза..."):
        res_vegf = run_neovascularization_sweep_oxygen(
            R_outer_microns=r_vegf,
            rho_million_per_ml=rho_vegf,
            D_oxygen_coefficient=1.5e-5,
            geometry=geo_vegf,
            L_fibrosis_microns=l_fib,
            D_fibrosis=1.0e-5,
            V_loaded_relative=v_loaded,
            k_clear_tissue=k_clear,
            beta_angiogenesis=beta_angio,
            K_vegf=0.1,
            p_base=30.0,
            p_max=p_max_angio,
            days=60 if species_selected == "Human" else 21,
            species=species_selected,
            av_loop_flow=av_loop_flow,
            crispr_hypoimmune=crispr_hypoimmune,
            cd47_overexpression=cd47_overexpression,
            tau_blood=tau_blood,
            anticoagulation=anticoagulation,
            pO2_pfc_saturation=pO2_pfc_saturation,
            turnover_rate=turnover_rate,
            complement_protection=complement_protection
        )
        
    # Динамические графики
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("📈 Концентрация VEGF и Давление O₂ на границе")
        st.markdown("График показывает, как VEGF диффундирует на границу раздела сред и как это стимулирует локальную оксигенацию со временем.")
        
        fig_timeline = gr_obj.Figure()
        fig_timeline.add_trace(
            gr_obj.Scatter(
                x=res_vegf["t"],
                y=res_vegf["C_interface"],
                name="VEGF на границе (отн. ед.)",
                line=dict(color="#10B981", width=3)
            )
        )
        fig_timeline.add_trace(
            gr_obj.Scatter(
                x=res_vegf["t"],
                y=res_vegf["p_boundary"],
                name="Локальное pO₂ на границе (mmHg)",
                line=dict(color="#60A5FA", width=3),
                yaxis="y2"
            )
        )
        
        fig_timeline.update_layout(
            template="plotly_dark",
            xaxis_title="Время (сутки)",
            yaxis_title="Концентрация VEGF",
            yaxis2=dict(
                title="Парциальное давление O₂ (mmHg)",
                overlaying="y",
                side="right",
                range=[25, p_max_angio + 5]
            ),
            legend=dict(x=0.05, y=0.95),
            margin=dict(l=20, r=20, t=20, b=20),
            height=380
        )
        st.plotly_chart(fig_timeline, use_container_width=True)
        
    with col2:
        st.subheader("💖 Выживаемость и Инсулиновая емкость")
        st.markdown("График показывает процесс терапевтического спасения клеток: рост жизнеспособности и секреции инсулина по мере прорастания капилляров.")
        
        fig_metrics = gr_obj.Figure()
        fig_metrics.add_trace(
            gr_obj.Scatter(
                x=res_vegf["t"],
                y=res_vegf["viability_over_time"],
                name="Жизнеспособность клеток (%)",
                line=dict(color="#EF4444", width=3)
            )
        )
        fig_metrics.add_trace(
            gr_obj.Scatter(
                x=res_vegf["t"],
                y=res_vegf["insulin_over_time"],
                name="Секреция инсулина (%)",
                line=dict(color="#F59E0B", width=2, dash="dash")
            )
        )
        
        fig_metrics.update_layout(
            template="plotly_dark",
            xaxis_title="Время (сутки)",
            yaxis_title="Доля (%)",
            yaxis=dict(range=[0, 105]),
            legend=dict(x=0.05, y=0.25),
            margin=dict(l=20, r=20, t=20, b=20),
            height=380
        )
        st.plotly_chart(fig_metrics, use_container_width=True)

        # Phase 9: CD47 & complement warnings in Mode 4
        if crispr_hypoimmune:
            if not cd47_overexpression:
                st.error("⚠️ **АТАКА NK-КИЛЛЕРОВ ('Missing Self')!** Полное отсутствие MHC-I на поверхности клеток активирует NK-клетки. Трансплантат будет уничтожен без гиперэкспрессии CD47.")
            if not complement_protection:
                st.error("⚠️ **ЛИЗИС СИСТЕМОЙ КОМПЛЕМЕНТА!** Без мембранных ингибиторов комплемента (CD55/CD59) белки крови пробивают поры в клетках за минуты.")
            if cd47_overexpression and complement_protection:
                st.success("✅ **Генетический щит активен!** Клетки защищены гиперэкспрессией CD47 от NK-клеток и CD55/CD59 от системы комплемента.")

        # Phase 9: Hyperoxia ROS stress warning in Mode 4
        if phi_pfc > 0.0 and pO2_pfc_saturation > 200.0:
            st.error(f"⚠️ **ГИПЕРОКСИЧЕСКИЙ ROS-ШОК!** Начальное насыщение PFC кислородом ({pO2_pfc_saturation:.1f} mmHg) превышает безопасные 200 mmHg. Бета-клетки страдают от мгновенного коллапса митохондрий под действием ROS.")

        # Phase 9: AV-loop thrombosis & hyperplasia warnings in Mode 4
        if av_loop_flow:
            if tau_blood < 1.5 or tau_blood > 8.0:
                st.error(f"⚠️ **ГЕМОДИНАМИЧЕСКИЙ КРИЗ (Тромбоз AV-петли)!** Сдвиговое напряжение ({tau_blood:.1f} Па) вне физиологической нормы [1.5, 8.0] Па. Поток крови замедляется тромбом.")
            st.warning(f"🩸 **Интимальная гиперплазия:** В месте анастомоза разрастается гладкая мускулатура (k = 0.005/день). Пропускная способность шунта снизится в долгосрочной перспективе.")
        
    # Пространственные профили VEGF по дням
    st.write("---")
    col3, col4 = st.columns([1, 2])
    
    with col3:
        st.subheader("⏱️ Интерактивный тайм-слайдер")
        st.markdown("Передвигайте ползунок, чтобы увидеть распределение концентраций внутри и снаружи капсулы в конкретные сутки.")
        max_days = 60 if species_selected == "Human" else 21
        default_selected_day = 14 if species_selected == "Human" else 7
        selected_day = st.slider("Выберите день для анализа:", 0, max_days, default_selected_day, step=1)
        
        # Вычисление метрик для выбранного дня
        day_idx = np.argmin(np.abs(res_vegf["t"] - selected_day))
        current_pO2_bound = res_vegf["p_boundary"][day_idx]
        current_viability = res_vegf["viability_over_time"][day_idx]
        current_insulin = res_vegf["insulin_over_time"][day_idx]
        
        st.metric("Давление O₂ на границе", f"{current_pO2_bound:.1f} mmHg")
        st.metric("Выживаемость клеток", f"{current_viability:.1f}%")
        st.metric("Секреция инсулина", f"{current_insulin:.1f}%")
        
    with col4:
        st.subheader("🔬 Пространственный профиль VEGF")
        st.markdown("Отображает профиль VEGF внутри полимера (0 - R) и его диффузию/затухание в окружающих тканях.")
        
        saved_profiles = res_vegf["saved_profiles"]
        r_coords = res_vegf["r_grid"]
        
        available_days = sorted(list(saved_profiles.keys()))
        nearest_day = available_days[np.argmin(np.abs(np.array(available_days) - selected_day))]
        profile_to_draw = saved_profiles[nearest_day]
        
        fig_profile = gr_obj.Figure()
        fig_profile.add_trace(
            gr_obj.Scatter(
                x=r_coords,
                y=profile_to_draw,
                name=f"VEGF Профиль (День {nearest_day})",
                line=dict(color="#10B981", width=4),
                fill="tozeroy",
                fillcolor="rgba(16, 185, 129, 0.1)"
            )
        )
        
        # Разделительная линия Capsule/Tissue
        fig_profile.add_vline(
            x=r_vegf,
            line_width=2,
            line_dash="dash",
            line_color="rgba(255,255,255,0.5)",
            annotation_text="Граница капсулы (R)",
            annotation_position="top left"
        )
        
        if l_fib > 0:
            fig_profile.add_vline(
                x=r_vegf + l_fib,
                line_width=2,
                line_dash="dot",
                line_color="rgba(239, 68, 68, 0.5)",
                annotation_text="Граница фиброза",
                annotation_position="top right"
            )
            fig_profile.add_vrect(
                x0=r_vegf, x1=r_vegf + l_fib,
                fillcolor="#EF4444", opacity=0.08, line_width=0
            )
            
        fig_profile.update_layout(
            template="plotly_dark",
            xaxis_title="Радиальная координата (мкм)",
            yaxis_title="Концентрация VEGF (отн. ед.)",
            yaxis=dict(range=[0, max(1.1 * v_loaded, 1.1)]),
            margin=dict(l=20, r=20, t=20, b=20),
            height=380
        )
        st.plotly_chart(fig_profile, use_container_width=True)

# ==============================================================================
# РЕЖИМ 5: МИНИ-ОРГАНОИДЫ (ФАЗА 10: BIOMIMESIS)
# ==============================================================================
elif app_mode == "🧫 Мини-органоиды (Фаза 10: Biomimesis)":
    st.header("🧫 Моделирование мини-органоидов (The Ultimate Biomimesis)")
    st.markdown(r"""
    В Фазе 10 мы полностью отказываемся от искусственных полимерных капсул в пользу **"голых" гипоиммунных мини-органоидов**, 
    трансплантируемых непосредственно во внутрипортальный кровоток печени. Это полностью решает проблему реакции на чужеродное тело (FBR), 
    устраняет диффузионный барьер полимерного геля и нормализует инсулиновый метаболизм за счет первого прохода через печень (First-Pass Hepatic Extraction).
    """)
    
    # Сайдбар настройки для Фазы 10
    st.sidebar.write("---")
    st.sidebar.header("🧫 Параметры мини-органоидов")
    
    N_0 = st.sidebar.slider(
        "Начальная популяция клеток N₀ (млн)",
        10.0, 500.0, 100.0,
        step=10.0,
        help="Исходное количество терапевтических бета-клеток в трансплантируемом объеме."
    )
    
    turnover_rate_10 = st.sidebar.slider(
        "Скорость апоптоза (Turnover, %/сутки)",
        0.1, 5.0, 1.0,
        step=0.1,
        help="Базовая скорость физиологической гибели бета-клеток."
    ) / 100.0
    
    phi_epc = st.sidebar.slider(
        "Доля эндотелиальных клеток EPCs (%)",
        0.0, 20.0, 10.0,
        step=1.0,
        help="Ко-инкапсулированные EPCs и MSCs стимулируют самосборку внутренней сосудистой сети органоида за первые дни."
    ) / 100.0
    
    st.sidebar.header("🧬 Конфигуратор CRISPR-редактирования")
    b2m_ko = st.sidebar.checkbox("Knock-out B2M (Удаление MHC Class I)", value=False, help="Защита от T-лимфоцитов. ВНИМАНИЕ: без CD47 активирует NK-клетки!")
    ciita_ko = st.sidebar.checkbox("Knock-out CIITA (Удаление MHC Class II)", value=False, help="Защита от хелперных T-клеток и активации макрофагов.")
    cd47_ki = st.sidebar.checkbox("Knock-in CD47 (Защитный сигнал)", value=False, help="Предотвращает атаку NK-клеток ('Missing Self') и фагоцитоз макрофагами.")
    cd55_cd59_ki = st.sidebar.checkbox("Knock-in CD55/CD59 (Ингибиторы комплемента)", value=False, help="Предотвращает лизис клеток белками системы комплемента плазмы.")
    pdl1_ki = st.sidebar.checkbox("Knock-in PD-L1 (Подавление T-клеток)", value=False, help="Локальное ингибирование аутореактивных T-лимфоцитов хозяина.")
    
    # Расчеты
    t_10_years = np.linspace(0, 3650, 500) # 10 лет
    pop_timeline = simulate_organoid_population(
        t_10_years,
        N_0=N_0,
        N_stem_fraction=0.05,
        r_proliferation=0.02,
        turnover_rate=turnover_rate_10,
        b2m_ko=b2m_ko,
        ciita_ko=ciita_ko,
        cd47_ki=cd47_ki,
        cd55_cd59_ki=cd55_cd59_ki,
        pdl1_ki=pdl1_ki
    )
    
    t_30_days = np.linspace(0, 30, 300) # 30 дней для васкуляризации
    pO2_timeline = simulate_organoid_oxygenation(t_30_days, phi_epc=phi_epc)
    
    # Расчет инсулина на первые 30 дней с учетом васкуляризации и выживания клеток
    N_cells_30 = simulate_organoid_population(
        t_30_days,
        N_0=N_0,
        N_stem_fraction=0.05,
        r_proliferation=0.02,
        turnover_rate=turnover_rate_10,
        b2m_ko=b2m_ko,
        ciita_ko=ciita_ko,
        cd47_ki=cd47_ki,
        cd55_cd59_ki=cd55_cd59_ki,
        pdl1_ki=pdl1_ki
    )
    ins_portal, ins_systemic = simulate_organoid_insulin(t_30_days, N_cells_30, pO2_timeline)
    
    # Иммунологический вердикт
    st.subheader("🛡️ Иммунологический статус трансплантата")
    k_leak = calculate_immune_leak(b2m_ko, ciita_ko, cd47_ki, cd55_cd59_ki, pdl1_ki)
    
    if k_leak == 0.0:
        st.success("🥇 **Идеальный статус (Off-the-shelf Universal Graft)!** Все векторы иммунного ответа (T-клетки, NK, комплемент, макрофаги) полностью заблокированы. Прогнозируется 10-летняя выживаемость органоида.")
    elif b2m_ko and not cd47_ki:
        st.error("🚨 **Патофизиологический кризис (Missing Self)!** Нокаут B2M спасает от T-клеток, но отсутствие MHC-I без оверэкспрессии CD47 вызывает мгновенный лизис натуральными киллерами (NK-клетки).")
    elif not b2m_ko and not pdl1_ki:
        st.error("⚠️ **Тканевая несовместимость (Аллореактивное отторжение)!** Клетки экспрессируют MHC I, но не защищены PD-L1. Иммунная система хозяина быстро уничтожит органоид с помощью цитотоксических T-лимфоцитов.")
    elif not cd55_cd59_ki:
        st.warning("⚠️ **Уязвимость к гуморальному иммунитету (Лизис комплементом)!** Отсутствие белков CD55/CD59 приведет к пробиванию пор в клеточной мембране белками сыворотки крови (MAC комплекс).")
    elif not ciita_ko:
        st.warning("⚠️ **Вялотекущее воспаление (Презентация антигена)!** Экспрессия MHC II (CIITA) приведет к постоянному привлечению и активации макрофагов.")
    else:
        st.info("ℹ️ **Частичная защита.** Органоид частично модифицирован, однако присутствует остаточная иммунная утечка.")
        
    # Колонки с графиками
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("📈 Выживаемость популяции клеток (Горизонт 10 лет)")
        st.markdown("Показывает численность популяции клеток (в миллионах) на долгосрочном интервале времени под влиянием естественного оборота и иммунного лизиса.")
        
        fig_pop = gr_obj.Figure()
        fig_pop.add_trace(
            gr_obj.Scatter(
                x=t_10_years / 365.0, # в годах
                y=pop_timeline,
                name="Популяция клеток (млн)",
                line=dict(color="#60A5FA", width=4),
                fill="tozeroy",
                fillcolor="rgba(96, 165, 250, 0.1)"
            )
        )
        fig_pop.update_layout(
            template="plotly_dark",
            xaxis_title="Время (годы)",
            yaxis_title="Количество клеток N (млн)",
            yaxis=dict(range=[0, N_0 * 1.5]),
            margin=dict(l=20, r=20, t=20, b=20),
            height=380
        )
        st.plotly_chart(fig_pop, use_container_width=True)
        
    with col2:
        st.subheader("🩸 Кинетика инсулина (Портальный vs Системный)")
        st.markdown("Сравнение концентрации инсулина в воротной вене печени и системном кровотоке. Демонстрирует 60% печеночную экстракцию для предотвращения инсулинорезистентности.")
        
        fig_ins = gr_obj.Figure()
        fig_ins.add_trace(
            gr_obj.Scatter(
                x=t_30_days,
                y=ins_portal,
                name="Портальный инсулин (в воротной вене)",
                line=dict(color="#F59E0B", width=3)
            )
        )
        fig_ins.add_trace(
            gr_obj.Scatter(
                x=t_30_days,
                y=ins_systemic,
                name="Системный инсулин (периферия)",
                line=dict(color="#EC4899", width=3, dash="dash")
            )
        )
        fig_ins.update_layout(
            template="plotly_dark",
            xaxis_title="Время (сутки)",
            yaxis_title="Концентрация инсулина (отн. ед.)",
            legend=dict(x=0.05, y=0.95),
            margin=dict(l=20, r=20, t=20, b=20),
            height=380
        )
        st.plotly_chart(fig_ins, use_container_width=True)
        
    # Внутренняя оксигенация (васкуляризация)
    st.write("---")
    col_ox_desc, col_ox_plot = st.columns([1, 2])
    with col_ox_desc:
        st.subheader("🫀 Самосборка микрососудистой сети")
        st.markdown(f"""
        Ко-инкапсулированные предшественники эндотелия (EPCs) формируют капилляры непосредственно **внутри** микро-органоида.
        Это обеспечивает быструю оксигенацию клеток с первых суток и исключает "окно смерти" гипоксии.
        
        *   **Доля EPCs:** {phi_epc * 100:.1f}%
        *   **Портальное pO₂:** {pO2_timeline[-1]:.1f} mmHg (насыщение к 30 дню)
        *   **Время полуваскуляризации:** {np.log(2) / (0.1 * (1.0 + 5.0 * phi_epc)):.1f} суток
        """)
    with col_ox_plot:
        fig_ox = gr_obj.Figure()
        fig_ox.add_trace(
            gr_obj.Scatter(
                x=t_30_days,
                y=pO2_timeline,
                name="Внутреннее pO₂ органоида (mmHg)",
                line=dict(color="#10B981", width=4),
                fill="tozeroy",
                fillcolor="rgba(16, 185, 129, 0.1)"
            )
        )
        fig_ox.update_layout(
            template="plotly_dark",
            xaxis_title="Время (сутки)",
            yaxis_title="Парциальное давление O₂ (mmHg)",
            yaxis=dict(range=[0, 55]),
            margin=dict(l=20, r=20, t=20, b=20),
            height=300
        )
        st.plotly_chart(fig_ox, use_container_width=True)
