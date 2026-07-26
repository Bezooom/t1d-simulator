"""
Panel 5: Mini-Organoids Phase 10 (Biomimesis).
Handles organoid population, genomic editing, IBMIR protection,
clinical calculator, OGTT, CGM, and clinical export.
"""
import streamlit as st
import numpy as np
import plotly.graph_objects as go

from t1d_simulator.organoid_simulator import (
    simulate_organoid_population,
    simulate_organoid_oxygenation,
    simulate_transplantation_site_comparison,
    simulate_ibmir_protection,
    simulate_base_editing_fidelity,
    calculate_immune_leak,
    calculate_patient_transplant_dose,
    simulate_ogtt_glycemic_control,
    simulate_cgm_30day_metrics,
    evaluate_patient_clinical_risk_profile,
)
from t1d_simulator.organoid_cad_exporter import (
    generate_patient_clinical_passport,
    generate_omental_scaffold_stl,
)


def render_organoid_panel():
    """Render the full Phase 10 Organoids panel."""

    st.header("🧫 Моделирование мини-органоидов (Фаза 10: Biomimesis)")
    st.markdown("""
    В Фазе 10 мы полностью отказываемся от искусственных полимерных капсул в пользу
    **голых гипоиммунных мини-органоидов**, трансплантируемых непосредственно на сальник
    или в портальный кровоток. Это устраняет диффузионный барьер полимерного геля,
    обеспечивает быструю внутреннюю оксигенацию и гарантирует 100% биологическую
    и онкогенную безопасность.
    """)

    # --- Sidebar parameters ---
    st.sidebar.write("---")
    st.sidebar.header("🧫 Параметры мини-органоидов")

    N_0 = st.sidebar.slider(
        "Начальная популяция клеток N₀ (млн)",
        10.0, 500.0, 100.0, step=10.0,
        help="Исходное количество терапевтических бета-клеток.",
    )
    turnover_rate_10 = st.sidebar.slider(
        "Скорость апоптоза (Turnover, %/сутки)",
        0.1, 5.0, 1.0, step=0.1,
        help="Базовая скорость физиологической гибели бета-клеток.",
    ) / 100.0
    phi_epc = st.sidebar.slider(
        "Доля эндотелиальных клеток EPCs (%)",
        0.0, 20.0, 10.0, step=1.0,
        help="Ко-инкапсулированные EPCs и MSCs стимулируют самосборку сосудистой сети.",
    ) / 100.0

    # --- Genomic editing ---
    st.sidebar.header("🧬 Геномное редактирование (CRISPR)")
    edit_method = st.sidebar.radio(
        "Технология редактирования ДНК:",
        ["Base Editing (CBE/ABE)", "Classic SpCas9 Nuclease"],
        help="Base Editing вводит точечные замены без DSB и транслокаций.",
    )

    b2m_ko = st.sidebar.checkbox(
        "Knock-out B2M (Удаление MHC Class I)", value=True,
        help="Защита от T-лимфоцитов. ВНИМАНИЕ: без CD47 активирует NK-клетки!",
    )
    ciita_ko = st.sidebar.checkbox(
        "Knock-out CIITA (Удаление MHC Class II)", value=True,
        help="Защита от хелперных T-клеток и активации макрофагов.",
    )
    cd47_ki = st.sidebar.checkbox(
        "Knock-in CD47 (Защитный сигнал)", value=True,
        help="Предотвращает атаку NK-клеток ('Missing Self').",
    )
    cd55_cd59_ki = st.sidebar.checkbox(
        "Knock-in CD55/CD59 (Ингибиторы комплемента)", value=True,
        help="Предотвращает лизис клеток белками системы комплемента.",
    )
    pdl1_ki = st.sidebar.checkbox(
        "Knock-in PD-L1 (Подавление T-клеток)", value=True,
        help="Локальное ингибирование аутореактивных T-лимфоцитов хозяина.",
    )

    # --- Nanochemistry ---
    st.sidebar.header("🛡️ Нанохимия & Защита от IBMIR")
    peg_lmwh_density = st.sidebar.slider(
        "Плотность Lipid-PEG-LMWH щита",
        0.0, 2.0, 1.0, step=0.1,
        help="Модификация мембраны нейтрализует Тканевой Фактор (CD142) и IBMIR.",
    )

    # --- Transplantation site ---
    st.sidebar.header("📍 Анатомический сайт трансплантации")
    transplant_site = st.sidebar.radio(
        "Зона пересадки:",
        ["Большой сальник (Omental Pouch)", "Воротная вена печени (Portal Vein)"],
        format_func=lambda x: x,
    )
    site_key_code = "omental_pouch" if "сальник" in transplant_site else "portal_vein"

    # --- Emergency switch ---
    st.sidebar.header("🚨 Система iCasp9 (Emergency Suicide Switch)")
    ap1903_conc = st.sidebar.slider(
        "Концентрация AP1903 / Rimiducid (нМ)",
        0.0, 50.0, 0.0, step=1.0,
        help="При выявлении тератокарцином запускает апоптоз каспазы-9.",
    )

    # =====================================================================
    # CALCULATIONS
    # =====================================================================
    t_10_years = np.linspace(0, 3650, 500)
    pop_timeline = simulate_organoid_population(
        t_10_years, N_0=N_0, N_stem_fraction=0.05,
        r_proliferation=0.02, turnover_rate=turnover_rate_10,
        b2m_ko=b2m_ko, ciita_ko=ciita_ko, cd47_ki=cd47_ki,
        cd55_cd59_ki=cd55_cd59_ki, pdl1_ki=pdl1_ki,
    )

    t_30_days = np.linspace(0, 30, 300)
    pO2_timeline = simulate_organoid_oxygenation(t_30_days, phi_epc=phi_epc)

    N_cells_30 = simulate_organoid_population(
        t_30_days, N_0=N_0, N_stem_fraction=0.05,
        r_proliferation=0.02, turnover_rate=turnover_rate_10,
        b2m_ko=b2m_ko, ciita_ko=ciita_ko, cd47_ki=cd47_ki,
        cd55_cd59_ki=cd55_cd59_ki, pdl1_ki=pdl1_ki,
    )

    site_res = simulate_transplantation_site_comparison(t_30_days, N_cells_30, site=site_key_code)
    ibmir_res = simulate_ibmir_protection(peg_lmwh_density=peg_lmwh_density)
    genomic_res = simulate_base_editing_fidelity(
        b2m_ko=b2m_ko, ciita_ko=ciita_ko, cd47_ki=cd47_ki,
        cd55_cd59_ki=cd55_cd59_ki, pdl1_ki=pdl1_ki,
        method="base_editing" if "Base" in edit_method else "spcas9",
    )

    # =====================================================================
    # METRICS
    # =====================================================================
    m1, m2, m3, m4 = st.columns(4)
    m1.metric(
        "48ч Выживание (IBMIR)",
        f"{ibmir_res['retention_48h_percent']:.1f}%",
        delta="Тромбин заблокирован" if peg_lmwh_density > 0.5 else "Высокий IBMIR",
    )
    m2.metric(
        "Геномная безопасность",
        f"{100.0 - genomic_res['translocation_risk_percent']:.1f}%",
        delta=f"Риск транслокаций: {genomic_res['translocation_risk_percent']:.1f}%",
    )
    m3.metric(
        "Риск стеатоза печени",
        f"{site_res['steatosis_risk_index']:.1f}%",
        delta="Низкий (Сальник)" if site_key_code == "omental_pouch" else "Высокий (Печень)",
        delta_color="normal" if site_key_code == "omental_pouch" else "inverse",
    )
    m4.metric(
        "Извлекаемость органоида",
        f"{site_res['retrievability_score']:.0f}%",
        delta="100% хирургическая" if site_key_code == "omental_pouch" else "Неизвлекаем",
    )

    # =====================================================================
    # IMMUNOLOGY STATUS
    # =====================================================================
    st.subheader("🛡️ Анализ безопасности и иммунологии трансплантата")
    k_leak = calculate_immune_leak(b2m_ko, ciita_ko, cd47_ki, cd55_cd59_ki, pdl1_ki)

    if ap1903_conc > 0.0:
        st.error(
            f"🚨 **Аварийный протокол апоптоза (AP1903 = {ap1903_conc:.1f} нМ)!** "
            f"Активирована каспаза-9. Запущен процесс уничтожения 100% трансплантата."
        )
    elif k_leak == 0.0 and peg_lmwh_density >= 0.8:
        st.success(
            "🥇 **Идеальный статус (Off-the-shelf Universal Graft)!** "
            "Полная иммунная маскировка + защиты от IBMIR гепариновым гликокаликсом."
        )
    elif b2m_ko and not cd47_ki:
        st.error(
            "🚨 **Патофизиологический кризис (Missing Self)!** "
            "Нокаут B2M спасает от T-клеток, но отсутствие MHC-I без CD47 вызывает "
            "мгновенный лизис натуральными киллерами."
        )
    elif not b2m_ko and not pdl1_ki:
        st.error(
            "⚠️ **Тканевая несовместимость (Аллореактивное отторжение)!** "
            "Клетки экспрессируют MHC I, но не защищены PD-L1."
        )
    elif peg_lmwh_density < 0.5:
        st.warning(
            "⚠️ **Риск ранней гибели от IBMIR!** Без PEG-LMWH тромбоциты и "
            "Тканевой Фактор (CD142) вызовут локальный тромбоз."
        )
    else:
        st.info("ℹ️ **Частичная защита.** Органоид модифицирован, но присутствуют остаточные риски.")

    # =====================================================================
    # POPULATION & INSULIN PLOTS
    # =====================================================================
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("📈 Выживаемость популяции клеток (Горизонт 10 лет)")
        fig_pop = go.Figure()
        fig_pop.add_trace(go.Scatter(
            x=t_10_years / 365.0, y=pop_timeline,
            name="Популяция клеток (млн)",
            line=dict(color="#60A5FA", width=4),
            fill="tozeroy", fillcolor="rgba(96, 165, 250, 0.1)",
        ))
        fig_pop.update_layout(
            template="plotly_dark",
            xaxis_title="Время (годы)", yaxis_title="Количество клеток N (млн)",
            height=350,
        )
        st.plotly_chart(fig_pop, use_container_width=True)

    with col2:
        st.subheader("🩸 Кинетика инсулина (Портальный vs Системный)")
        fig_ins = go.Figure()
        fig_ins.add_trace(go.Scatter(
            x=t_30_days, y=site_res["ins_portal"],
            name="Портальный инсулин (печень)",
            line=dict(color="#F59E0B", width=3),
        ))
        fig_ins.add_trace(go.Scatter(
            x=t_30_days, y=site_res["ins_systemic"],
            name="Системный инсулин (периферия)",
            line=dict(color="#EC4899", width=3, dash="dash"),
        ))
        fig_ins.update_layout(
            template="plotly_dark",
            xaxis_title="Время (сутки)", yaxis_title="Инсулин (отн. ед.)",
            height=350,
        )
        st.plotly_chart(fig_ins, use_container_width=True)

    # =====================================================================
    # APOPTOSIS & VASCULARIZATION
    # =====================================================================
    if ap1903_conc > 0.0:
        st.subheader("⚡ Динамика аварийного апоптоза (iCasp9 / AP1903)")
        from organoid_simulator import simulate_icasp9_ap1903_apoptosis
        t_apop_hours = np.linspace(0, 6, 100)
        surv_cells, elim_pct = simulate_icasp9_ap1903_apoptosis(
            t_apop_hours, ap1903_conc_nM=ap1903_conc, N_0=N_0,
        )
        fig_apop = go.Figure()
        fig_apop.add_trace(go.Scatter(
            x=t_apop_hours, y=surv_cells,
            name="Остаточные клетки (млн)",
            line=dict(color="#EF4444", width=4),
        ))
        fig_apop.update_layout(
            template="plotly_dark",
            xaxis_title="Время после введения AP1903 (часы)",
            yaxis_title="Популяция клеток (млн)",
            height=320,
        )
        st.plotly_chart(fig_apop, use_container_width=True)

    st.write("---")
    col_ox_desc, col_ox_plot = st.columns([1, 2])
    with col_ox_desc:
        st.subheader("🫀 Самосборка микрососудистой сети")
        st.markdown(f"""
        Ко-инкапсулированные EPCs и MSCs формируют капилляры непосредственно **внутри** органоида.

        *   **Доля EPCs:** {phi_epc * 100:.1f}%
        *   **Внутреннее pO₂:** {pO2_timeline[-1]:.1f} mmHg
        *   **Время васкуляризации:** {np.log(2) / (0.1 * (1.0 + 5.0 * phi_epc)):.1f} суток
        """)
    with col_ox_plot:
        fig_ox = go.Figure()
        fig_ox.add_trace(go.Scatter(
            x=t_30_days, y=pO2_timeline,
            name="pO₂ (mmHg)",
            line=dict(color="#10B981", width=4),
            fill="tozeroy", fillcolor="rgba(16, 185, 129, 0.1)",
        ))
        fig_ox.update_layout(
            template="plotly_dark",
            xaxis_title="Время (сутки)", yaxis_title="pO₂ (mmHg)",
            height=280,
        )
        st.plotly_chart(fig_ox, use_container_width=True)

    # =====================================================================
    # CLINICAL CALCULATOR
    # =====================================================================
    st.write("---")
    st.header("🏥 Клинический трансплантационный калькулятор для пациента СД1")
    st.markdown("Калькулятор преобразует физиометрические данные больного СД1 в точный расчет дозы IEQ, числа органоидов и площади хирургического сальника.")

    col_pat1, col_pat2, col_pat3 = st.columns(3)
    with col_pat1:
        p_weight = st.slider("Масса тела пациента (кг)", 40.0, 120.0, 70.0, step=1.0)
    with col_pat2:
        p_tdi = st.slider("Суточная доза инсулина TDI (ЕД/сутки)", 10.0, 100.0, 45.0, step=1.0)
    with col_pat3:
        p_cpeptide = st.slider(
            "Базальный C-пептид (пмоль/л)", 0.0, 300.0, 10.0, step=10.0,
            help="< 30 пмоль/л означает абсолютную бета-клеточную недостаточность.",
        )

    p_dose = calculate_patient_transplant_dose(
        weight_kg=p_weight, tdi_units=p_tdi,
        c_peptide_pmol_l=p_cpeptide,
    )

    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Целевая доза (IEQ)", f"{p_dose['total_ieq']:,.0f}",
              delta=f"{p_dose['target_ieq_per_kg']:.0f} IEQ/кг")
    c2.metric("Количество β-клеток", f"{p_dose['total_cells_millions']:.1f} млн",
              delta="1,560 клеток / IEQ")
    c3.metric("Число органоидов", f"{p_dose['total_organoids_count']:,} шт",
              delta="R = 125 мкм")
    c4.metric("Объем матрикса сальника", f"{p_dose['matrix_volume_ml']:.1f} мл",
              delta="Фибрин/dECM")
    c5.metric(
        "Прогноз инсулинонезависимости",
        f"{p_dose['insulin_independence_forecast']:.0f}%",
        delta="Полный отказ от инъекций" if p_dose['insulin_independence_forecast'] >= 100 else "Частичный",
    )

    # OGTT
    st.subheader("🍽️ Динамика гликемического контроля при приеме пищи (OGTT / 4 часа)")
    carb_load = st.slider("Углеводная нагрузка (грамм углеводов)", 20.0, 100.0, 50.0, step=5.0)
    t_ogtt_hours = np.linspace(0, 4, 200)
    g_pre = simulate_ogtt_glycemic_control(t_ogtt_hours, meal_carbs_g=carb_load, is_transplanted=False)
    g_post = simulate_ogtt_glycemic_control(t_ogtt_hours, meal_carbs_g=carb_load, is_transplanted=True)

    fig_ogtt = go.Figure()
    fig_ogtt.add_trace(go.Scatter(
        x=t_ogtt_hours * 60.0, y=g_pre,
        name="ДО трансплантации (СД1 без инсулина)",
        line=dict(color="#EF4444", width=3, dash="dash"),
    ))
    fig_ogtt.add_trace(go.Scatter(
        x=t_ogtt_hours * 60.0, y=g_post,
        name="ПОСЛЕ трансплантации органоидов",
        line=dict(color="#10B981", width=4),
    ))
    fig_ogtt.add_hrect(
        y0=3.9, y1=7.8,
        fillcolor="#10B981", opacity=0.1, line_width=0,
        annotation_text="Целевой физиологический диапазон (3.9 - 7.8 ммоль/л)",
        annotation_position="top left",
    )
    fig_ogtt.add_hrect(
        y0=11.0, y1=25.0,
        fillcolor="#EF4444", opacity=0.1, line_width=0,
        annotation_text="Зона тяжелой гипергликемии",
        annotation_position="top left",
    )
    fig_ogtt.update_layout(
        template="plotly_dark",
        xaxis_title="Время после приема пищи (минуты)",
        yaxis_title="Глюкоза крови (ммоль/л)",
        yaxis=dict(range=[2.0, 24.0]),
        height=400,
        legend=dict(x=0.02, y=0.98),
    )
    st.plotly_chart(fig_ogtt, use_container_width=True)

    # =====================================================================
    # EXPORT
    # =====================================================================
    st.write("---")
    st.header("📄 Экспорт документов и 3D-моделей для клиники")
    col_exp1, col_exp2 = st.columns(2)

    with col_exp1:
        st.subheader("📋 Паспорт трансплантата пациента")
        pat_info = {"weight_kg": p_weight, "tdi_units": p_tdi, "c_peptide_pmol_l": p_cpeptide}
        passport_text = generate_patient_clinical_passport(pat_info, p_dose, {})
        st.download_button(
            label="💾 Скачать Клинический Паспорт Пациента (.md)",
            data=passport_text,
            file_name=f"Clinical_Passport_T1D_Patient_{int(p_weight)}kg.md",
            mime="text/markdown",
        )

    with col_exp2:
        st.subheader("🖨️ 3D-Модель скаффолда для биопринтера")
        stl_content = generate_omental_scaffold_stl(
            area_cm2=p_dose['omental_area_coverage_cm2'],
            thickness_mm=0.5,
        )
        st.download_button(
            label="📐 Скачать 3D-модель скаффолда сальника (.stl)",
            data=stl_content,
            file_name=f"Omental_Scaffold_{int(p_dose['omental_area_coverage_cm2'])}sqcm.stl",
            mime="model/stl",
        )

    # =====================================================================
    # AI RISK STRATIFIER
    # =====================================================================
    st.write("---")
    st.header("🧠 ИИ-Стратификатор клинических рисков и 30-дневный CGM-профиль")
    st.markdown(
        "ИИ-анализатор рассчитывает индекс лабильности гликемии, риски гипогликемических "
        "состояний и прогнозирует 30-дневные стандартизированные метрики CGM."
    )

    col_r1, col_r2 = st.columns([1, 2])

    with col_r1:
        st.subheader("🎯 Приоритетность трансплантации")
        p_hba1c = st.slider("Текущий HbA1c (%)", 6.0, 14.0, 8.8, step=0.1)
        p_hypo_events = st.slider("Тяжелых гипогликемий в месяц", 0, 20, 6, step=1)

        risk_profile = evaluate_patient_clinical_risk_profile(
            tdi_units=p_tdi, c_peptide_pmol_l=p_cpeptide,
            hba1c_percent=p_hba1c, hypo_events_per_month=p_hypo_events,
        )

        st.metric(
            "Индекс приоритетности",
            f"{risk_profile['priority_score']:.0f} / 100",
            delta="Критический" if risk_profile['priority_score'] > 75 else "Высокий",
        )
        st.info(f"**Рекомендация ИИ:** {risk_profile['recommendation']}")
        st.caption(f"**Обоснование:** {risk_profile['rationale']}")

    with col_r2:
        st.subheader("📊 Прогноз 30-дневных показателей CGM (TIR / TBR / TAR)")
        cgm_pre = simulate_cgm_30day_metrics(is_transplanted=False)
        cgm_post = simulate_cgm_30day_metrics(is_transplanted=True)

        cm1, cm2, cm3, cm4 = st.columns(4)
        cm1.metric(
            "Время в норме TIR (3.9-10.0)",
            f"{cgm_post['TIR_percent']:.1f}%",
            delta=f"+{cgm_post['TIR_percent'] - cgm_pre['TIR_percent']:.1f}% (Цель >70%)",
        )
        cm2.metric(
            "Гипогликемии TBR (<3.9)",
            f"{cgm_post['TBR_percent']:.1f}%",
            delta=f"{cgm_post['TBR_percent'] - cgm_pre['TBR_percent']:.1f}% (Цель <4%)",
            delta_color="normal",
        )
        cm3.metric(
            "Гипергликемии TAR (>10.0)",
            f"{cgm_post['TAR_percent']:.1f}%",
            delta=f"{cgm_post['TAR_percent'] - cgm_pre['TAR_percent']:.1f}% (Цель <25%)",
            delta_color="normal",
        )
        cm4.metric(
            "Гликированный HbA1c",
            f"{cgm_post['gmi_hba1c_percent']:.1f}%",
            delta=f"{cgm_post['gmi_hba1c_percent'] - p_hba1c:.1f}% (Норма <6.0%)",
            delta_color="normal",
        )

        # Comparative bar chart
        fig_cgm = go.Figure()
        fig_cgm.add_trace(go.Bar(
            name="ДО трансплантации (Инсулинотерапия)",
            x=["TIR (Норма)", "TBR (Гипо)", "TAR (Гипер)"],
            y=[cgm_pre['TIR_percent'], cgm_pre['TBR_percent'], cgm_pre['TAR_percent']],
            marker_color=["#3B82F6", "#EF4444", "#F59E0B"],
        ))
        fig_cgm.add_trace(go.Bar(
            name="ПОСЛЕ трансплантации органоидов",
            x=["TIR (Норма)", "TBR (Гипо)", "TAR (Гипер)"],
            y=[cgm_post['TIR_percent'], cgm_post['TBR_percent'], cgm_post['TAR_percent']],
            marker_color=["#10B981", "#10B981", "#10B981"],
        ))
        fig_cgm.update_layout(
            template="plotly_dark",
            barmode="group",
            yaxis_title="Процент времени (%)",
            height=320,
            margin=dict(l=20, r=20, t=20, b=20),
        )
        st.plotly_chart(fig_cgm, use_container_width=True)
