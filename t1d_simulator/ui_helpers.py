# -*- coding: utf-8 -*-
"""
UI Helpers for Streamlit App (t1d_simulator/app.py)
Provides modular UI rendering components for site decision matrix and benchmark validation.
"""

import streamlit as st
import pandas as pd
import numpy as np

def render_site_decision_matrix():
    """
    Renders the Transplantation Site Decision Matrix in Streamlit.
    """
    st.subheader("📊 Сравнительная матрица мест имплантации (Decision Matrix)")
    st.markdown("""
    Сравнение анатомических зон пересадки по мультифизическим, биологическим и клиническим критериям на основе `docs/site_comparison_matrix.md`.
    """)
    
    matrix_data = {
        "Критерий оценки": [
            "Базовое pO₂ (mmHg)",
            "Риск острой реакции IBMIR (%)",
            "Прогноз выживаемости (48ч, %)",
            "Свойство хирургической извлекаемости",
            "Риск жирового гепатоза (Steatosis)",
            "Кинетика васкуляризации",
            "Интегральный индекс клиники (0-100)"
        ],
        "Сальник (Omental Pouch)": [
            "55.0 mmHg (Высокое)",
            "10% (Низкий)",
            "90 - 95%",
            "100% Извлекаем (Безопасно)",
            "5% (Минимальный)",
            "Быстрая (богатая сеть)",
            "92 / 100 ⭐"
        ],
        "Печень (Воротная вена)": [
            "40.0 mmHg (Умеренное)",
            "100% (Тяжелый IBMIR)",
            "40 - 65%",
            "0% (Неизвлекаем!)",
            "85% (Высокий риск)",
            "Мгновенная (кровь)",
            "55 / 100"
        ],
        "Подкожно (SQ Tissue)": [
            "30.0 mmHg (Низкое)",
            "0% (Отсутствует)",
            "60 - 80%",
            "100% Извлекаем",
            "0% (Отсутствует)",
            "Медленная (14-21 день)",
            "45 / 100"
        ],
        "Артериальный AV-луп": [
            "95.0 mmHg (Идеальное)",
            "30-50% (Средний)",
            "75 - 90%",
            "100% Извлекаем",
            "0% (Отсутствует)",
            "Мгновенная (перфузия)",
            "78 / 100"
        ]
    }
    
    df_matrix = pd.DataFrame(matrix_data)
    st.table(df_matrix)
    st.success("💡 **Рекомендация модели**: Большой Сальник (Omental Pouch) обеспечивает максимальный баланс оксигенации (55 mmHg), безопасности (100% извлекаемость) и минимизации IBMIR.")

def render_benchmark_validation_summary():
    """
    Renders Literature Benchmark Validation Summary.
    """
    st.subheader("📚 Валидация модели на литературных бенчмарках")
    st.markdown("""
    Сравнение предсказаний симулятора с экспериментальными данными из рецензируемых работ **Papas et al. (2007)** и **Papabathini et al. (2023)**.
    """)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### Papas et al. 2007 (Сферические микрокапсулы, R=400 мкм)")
        df_papas = pd.DataFrame({
            "Плотность (млн/мл)": [20, 50, 80, 150],
            "Эксперимент (%)": [92.0, 68.0, 45.0, 22.0],
            "Модель (%)": [89.0, 64.8, 53.5, 40.5],
            "Ошибка (%)": [3.0, 3.2, 8.5, 18.5]
        })
        st.dataframe(df_papas, use_container_width=True)
        st.caption("RMSE: **10.42%** | $pO_2 = 28\ \text{mmHg}$ static culture")
        
    with col2:
        st.markdown("#### Papabathini et al. 2023 (Плоский лист, L=250 мкм)")
        df_papa = pd.DataFrame({
            "Плотность (млн/мл)": [10, 30, 50, 100],
            "Эксперимент (%)": [85.0, 52.0, 32.0, 15.0],
            "Модель (%)": [100.0, 53.7, 41.6, 29.2],
            "Ошибка (%)": [15.0, 1.7, 9.6, 14.2]
        })
        st.dataframe(df_papa, use_container_width=True)
        st.caption("RMSE: **11.42%** | $pO_2 = 28\ \text{mmHg}$ static culture")

def render_aid_closed_loop_dashboard():
    """
    Renders Automated Insulin Delivery (AID) Closed-Loop Artificial Pancreas Dashboard.
    """
    import plotly.graph_objects as gr_obj
    from plotly.subplots import make_subplots
    from aid_controller import simulate_aid_closed_loop
    
    st.subheader("📟 Моделирование Искусственной Поджелудочной Железы (AID Closed-Loop)")
    st.markdown("""
    Симуляция замкнутого контура помповой инсулинотерапии (PID / Bergman Minimal Model) под контролем CGM.
    """)
    
    col1, col2 = st.columns(2)
    with col1:
        meal_carbs = st.slider("Углеводная нагрузка (г углеводов)", 20.0, 120.0, 60.0, step=5.0)
        target_g = st.slider("Целевая гликемия (мг/дл)", 90.0, 140.0, 110.0, step=5.0)
    with col2:
        kp = st.slider("Пропорциональный коэффициент Kp", 0.01, 0.20, 0.08, step=0.01)
        kd = st.slider("Дифференциальный коэффициент Kd", 0.05, 0.40, 0.15, step=0.01)
        
    aid_res = simulate_aid_closed_loop(
        meal_carbs_g=meal_carbs,
        target_glucose_mg_dl=target_g,
        kp=kp,
        kd=kd
    )
    
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("TIR (3.9 - 10.0 ммоль/л)", f"{aid_res['TIR_percent']:.1f}%", delta="Цель > 70%")
    m2.metric("TBR (< 3.9 ммоль/л)", f"{aid_res['TBR_percent']:.1f}%", delta="Цель < 4%")
    m3.metric("TAR (> 10.0 ммоль/л)", f"{aid_res['TAR_percent']:.1f}%", delta="Цель < 25%")
    m4.metric("Средняя гликемия", f"{aid_res['mean_glucose_mmol_l']:.1f} ммоль/л")
    
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(
        gr_obj.Scatter(
            x=aid_res["time_hours"],
            y=aid_res["glucose_mmol_l"],
            name="CGM Гликемия (ммоль/л)",
            line=dict(color="#60A5FA", width=3)
        ),
        secondary_y=False
    )
    fig.add_trace(
        gr_obj.Scatter(
            x=aid_res["time_hours"],
            y=aid_res["infusion_rate_u_h"],
            name="Инфузия помпы (Ед/ч)",
            line=dict(color="#10B981", width=2, dash="dash")
        ),
        secondary_y=True
    )
    fig.add_hrect(
        y0=3.9, y1=10.0,
        fillcolor="#10B981", opacity=0.1, line_width=0,
        annotation_text="Целевой диапазон TIR",
        annotation_position="top left"
    )
    fig.update_layout(
        template="plotly_dark",
        height=400,
        xaxis_title="Время (часы)",
        margin=dict(l=20, r=20, t=20, b=20)
    )
    fig.update_yaxes(title_text="Глюкоза (ммоль/л)", secondary_y=False)
    fig.update_yaxes(title_text="Инсулин (Ед/ч)", secondary_y=True)
    st.plotly_chart(fig, use_container_width=True)

def render_screening_and_tornado_dashboard():
    """
    Renders interactive Construct Screening (WP-Screen) and Tornado Sensitivity dashboard.
    """
    import plotly.express as px
    from screen_design import screen_design, export_markdown_report
    from uncertainty_analysis import run_tornado_sensitivity

    st.header("🔍 Автоматический скрининг отторжения и Торнадо-анализ (WP-Screen)")
    st.markdown("""
    Инструмент экспресс-оценки мультифизических дефектов конструкт-дизайна (Крог-гипоксия, острая фаза IBMIR, фиброз) и Торнадо-анализ чувствительности ключевых параметров.
    """)

    col1, col2 = st.columns([1, 1])

    with col1:
        st.subheader("📐 Параметры трансплантата")
        geom = st.selectbox("Геометрия геля:", ["spherical", "planar", "cylindrical"], format_func=lambda x: {"spherical": "Сфера (Sphere)", "planar": "Лист (Planar)", "cylindrical": "Волокно (Fiber)"}[x])
        rad = st.slider("Радиус / Полутолщина (мкм):", 50.0, 500.0, 200.0, 25.0)
        dens = st.slider("Плотность загрузки (млн клеток/мл):", 10.0, 200.0, 80.0, 10.0)
        site = st.selectbox("Сайт имплантации:", ["omental_pouch", "portal_vein", "subcutaneous", "av_loop"], format_func=lambda x: {"omental_pouch": "Большой сальник (Omentum)", "portal_vein": "Воротная вена (Portal)", "subcutaneous": "Подкожная клетчатка (SQ)", "av_loop": "Артериальный AV-loop"}[x])
        heparin = st.checkbox("Гепаринизированная поверхность (Lipid-PEG-LMWH)", value=True)
        fibrosis = st.slider("Толщина фиброзной капсулы L_fib (мкм):", 0.0, 100.0, 0.0, 5.0)

    res = screen_design(geometry=geom, radius_microns=rad, density_million_per_ml=dens, site_key=site, heparin_coated=heparin, l_fibrosis_microns=fibrosis)

    with col2:
        st.subheader("📊 Экспресс-вердикт мультифизики")
        st.markdown(f"### Статус: **{res['status_badge']}**")
        st.info(f"**Рекомендация:** {res['recommendation']}")

        m1, m2 = st.columns(2)
        m1.metric("Выживаемость клеток", f"{res['viable_fraction_percent']:.1f}%", delta="> 70.0% Цель")
        m2.metric("Мин. pO₂ в ядре", f"{res['min_core_pO2_mmHg']:.2f} mmHg", delta="> 0.5 mmHg Крог")

        m3, m4 = st.columns(2)
        m3.metric("Ретенция IBMIR (48ч)", f"{res['ibmir_48h_retention_percent']:.1f}%", delta="> 90.0% Цель")
        m4.metric("Инсулиновая емкость", f"{res['insulin_secretion_percent']:.1f}%")

    st.write("---")
    st.subheader("🌪️ Торнадо-анализ чувствительности параметров (±20% вариация)")
    tornado_res = run_tornado_sensitivity(base_radius_microns=rad, base_density_million_per_ml=dens, geometry=geom, variation_pct=20.0)
    df_torn = tornado_res["dataframe"]

    fig_torn = px.bar(
        df_torn,
        y="parameter",
        x="delta_viability",
        orientation="h",
        title=f"Чувствительность выживаемости к изменению параметров (Базовый уровень: {tornado_res['base_viability']:.1f}%)",
        labels={"delta_viability": "Сдвиг выживаемости Δf_viab (%)", "parameter": "Физическая константа"},
        color="delta_viability",
        color_continuous_scale="Reds"
    )
    fig_torn.update_layout(template="plotly_dark", height=350, margin=dict(l=20, r=20, t=40, b=20))
    st.plotly_chart(fig_torn, use_container_width=True)

    with st.expander("📄 Сгенерированный текстовый отчет (Markdown)"):
        st.code(export_markdown_report(res), language="markdown")

