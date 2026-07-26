"""
Main Streamlit application for T1D Simulator.
Orchestrates all 6 simulation modes using modular panels.
"""
import streamlit as st
import numpy as np

from t1d_simulator.ui.panels import (
    render_oxygen_panel,
    render_cad_panel,
    render_gnn_panel,
    render_angiogenesis_panel,
    render_organoid_panel,
    render_pinn_panel,
)

# Настройка страницы
st.set_page_config(
    page_title="T1D Simulator - Beta-Cell Encapsulation Twin",
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
    .main-header {
        font-size: 2.5rem;
        font-weight: 700;
        margin-bottom: 0.5rem;
    }
    .main-subheader {
        font-size: 1.1rem;
        color: #9CA3AF;
        margin-bottom: 2rem;
    }
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-header">🔮 Цифровой двойник капсулы β-клеток</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="main-subheader">Simulator solves nonlinear oxygen reaction-diffusion equation '
    'in various geometries and procedurally generates 3D TPMS scaffolds for bioprinting optimization.</div>',
    unsafe_allow_html=True
)

# --- САЙДБАР: РЕЖИМ РАБОТЫ ---
st.sidebar.title("🛠️ Select Mode")
app_mode = st.sidebar.radio(
    "Mode:",
    [
        "1️⃣ 1D Oxygen Diffusion",
        "🔮 3D TPMS Design",
        "🧪 GNN Drug Screening",
        "🩸 Angiogenesis",
        "🧫 Mini-Organoids (Phase 10)",
        "📊 PINN Comparison",
    ],
    index=0,
)

# ==============================================================================
# RENDER PANELS BASED ON MODE
# ==============================================================================
if app_mode == "1️⃣ 1D Oxygen Diffusion":
    render_oxygen_panel()
elif app_mode == "🔮 3D TPMS Design":
    render_cad_panel()
elif app_mode == "🧪 GNN Drug Screening":
    render_gnn_panel()
elif app_mode == "🩸 Angiogenesis":
    render_angiogenesis_panel()
elif app_mode == "🧫 Mini-Organoids (Phase 10)":
    render_organoid_panel()
elif app_mode == "📊 PINN Comparison":
    render_pinn_panel()

# --- Footer ---
st.markdown("---")
st.caption("T1D Simulator v1.0 | Beta-Cell Encapsulation Twin | Physics-Informed Modeling")
