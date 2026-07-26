"""
Page setup, sidebar configuration, and mode selection for Streamlit app.
"""
import streamlit as st

from t1d_simulator.simulator import HYDROGELS, IMPLANTATION_SITES


# Mode selection labels (matching the original app.py)
MODE_LABELS = [
    "1D Симуляция диффузии O₂",
    "🔮 Генеративный 3D-дизайн (TPMS)",
    "🧪 ML-подбор антифиброзных покрытий (GNN)",
    "🩸 Неоваскуляризация (VEGF / Ангиогенез)",
    "🧫 Мини-органоиды (Фаза 10: Biomimesis)",
    "📍 Сравнение мест пересадки & Бенчмарки",
]

GEOMETRY_FORMAT = {
    "planar": "Плоская пластина (Slab / Лист)",
    "cylindrical": "Цилиндрическая нить (Fiber / Волокно)",
    "spherical": "Сферическая микрокапсула (Microsphere)",
}

TPMS_TYPE_FORMAT = {
    "gyroid": "Гироид (Gyroid)",
    "schwarz_p": "Поверхность Шварца P (Schwarz P)",
}


def setup_page():
    """Configure Streamlit page settings and custom CSS."""
    st.set_page_config(
        page_title="In Silico Beta-Cell Encapsulation Twin",
        page_icon="🔮",
        layout="wide",
        initial_sidebar_state="expanded",
    )

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


def build_sidebar_mode():
    """Render sidebar mode radio and return selected mode string."""
    st.sidebar.title("🛠️ Выберите режим")
    return st.sidebar.radio(
        "Режим работы:",
        MODE_LABELS,
    )


def format_geometry(x):
    return GEOMETRY_FORMAT.get(x, x)


def format_tpms_type(x):
    return TPMS_TYPE_FORMAT.get(x, x)
