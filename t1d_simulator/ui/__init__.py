"""UI package for t1d_simulator."""
from t1d_simulator.ui.layout import setup_page, build_sidebar_mode
from t1d_simulator.ui.helpers import (
    render_site_decision_matrix,
    render_benchmark_validation_summary,
    render_aid_closed_loop_dashboard,
)

__all__ = [
    "setup_page",
    "build_sidebar_mode",
    "render_site_decision_matrix",
    "render_benchmark_validation_summary",
    "render_aid_closed_loop_dashboard",
]
