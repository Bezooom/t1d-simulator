"""UI panels package for t1d_simulator."""
from t1d_simulator.ui.panels.oxygen import render_oxygen_panel
from t1d_simulator.ui.panels.organoid import render_organoid_panel
from t1d_simulator.ui.panels.gnn import render_gnn_panel
from t1d_simulator.ui.panels.pinn import render_pinn_panel
from t1d_simulator.ui.panels.cad import render_cad_panel
from t1d_simulator.ui.panels.angiogenesis import render_angiogenesis_panel

__all__ = [
    "render_oxygen_panel",
    "render_organoid_panel",
    "render_gnn_panel",
    "render_pinn_panel",
    "render_cad_panel",
    "render_angiogenesis_panel",
]
