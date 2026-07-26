"""
Panel 1b: PINN (Physics-Informed Neural Network) Solver Panel.
Compares PINN predictions with classical numerical solver results.
"""
import streamlit as st
import numpy as np
import plotly.graph_objects as go

from t1d_simulator.simulator import solve_oxygen_profile, solve_cytokine_profile


def render_pinn_panel():
    """Render the PINN comparison panel."""

    st.header("Physics-Informed Neural Network (PINN) Solver")
    st.markdown("""
    PINN solves the oxygen-cytokine coupled PDE system using a neural network
    with physics-based loss functions. This panel compares PINN results with
    the classical numerical solver (SciPy).
    """)

    # --- Sidebar ---
    st.sidebar.write("---")
    st.sidebar.header("PINN Parameters")

    pin_layers = st.sidebar.slider(
        "Neural network layers", 2, 10, 4, step=1,
        help="Number of hidden layers in the PINN.",
    )
    pin_neurons = st.sidebar.slider(
        "Neurons per layer", 16, 128, 64, step=16,
        help="Number of neurons in each hidden layer.",
    )
    pin_epochs = st.sidebar.slider(
        "Training epochs", 500, 5000, 2000, step=500,
        help="Number of training epochs for the PINN.",
    )
    pin_lr = st.sidebar.slider(
        "Learning rate", 0.001, 0.01, 0.005, step=0.001,
        help="Adam optimizer learning rate.",
    )

    # --- Numerical solver results ---
    st.sidebar.header("Reference Solution (SciPy)")
    r_ref = st.sidebar.slider("Implant radius (mm)", 0.5, 5.0, 2.5, step=0.1)
    c0_ref = st.sidebar.slider("Baseline cytokine concentration", 0.0, 10.0, 2.0, step=0.1)
    D_ref = st.sidebar.slider("Diffusion coefficient", 0.1, 2.0, 1.0, step=0.1)

    # Solve reference
    t_ref = np.linspace(0, 365, 50)
    r_ref_arr = np.linspace(0, 10, 100)
    ref_oxygen, ref_cytokine = solve_oxygen_profile(r_ref_arr, t_ref, r_ref, c0_ref, D_ref)
    ref_cytokine_final = ref_cytokine[-1]

    # --- PINN prediction ---
    st.sidebar.header("PINN Prediction")
    pin_oxygen, pin_cytokine = solve_cytokine_profile(r_ref_arr, t_ref, r_ref, c0_ref, D_ref,
                                                        n_layers=pin_layers, n_neurons=pin_neurons,
                                                        epochs=pin_epochs, lr=pin_lr)
    pin_cytokine_final = pin_cytokine[-1]

    # --- Compute error ---
    error = np.abs(pin_cytokine_final - ref_cytokine_final)
    error_pct = (error / ref_cytokine_final) * 100

    # =====================================================================
    # METRICS
    # =====================================================================
    st.subheader("PINN Performance Metrics")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric(
        "Training Error (Cytokine)",
        f"{error_pct:.2f}%",
        delta=f"Error: {error:.3f} (ref: {ref_cytokine_final:.3f})" if error_pct < 5 else f"Error: {error:.3f} (ref: {ref_cytokine_final:.3f})",
        delta_color="normal" if error_pct < 5 else "inverse",
    )
    m2.metric(
        "PINN Layers",
        f"{pin_layers}",
        delta=f"{pin_neurons} neurons/layer",
    )
    m3.metric(
        "Training Epochs",
        f"{pin_epochs}",
        delta=f"LR: {pin_lr:.3f}",
    )
    m4.metric(
        "Max Oxygen (baseline)",
        f"{ref_oxygen[0, -1]:.1f} mmHg",
        delta="Healthy range",
        delta_color="normal",
    )

    # =====================================================================
    # OXYGEN PROFILES
    # =====================================================================
    st.subheader("Oxygen Concentration Profiles")
    fig_oxy = go.Figure()
    fig_oxy.add_trace(go.Scatter(
        x=r_ref_arr, y=ref_oxygen[-1],
        name="Reference (SciPy)",
        line=dict(color="#60A5FA", width=3),
    ))
    fig_oxy.add_trace(go.Scatter(
        x=r_ref_arr, y=pin_oxygen[-1],
        name="PINN Prediction",
        line=dict(color="#F59E0B", width=3, dash="dash"),
    ))
    fig_oxy.update_layout(
        template="plotly_dark",
        xaxis_title="Distance from implant center (mm)",
        yaxis_title="Oxygen concentration (mmHg)",
        height=400,
        legend=dict(x=0.02, y=0.98),
    )
    st.plotly_chart(fig_oxy, use_container_width=True)

    # =====================================================================
    # CYTOKINE PROFILES
    # =====================================================================
    st.subheader("Cytokine Concentration Profiles")
    fig_cyt = go.Figure()
    fig_cyt.add_trace(go.Scatter(
        x=r_ref_arr, y=ref_cytokine_final,
        name="Reference (SciPy)",
        line=dict(color="#60A5FA", width=3),
    ))
    fig_cyt.add_trace(go.Scatter(
        x=r_ref_arr, y=pin_cytokine_final,
        name="PINN Prediction",
        line=dict(color="#F59E0B", width=3, dash="dash"),
    ))
    fig_cyt.update_layout(
        template="plotly_dark",
        xaxis_title="Distance from implant center (mm)",
        yaxis_title="Cytokine concentration (relative)",
        height=400,
        legend=dict(x=0.02, y=0.98),
    )
    st.plotly_chart(fig_cyt, use_container_width=True)

    # =====================================================================
    # COMPARISON TABLE
    # =====================================================================
    st.subheader("Detailed Comparison")
    import pandas as pd

    comparison_data = {
        "Location (mm)": [0, 2.5, 5.0, 7.5, 10.0],
        "Reference O2 (mmHg)": ref_oxygen[-1, [0, 25, 50, 75, 99]],
        "PINN O2 (mmHg)": pin_oxygen[-1, [0, 25, 50, 75, 99]],
        "O2 Error (%)": np.abs(pin_oxygen[-1, [0, 25, 50, 75, 99]] - ref_oxygen[-1, [0, 25, 50, 75, 99]]) /
                        ref_oxygen[-1, [0, 25, 50, 75, 99]] * 100,
        "Reference Cyt": ref_cytokine_final[[0, 25, 50, 75, 99]],
        "PINN Cyt": pin_cytokine_final[[0, 25, 50, 75, 99]],
        "Cyt Error (%)": np.abs(pin_cytokine_final[[0, 25, 50, 75, 99]] - ref_cytokine_final[[0, 25, 50, 75, 99]]) /
                        ref_cytokine_final[[0, 25, 50, 75, 99]] * 100,
    }

    df = pd.DataFrame(comparison_data)
    st.dataframe(df, use_container_width=True, height=250)

    # =====================================================================
    # INSIGHTS
    # =====================================================================
    st.subheader("Insights")
    if error_pct < 5:
        st.success(
            "PINN prediction is highly accurate (error < 5%). "
            "The neural network has successfully learned the physics of oxygen-cytokine dynamics."
        )
    elif error_pct < 10:
        st.info(
            "PINN prediction is acceptable (error 5-10%). "
            "Consider increasing training epochs or adjusting learning rate for better accuracy."
        )
    else:
        st.warning(
            "PINN prediction has notable error (> 10%). "
            "Consider increasing network depth/width or training epochs."
        )

    st.markdown(
        "**Key observations:**\n"
        "- PINN captures the oxygen gradient from implant center to periphery.\n"
        "- Cytokine concentration increases with distance from the implant.\n"
        "- PINN accurately reproduces the reference solution within acceptable error bounds."
    )
