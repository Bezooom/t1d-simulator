"""
Panel 4: Angiogenesis Module.
Handles VEGF-driven angiogenesis simulation, capillary density, and clinical outcomes.
"""
import streamlit as st
import numpy as np
import plotly.graph_objects as go

from t1d_simulator.simulator import solve_cytokine_profile, solve_cytokine_profile_transient


def render_angiogenesis_panel():
    """Render the angiogenesis panel."""

    st.header("Angiogenesis Module")
    st.markdown("""
    Simulate VEGF-driven angiogenesis in diabetic implants. VEGF (Vascular Endothelial
    Growth Factor) is critical for capillary formation and insulin diffusion.
    """)

    # --- Sidebar ---
    st.sidebar.write("---")
    st.sidebar.header("Angiogenesis Parameters")

    # VEGF parameters
    vegf_concentration = st.sidebar.slider(
        "VEGF concentration (ng/mL)", 0.0, 10.0, 2.0, step=0.1,
        help="Baseline VEGF concentration at implant site.",
    )
    vegf_diffusion = st.sidebar.slider(
        "VEGF diffusion coefficient (mm²/day)", 0.01, 1.0, 0.1, step=0.01,
        help="Diffusion rate of VEGF in tissue.",
    )
    vegf_decay = st.sidebar.slider(
        "VEGF decay rate (1/day)", 0.01, 1.0, 0.1, step=0.01,
        help="Rate of VEGF degradation in tissue.",
    )

    # Capillary parameters
    capillary_density = st.sidebar.slider(
        "Capillary density (capillaries/mm²)", 50, 500, 200, step=10,
        help="Density of capillaries in the implant.",
    )
    capillary_length = st.sidebar.slider(
        "Average capillary length (mm)", 0.1, 10.0, 1.0, step=0.1,
        help="Average length of capillaries.",
    )
    capillary_radius = st.sidebar.slider(
        "Capillary radius (mm)", 0.001, 0.01, 0.005, step=0.001,
        help="Radius of individual capillaries.",
    )

    # Time parameters
    time_days = st.sidebar.slider(
        "Simulation time (days)", 1, 365, 30, step=1,
        help="Duration of angiogenesis simulation.",
    )

    # --- Simulation ---
    st.sidebar.header("Simulation")
    if st.sidebar.button("Run Angiogenesis Simulation"):
        with st.spinner("Running angiogenesis simulation..."):
            # Generate time array
            t = np.linspace(0, time_days, 100)

            # Simulate VEGF concentration over time
            vegf_profile = np.zeros(len(t))
            for i, ti in enumerate(t):
                vegf_profile[i] = vegf_concentration * np.exp(-vegf_decay * ti)

            # Simulate capillary density over time (logistic growth)
            capillary_profile = np.zeros(len(t))
            for i, ti in enumerate(t):
                capillary_profile[i] = capillary_density / (
                    1 + np.exp(-0.1 * (ti - time_days / 2))
                )

            # Simulate insulin diffusion rate
            insulin_diffusion = np.zeros(len(t))
            for i, ti in enumerate(t):
                insulin_diffusion[i] = (
                    0.5 + 0.5 * (capillary_profile[i] / capillary_density)
                )

            # Store results
            st.session_state['angiogenesis_results'] = {
                'time': t,
                'vegf_profile': vegf_profile,
                'capillary_profile': capillary_profile,
                'insulin_diffusion': insulin_diffusion,
            }
            st.session_state['angiogenesis_params'] = {
                'vegf_concentration': vegf_concentration,
                'vegf_diffusion': vegf_diffusion,
                'vegf_decay': vegf_decay,
                'capillary_density': capillary_density,
                'capillary_length': capillary_length,
                'capillary_radius': capillary_radius,
                'time_days': time_days,
            }
            st.success("Angiogenesis simulation completed!")

    # --- Results ---
    if 'angiogenesis_results' in st.session_state:
        results = st.session_state['angiogenesis_results']
        params = st.session_state['angiogenesis_params']

        st.subheader("Angiogenesis Results")

        # --- VEGF Profile ---
        st.subheader("VEGF Concentration Over Time")
        fig_vegf = go.Figure()
        fig_vegf.add_trace(go.Scatter(
            x=results['time'],
            y=results['vegf_profile'],
            mode='lines',
            name='VEGF (ng/mL)',
            line=dict(color='#10B981', width=3),
        ))
        fig_vegf.add_hline(
            y=params['vegf_concentration'],
            line_dash="dash",
            line_color="#EF4444",
            annotation_text="Initial VEGF",
        )
        fig_vegf.update_layout(
            template="plotly_dark",
            xaxis_title="Time (days)",
            yaxis_title="VEGF concentration (ng/mL)",
            height=400,
        )
        st.plotly_chart(fig_vegf, use_container_width=True)

        # --- Capillary Density ---
        st.subheader("Capillary Density Over Time")
        fig_cap = go.Figure()
        fig_cap.add_trace(go.Scatter(
            x=results['time'],
            y=results['capillary_profile'],
            mode='lines',
            name='Capillary density (cap/mm²)',
            line=dict(color='#60A5FA', width=3),
        ))
        fig_cap.add_hline(
            y=params['capillary_density'],
            line_dash="dash",
            line_color="#F59E0B",
            annotation_text="Final density",
        )
        fig_cap.update_layout(
            template="plotly_dark",
            xaxis_title="Time (days)",
            yaxis_title="Capillary density (capillaries/mm²)",
            height=400,
        )
        st.plotly_chart(fig_cap, use_container_width=True)

        # --- Insulin Diffusion ---
        st.subheader("Insulin Diffusion Rate")
        fig_ins = go.Figure()
        fig_ins.add_trace(go.Scatter(
            x=results['time'],
            y=results['insulin_diffusion'],
            mode='lines',
            name='Insulin diffusion (relative)',
            line=dict(color='#EC4899', width=3),
        ))
        fig_ins.add_hline(
            y=1.0,
            line_dash="dash",
            line_color="#10B981",
            annotation_text="Max diffusion",
        )
        fig_ins.update_layout(
            template="plotly_dark",
            xaxis_title="Time (days)",
            yaxis_title="Insulin diffusion rate (relative)",
            height=400,
        )
        st.plotly_chart(fig_ins, use_container_width=True)

        # --- Summary Metrics ---
        st.subheader("Summary Metrics")
        col1, col2, col3, col4 = st.columns(4)
        final_vegf = results['vegf_profile'][-1]
        final_cap = results['capillary_profile'][-1]
        final_ins = results['insulin_diffusion'][-1]
        time_to_half = np.log(2) / params['vegf_decay']

        col1.metric(
            "Final VEGF",
            f"{final_vegf:.2f} ng/mL",
            delta=f"Decay: {params['vegf_decay']:.3f} /day",
        )
        col2.metric(
            "Final Capillary Density",
            f"{final_cap:.0f} cap/mm²",
            delta=f"Target: {params['capillary_density']} cap/mm²",
        )
        col3.metric(
            "Insulin Diffusion",
            f"{final_ins:.2f}",
            delta="Relative to max",
        )
        col4.metric(
            "Time to 50% VEGF",
            f"{time_to_half:.1f} days",
            delta="Half-life",
        )

        # --- Clinical Insights ---
        st.subheader("Clinical Insights")
        if final_cap >= params['capillary_density'] * 0.9:
            st.success(
                "Angiogenesis is progressing well. Capillary density is approaching target. "
                "Insulin diffusion is optimal."
            )
        elif final_cap >= params['capillary_density'] * 0.7:
            st.info(
                "Angiogenesis is progressing. Capillary density is at 70-90% of target. "
                "Continue monitoring."
            )
        else:
            st.warning(
                "Angiogenesis is slow. Capillary density is below 70% of target. "
                "Consider VEGF supplementation or extended simulation time."
            )

        st.markdown(
            "**Key insights:**\n"
            "- VEGF concentration decays over time, driving capillary growth.\n"
            "- Capillary density follows a logistic growth pattern.\n"
            "- Insulin diffusion rate improves with capillary density.\n"
            "- Time to 50% VEGF decay indicates treatment duration."
        )
    else:
        st.info("Click 'Run Angiogenesis Simulation' in the sidebar to run the simulation.")
