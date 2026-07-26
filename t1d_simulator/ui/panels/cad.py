"""
Panel 2: 3D TPMS Design Module.
Handles TPMS surface generation, meshing, and 3D visualization.
"""
import streamlit as st
import numpy as np
import plotly.graph_objects as go
import mesh_generator

from t1d_simulator.simulator import HYDROGELS, IMPLANTATION_SITES, K_M, K_M_INSULIN, V_MAX


def render_cad_panel():
    """Render the 3D TPMS design panel."""

    st.header("3D TPMS Design Module")
    st.markdown("""
    Generate and visualize Triply Periodic Minimal Surface (TPMS) structures
    for diabetic implant encapsulation. TPMS surfaces provide high surface-area
    to volume ratios for cell attachment and nutrient diffusion.
    """)

    # --- Sidebar ---
    st.sidebar.write("---")
    st.sidebar.header("TPMS Parameters")

    # Material selection
    hydrogel = st.sidebar.selectbox(
        "Hydrogel type",
        list(HYDROGELS.keys()),
        index=0,
    )
    hydrogel_props = HYDROGELS[hydrogel]

    # Geometry parameters
    width = st.sidebar.slider("Width (mm)", 5.0, 50.0, 20.0, step=1.0)
    height = st.sidebar.slider("Height (mm)", 5.0, 50.0, 20.0, step=1.0)
    depth = st.sidebar.slider("Depth (mm)", 5.0, 50.0, 20.0, step=1.0)

    # TPMS type
    tpms_type = st.sidebar.selectbox(
        "TPMS surface type",
        ["Schwarz P", "Schwarz D", "Gyroid", "Diamond"],
    )

    # Resolution
    resolution = st.sidebar.slider("Mesh resolution", 20, 100, 50, step=10)

    # --- Generate TPMS mesh ---
    st.sidebar.header("Mesh Generation")
    if st.sidebar.button("Generate TPMS Mesh"):
        with st.spinner("Generating TPMS mesh..."):
            tpms_mesh = mesh_generator.generate_tpms(
                width=width,
                height=height,
                depth=depth,
                tpms_type=tpms_type,
                resolution=resolution,
            )

            # Store in session state
            st.session_state['tpms_mesh'] = tpms_mesh
            st.session_state['tpms_params'] = {
                'hydrogel': hydrogel,
                'width': width,
                'height': height,
                'depth': depth,
                'tpms_type': tpms_type,
                'resolution': resolution,
            }
            st.success("TPMS mesh generated successfully!")

    # --- 3D Visualization ---
    if 'tpms_mesh' in st.session_state:
        tpms_mesh = st.session_state['tpms_mesh']
        tpms_params = st.session_state['tpms_params']

        st.subheader("3D TPMS Structure")

        # Create mesh plot
        fig = go.Figure(data=[
            go.Mesh3d(
                x=tpms_mesh['vertices'][:, 0],
                y=tpms_mesh['vertices'][:, 1],
                z=tpms_mesh['vertices'][:, 2],
                i=tpms_mesh['faces'][:, 0],
                j=tpms_mesh['faces'][:, 1],
                k=tpms_mesh['faces'][:, 2],
                color='#60A5FA',
                opacity=0.6,
                name='TPMS Surface',
                showscale=False,
            )
        ])

        fig.update_layout(
            scene=dict(
                xaxis_title='X (mm)',
                yaxis_title='Y (mm)',
                zaxis_title='Z (mm)',
                aspectmode='data',
            ),
            template="plotly_dark",
            height=600,
        )

        st.plotly_chart(fig, use_container_width=True)

        # --- Mesh Statistics ---
        st.subheader("Mesh Statistics")
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Vertices", f"{tpms_mesh['n_vertices']:,}")
        col2.metric("Faces", f"{tpms_mesh['n_faces']:,}")
        col3.metric("Surface Area (mm²)", f"{tpms_mesh['surface_area']:.1f}")
        col4.metric("Volume (mm³)", f"{tpms_mesh['volume']:.1f}")

        # --- Material Properties ---
        st.subheader("Material Properties")
        col5, col6, col7, col8 = st.columns(4)
        col5.metric("Elastic Modulus", f"{hydrogel_props['elastic_modulus']:.2f} MPa")
        col6.metric("Porosity", f"{hydrogel_props['porosity']:.2%}")
        col7.metric("Swelling Ratio", f"{hydrogel_props['swelling_ratio']:.2f}")
        col8.metric("Biodegradation (days)", f"{hydrogel_props['biodegradation_time']:.0f}")

        # --- Export Options ---
        st.subheader("Export Options")
        col_exp1, col_exp2, col_exp3 = st.columns(3)
        with col_exp1:
            st.download_button(
                label="Download STL",
                data=tpms_mesh['stl_data'],
                file_name=f"tpms_{tpms_params['tpms_type'].lower()}_{width}x{height}x{depth}.stl",
                mime="model/stl",
            )
        with col_exp2:
            st.download_button(
                label="Download JSON",
                data=str(tpms_mesh),
                file_name=f"tpms_{tpms_params['tpms_type'].lower()}.json",
                mime="application/json",
            )
        with col_exp3:
            st.download_button(
                label="Download CSV",
                data=tpms_mesh['csv_data'],
                file_name=f"tpms_vertices_{width}x{height}x{depth}.csv",
                mime="text/csv",
            )
    else:
        st.info("Click 'Generate TPMS Mesh' in the sidebar to create a TPMS structure.")

    # --- Design Guidelines ---
    st.subheader("Design Guidelines")
    st.markdown("""
    **TPMS Surface Types:**
    - **Schwarz P**: Primitive surface, high connectivity, good for diffusion
    - **Schwarz D**: Diamond surface, complex topology, excellent for cell attachment
    - **Gyroid**: Minimal surface, triply periodic, balanced properties
    - **Diamond**: High surface area, ideal for tissue engineering

    **Key Design Parameters:**
    - **Resolution**: Higher resolution = smoother surface but more computational cost
    - **Dimensions**: Match to implant site (pancreas, omentum, etc.)
    - **Material**: Choose hydrogel based on mechanical and biological requirements
    """)
