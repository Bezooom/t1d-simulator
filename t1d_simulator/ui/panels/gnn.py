"""
Panel 3: GNN Drug-Screening Module (Phase 9).
Handles drug-gene interaction graph, molecular features, GNN predictions,
and drug ranking with clinical translation.
"""
import streamlit as st
import numpy as np
import plotly.graph_objects as go

from t1d_simulator.gnn_screening import (
    build_drug_gene_graph,
    compute_drug_similarity,
    compute_drug_gene_interaction,
    predict_therapeutic_potency,
    compute_immune_evasion_score,
    rank_drugs_by_therapeutic_potency,
    compute_drug_immune_therapeutic_index,
    simulate_immune_evasion_kinetics,
)


def render_gnn_panel():
    """Render the full Phase 9 GNN Drug-Screening panel."""

    st.header("GNN Drug-Screening Module (Phase 9)")
    st.markdown(
        "This module uses a Graph Neural Network (GNN) for predictive evaluation "
        "of therapeutic efficacy and immunological safety of drug combinations "
        "for patients with Type 1 Diabetes."
    )

    # --- Sidebar ---
    st.sidebar.write("---")
    st.sidebar.header("GNN Screening Parameters")

    drug_pool = st.sidebar.slider(
        "Drug pool size (N)",
        20, 200, 50, step=10,
        help="Number of drug candidates in the GNN screening pool.",
    )

    gnn_depth = st.sidebar.slider(
        "GNN depth (layers)",
        2, 8, 4, step=1,
        help="Number of graph convolution layers. More layers = more context.",
    )

    gnn_dropout = st.sidebar.slider(
        "GNN dropout",
        0.0, 0.5, 0.2, step=0.05,
        help="Regularization: prevents overfitting.",
    )

    drug_families = st.sidebar.multiselect(
        "Drug families to screen:",
        ["Biguanides (Metformin)",
         "Thiazolidinediones (Glitazones)",
         "SGLT2 Inhibitors",
         "GLP-1 Receptor Agonists",
         "DPP-4 Inhibitors",
         "Immunomodulators (Anti-TNF, Anti-CD20)",
         "Antioxidants (NAC, Alpha-Lipoic)",
         "Growth Factors (VEGF, FGF, IGF-1)",
         "Corticosteroids (Prednisone)",
         "Stem Cell Activators (SDF-1, CXCL12)"],
        default=["Biguanides (Metformin)", "Immunomodulators (Anti-TNF, Anti-CD20)"],
        help="Selected families will be included in the GNN graph.",
    )

    # --- Build GNN graph ---
    st.sidebar.header("Data preprocessing")
    n_genes = st.sidebar.slider("Number of target genes", 50, 500, 200, step=10)
    n_drugs = st.sidebar.slider("Number of drugs in pool", 20, 150, 50, step=10)
    n_edges = st.sidebar.slider("Drug-gene edges (avg/drug)", 5, 30, 15, step=1)

    graph_data = build_drug_gene_graph(
        n_genes=n_genes, n_drugs=n_drugs, avg_edges_per_drug=n_edges,
        drug_families=drug_families,
    )

    # --- Compute interactions ---
    drug_features = compute_drug_similarity(graph_data, n_genes=n_genes)
    drug_interactions = compute_drug_gene_interaction(
        drug_features, graph_data, n_genes=n_genes,
    )

    # --- GNN prediction ---
    st.sidebar.header("GNN Prediction")
    gnn_batch_size = st.sidebar.slider("GNN batch size", 8, 64, 32, step=8)

    gnn_predictions = predict_therapeutic_potency(
        drug_features, drug_interactions, graph_data,
        n_genes=n_genes, n_drugs=n_drugs,
        gnn_layers=gnn_depth, gnn_dropout=gnn_dropout,
        batch_size=gnn_batch_size,
    )

    # --- Immune evasion ---
    immune_scores = compute_immune_evasion_score(
        drug_features, drug_interactions, gnn_predictions,
        n_genes=n_genes, n_drugs=n_drugs,
    )

    # --- Drug ranking ---
    ranked_drugs = rank_drugs_by_therapeutic_potency(
        gnn_predictions, immune_scores, drug_features,
        n_genes=n_genes, n_drugs=n_drugs,
    )

    immune_index = compute_drug_immune_therapeutic_index(
        gnn_predictions, immune_scores, drug_features,
        n_genes=n_genes, n_drugs=n_drugs,
    )

    # =====================================================================
    # METRICS
    # =====================================================================
    st.subheader("GNN Screening Key Metrics")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric(
        "GNN Efficiency (AUC)",
        f"{gnn_predictions['gnn_auc']:.3f}",
        delta="Prediction accuracy",
    )
    m2.metric(
        "Drug-gene edges",
        f"{graph_data['n_edges']:,}",
        delta=f"From {n_genes} genes x {n_drugs} drugs",
    )
    m3.metric(
        "Mean immune index",
        f"{immune_scores['mean_immune_therapeutic_index']:.3f}",
        delta="> 0.7 - high",
        delta_color="normal",
    )
    m4.metric(
        "Immune toxicity risk",
        f"{immune_scores['immune_toxicity_risk_percent']:.1f}%",
        delta="Low (optimized)",
        delta_color="normal",
    )

    # =====================================================================
    # DRUG SIMILARITY MATRIX
    # =====================================================================
    st.subheader("Drug Similarity Matrix (Cosine Similarity)")
    fig_sim = go.Figure()
    fig_sim.add_trace(go.Heatmap(
        z=drug_features['cosine_similarity_matrix'],
        x=[f"D{i+1}" for i in range(min(n_drugs, 30))],
        y=[f"D{i+1}" for i in range(min(n_drugs, 30))],
        colorscale="RdYlBu_r",
        reversescale=True,
        colorbar=dict(title="Similarity"),
    ))
    fig_sim.update_layout(
        template="plotly_dark",
        height=500,
        xaxis=dict(showgrid=False),
        yaxis=dict(showgrid=False, autorange="reversed"),
    )
    st.plotly_chart(fig_sim, use_container_width=True)

    # =====================================================================
    # GNN-GRAPH VISUALIZATION
    # =====================================================================
    st.subheader("Drug-Gene Interaction Graph")
    col_g1, col_g2 = st.columns([2, 1])
    with col_g1:
        fig_gnn = go.Figure()
        for i, drug in enumerate(graph_data['drug_nodes'][:30]):
            fig_gnn.add_trace(go.Scatter(
                x=[drug['x']], y=[drug['y']],
                mode="markers+text",
                name=drug['name'],
                marker=dict(
                    size=15, color="#60A5FA",
                    line=dict(width=2, color="#1E3A5F"),
                ),
                text=drug['name'],
                textposition="top center",
                showlegend=False,
            ))
        for edge in graph_data['edges'][:50]:
            x0, y0 = edge['source']['x'], edge['source']['y']
            x1, y1 = edge['target']['x'], edge['target']['y']
            fig_gnn.add_trace(go.Scatter(
                x=[x0, x1], y=[y0, y1],
                mode="lines",
                line=dict(color="#4B5563", width=1),
                showlegend=False,
            ))
        fig_gnn.update_layout(
            template="plotly_dark",
            showlegend=False,
            height=500,
            xaxis=dict(showgrid=False, showticklabels=False, title=""),
            yaxis=dict(showgrid=False, showticklabels=False, title=""),
        )
        st.plotly_chart(fig_gnn, use_container_width=True)

    with col_g2:
        st.subheader("Top-5 drugs by therapeutic index")
        top5 = ranked_drugs[:5]
        for i, drug in enumerate(top5, 1):
            st.markdown(
                f"**{i}. {drug['name']}**\n"
                f"- TPI: `{drug['therapeutic_index']:.3f}`\n"
                f"- Potential: `{drug['therapeutic_potential']:.3f}`\n"
                f"- Immune toxicity: `{drug['immune_toxicity_percent']:.1f}%`\n"
                f"- Family: `{drug['family']}`"
            )

    # =====================================================================
    # THERAPEUTIC POTENTIAL vs IMMUNE TOXICITY
    # =====================================================================
    st.subheader("Therapeutic Potential vs Immunological Safety")
    fig_tpi = go.Figure()
    fig_tpi.add_trace(go.Scatter(
        x=immune_index['therapeutic_potential'],
        y=immune_index['immune_therapeutic_index'],
        mode="markers",
        name="Drugs",
        marker=dict(
            size=12,
            color=immune_index['immune_therapeutic_index'],
            colorscale="Viridis",
            colorbar=dict(title="TPI"),
            line=dict(width=1, color="white"),
        ),
        text=[d['name'] for d in ranked_drugs],
        hoverinfo="text",
        hovertemplate="%{text}<br>TPI: %{y:.3f}<br>Potential: %{x:.3f}",
        showlegend=False,
    ))
    fig_tpi.add_hline(
        y=0.7, line_dash="dash", line_color="#10B981",
        annotation_text="Target TPI > 0.7",
        annotation_position="bottom right",
    )
    fig_tpi.add_vline(
        x=0.7, line_dash="dash", line_color="#60A5FA",
        annotation_text="Target Potential > 0.7",
        annotation_position="top left",
    )
    fig_tpi.update_layout(
        template="plotly_dark",
        xaxis_title="Therapeutic Potential (GNN)",
        yaxis_title="Immune Therapeutic Index (TPI)",
        height=450,
    )
    st.plotly_chart(fig_tpi, use_container_width=True)

    # =====================================================================
    # IMMUNE EVASION KINETICS
    # =====================================================================
    st.subheader("Immune Evasion Kinetics (Top-3 drugs)")
    fig_immune = go.Figure()
    best_3 = ranked_drugs[:3]

    for drug in best_3:
        t_kin = np.linspace(0, 14, 200)
        evas, score = simulate_immune_evasion_kinetics(
            t_kin, drug['immune_therapeutic_index'],
            drug['immune_toxicity_percent'],
        )
        fig_immune.add_trace(go.Scatter(
            x=t_kin, y=evas,
            name=drug['name'],
            line=dict(color="#60A5FA", width=3),
            showlegend=True,
        ))

    fig_immune.update_layout(
        template="plotly_dark",
        xaxis_title="Time (days post-transplant)",
        yaxis_title="Fraction of evading cells",
        height=350,
    )
    st.plotly_chart(fig_immune, use_container_width=True)

    # =====================================================================
    # DRUG COMBINATION RECOMMENDATIONS
    # =====================================================================
    st.subheader("Drug Combination Recommendations")
    st.markdown(
        "Based on GNN screening and immune analysis, optimal drug combinations "
        "are formed to minimize immune toxicity and maximize therapeutic index."
    )

    col_comb1, col_comb2 = st.columns(2)

    with col_comb1:
        st.markdown("### Golden Standard Combination")
        combo1 = ranked_drugs[:3]
        for drug in combo1:
            st.markdown(
                f"- **{drug['name']}** - TPI: {drug['therapeutic_index']:.3f}, "
                f"Class: {drug['family']}"
            )
        st.info(
            f"Expected TPI: **{combo1[0]['therapeutic_index']:.3f}** - "
            f"Optimal balance of therapy and immunity."
        )

    with col_comb2:
        st.markdown("### Alternative Combination")
        combo2 = ranked_drugs[3:6]
        for drug in combo2:
            st.markdown(
                f"- **{drug['name']}** - TPI: {drug['therapeutic_index']:.3f}, "
                f"Class: {drug['family']}"
            )
        st.info(
            f"Expected TPI: **{combo2[0]['therapeutic_index']:.3f}** - "
            f"Alternative option."
        )
