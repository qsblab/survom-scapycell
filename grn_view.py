import pandas as pd
import celloracle as co
import plotly.express as px
import plotly.graph_objects as go
import plotly.express as px

from plotly.subplots import make_subplots


def plot_tf_centrality_across_clusters(
    merged_df: pd.DataFrame,
    tf_name: str = "GATA2",
    height: int = 500,
    marker_size: int = 14
):
    """
    Plot 3 compact scatter subplots showing centrality scores for a TF across clusters.
    Each subplot has its own small, responsive colorbar.

    Parameters:
    - merged_df: CellOracle merged_score DataFrame
    - tf_name: TF to plot (default: 'GATA2')
    - height: Total figure height
    - marker_size: Size of the marker points

    Returns:
    - fig: Plotly Figure with 3 subplots
    """

    metrics = [
        "degree_centrality_all",
        "betweenness_centrality",
        "eigenvector_centrality"
    ]

    # Prep and filter
    df = merged_df.copy()
    if "tf" not in df.columns:
        df["tf"] = df.index
    df["cluster"] = df["cluster"].astype(str)

    tf_df = df[df["tf"] == tf_name]
    if tf_df.empty:
        raise ValueError(f"❌ TF '{tf_name}' not found in merged_df.")

    tf_df = tf_df.sort_values("cluster")

    # Create subplots
    fig = make_subplots(
        rows=1,
        cols=3,
        shared_yaxes=True,
        horizontal_spacing=0.1,
        subplot_titles=[m.replace("_", " ").title() for m in metrics]
    )

    for i, metric in enumerate(metrics):
        if metric not in tf_df.columns:
            continue

        fig.add_trace(
            go.Scatter(
                x=tf_df[metric],
                y=tf_df["cluster"],
                mode="markers",
                marker=dict(
                    size=marker_size,
                    color=tf_df[metric],
                    colorscale="Viridis",
                    showscale=True,
                    colorbar=dict(
                        title="",  # Optional: remove title for compactness
                        thickness=8,
                        len=0.5,
                        x=0.31 * (i+1),  # tightly aligned with subplot
                        xanchor="left",
                        xpad=2
                    ),
                    line=dict(width=1, color="DarkSlateGrey")
                ),
                text=tf_df[metric].round(3),
                hovertemplate=f"Cluster: %{{y}}<br>{metric}: %{{x}}<extra></extra>",
                name=metric
            ),
            row=1,
            col=i+1
        )

        fig.update_xaxes(title_text="Score", row=1, col=i+1)

    # Layout
    fig.update_layout(
        title_text=f"Centrality Scores Across Clusters for TF: {tf_name}",
        margin=dict(l=60, r=20, t=50, b=50),
        template="plotly_white",
        hovermode="closest",
        transition_duration=300,
        showlegend=False,
        autosize=True
    )

    fig.update_yaxes(
        title_text="Cluster",
        row=1, col=1,
        automargin=True,
        autorange="reversed"
    )

    return fig



def plot_cluster_centrality_from_df(
    merged_df: pd.DataFrame,
    cluster_id: str = "9",
    centrality_metric: str = "degree_centrality_all",
    top_n: int = 20,
    height: int = 500
):
    """
    Create a responsive Plotly figure of top N TFs in a cluster based on a centrality metric.
    
    Parameters:
    - merged_df: DataFrame from CellOracle's `links.merged_score`
    - cluster_id: cluster to filter by (as str or int)
    - centrality_metric: metric column (e.g. 'eigenvector_centrality')
    - top_n: number of top TFs to plot
    - height: plot height in pixels
    
    Returns:
    - fig: Plotly Figure (responsive layout, ready for Dash or browser)
    """

    # Ensure correct types and columns
    df = merged_df.copy()
    if "tf" not in df.columns:
        df["tf"] = df.index

    df["cluster"] = df["cluster"].astype(str)
    cluster_id = str(cluster_id)

    if centrality_metric not in df.columns:
        raise ValueError(f"❌ Centrality metric '{centrality_metric}' not found in DataFrame.")

    # Filter and sort
    cluster_df = df[df["cluster"] == cluster_id]
    if cluster_df.empty:
        raise ValueError(f"❌ No TFs found for cluster '{cluster_id}'.")

    top_df = cluster_df.sort_values(by=centrality_metric, ascending=False).head(top_n)

    # Plot
    fig = px.scatter(
        top_df,
        x=centrality_metric,
        y="tf",
        text="tf",
        color=centrality_metric,
        size=centrality_metric,
        color_continuous_scale="Viridis",
        labels={
            "tf": "Genes",
            centrality_metric: "Centrality Score"
        },
        title=f"Top {top_n} Genes in Cluster {cluster_id} by {centrality_metric.replace('_', ' ').title()}",
        height=height
    )

    # Layout tweaks
    fig.update_layout(
        yaxis=dict(autorange="reversed", automargin=True),
        xaxis=dict(automargin=True),
        margin=dict(l=80, r=30, t=60, b=50),
        showlegend=False,
        hovermode="closest",
        transition_duration=300,
        template="plotly_white",
    )

    fig.update_traces(
        marker=dict(line=dict(width=1, color="DarkSlateGrey")),
        textposition="top right"
    )

    return fig
