import plotly.io as pio
pio.renderers.default = "browser"  # 👈 REQUIRED
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
    Plot 3 scatter subplots (1 row, 3 cols) showing centrality scores for a TF across clusters.
    Color and X-position represent score. Y-axis is cluster.

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

    # Validate and prep
    df = merged_df.copy()
    if "tf" not in df.columns:
        df["tf"] = df.index
    df["cluster"] = df["cluster"].astype(str)

    tf_df = df[df["tf"] == tf_name]
    if tf_df.empty:
        raise ValueError(f"❌ TF '{tf_name}' not found in merged_df.")

    tf_df = tf_df.sort_values("cluster")

    fig = make_subplots(
        rows=1,
        cols=3,
        shared_yaxes=True,
        horizontal_spacing=0.1,
        subplot_titles=[m.replace("_", " ").title() for m in metrics]
    )

    for i, metric in enumerate(metrics):
        fig.add_trace(
            go.Scatter(
                x=tf_df[metric],
                y=tf_df["cluster"],
                mode="markers",
                marker=dict(
                    size=marker_size,
                    color=tf_df[metric],
                    colorscale="Viridis",
                    colorbar=dict(title="Score") if i == 2 else None,
                ),
                text=tf_df[metric].round(3),
                hovertemplate=f"Cluster: %{{y}}<br>{metric}: %{{x}}<extra></extra>",
                name=metric
            ),
            row=1,
            col=i+1
        )

        fig.update_xaxes(title_text="Score", row=1, col=i+1)

    fig.update_layout(
        title_text=f"Centrality Scores Across Clusters for TF: {tf_name}",
        height=height,
        width=1100,
        template="plotly_white",
        margin=dict(t=60, l=60, r=20, b=40),
        showlegend=False
    )

    fig.update_yaxes(title_text="Cluster", row=1, col=1)

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
            "tf": "Transcription Factor",
            centrality_metric: "Centrality Score"
        },
        title=f"Top {top_n} TFs in Cluster {cluster_id} by {centrality_metric.replace('_', ' ').title()}",
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


# === CONFIGURATION ===
LINKS_FILE = "/data/users/vikash/scom/tamalika_grn/links.celloracle.links"
CLUSTER_NAME = "9"       # Cluster to analyze
TF_NAME = "GATA2"        # Transcription factor to filter
TOP_N = 15               # Number of top TFs per cluster for centrality plot


# === LOAD LINKS ===
print("📥 Loading links object...")
links = co.load_hdf5(LINKS_FILE)


# === PART 1: FILTERED GRN LINKS TABLE ===
print(f"📊 Extracting filtered links for cluster {CLUSTER_NAME} and TF {TF_NAME}")
if CLUSTER_NAME not in links.filtered_links:
    raise ValueError(f"❌ Cluster '{CLUSTER_NAME}' not found in links.")

filtered_links_df = links.filtered_links[CLUSTER_NAME]
tf_links_df = filtered_links_df[filtered_links_df["source"] == TF_NAME]

print(f"✅ Found {len(tf_links_df)} links where source == '{TF_NAME}' in cluster {CLUSTER_NAME}\n")
print(tf_links_df.head())

# === OPTIONAL: Save to CSV ===
tf_links_df.to_csv(f"links_from_{TF_NAME}_cluster_{CLUSTER_NAME}.csv", index=False)
print(f"💾 Saved filtered TF links to: links_from_{TF_NAME}_cluster_{CLUSTER_NAME}.csv")


# === PART 2: CENTRALITY PLOT FOR TOP TFs ACROSS CLUSTERS ===
print("\n📈 Generating centrality comparison plot...")

# Ensure scores are computed
links.get_network_score()
merged = links.merged_score.copy()

# Add TF name as a column (from index)
merged["tf"] = merged.index

fig = plot_tf_centrality_across_clusters(merged_df=merged, tf_name="GATA2")
fig.show()


fig = plot_cluster_centrality_from_df(
    merged_df=merged,
    cluster_id="9",
    centrality_metric="eigenvector_centrality",  # or degree_centrality_all
    top_n=20
)

fig.show()



# Get top TFs by eigenvector_centrality per cluster
top_genes = (
    merged.groupby("cluster")
    .apply(lambda df: df.nlargest(TOP_N, "eigenvector_centrality"))
    .reset_index(drop=True)
)

# Create long-form DataFrame for plotting
plot_df = top_genes[["cluster", "tf", "degree_all", "betweenness_centrality", "eigenvector_centrality"]].melt(
    id_vars=["cluster", "tf"],
    var_name="centrality_metric",
    value_name="score"
)

# Plot eigenvector centrality only for now
ev_df = plot_df[plot_df["centrality_metric"] == "eigenvector_centrality"]

fig = px.bar(
    ev_df,
    x="score",
    y="tf",
    color="cluster",
    orientation="h",
    title=f"Top {TOP_N} TFs by Eigenvector Centrality per Cluster",
    labels={"score": "Eigenvector Centrality", "tf": "Transcription Factor"}
)

fig.update_layout(height=700, yaxis={'categoryorder':'total ascending'})
fig.show()
