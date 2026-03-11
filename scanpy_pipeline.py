import scanpy as sc
import shutil
import uuid
import time
import threading
import os, base64, uuid
from datetime import datetime, timedelta
import logging
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np 
UPLOAD_ROOT = "uploaded_data"  # Adjust to match your upload folder
px_size=1000

# Basic logger setup
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),  # logs to console
        logging.FileHandler("scanpy.log")  # logs to a file
    ]
)

logger = logging.getLogger(__name__)





def cleanup_folders(base_folder, interval_minutes=30, expiry_minutes=60):
    """Generic cleaner for folders under a given base directory."""
    while True:
        os.makedirs(base_folder, exist_ok=True)

        now = datetime.now()
        cutoff = now - timedelta(minutes=expiry_minutes)

        for folder in os.listdir(base_folder):
            folder_path = os.path.join(base_folder, folder)
            try:
                if os.path.isdir(folder_path):
                    folder_time = datetime.fromtimestamp(os.path.getctime(folder_path))
                    if folder_time < cutoff:
                        shutil.rmtree(folder_path)
                        print(f"[CLEANUP] Deleted: {folder_path}")
            except Exception as e:
                print(f"[CLEANUP ERROR] Failed to delete {folder_path}: {e}")

        time.sleep(interval_minutes * 60)


def plot_pca_variance(adata, n_pcs):

    variance_ratio = adata.uns['pca']['variance_ratio'][:n_pcs]
    pcs = [f'PC{i+1}' for i in range(len(variance_ratio))]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=pcs, y=variance_ratio, mode='lines+markers'))

    fig.update_layout(
        title='PCA Variance Ratio',
        xaxis_title='Principal Component',
        yaxis_title='Variance Ratio',
        width=px_size, height=px_size
    )
    return fig





def plot_umap(adata, color_by="louvain"):
    """
    Plot UMAP colored either by cluster (categorical) or gene expression (continuous).

    Parameters:
    - adata: AnnData object
    - color_by: str, cluster key or gene name

    Returns:
    - Plotly figure
    """
    import plotly.graph_objs as go
    import plotly.express as px
    import pandas as pd
    import numpy as np

    if "X_umap" not in adata.obsm:
        raise ValueError("UMAP coordinates not found in adata.obsm")

    umap_1 = adata.obsm['X_umap'][:, 0]
    umap_2 = adata.obsm['X_umap'][:, 1]

    # Determine whether color_by is cluster (obs) or gene (var)
    if color_by in adata.obs.columns:
        # Categorical coloring
        cluster_data = adata.obs[color_by].astype(str)
        unique_categories = sorted(pd.unique(cluster_data), key=lambda x: int(x) if x.isdigit() else x)
        color_map = px.colors.qualitative.Safe
        color_dict = {category: color_map[i % len(color_map)] for i, category in enumerate(unique_categories)}

        fig = go.Figure()
        for category in unique_categories:
            indices = cluster_data == category
            fig.add_trace(go.Scatter(
                x=umap_1[indices],
                y=umap_2[indices],
                mode='markers',
                marker=dict(
                    size=6,
                    color=color_dict[category],
                    line=dict(width=0.3, color='black')
                ),
                name=f"Cluster {category}"
            ))

        fig.update_layout(title=f"UMAP colored by cluster: '{color_by}'",width=px_size, height=px_size)

    elif color_by in adata.var_names:
        # Continuous gene expression coloring
        expr = adata[:, color_by].X.toarray().flatten() if hasattr(adata[:, color_by].X, "toarray") else adata[:, color_by].X

        fig = go.Figure(data=go.Scatter(
            x=umap_1,
            y=umap_2,
            mode='markers',
            marker=dict(
                size=6,
                color=expr,
                colorscale='Viridis',
                colorbar=dict(title=color_by),
                showscale=True
            ),
            text=[f"{color_by}: {val:.2f}" for val in expr]
        ))

        fig.update_layout(title=f"UMAP colored by expression of '{color_by}'")

    else:
        raise ValueError(f"'{color_by}' not found in adata.obs or adata.var_names")

    fig.update_layout(
        xaxis_title='UMAP1',
        yaxis_title='UMAP2',
        paper_bgcolor='white',
        plot_bgcolor='white',
        width=px_size, height=px_size
    )


    # 2) Enforce a 1:1 aspect ratio (one unit in x equals one unit in y)
    # fig.update_yaxes(scaleanchor="x", scaleratio=1)
    # fig.update_xaxes(constrain='domain')

    return fig



def plot_umap_clusters2(adata, cluster_key="louvain"):
    """
    Plot UMAP scatter plot with clusters using different colors.
    
    Parameters:
    - adata: AnnData object containing UMAP and clustering information
    - cluster_key: Key in adata.obs to use for coloring clusters (e.g., 'louvain', 'leiden')
    
    Returns:
    - fig: Plotly Figure object
    """
    # Ensure UMAP has been computed
    if "X_umap" not in adata.obsm:
        raise ValueError("UMAP coordinates not found in `adata.obsm['X_umap']`")

    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")

    # Extract UMAP coordinates
    umap_1 = adata.obsm['X_umap'][:, 0]
    umap_2 = adata.obsm['X_umap'][:, 1]

    # Get cluster data
    cluster_data = adata.obs[cluster_key].astype(str)

    # Define unique clusters and color map
    unique_categories = sorted(pd.unique(cluster_data), key=lambda x: int(x) if x.isdigit() else x)
    color_map = px.colors.qualitative.Safe
    color_dict = {category: color_map[i % len(color_map)] for i, category in enumerate(unique_categories)}

    # Initialize plot
    fig = go.Figure()

    # Add trace for each cluster
    for category in unique_categories:
        indices = cluster_data == category
        fig.add_trace(go.Scatter(
            x=umap_1[indices],
            y=umap_2[indices],
            mode='markers',
            marker=dict(
                size=6,
                color=color_dict[category],
                line=dict(width=0.5, color='black')
            ),
            name=f"Cluster {category}"
        ))

    # Update layout
    fig.update_layout(
        title=f"UMAP Clustering by '{cluster_key}'",
        xaxis_title="UMAP1",
        yaxis_title="UMAP2",
        height=600,
        paper_bgcolor='white',
        plot_bgcolor='white'
    )

    return fig





def save_uploaded_files(contents_list, filenames_list, debug=False):
    if debug:
        upload_dir = "uploaded_data/debug_session"
    else:
        session_id = str(uuid.uuid4())[:8]
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        upload_dir = f"uploaded_data/{timestamp}_{session_id}"
        
    os.makedirs(upload_dir, exist_ok=True)
    saved_paths = []

    for content, filename in zip(contents_list, filenames_list):
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)
        filepath = os.path.join(upload_dir, filename)
        with open(filepath, 'wb') as f:
            f.write(decoded)
        saved_paths.append(filepath)

    return saved_paths, upload_dir






def plot_volcano(df: pd.DataFrame, cluster_id: str, lfc_cutoff=1.0, pval_cutoff=0.05):
    """
    Create a volcano plot from DE results DataFrame.
    Required columns: ['gene', 'log2FC', 'pvals_adj']
    """
    df = df.copy()

    if "pvals_adj" not in df.columns:
        raise ValueError("Missing 'pvals_adj' column in DE results.")

    df = df.rename(columns={"pvals_adj": "padj"})  # Normalize column name

    df = df.dropna(subset=["log2FC", "padj", "gene"])
    df["-log10(padj)"] = -np.log10(df["padj"] + 1e-300)

    def classify(row):
        if row["padj"] < pval_cutoff:
            if row["log2FC"] > lfc_cutoff:
                return "Up"
            elif row["log2FC"] < -lfc_cutoff:
                return "Down"
        return "Not Sig"

    df["significance"] = df.apply(classify, axis=1)

    color_map = {
        "Up": "red",
        "Down": "blue",
        "Not Sig": "gray"
    }

    fig = px.scatter(
        df,
        x="log2FC",
        y="-log10(padj)",
        color="significance",
        color_discrete_map=color_map,
        hover_name="gene",
        title=f"Volcano Plot – Cluster {cluster_id}",
        labels={"-log10(padj)": "-log10(p-adjusted)", "log2FC": "log2 Fold Change"},
    )

    fig.update_traces(marker=dict(size=6, opacity=0.7, line=dict(width=0.3, color='black')))
    fig.update_layout(
        margin=dict(l=10, r=10, t=40, b=10),
        paper_bgcolor='white',
        plot_bgcolor='white',
        legend_title_text="Significance",width=px_size, height=px_size
    )

    return fig




def get_volcano_data(adata, cluster_id):
    """
    Extract volcano plot data for a given cluster from AnnData DE results.
    Returns a DataFrame with columns: ['gene', 'log2FC', 'pvals_adj']
    """
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    if cluster_id not in groups:
        raise ValueError(f"Cluster {cluster_id} not found in DE results.")

    names = result["names"][cluster_id]
    logfc = result["logfoldchanges"][cluster_id]
    pvals_adj = result["pvals_adj"][cluster_id]

    df = pd.DataFrame({
        "gene": names,
        "log2FC": logfc,
        "pvals_adj": pvals_adj
    })

    return df




def run_scanpy_analysis(
    dataset_path: str,
    min_genes_per_cell: int,
    min_cells_per_gene: int,
    norm_method: str,
    n_variable_genes: int,
    feature_selection_method: str,
    n_pcs: int,
    umap_min_dist: float,
    umap_n_neighbors: int,
    clustering_resolution: float,
    clustering_method: str,
    ranking_method: str,
    tool_choice: str = "scanpy",
    delete_after_hours: int = 10,
    debug: bool = False,
    return_figs=True,
     **kwargs
) -> str:
    """
    Full Scanpy analysis with UI-driven parameters.
    """

    mt_enabled = kwargs.get("mt_filter_enabled", False)
    mt_lower   = kwargs.get("mt_lower")
    mt_upper   = kwargs.get("mt_upper")
    mt_prefix  = (kwargs.get("mt_prefix") or "MT-").strip()

    if not os.path.exists(dataset_path):
        logger.error(f"Input dataset folder not found: {dataset_path}")
        raise FileNotFoundError("Input dataset folder not found.")

    # Create unique output directory
    base_dir = "output"
    os.makedirs(base_dir, exist_ok=True)

    if debug:
        output_dir = os.path.join(base_dir, "debug_output")
    else:
        output_dir = os.path.join(base_dir, f"results_{uuid.uuid4().hex}")

    os.makedirs(output_dir, exist_ok=True)

    try:
        # Load data
        logger.info("Loading 10X matrix from: %s", dataset_path)
        try:
            is_compressed = all(
                os.path.exists(os.path.join(dataset_path, f))
                for f in ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]
            )
            adata = sc.read_10x_mtx(
                dataset_path,
             # helps with speed
            )
            logger.info(f"✅ Loaded {adata.shape[0]} cells and {adata.shape[1]} genes.")
            
        except Exception as e:
            logger.error(f"❌ Failed to load 10X data: {e}")
            raise
        adata.layers["counts"] = adata.X.copy()
        adata.raw = adata
        # Quality Control
        logger.info("Filtering cells and genes")
        sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)
        sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
        # Case-insensitive startswith using the provided prefix
        names = (adata.var.get("gene_symbols") or adata.var_names).astype(str)
        prefix = mt_prefix.upper()
        mt_mask = pd.Index(names).str.upper().str.startswith(prefix)
        adata.var["mt"] = np.asarray(mt_mask)
        
        logger.info(f"before mitochondria: {adata.shape}")

        # Compute QC metrics (adds 'pct_counts_mt')
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

        # Apply filter only if enabled and bounds provided
        if mt_enabled and mt_lower is not None and mt_upper is not None:
            mask = (adata.obs["pct_counts_mt"] >= float(mt_lower)) & (adata.obs["pct_counts_mt"] <= float(mt_upper))
            adata = adata[mask].copy()
            logger.info(f"after mitochondria: {adata.shape}")

        # Normalization
        logger.info(f"Running normalization method: {norm_method}")
        if norm_method.lower() == "median":
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
        elif norm_method.lower() == "cpm":
            sc.pp.normalize_total(adata,target_sum=1e6)
            sc.pp.log1p(adata)
        # Feature selection
        flavor_map = {
            "vst": "seurat_v3",
            "mean.var.plot": "cell_ranger",
            "dispersion": "cell_ranger"
        }
        flavor = flavor_map.get(feature_selection_method.lower(), "seurat_v3")
        logger.info(f"Selecting {n_variable_genes} variable genes using {flavor}")
        sc.pp.highly_variable_genes(adata, n_top_genes=n_variable_genes, flavor=flavor)
        adata = adata[:, adata.var.highly_variable]

        # PCA
        logger.info(f"Running PCA with {n_pcs} components")
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)

        # Neighbors & UMAP
        logger.info("Computing neighbors and UMAP")
        sc.pp.neighbors(adata, n_neighbors=umap_n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(adata, min_dist=umap_min_dist)

        # Clustering
        logger.info(f"Running clustering using {clustering_method}")
        if clustering_method == "louvain":
            sc.tl.louvain(adata, resolution=clustering_resolution)
            cluster_key = "louvain"
        else:
            sc.tl.leiden(adata, resolution=clustering_resolution)
            cluster_key = "leiden"

        # Marker gene identification
        adata.layers["log_counts"] = sc.pp.log1p(adata.layers["counts"].copy())

        logger.info(f"Ranking genes by: {ranking_method}")
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=ranking_method,use_raw=False,
        layer="log_counts")

        # Save DE results to CSV
        if "rank_genes_groups" in adata.uns:
            try:
                result = adata.uns["rank_genes_groups"]
                groups = result["names"].dtype.names

                for group in groups:
                    df = pd.DataFrame({
                        "gene": result["names"][group],
                        "log2FC": result["logfoldchanges"][group],
                        "pval": result["pvals"][group],
                        "pval_adj": result["pvals_adj"][group]
                    })
                    df.to_csv(os.path.join(output_dir, f"de_cluster_{group}.csv"), index=False)
            except Exception as e:
                logger.warning(f"⚠️ Failed to export DE CSVs: {e}")



        # Save results
        adata_path = os.path.join(output_dir, "adata_processed.h5ad")
        adata.write(adata_path)
        

    
        # Generate figures (using modular functions)
        pca_fig = plot_pca_variance(adata, n_pcs)
        umap_fig = plot_umap(adata, cluster_key)

        logger.info(f"Saved AnnData object: {adata_path}")


        ### plot umap


        
        return (output_dir, pca_fig, umap_fig,adata_path) if return_figs else output_dir

    except Exception as e:
        logger.exception("❌ Error during Scanpy analysis:")
        raise e
