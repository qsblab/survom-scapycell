#!/usr/bin/env python
# coding: utf-8

import logging
import warnings
from pathlib import Path
from typing import Dict, Any, List, Optional, Union

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
import scipy.io

# --- Setup basic logging ---
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()],
)

# --- Default parameters for the analysis workflow ---
DEFAULT_PARAMS = {
    'input_type': 'mtx',
    'data_path': '.',
    'output_dir': 'sc_results',
    'min_genes_per_cell': 200,
    'min_counts_per_cell': 500,
    'min_cells_per_gene': 3,
    'max_pct_mt': 10.0,
    'n_top_genes': 2000,
    'regress_vars': ['total_counts', 'pct_counts_mt'],
    'n_neighbors': 15,
    'n_pcs': 50,
    'umap_min_dist': 0.5,
    'clustering_resolution': 0.5,
    'random_state': 42,
    'show_plots': True,
    'save_plots': True,
    'save_adata_path': 'processed_adata.h5ad',
}


def find_existing_file(base_path: Path, potential_filenames: List[str]) -> Optional[Path]:
    """Checks a list of potential filenames in a directory and returns the first one found."""
    for fname in potential_filenames:
        fpath = base_path / fname
        if fpath.exists():
            return fpath
    return None


def run_single_cell_analysis(params: Dict[str, Any]) -> Optional[sc.AnnData]:
    """
    Runs a standard single-cell RNA-seq analysis workflow based on Scanpy best practices.

    This workflow includes data loading (MTX or 10x H5), quality control, filtering,
    normalization, feature selection, dimensionality reduction (PCA, UMAP), and clustering.
    Plots are generated using Plotly and can be optionally saved as HTML files.

    Args:
        params (dict): A dictionary containing analysis parameters. Missing parameters
                       will be replaced by values from `DEFAULT_PARAMS`.

            - input_type (str): 'mtx' or '10x_h5'.
            - data_path (str): Path to the directory containing input files.
            - output_dir (str): Directory to save results (plots, AnnData object).
            - min_genes_per_cell (int): Min genes per cell for filtering.
            - min_counts_per_cell (int): Min counts (UMIs) per cell for filtering.
            - min_cells_per_gene (int): Min cells a gene must be in for filtering.
            - max_pct_mt (float): Max mitochondrial percentage per cell.
            - n_top_genes (int): Number of highly variable genes to select.
            - regress_vars (Optional[List[str]]): Variables to regress out.
            - n_neighbors (int): Number of neighbors for KNN graph.
            - n_pcs (int): Number of principal components to compute.
            - umap_min_dist (float): UMAP `min_dist` parameter.
            - clustering_resolution (float): Leiden clustering resolution.
            - random_state (int): Seed for reproducible PCA and UMAP.
            - show_plots (bool): Whether to display interactive plots.
            - save_plots (bool): Whether to save plots as HTML files.
            - save_adata_path (Optional[str]): Filename for the final AnnData object.

    Returns:
        Optional[sc.AnnData]: The processed AnnData object, or None if a critical error occurs.
    """
    # --- Parameter Merging and Validation ---
    config = DEFAULT_PARAMS.copy()
    config.update(params)

    input_type = config['input_type']
    data_path = Path(config['data_path'])
    output_dir = Path(config['output_dir'])
    min_genes_per_cell = config['min_genes_per_cell']
    min_counts_per_cell = config['min_counts_per_cell']
    min_cells_per_gene = config['min_cells_per_gene']
    max_pct_mt = config['max_pct_mt']
    n_top_genes = config['n_top_genes']
    regress_vars = config['regress_vars']
    n_neighbors = config['n_neighbors']
    n_pcs = config['n_pcs']
    umap_min_dist = config['umap_min_dist']
    clustering_resolution = config['clustering_resolution']
    random_state = config['random_state']
    show_plots = config['show_plots']
    save_plots = config['save_plots']
    save_adata_path = config['save_adata_path']

    logging.info("--- Starting Single-Cell Analysis ---")
    logging.info(f"Using parameters: {config}")

    # --- Setup Output Directory and Reproducibility ---
    output_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.set_figure_params(dpi=80, facecolor='white')
    sc.settings.seed = random_state

    if not data_path.exists():
        logging.error(f"Input data path does not exist: {data_path}")
        raise FileNotFoundError(f"Input data path does not exist: {data_path}")

    # --- 1. Read Data ---
    if input_type == 'mtx':
        logging.info(f"Reading MTX files from {data_path}")
        matrix_file = find_existing_file(data_path, ["matrix.mtx.gz", "matrix.mtx"])
        feature_file = find_existing_file(data_path, ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"])
        barcode_file = find_existing_file(data_path, ["barcodes.tsv.gz", "barcodes.tsv"])

        if not matrix_file: raise FileNotFoundError(f"Matrix file not found in {data_path}")
        if not feature_file: raise FileNotFoundError(f"Features/genes file not found in {data_path}")
        if not barcode_file: raise FileNotFoundError(f"Barcodes file not found in {data_path}")

        adata = sc.read_mtx(matrix_file).T
        genes = pd.read_csv(feature_file, header=None, sep='\t')
        barcodes = pd.read_csv(barcode_file, header=None, sep='\t')

        adata.var_names = genes.iloc[:, 0].astype(str)
        adata.obs_names = barcodes.iloc[:, 0].astype(str)
        
        if genes.shape[1] > 1:
            adata.var['gene_symbols'] = genes.iloc[:, 1].values
            logging.info("Found gene symbols in the second column of features file.")
        else:
            logging.warning("Only one column in features file. Using gene IDs as gene_symbols.")
            adata.var['gene_symbols'] = adata.var_names.values

    elif input_type == '10x_h5':
        h5_file = find_existing_file(data_path, ["filtered_feature_bc_matrix.h5", "filtered_gene_bc_matrix.h5"])
        if not h5_file: raise FileNotFoundError(f"No 10x HDF5 file found in {data_path}")
        logging.info(f"Reading 10x HDF5 file: {h5_file.name}")
        adata = sc.read_10x_h5(h5_file)
        adata.var['gene_symbols'] = adata.var_names.values # Use original names as symbols

    else:
        raise ValueError(f"Unsupported input_type: '{input_type}'. Use 'mtx' or '10x_h5'.")

    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    logging.info(f"Raw data dimensions: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # --- 2. Basic Filtering (Genes) ---
    logging.info(f"Filtering genes present in < {min_cells_per_gene} cells.")
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    logging.info(f"Dimensions after gene filtering: {adata.shape}")

    # --- 3. Quality Control (QC) ---
    logging.info("Calculating QC metrics...")
    # Use gene symbols for MT check, as var_names might be Ensembl IDs
    adata.var['mt'] = adata.var['gene_symbols'].str.upper().str.startswith(('MT-', 'MT.'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    if show_plots or save_plots:
        logging.info("Generating QC plots...")
        qc_df = adata.obs.copy()
        plot_configs = {
            'n_genes_by_counts': "QC: Genes per Cell (before filtering)",
            'total_counts': "QC: Total Counts per Cell (before filtering)",
            'pct_counts_mt': "QC: Mitochondrial % per Cell (before filtering)"
        }
        for y_val, title in plot_configs.items():
            fig = px.violin(qc_df, y=y_val, box=True, points="all", title=title)
            if show_plots: fig.show()
            if save_plots: fig.write_html(output_dir / f"{y_val}_violin.html")

        fig_scatter_qc = px.scatter(
            qc_df, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt',
            title="QC: Counts vs Genes (color = MT%)",
            labels={'total_counts': 'Total Counts', 'n_genes_by_counts': 'N Genes', 'pct_counts_mt': 'Mitochondrial %'}
        )
        if show_plots: fig_scatter_qc.show()
        if save_plots: fig_scatter_qc.write_html(output_dir / "qc_scatter.html")

    # --- 4. Filtering Cells based on QC ---
    logging.info("Applying cell filters based on QC metrics...")
    cells_before = adata.shape[0]
    sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)
    sc.pp.filter_cells(adata, min_counts=min_counts_per_cell)
    adata = adata[adata.obs.pct_counts_mt < max_pct_mt, :]
    logging.info(f"Filtered out {cells_before - adata.shape[0]} cells.")
    logging.info(f"Dimensions after cell filtering: {adata.shape}")

    # --- 5. Normalization & Transformation ---
    logging.info("Normalizing counts and log-transforming data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # --- 6. Feature Selection (Highly Variable Genes - HVGs) ---
    logging.info(f"Identifying top {n_top_genes} highly variable genes...")
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=False, flavor='seurat_v3')
    except ValueError as e:
        logging.warning(f"Could not compute highly variable genes. Low cell/gene count? Error: {e}")
        logging.warning("Skipping HVG selection and downstream steps.")
        return adata

    if show_plots or save_plots:
        if 'highly_variable' in adata.var:
            logging.info("Generating HVG plot...")
            hvg_df = adata.var.copy()
            fig_hvg = px.scatter(
                hvg_df.dropna(subset=['means', 'variances_norm']),
                x='means', y='variances_norm', color='highly_variable',
                hover_data=['gene_symbols'], title="Highly Variable Genes",
                labels={'means': 'Mean Expression', 'variances_norm': 'Normalized Variance'}
            )
            if show_plots: fig_hvg.show()
            if save_plots: fig_hvg.write_html(output_dir / "hvg_scatter.html")
        else:
            logging.warning("Could not generate HVG plot. 'highly_variable' not in adata.var.")

    logging.info("Subsetting data to HVGs...")
    adata.raw = adata
    if not adata.var['highly_variable'].any():
        logging.warning("No highly variable genes were identified. Cannot proceed with dimensionality reduction.")
        return adata
    adata = adata[:, adata.var.highly_variable]
    logging.info(f"Dimensions after HVG selection: {adata.shape}")

    # --- 7. Regress out unwanted variation & Scale data ---
    if regress_vars and len(regress_vars) > 0:
        logging.info(f"Regressing out variables: {regress_vars}...")
        sc.pp.regress_out(adata, keys=regress_vars)
    logging.info("Scaling data to unit variance...")
    sc.pp.scale(adata, max_value=10)

    # --- 8. Principal Component Analysis (PCA) ---
    logging.info("Running PCA...")
    n_comps = min(n_pcs, adata.shape[0] - 1, adata.shape[1] - 1)
    if n_comps < n_pcs:
        logging.warning(f"Requested {n_pcs} PCs, but reduced to {n_comps} due to data dimensions.")
    if n_comps < 2:
        logging.error(f"Cannot run PCA with fewer than 2 components ({n_comps}). Check filtering steps.")
        return adata
    sc.pp.pca(adata, n_comps=n_comps, random_state=random_state)

    if (show_plots or save_plots) and 'pca' in adata.uns:
        logging.info("Generating PCA variance plot...")
        fig_pca_var = px.line(
            y=np.cumsum(adata.uns['pca']['variance_ratio']),
            x=range(1, len(adata.uns['pca']['variance_ratio']) + 1),
            title='PCA Cumulative Variance Explained', markers=True, range_y=[0, 1.05],
            labels={'x': 'Principal Component', 'y': 'Cumulative Variance Ratio'}
        )
        if show_plots: fig_pca_var.show()
        if save_plots: fig_pca_var.write_html(output_dir / "pca_variance.html")

    # --- 9. Neighborhood Graph, UMAP, and Clustering ---
    effective_pcs = adata.obsm['X_pca'].shape[1]
    logging.info(f"Computing neighborhood graph (k={n_neighbors}, PCs={effective_pcs})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=random_state) # Use all computed PCs

    logging.info(f"Running UMAP (min_dist={umap_min_dist})...")
    sc.tl.umap(adata, min_dist=umap_min_dist, random_state=random_state)

    logging.info(f"Running Leiden clustering (resolution={clustering_resolution})...")
    sc.tl.leiden(adata, resolution=clustering_resolution, random_state=random_state, key_added='leiden')

    # --- 10. Visualize Clusters on UMAP ---
    if (show_plots or save_plots) and 'X_umap' in adata.obsm:
        logging.info("Generating UMAP plot colored by Leiden clusters...")
        umap_df = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names, columns=['UMAP1', 'UMAP2'])
        umap_df['Cluster'] = adata.obs['leiden'].astype('category')
        
        # Sort clusters numerically for consistent plotting
        cluster_order = sorted(umap_df['Cluster'].cat.categories, key=lambda x: int(x))

        fig_umap_clusters = px.scatter(
            umap_df, x='UMAP1', y='UMAP2', color='Cluster',
            title=f"UMAP colored by Leiden Cluster (res={clustering_resolution})",
            labels={'Cluster': 'Leiden Cluster'}, category_orders={"Cluster": cluster_order}
        )
        fig_umap_clusters.update_traces(marker=dict(size=3, opacity=0.8))
        fig_umap_clusters.update_layout(legend_title_text='Cluster')

        if show_plots: fig_umap_clusters.show()
        if save_plots: fig_umap_clusters.write_html(output_dir / "umap_leiden_clusters.html")

    # --- 11. Save Final AnnData Object ---
    if save_adata_path:
        final_path = output_dir / save_adata_path
        try:
            logging.info(f"Saving final AnnData object to {final_path}")
            if 'leiden' in adata.obs:
                sc.pl._utils.add_colors_for_categorical_sample_annotation(adata, 'leiden')
            adata.write_h5ad(final_path, compression="gzip")
        except Exception as e:
            logging.error(f"Could not save AnnData object to {final_path}. Error: {e}")

    logging.info("--- Analysis Complete ---")
    return adata


# --- Example Usage ---
if __name__ == "__main__":
    # --- IMPORTANT: Change data_path to your actual data location ---
    # This example assumes the public 1k PBMC dataset from 10x Genomics
    # Download from: https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5
    
    # Example for HDF5 file
    params_h5 = {
        'input_type': '10x_h5',
        'data_path': './pbmc1k_data/', # <<< --- CHANGE THIS PATH --- <<<
        'output_dir': 'results_pbmc1k_h5',
        'max_pct_mt': 5.0, # PBMCs have low MT content
        'n_pcs': 30,
        'clustering_resolution': 0.4,
    }

    # Example for MTX files (using the path from your original script)
    params_mtx = {
        'input_type': 'mtx',
        'data_path': '/mnt/c/Users/TAMALIKA/Onedrive/Desktop/survom/filtered_matrices_mex/hg19/', # <<< --- CHANGE THIS PATH --- <<<
        'output_dir': 'results_survom_mtx',
        'show_plots': True,
        'save_plots': True,
    }

    # Choose which parameters to run
    # To run, uncomment the desired line and ensure 'data_path' is correct.
    # chosen_params = params_h5 
    chosen_params = params_mtx

    try:
        adata_processed = run_single_cell_analysis(chosen_params)
        if adata_processed:
            print("\n--- SCRIPT FINISHED ---")
            print("Processed AnnData object info:")
            print(adata_processed)
            if 'leiden' in adata_processed.obs:
                print("\nCluster counts:")
                print(adata_processed.obs['leiden'].value_counts().sort_index())
                
    except FileNotFoundError as e:
        logging.error(f"Fatal Error: Input data not found. {e}")
        logging.error("Please ensure the 'data_path' parameter is correct and the files exist.")
    except Exception as e:
        logging.error(f"An unexpected error occurred during the analysis: {e}", exc_info=True)