import os
# Set the cache path before importing celloracle
os.environ["CELLORACLE_DATA_DIR"] = os.path.join(os.getcwd(), "celloracle_cache")

import scanpy as sc
import celloracle as co
import pandas as pd
import logging
import time

logger = logging.getLogger(__name__)

def load_cached_base_grn(species_key: str, cache_dir="celloracle_cache/promoter_base_GRN"):
    """
    Load pre-cached .parquet GRN file for the selected species, with detailed timing and file info.
    """
    species_to_filename = {
        "human": "hg19_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet",
        "mouse": "mm10_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet",
        "pig": "Sscrofa11.1_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20230218.parquet",
        "rat": "rn6_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet",
        "zebrafish": "danRer11_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
        "chicken": "galGal6_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
        "celegans": "ce10_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
        "drosophila": "dm6_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
        "scerevisiae": "sacCer3_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
        "arabidopsis": "TAIR10_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
        "xenopus_laevis": "Xenopus_laevis_v10.1_TFinfo_dataframe_CisBPv2_Xenopus_laevis_fpr2_threshold_10_20221228.parquet",
        "xenopus_tropicalis": "xenTro3_TFinfo_dataframe_CisBPv2_Xenopus_tropicalis_fpr2_threshold_10_20221231.parquet"
    }

    if species_key not in species_to_filename:
        raise ValueError(f"❌ Unknown species: {species_key}")

    filepath = os.path.join(cache_dir, species_to_filename[species_key])
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"❌ GRN file for {species_key} not found at: {filepath}")

    file_size_mb = os.path.getsize(filepath) / (1024 * 1024)
    logger.info(f"📂 Loading GRN from: {filepath}")
    logger.info(f"📦 File size: {file_size_mb:.2f} MB")

    start_time = time.time()
    tf_info_df = pd.read_parquet(filepath)
    elapsed = time.time() - start_time

    logger.info(f"✅ Loaded GRN for '{species_key}' in {elapsed:.2f} seconds")
    return tf_info_df

# Mapping species to CellOracle GRN loader
SPECIES_TO_GRN_LOADER = {
    "human": co.data.load_human_promoter_base_GRN,
    "mouse": co.data.load_mouse_promoter_base_GRN,
    "pig": co.data.load_Pig_promoter_base_GRN,
    "rat": co.data.load_rat_promoter_base_GRN,
    "zebrafish": co.data.load_zebrafish_promoter_base_GRN,
    "chicken": co.data.load_chicken_promoter_base_GRN,
    "celegans": co.data.load_Celegans_promoter_base_GRN,
    "drosophila": co.data.load_drosophila_promoter_base_GRN,
    "scerevisiae": co.data.load_Scerevisiae_promoter_base_GRN,
    "arabidopsis": co.data.load_arabidopsis_promoter_base_GRN,
    "xenopus_laevis": co.data.load_xenopus_laevis_promoter_base_GRN,
    "xenopus_tropicalis": co.data.load_xenopus_tropicalis_promoter_base_GRN,
}


def run_grn_analysis(h5ad_path, species, clustering_method, alpha, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"🚀 GRN analysis started for species={species}, clustering={clustering_method}, alpha={alpha}")
    
    adata = sc.read(h5ad_path)

    # --- 1. Identify the Correct Column ---
    if clustering_method not in adata.obs:
        # If user asked for 'annotation' but it's missing, try to find a fallback
        found = False
        # Common names for annotation columns
        for fallback in ["cell_type", "annotation", "louvain_labels", "final_annotation"]:
            if fallback in adata.obs:
                logger.info(f"⚠️ '{clustering_method}' not found. Switching to '{fallback}'")
                clustering_method = fallback
                found = True
                break
        
        if not found:
            raise ValueError(f"❌ Column '{clustering_method}' not found in AnnData. Did you run the annotation step?")

    cluster_key = clustering_method

    # --- 2. CRITICAL FIX: Remove Cells with Missing Labels (NaN) ---
    # This prevents the "Task errored: nan" crash
    initial_cells = adata.n_obs
    
    # Check for NaNs
    if adata.obs[cluster_key].isnull().any():
        n_nan = adata.obs[cluster_key].isnull().sum()
        logger.warning(f"⚠️ Found {n_nan} cells with 'NaN' labels in '{cluster_key}'. Removing them...")
        adata = adata[~adata.obs[cluster_key].isnull()].copy()
    
    # Check for empty strings if the column is string/object type
    if adata.obs[cluster_key].dtype == "object":
        empty_mask = adata.obs[cluster_key] == ""
        if empty_mask.any():
            n_empty = empty_mask.sum()
            logger.warning(f"⚠️ Found {n_empty} cells with empty labels. Removing them...")
            adata = adata[~empty_mask].copy()

    if adata.n_obs < initial_cells:
        logger.info(f"📉 Reduced cell count from {initial_cells} to {adata.n_obs} after cleaning labels.")

    if adata.n_obs == 0:
        raise ValueError("❌ All cells were removed because they had missing labels! Please check your annotation file.")

    # --- 3. Format as Category ---
    if adata.obs[cluster_key].dtype.name != 'category':
        logger.info(f"⚠️ Converting '{cluster_key}' to category")
        adata.obs[cluster_key] = adata.obs[cluster_key].astype('category')
    
    # ---------------------------------------------------------------

    # Load base GRN
    if species not in SPECIES_TO_GRN_LOADER:
        raise ValueError(f"❌ Unsupported species '{species}'")

    downsample_to=30000
    if adata.shape[0] > downsample_to:
        print(f"Downsampling to {downsample_to} cells.")
        sc.pp.subsample(adata, n_obs=downsample_to, random_state=123)

    logger.info(f"📦 Loading base GRN for {species}")
    base_grn = load_cached_base_grn(species, cache_dir="celloracle_cache/promoter_base_GRN")
    
    if "counts" not in adata.layers:
        raise ValueError("❌ 'counts' layer not found. Please ensure you saved raw counts.")

    # FORCE Raw Counts: Overwrite X with the raw counts layer to avoid log-transform warnings
    adata.X = adata.layers["counts"].copy()

    # Build Oracle object
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata,
        cluster_column_name=cluster_key,
        embedding_name="X_umap"
    )
    oracle.import_TF_data(TF_info_matrix=base_grn)

    # PCA & Imputation
    oracle.perform_PCA()
    n_comps = min(50, 20)
    k = int(adata.n_obs * 0.025)
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k * 8, b_maxl=k * 4, n_jobs=4)

    # Save oracle object
    oracle_path = os.path.join(output_dir, "oracle.celloracle.oracle")
    oracle.to_hdf5(oracle_path)

    # GRN Inference
    links = oracle.get_links(cluster_name_for_GRN_unit=cluster_key, alpha=alpha)

    # Save raw/unfiltered links
    raw_links_path = os.path.join(output_dir, "raw_links.celloracle.links")
    links.to_hdf5(file_path=raw_links_path)

    links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

    # Save filtered links
    filtered_links_path = os.path.join(output_dir, "filtered_links.celloracle.links")
    links.to_hdf5(file_path=filtered_links_path)
    
    summary = {
        "oracle_path": oracle_path,
        "raw_links_path": raw_links_path,
        "filtered_links_path": filtered_links_path,
        "clusters": list(map(str, links.links_dict.keys()))
    }

    logger.info("✅ GRN analysis complete.")
    return summary

# if __name__ == "__main__":
#     import argparse

#     parser = argparse.ArgumentParser(description="Run CellOracle GRN analysis.")
#     parser.add_argument("--h5ad", required=True, help="Path to input .h5ad file")
#     parser.add_argument("--species", required=True, help="Species key (e.g., human, mouse, pig)")
#     parser.add_argument("--cluster", required=True, help="Clustering method (column in .obs)")
#     parser.add_argument("--alpha", type=float, default=10.0, help="Alpha value for GRN inference")
#     parser.add_argument("--outdir", default="celloracle_output", help="Directory to save output")

#     args = parser.parse_args()

#     run_grn_analysis(
#         h5ad_path=args.h5ad,
#         species=args.species,
#         clustering_method=args.cluster,
#         alpha=args.alpha,
#         output_dir=args.outdir
#     )