# Survom: End-to-End Single-Cell & GRN Tutorial

Welcome to Survom! This tutorial is designed for first-time users and covers all topics from loading datasets, performing single-cell RNA-seq preprocessing, to the visualization and biological interpretation of Gene Regulatory Networks (GRNs). 

This tutorial assumes a basic knowledge of:
* The scRNA-seq protocol
* Transcriptomic data analysis
* Basics of gene expression and its relevance to cell identity

---

## Phase 1: Single-Cell Preprocessing & Clustering

This phase cleans your raw data, corrects for technical noise, and groups similar cells together so you can identify distinct cell populations.

### 1. Upload Your Dataset
Load your expression data into the tool. Survom supports standard 10x Genomics output files (`matrix.mtx`, `barcodes.tsv`, and `features.tsv`). In these matrices, genes should be in rows and samples (cells) should be in columns. Ensure your files are uncompressed or standard gzipped formats before uploading.

### 2. Quality Control (QC) & Filtering
Apply threshold filters to exclude cells with too few or too many genes detected (indicating empty droplets or multiplet/doublet artifacts).
* **Min. Genes per Cell:** Set to 200 as a starting point to remove empty droplets. 
* **Max Genes (Implicit):** While our tool focuses on minimums, highly extreme gene counts naturally drop out during normalization. 
* **Min. Cells per Gene:** Set to 3 to remove uninformative genes and speed up computation.
* **Mitochondrial Filter (% mt):** Dying cells leak standard RNA but retain mitochondrial RNA. Enable the filter, input the MT prefix (e.g., `MT-`), and set an upper limit of 5% to 10% to ensure you are only analyzing healthy, viable cells.

### 3. Normalization & Feature Selection
* **Normalization:** Raw counts must be corrected for differences in sequencing depth. Median Normalization or Counts Per Million (CPM) are excellent defaults. The data is automatically log-transformed to stabilize variance.
* **Feature Selection:** Select "VST" (Variance Stabilizing Transformation) and retain the top 2,000 Highly Variable Genes (HVGs). This isolates the genes that capture true biological variation rather than baseline "housekeeping" functions.

### 4. Dimensionality Reduction (PCA & UMAP)
* **PCA:** Run Principal Component Analysis to compress the dataset. Retaining the first 30 principal components is a balanced choice for most datasets.
* **UMAP:** Set the minimum distance (`min_dist`) to 0.3 for a balance between global and local structure, and use 15 to 30 neighbors. This will flatten your complex data into a readable 2D map. 


### 5. Clustering & Marker Identification
* **Clustering:** Use the Louvain algorithm. A resolution of 0.5 to 0.8 works perfectly for datasets of moderate complexity (5,000–20,000 cells). Increase this value if you want finer, more specific sub-clusters.
* **Differential Expression:** Choose the Wilcoxon rank-sum test to identify "Marker Genes" that are significantly enriched in each cluster. 
* *Downstream Step:* You can export these marker genes to perform Cell Type Annotation (matching signatures with public databases) and Functional Enrichment (e.g., KEGG, GO pathways) offline.

---

## Phase 2: Gene Regulatory Network (GRN) Inference

Once your cells are clustered, use this phase to discover the Transcription Factors (the "master switches") controlling those specific cell identities.

### 1. Initialization & Target Selection
* Select your Organism (e.g., Mouse or Human) so Survom can load the correct Base GRN containing known DNA binding sites.
* Choose to build the network based on your Computed Clusters (from Phase 1) or an Uploaded Annotation column, and specify your target lineage (e.g., "Erythroids").

### 2. Network Inference Parameters
* **Sparsity (Alpha):** Set Alpha to 10. Higher values apply a strict mathematical penalty, forcing the algorithm to delete weak connections and only retain the highest-confidence regulatory links.
* **Link Filtering:** Keep only the top 2,000 links to prevent the final network from looking like an unreadable hairball.

---

## Phase 3: Image Analysis & Interactive Results

Survom generates highly interactive plots and downloadable tables for publication and further research. 

### 1. Interactive Image Analysis
Every plot generated in Survom (UMAP, Volcano Plots, Centrality Scatters) can be zoomed, panned, and adjusted. Hover over points to reveal detailed statistical information. On the top right corner of all result images, you will find a toolbar with the following options:
* **Download:** Download the plot as a high-resolution PNG or SVG format for publication.
* **Zoom:** Click and drag to zoom into specific subgraphs or dense clusters.
* **Pan:** Move seamlessly across the image.
* **Select:** Use the Box or Lasso select tools to highlight specific areas of the plot.
* **Reset:** Autoscale or Reset the axes of the image to return to the default view.

### 2. Network Visualization
* **Centrality Plot:** Identifies the master regulators of your cluster. Genes with an Eigenvector score near 1.0 are the most critical drivers of that cell identity.
* **Cytoscape Viewer:** After extracting a subnetwork for specific target genes, you can interact with the Cytoscape graph. Drag nodes to rearrange them for clarity, and click "Download Network Image" to save the topology. 


### 3. Data Export
When you are satisfied with your analysis, use the "Export Results" section. You can download your processed Single-Cell data (`.h5ad`), Marker Gene lists (`.csv`), and Full Network Tables (`.csv`), or simply click **"Download All Results"** to get a comprehensive ZIP archive of your entire session.