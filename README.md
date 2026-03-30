# Survom-scapycell

Scapycell is an interactive web platform for performing end-to-end single-cell RNA-seq analysis and Gene Regulatory Network (GRN) Inference. It is designed to identify master transcription factors and visualize complex cell populations. 

### Requirements

This module was developed and optimized for **Python 3.10**, and it is strictly recommended to run Python 3.10 to ensure compatibility with complex bioinformatics dependencies such as `velocyto` and `gimmemotifs`. 

Scapycell depends on several core Python packages including:
* Scanpy
* Celloracle
* Dash & Plotly
* Cytoscape.js
* Velocyto
* Gimmemotifs

The installation process should handle these dependencies for you. However, because `velocyto` requires pre-compiled binaries, we highly recommend using `conda` alongside `pip` for the setup.

### Setup

You will need to install Conda/Miniconda and create an environment using Python 3.10:


```conda create -n survom python=3.10```
```conda activate survom```

Installing from Source (Recommended)
Inside your activated survom environment, download the source code from GitHub and install the specific dependencies:

Bash
# Clone the repository
```git clone [https://github.com/qsblab/survom-scapycell.git](https://github.com/qsblab/survom-scapycell.git)```
```cd survom-scapycell```

# Install heavy dependencies via Conda to avoid build errors
```conda install -c conda-forge -c bioconda numpy=1.26.4 velocyto.py -y```

# Ensure compatible setuptools for legacy packages
```pip install "setuptools<70"```

# Install remaining requirements
```pip install -r requirements.txt```

With your environment activated and dependencies installed, you can immediately launch the interactive application:

`python app.py`

The terminal will generate a local web address (typically http://127.0.0.1:8050/scom/). Open this link in your web browser to access the graphical interface.

# Survom Documentation & Tutorials
We have provided comprehensive documentation to guide users step-by-step through running Survom for different case studies.
* Main Tutorial: Explains how to upload standard 10x Genomics matrices and seamlessly perform QC filtering, normalization, and highly variable gene (HVG) selection.

* GRN Inference: Explains how to build and visualize complex regulatory networks.

# Usage
First, launch the app via app.py. Ensure you have your single-cell matrices ready (e.g., .mtx, .tsv, or .h5ad formats).

Upload your data via the interface. If you are computing differential expression to identify cluster-specific marker genes, make sure to configure the thresholds in the QC tab before running the Scanpy pipeline. Finally, navigate to the Network Visualization tab to explore interactive node-and-edge graphs.

# License
The software in this repository is licensed under the Apache-2.0 License. Please see the LICENSE file for more details.


