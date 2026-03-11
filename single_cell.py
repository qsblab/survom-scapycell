import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated.*")

import dash
import dash_bootstrap_components as dbc
from dash import html, dcc, Input, Output, State
from scanpy_pipeline import run_scanpy_analysis, save_uploaded_files, cleanup_folders
import base64, os, io
from datetime import datetime
import threading
from flask_session import Session
import logging
import uuid
from flask import Flask, request, jsonify

def single_cell_layout():
    
    # --- 1. SETUP FORM CARD ---
    setup_card = dbc.Card(
        dbc.CardBody([
            # Top Header & Action Buttons
            dbc.Row([
                dbc.Col(html.H4("Analysis Setup", className="fw-bold mb-0"), md=6),
                dbc.Col([
                    dbc.Button(
                        [html.I(className="bi bi-download me-2"), "Load Example Data"],
                        id="load-example-singlecell-button",
                        color="outline-secondary",
                        size="sm",
                        className="me-2 shadow-sm"
                    ),
                    dbc.Button("🔁 Fresh Start", id="fresh-start-button", color="warning", size="sm", className="shadow-sm"),
                ], md=6, className="text-end")
            ], className="mb-2 align-items-center"),
            
            html.Div(id="fresh-start-status", className="text-muted text-end small mb-4"),

            # Step 1: Upload
            html.Div([
                dbc.Label("1. Upload Dataset", className="fw-semibold text-secondary"),
                html.Div(id="upload-status", className="mb-2"),
                html.Div([
                    html.Iframe(
                        srcDoc=open("dropzone.html").read(),  # dropzone.html is the custom uploader
                        style={"width": "100%", "height": "130px", "border": "none"}
                    ),
                ], className="border border-2 border-secondary border-dashed rounded bg-light p-2 mb-1"),
                html.Small("Format: 10x Matrix (matrix.mtx, features.tsv, barcodes.tsv — gzipped or not)", className="text-muted d-block mb-4"),
            ]),

            # Step 2: Quality Control
            html.Div([
                dbc.Label("2. Quality Control", className="fw-semibold text-secondary"),
                dbc.Row([
                    dbc.Col(dbc.InputGroup([
                        dbc.InputGroupText("Min. Genes/Cell"),
                        dbc.Input(id="min-genes", type="number", value=200, className="shadow-sm")
                    ]), md=6, className="mb-3 mb-md-0"),
                    dbc.Col(dbc.InputGroup([
                        dbc.InputGroupText("Min. Cells/Gene"),
                        dbc.Input(id="min-cells", type="number", value=3, className="shadow-sm")
                    ]), md=6),
                ], className="mb-3"),

                # MT Filter Section (Boxed for visual clarity)
                html.Div([
                    dbc.Row([
                        dbc.Col(
                            dbc.Checklist(
                                options=[{"label": "Enable % mt filter", "value": "on"}],
                                value=[],  # off by default
                                id="mt-filter-enabled",
                                switch=True,
                                className="fw-bold text-primary mt-2"
                            ), width="auto", className="pe-4"
                        ),
                        dbc.Col(dbc.InputGroup([
                            dbc.InputGroupText("Prefix"),
                            dbc.Input(id="mt-prefix", type="text", value="MT-", placeholder="e.g. MT-")
                        ]), md=3),
                        dbc.Col(dbc.InputGroup([
                            dbc.InputGroupText("Lower %"),
                            dbc.Input(id="mt-lower", type="number", min=0, max=100, step=0.1, placeholder="0.0", disabled=True)
                        ]), md=3),
                        dbc.Col(dbc.InputGroup([
                            dbc.InputGroupText("Upper %"),
                            dbc.Input(id="mt-upper", type="number", min=0, max=100, step=0.1, placeholder="20.0", disabled=True)
                        ]), md=3),
                    ], className="align-items-center"),
                    html.Small("Prefix marks mitochondrial genes. Toggle switch to apply filter.", className="text-muted mt-2 d-block")
                ], className="p-3 bg-light border rounded mb-4")
            ]),

            # Step 3 & 4: Normalization & Feature Selection
            dbc.Row([
                dbc.Col([
                    dbc.Label("3. Normalization", className="fw-semibold text-secondary"),
                    dbc.Select(
                        id="norm-method",
                        options=[
                            {"label": "Median", "value": "median"},
                            {"label": "CPM", "value": "cpm"},
                        ],
                        value="median", # Changed from "" to "median" for better UX default
                        className="shadow-sm mb-1"
                    ),
                    html.Small("Scales library sizes across cells.", className="text-muted d-block mb-4"),
                ], md=4),

                dbc.Col([
                    dbc.Label("4. Feature Selection", className="fw-semibold text-secondary"),
                    dbc.InputGroup([
                        dbc.InputGroupText("# Var Genes"),
                        dbc.Input(id="n-features", type="number", value=2000, className="shadow-sm")
                    ], className="mb-2"),
                    dbc.Select(
                        id="feature-method",
                        options=[
                            {"label": "VST (Recommended)", "value": "vst"},
                            {"label": "MeanVarPlot", "value": "mean.var.plot"},
                            {"label": "Dispersion", "value": "dispersion"}
                        ],
                        value="vst",
                        className="shadow-sm mb-1"
                    ),
                    html.Small("Finds high cell-to-cell variation.", className="text-muted d-block mb-4"),
                ], md=8),
            ]),

            # Step 5 & 6: PCA & UMAP
            dbc.Row([
                dbc.Col([
                    dbc.Label("5. PCA Components", className="fw-semibold text-secondary"),
                    dbc.InputGroup([
                        dbc.InputGroupText("# PCs"),
                        dbc.Input(id="n-pcs", type="number", value=30, className="shadow-sm")
                    ], className="mb-1"),
                    html.Small("Used for UMAP & Clustering.", className="text-muted d-block mb-4"),
                ], md=4),

                dbc.Col([
                    dbc.Label("6. Visualization (UMAP)", className="fw-semibold text-secondary"),
                    dbc.Row([
                        dbc.Col(dbc.InputGroup([
                            dbc.InputGroupText("Min. Dist"),
                            dbc.Input(id="umap-dist", type="number", value=0.3, step=0.1, className="shadow-sm")
                        ])),
                        dbc.Col(dbc.InputGroup([
                            dbc.InputGroupText("Neighbors"),
                            dbc.Input(id="umap-neighbors", type="number", value=30, className="shadow-sm")
                        ])),
                    ], className="mb-1"),
                    html.Small("Manifold approximation parameters.", className="text-muted d-block mb-4"),
                ], md=8),
            ]),

            # Step 7 & 8: Clustering & DE
            dbc.Row([
                dbc.Col([
                    dbc.Label("7. Clustering", className="fw-semibold text-secondary"),
                    dbc.InputGroup([
                        dbc.InputGroupText("Resolution"),
                        dbc.Input(id="resolution", type="number", value=0.8, step=0.1, className="shadow-sm")
                    ], className="mb-2"),
                    dbc.Select(
                        id="clustering-method",
                        options=[
                            {"label": "Graph-based (Louvain)", "value": "louvain"},
                            {"label": "Leiden", "value": "leiden"},
                        ],
                        value="louvain",
                        className="shadow-sm mb-4"
                    ),
                ], md=6),

                dbc.Col([
                    dbc.Label("8. Differential Expression", className="fw-semibold text-secondary"),
                    dbc.Select(
                        id="ranking-method",
                        options=[
                            {"label": "Wilcoxon Rank Sum", "value": "wilcoxon"},
                            {"label": "t-test", "value": "ttest"},
                            {"label": "Logistic Regression", "value": "logreg"},
                        ],
                        value="wilcoxon",
                        className="shadow-sm mb-4"
                    ),
                ], md=6),
            ]),

            html.Hr(className="my-3"),

            # Run Section
            dbc.Row([
                dbc.Col([
                    html.H6("Backend Engine", className="fw-bold mb-2 text-secondary"),
                    dbc.RadioItems(
                        id="tool-choice",
                        options=[{"label": " Scanpy", "value": "scanpy"}],
                        value="scanpy",
                        inline=True,
                        className="fw-bold"
                    )
                ], width="auto", className="pe-4 align-self-center"),
                
                dbc.Col([
                    dbc.Button(
                        "RUN ANALYSIS PIPELINE", 
                        id="run-button", 
                        color="primary", 
                        size="lg",
                        className="w-100 shadow border-0 rounded-3 fw-bold",
                        style={"letterSpacing": "1px"}
                    )
                ])
            ], className="mt-2 mb-2"),
            
            html.Div(id="run-status-msg", className="text-center text-info fw-bold mt-2"),

        ]),
        className="shadow-lg border-0 rounded-4 p-3 bg-white mb-5"
    )

    # --- 2. RESULTS PANEL CARD ---
    results_card = dbc.Card(
        dbc.CardBody([
            html.H4("Analysis Results", className="fw-bold mb-4"),
            html.Hr(),
            
            # Main UMAP/Plots Output
            dcc.Loading(
                id="loading-output",
                type="circle",
                children=html.Div(id="analysis-results", style={"minHeight": "200px"})
            ),
            
            # Hidden Placeholders (Kept intact for your callbacks)
            dcc.Dropdown(id="gene-dropdown", options=[], style={"display": "none"}),
            dcc.Dropdown(id="metadata-field-dropdown", options=[], style={"display": "none"}),

            html.Hr(className="my-5"),

            # Volcano Plot Section
            html.Div([
                html.H5("Differential Expression – Volcano Plot", className="fw-bold mb-3"),
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Select Cluster to View", className="text-secondary"),
                        dcc.Dropdown(id="volcano-cluster-dropdown", options=[], className="shadow-sm")
                    ], md=6)
                ], className="mb-4"),
                dcc.Graph(id="volcano-plot")
            ], id="volcano-section", style={"display": "none"}, className="p-4 bg-light rounded border mb-5"),

            # Download Section
            html.Div([
                html.H5("Export Results", className="fw-bold mb-3"),
                dbc.Row([
                    dbc.Col([
                        dcc.Dropdown(
                            id="file-download-dropdown",
                            placeholder="Select a generated result file...",
                            className="shadow-sm"
                        )
                    ], md=8),
                    dbc.Col([
                        dbc.Button("⬇️ Download File", id="file-download-button", color="info", className="w-100 shadow-sm text-white fw-bold"),
                        dcc.Download(id="file-download")
                    ], md=4),
                ], className="align-items-center p-4 bg-light rounded border")
            ]),
        ]),
        className="shadow-lg border-0 rounded-4 p-4 bg-white"
    )

    # --- 3. MAIN CONTAINER ASSEMBLY ---
    return dbc.Container([
        
        # Hidden Utilities (State Management)
        dcc.Interval(id="poll-interval", interval=3000, n_intervals=0),
        dcc.Store(id="upload-folder", storage_type="local"),
        dcc.Store(id="umap-mode", data="gene"),
        dcc.Store(id="run-started", data=False),

        # Header Text
        dbc.Row(dbc.Col(html.P("Upload 10x matrix files and run standard preprocessing, clustering, and visualization.", className="text-center text-muted fs-5 mt-2 mb-5"))),

        # Form
        dbc.Row(dbc.Col(setup_card, xs=12, md=10, lg=8, xl=8), className="justify-content-center"),

        # Results
        dbc.Row(dbc.Col(results_card, width=12), className="mt-2 mb-5")

    ], fluid=True, className="px-4")