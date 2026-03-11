import dash_bootstrap_components as dbc
from dash import html, dcc
import dash_cytoscape as cyto
from dash import dash_table

def grn_layout():
    species_options = [
        {"label": "Human", "value": "human"},
        {"label": "Mouse", "value": "mouse"},
        {"label": "Pig", "value": "pig"},
        {"label": "Rat", "value": "rat"},
        {"label": "Zebrafish", "value": "zebrafish"},
        {"label": "Chicken", "value": "chicken"},
        {"label": "C. elegans", "value": "celegans"},
        {"label": "Drosophila", "value": "drosophila"},
        {"label": "Yeast", "value": "scerevisiae"},
        {"label": "Arabidopsis", "value": "arabidopsis"},
        {"label": "Xenopus laevis", "value": "xenopus_laevis"},
        {"label": "Xenopus tropicalis", "value": "xenopus_tropicalis"},
    ]

    # --- 1. SETUP FORM CARD ---
    setup_card = dbc.Card(
        dbc.CardBody([
            dbc.Row([
                dbc.Col(html.H4("GRN Inference Setup", className="fw-bold mb-0"), width=8),
                dbc.Col(dbc.Button("🔍 Load Example", id="load-example-button", color="outline-secondary", size="sm", className="float-end"), width=4)
            ], className="mb-4 align-items-center"),

            # Step 1
            html.Div([
                dbc.Label("1. Select Organism", className="fw-semibold text-secondary"),
                dbc.Select(
                    id="grn-species-dropdown",
                    options=species_options,
                    value=None,
                    placeholder="Choose organism for base GRN...",
                    className="shadow-sm mb-4"
                )
            ]),

            # Step 2
            html.Div([
                dbc.Label("2. Define Cell Groups", className="fw-semibold text-secondary"),
                dbc.RadioItems(
                    id="grn-grouping-choice",
                    options=[
                        {"label": " Use Computed Clusters (e.g., Louvain)", "value": "numeric"},
                        {"label": " Use Uploaded Annotation (cell_type)", "value": "annotation_col"},
                    ],
                    value="numeric",
                    className="mb-2"
                ),
                # Sub-option for clustering
                html.Div([
                    html.Small("Select Algorithm (only if using clusters):", className="text-muted d-block mb-1"),
                    dbc.RadioItems(
                        id="grn-clustering-method",
                        options=[
                            {"label": "Louvain", "value": "louvain"},
                            {"label": "Leiden", "value": "leiden"},
                        ],
                        value="louvain",
                        inline=True,
                    ),
                ], className="p-3 bg-light border rounded mb-4"),
            ]),

            # Step 3
            html.Div([
                dbc.Label("3. GRN Sparsity (Alpha)", className="fw-semibold text-secondary"),
                dbc.InputGroup([
                    dbc.InputGroupText("Alpha"),
                    dbc.Input(id="grn-alpha", type="number", min=1, value=10, step=1, className="shadow-sm"),
                ], className="mb-1"),
                html.Small("Higher alpha increases sparsity (fewer edges retained).", className="text-muted d-block mb-4"),
            ]),

            # Step 4
            html.Div([
                dbc.Label("4. Link Filtering Threshold", className="fw-semibold text-secondary"),
                dbc.InputGroup([
                    dbc.InputGroupText("Top N Links"),
                    dbc.Input(id="grn-threshold-number", type="number", min=100, step=100, value=2000, className="shadow-sm"),
                ], className="mb-4"),
            ]),

            # Warning
            dbc.Alert("Note: Requires a preprocessed .h5ad file. Run Single-Cell Analysis first.", color="warning", className="small py-2 mb-4"),

            # Submit
            dbc.Button("Run GRN Analysis", id="grn-submit-button", color="primary", size="lg", className="w-100 shadow border-0 rounded-3 fw-bold")
        ]),
        className="shadow-lg border-0 rounded-4 p-2 bg-white mb-5"
    )

    # --- 2. RESULTS PANEL CARD ---
    results_card = html.Div(
        id="grn-results-panel",
        style={"display": "none"}, # Hidden by default
        children=[
            dbc.Card(
                dbc.CardBody([
                    
                    # Top Bar: Download Global Table
                    dbc.Row([
                        dbc.Col(html.H4("GRN Analysis Results", className="fw-bold mb-0"), width=8),
                        dbc.Col([
                            dbc.Button("Download Full Network (CSV)", id="export-filtered-button", color="info", className="float-end shadow-sm"),
                            dcc.Download(id="download-filtered-table")
                        ], width=4)
                    ], className="mb-4 align-items-center"),
                    html.Hr(),

                    # Centrality Plots Section
                    html.H5("Network Centrality Analysis", className="fw-bold mt-4 mb-3"),
                    dbc.Row([
                        dbc.Col([dbc.Label("Select Cluster", className="text-secondary"), dcc.Dropdown(id="grn-cluster-select", placeholder="Select cluster...", clearable=False, className="shadow-sm")]),
                        dbc.Col([dbc.Label("Centrality Metric", className="text-secondary"), dcc.Dropdown(id="grn-centrality-metric", options=[{"label": "Degree", "value": "degree_centrality_all"}, {"label": "Betweenness", "value": "betweenness_centrality"}, {"label": "Eigenvector", "value": "eigenvector_centrality"}], value="degree_centrality_all", className="shadow-sm")]),
                        dbc.Col([dbc.Label("Top N Genes", className="text-secondary"), dcc.Slider(id="grn-top-n-tfs", min=5, max=30, step=5, value=20, marks={i: str(i) for i in range(5, 35, 5)})]),
                    ], className="mb-4 bg-light p-3 rounded border"),
                    
                    dcc.Loading(dcc.Graph(id="grn-centrality-plot", style={"width": "100%", "height": "500px"})),

                    html.Hr(className="my-5"),

                    # TF Subnetwork Section
                    html.H5("Subnetwork Discovery (Target TFs)", className="fw-bold mb-3"),
                    dbc.Row([
                        dbc.Col([
                            dbc.Label("1. Input Target Genes", className="fw-semibold text-secondary"),
                            dcc.Textarea(
                                id="target-genes-input",
                                style={"width": "100%", "height": "120px", "borderRadius": "8px"},
                                className="p-2 border shadow-sm",
                                placeholder="Paste gene symbols here (one per line or comma-separated)...",
                            ),
                            dcc.Interval(id="grn-auto-trigger", interval=300, n_intervals=0, max_intervals=1),
                            dbc.Button("Extract Subnetwork", id="analyze-grn-button", color="primary", className="mt-3 shadow-sm w-100 fw-bold"),
                        ], md=4),
                        
                        dbc.Col([
                            dbc.Label("2. Top Regulating TFs", className="fw-semibold text-secondary"),
                            html.Div([
                                dash_table.DataTable(
                                    id="grn-tf-table",
                                    columns=[
                                        {"name": "TF", "id": "TF"},
                                        {"name": "#Edges", "id": "edge_count"},
                                        {"name": "%Edges", "id": "edge_pct"},
                                        {"name": "Max -logp", "id": "max_logp"},
                                    ],
                                    row_selectable="multi",
                                    style_table={"overflowX": "auto", "border": "1px solid #dee2e6", "borderRadius": "5px"},
                                    page_size=5,
                                    style_cell={"textAlign": "left", "padding": "10px", "fontFamily": "sans-serif"},
                                    style_header={"fontWeight": "bold", "backgroundColor": "#f8f9fa"},
                                )
                            ], className="shadow-sm")
                        ], md=8),
                    ], className="mb-5"),

                    # Cytoscape Viewer
                    html.H5("Interactive Subnetwork Viewer", className="fw-bold mb-3"),
                    html.Div(id="cyto-hint-msg", className="text-muted mb-2"),
                    dbc.Row([
                        dbc.Col([
                            html.Div(
                                id="edge-hover-info",
                                className="p-3 bg-light border rounded shadow-sm mb-3",
                                style={"minHeight": "60px"}
                            ),
                            html.Div(
                                cyto.Cytoscape(
                                    id="grn-network-cytoscape",
                                    layout={"name": "cose"},
                                    style={"height": "600px", "width": "100%", "backgroundColor": "#ffffff"},
                                    userZoomingEnabled=True,
                                    stylesheet=[
                                        {"selector": ".tf-node", "style": {"shape": "ellipse", "width": "50px", "height": "50px", "background-color": "#e74c3c", "label": "data(label)", "color": "#333", "font-size": "14px", "text-valign": "center"}},
                                        {"selector": ".target-node", "style": {"shape": "ellipse", "width": "40px", "height": "40px", "background-color": "#95a5a6", "label": "data(label)", "color": "#333", "font-size": "12px", "text-valign": "center"}},
                                        {"selector": "edge", "style": {"line-color": "#bdc3c7", "width": 2, "target-arrow-shape": "triangle", "curve-style": "bezier"}}
                                    ]
                                ),
                                className="border rounded shadow-sm p-2 bg-white"
                            ),
                            dbc.Button("Download Network Image (PNG)", id="download-png-button", color="secondary", outline=True, className="mt-3"),
                            dcc.Store(id="filtered-edges-store"),
                            html.Div(id="dummy-output", style={"display": "none"}),
                        ])
                    ]),

                    html.Hr(className="my-5"),

                    # TF Activity Across Clusters
                    html.H5("Inspect Specific TF Activity", className="fw-bold mb-3"),
                    dbc.Row([
                        dbc.Col([
                            dcc.Dropdown(id="grn-tf-input", placeholder="Select or type a TF to inspect across all clusters...", searchable=True, clearable=True, className="shadow-sm mb-3")
                        ], md=6)
                    ]),
                    dcc.Loading(dcc.Graph(id="grn-top-tfs-plot", style={"width": "100%", "height": "400px"})),

                    # Bottom Download
                    html.Div([
                        dbc.Button("Download All Analysis Results (ZIP)", id="grn-download-button", color="success", size="lg", className="mt-5 shadow"), 
                        dcc.Download(id="grn-download-all")
                    ], className="text-center")

                ]),
                className="shadow-lg border-0 rounded-4 p-4 bg-white"
            )
        ]
    )

    # --- 3. MAIN CONTAINER ASSEMBLY ---
    return dbc.Container([
        # Hidden Utilities
        html.Pre(id="debug-task-id", style={"display": "none"}),
        dcc.Interval(id="grn-status-poll", interval=3000, n_intervals=0, disabled=True),
        dcc.Store(id="grn-task-id"),
        dcc.Store(id="grn-species"),

        # Header
        dbc.Row(dbc.Col(html.P("Predict and visualize cell-type specific gene regulatory networks using CellOracle.", className="text-center text-muted fs-5 mt-2 mb-5"))),

        # Form
        dbc.Row(dbc.Col(setup_card, xs=12, md=10, lg=8, xl=6), className="justify-content-center"),

        # Status Message / Loader
        dbc.Row(dbc.Col([
            html.Div(id="grn-status-msg", className="text-center fs-5 text-primary fw-bold mb-3"),
            dcc.Loading(id="grn-loading", type="circle", children=html.Div(id="grn-analysis-results"))
        ])),

        # Results
        dbc.Row(dbc.Col(results_card, width=12), className="mt-2 mb-5")

    ], fluid=True, className="px-4")