import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated.*")

import os
import logging
from datetime import datetime
import threading
import uuid

import dash
import dash_bootstrap_components as dbc
from dash import html, dcc, Input, Output, State, no_update
from flask import Flask, request, jsonify
from flask_session import Session

from grn_layout import grn_layout
from single_cell import single_cell_layout
from scanpy_pipeline import cleanup_folders

# Externalized callbacks
from grn_callbacks import register_grn_callbacks
from sc_callbacks import register_single_cell_callbacks

# Basic logger setup
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("pipeline.log")
    ]
)
logger = logging.getLogger(__name__)

# Create Flask server
server = Flask(__name__)
server.secret_key = "123456"
server.config["SESSION_TYPE"] = "filesystem"
Session(server)

UPLOAD_FOLDER = "uploaded_data"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
logger.info("🚀 Flask/Dash app started.")

# Flask route for file upload
@server.route("/scom/upload", methods=["POST"])
def upload_files():
    try:
        all_files = [file for key, file in request.files.items() if key.startswith("file[")]
        if not all_files:
            return "No files uploaded", 400

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        folder = os.path.join("uploaded_data", timestamp)
        os.makedirs(folder, exist_ok=True)

        saved = []
        for file in all_files:
            if file and file.filename:
                path = os.path.join(folder, file.filename)
                file.save(path)
                saved.append(file.filename)

        return jsonify({"message": "Upload successful", "saved_files": saved, "session_folder": folder}), 200

    except Exception as e:
        return f"❌ Upload failed: {e}", 500

# Initialize Dash app
app = dash.Dash(
    __name__,
    server=server,
    external_stylesheets=[dbc.themes.FLATLY], # FLATLY gives it a super clean, modern font/color palette
    url_base_pathname='/scom/',
    suppress_callback_exceptions=True
)
app.title = "Survom - Single-Cell Analysis"

# Layout
app.layout = html.Div([
    # Hidden Stores for State Management
    dcc.Store(id="adata-path"),
    dcc.Store(id='grn-summary-store'),
    dcc.Store(id="grn-analysis-complete"),
    dcc.Store(id="grn-merged-score"),
    dcc.Store(id="grn-links-store"),

    # --- 1. Top Navigation/Header Bar ---
    dbc.Navbar(
        dbc.Container([
            dbc.Row([
                dbc.Col(html.Img(src="assets/bhu_logo.png", height="60px"), width="auto"),
                dbc.Col(
                    html.H2("Survom Analysis Pipeline", className="text-white fw-bold mb-0 mx-4"), 
                    width="auto", className="d-flex align-items-center"
                ),
                dbc.Col(html.Img(src="assets/umea_logo.png", height="60px"), width="auto"),
            ], align="center", className="g-0 w-100 justify-content-center"),
        ], fluid=True),
        color="primary", # Uses the theme's main color
        dark=True,       # Makes the text/icons white to contrast the background
        className="shadow-sm mb-4 py-3" # Adds a drop shadow and vertical padding
    ),

    # --- 2. Main Content Container ---
    dbc.Container([
        
        # Modern Bootstrap Tabs (Pills)
        dbc.Tabs(
            [
                dbc.Tab(label="Single Cell Analysis", tab_id="single-cell", label_class_name="fw-bold px-4"),
                dbc.Tab(label="Gene Regulatory Network", tab_id="grn", label_class_name="fw-bold px-4"),
            ],
            id="analysis-tabs",
            active_tab="single-cell", 
            className="mb-4 justify-content-center nav-pills" # nav-pills makes them look like modern buttons
        ),

        # Content area where the layouts will load
        html.Div(id="tab-content", className="p-4 bg-white shadow-sm rounded border")
        
    ], fluid=True, className="pb-5") # fluid=True uses full screen width; pb-5 adds padding at the bottom
    
], className="bg-light", style={"minHeight": "100vh"}) # Sets the entire background to a subtle gray, making the white cards pop

# Tab content switcher
@app.callback(
    Output("tab-content", "children"),
    Input("analysis-tabs", "active_tab") # CRITICAL: dbc.Tabs uses 'active_tab', not 'value'
)
def display_tab_content(tab_name):
    if tab_name == "grn":
        return grn_layout()
    elif tab_name == "single-cell":
        return single_cell_layout()
    return html.Div(
        [html.H4("🔍 Unknown tab selected", className="text-muted")], 
        className="text-center mt-5"
    )

# Register modular callbacks
def register_all_callbacks(app):
    register_single_cell_callbacks(app)
    register_grn_callbacks(app)

# Trigger cleanup threads on startup
if __name__ == "__main__":
    register_all_callbacks(app)

    threading.Thread(target=cleanup_folders, kwargs={
        "base_folder": "uploaded_data", "interval_minutes": 10, "expiry_minutes": 30
    }, daemon=True).start()

    threading.Thread(target=cleanup_folders, kwargs={
        "base_folder": "output", "interval_minutes": 10, "expiry_minutes": 30
    }, daemon=True).start()

    app.run(host="0.0.0.0", debug=True, port=8050, use_reloader=False)