# sc_callbacks.py

import threading
import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, no_update,dcc,ctx
from scanpy_pipeline import run_scanpy_analysis, plot_umap,get_volcano_data, plot_volcano
import uuid
import scanpy as sc
import os
from globals import task_status,EXAMPLE_FOLDER_PATH,UPLOAD_BASE_DIR 
import shutil
import zipfile

from dash.exceptions import PreventUpdate

import base64, io, pandas as pd

import logging
logging.basicConfig(level=logging.DEBUG)



def register_single_cell_callbacks(app):


    @app.callback(
    Output("metadata-upload-status", "children"),
    Output("metadata-field-dropdown", "options"),
    Input("upload-metadata-after", "contents"),
    State("upload-metadata-after", "filename"),
    State("adata-path", "data"),
    prevent_initial_call=True
    )
    def add_metadata_after(contents, filename, adata_path):
        from dash.exceptions import PreventUpdate
        import base64, io, os
        import pandas as pd
        import scanpy as sc

        if not contents or not adata_path:
            raise PreventUpdate

        # --- read CSV/TSV; parse common "NA" strings as missing ---
        try:
            _, content_string = contents.split(",", 1)
            decoded = base64.b64decode(content_string)
            meta = pd.read_csv(
                io.StringIO(decoded.decode("utf-8")),
                sep=None, engine="python",
                na_values=["NA", "na", "Na", "N/A", ""]
            )
        except Exception as e:
            return f"❌ Failed to read {filename}: {e}", []

        if meta.empty:
            return "❌ Metadata file appears to be empty.", []

        # --- pick the key column: 'barcode' → common aliases → first non-unnamed column ---
        preferred = ["barcode", "BARCODE", "Barcode"]
        aliases   = ["CELL_ID", "cell_id", "Cell_ID", "CellId", "cell", "barcodes", "CB", "cell_barcode"]

        cols = list(meta.columns)
        key_col = next((c for c in preferred if c in cols), None)
        used_fallback = False
        alias_used = None

        if key_col is None:
            alias_used = next((c for c in aliases if c in cols), None)
            if alias_used is not None:
                key_col = alias_used
            else:
                # first non-index-like column
                first_cols = [c for c in cols if not str(c).lower().startswith("unnamed")]
                key_col = first_cols[0]
                used_fallback = True

        # normalize to 'barcode'
        if key_col != "barcode":
            meta.rename(columns={key_col: "barcode"}, inplace=True)

        # ensure string keys
        meta["barcode"] = meta["barcode"].astype(str).str.strip()

        # --- read adata and ensure keys are comparable ---
        if not os.path.exists(adata_path):
            return "❌ Run the analysis first so 'adata.h5ad' exists.", []
        ad = sc.read(adata_path)
        ad.obs_names = ad.obs_names.astype(str)

        # (optional) normalize suffixes if your IDs differ (uncomment if needed)
        # meta["barcode"] = meta["barcode"].str.replace(r"-1$", "", regex=True)
        # ad.obs_names = ad.obs_names.str.replace(r"-1$", "", regex=True)

        # quick sanity: any overlap?
        overlap = ad.obs_names.isin(meta["barcode"]).sum()
        if overlap == 0:
            msg = f"❌ No matching cell IDs between AnnData (n={ad.n_obs}) and metadata using '{key_col}'."
            if used_fallback:
                msg += " (Used first column as IDs.)"
            elif alias_used:
                msg += f" (Using alias column '{alias_used}'.)"
            return msg, []

        # avoid name collisions
        collisions = [c for c in meta.columns if c != "barcode" and c in ad.obs.columns]
        if collisions:
            meta = meta.rename(columns={c: f"{c}_meta" for c in collisions})

        # merge by barcode
        ad.obs["__barcode_tmp__"] = ad.obs_names
        ad.obs = ad.obs.merge(
            meta.set_index("barcode"),
            left_on="__barcode_tmp__", right_index=True, how="left"
        )
        ad.obs.drop(columns=["__barcode_tmp__"], inplace=True)

        # tidy dtypes: numeric stays numeric; other objects → category
        for c in ad.obs.columns:
            newcol = pd.to_numeric(ad.obs[c], errors="ignore")
            if newcol.dtype == "object":
                newcol = newcol.astype("category")
            ad.obs[c] = newcol

        ad.write(adata_path)

        # build dropdown options (everything except 'barcode')
        fields = [c for c in meta.columns if c != "barcode"]
        opts = [{"label": f, "value": f} for f in fields]

        # status message
        if used_fallback:
            status = f"✅ Loaded metadata from {filename}. Using FIRST column '{key_col}' as barcode. " \
                    f"Matched {overlap} cells. " \
                    + (f"(Renamed collisions: {', '.join(collisions)} → *_meta*)" if collisions else "")
        elif alias_used:
            status = f"✅ Loaded metadata from {filename}. Using alias column '{alias_used}' as barcode. " \
                    f"Matched {overlap} cells. " \
                    + (f"(Renamed collisions: {', '.join(collisions)} → *_meta*)" if collisions else "")
        else:
            status = f"✅ Loaded metadata from {filename}. Using 'barcode'. Matched {overlap} cells. " \
                    + (f"(Renamed collisions: {', '.join(collisions)} → *_meta*)" if collisions else "")

        return status, opts



    def background_scanpy_runner(task_id, kwargs):
        try:
            result_dir, pca_fig, umap_fig, adata_path = run_scanpy_analysis(**kwargs)

            adata = sc.read(adata_path)
            gene_options = [{"label": g, "value": g} for g in adata.var_names]
            # After adata is read
            rank_info = adata.uns.get("rank_genes_groups", {})
            cluster_options = []

            if rank_info:
                cluster_names = rank_info["names"].dtype.names  # tuple of group names
                cluster_options = [{"label": f"Cluster {g}", "value": g} for g in cluster_names]


            task_status[task_id] = {
                "status": "done",
                "result_dir": result_dir,
                "pca_fig": pca_fig,
                "umap_fig": umap_fig,
                "adata_path": adata_path,
                "gene_options": gene_options,
                "cluster_options": cluster_options,
            }
        except Exception as e:
            task_status[task_id] = {
                "status": "error",
                "message": str(e)
            }

    @app.callback(Output("upload-status", "children"), Input("upload-folder", "data"))
    def check_folder_valid(folder):
        if not folder:
            return "📂 No folder uploaded yet."
        elif not os.path.exists(folder):
            return f"❌ Folder not found: {folder}. Please re-upload."
        return f"✅ Stored upload folder: {folder}"


    @app.callback(
    Output("mt-lower", "disabled"),
    Output("mt-upper", "disabled"),
    Input("mt-filter-enabled", "value"),
    )
    def _toggle_mt_inputs(enabled_vals):
        enabled = "on" in (enabled_vals or [])
        return (not enabled, not enabled)

    @app.callback(
        Output("run-status-msg", "children"),
        Output("run-started", "data"),
        Input("run-button", "n_clicks"),
        State("upload-folder", "data"),
        State("min-genes", "value"),
        State("min-cells", "value"),
        State("norm-method", "value"),
        State("n-features", "value"),
        State("feature-method", "value"),
        State("n-pcs", "value"),
        State("umap-dist", "value"),
        State("umap-neighbors", "value"),
        State("resolution", "value"),
        State("clustering-method", "value"),
        State("ranking-method", "value"),
        State("tool-choice", "value"),
        State("mt-filter-enabled", "value"),
        State("mt-lower", "value"),
        State("mt-upper", "value"),
        State("mt-prefix", "value"),
        prevent_initial_call=True
    )
    def start_background_job(n_clicks, folder, min_genes, min_cells, norm_method,
                             n_feats, feat_method, n_pcs, umap_dist, umap_neighbors, 
                             resolution, cluster_method, rank_method, tool_choice,
                             mt_filter_enabled_vals, mt_lower, mt_upper, mt_prefix):

        if not n_clicks or not folder:
            raise dash.exceptions.PreventUpdate

        task_id = str(uuid.uuid4())
        task_status[task_id] = {"status": "running"}
        mt_enabled = "on" in (mt_filter_enabled_vals or [])

        thread_args = {
            "dataset_path": folder,
            "min_genes_per_cell": min_genes,
            "min_cells_per_gene": min_cells,
            "norm_method": norm_method,
            "n_variable_genes": n_feats,
            "feature_selection_method": feat_method,
            "n_pcs": n_pcs,
            "umap_min_dist": umap_dist,
            "umap_n_neighbors": umap_neighbors,
            "clustering_resolution": resolution,
            "clustering_method": cluster_method,
            "ranking_method": rank_method,
            "tool_choice": tool_choice,
            "debug": False,
            "mt_filter_enabled": mt_enabled,
            "mt_lower": mt_lower,        # may be None if disabled
            "mt_upper": mt_upper,        # may be None if disabled
            "mt_prefix": mt_prefix or "MT-",
        }

        threading.Thread(
            target=background_scanpy_runner,
            kwargs={"task_id": task_id, "kwargs": thread_args},
            daemon=True
        ).start()

        return f"🔄 Analysis started. Task ID: {task_id}", task_id

    @app.callback(
        Output("run-button", "children"),
        Output("analysis-results", "children"),
        Output("gene-dropdown", "options"),
        Output("adata-path", "data"),
        Output("poll-interval", "disabled"),
        Output("volcano-cluster-dropdown", "options"),
        Output("volcano-section", "style"),
        Output("file-download-dropdown", "options"),
        Input("run-started", "data"),
        Input("poll-interval", "n_intervals"),
        State("run-button", "n_clicks"),
        State("upload-folder", "data"),
        State("min-genes", "value"),
        State("min-cells", "value"),
        State("norm-method", "value"),
        State("n-features", "value"),
        State("feature-method", "value"),
        State("n-pcs", "value"),
        State("umap-dist", "value"),
        State("umap-neighbors", "value"),
        State("resolution", "value"),
        State("clustering-method", "value"),
        State("ranking-method", "value"),
        State("tool-choice", "value"),
        prevent_initial_call=True,
        allow_duplicate=True
    )
    def on_run_click(started, n_intervals, n_clicks, folder, min_genes, min_cells, norm_method,
                     n_feats, feat_method, n_pcs, umap_dist, umap_neighbors, 
                     resolution, cluster_method, rank_method, tool_choice):

        if not started or not n_clicks or not folder:
            raise dash.exceptions.PreventUpdate

        task_id = started
        result = task_status.get(task_id)

        if not result:
            return ..., no_update, no_update, no_update, False, [], {"display": "none"},[]


        if result.get("status") == "error":
            return (
                f"❌ Error: {result['message']}",
                f"<div style='color:red;'>Error: {result['message']}</div>",
                [],
                None,
                True,
                [], {"display": "none"},
                []
            )

        if result.get("status") != "done":
            return "⏳ Still processing...", no_update, no_update, no_update, False, no_update, no_update,[]

        try:
            file_options = [
                {"label": f, "value": f}
                for f in os.listdir(result["result_dir"])
                if not f.endswith(".log")
            ]
        except Exception as e:
            print(f"❌ Failed to read result_dir: {e}")
            file_options = []
            

        return (
            f"✅ Done — Results saved to {result['result_dir']}",
            [
                dcc.Graph(figure=result["pca_fig"], id="pca-plot"),
                html.Hr(),
                html.H5("UMAP Visualization"),

                # --- New: metadata upload + instructions ---
                html.Div([
                    html.H6("Add/Update Cell Metadata"),
                    dcc.Upload(
                        id="upload-metadata-after",
                        children=html.Div(["Drag & drop or ", html.A("select a CSV/TSV")]),
                        accept=".csv,.tsv,.txt",
                        multiple=False,
                        className="border rounded p-3 text-center mb-2"
                    ),
                    html.Small(
                        "To color cells by metadata (categories or numeric values), upload a file with a 'barcode' column "
                        "matching barcodes.tsv. Then pick a field below.",
                        className="text-muted"
                    ),
                    html.Div(id="metadata-upload-status", className="text-success mt-1 mb-3"),
                ]),

                dbc.Row([
                    dbc.Col([
                        html.Label("Select a gene"),
                        dcc.Dropdown(id="gene-dropdown", placeholder="Select a gene...", className="mb-3", options=result["gene_options"])
                    ], width=8),
                    dbc.Col([
                        html.Label("Toggle UMAP View"),
                        dbc.Button("Switch to Cluster View", id="umap-toggle-button", color="secondary", className="w-100", n_clicks=0)
                    ], width=4),

                    html.Div([
                        html.Small("Optional: color by metadata — upload file above first.", className="text-muted"),
                        dcc.Dropdown(id="metadata-field-dropdown", options=[], placeholder="Choose a metadata field…", className="mb-2")
                    ]),
                ]),

                dcc.Graph(id="umap-plot")
            ],
            result["gene_options"],
            result["adata_path"],
            True,
            result.get("cluster_options", []),
            {"display": "block"},
            file_options
        )

    @app.callback(
    Output("umap-plot", "figure"),
    Input("gene-dropdown", "value"),
    Input("umap-mode", "data"),
    Input("metadata-field-dropdown", "value"),   # ← add this line
    State("adata-path", "data")
    )

    def update_umap_plot_view(gene_name, mode, meta_field, adata_path):
        if not adata_path:
            raise dash.exceptions.PreventUpdate

        
        if not adata_path:
            raise dash.exceptions.PreventUpdate

        adata = sc.read(adata_path)

        # ✅ metadata wins if selected
        if meta_field:
            return plot_umap(adata, color_by=meta_field)


        adata = sc.read(adata_path)
        if mode == "gene":
            if not gene_name:
                raise dash.exceptions.PreventUpdate
            return plot_umap(adata, color_by=gene_name)
        else:
            cluster_key = "leiden" if "leiden" in adata.obs else adata.obs.select_dtypes("category").columns[0]
            return plot_umap(adata, color_by=cluster_key)


    @app.callback(
        Output("volcano-plot", "figure"),
        Input("volcano-cluster-dropdown", "value"),
        State("adata-path", "data"),
        prevent_initial_call=True
    )
    def update_volcano_plot_view(cluster_id, adata_path):
        if not adata_path or not cluster_id:
            raise dash.exceptions.PreventUpdate

        adata = sc.read(adata_path)
        df = get_volcano_data(adata, cluster_id)
        return plot_volcano(df, cluster_id)



    @app.callback(
    Output("upload-folder", "data", allow_duplicate=True),
    Output("run-started", "data", allow_duplicate=True),
    Input("load-example-singlecell-button", "n_clicks"),
    prevent_initial_call=True
    )
    def trigger_example_load(n_clicks):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate

        # Use configured path
        folder_path = EXAMPLE_FOLDER_PATH

        task_id = str(uuid.uuid4())
        task_status[task_id] = {"status": "running"}

        thread_args = {
            "dataset_path": folder_path,
            "min_genes_per_cell": 200,
            "min_cells_per_gene": 3,
            "norm_method": "log1p",
            "n_variable_genes": 2000,
            "feature_selection_method": "seurat",
            "n_pcs": 20,
            "umap_min_dist": 0.3,
            "umap_n_neighbors": 15,
            "clustering_resolution": 0.5,
            "clustering_method": "leiden",
            "ranking_method": "wilcoxon",
            "tool_choice": "scanpy",
            "debug": False
        }

        threading.Thread(
            target=background_scanpy_runner,
            kwargs={"task_id": task_id, "kwargs": thread_args},
            daemon=True
        ).start()

        return folder_path, task_id





    @app.callback(
    Output("file-download-dropdown", "options", allow_duplicate=True),  # 👈 allow duplicates
    Input("run-started", "data"),
    prevent_initial_call=True  # 👈 also helps avoid early triggering
    )
    def list_downloadable_files(task_id):
        print(f"📂 Populating dropdown for task_id: {task_id}")

        if not task_id:
            print("❌ No task_id")
            return []

        task = task_status.get(task_id)
        if not task or task.get("status") != "done":
            print("⏳ Task not done yet")
            return []

        result_dir = task.get("result_dir")
        print(f"📁 Scanning directory: {result_dir}")

        try:
            files = os.listdir(result_dir)
            csv_files = [f for f in files if f.endswith(".csv")]
            options = [{"label": f, "value": f} for f in csv_files]
            print(f"✅ Found files: {files}")
            return options
        except Exception as e:
            print(f"❌ Error reading result_dir: {e}")
            return []



    @app.callback(
    Output("file-download", "data"),
    Input("file-download-button", "n_clicks"),
    State("file-download-dropdown", "value"),
    State("run-started", "data"),
    prevent_initial_call=True
    )
    def serve_selected_file(n_clicks, selected_file, task_id):
        if not selected_file or not task_id:
            raise PreventUpdate

        task = task_status.get(task_id)
        if not task:
            raise PreventUpdate

        result_dir = task.get("result_dir")
        file_path = os.path.join(result_dir, selected_file)

        if os.path.exists(file_path):
            print(f"🚀 Serving file: {file_path}")
            return dcc.send_file(file_path)
        else:
            print("❌ File not found")
            raise PreventUpdate




    @app.callback(
    Output("sc-download-container", "style"),
    Input("run-started", "data"),
    )
    def show_download_button(task_id):
        if task_id:
            return {"display": "block"}
        return {"display": "none"}









    @app.callback(
        Output("fresh-start-status", "children"),
        Input("fresh-start-button", "n_clicks"),
        State("upload-folder", "data"),  # You must have this store set somewhere
        prevent_initial_call=True
    )
    def cleanup_uploaded_data(n_clicks, session_folder):
        if not session_folder:
            return "⚠️ No session folder available."


        if session_folder == EXAMPLE_FOLDER_PATH:
            return "ℹ️ Example data is not deleted."

        if not os.path.exists(session_folder):
            return f"⚠️ Folder not found: {session_folder}"

        try:
            shutil.rmtree(session_folder)
            return f"✅ Deleted folder: {session_folder}"
        except Exception as e:
            return f"❌ Failed to delete folder: {e}"


    @app.callback(
        Output("umap-mode", "data"),
        Output("umap-toggle-button", "children"),
        Input("umap-toggle-button", "n_clicks"),
        State("umap-mode", "data"),
        prevent_initial_call=True
    )
    def toggle_umap_mode(n_clicks, current_mode):
        if current_mode == "gene":
            return "cluster", "Switch to Gene Expression View"
        return "gene", "Switch to Cluster View"

