import uuid
import dash
import threading
from dash import Input, Output, State, html, no_update,dcc,ctx
from run_grn_analysis import run_grn_analysis
from globals import grn_tasks  # if you move this to globals.py or similar
import plotly.graph_objects as go
import pandas as pd 
from dash import no_update
from grn_view import (
    plot_tf_centrality_across_clusters,
    plot_cluster_centrality_from_df
)
import dash_table
from dash import callback
import numpy as np
import os 
import zipfile
# PRECOMPUTED_GRN_SUMMARY = {
#     "oracle_path": "/data/users/vikash/scom/tamalika_grn/a6c64529-9604-4106-bcc0-ccd93c2a745d/oracle.celloracle.oracle",
#     "raw_links_path": "/data/users/vikash/scom/tamalika_grn/a6c64529-9604-4106-bcc0-ccd93c2a745d/raw_links.celloracle.links",
#     "filtered_links_path": "/data/users/vikash/scom/tamalika_grn/a6c64529-9604-4106-bcc0-ccd93c2a745d/filtered_links.celloracle.links",
# }

PRECOMPUTED_GRN_SUMMARY = {
    "oracle_path": "/app/tamalika_grn/a6c64529-9604-4106-bcc0-ccd93c2a745d/oracle.celloracle.oracle",
    "raw_links_path": "/app/tamalika_grn/a6c64529-9604-4106-bcc0-ccd93c2a745d/raw_links.celloracle.links",
    "filtered_links_path": "/app/tamalika_grn/a6c64529-9604-4106-bcc0-ccd93c2a745d/filtered_links.celloracle.links",
}

from dash.exceptions import PreventUpdate

import celloracle as co

# Failsafe: Try to load precomputed clusters, but don't crash if file missing
try:
    if os.path.exists(PRECOMPUTED_GRN_SUMMARY["filtered_links_path"]):
        links = co.load_hdf5(PRECOMPUTED_GRN_SUMMARY["filtered_links_path"])
        # print(list(links.links_dict.keys()))
        PRECOMPUTED_GRN_SUMMARY["clusters"] = list(map(str, links.links_dict.keys()))
    else:
        PRECOMPUTED_GRN_SUMMARY["clusters"] = []
except Exception as e:
    print(f"⚠️ Could not preload example clusters: {e}")
    PRECOMPUTED_GRN_SUMMARY["clusters"] = []

from dash import ClientsideFunction

import logging
logger = logging.getLogger(__name__)




def register_grn_callbacks(app):
    def grn_background_runner(task_id, h5ad_path, species, clustering_method, alpha):
        try:
            output_dir = f"output/{task_id}"
            summary = run_grn_analysis(h5ad_path, species, clustering_method, alpha, output_dir)
            logger.info(f"✅ grn_background_runner: {h5ad_path}")
            
            grn_tasks[task_id] = {"status": "done", "output": summary}
        except Exception as e:
            grn_tasks[task_id] = {"status": "error", "error": str(e)}

    # @app.callback(
    #     Output("grn-status-msg", "children"),
    #     Output("grn-task-id", "data"),
    #     Input("grn-submit-button", "n_clicks"),
    #     State("grn-species-dropdown", "value"),
    #     State("grn-clustering-method", "value"),
    #     State("grn-alpha", "value"),
    #     State("adata-path", "data"),
    #     prevent_initial_call=True
    # )
    # def submit_grn_job(n_clicks, species, clustering_method, alpha, h5ad_path):
    #     if not h5ad_path or not species or not clustering_method:
    #         return "❌ Missing input parameters", no_update

    #     task_id = str(uuid.uuid4())
    #     grn_tasks[task_id] = {"status": "running"}

    #     threading.Thread(
    #         target=grn_background_runner,
    #         args=(task_id, h5ad_path, species, clustering_method, alpha),
    #         daemon=True
    #     ).start()

    #     return f"🔄 GRN job started (Task ID: {task_id})", task_id


    # @app.callback(
    # Output("grn-status-msg", "children"),       # This will be set to an empty string
    # Output("grn-task-id", "data"),              # Still sets the task ID for internal use
    # Input("grn-submit-button", "n_clicks"),
    # prevent_initial_call=True
    # )
    # def submit_grn_job(n_clicks):
    #     task_id = "a6c64529-9604-4106-bcc0-ccd93c2a745d"

    #     # Set task as completed using precomputed summary
    #     grn_tasks[task_id] = {
    #         "status": "done",
    #         "output": PRECOMPUTED_GRN_SUMMARY
    #     }

    #     return "", task_id  # ❌ No visible message, but ✅ task_id still stored

    @app.callback(
    Output("debug-task-id", "children"),
    Input("grn-task-id", "data"),
    prevent_initial_call=True
    )
    def show_task_id(task_id):
        print("🧠 Task ID updated in UI:", task_id)
        return f"Current task ID: {task_id}"

    @app.callback(
    Output("grn-status-msg", "children",allow_duplicate=True),
    Output("grn-task-id", "data",allow_duplicate=True),
    Output("grn-status-poll", "disabled", allow_duplicate=True),
    Input("grn-submit-button", "n_clicks"),
    Input("load-example-button", "n_clicks"),
    State("grn-species-dropdown", "value"),
    State("grn-clustering-method", "value"),
    State("grn-grouping-choice", "value"),   # <--- MODIFIED: Added this State
    State("grn-alpha", "value"),
    State("adata-path", "data"),
    prevent_initial_call=True
    )
    def submit_or_load_grn(submit_clicks, example_clicks, species, clustering_method, grouping_choice, alpha, h5ad_path):
        from dash import ctx

        triggered = ctx.triggered_id
        if triggered == "load-example-button":
            # Load precomputed task for example
            task_id = "a6c64529-9604-4106-bcc0-ccd93c2a745d"
            grn_tasks[task_id] = {
                "status": "done",
                "output": PRECOMPUTED_GRN_SUMMARY
            }
            return "", task_id,False

        elif triggered == "grn-submit-button":
            # Validate input
            if not h5ad_path or not species:
                return "❌ Missing input parameters", dash.no_update,True

            # --- MODIFIED SECTION START ---
            # Determine target column based on user toggle
            if grouping_choice == "annotation_col":
                target_column = "annotation"  # Hardcoded standard name
                print("👉 User selected Annotation mode. Using 'annotation' column.")
            else:
                if not clustering_method:
                    return "❌ Select a clustering method (Leiden/Louvain)", dash.no_update, True
                target_column = clustering_method
                print(f"👉 User selected Cluster mode. Using '{target_column}' column.")
            # --- MODIFIED SECTION END ---

            
            print("👀 Clicks:", submit_clicks, example_clicks)
            print("📦 h5ad_path:", h5ad_path)
            print("📦 species:", species)
            print("📦 target_column:", target_column) # Changed from clustering_method
            
            # Start new GRN analysis
            task_id = str(uuid.uuid4())
            grn_tasks[task_id] = {"status": "running"}

            threading.Thread(
                target=grn_background_runner,
                args=(task_id, h5ad_path, species, target_column, alpha), # Passed target_column
                daemon=True
            ).start()

            return f"🔄 GRN job started using '{target_column}' (Task ID: {task_id})\nThis may take up to 1 hour. Please keep the window open.", task_id,False

        # No known trigger
        raise dash.exceptions.PreventUpdate


    @app.callback(
    Output("grn-summary-store", "data", allow_duplicate=True),
    Output("grn-analysis-complete", "data", allow_duplicate=True),
    Output("grn-merged-score", "data", allow_duplicate=True),
    Output("grn-status-poll", "disabled"),  # to stop polling
    Input("grn-status-poll", "n_intervals"),
    State("grn-task-id", "data"),
    prevent_initial_call=True
    )
    def poll_grn_status(n_intervals, task_id):
        print(f"📡 Polling @ interval {n_intervals} for task ID: {task_id}")

        result = grn_tasks.get(task_id)
        # print("📦 Task lookup result:", result)

        if not result:
            print("❌ No task found for task_id")
            return no_update, False, no_update, False

        if result["status"] == "done":
            summary = result["output"]
            print("✅ Task is done. Summary before enrichment:", summary)

            try:
                links = co.load_hdf5(summary["filtered_links_path"])
                links.get_network_score()
            except Exception as e:
                print(f"❌ Error loading filtered links: {e}")
                return {"error": f"Failed to load links: {e}"}, False, no_update, True

            # Add cluster information
            try:
                summary["clusters"] = list(map(str, links.links_dict.keys()))
                print("📊 Clusters extracted and added to summary:", summary["clusters"])
            except Exception as e:
                print(f"⚠️ Failed to extract clusters: {e}")
                summary["clusters"] = []

            # Prepare merged score
            try:
                merged_df = links.merged_score.copy()
                merged_df["tf"] = merged_df.index
                merged_json = merged_df.to_json(date_format="iso", orient="split")
                # print("📈 Merged score shape:", merged_df.shape)
            except Exception as e:
                print(f"❌ Failed to compute merged score: {e}")
                return {"error": f"Failed to compute score: {e}"}, False, no_update, True

            print("🚀 Returning updated summary and merged data")
            return summary, True, merged_json, True  # ✅ Stop polling

        elif result["status"] == "error":
            print("❌ Task errored:", result.get("error"))
            return {"error": result["error"]}, False, no_update, True

        print("⏳ Task still running")
        return no_update, False, no_update, False


#     @app.callback(
#     Output("grn-status-msg", "children"),
#     Output("grn-task-id", "data"),
#     Input("grn-submit-button", "n_clicks"),
#     Input("load-example-button", "n_clicks"),
#     prevent_initial_call=True
# )
#     def submit_grn_job(submit_clicks, example_clicks):
#         from dash import ctx

#         if ctx.triggered_id == "load-example-button":
#             # Use precomputed task & summary
#             task_id = "a6c64529-9604-4106-bcc0-ccd93c2a745d"
#             grn_tasks[task_id] = {
#                 "status": "done",
#                 "output": PRECOMPUTED_GRN_SUMMARY
#             }
#             return "", task_id

#         # Otherwise, handle real GRN run (if needed)
#         raise dash.exceptions.PreventUpdate



    @app.callback(
    Output("target-genes-input", "value", allow_duplicate=True),
    Output("grn-auto-trigger", "n_intervals", allow_duplicate=True),
    Input("load-example-button", "n_clicks"),
    prevent_initial_call=True
    )
    def load_example_and_trigger_analysis(n_clicks):
        default_genes = """MAFF
            KLF7
            NR4A1
            FLI1
            RORA
            SPI1
            NR4A2
            A1BG"""

        return default_genes, 1



        # @app.callback(
    #     Output("debug-output", "children"),
    #     Input("grn-summary-store", "data")  # <-- changed
    # )
    # def show_debug(data):
    #     import json
    #     return json.dumps(data, indent=2)

        
    

    # @app.callback(
    # Output("grn-analysis-results", "children", allow_duplicate=True),
    # Output("grn-summary-store", "data", allow_duplicate=True),
    # Output("grn-analysis-complete", "data", allow_duplicate=True),
    # Output("grn-merged-score", "data", allow_duplicate=True),
    # Input("grn-task-id", "data"),                     # for example
    # Input("grn-analysis-complete", "data"),           # for real run
    # prevent_initial_call=True
    # )
    # def display_grn_results(task_id, is_complete):
    #     from dash import ctx
    #     trigger = ctx.triggered_id
    #     print(f"🧪 display_grn_results triggered via {trigger} | Task ID: {task_id} | Complete: {is_complete}")

    #     if not task_id:
    #         raise dash.exceptions.PreventUpdate

    #     result = grn_tasks.get(task_id)
    #     if not result:
    #         raise dash.exceptions.PreventUpdate

    #     if trigger == "grn-analysis-complete" and not is_complete:
    #         raise dash.exceptions.PreventUpdate

    #     if result.get("status") != "done":
    #         raise dash.exceptions.PreventUpdate

    #     summary = result["output"]

    #     try:
    #         links = co.load_hdf5(summary["filtered_links_path"])
    #         links.get_network_score()
    #         merged_df = links.merged_score.copy()
    #         merged_df["tf"] = merged_df.index
    #         merged_json = merged_df.to_json(date_format="iso", orient="split")
    #     except Exception as e:
    #         return html.Div(f"❌ Failed to load results: {e}"), dash.no_update, False, dash.no_update

    #     return html.Div("✅ GRN analysis complete."), summary, True, merged_json

    @app.callback(
    Output("grn-analysis-results", "children", allow_duplicate=True),
    Input("grn-task-id", "data"),                    # triggered by example load
    Input("grn-analysis-complete", "data"),          # triggered after real run
    prevent_initial_call=True
    )
    def display_grn_results(task_id, is_complete):
        from dash import ctx
        trigger = ctx.triggered_id
        # print(f"🧪 display_grn_results triggered via: {trigger} | Task ID: {task_id} | Complete flag: {is_complete}")

        if not task_id:
            raise dash.exceptions.PreventUpdate

        result = grn_tasks.get(task_id)
        if not result:
            # print(f"⚠️ No result found for Task ID: {task_id}")
            raise dash.exceptions.PreventUpdate

        if trigger == "grn-analysis-complete" and not is_complete:
            # print("⚠️ Analysis not marked complete yet.")
            raise dash.exceptions.PreventUpdate

        if result.get("status") != "done":
            # print("⚠️ Task is not yet completed.")
            raise dash.exceptions.PreventUpdate

        # print("✅ GRN analysis result is ready.")
        return html.Div("✅ GRN analysis complete.")


    @app.callback(
    Output("grn-results-panel", "style"),
    Input("grn-task-id", "data"),
    Input("grn-analysis-complete", "data"),  # 👈 added
    prevent_initial_call=True
    )
    def show_grn_results_panel(task_id, is_complete):
        # print("🧪 show_grn_results_panel | Task ID:", task_id)
        # print("📦 is_complete:", is_complete)

        result = grn_tasks.get(task_id)
        if result and result.get("status") == "done" and is_complete:
            return {"display": "block"}

        return {"display": "none"}



    # @app.callback(
    # Output("grn-links-table", "data"),
    # Output("grn-links-table", "columns"),
    # Input("grn-summary-store", "data"),  # <-- changed
    # Input("grn-cluster-select", "value"),
    # Input("grn-tf-filter", "value"),
    # prevent_initial_call=True
    # )
    # def update_links_table(summary, selected_cluster, tf_filter):
    #     if not summary or not selected_cluster:
    #         return [], []

    #     selected_cluster = str(selected_cluster)  # ensure key matches
    #     links = co.load_hdf5(summary["filtered_links_path"])
    #     df = links.links_dict[selected_cluster].reset_index()
    #     if tf_filter and "source" in df.columns:
    #         df = df[df["source"].str.contains(tf_filter, case=False, na=False)]

    #     columns = [{"name": col, "id": col} for col in df.columns]
    #     data = df.to_dict("records")

    #     return data, columns




    @app.callback(
    Output("grn-cluster-select", "value"),
    Input("grn-cluster-select", "options"),
    prevent_initial_call=True
    )
    def set_default_cluster(options):
        if options:
            return options[0]["value"]
        return None

    @app.callback(
    Output("grn-cluster-select", "options"),
    Input("grn-summary-store", "data"),
    prevent_initial_call=True
    )
    def populate_cluster_options(summary):
        # print("📦 summary received in cluster options:", summary)
        if not summary or "clusters" not in summary:
            return []

        return [{"label": str(c), "value": str(c)} for c in summary["clusters"]]




    @app.callback(
        Output("grn-download-all", "data"),
        Input("grn-download-button", "n_clicks"),
        State("grn-task-id", "data"),
        prevent_initial_call=True
    )
    def download_grn_results(n_clicks, task_id):
        if not n_clicks or not task_id:
            raise PreventUpdate  # prevents download unless button was clicked

        task = grn_tasks.get(task_id)
        if not task or "output" not in task:
            raise PreventUpdate

        summary = task["output"]
        output_dir = os.path.dirname(summary["oracle_path"])
        zip_path = os.path.join(output_dir, f"grn_output_{task_id}.zip")

        if not os.path.exists(zip_path):
            with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
                for root, _, files in os.walk(output_dir):
                    for file in files:
                        file_path = os.path.join(root, file)
                        arcname = os.path.relpath(file_path, start=output_dir)
                        zipf.write(file_path, arcname)

        return dcc.send_file(zip_path)




    

    @app.callback(
        Output("grn-tf-table", "data"),
        Output("filtered-edges-store", "data"),
        Output("cyto-hint-msg", "children"),
        Input("analyze-grn-button", "n_clicks"),
        Input("grn-auto-trigger", "n_intervals"),  # 👈 Add this line
        State("target-genes-input", "value"),
        State("grn-links-store", "data"),
        prevent_initial_call=True
    )
    def analyze_grn(n_clicks, auto_trigger, gene_text, links_json):
        # Optional: restrict to just valid trigger
        if ctx.triggered_id not in ["analyze-grn-button", "grn-auto-trigger"]:
            raise dash.exceptions.PreventUpdate

        if not gene_text or not links_json:
            return [], None, "ℹ️ Please input genes and ensure GRN data is loaded."

        target_genes = [g.strip() for g in gene_text.strip().splitlines() if g.strip()]
        if not target_genes:
            return [], None, "⚠️ No valid genes provided."

        try:
            df = pd.read_json(links_json, orient="split")
        except Exception as e:
            return [], None, f"⚠️ Error reading GRN data: {e}"

        required_cols = {"source", "target", "-logp", "p", "coef_mean", "coef_abs"}
        if not required_cols.issubset(df.columns):
            return [], None, f"⚠️ GRN data must include {', '.join(required_cols)}."

        # Step 1: Filter for target genes
        filtered = df[df["target"].isin(target_genes)]
        if filtered.empty:
            return [], None, "⚠️ No interactions found for these genes."

        # Step 1.5: Remove zero-weight edges
        filtered = filtered[filtered["coef_abs"] > 1e-6]

        # Step 2: Keep top 2000
        filtered = filtered.sort_values(by="-logp", ascending=False).head(2000)

        # Step 3: Rank TFs
        tf_stats = (
            filtered.groupby("source")
            .agg(
                edge_count=("target", "count"),
                max_logp=("-logp", "max"),
                min_p=("p", "min"),
                median_coef=("coef_mean", "median")
            )
            .reset_index()
            .rename(columns={"source": "TF"})
        )

        tf_stats["edge_pct"] = 100 * tf_stats["edge_count"] / len(filtered)
        tf_stats["edge_pct"] = tf_stats["edge_pct"].round(1)

        tf_stats = tf_stats.sort_values(by="edge_pct", ascending=False)

        return (
            tf_stats.to_dict("records"),
            filtered.to_json(date_format="iso", orient="split"),
            f"✅ Found {len(filtered)} edges. Select up to 5 TFs to visualize."
        )



    @app.callback(
    Output("grn-network-cytoscape", "elements"),
    Output("cyto-hint-msg", "children", allow_duplicate=True),
    Input("grn-tf-table", "derived_virtual_selected_rows"),
    State("grn-tf-table", "data"),
    State("filtered-edges-store", "data"),
    prevent_initial_call=True
)
    def build_network_from_tfs(selected_rows, tf_table_data, filtered_edges_json):
        if not selected_rows or not tf_table_data or not filtered_edges_json:
            return [], "ℹ️ Select up to 3 TFs to display their network."

        selected_tfs = [tf_table_data[i]["TF"] for i in selected_rows][:5]

        df = pd.read_json(filtered_edges_json, orient="split")
        sub_df = df[df["source"].isin(selected_tfs)]

        # Limit to 10 edges per TF by -logp
        top_edges = (
            sub_df.sort_values(by=["source", "-logp"], ascending=[True, False])
            .groupby("source")
            .head(10)
        )

        nodes_set = set(top_edges["source"]).union(top_edges["target"])

        node_elements = []
        for node in nodes_set:
            is_tf = node in selected_tfs
            node_elements.append({
                "data": {"id": node, "label": node},
                "classes": "tf-node" if is_tf else "target-node"
            })

        edge_elements = []
        for _, row in top_edges.iterrows():
            edge_elements.append({
                "data": {
                    "source": row["source"],
                    "target": row["target"],
                    "title": f"{row['source']} → {row['target']}\n-logp: {row['-logp']:.2f}\np: {row['p']:.2e}\ncoef: {row['coef_mean']:.2e}",
                    "weight": f"{row['-logp']:.2f}",
                    "coef_mean": f"{row['coef_mean']:.2e}",
                    "coef_abs": f"{row['coef_abs']:.2e}",
                    "pval": f"{row['p']:.2e}",
                    "-logp": f"{row['-logp']:.2f}"
                }
            })

        return node_elements + edge_elements, f"✅ Showing {len(top_edges)} edges for {len(selected_tfs)} TF(s)."



    @app.callback(
    Output("edge-hover-info", "children"),
    Input("grn-network-cytoscape", "mouseoverEdgeData"),
    prevent_initial_call=True
    )
    def show_edge_hover(data):
        if not data:
            return ""

        return html.Div([
            html.H6(f"{data['source']} → {data['target']}", style={"marginBottom": "8px"}),
            html.Ul([
                html.Li(f"-logp: {data.get('-logp', 'N/A')}"),
                html.Li(f"p-value: {data.get('pval', 'N/A')}"),
                html.Li(f"coef: {data.get('coef_mean', 'N/A')}"),
                html.Li(f"|coef|: {data.get('coef_abs', 'N/A')}"),
            ], style={"paddingLeft": "18px", "margin": 0})
        ])


    
    
    @app.callback(
    Output("grn-links-store", "data"),  # Store full df
    Input("grn-summary-store", "data"),
    Input("grn-cluster-select", "value"),
    prevent_initial_call=True
    )
    def store_full_links(summary, selected_cluster):
        if not summary or not selected_cluster:
            raise PreventUpdate

        selected_cluster = str(selected_cluster)
        links = co.load_hdf5(summary["filtered_links_path"])
        df = links.links_dict[selected_cluster].reset_index()

        return df.to_json(orient="split")


    @app.callback(
    Output("download-filtered-table", "data"),
    Input("export-filtered-button", "n_clicks"),
    State("grn-links-store", "data"),
    prevent_initial_call=True
    )
    def export_filtered_links(n_clicks, df_json):
        if not df_json:
            raise PreventUpdate

        df = pd.read_json(df_json, orient="split")
        return dcc.send_string(df.to_csv(index=False), filename="filtered_links.csv")



    

    @app.callback(
    Output("grn-tf-input", "options"),
    Output("grn-tf-input", "value"),
    Input("grn-merged-score", "data"),       # <- now using merged centrality data
    Input("grn-centrality-metric", "value"), # centrality metric
    prevent_initial_call=True
    )
    def populate_tf_dropdown(merged_json, metric):
        if not merged_json or not metric:
            raise PreventUpdate

        try:
            df = pd.read_json(merged_json, orient="split")
        except Exception as e:
            print(f"Error loading merged_df: {e}")
            raise PreventUpdate

        if "tf" not in df.columns:
            print("Missing 'tf' column in merged_df.")
            raise PreventUpdate

        tf_list = sorted(df["tf"].dropna().unique())
        options = [{"label": tf, "value": tf} for tf in tf_list]

        # Suggest top TF by selected centrality
        suggested_tf = None
        if metric in df.columns:
            top_row = df.sort_values(by=metric, ascending=False).dropna(subset=[metric])
            if not top_row.empty:
                suggested_tf = top_row.iloc[0]["tf"]

        return options, suggested_tf or (tf_list[0] if tf_list else None)






    app.clientside_callback(
        ClientsideFunction(
            namespace="client",
            function_name="download_cytoscape_png"
        ),
        Output("dummy-output", "children"),
        Input("download-png-button", "n_clicks"),
        prevent_initial_call=True
    )




    @app.callback(
    Output("grn-cluster-select", "value",allow_duplicate=True),
    Output("grn-centrality-metric", "value",allow_duplicate=True),
    Output("grn-tf-input", "value",allow_duplicate=True),
    Input("grn-merged-score", "data"),
    State("grn-summary-store", "data"),
    State("grn-merged-score", "data"),
    prevent_initial_call=True
    )
    def set_initial_cluster_and_tf(is_complete, summary, merged_json):
        # print('this is used -----')
        if not is_complete or not summary or not merged_json:
            raise dash.exceptions.PreventUpdate

        try:
            df = pd.read_json(merged_json, orient="split")
        except Exception as e:
            print(f"❌ Error reading merged score: {e}")
            raise dash.exceptions.PreventUpdate

        # Default to first available cluster
        default_cluster = summary["clusters"][0] if summary["clusters"] else None

        # Default to top TF by centrality
        default_metric = "degree_centrality_all"
        top_tf = None
        if default_metric in df.columns:
            df_cluster = df[df["cluster"] == default_cluster]
            top = df_cluster.sort_values(by=default_metric, ascending=False).dropna(subset=[default_metric])
            if not top.empty:
                top_tf = top.iloc[0]["tf"]

        return default_cluster, default_metric, top_tf



    @app.callback(
    Output("grn-centrality-plot", "figure"),
    Output("grn-top-tfs-plot", "figure"),
    Input("grn-merged-score", "data"),
    Input("grn-cluster-select", "value"),
    Input("grn-centrality-metric", "value"),
    Input("grn-top-n-tfs", "value"),
    Input("grn-tf-input", "value"),
    prevent_initial_call=True
    )
    def update_grn_plots(merged_json, selected_cluster, metric, top_n, tf_name):


        if isinstance(selected_cluster, dict):
            selected_cluster = selected_cluster.get("value")

        if not merged_json or not selected_cluster or not metric:
            return go.Figure(), go.Figure()

        df = pd.read_json(merged_json, orient="split")
        # print("📦 fig1  fig2:", df.head())
        # Top TFs (fig2)
        fig2 = go.Figure()
        try:
            fig2 = plot_cluster_centrality_from_df(
                df,
                cluster_id=selected_cluster,
                centrality_metric=metric,
                top_n=top_n
            )
        except Exception as e:
            fig2 = go.Figure(layout={"title": f"Top TF error: {e}"})

        # TF plot (fig1)
        fig1 = go.Figure()
        if tf_name:
            try:
                fig1 = plot_tf_centrality_across_clusters(df, tf_name=tf_name)
            except Exception as e:
                fig1 = go.Figure(layout={"title": f"TF plot error: {e}"})

        return fig2, fig1