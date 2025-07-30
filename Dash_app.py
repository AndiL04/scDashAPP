import scanpy as sc
import pandas as pd
import numpy as np
import plotly.express as px
from dash import Dash, dcc, html, Input, Output
import subprocess
import atexit
import webbrowser
import socket

# Load the h5ad file
adata = sc.read_h5ad("/home/aliu/projects/vascular/cxg/human_brain_vasculature_atlas_normalized.h5ad")

# Precompute: list of all genes
all_genes = adata.var_names.tolist()

# Define colors for cell.class
cell_class_colors = {
    "Astrocyte": "#F06719",
    "Ependymal_Epithelial": "#23767C",
    "Oligodendrocytes": "#00BFC4",
    "Neuron": "#08415C",
    "Microglia_Macrophage_T": "#DC143C",
    "Mural_Cell": "#A26DC2",
    "Endothelial": "#fcbe05",
    "OPC": "#0072B2",
    "Fibroblast": "#5b844d",
}

def is_port_in_use(port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) == 0

cb_process = None
if not is_port_in_use(8888):
    cb_process = subprocess.Popen(["cbBuild", "-o", "./public_html/", "-p", "8888"])
    atexit.register(cb_process.terminate)
else:
    print("✅ UCSC Cell Browser is already running at http://127.0.0.1:8888/")
    
# # Ensure it terminates when Dash app stops
# atexit.register(cb_process.terminate)

# Dash app setup
app = Dash(__name__)
app.title = "Gene Expression Explorer"

app.layout = html.Div([

    # Header with logo and title
    html.Div([
        html.Img(
            src="/assets/logo.png",
            style={
                "height": "60px",
                "marginRight": "25px"
            }
        ),
        html.H1("Human Brain Vasculature Single-Cell Explorer", style={
            "fontFamily": "Alfabet",
            'color': '#31175A',
            "fontWeight": "bold",
            "fontSize": "28px",
            "textAlign": "center",
            "margin": "0",
        })
    ], style={
        "display": "flex",
        "alignItems": "center",
        "paddingLeft": "20px",
        "marginBottom": "20px"
    }),

    # Tabs section
    dcc.Tabs([
        dcc.Tab(label='Gene Expression',
                style={'fontFamily': 'Arial', 'fontSize': '15px'},
                selected_style={
                    'fontFamily': 'Alfabet', 
                    'fontWeight': 'bold',
                    'color': '#31175A',
                    'borderTop': '2px solid #31175A',
                },
                children=[
            html.Div([
                # Left control panel
                html.Div([
                    html.Label("Select a gene:", style={
                        'fontWeight': 'bold',
                        "fontFamily": "Alfabet",
                        "fontSize": "15px"
                    }),
                    dcc.Dropdown(
                        id='gene-dropdown',
                        options=[{'label': gene, 'value': gene} for gene in all_genes],
                        value='CLDN5',
                        searchable=True,
                        style={
                            'width': '250px',
                            "fontFamily": "Alfabet",
                            "fontSize": "14px"
                        }
                    ),

                    html.Label("Group by:", style={
                        'fontWeight': 'bold',
                        'marginTop': '30px',
                        "fontFamily": "Alfabet",
                        "fontSize": "15px"
                    }),
                    dcc.RadioItems(
                        id='groupby-toggle',
                        options=[
                            {'label': 'Cell Class', 'value': 'cell.class'},
                            {'label': 'Cell Type Middle', 'value': 'cell.type.middle'},
                            {'label': 'Cell Type', 'value': 'cell.type'}
                        ],
                        value='cell.class',
                        labelStyle={
                            'display': 'block',
                            'marginTop': '5px',
                            "fontFamily": "Alfabet",
                            "fontSize": "14px"
                        },
                        style={"fontFamily": "Alfabet"}
                    )
                ], style={
                    'width': '250px',
                    'paddingRight': '30px',
                    "fontFamily": "Alfabet"
                }),

                # Right panel with plot
                html.Div([
                    dcc.Graph(
                        id='expression-barplot',
                        config={'responsive': True},
                        style={"fontFamily": "Alfabet"}
                    )
                ], style={'flexGrow': 1, 'minWidth': '0'})
            ], style={'display': 'flex', 'flexDirection': 'row'})
        ]),

        dcc.Tab(label='Cell Browser',
                style={'fontFamily': 'Alfabet', 'fontSize': '15px'},
                selected_style={
                    'fontFamily': 'Alfabet', 
                    'fontWeight': 'bold',
                    'color': '#31175A',
                    'borderTop': '2px solid #31175A',
                },
                children=[
            html.Div([
                html.Iframe(
                    src="http://127.0.0.1:8888/",
                    style={'width': '100%', 'height': '85vh', 'border': 'none'}
                )
            ], style={"fontFamily": "Alfabet"})
        ])
    ])
], style={"fontFamily": "Alfabet"})  # Apply Alfabet globally to the whole page


@app.callback(
    Output('expression-barplot', 'figure'),
    Input('gene-dropdown', 'value'),
    Input('groupby-toggle', 'value')
)
def update_figures(gene, groupby_col):
    if gene not in adata.var_names:
        return px.bar(title="Gene not found")

    x = adata[:, gene].X
    gene_expr = x.toarray().flatten() if hasattr(x, "toarray") else x

    df = pd.DataFrame({
        'group': adata.obs[groupby_col],
        'expression': gene_expr
    })

    df_grouped = df.groupby('group').agg(
        expression_mean=('expression', 'mean'),
        expression_sem=('expression', 'sem')
    ).reset_index()
    df_grouped.rename(columns={'expression_mean': 'expression', 'expression_sem': 'error'}, inplace=True)

    if groupby_col == 'cell.class':
        df_grouped["color"] = df_grouped["group"].map(cell_class_colors)
        color_map = cell_class_colors
    else:
        color_map = {}

    fig = px.bar(
        df_grouped,
        x='group',
        y='expression',
        color='group',
        error_y='error',
        color_discrete_map=color_map,
        title=f"Average expression of {gene}",
        labels={'expression': 'Avg. Expression', 'group': groupby_col.replace('.', ' ').title()},
        hover_data={'expression': ':.2f', 'error': ':.2f'}
    )

    fig.update_traces(marker_line_width=1, marker_line_color='black')
    fig.update_xaxes(tickangle=45, tickfont=dict(size=12))
    fig.update_yaxes(ticklabelposition='outside', tickfont=dict(size=12), ticks="outside")
    fig.update_layout(
        autosize=True,
        paper_bgcolor='white',
        plot_bgcolor='white',
        margin=dict(l=70, r=40, t=60, b=120),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=-2,
            xanchor="center",
            x=0.5
        ),
        shapes=[{
            'type': 'rect',
            'xref': 'paper', 'yref': 'paper',
            'x0': 0.00, 'y0': -0.01, 'x1': 1.005, 'y1': 1.005,
            'line': {'color': 'black', 'width': 1},
            'layer': 'below'
        }]
    )

    return fig


if __name__ == '__main__':
    # Optional: open the app in a web browser automatically
    webbrowser.open("http://127.0.0.1:8050/")
    app.run_server(debug=True)
