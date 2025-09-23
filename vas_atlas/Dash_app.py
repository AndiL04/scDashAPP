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
adata = sc.read_h5ad("../merged_normalized_sub.h5ad")

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

# Function to check if a port is in use
def is_port_in_use(port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) == 0

cb_process = None
if not is_port_in_use(8888):
    cb_process = subprocess.Popen(["cbBuild", "-o", "./public_html/", "-p", "8888"])
    atexit.register(cb_process.terminate)
else:
    print("✅ UCSC Cell Browser is already running at http://127.0.0.1:8888/")
    
# Dash app setup
app = Dash(__name__)
app.title = "Gene Expression Explorer"

app.layout = html.Div([
    # Header with logo and title
    html.Div([
        html.Img(src="/assets/logo.png", style={"height": "60px","marginRight": "25px"}),
        html.H1("Human Brain Vasculature Single-Cell Explorer", style={
            "fontFamily": "Alfabet",
            'color': '#31175A',
            "fontWeight": "bold",
            "fontSize": "28px",
            "textAlign": "center",
            "margin": "0",
        })
    ], style={"display": "flex", "alignItems": "center", "paddingLeft": "20px","marginBottom": "20px"}),

    # Tabs section
    dcc.Tabs([
        dcc.Tab(label='Gene Expression',
                style={'fontFamily': 'Arial', 'fontSize': '15px'},
                selected_style={'fontFamily': 'Alfabet', 'fontWeight': 'bold','color': '#31175A','borderTop': '2px solid #31175A',},
                children=[
            html.Div([
                # Left control panel
                html.Div([
                    html.Label("Select a gene:", style={'fontWeight': 'bold',"fontFamily": "Alfabet","fontSize": "15px"}),
                    dcc.Dropdown(
                        id='gene-dropdown',
                        options=[{'label': gene, 'value': gene} for gene in all_genes],
                        value='CLDN5',
                        searchable=True,
                        style={'width': '250px',"fontFamily": "Alfabet","fontSize": "14px"}
                    ),
                    html.Label("Group by:", style={'fontWeight': 'bold','marginTop': '30px',"fontFamily": "Alfabet","fontSize": "15px"}),
                    dcc.RadioItems(
                        id='groupby-toggle',
                        options=[
                            {'label': 'Cell Class', 'value': 'cell.class'},
                            {'label': 'Cell Type', 'value': 'cell.type.middle'},
                            {'label': 'Sub Cell Type', 'value': 'cell.type'}
                        ],
                        value='cell.class',
                        labelStyle={'display': 'block','marginTop': '5px',"fontFamily": "Alfabet","fontSize": "14px"},
                        style={"fontFamily": "Alfabet"})
                ], style={'width': '250px','paddingRight': '30px',"fontFamily": "Alfabet"}),
                # Right panel with plot
                html.Div([
                    dcc.Graph(id='expression-barplot',config={'responsive': True},style={"fontFamily": "Alfabet"}),
                    dcc.Graph(id='dot-plot', config={'responsive': True}, style={"fontFamily": "Alfabet", "marginTop": "20px"})
                ], style={'flexGrow': 1, 'minWidth': '0'})
            ], style={'display': 'flex', 'flexDirection': 'row'})
        ]),

        dcc.Tab(label='Cell Browser',
                style={'fontFamily': 'Alfabet', 'fontSize': '15px'},
                selected_style={'fontFamily': 'Alfabet', 'fontWeight': 'bold','color': '#31175A','borderTop': '2px solid #31175A'},
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
    Output('dot-plot', 'figure'),
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

    df_grouped = df.groupby('group').agg(expression_mean=('expression', 'mean'),expression_sem=('expression', 'sem')).reset_index()
    df_grouped.rename(columns={'expression_mean': 'expression', 'expression_sem': 'error'}, inplace=True)

    # set colors based on groupby_col
    color_map = cell_class_colors if groupby_col == 'cell.class' else {}

    bar_fig = px.bar(
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

    bar_fig.update_traces(marker_line_width=1, marker_line_color='black')
    bar_fig.update_xaxes(tickangle=45, tickfont=dict(size=12))
    bar_fig.update_yaxes(ticklabelposition='outside', tickfont=dict(size=12), ticks="outside")
    bar_fig.update_layout(
        autosize=True,paper_bgcolor='white', plot_bgcolor='white',
        margin=dict(l=70, r=40, t=60, b=120),height=600,
        legend=dict(orientation="h", yanchor="bottom",y=-0.75,xanchor="center",x=0.5),
        shapes=[{
            'type': 'rect',
            'xref': 'paper', 'yref': 'paper',
            'x0': 0.00, 'y0': -0.01, 'x1': 1.005, 'y1': 1.005,
            'line': {'color': 'black', 'width': 1},
            'layer': 'below'
        }]
    )

    ## Dot plot
    obs = adata.obs.copy()
    df_dot = pd.DataFrame({
        'region': obs['region_name'],
        'group': obs[groupby_col],
        'expression': gene_expr
    })
    df_dot['is_expressed'] = df_dot['expression'] > 0
    grouped_dot = df_dot.groupby(['region', 'group']).agg(
        avg_expr=('expression', 'mean'),
        pct_expressed=('is_expressed', 'mean')
    ).reset_index()
    grouped_dot['pct_expressed'] *= 100
    grouped_dot['pct_expressed'] = grouped_dot['pct_expressed'].fillna(0)

    dot_fig = px.scatter(
        grouped_dot, x='group', y='region', size='pct_expressed', color='avg_expr',
        color_continuous_scale='magma',
        labels={'group': groupby_col.replace('.', ' ').title(), 'region': 'Brain Region', 'avg_expr': 'Avg. Expression', 'pct_expressed': '% Expressing Cells'},
        title=f"{gene} expression across brain regions and {groupby_col.replace('.', ' ')}"
    )
    dot_fig.update_traces(marker=dict(line=dict(width=0.5, color='gray')))
    dot_fig.update_layout(
        autosize=True, paper_bgcolor='white', #plot_bgcolor='white',
        margin=dict(l=70, r=40, t=60, b=80), height=950,
        coloraxis_colorbar=dict(title='Avg Expr')
    )    

    return bar_fig, dot_fig


if __name__ == '__main__':
    # Optional: open the app in a web browser automatically
    webbrowser.open("http://127.0.0.1:8050/")
    app.run_server(debug=True)
