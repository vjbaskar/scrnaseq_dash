## 1st Example Application: Basic
import dash
#import dash_core_components as dcc
#import dash_html_components as html
import scanpy as sc
import pandas as pd
import anndata
from dash import Dash, html, dash_table, dcc
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template
import plotly.express as px
import plotly.graph_objects as go


# Styling

SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
}

CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}


# Read data table
df = pd.read_csv('data.csv')
df.index = list(df.DataId)

dt = dash_table.DataTable(
    id='project_dt',
    data=df.to_dict('records'),
    row_selectable='single',
    columns=[{"name": i, "id": i} for i in df.columns],
    style_as_list_view=True,  # No vertical lines
    style_cell={'padding': '5px'},
    style_header={
        'backgroundColor': 'Orange',
        'fontWeight': 'bold'
    }
)



#adatafile = 'ref_landscape.h5ad'
def read_adata(adatafile):
    adata = sc.read_h5ad(adatafile)
    return adata

def get_umap(adata):
    return pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index = adata.obs_names)

#adata = read_adata(adatafile)

# Start Application
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])
# Define application window name/ title
app.title = 'scrnaseq'
load_figure_template('BOOTSTRAP')



#umap_base = dcc.Graph(id = 'umap',
#                      figure = px.scatter(data_frame=umap, x = 'UMAP1', y = 'UMAP2', width=1000, height=800)
#                      )
#dt = dash_table.DataTable(data=adata.obs.head(10).to_dict('records'), page_size=10)

def blank_figure():
    fig = go.Figure(go.Scatter(x=[], y=[]))
    fig.update_layout(template=None)
    fig.update_xaxes(showgrid=False, showticklabels=False, zeroline=False)
    fig.update_yaxes(showgrid=False, showticklabels=False, zeroline=False)

    return fig


# Application Layout
#  Drop down menu


# sidebar_layout = html.Div(
#     [
#         html.H4('Metadata'),
#         obs_dropdown,
#         html.Br(),
#         html.H4('Genes'),
#         var_dropdown,
#         html.Br(),
#         html.H4('Point size'),
#         point_sizes
#     ],
#     style=SIDEBAR_STYLE
# )

sidebar_layout = html.Div([
    html.Div(
        id = 'umap_options',
    ),
    html.Div(
        [
            dcc.Loading(
                id="loading_adata",
                type="default",
                children=html.Div(id="umap_options"),
                fullscreen=True
            )
        ]
    ),
],
style = SIDEBAR_STYLE
)


content_layout = html.Div([
    html.Div([
        html.H1(children='UMAP'),
        html.Br(),
        html.Div(dt)
    ]),
    html.Br(),
    #html.Button('Run', id='expt_id', n_clicks=0),
    # html.Div(
    #     [
    #         html.P(id = 'adata'),
    #         dcc.Loading(
    #             id="loading_adata",
    #             type="default",
    #             children=html.Div(id="adata")
    #         )
    #     ]
    # ),
    html.Div(
        [
            html.Div([
                dcc.Graph(id='umap', figure=blank_figure())
            ], id = 'centre'),
            html.Div(
                [
                    dcc.Loading(
                        id="loading_adata1",
                        type="default",
                        children=html.Div(id="centre"),
                        fullscreen=False
                    )
                ]
            )

        ]
    )
    #html.Div(id = 'umap')
],
    style=CONTENT_STYLE
)

app.layout = html.Div([sidebar_layout, content_layout])

# Trigger when data table is selected
@app.callback(
    Output('umap_options', 'children'),
    Input('project_dt', 'selected_rows'),
    prevent_initial_call=True
)
def get_adata(selected_rows):
    print(selected_rows)
    selected_adata = df.iloc[selected_rows].File
    expt = df.iloc[selected_rows].ExptName
    global dataid
    dataid = df.iloc[selected_rows].DataId
    filename = df.iloc[selected_rows].File.values[0]
    global adata
    print(filename)
    adata = sc.read_h5ad(filename)
    print("adata read successfully")
    print(adata.obs.columns)
    obs_dropdown = dcc.Dropdown(adata.obs.columns, id='obs_names')
    var_dropdown = dcc.Dropdown(adata.var_names, id='var_names')
    point_sizes = dcc.Dropdown(list(range(20)), id='point_size', value=3)

    send_div = [
            html.H4('Metadata'),
            obs_dropdown,
            html.Br(),
            html.H4('Genes'),
            var_dropdown,
            html.Br(),
            html.H4('Point size'),
            point_sizes

    ]

    return send_div

# Plots obs columns
@app.callback(
    Output('umap', 'figure', allow_duplicate=True),
    Input('obs_names', 'value'),
    Input('point_size', 'value'),
    prevent_initial_call=True,
)
def plot_umap(value, ps):
    umap = get_umap(adata)
    umap = pd.merge(umap, adata.obs, left_index=True, right_index=True)
    fig = px.scatter(data_frame=umap, x='UMAP1', y='UMAP2', width=1000, height=800, color = value)
    fig.update_traces(marker_size=ps)
    return fig

# Plot gene expression
@app.callback(
    Output('umap', 'figure' ),
    Input('var_names', 'value'),
    Input('point_size', 'value'),
    prevent_initial_call=True
)
def plot_umap_gene_exp(value, ps):
    umap = get_umap(adata)
    if isinstance(adata[:,value].X, anndata._core.views.SparseCSRView):
        x = adata[:, value].X.toarray()
    else:
        x = adata[:, value].X.toarray()

    t = pd.DataFrame(x, index=adata.obs.index, columns=[value])

    umap = pd.merge(umap, t, left_index=True, right_index=True)
    fig = px.scatter(data_frame=umap, x='UMAP1', y='UMAP2', width=1000, height=800, color = value)
    fig.update_traces(marker_size=ps)
    return fig


#if __name__ == '__main__':
app.run_server(debug=True)
