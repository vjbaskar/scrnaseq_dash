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
    "width": "20rem",
    "padding": "2rem 1rem",
    "background-color": "#9C0F0F",
}

CONTENT_STYLE = {
    "margin-left": "16rem",
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
    style_as_list_view=False,  # No vertical lines
    filter_action='native',
    style_cell={'padding': '5px'},
    style_header={
        'backgroundColor': 'rgb(156, 15, 15)',
        'color': 'white'
    },
    style_data={
        'backgroundColor': 'rgb(50, 50, 50)',
        'color': 'white'
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
app = dash.Dash(external_stylesheets=[dbc.themes.LUX])
# Define application window name/ title
app.title = 'scrnaseq'
load_figure_template('LUX')



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
                parent_className='loading_wrapper',
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
            html.Div(id = 'umap'),
            html.Div(
                [
                    dcc.Loading(
                        id="loading_adata1",
                        type="default",
                        children=html.Div(id="centre"),
                        fullscreen=False,
                        className='loading_wrapper'
                    )
                ],

            )

        ]
    )
    #html.Div(id = 'umap')
],
    style=CONTENT_STYLE
)

app.layout = dbc.Container(html.Div([sidebar_layout, content_layout]))

def anndata_sidebar(adata):
    """
    Returns the sidebar containing drop down menus in the left bar
    :param adata: anndata obj
    :return: html.Div()
    """
    obs_dropdown = dcc.Dropdown(adata.obs.columns, id='obs_names', value=adata.obs.columns[0])
    var_dropdown = dcc.Dropdown(adata.var_names, id='var_names', value = adata.var.index[0])
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


# Trigger when data table is selected
@app.callback(
    Output('umap_options', 'children'),
    Input('project_dt', 'selected_rows'),
    prevent_initial_call=True
)
def get_adata(selected_rows):
    """
    From table, send in selected_rows to get adata
    :param selected_rows: Returns a row data from 'dt'
    :return: html.Div (left drop down menus)
    """
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
    return anndata_sidebar(adata)

# Plots obs columns
@app.callback(
    Output('umap', 'children', allow_duplicate=True),
    Input('obs_names', 'value'),
    Input('point_size', 'value'),
    prevent_initial_call=True,
)
def plot_umap(value, ps):
    umap = get_umap(adata)
    umap = pd.merge(umap, adata.obs, left_index=True, right_index=True)
    fig = px.scatter(data_frame=umap, x='UMAP1', y='UMAP2', width=1000, height=800, color = value)
    fig.update_traces(marker_size=ps)
    fig.update_layout(legend={'itemsizing': 'constant'})
    graph = dcc.Graph(figure=fig)
    return graph



# Expression sliders
def sliders(fig):
    steps = []
    for i in range(5):
        step = dict(
            method="update",
            args=[{"visible": [False] * 5},
                  {"title": "Slider switched to step: " + str(i)}],  # layout attribute
        )
        step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
        steps.append(step)

    sliders = [dict(
        active=10,
        currentvalue={"prefix": "Frequency: "},
        pad={"t": 50},
        steps=steps
    )]

    fig.update_layout(
        sliders=sliders
    )
    return fig


# Plot gene expression
@app.callback(
    Output('umap', 'children' ),
    Input('var_names', 'value'),
    Input('point_size', 'value'),
    Input('range-slider', 'value'),
    prevent_initial_call=True
)
def plot_umap_gene_exp(value, ps, range_exp):
    umap = get_umap(adata)
    if isinstance(adata[:,value].X, anndata._core.views.SparseCSRView):
        x = adata[:, value].X.toarray()
    else:
        x = adata[:, value].X.toarray()

    t = pd.DataFrame(x, index=adata.obs.index, columns=[value])
    print(range_exp)
    umap = pd.merge(umap, t, left_index=True, right_index=True)
    fig = px.scatter(data_frame=umap, x='UMAP1', y='UMAP2', width=1000, height=800, color = value)
    fig.update_traces(marker_size=ps)
    fig.update_layout(legend={'itemsizing': 'constant'})
    #fig = sliders(fig)

    rs = dcc.RangeSlider(
        id='range-slider',
        min=0, max=2.5, step=0.1,
        marks={0: '0', 2.5: '2.5'},
        value=[0.5, 2]
    )

    graph = [rs, dcc.Graph(figure=fig)]
    return graph


#if __name__ == '__main__':
app.run_server(debug=False)
