#!/usr/bin/env python3

import sys
import os
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
#import muon as mu
TITLE = 'scRepo'

if len(sys.argv) < 3:
    print("Usage: screpo mode database.csv port")
    print("mode: local or remote")
    print("local: run in your local machine")
    print("remote: run in remote machine")
    exit(0)


# Styling
#dbc_css = "https://cdn.jsdelivr.net/gh/AnnMarieW/dash-bootstrap-templates/dbc.min.css"
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
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
    "margin-left": "2rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}

nav = dbc.Nav(
    [
        dbc.NavLink("scRepo", active=True, href="#"),
        dbc.NavLink("Usage", disabled=True, href="#"),
        dbc.NavLink("Contact", disabled=True,  href="#")
    ],
    pills=True
)


navbar = dbc.NavbarSimple(
    children=[
        dbc.NavItem(dbc.NavLink("Page 1", href="#")),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("More pages", header=True),
                dbc.DropdownMenuItem("Page 2", href="#"),
                dbc.DropdownMenuItem("Page 3", href="#"),
            ],
            nav=True,
            in_navbar=True,
            label="More",
        ),
    ],
    brand="NavbarSimple",
    brand_href="#",
    color="primary",
    dark=True,
)


def get_data_table(df, page_size = 5):

    """
    Reads in a pandas dataframd and outputs a dash_table.Datatable obj
    """

    df.index = list(df.DataId)
    dt = dash_table.DataTable(
        id='project_dt',
        data=df.to_dict('records'),
        #columns=[{"name": i, "id": i, 'disable_sort': True } if i == "File"  else {"name": i, "id": i } for i in df.columns ],
        columns=[ {"name": i, "id": i, "hideable": True} for i in df.columns],
        style_cell={
            #'overflow': 'hidden',
            #'textOverflow': 'ellipsis',
            #'maxWidth': 4,
            'font-family': 'sans-serif'
        },
        style_header={
            'backgroundColor': 'rgb(156, 15, 15)',
            'color': 'white'
        },
        style_data={

            'whiteSpace': 'normal',
            'height': 'auto',
        },
        style_as_list_view=False,  # No vertical lines
        filter_action='native',
        row_selectable='single',
        page_size= page_size,
        hidden_columns = ['File']
    )
    return dt
#adatafile = 'ref_landscape.h5ad'
def read_adata(adatafile):
    adata = sc.read_h5ad(adatafile)
    return adata

def get_umap(adata):
    return pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index = adata.obs_names)

def blank_figure():
    fig = go.Figure(go.Scatter(x=[], y=[]))
    fig.update_layout(template=None)
    fig.update_xaxes(showgrid=False, showticklabels=False, zeroline=False)
    fig.update_yaxes(showgrid=False, showticklabels=False, zeroline=False)

    return fig


def populate_sidebar(adata, filename):
    """
    Returns the sidebar containing drop down menus in the left bar
    :param adata: anndata obj
    :return: html.Div()
    """
    cols = adata.obs.columns
    vars = adata.var_names
    #Dim reduction obsm
    dim_reds = list(adata.obsm.keys())
    dim_reds = [ i.replace("X_", "") for i in dim_reds ]
    if "umap" in dim_reds:
        dim_reds_default = "umap"
    else:
        dim_reds_default = dim_reds[0]

    obs_dropdown = dcc.Dropdown(cols, id='obs_names', value=cols[len(cols)-1])
    var_dropdown = dcc.Dropdown(vars, id='var_names')
    dims_dropdown = dcc.Dropdown(dim_reds, id='dim_reds', value=dim_reds_default)
    point_sizes = dcc.Dropdown(list(range(20)), id='point_size', value=3)

    # Download link
    dl = "rsync -avzP user@login-cpu.hpc.cam.ac.uk:" + os.path.abspath(filename) + " ./"

    download_link = dcc.ConfirmDialogProvider(
        children=html.Button('Download'),
        id='wget',
        message="Download data\n " + dl
    )

    send_div = [
        html.H4('Reduction maps', style={'color': 'white'}),
        dims_dropdown,
        html.Br(),
        html.H4('Metadata', style={'color': 'white'}),
        obs_dropdown,
        html.Br(),
        html.H4('Genes', style={'color': 'white'}),
        var_dropdown,
        html.Br(),
        html.H4('Point size', style={'color': 'white'}),
        point_sizes,
        html.Br(),
        html.H4('Download File', style={'color': 'white'}),
        download_link,
        #html.Button("Download File2", id="btn_file"),
        #dcc.Download(id="download-image")
    ]
    return send_div

def get_dimred_pd(adata, dimred='umap'):
    obsm_col = "X_" + dimred
    print(obsm_col)
    colnames = [ dimred+i for i in ["1","2"] ]
    matx = adata.obsm[obsm_col]
    print(matx.shape)
    matx = matx[:,[0,1]]
    df  = pd.DataFrame(
        matx,
        columns = colnames,
        index = adata.obs_names
    )
    return df

def dimred_plot(colname, ps, dimred="umap", obs_or_var = 'obs' ):
    umap = get_dimred_pd(adata, dimred = dimred)
    colnames = list(umap.columns)
    if obs_or_var == 'obs':
        umap = pd.merge(umap, adata.obs, left_index=True, right_index=True)
        fig = px.scatter(data_frame=umap,
                         x=colnames[0],
                         y=colnames[1],
                         color = colname,
                         width=1000, height=800)
    if obs_or_var == 'var':
        if colname == None:
            return None
        if isinstance(adata[:, colname].X, anndata._core.views.SparseCSRView):
            x = adata[:, colname].X.toarray()
        else:
            x = adata[:, colname].X.toarray()
        t = pd.DataFrame(x, index=adata.obs.index, columns=[colname])
        umap = pd.merge(umap, t, left_index=True, right_index=True)
        print(umap.head())
        fig = px.scatter(data_frame=umap,
                         x=colnames[0],
                         y=colnames[1],
                         color=colname,
                         width=1000, height=800,
                         color_continuous_scale='brwnyl')
    fig.update_traces(marker_size=ps)
    fig.update_layout(legend={'itemsizing': 'constant'})
    graph = dcc.Graph(figure=fig)
    return graph

# Start Application
#app = dash.Dash(external_stylesheets=[dbc.themes.PULSE, dbc_css])
app = dash.Dash(
    external_stylesheets=[dbc.themes.ZEPHYR],
    suppress_callback_exceptions=True
)
app.title = 'scRepo'
load_figure_template('LUX')


# Read data table
if len(sys.argv) == 4:
    csvfile = 'data.csv'
    csvfile = sys.argv[2]
    print("Reading from ", csvfile)
    df = pd.read_csv(csvfile)
    dt = get_data_table(df, page_size=10)

# Application Layout
# Side bar layout details and vars
sidebar_layout = html.Div([
    html.H1("scRepo", style = {'color': 'white'}),
    html.Hr(
        style={
            "borderWidth": "0.2vh",
            "width": "100%",
            "borderColor": "white",
            "opacity": "unset",
        }
    ),

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
        #html.H1(children=TITLE),
        html.Br(),
        html.Div(dt , className = "dbc dbc-row-selectable", style={'justify': 'left'})
    ]),
    html.Br(),
    html.Div( #--> UMAP
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

import dash_bootstrap_components as dbc
app.layout = dbc.Container(
   [
    dbc.Row([sidebar_layout, content_layout])
   ]
)

#temp = dict(adata = None)

#adata = None
# Trigger when data table is selected
@app.callback(
    Output('umap_options', 'children'),
    Input('project_dt', 'selected_rows'),
    prevent_initial_call=True,
    suppress_callback_exceptions=True
)

# From a unique row selected, read adata and also populate the sidebar
def get_adata(selected_rows):
    """
    From table, send in selected_rows to get adata
    This also populates the sidebar
    :param selected_rows: Returns a row data from 'dt'
    :return: html.Div (left drop down menus)

    """
    print(selected_rows)
    selected_adata = df.iloc[selected_rows].File
    expt = df.iloc[selected_rows].ExptName
    global dataid # A potential bug here. The dash servers may be sharing instances in which case adata will be overwritten by other instances. Rare but possible.
    global adata

    dataid = df.iloc[selected_rows].DataId
    filename = df.iloc[selected_rows].File.values[0]
    modality = df.iloc[selected_rows].Modality.values[0]
    #global adata # A potential bug here. The dash servers may be sharing instances in which case adata will be overwritten by other instances. Rare but possible.
    if modality == 'Unimodal':
        print("Loading unimodal anndata")
        adata = sc.read_h5ad(filename)

    elif modality == 'Multimodal':
        """
        Have not used multimodal because muon is ridculously slow in loading data. 
        Instead I have used a hack to add multimodal in the same adata.X with varying vars. 
        The only downside is that this requires the same set of cells across all modalities.
        """
        print("Loading multimodal mudata")
        #adata = mu.read_h5mu(filename)
        adata = sc.read_h5ad(filename) # For now just read multimodal as unimodal files
        #temp['adata'] = adata
        print("Successful")
    print("adata read successfully")

    #print(adata.obs.columns)
    return populate_sidebar(adata, filename)


# Plot gene expression values
@app.callback(
    Output('umap', 'children', allow_duplicate=True),
    Input('var_names', 'value'),
    Input('point_size', 'value'),
    Input('dim_reds', 'value'),
    prevent_initial_call=True,
    suppress_callback_exceptions=True
)
def plot_var(varname, ps, dimred):
    g = dimred_plot(varname, ps, dimred, obs_or_var='var')
    return g

# Plots obs columns
@app.callback(
    Output('umap', 'children', allow_duplicate=True),
    Input('obs_names', 'value'),
    Input('point_size', 'value'),
    Input('dim_reds', 'value'),
    prevent_initial_call=True,
)

def plot_obs(obs_col, ps, dimred):
    g = dimred_plot(obs_col, ps, dimred, obs_or_var='obs')
    return g

"""
@app.callback(
    Output("download-image", "data"),
    Input("btn_file", "n_clicks"),
    Input('project_dt', 'selected_rows'),
    prevent_initial_call=True
)
def get_download(n_clicks, selected_rows):
    print(selected_rows)
    selected_adata = df.iloc[selected_rows].File
    filename = df.iloc[selected_rows].File.values[0]
    return dcc.send_file(
        filename
    )


"""

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: screpo mode database.csv port")
        print("mode: local or remote")
        print("local: run in your local machine")
        print("remote: run in remote machine")
        exit(0)
    port=sys.argv[3]
    runmode = sys.argv[1]


    if runmode == 'local':
        app.run_server(debug=True, port = port) # For internal testing on local server
    if rumode == 'remote':
        app.run_server(debug=True, port = port, host = '0.0.0.0') # For external deployment
