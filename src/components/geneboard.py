# Gene Exploration. In common-gene 2 gene. 
import base64
import datetime
import time
import io
import os 
import pandas as pd 
import itertools
from collections import defaultdict
import networkx as nx 

import plotly.graph_objs as go
import dash
from dash import dash_table
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash import dcc
from dash import html

from app import app
from src.components.header import headerComponent_geneboard


def get_genes_relations(dfgenes):
    
    """ From the table of articles per gene term, obtain the number of articles
    associated with a gene and the number of common articles between genes
    
    Arguments:
        dfgenes {[dataframe]} -- dfgenes: data frame result of the IdMiner run
    
    Returns:
        [dict] -- [article_by_gene: key is gene, value is the list of items
                                    related to the gene]
        [dict] -- [gene_common: key is gene, value is the dictionary list with
                                genes as keys (excluding the first) and values
                                are items the genes have in common]
        [dict] -- [gene_common_count: key is gene, value is the dictionary list
                                with genes as keys (excluding the first) and 
                                values are the number of items the genes have
                                in common]
    """
    print("****** get_genes_relations ******")

    genes = dfgenes.columns[0:-1]
    article_by_gene = {}
    for gene in genes:
        article_by_gene[gene] = list(set(itertools.chain.from_iterable(dfgenes[gene].dropna().astype(str).str.split(",").tolist())))

    gene_common = defaultdict(dict)
    gene_common_count = defaultdict(dict)

    for gene in genes:
        for gene2 in genes:
            if gene != gene2:
                in_common = list(set(article_by_gene[gene]) & set(article_by_gene[gene2]))
                if len(in_common) > 0:
                    gene_common[gene][gene2] = in_common
                    gene_common_count[gene][gene2] = len(in_common)

    return article_by_gene,gene_common,gene_common_count



def generate_dataframe(article_by_gene,gene_common_count):
    """[summary]
    
    Arguments:
        article_by_gene {[dict]} -- key is gene, value is the list of items 
                                    related to the gene
        gene_common_count {[dict]} -- key is gene, value is the dictionary with
                                      keys as genes (excluding the first) and 
                                      values as the number of items the genes
                                      have in common
    
    Returns:
        [dataframe] -- data frame showing the number of articles in common
                       between 2 genes
    """
    print("****** generate_dataframe ******")

    genes = list(article_by_gene.keys()) #Define gene names as the columns
    records = []
    start_query_value = 0
    start_query_max = 0
    start_query = ""
    for first_gene in genes: #Obtain relationships between genes
        row = []
        for second_gene in genes:
            if first_gene == second_gene:
                own_articles = len(article_by_gene[first_gene])
                row.append(own_articles)
            else:
                if second_gene in gene_common_count[first_gene].keys():
                    row.append(gene_common_count[first_gene][second_gene])
                else:
                    row.append(0)
        start_query_value = sum(row) - own_articles
        if start_query_value > start_query_max:
            start_query_max = start_query_value 
            start_query = first_gene
        row.insert(0,len(gene_common_count[first_gene]))
        row.insert(0, first_gene)            
        records.append(row)
    genes.insert(0,"In-common")    
    genes.insert(0,"Gene")
    df = pd.DataFrame(records, columns = genes)
    return df,start_query


def get_gene_edges(query_gene,gene_common_count):
    """
    Function that creates the structure of tuples (list of tuples) necessary to
    establish the edges of the graph. Each tuple consists of query_gene, 
    subject_gene, and the number of items common to these two. 
    
    Arguments:
        query_gene {[str]} -- name of the gene to establish relationships with
        gene_common_count {[dict]} -- keys are genes, values are dictionary 
                                      lists with keys as genes (excluding the
                                      first) and values as the number of items
                                      the genes have in common
    
    Returns:
        [list] -- [List of 3-valued tuples: query_gene (node 1); subject_gene 
                                            (node 2) and number of items in
                                            common (weight of interaction)]
    """
    print("****** get_gene_edges ******")

    gene_dict = gene_common_count[query_gene] #Obtain dictionary for gene of interest
    edge_list = [] # Empty list of future vertices (term, gene, number of articles
                   # where the term and gene are mentioned)

    # Iterate through dictionary to establish relationships with the rest
    # (subjects) of the genes
    for subject_gene,num_articles in gene_dict.items(): 
        edge_list.append((query_gene, subject_gene, num_articles))
    return edge_list


def create_gene_node_trace(network,query_gene,edge_list,gene_common_count):
    """ Creates a layout of the scatter type (plotly) to visualize the network.
    In this case, the place where the nodes will be in an object of the scatter
    type have been defined. 
    
    Returns:
        [go] -- Returns a plotly scatter plot object. 
    """
    print("****** create_gene_node_trace ******")

    pos = nx.fruchterman_reingold_layout(network) # Defines a layout type of reingold
    node_trace = go.Scatter(
        x=[], # x-axis position of the node
        y=[],# y-axis position of the node
        text=[],
        mode='markers', # Markers
        name='ntw',
        hoverinfo='text', # Displays information when hovering over the point
        marker=dict(symbol='circle-dot',
            showscale=True, # Shows the scale
            colorscale='YlOrRd', # Color scale: yellow (low) to red (high)
            reversescale=True, 
            color=[], # List where the colors of each node will be encoded
            size=[], # List where the size of each node will be encoded
            colorbar=dict( # Color bar to show the color scale 
                thickness=30,
                title='Articles in common', 
                xanchor='left',
                titleside='right'),
            )
    )
    for node in network.nodes(): # Iterate through nodes of the graph
        x, y = pos[node] # Obtain positions from the map generated by reingold
        node_trace['x'] += tuple([x]) # x-axis position of the node
        node_trace['y'] += tuple([y]) # y-axis position of the node
    
    # Creates a dictionary to determine the number of genes in common between 
    # the gene subjects and the gene query, with color as an indicator
    nodecolor = dict([(gene, num_articles["weight"]) for term, gene, num_articles in network.edges(data=True)])
    nodecolor[query_gene] = 0 # Query_gene has no color (will be excluded from the graph)
    maxinter = max(nodecolor.values()) # Maximum scale value for query_gene: subject_gene that has more articles in common
    for node in nodecolor: # Normalize by the maximum scale value
        nodecolor[node] = int((nodecolor[node]/maxinter)*100)
    for node in network.nodes(): # Iterate through the nodes of the network (graph)
        if node != query_gene: # In the case that the node is different from the query_gene
            node_info = node + ": " + str(gene_common_count[query_gene][node]) # Information that appears in the node: number of articles per gene.
            node_trace['text'] += tuple([node_info]) # Put the information in the key text
            node_trace['marker']['color'] += tuple([int(nodecolor[node])]) # Escalation information - colors by number of interactions. Max is set to 100 for the most interactions. 
            node_trace['marker']['size'] += tuple([70]) # Fixed node size at 70.
        else: # When the node is the query_gene.
            node_trace['marker']['color'] += tuple([0]) # Set the color values 
            node_trace['marker']['size'] += tuple([0]) # Set size to 0
            node_trace['text'] += tuple([""]) # Do not show information
    return node_trace




def create_gene_name_trace(network,query_gene,node_trace,gene_common):
    """[summary]
    
    Arguments:
        network {[networx object]} -- Graph created from networkx
        query_gene {[str]} -- Gene of interest
        node_trace {[dict]} -- Dictionary that contains information of the nodes
                               in an object of the scatter type of plotly
        gene_common {[dict]} -- Dictionary: keys are genes, values are dictionary 
                                lists with keys as genes (excluding the first)
                                and values as items the genes have in common
     
    Returns:
        [dict] -- [names_nodes: dictionary that contains information about the 
                                names of the nodes and articles associated with
                                each node (gene)]
    """
    print("****** create_gene_name_trace ******")

    # Determines the names of nodes. The name of the node is composed of the name 
    # of the subject gene and contains a link to the articles it has in common 
    # with the query_gene
    names_nodes = ["<a href=" + "'https://www.ncbi.nlm.nih.gov/pubmed/{0}'".format(",".join(gene_common[query_gene][gene]).replace(".0","")) + 
                  'style="color: #000000">' + gene + 
                  "</a>" if gene != query_gene 
                  else "" for gene in list(network.nodes())]
    names_trace = go.Scatter( # Creates a scatter object
    x=node_trace["x"], # Defines the position of the name
    y=node_trace["y"], # Defines the position of the names in the place where 
                       # the nodes were previously determined (node_trace)
    text=names_nodes, # List containing the names of the nodes
    hoverinfo='none', # Turn off information display when hovering
    textposition='middle center', # Centers the text
    textfont=dict( # Deterimines the text font
        family='arial',
        size=20,
        color='#000000'
    ),
    mode='text')
    return names_trace

def gene_network_layout(query_gene,gene_common):
# def gene_network_layout(query_gene,gene_common,articles_by_gene):
    """Generate the layout of the graph.
    
    Arguments:
        query_gene {[str]} --  Name of the gene to establish relationships with
        articles_by_gene {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """
    print("****** gene_network_layout ******")

    # All items related to the query, in a single string. Remove duplicates 
    # from the set. With itertools chains from iterable, list of lists is 
    # transformed to a single list
    pubmed = ",".join(article_by_gene[query_gene]).replace(".0","")
    link = "<a href=" + "'https://www.ncbi.nlm.nih.gov/pubmed/{0}'".format(pubmed) + '>'+ query_gene +'</a>' # Links to all articles
    title = link + ": # Genes: " + str(len(gene_common[query_gene]))  + "; # Articles: " + str(len(article_by_gene[query_gene])) # Determine title
    axis = dict(showline=False,  # Hide axis line, grid, ticklabels and title
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        title=''
        )
    layout = go.Layout(title=title, # Set title style
        titlefont=dict(
            family='Gadugi',
            size=25,
            color='black'
        ),
        font=dict(size=15), # Set graph style
        plot_bgcolor='#EDEEF0',
        showlegend=False,
        autosize=True,
        height=800,
        xaxis=go.layout.XAxis(axis),
        yaxis=go.layout.YAxis(axis),
        margin=go.layout.Margin(
            l=100,
            r=100,
            b=100,
            t=100,
        ),
        annotations=[  # Set annotation style
            dict(
                showarrow=False,
                text=query_gene,
                xref='paper',
                yref='paper',
                x=0,
                y=-0.1,
                xanchor='left',
                yanchor='bottom',
                font=dict(
                    size=20
                )
            )
        ]
        )
    return layout


def create_gene_network(query_gene,article_by_gene,gene_common,gene_common_count):
    """Creacion del grafo. Centrado en un gen y la realcion con el resto de los genes (articulos en comun)
    
    Arguments:
        query_gene {[str]} --  Nombre del gen del cual quiero establecer relaciones
        article_by_gene {[dict]} -- Diccionario donde el key es el gen y el value es la lista de articulos relacionados al gen
        gene_common {[dict]} -- Diccionario donde el key es el gen y el value es la lista de diccionario donde los key son los genes (distintos al primero) y los values son los articulos que tienen los genes en comun
        gene_common_count {[dict]} -- Diccionario donde el key es el gen y el value es la lista de diccionario donde los key son los genes (distintos al primero) y los values son la cantidad de articulos que tienen los genes en comun
    
    Returns:
        [dict] -- fig: es un diccionario que contiene la figura del grafo, que es ploteado por dassh-plotly
    """
    print("****** create_gene_network ******")

    edge_list = get_gene_edges(query_gene,gene_common_count) # Creo la estructura de EDGES.
    if len(edge_list) == 0: #Cuando no hay resultados en una union o interseccion.
        return False
    network = nx.Graph() # genero el objeto de grafo vacio.
    node_list = [gene[1] for gene in edge_list] # nombre de los genes como nodos
    network.add_weighted_edges_from(edge_list) # Creo lista de arcos. edges.
    node_trace = create_gene_node_trace(network,query_gene,edge_list,gene_common_count)
    name_trace = create_gene_name_trace(network,query_gene,node_trace,gene_common)
    # layout = gene_network_layout(query_gene,gene_common,article_by_gene)
    layout = gene_network_layout(query_gene, gene_common)
    annot = "<a href='http://www.genomica.weebly.com'>IdMiner: Departamento de Genomica - IIBCE</a>"
    data = [node_trace, name_trace]
    fig = go.Figure(data=data, layout=layout)
    fig['layout']['annotations'][0]['text'] = annot
    return fig


def get_genes_df(contents, filename, date):
    print("****** get_genes_df ******")
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    global df_genes_relation
    global start_gene_query
    global article_by_gene
    global gene_common
    global gene_common_count
    global dfgenes
    try:
         # Global setting as it is a variable I will need, there are no 
         # problems. 
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            print("****** Reading csv file ******")
            dfgenes = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            dfgenes = pd.read_excel(io.BytesIO(decoded))
        article_by_gene,gene_common,gene_common_count = get_genes_relations(dfgenes)
        df_genes_relation, start_gene_query = generate_dataframe(article_by_gene,gene_common_count)
        return None
    except Exception as e:
        print(e,filename,"no correct format")



def parse_contents(df_genes_relation):
    print("****** parse_contents ******")
    try:
        col = df_genes_relation.columns
        return html.Div(children = [
            html.Hr(),
            dash_table.DataTable(
                id='gene_table-sorting-filtering',
                columns=[{"name": i, "id": i, 'deletable': True}
                        for i in df_genes_relation[col].columns],
                page_current = 0,
                page_size = 10,
                page_action='custom',
                sort_action='custom',
                sort_mode='multi',
                sort_by=[],
                filter_action='custom',
                filter_query='',
                style_table={'overflowX': 'scroll'},
                style_header={
                    'backgroundColor': '#91B9E5',
                    'minWidth': '0px', 'maxWidth': '800px',
                    'fontWeight': 'bold',
                    'font-family': 'inherit',
                    'textAlign': 'center',
                    'padding-right': '20px',
                    'padding-left': '20px',
                    'font-size': '15'
                },
                style_cell={
                    'backgroundColor': '#FAFAFA',
                    'minWidth': '0px', 'maxWidth': '800px',
                    'whiteSpace': 'no-wrap',
                    'overflow': 'hidden',
                    'textAlign': 'center',
                    'font-size': '15'
                }
            ),
            html.Hr(),
            dcc.Markdown('''#### Select Gene:'''),
            html.Div(
                id='gene-dropdown-container',
                children=[
                    dcc.Dropdown(
                        id='query-gene-dropdown',
                        options=[
                            {'label': i.title(), 'value': i} for i in sorted(dfgenes.columns[:-1])
                        ],
                        multi=False,
                        placeholder="dropdown menu please work",
                        disabled=False,
                        value=start_gene_query #Start query value
                    )
                ]
            ),
            html.Hr(),
            dcc.Graph(id='gene_net_graph')
        ]
        )
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])





#MAIN


print("****** uploadOrLoadSample ******")
uploadOrLoadSample = html.Div(
    className='flex-container',
    children=[
        dcc.Upload(
            id="upload-gene-data",
            className='dashed-file-upload',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Files')
                ],
                id='uploadedGenes'
            )),
        html.Div(id='output-gene-data-upload')
    ]
)

print("****** Callback to update filename ******")
@app.callback(Output('uploadedGenes','children'),
                [Input('upload-gene-data', 'contents')],
                [State('upload-gene-data', 'filename')])
def update_filename(contents, fname):
    if contents == None:
        raise PreventUpdate
    else:
        return fname


print("****** output-gene-data-upload layout ******")
layout = html.Div(
    children=[
        headerComponent_geneboard,
        html.Div(
            id="configuration-form-container",
        children=[
            html.H4('EXPLORING GENES',className='configuration-subsection'),
            uploadOrLoadSample,
            html.Div(id='output-gene-data-upload')
        ]
    )
    ]
)

print("****** Callback to update_output ******")
@app.callback(Output('output-gene-data-upload', 'children'),
              [Input('upload-gene-data', 'contents')],
              [State('upload-gene-data', 'filename'),
               State('upload-gene-data', 'last_modified')])
def update_output(file_content, file_name, file_date):
    if file_content is not None:
        if file_name.split("_")[-1] == "IDMiner-Genes.csv":
            get_genes_df(file_content, file_name, file_date)
            children = [parse_contents(df_genes_relation)]
        else:
             return html.Div(['There was an error processing this file. You need to upload this file, (yourname_IDMiner-Genes.csv)'])
        print("****** children returned ******")
        return children
        

print("****** Callback to update_graph ******")
@app.callback(
    Output('gene_table-sorting-filtering', 'data'),
    [Input('gene_table-sorting-filtering', 'page_current'),
     Input('gene_table-sorting-filtering', 'page_size'),
     Input('gene_table-sorting-filtering', 'sort_by'),
     Input('gene_table-sorting-filtering', 'filter_query')])
def update_graph(page_current, page_size, sort_by, filter_query):
    filtering_expressions = filter_query.split(' && ')
    dff = df_genes_relation
    for filter in filtering_expressions:
        if ' eq ' in filter:
            col_name = filter.split(' eq ')[0]
            filter_value = filter.split(' eq ')[1]
            dff = dff.loc[dff[col_name] == filter_value]
        if ' > ' in filter:
            col_name = filter.split(' > ')[0]
            filter_value = float(filter.split(' > ')[1])
            dff = dff.loc[dff[col_name] > filter_value]
        if ' < ' in filter:
            col_name = filter.split(' < ')[0]
            filter_value = float(filter.split(' < ')[1])
            dff = dff.loc[dff[col_name] < filter_value]

    if len(sort_by):
        dff = dff.sort_values(
            [col['column_id'] for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )

    return dff.iloc[page_current * page_size:
        (page_current + 1) * page_size
    ].to_dict('records')


print("****** app callback to load_graph ******")
@app.callback(
    Output('gene_net_graph', 'figure'), 
    [Input('gene-dropdown-container', 'children')])
def load_graph(selected_dropdown_value):
    print("Dropdown type: ", type(selected_dropdown_value))
    print("Dropdown value: ", selected_dropdown_value)
    graph = create_gene_network(selected_dropdown_value,article_by_gene,gene_common,gene_common_count)
    if graph:
        return graph
    else:
        return {}