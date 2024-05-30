# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 09:21:41 2023

@author: sleclerc

This program perform various data processing and analysis tasks related to gene expression data. Here is an overview of what the program does:

1. Imports necessary libraries: `gzip`, `json`, `pandas`, `networkx`, `matplotlib`, `plotly`, and `numpy`.
2. Sets global parameters such as the value of `fc` and `folder`.
3. Defines several functions:
   - `grab_merge_data`: Extracts specific data from a DataFrame based on gene names or gene synonyms.
   - `hist_from_dict_value`: Generates a histogram from values in a dictionary.
   - `select_junction`: Filters a DataFrame based on a list of protein names.
   - `grab_data`: Extracts specific data from a protein atlas based on gene names or gene synonyms.
   - `remove_link`: Removes links from a networkx graph based on a DataFrame and a specified value.
4. Opens a protein atlas file in JSON format and reads the data into a variable called `data`.
5. Reads gene expression data from an Excel file into a DataFrame called `groseq_data`.
6. Filters the `groseq_data` DataFrame based on certain conditions and assigns the filtered data to variables `up_4h`, `down_4h`, `up_8h`, and `down_8h`.
7. Defines a list of data fields to extract from the protein atlas.
8. Calls the `grab_merge_data` function to extract information for the up- and down-regulated proteins at 4h and 8h time points.
9. Combines the extracted information for up- and down-regulated proteins into dictionaries (`up_4h_info`, `down_4h_info`, `up_8h_info`, and `down_8h_info`).
10. Combines the information for all regulated proteins into dictionaries (`eight_h_info`, `for_h_info`, and `all_info`).
11. Creates dictionaries (`unique_4h_info` and `unique_8h_info`) containing information for uniquely regulated proteins at 4h and 8h time points.
12. Stores the combined information in a DataFrame called `data_df` and saves it as a pickle file.
13. Generates histogram figures and saves them as images.
14. Creates a Venn diagram to compare gene expression at different logFC values and saves it as an image.
15. Performs additional processing steps involving protein name mapping.

Require files (can be downloaded):
- 9606.protein.info.v11.5.txt.gz
- 9606.protein.links.full.v11.5.txt.gz
- proteinatlas.json.gz

Gro-Seq dataset:
- getDiffExpression_VERO_ChlSab1.1_noadj_rpkmAdded_GeneInfoAdded_EXCEL2.xlsx
"""
import gzip
import json
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import plotly.io as pio
pio.renderers.default='browser'
from matplotlib_venn import venn2
#global parameters
fc = 1 #can be 0.585, 1 or 2 for corresponding FC
folder = r'GRO-Seq image//'

#functions
def grab_merge_data(df, gene_name, to_extract):
    info = {}
    for gene in list(df['Gene name']):
        if not isinstance(gene, str):
            print('Not a named gene.')
            continue
        found = False
        if gene in gene_name:
            idx = gene_name.index(gene)
            sub_data = {k: data[idx][k] for k in to_extract}
            print('Found as '+gene)
            found = True
        else:
            for i in data:
                if len(i['Gene synonym']) > 0:
                    if gene in i['Gene synonym']:
                        sub_data = {k: i[k] for k in to_extract}
                        print('Found as '+i['Gene']+', not as '+gene)
                        found = True
                        break
                    
        if found == True:
            sel = df[df['Gene name'] == gene].squeeze()
            if isinstance(sel, pd.Series):
                sub_data['logFC_4h'] = sel['mock vs. HSV1_4h logFC']
                sub_data['logFC_8h'] = sel['mock vs. HSV1_8h logFC']
                sub_data['PValue_4h'] = sel['mock vs. HSV1_4h PValue']
                sub_data['PValue_8h'] = sel['mock vs. HSV1_8h PValue']
                sub_data['chr'] = sel['chr']
                sub_data['Gene_ID'] = sel['Gene stable ID']
                info[gene] = sub_data
            elif isinstance(sel, pd.DataFrame) and len(sel) > 1: #security, should not happen!
                print(gene+' has multiple entry in the GRO-seq database! '+str(len(sel)))
                sub_data['logFC_4h'] = sel['mock vs. HSV1_4h logFC'].max()
                sub_data['logFC_8h'] = sel['mock vs. HSV1_8h logFC'].max()
                sub_data['PValue_4h'] = sel['mock vs. HSV1_4h PValue'].max()
                sub_data['PValue_8h'] = sel['mock vs. HSV1_8h PValue'].max()
                sub_data['chr'] = sel['chr']
                sub_data['Gene_ID'] = sel['Gene stable ID']
                info[gene] = sub_data
    return info

def hist_from_dict_value(d, key, value=6):
    #d is dictionary, key is the dict key, value is the top most (number)
    ex = [x for xs in [d[x][key] for x in d if d[x][key] is not None] for x in xs]
    ex = np.unique(ex, return_counts=True)
    name = [str(x) for x in ex[0][ex[1].argsort()][-value:]]
    value = ex[1][ex[1].argsort()[-value:]]
    return name, value

def select_junction(new_df, selection_list):
    #new_df is the interactome, selection list is the list to isolate
    #return a df similar to new_df, but containing only the protein of selection
    selection = []
    for idx, row in new_df.iterrows():
        if row['protein1'] in selection_list:
            selection.append(idx)
            continue
        if row['protein2'] in selection_list:
            selection.append(idx)
    df = new_df.loc[selection]
    df = df.drop_duplicates()
    return df

def grab_data(gene_list, to_extract, data, verbose=True):
    #gene_list and to_extract should be a list of string
    #data should be the protein atlas
    gene_name = [x['Gene'] for x in data] #get the list of gene in the atlas
    info = {}
    for gene in gene_list:
        if not isinstance(gene, str): #check is entry is a string
            if verbose: print('Not a named gene.')
            continue
        found = False
        if gene in gene_name: #corresponding name for the protein
            idx = gene_name.index(gene)
            sub_data = {k: data[idx][k] for k in to_extract}
            if verbose: print('Found as '+gene)
            found = True
        else:
            for i in data: #look at the synonym
                if len(i['Gene synonym']) > 0:
                    if gene in i['Gene synonym']:
                        sub_data = {k: i[k] for k in to_extract}
                        if verbose: print('Found as '+i['Gene']+', not as '+gene)
                        found = True
                        break
        if found == True:
            info[gene] = sub_data
            
    return info

def remove_link(G, df, value):
    #G is a networkx object, df is the connectome dataframe, value is the number of connection (int)
    #return a dataframe with only proteins with the required number of connections
    dict_adjacencies = {}
    for i, adjacencies in enumerate(G.adjacency()):
        dict_adjacencies[adjacencies[0]] = len(adjacencies[1])
        
    pd_adj = pd.DataFrame.from_dict(dict_adjacencies, orient='index', columns=['connections'])
    to_remove = list(pd_adj[pd_adj.connections <= value].index)

    prots = []
    for idx, prot in df.iterrows():
        if prot['protein1'] in to_remove or prot['protein2'] in to_remove:
            continue
        prots.append(prot)
    filtered = pd.concat(prots, axis=1, ignore_index=True)
    return filtered.T


#%% open the protein atlas
with gzip.open('proteinatlas.json.gz', 'r') as fin:
    data = json.loads(fin.read().decode('utf-8'))

# get the expression data
groseq_data = pd.read_excel('getDiffExpression_VERO_ChlSab1.1_noadj_rpkmAdded_GeneInfoAdded_EXCEL2.xlsx')

up_4h = groseq_data[(groseq_data['mock vs. HSV1_4h logFC']>fc)&(groseq_data['mock vs. HSV1_4h FDR']<0.05)&(groseq_data['MAX_RPKM']>0.5)]
down_4h = groseq_data[(groseq_data['mock vs. HSV1_4h logFC']<-fc)&(groseq_data['mock vs. HSV1_4h FDR']<0.05)&(groseq_data['MAX_RPKM']>0.5)]

up_8h = groseq_data[(groseq_data['mock vs. HSV1_8h logFC']>fc)&(groseq_data['mock vs. HSV1_8h FDR']<0.05)&(groseq_data['MAX_RPKM']>0.5)]
down_8h = groseq_data[(groseq_data['mock vs. HSV1_8h logFC']<-fc)&(groseq_data['mock vs. HSV1_8h FDR']<0.05)&(groseq_data['MAX_RPKM']>0.5)]

#Data to extract
gene_name = [x['Gene'] for x in data] #get the list of gene
to_extract = ('Gene', 'Molecular function', 'Biological process', 'Subcellular location', 'Uniprot')

#up and down regulated protein information
up_4h_info = grab_merge_data(up_4h, gene_name, to_extract)
down_4h_info = grab_merge_data(down_4h, gene_name, to_extract)
up_8h_info = grab_merge_data(up_8h, gene_name, to_extract)
down_8h_info = grab_merge_data(down_8h, gene_name, to_extract)

#All regulated protein information (fusion of up and down)
eight_h_info = {**up_8h_info, **down_8h_info}
for_h_info = {**up_4h_info, **down_4h_info}
all_info = {**eight_h_info, **for_h_info} #fused

#unique protein found at a specific time
unique_4h_info = {x: for_h_info[x] for x in for_h_info if x not in eight_h_info.keys()}
unique_8h_info = {x: eight_h_info[x] for x in eight_h_info if x not in for_h_info.keys()}

#%% Make the database to save everything
ex = []
for k in all_info:
    if 'logFC_4h' in all_info[k].keys():
        if all_info[k]['logFC_4h'] >= fc or all_info[k]['logFC_4h'] <= -fc:
            info = [all_info[k]['chr'], all_info[k]['Gene_ID'],
                    all_info[k]['Gene'], all_info[k]['logFC_4h'], '4h']
            ex.append(info)
    if 'logFC_8h' in all_info[k].keys():
        if all_info[k]['logFC_8h'] >= fc or all_info[k]['logFC_8h'] <= -fc:
            info = [all_info[k]['chr'], all_info[k]['Gene_ID'],
                    all_info[k]['Gene'], all_info[k]['logFC_8h'], '8h']
            ex.append(info)

data_df = pd.DataFrame(data=ex, columns=['Chr', 'Nearest Ensembl', 'gene name', 'logFC', 'condition'])

data_df.to_pickle('GRO_data.pckl')

#%% histogram figures!
fig, ax = plt.subplots(ncols=2, nrows=3, figsize=(10,15))

#4h
_proc, _c = hist_from_dict_value(unique_4h_info, 'Biological process')
ax[0,0].bar(_proc, _c)
ax[0,0].set_title('4h unique - Biological process')
ax[0,0].set_xticklabels(_proc, rotation=45, ha='right')
_proc, _c = hist_from_dict_value(unique_4h_info, 'Molecular function')
ax[1,0].bar(_proc, _c)
ax[1,0].set_title('4h unique - Molecular function')
ax[1,0].set_xticklabels(_proc, rotation=45, ha='right')
_proc, _c = hist_from_dict_value(unique_4h_info, 'Subcellular location')
ax[2,0].bar(_proc, _c)
ax[2,0].set_title('4h unique - Subcellular location')
ax[2,0].set_xticklabels(_proc, rotation=45, ha='right')

#8h
_proc, _c = hist_from_dict_value(unique_8h_info, 'Biological process')
ax[0,1].bar(_proc, _c)
ax[0,1].set_title('8h unique - Biological process')
ax[0,1].set_xticklabels(_proc, rotation=45, ha='right')
_proc, _c = hist_from_dict_value(unique_8h_info, 'Molecular function')
ax[1,1].bar(_proc, _c)
ax[1,1].set_title('8h unique - Molecular function')
ax[1,1].set_xticklabels(_proc, rotation=45, ha='right')
_proc, _c = hist_from_dict_value(unique_8h_info, 'Subcellular location')
ax[2,1].bar(_proc, _c)
ax[2,1].set_title('8h unique - Subcellular location')
ax[2,1].set_xticklabels(_proc, rotation=45, ha='right')

plt.tight_layout()
plt.savefig(folder+'unique_protein_hist.png', dpi=300)
plt.show()

fig, ax = plt.subplots(ncols=2, nrows=3, figsize=(10,15))
#4h
_proc, _c = hist_from_dict_value(for_h_info, 'Biological process')
ax[0,0].bar(_proc, _c)
ax[0,0].set_title('All 4h proteins - Biological process')
ax[0,0].set_xticklabels(_proc, rotation=45, ha='right')
_proc, _c = hist_from_dict_value(for_h_info, 'Molecular function')
ax[1,0].bar(_proc, _c)
ax[1,0].set_title('All 4h proteins - Molecular function')
ax[1,0].set_xticklabels(_proc, rotation=45, ha='right')
_proc, _c = hist_from_dict_value(for_h_info, 'Subcellular location')
ax[2,0].bar(_proc, _c)
ax[2,0].set_title('All 4h proteins - Subcellular location')
ax[2,0].set_xticklabels(_proc, rotation=45, ha='right')

#8h
_proc, _c = hist_from_dict_value(eight_h_info, 'Biological process')
ax[0,1].bar(_proc, _c)
ax[0,1].set_title('All 8h proteins - Biological process')
ax[0,1].set_xticklabels(_proc, rotation=45, ha='right')
_proc, _c = hist_from_dict_value(eight_h_info, 'Molecular function')
ax[1,1].bar(_proc, _c)
ax[1,1].set_title('All 8h proteins - Molecular function')
ax[1,1].set_xticklabels(_proc, rotation=45, ha='right')
_proc, _c = hist_from_dict_value(eight_h_info, 'Subcellular location')
ax[2,1].bar(_proc, _c)
ax[2,1].set_title('All 8h proteins - Subcellular location')
ax[2,1].set_xticklabels(_proc, rotation=45, ha='right')

plt.tight_layout()
plt.savefig(folder+'all_protein_hist.png', dpi=300)
plt.show()

#%% Venn diagram  
plt.figure(figsize=(5,5))
v = venn2([for_h_info.keys(),eight_h_info.keys()], ['4 hpi', '8 hpi'])
plt.title('Compare genes expression at LogFC='+str(fc))
plt.savefig(folder+'Venn Comparison.png')
plt.show()


#%% Get the link between ENPS format to classical protein name
prot_name = '9606.protein.info.v11.5.txt.gz'
prot_name = pd.read_csv(prot_name, sep='\t')

enps_prot = {}
for key in all_info:
    sel = prot_name[prot_name['preferred_name']==all_info[key]['Gene']].squeeze()
    if len(sel) == 0:
        continue
    
    prot = sel['#string_protein_id']
    if prot not in enps_prot: #security, should not happen since working with dict
        enps_prot[prot] = key

#%% read the human interactome
interactome_name = '9606.protein.links.full.v11.5.txt.gz'
df = pd.read_csv(interactome_name, sep=' ')

#select only the interacting protein. Slow!
enps_prot_list = enps_prot.keys()
selection = {}
for idx, row in df.iterrows():
    if row.experiments > 500: #value to change!
        if row.protein1 in enps_prot_list:
            selection[idx] = [enps_prot[row.protein1],
                              prot_name[prot_name['#string_protein_id']==row.protein2].squeeze()['preferred_name'],
                              row.experiments]
            continue
        if row.protein2 in enps_prot_list:
            selection[idx] = [prot_name[prot_name['#string_protein_id']==row.protein1].squeeze()['preferred_name'],
                              enps_prot[row.protein2],
                              row.experiments]
            continue
        
new_df = pd.DataFrame.from_dict(selection).T
new_df = new_df.rename(columns={0:'protein1', 1:'protein2', 2:'score'})
print('Finish to read the human interactome')
#%% Network generation
def mpl_graph_draw(df, prot_list, time, title, largest=True):
    """
    Draw a Matplotlib/Networkx plot. Difficult to read if there are a lot of proteins.
    The thickness of the edge is proportional to the protein-protein connection
    The size of the node is proportional to its LogFC score
    The color of the node indicate:
        if white, a middle man protein, not part of the prot_list
        if blue, part of the prot_list
        if yellow, part of the prot_list AND in the location of interest
        if pale yellow, middle man protein at the location of interest
    Parameters
    ----------
    df : DataFrame of connecting edges (protein1, protein2, score)
    prot_list : list of protein (str) being the main actors
        Will be draw light blue node
    time : str
        The LogFC to focus. either 4h or 8h
    title : str
        The title.
    largest : bool, optional
        If selecting the largest subnetwork or not. The default is True.

    """
    G_main = nx.Graph() #main graph
    for idx, row in df.iterrows():
        G_main.add_edge(row.protein1, row.protein2, weight=row.score)
    
    #Select only the largest subgraph
    sub_network = sorted(nx.connected_components(G_main), key=len, reverse=True)
    if largest:
        sub_network = [sub_network[0]] #new graph
    
    blue_cmap = cm.get_cmap('Blues') #Down cell
    reds_cmap = cm.get_cmap('Purples') #Up cell
    oran_cmap = cm.get_cmap('Oranges') #Up in mito only
    gree_cmap = cm.get_cmap('Greens') #down in mito only
    #middle man white, silver if in mito
    cmap_div = [5, 5] #first is higher than 0, second is lower than 0
    idx = 0
    for network in sub_network: #for each subgraph
        plt.figure(figsize=(40,40))
        G = G_main.subgraph(network)
        if G.number_of_nodes() < 4:
            print('No more subgraph')
            break
        #add color to graph based on their owner or properties
        #Node size is function of LogFC
        color_map = []
        node_size = []
        for node in G.nodes:
            if node in sele_func[i].keys(): #change here!
                c = 'skyblue'
                n = all_prot[node]
                s = max(1000, abs(int(1000*n['logFC_'+time])))
                if n['logFC_'+time] < 0: c = blue_cmap(min(1, abs(n['logFC_'+time])/cmap_div[1]))
                if n['logFC_'+time] > 0: c = reds_cmap(min(1, n['logFC_'+time]/cmap_div[0]))
                if isinstance(n['Subcellular location'], list): #to not grab a NaN
                    if focus_loc in n['Subcellular location']:
                        c = 'orange'
                        if n['logFC_'+time] < 0: c = gree_cmap(min(1, abs(n['logFC_'+time])/cmap_div[1]))
                        if n['logFC_'+time] > 0: c = oran_cmap(min(1, n['logFC_'+time]/cmap_div[0]))

                        
            else: #a middle man
                c = 'snow'
                s = 1000
                if node in interact_prot_info.keys():
                    n = interact_prot_info[node]
                    if isinstance(n['Subcellular location'], list): #to not grab a NaN
                        if focus_loc in n['Subcellular location']:
                            c = 'silver'
            color_map.append(c)
            node_size.append(s)
        
        #need to choose between the two distribution
        pos = nx.spring_layout(G, k=0.3*1/np.sqrt(len(G.nodes())))
        pos = nx.spring_layout(G)
        for k in pos.keys(): #to obtain positive only node position, allowing to square it
            pos[k] = (pos[k]+3)**2
        
        print('Start to draw function')
        nx.draw_networkx_nodes(G, pos, node_color=color_map, node_size=node_size,
                               edgecolors='black')
        #connection thickness between nodes indicate the score
        weights = np.array([G[u][v]['weight'] for u,v in G.edges()], dtype=int)
        weights = (((weights-np.min(weights))/(np.max(weights)-np.min(weights)))+0.5)
        nx.draw_networkx_edges(G, pos, width=weights*3)
        
        #Write the protein name in the node
        prot = {k: v for k, v in pos.items() if k in prot_list}
        not_prot = {k: v for k, v in pos.items() if k not in prot_list}
        nx.draw_networkx_labels(G.subgraph(prot), prot, font_weight='bold')
        nx.draw_networkx_labels(G.subgraph(not_prot), not_prot)
        
        plt.title(title)
        plt.tight_layout()
        plt.savefig(folder+'protein-subnetwork_'+str(idx)+'_'+title+'.svg')
        plt.savefig(folder+'protein-subnetwork_'+str(idx)+'_'+title+'.png', dpi=75)
        plt.show()
        idx += 1


#%%MATPLOTLIB!
sele_func = [for_h_info, eight_h_info] #data to use
time = ['4h', '8h'] #time point for the LogFC
kind = ['short', 'long']

focus_loc = 'Mitochondria' #localization to focus on
connection_n = 4 #number of connection for middle man to exist - LONG only


for i, t in enumerate(time):
    for ii, k in enumerate(kind):
        print(t, k)
    
        #Data selection
        sele_df = select_junction(new_df, sele_func[i].keys())
        #grab all interactors NOT in the dataset
        all_prot = set(list(sele_df['protein1']) + list(sele_df['protein2']))
        interact_prot = all_prot - set(sele_func[i].keys())
        #get their information as a dict
        interact_prot_info = grab_data(list(interact_prot), to_extract, data)
        all_prot = {**sele_func[i], **interact_prot_info}
        
        if k == 'long':
            #initial network, allow to calculate and drop low interaction link
            G = nx.from_pandas_edgelist(sele_df, 'protein1', 'protein2')
            filtered = remove_link(G, sele_df, connection_n)
            mpl_graph_draw(filtered, sele_func[i].keys(), t, t+'_'+k, largest=False)
            
        if k == 'short':
            #short network (without any middle man)
            short = []
            for idx, row in sele_df.iterrows():
                if row['protein1'] in sele_func[i].keys() and row['protein2'] in sele_func[i].keys():
                    short.append(idx)
            short = sele_df.loc[short]
            short = short.drop_duplicates()
            G = nx.from_pandas_edgelist(short, 'protein1', 'protein2')
            filtered = remove_link(G, short, 0)
            mpl_graph_draw(filtered, sele_func[i].keys(), t, t+'_'+k, largest=False)
        
        
#%%MATPLOTLIB Mitochondria focused
sele_func = [for_h_info, eight_h_info] #data to use
time = ['4h', '8h'] #time point for the LogFC

focus_loc = 'Mitochondria' #localization to focus on
connection_n = 0 #number of connection for middle man to exist - LONG only


for i, t in enumerate(time):
    print(t)
    #select protein with Mitochondria in their location
    selected_prot = {}
    for prot in sele_func[i].keys():
        if sele_func[i][prot]['Subcellular location'] is not None:
            if 'Mitochondria' in sele_func[i][prot]['Subcellular location']:
                selected_prot[prot] = sele_func[i][prot]
    #Data selection
    sele_df = select_junction(new_df, selected_prot.keys())
    #grab all interactors NOT in the dataset
    all_prot = set(list(sele_df['protein1']) + list(sele_df['protein2']))
    interact_prot = all_prot - set(sele_func[i].keys())
    #get their information as a dict
    interact_prot_info = grab_data(list(interact_prot), to_extract, data, verbose=False)
    all_prot = {**sele_func[i], **interact_prot_info}
    
    #initial network, allow to calculate and drop low interaction link
    G = nx.from_pandas_edgelist(sele_df, 'protein1', 'protein2')
    filtered = remove_link(G, sele_df, connection_n)
    mpl_graph_draw(filtered, selected_prot.keys(), t, 'Mito_focus at '+t+'_'+k, largest=False)
    



#%%
from PIL import Image

blue_cmap = cm.get_cmap('Purples') #Down cell
reds_cmap = cm.get_cmap('Reds') #Up cell
oran_cmap = cm.get_cmap('Oranges') #Up in mito only
gree_cmap = cm.get_cmap('Greens') #down in mito only



v4 = np.array([abs(for_h_info[x]['logFC_4h']) for x in for_h_info])
v8 = np.array([abs(eight_h_info[x]['logFC_8h']) for x in eight_h_info])

new_x = np.round(np.histogram(v8, bins=50)[1][np.histogram(v8, bins=50)[1] < 5],
                 decimals=2)

# new_x = new_x[::-1]

im = Image.fromarray(np.expand_dims(np.uint8(blue_cmap(new_x/5)*255), axis=1))


plt.imshow(im, aspect=0.2)
x = plt.yticks()[0]

plt.yticks(np.arange(0, max(x), max(x)/len(x)),
           np.round(np.arange(0, max(new_x), max(new_x)/len(x)), decimals=3)+1)

plt.xticks(np.arange(0))
plt.savefig('Purple_8h_cmap.svg')
plt.show()





