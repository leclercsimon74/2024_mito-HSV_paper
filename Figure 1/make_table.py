# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 09:03:23 2023

@author: sleclerc
"""

import gzip
import json
import pandas as pd
import numpy as np
#global parameters
fc = 1 #can be 0.585, 1 or 2 for corresponding FC
folder = r'table//'


#functions
def grab_merge_data(df, gene_name, to_extract):
    info = {}
    for gene in list(df['Gene name']):
        new_name = False
        if not isinstance(gene, str):
            #print('Not a named gene.')
            continue
        found = False
        if gene in gene_name:
            idx = gene_name.index(gene)
            sub_data = {k: data[idx][k] for k in to_extract}
            #print('Found as '+gene)
            found = True
        else:
            for i in data:
                if len(i['Gene synonym']) > 0:
                    if gene in i['Gene synonym']:
                        sub_data = {k: i[k] for k in to_extract}
                        print('Found as '+i['Gene']+', not as '+gene)
                        found = True
                        new_name = True
                        break
        if found == False:
            continue
        
        sel = df[df['Gene name'] == gene].squeeze()
        
        if new_name: gene = i['Gene']
            
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

full_df = pd.DataFrame.from_dict(all_info)
full_df = full_df.T
full_df.to_csv(folder+'full_table.csv')

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

#%% count the fav localization
fav = 'Mitochondria'
up4 = full_df[(full_df['logFC_4h']>=1) &
              (full_df['PValue_4h']<=0.05) &
              (full_df['Subcellular location'].astype(str).str.contains(fav))]

down4 = full_df[(full_df['logFC_4h']<=-1) &
                (full_df['PValue_4h']<=0.05) &
                (full_df['Subcellular location'].astype(str).str.contains(fav))]

up8 = full_df[(full_df['logFC_8h']>=1) &
              (full_df['PValue_8h']<=0.05) &
              (full_df['Subcellular location'].astype(str).str.contains(fav))]

down8 = full_df[(full_df['logFC_8h']<=-1) &
                (full_df['PValue_8h']<=0.05) &
                (full_df['Subcellular location'].astype(str).str.contains(fav))]

unique_up = len(set(list(up4['Gene'])+list(up8['Gene'])))
unique_down = len(set(list(down4['Gene'])+list(down8['Gene'])))


table = [['','Upregulated', 'Downregulated', 'Total'],
         ['4 hpi',len(up4), len(down4), len(up4)+len(down4)],
         ['8 hpi', len(up8), len(down8), len(up8)+len(down8)],
         ['Total', unique_up, unique_down, unique_up+unique_down]]
for row in table:
    print('| {:1} | {:1} | {:1} | {:1} |'.format(*row))

unique = list(set(list(up4['Gene'])+list(up8['Gene']))) + list(set(list(down4['Gene'])+list(down8['Gene'])))

fav_df = full_df[full_df.Gene.isin(unique)]
fav_df.to_csv(folder+fav+'_table.csv')

# #%% read the human interactome
# interactome_name = '9606.protein.links.full.v11.5.txt.gz'
# df = pd.read_csv(interactome_name, sep=' ')

# #select only the interacting protein. Slow!
# enps_prot_list = enps_prot.keys()
# selection = {}
# for idx, row in df.iterrows():
#     if row.experiments > 500: #value to change!
#         if row.protein1 in enps_prot_list:
#             selection[idx] = [enps_prot[row.protein1],
#                               prot_name[prot_name['#string_protein_id']==row.protein2].squeeze()['preferred_name'],
#                               row.experiments]
#             continue
#         if row.protein2 in enps_prot_list:
#             selection[idx] = [prot_name[prot_name['#string_protein_id']==row.protein1].squeeze()['preferred_name'],
#                               enps_prot[row.protein2],
#                               row.experiments]
#             continue
        
# new_df = pd.DataFrame.from_dict(selection).T
# new_df = new_df.rename(columns={0:'protein1', 1:'protein2', 2:'score'})
# print('Finish to read the human interactome')
# new_df.to_csv(folder+'Interactome.csv')

#%%Venn diagram full cell + mito focus
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

#full cell
plt.figure(figsize=(10,10))
v = venn2([for_h_info.keys(),eight_h_info.keys()], ['4 hpi', '8 hpi'])
plt.title('Compare genes expression at LogFC='+str(fc))
plt.savefig(folder+'Venn Comparison full data.png', dpi=600)
plt.show()

#mito focus

def extract_loc_dict(dictionnary, location):
    in_ = {}
    for item in dictionnary:
        t = dictionnary[item]['Subcellular location']
        if isinstance(t, list):
            if location in t:
                in_[item] = dictionnary[item]
    return in_

for_h_mito_info = extract_loc_dict(for_h_info, 'Mitochondria')
eight_h_mito_info = extract_loc_dict(eight_h_info, 'Mitochondria')

plt.figure(figsize=(10,10))
v = venn2([for_h_mito_info.keys(), eight_h_mito_info.keys()],
          ['4 hpi', '8 hpi'])
plt.title('Compare genes expression in mitochondria at LogFC='+str(fc))
plt.savefig(folder+'Venn Comparison mitochondria focus.png', dpi=600)
plt.show()


#%% Dual direction bar plot for mito only
# the idea is to categorize the data by their biological process
# Then up and down regulated count (in corresponding direction)
# Then split the 4 and 8 hpi (color)

graph_value = {}
key = 'Biological process'
d = {**for_h_mito_info, **eight_h_mito_info}

for item in d:
    t = d[item][key]
    if isinstance(t, list):
        for i in t:
            if i not in graph_value:
                dd = {'up_4h':0,'down_4h':0,
                      'up_8h':0,'down_8h':0,
                      'count':0}
                graph_value[i] = dd
            
            if d[item]['logFC_4h'] > 1:
                graph_value[i]['up_4h'] += 1
            elif d[item]['logFC_4h'] < -1:
                graph_value[i]['down_4h'] -= 1
            if d[item]['logFC_4h'] > 1 or d[item]['logFC_4h'] < -1:
                graph_value[i]['count'] +=1
                
            if d[item]['logFC_8h'] > 1:
                graph_value[i]['up_8h'] += 1
            elif d[item]['logFC_8h'] < -1:
                graph_value[i]['down_8h'] -= 1
            if d[item]['logFC_8h'] > 1 or d[item]['logFC_8h'] < -1:
                graph_value[i]['count'] +=1      


clean_graph_value = {}
for item in graph_value:
    if graph_value[item]['count'] >= 5:
        clean_graph_value[item] = graph_value[item]

df = pd.DataFrame.from_dict(clean_graph_value)
df = df.T
df = df.sort_index()


name = list(df.index)
x = np.arange(len(name))  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots(layout='constrained', figsize=(10,10))
ax.grid(True, linestyle='--', color='#D3D3D3', zorder=0)

ax.bar(x, np.array(df['up_4h']), width, label='4hpi upregulated', color='#A4D8D8', zorder=3)
ax.bar(x, np.array(df['down_4h']), width, label='4hpi downregulated', color='#F1B3A9', zorder=3)   
ax.bar(x + width, np.array(df['up_8h']), width, label='8hpi upregulated', color='#047c99', zorder=3)
ax.bar(x + width, np.array(df['down_8h']), width, label='8hpi downregulated', color='#C3B1E1', zorder=3)  
ax.axhline(y=0,  color='k', zorder=5, linewidth=0.5)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Count')
ax.set_title('Gro-Seq protein biological process in Mitochondria')
ax.set_xticks(x + width, name, rotation=45, ha='right')
ticks =  ax.get_yticks()
ax.set_yticklabels([int(abs(tick)) for tick in ticks])
ax.legend(loc='upper left', ncols=2)
ax.grid(True)
plt.savefig(folder+'Gro-Seq protein biological process in Mitochondria.png', dpi=600)
plt.show()



#%% Same on the full dataset
graph_value = {}
key = 'Biological process'
d = {**for_h_info, **eight_h_info}

for item in d:
    t = d[item][key]
    if isinstance(t, list):
        for i in t:
            if i not in graph_value:
                dd = {'up_4h':0,'down_4h':0,
                      'up_8h':0,'down_8h':0,
                      'count':0}
                graph_value[i] = dd
            
            if d[item]['logFC_4h'] > 1:
                graph_value[i]['up_4h'] += 1
            elif d[item]['logFC_4h'] < -1:
                graph_value[i]['down_4h'] -= 1
            if d[item]['logFC_4h'] > 1 or d[item]['logFC_4h'] < -1:
                graph_value[i]['count'] +=1
                
            if d[item]['logFC_8h'] > 1:
                graph_value[i]['up_8h'] += 1
            elif d[item]['logFC_8h'] < -1:
                graph_value[i]['down_8h'] -= 1
            if d[item]['logFC_8h'] > 1 or d[item]['logFC_8h'] < -1:
                graph_value[i]['count'] +=1            

clean_graph_value = {}
for item in graph_value:
    if graph_value[item]['count'] >= 50:
        clean_graph_value[item] = graph_value[item]

df = pd.DataFrame.from_dict(clean_graph_value)
df = df.T
df = df.sort_index()


name = list(df.index)
x = np.arange(len(name))  # the label locations
width = 0.25  # the width of the bars

# c = u.get_color(4, 'bar')

fig, ax = plt.subplots(layout='constrained', figsize=(10,10))
ax.grid(True, linestyle='--', color='#D3D3D3', zorder=0)

ax.bar(x, np.array(df['up_4h']), width, label='4hpi upregulated', color='#A4D8D8', zorder=3)
ax.bar(x, np.array(df['down_4h']), width, label='4hpi downregulated', color='#F1B3A9', zorder=3)   
ax.bar(x + width, np.array(df['up_8h']), width, label='8hpi upregulated', color='#047c99', zorder=3)
ax.bar(x + width, np.array(df['down_8h']), width, label='8hpi downregulated', color='#C3B1E1', zorder=3)  
ax.axhline(y=0,  color='k', zorder=5, linewidth=0.5)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Count')
ax.set_title('Gro-Seq protein biological process')
ax.set_xticks(x + width, name, rotation=45, ha='right')
ticks =  ax.get_yticks()
ax.set_yticklabels([int(abs(tick)) for tick in ticks])
ax.legend(loc='upper left', ncols=2)

plt.savefig(folder+'Gro-Seq protein biological process.png', dpi=600)
plt.show()

#%% Own annotation based on Uniprot summary
#load the data
df = pd.read_csv(folder+'Mitochondria_table - Modify.csv')
df['Mitochondrial function'] = df['Mitochondrial function'].tolist()

graph_value = {}
key = 'Mitochondrial function'
d = df.T.to_dict()

for item in d:
    t = d[item][key]
    if not isinstance(t, float):
        t = [x for x in t.split('\'') if x not in '[], ']
        for i in t:
            if i not in graph_value:
                dd = {'up_4h':0,'down_4h':0,
                      'up_8h':0,'down_8h':0,
                      'count':0}
                graph_value[i] = dd
            
            if d[item]['logFC_4h'] > 1:
                graph_value[i]['up_4h'] += 1
            elif d[item]['logFC_4h'] < -1:
                graph_value[i]['down_4h'] -= 1
            if d[item]['logFC_4h'] > 1 or d[item]['logFC_4h'] < -1:
                graph_value[i]['count'] +=1
                
            if d[item]['logFC_8h'] > 1:
                graph_value[i]['up_8h'] += 1
            elif d[item]['logFC_8h'] < -1:
                graph_value[i]['down_8h'] -= 1
            if d[item]['logFC_8h'] > 1 or d[item]['logFC_8h'] < -1:
                graph_value[i]['count'] +=1            

clean_graph_value = {}
for item in graph_value:
    if graph_value[item]['count'] >= 5:
        clean_graph_value[item] = graph_value[item]

df = pd.DataFrame.from_dict(clean_graph_value)
df = df.T
df = df.sort_index()


name = list(df.index)
x = np.arange(len(name))  # the label locations
width = 0.25  # the width of the bars

# c = u.get_color(4, 'bar')

fig, ax = plt.subplots(layout='constrained', figsize=(10,10))
ax.grid(True, linestyle='--', color='#D3D3D3', zorder=0)

ax.bar(x, np.array(df['up_4h']), width, label='4hpi upregulated', color='#A4D8D8', zorder=3)
ax.bar(x, np.array(df['down_4h']), width, label='4hpi downregulated', color='#F1B3A9', zorder=3)   
ax.bar(x + width, np.array(df['up_8h']), width, label='8hpi upregulated', color='#047c99', zorder=3)
ax.bar(x + width, np.array(df['down_8h']), width, label='8hpi downregulated', color='#C3B1E1', zorder=3)  
ax.axhline(y=0,  color='k', zorder=5, linewidth=0.5)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Count')
ax.set_title('Gro-Seq protein Mitochondria function')
ax.set_xticks(x + width, name, rotation=45, ha='right')
ticks =  ax.get_yticks()
ax.set_yticklabels([int(abs(tick)) for tick in ticks])
ax.legend(loc='upper left', ncols=2)

plt.savefig(folder+'Gro-Seq protein mitochondria function.png', dpi=600)
plt.show()


#%% GO style
import json
import os

#extract value from GO based and classify based on the protein name
with open('GO'+os.sep+'analysis.json') as file:
    js = json.load(file)
    js = js['overrepresentation']

data = js['group']
protein_GO = {}
level = 0 #0 to 5

for d in data:
    sub = d['result']
    
    if isinstance(sub, list):
        for s in sub:
            inp = s['input_list']
            t = s['term']
            if t['level'] != level:
                continue
            l = t['label']
            i = t['id']        
            for prot in inp['mapped_id_list']['mapped_id']:
                if prot in protein_GO.keys():
                    protein_GO[prot]['label'].append(l)
                    protein_GO[prot]['id'].append(i)
                else:
                    dd = {'label':[l],
                          'id':[i]}
                    protein_GO[prot] = dd
                    
    if isinstance(sub, dict):
        inp = sub['input_list']
        t = sub['term']
        if t['level'] > level:
            continue
        l = t['label']
        i = t['id']        
        for prot in inp['mapped_id_list']['mapped_id']:
            if prot in protein_GO.keys():
                protein_GO[prot]['label'].append(l)
                protein_GO[prot]['id'].append(i)
            else:
                dd = {'label':[l],
                      'id':[i]}
                protein_GO[prot] = dd        
                
        
graph_value = {}

d = {**for_h_mito_info, **eight_h_mito_info}

        
for item in protein_GO:
    t = protein_GO[item]['label']
    if isinstance(t, list):
        for i in t:
            if i not in graph_value:
                dd = {'up_4h':0,'down_4h':0,
                      'up_8h':0,'down_8h':0,
                      'count':0}
                graph_value[i] = dd
            
            if d[item]['logFC_4h'] > 1:
                graph_value[i]['up_4h'] += 1
            elif d[item]['logFC_4h'] < -1:
                graph_value[i]['down_4h'] -= 1
            if d[item]['logFC_4h'] > 1 or d[item]['logFC_4h'] < -1:
                graph_value[i]['count'] +=1
                
            if d[item]['logFC_8h'] > 1:
                graph_value[i]['up_8h'] += 1
            elif d[item]['logFC_8h'] < -1:
                graph_value[i]['down_8h'] -= 1
            if d[item]['logFC_8h'] > 1 or d[item]['logFC_8h'] < -1:
                graph_value[i]['count'] +=1            

clean_graph_value = {}
for item in graph_value:
    if graph_value[item]['count'] >= 2:
        clean_graph_value[item] = graph_value[item]

df = pd.DataFrame.from_dict(clean_graph_value)
df = df.T
df = df.sort_index()


name = list(df.index)
x = np.arange(len(name))  # the label locations
width = 0.25  # the width of the bars

# c = u.get_color(4, 'bar')

fig, ax = plt.subplots(layout='constrained', figsize=(10,10))
ax.grid(True, linestyle='--', color='#D3D3D3', zorder=0)

ax.bar(x, np.array(df['up_4h']), width, label='4hpi upregulated', color='#A4D8D8', zorder=3)
ax.bar(x, np.array(df['down_4h']), width, label='4hpi downregulated', color='#F1B3A9', zorder=3)   
ax.bar(x + width, np.array(df['up_8h']), width, label='8hpi upregulated', color='#047c99', zorder=3)
ax.bar(x + width, np.array(df['down_8h']), width, label='8hpi downregulated', color='#C3B1E1', zorder=3)  
ax.axhline(y=0,  color='k', zorder=5, linewidth=0.5)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Count')
ax.set_title('Gro-Seq protein Mitochondria function (GO) at level '+str(level))
ax.set_xticks(x + width, name, rotation=45, ha='right')
ticks =  ax.get_yticks()
ax.set_yticklabels([int(abs(tick)) for tick in ticks])
ax.legend(loc='upper left', ncols=2)

plt.savefig(folder+'Gro-Seq protein mitochondria function from GO at level='+str(level)+'.png', dpi=600)
plt.show()

#%% Sankey diagram

import random
import numpy as np

#extract value from GO based and classify based on the protein name
with open('GO'+os.sep+'analysis.json') as file:
    js = json.load(file)
    js = js['overrepresentation']

data = js['group']
table = []
max_lvl = 15

for d in data:
    sub = d['result']
    hierarchy = {x:[] for x in range(max_lvl)}
    go_hie = {x:[] for x in range(max_lvl)}
    color = '#%02X%02X%02X' % (random.randint(0,255),
                               random.randint(0,255),
                               random.randint(0,255))
    
    if isinstance(sub, dict):
        inp = sub['input_list']
        t = sub['term']
        hierarchy[t['level']].append(inp['mapped_id_list']['mapped_id'])
        go_hie[t['level']].append(t['label'])        
    else:
        for s in sub:
            inp = s['input_list']
            t = s['term']
            hierarchy[t['level']].append(inp['mapped_id_list']['mapped_id'])
            go_hie[t['level']].append(t['label'])
    
    sub_hie = {x:set() for x in range(max_lvl)}
    for prot_sublist in hierarchy:

        for i, prot_list in enumerate(hierarchy[prot_sublist]):
            prots = set(prot_list)

            value = len(sub_hie[i].symmetric_difference(prots))
            go = go_hie[prot_sublist][i]


            re_go = {x:'' for x in range(max_lvl)}
            for prot_sublist2 in hierarchy:
                for j, prot_list2 in enumerate(hierarchy[prot_sublist2]):
                    #if prot_list == prot_list2: #same hierarchy, ignore
                    if prots.issubset(prot_list2):
                        re_go[prot_sublist2] += go_hie[prot_sublist2][j]+'\n'
                        sub_hie[prot_sublist2].update(prots)
              
                        

            if len(re_go[0]) > 0 or len(re_go[3]) > 5:

                re_go['value'] = value
                re_go['color'] = color
                table.append(re_go)
                

    

wished_org = ['level5', 'level4', 'level3', 'level2', 'level1', 'level0', 'value', 'color']
df_go = pd.DataFrame(table)
df_go = df_go.rename(columns={0:'level0', 1:'level1', 2:'level2', 3:'level3', 4:'level4', 5:'level5'})
df_go = df_go[wished_org]
df_go[df_go == ''] = np.nan
df_go = df_go[df_go.value>0]
df_go = df_go.reset_index()

#%% Own Sankey
#with fix y position for the node - not done
#1 color by network - DONE
#be sure to separate the mother node - DONE
# trouble with too long name
# trouble with the level position not respected
import plotly
import plotly.graph_objects as go

df = df_go
value_cols='value'
cat_cols=['level5', 'level4', 'level3', 'level2', 'level1', 'level0']

labelList = []
for catCol in cat_cols:
    labelListTemp =  list(set(df[catCol].values))
    labelList = labelList + labelListTemp
    
for i in range(len(cat_cols)-1):
    if i==0:
        sourceTargetDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
        sourceTargetDf.columns = ['source','target','count']
    else:
        tempDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
        tempDf.columns = ['source','target','count']
        sourceTargetDf = pd.concat([sourceTargetDf,tempDf])
    sourceTargetDf = sourceTargetDf.groupby(['source','target']).agg({'count':'sum'}).reset_index()

#labelList = list(set(sourceTargetDf.source) | (set(sourceTargetDf.target)))
sourceTargetDf['sourceID'] = sourceTargetDf['source'].apply(lambda x: labelList.index(x))
sourceTargetDf['targetID'] = sourceTargetDf['target'].apply(lambda x: labelList.index(x))

x = [] #fix the axis
c = []
for l in labelList:
    if isinstance(l, str):
        col_idx = np.where(df == l)[1][0] #col index
        x.append(col_idx)
        c.append(df.loc[np.where(df == l)[0][0], 'color'])
    else:
        x.append([])
        c.append([])


fig = go.Figure(go.Sankey(
    arrangement = "freeform",
    node = {
        "label": labelList,
        #"x":x,
        "y":x,
        'color': c,
        'pad':10},  # 10 Pixels
    link = {
      'source' : sourceTargetDf['sourceID'],
      'target' : sourceTargetDf['targetID'],
      'value' : sourceTargetDf['count']}
    ))

plotly.offline.plot(fig, validate=False)


#%% False heatmap with the LogC at 4 and 8 hpi for the best GO term at lvl 0
# map the value to a heatmap (bidirectional)
# create 2 boxes by protein and plot them
# save as SVG

#extract value from GO based and found the lvl0 of interest, extract protein list from them
with open('GO'+os.sep+'analysis.json') as file:
    js = json.load(file)
    js = js['overrepresentation']

data = js['group']

interest_lvl0 = ['apoptotic mitochondrial changes',
                 'cristae formation',
                 'regulation of mitochondrial membrane permeability',
                 'negative regulation of mitochondrion organization',
                 'mitochondrial fission',
                 'mitochondrial respiratory chain complex I assembly',
                 'mitochondrial cytochrome c oxidase assembly',
                 'mitochondrial respiratory chain complex III assembly',
                 'proton motive force-driven mitochondrial ATP synthesis']

poi = {} #protein of interest


for d in data:
    sub = d['result']
    if isinstance(sub, dict):
        inp = sub['input_list']
        t = sub['term']
        if t['level'] == 0: #if level 0
            if t['label'] in interest_lvl0: #and in the list
                dd = {'label':t['label'],
                      'id':t['id'],
                      'pValue':inp['pValue'],
                      'enrichment':inp['fold_enrichment'],
                      'protein':inp['mapped_id_list']['mapped_id']}
                poi[t['label']] = dd
                
    else:
        for s in sub:
            inp = s['input_list']
            t = s['term']
            if t['level'] == 0: #if level 0
                if t['label'] in interest_lvl0: #and in the list
                    dd = {'label':t['label'],
                          'id':t['id'],
                          'pValue':inp['pValue'],
                          'enrichment':inp['fold_enrichment'],
                          'protein':inp['mapped_id_list']['mapped_id']}
                    poi[t['label']] = dd


# make the graph

#need to find the highest (max) and lowest (min) logFC in GRO-Seq for the colormap
log = []
for p in poi:
    for prot in poi[p]['protein']:
        if all_info[prot]['PValue_4h'] < 0.05:
            log.append(all_info[prot]['logFC_4h'])
        if all_info[prot]['PValue_8h'] < 0.05:
            log.append(all_info[prot]['logFC_8h'])

mi, ma = np.percentile(log, [5,95])

cm_value = max(abs(mi), abs(ma))

import matplotlib
from matplotlib.patches import Rectangle
import math

cmap = matplotlib.colormaps.get_cmap('Spectral_r') #from 0 to 1!
norm = matplotlib.colors.Normalize(vmin=-cm_value, vmax=cm_value) #normalize


for p in poi:
    names = sorted(poi[p]['protein']) #sort alphabeticelly
    prot_logfc_4h = []
    prot_signi_4h = []
    prot_logfc_8h = []
    prot_signi_8h = []
    for prot in names:
        prot_signi_4h.append(all_info[prot]['PValue_4h'])
        prot_logfc_4h.append(all_info[prot]['logFC_4h'])
        prot_signi_8h.append(all_info[prot]['PValue_8h'])
        prot_logfc_8h.append(all_info[prot]['logFC_8h'])      
    
    n = len(names)
    nrows = math.ceil(n / 1)
    
    cell_width = 50
    cell_height = 22
    swatch_width = 48
    margin = 12
    
    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + 2 * margin
    dpi = 72
    
    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-margin)/height)
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    ax.set_title(p+' ('+poi[p]['id']+')')
    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height
    
        swatch_start_x = cell_width * col + 75
        text_pos_x = cell_width * col #+ swatch_width + 7
    
        ax.text(text_pos_x, y, name, fontsize=10,
                horizontalalignment='left',
                verticalalignment='center')
        
        rec = Rectangle(xy=(swatch_start_x, y-9), width=swatch_width,
                  height=18, facecolor=cmap(norm(prot_logfc_4h[i])), edgecolor='0') 
        if prot_signi_4h[i] > 0.05:
            rec.set_hatch('++')
        if prot_logfc_4h[i] > -1 and prot_logfc_4h[i] < 1:
            rec.set_hatch('xx')
        ax.add_patch(rec)
        
        rec = Rectangle(xy=(swatch_start_x+50, y-9), width=swatch_width,
                  height=18, facecolor=cmap(norm(prot_logfc_8h[i])), edgecolor='0')
        if prot_signi_8h[i] > 0.05:
            rec.set_hatch('++')
        if prot_logfc_8h[i] > -1 and prot_logfc_8h[i] < 1:
            rec.set_hatch('xx')
        ax.add_patch(rec)    
    
    plt.savefig('GO'+os.sep+p+'.svg')
    plt.show()


#cbar
fig,ax = plt.subplots(figsize=(3,0.3))
cbar = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                                        orientation='horizontal')
plt.xticks(np.arange(round(-cm_value, 0), round(cm_value, 0)+1, 1))
plt.savefig('GO'+os.sep+"colorbar.svg", bbox_inches='tight')
plt.show

# save the GO data used for the heatmap 
table = [['','Upregulated', 'Downregulated', 'Total'],
         ['4 hpi',len(up4), len(down4), len(up4)+len(down4)],
         ['8 hpi', len(up8), len(down8), len(up8)+len(down8)],
         ['Total', unique_up, unique_down, unique_up+unique_down]]

table = [['Name','GO id', 'pValue', 'fold enrichment', 'protein list']]
for p in poi:
    pp = poi[p]
    table.append([pp['label'], pp['id'], pp['pValue'],
                  pp['enrichment'],pp['protein']])

for row in table:
    print('| {:1} | {:1} | {:1} | {:1} |'.format(*row))

import csv

with open("GO"+os.sep+"GO_table.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(table)