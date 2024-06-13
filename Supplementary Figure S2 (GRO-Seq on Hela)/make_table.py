# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 09:03:23 2023

@author: sleclerc
"""

import gzip
import json
import pandas as pd
import numpy as np
import matplotlib
from matplotlib.patches import Rectangle
import math
import os

#global parameters
fc = 1 #can be 0.585, 1 or 2 for corresponding FC
folder = r'Table//'
r_folder = r'Result//'


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
            sub_data['logFC_12h'] = sel['mock vs. HSV1_12h logFC']
            sub_data['PValue_4h'] = sel['mock vs. HSV1_4h PValue']
            sub_data['PValue_8h'] = sel['mock vs. HSV1_8h PValue']
            sub_data['PValue_12h'] = sel['mock vs. HSV1_12h PValue']
            sub_data['chr'] = sel['chr']
            info[gene] = sub_data
        elif isinstance(sel, pd.DataFrame) and len(sel) > 1: #security, should not happen!
            print(gene+' has multiple entry in the GRO-seq database! '+str(len(sel)))
            sub_data['logFC_4h'] = sel['mock vs. HSV1_4h logFC'].max()
            sub_data['logFC_8h'] = sel['mock vs. HSV1_8h logFC'].max()
            sub_data['logFC_12h'] = sel['mock vs. HSV1_12h logFC'].max()
            sub_data['PValue_4h'] = sel['mock vs. HSV1_4h PValue'].max()
            sub_data['PValue_8h'] = sel['mock vs. HSV1_8h PValue'].max()
            sub_data['PValue_12h'] = sel['mock vs. HSV1_12h PValue'].max()
            sub_data['chr'] = sel['chr']
            info[gene] = sub_data

    return info



with gzip.open('proteinatlas.json.gz', 'r') as fin:
    data = json.loads(fin.read().decode('utf-8'))

# get the expression data
groseq_data = pd.read_excel('getDiffExpression_JKL_Hela112018_HSV1_0_12h_rpkm_added_condensed_EXCEL_sorting.xlsx')

#no Gene name column. Instead, called Annotation/Divergence
groseq_data['Gene name'] = groseq_data['Annotation/Divergence'].str.split('|', expand=True, n=1)[0]


up_4h = groseq_data[(groseq_data['mock vs. HSV1_4h logFC']>fc)&(groseq_data['mock vs. HSV1_4h FDR']<0.05)&(groseq_data['MAX_RPKM']>0.5)]
down_4h = groseq_data[(groseq_data['mock vs. HSV1_4h logFC']<-fc)&(groseq_data['mock vs. HSV1_4h FDR']<0.05)&(groseq_data['MAX_RPKM']>0.5)]

up_8h = groseq_data[(groseq_data['mock vs. HSV1_8h logFC']>fc)&(groseq_data['mock vs. HSV1_8h FDR']<0.05)&(groseq_data['MAX_RPKM']>0.5)]
down_8h = groseq_data[(groseq_data['mock vs. HSV1_8h logFC']<-fc)&(groseq_data['mock vs. HSV1_8h FDR']<0.05)&(groseq_data['MAX_RPKM']>0.5)]

up_12h = groseq_data[(groseq_data['mock vs. HSV1_12h logFC']>fc)&(groseq_data['mock vs. HSV1_8h FDR']<0.05)&(groseq_data['MAX_RPKM']>0.5)]
down_12h = groseq_data[(groseq_data['mock vs. HSV1_12h logFC']<-fc)&(groseq_data['mock vs. HSV1_8h FDR']<0.05)&(groseq_data['MAX_RPKM']>0.5)]


#Data to extract
gene_name = [x['Gene'] for x in data] #get the list of gene
to_extract = ('Gene', 'Molecular function', 'Biological process', 'Subcellular location', 'Uniprot')

#up and down regulated protein information
up_4h_info = grab_merge_data(up_4h, gene_name, to_extract)
down_4h_info = grab_merge_data(down_4h, gene_name, to_extract)
up_8h_info = grab_merge_data(up_8h, gene_name, to_extract)
down_8h_info = grab_merge_data(down_8h, gene_name, to_extract)
up_12h_info = grab_merge_data(up_12h, gene_name, to_extract)
down_12h_info = grab_merge_data(down_12h, gene_name, to_extract)

#All regulated protein information (fusion of up and down)
twelve_h_info = {**up_12h_info, **down_12h_info}
eight_h_info = {**up_8h_info, **down_8h_info}
for_h_info = {**up_4h_info, **down_4h_info}
all_info = {**twelve_h_info, **eight_h_info, **for_h_info} #fused

full_df = pd.DataFrame.from_dict(all_info)
full_df = full_df.T
full_df.to_csv(folder+'full_table.csv')

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

up12 = full_df[(full_df['logFC_12h']>=1) &
              (full_df['PValue_12h']<=0.05) &
              (full_df['Subcellular location'].astype(str).str.contains(fav))]

down12 = full_df[(full_df['logFC_12h']<=-1) &
                (full_df['PValue_12h']<=0.05) &
                (full_df['Subcellular location'].astype(str).str.contains(fav))]


unique_up = len(set(list(up4['Gene'])+list(up8['Gene'])+list(up12['Gene'])))
unique_down = len(set(list(down4['Gene'])+list(down8['Gene'])+list(down12['Gene'])))


table = [['','Upregulated', 'Downregulated', 'Total'],
         ['4 hpi',len(up4), len(down4), len(up4)+len(down4)],
         ['8 hpi', len(up8), len(down8), len(up8)+len(down8)],
         ['12 hpi', len(up12), len(down12), len(up12)+len(down12)],
         ['Total', unique_up, unique_down, unique_up+unique_down]]
for row in table:
    print('| {:1} | {:1} | {:1} | {:1} |'.format(*row))

unique = list(set(list(up4['Gene'])+list(up8['Gene'])+list(up12['Gene']))) + list(set(list(down4['Gene'])+list(down8['Gene'])+list(down12['Gene'])))

fav_df = full_df[full_df.Gene.isin(unique)]
fav_df.to_csv(folder+fav+'_table.csv')

#%%Venn diagram full cell + mito focus
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

#full cell
plt.figure(figsize=(10,10))
v = venn2([for_h_info.keys(),eight_h_info.keys(),
           #twelve_h_info.keys()
           ], ['4 hpi', '8 hpi',
               #'12 hpi'
               ])
plt.title('Compare genes expression at LogFC='+str(fc))
plt.savefig(r_folder+'Venn Comparison full data.png', dpi=600)
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
twelve_h_mito_info = extract_loc_dict(twelve_h_info, 'Mitochondria')

plt.figure(figsize=(10,10))
v = venn2([for_h_mito_info.keys(),
           eight_h_mito_info.keys(), 
           #twelve_h_mito_info.keys()
           ],
          ['4 hpi', '8 hpi',
           #' 12hpi'
           ])
plt.title('Compare genes expression in mitochondria at LogFC='+str(fc))
plt.savefig(r_folder+'Venn Comparison mitochondria focus.png', dpi=600)
plt.savefig(r_folder+'Venn Comparison mitochondria focus.svg')
plt.show()


#%% Dual direction bar plot for mito only
# the idea is to categorize the data by their biological process
# Then up and down regulated count (in corresponding direction)
# Then split the 4 and 8 hpi (color)

graph_value = {}
key = 'Biological process'
d = {**for_h_mito_info, **eight_h_mito_info, **twelve_h_mito_info}

for item in d:
    t = d[item][key]
    if isinstance(t, list):
        for i in t:
            if i not in graph_value:
                dd = {'up_4h':0,'down_4h':0,
                      'up_8h':0,'down_8h':0,
                      'up_12h':0, 'down_12h':0,
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
    
            if d[item]['logFC_12h'] > 1:
                graph_value[i]['up_12h'] += 1
            elif d[item]['logFC_12h'] < -1:
                graph_value[i]['down_12h'] -= 1
            if d[item]['logFC_12h'] > 1 or d[item]['logFC_12h'] < -1:
                graph_value[i]['count'] +=1   


clean_graph_value = {}
for item in graph_value:
    if graph_value[item]['count'] >= 5:
        clean_graph_value[item] = graph_value[item]

df = pd.DataFrame.from_dict(clean_graph_value)
df = df.T
df = df.sort_index()
df.to_csv('Result'+os.sep+'Gro-Seq protein biological process in Mitochondria.csv')

name = list(df.index)
x = np.arange(len(name))  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots(layout='constrained', figsize=(10,10))
ax.grid(True, linestyle='--', color='#D3D3D3', zorder=0)

ax.bar(x, np.array(df['up_4h']), width, label='4hpi upregulated', color='#A4D8D8', zorder=3)
ax.bar(x, np.array(df['down_4h']), width, label='4hpi downregulated', color='#F1B3A9', zorder=3)   
ax.bar(x + width, np.array(df['up_8h']), width, label='8hpi upregulated', color='#047c99', zorder=3)
ax.bar(x + width, np.array(df['down_8h']), width, label='8hpi downregulated', color='#C3B1E1', zorder=3)  
# ax.bar(x + width*2, np.array(df['up_12h']), width, label='12hpi upregulated', color='#011f4b', zorder=3)
# ax.bar(x + width*2, np.array(df['down_12h']), width, label='12hpi downregulated', color='#d291ff', zorder=3)  
ax.axhline(y=0,  color='k', zorder=5, linewidth=0.5)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Count')
ax.set_title('Gro-Seq protein biological process in Mitochondria')
ax.set_xticks(x + width, name, rotation=45, ha='right')
ticks =  ax.get_yticks()
ax.set_yticklabels([int(abs(tick)) for tick in ticks])
ax.legend(loc='upper left', ncols=2)
ax.grid(True)
plt.savefig(r_folder+'Gro-Seq protein biological process in Mitochondria.png', dpi=600)
plt.savefig(r_folder+'Gro-Seq protein biological process in Mitochondria.svg')
plt.show()



#%% Same on the full dataset
graph_value = {}
key = 'Biological process'
d = {**for_h_info, **eight_h_info, **twelve_h_info}

for item in d:
    t = d[item][key]
    if isinstance(t, list):
        for i in t:
            if i not in graph_value:
                dd = {'up_4h':0,'down_4h':0,
                      'up_8h':0,'down_8h':0,
                      'up_12h':0,'down_12h':0,
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

            if d[item]['logFC_12h'] > 1:
                graph_value[i]['up_12h'] += 1
            elif d[item]['logFC_12h'] < -1:
                graph_value[i]['down_12h'] -= 1
            if d[item]['logFC_12h'] > 1 or d[item]['logFC_12h'] < -1:
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
# ax.bar(x + width*2, np.array(df['up_12h']), width, label='12hpi upregulated', color='#011f4b', zorder=3)
# ax.bar(x + width*2, np.array(df['down_12h']), width, label='12hpi downregulated', color='#d291ff', zorder=3)  
ax.axhline(y=0,  color='k', zorder=5, linewidth=0.5)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Count')
ax.set_title('Gro-Seq protein biological process')
ax.set_xticks(x + width, name, rotation=45, ha='right')
ticks =  ax.get_yticks()
ax.set_yticklabels([int(abs(tick)) for tick in ticks])
ax.legend(loc='upper left', ncols=2)

plt.savefig(r_folder+'Gro-Seq protein biological process.png', dpi=600)
plt.show()

#%% GO style TODO! Generate the analysis file
import json
import os

#extract value from GO based and classify based on the protein name
with open('analysis.json') as file:
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

d = {**for_h_mito_info, **eight_h_mito_info, **twelve_h_mito_info}

        
for item in protein_GO:
    t = protein_GO[item]['label']
    if isinstance(t, list):
        for i in t:
            if i not in graph_value:
                dd = {'up_4h':0,'down_4h':0,
                      'up_8h':0,'down_8h':0,
                      'up_12h':0,'down_12h':0,
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

            if d[item]['logFC_12h'] > 1:
                graph_value[i]['up_12h'] += 1
            elif d[item]['logFC_12h'] < -1:
                graph_value[i]['down_12h'] -= 1
            if d[item]['logFC_12h'] > 1 or d[item]['logFC_12h'] < -1:
                graph_value[i]['count'] +=1  

clean_graph_value = {}
for item in graph_value:
    if graph_value[item]['count'] >= 0:
        clean_graph_value[item] = graph_value[item]

df = pd.DataFrame.from_dict(clean_graph_value)
df = df.T
df = df.sort_index()


name = list(df.index)
name_GO = name
x = np.arange(len(name))  # the label locations
width = 0.25  # the width of the bars

# c = u.get_color(4, 'bar')

fig, ax = plt.subplots(layout='constrained', figsize=(10,10))
ax.grid(True, linestyle='--', color='#D3D3D3', zorder=0)

ax.bar(x, np.array(df['up_4h']), width, label='4hpi upregulated', color='#A4D8D8', zorder=3)
ax.bar(x, np.array(df['down_4h']), width, label='4hpi downregulated', color='#F1B3A9', zorder=3)   
ax.bar(x + width, np.array(df['up_8h']), width, label='8hpi upregulated', color='#047c99', zorder=3)
ax.bar(x + width, np.array(df['down_8h']), width, label='8hpi downregulated', color='#C3B1E1', zorder=3)  
# ax.bar(x + width*2, np.array(df['up_12h']), width, label='12hpi upregulated', color='#011f4b', zorder=3)
# ax.bar(x + width*2, np.array(df['down_12h']), width, label='12hpi downregulated', color='#d291ff', zorder=3)  
ax.axhline(y=0,  color='k', zorder=5, linewidth=0.5)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Count')
ax.set_title('Gro-Seq protein Mitochondria function (GO) at level '+str(level))
ax.set_xticks(x + width, name, rotation=45, ha='right')
ticks =  ax.get_yticks()
ax.set_yticklabels([int(abs(tick)) for tick in ticks])
ax.legend(loc='upper left', ncols=2)

plt.savefig(r_folder+'Gro-Seq protein mitochondria function from GO at level='+str(level)+'.png', dpi=600)
plt.savefig(r_folder+'Gro-Seq protein mitochondria function from GO at level='+str(level)+'.svg')
plt.show()


#%% False heatmap with the LogC at 4 and 8 hpi for the best GO term at lvl 0
# map the value to a heatmap (bidirectional)
# create 2 boxes by protein and plot them
# save as SVG

#extract value from GO based and found the lvl0 of interest, extract protein list from them
with open('GO'+os.sep+'analysis.json') as file:
    js = json.load(file)
    js = js['overrepresentation']

data = js['group']

interest_lvl0 = name_GO

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
        if all_info[prot]['PValue_12h'] < 0.05:
            log.append(all_info[prot]['logFC_12h'])
            
mi, ma = np.percentile(log, [5,95])

cm_value = max(abs(mi), abs(ma))



cmap = matplotlib.colormaps.get_cmap('Spectral_r') #from 0 to 1!
norm = matplotlib.colors.Normalize(vmin=-cm_value, vmax=cm_value) #normalize


for p in poi:
    names = sorted(poi[p]['protein']) #sort alphabeticelly
    prot_logfc_4h = []
    prot_signi_4h = []
    prot_logfc_8h = []
    prot_signi_8h = []
    prot_logfc_12h = []
    prot_signi_12h = []
    for prot in names:
        prot_signi_4h.append(all_info[prot]['PValue_4h'])
        prot_logfc_4h.append(all_info[prot]['logFC_4h'])
        prot_signi_8h.append(all_info[prot]['PValue_8h'])
        prot_logfc_8h.append(all_info[prot]['logFC_8h'])      
        prot_signi_12h.append(all_info[prot]['PValue_12h'])
        prot_logfc_12h.append(all_info[prot]['logFC_12h'])  
        
    n = len(names)
    nrows = math.ceil(n / 1)
    
    cell_width = 60
    cell_height = 22
    swatch_width = 48
    margin = 10
    
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
        
        rec = Rectangle(xy=(swatch_start_x+100, y-9), width=swatch_width,
                  height=18, facecolor=cmap(norm(prot_logfc_12h[i])), edgecolor='0')
        if prot_signi_12h[i] > 0.05:
            rec.set_hatch('++')
        if prot_logfc_12h[i] > -1 and prot_logfc_12h[i] < 1:
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
