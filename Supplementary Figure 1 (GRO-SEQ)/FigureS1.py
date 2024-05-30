# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 10:12:47 2023

@author: sleclerc
"""

import utils as u

import gzip
import json
import pandas as pd
import numpy as np

fc = 1 #can be 0.585, 1 or 2 for corresponding FC
folder = r'Data//'

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

with gzip.open(folder+'proteinatlas.json.gz', 'r') as fin:
    data = json.loads(fin.read().decode('utf-8'))

# get the expression data
groseq_data = pd.read_excel(folder+'getDiffExpression_VERO_ChlSab1.1_noadj_rpkmAdded_GeneInfoAdded_EXCEL2.xlsx')

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

prot_name = '9606.protein.info.v11.5.txt.gz'
prot_name = pd.read_csv(folder+prot_name, sep='\t')

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
fav_df.to_csv(folder+fav+'_table.csv') #Supplementary Table 1

#%%
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

for_h_mito_info = extract_loc_dict(for_h_info, fav)
eight_h_mito_info = extract_loc_dict(eight_h_info, fav)

plt.figure(figsize=(5,5))
c = u.get_color(4, 'original')
v = venn2([for_h_mito_info.keys(), eight_h_mito_info.keys()],
          ['4 hpi', '8 hpi'])

v.get_patch_by_id('01').set_color(c[0])
v.get_patch_by_id('10').set_color(c[1])
v.get_patch_by_id('11').set_color(c[2])

#plt.title('Compare genes expression in '+fav+' at LogFC='+str(fc))
plt.savefig(folder+'Venn Comparison '+fav+' focus.svg')
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
df.to_csv(folder+'Biological process in mitochondria.csv')

c = u.get_color(4, 'original')

name = list(df.index)
x = np.arange(len(name))  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots(layout='constrained', figsize=(7,7))
ax.grid(True, linestyle='--', color='#D3D3D3', zorder=0)

ax.bar(x, np.array(df['up_4h']), width, label='4hpi upregulated', color='#A4D8D8', zorder=3)
ax.bar(x, np.array(df['down_4h']), width, label='4hpi downregulated', color='#C3B1E1', zorder=3)   
ax.bar(x + width, np.array(df['up_8h']), width, label='8hpi upregulated', color='#047c99', zorder=3)
ax.bar(x + width, np.array(df['down_8h']), width, label='8hpi downregulated', color='#9270C9', zorder=3)  
ax.axhline(y=0,  color='k', zorder=5, linewidth=0.5)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Count')
ax.set_title('Gro-Seq protein biological process in '+fav)
ax.set_xticks(x + width, name, rotation=45, ha='right')
ticks =  ax.get_yticks()
ax.set_yticklabels([int(abs(tick)) for tick in ticks])
ax.legend(loc='lower left', ncols=2)
ax.grid(True)
plt.savefig(folder+'Gro-Seq protein biological process in '+fav+'.svg')
plt.show()

