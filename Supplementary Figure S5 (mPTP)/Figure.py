# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 14:58:41 2024

@author: LECLERC
"""

import numpy as np
import tifffile
import os
from scipy import stats, ndimage
import seaborn as sns
import pandas as pd
from statannotations.Annotator import Annotator
import matplotlib.pyplot as plt
import matplotlib as mpl
from skimage.transform import resize
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
mpl.rcParams['mathtext.default'] = 'regular'

import utils as u

color = u.get_color(4, 'original')
sns.set_palette(sns.color_palette(color, len(color)))

def make_swarmplot(data, ax, dotsize):
    sns.boxplot(data=data, ax=ax, 
                      showcaps=False, boxprops={'facecolor':'None'},
                      showfliers=False, whiskerprops={'linewidth':0},
                      medianprops={'linewidth':0}, showmeans=True,
                      meanline=True, meanprops={'color':'black'})
    sns.swarmplot(data=data, size=dotsize, ax=ax, zorder=0)
    annotator = Annotator(ax, pairs, data=data)
    annotator.configure(test=test, verbose=True).apply_and_annotate()


def add_scalebar(pxl_size, ax, scale_bar_size=10):
    return AnchoredSizeBar(ax.transData,
                           scale_bar_size/pxl_size, '', 'lower right', 
                           borderpad = 5,
                           color='white',
                           frameon=False,
                           size_vertical=3,
                           fontproperties={'size':0})


def pick_cell(data, idx, interval_incr): #grab an average cell
    average_cell_value = data.mean()
    interval = 0.0
    while True:
        interval += interval_incr #smaller lead to more precise cell selection
        min_cell = average_cell_value - average_cell_value * interval
        max_cell = average_cell_value + average_cell_value * interval
        #select average ells
        sub = data[(data.nucleus_area > min_cell.nucleus_area) & (data.nucleus_area < max_cell.nucleus_area)]
        sub = sub[(sub.cyto_area > min_cell.cyto_area) & (sub.cyto_area < max_cell.cyto_area)]
        sub = sub[(sub.Infection_intensity > min_cell.Infection_intensity) & (sub.Infection_intensity < max_cell.Infection_intensity)]
        sub = sub[(sub.AVG_mPTP_IMAGE3_intensity > min_cell.AVG_mPTP_IMAGE3_intensity) & (sub.AVG_mPTP_IMAGE3_intensity < max_cell.AVG_mPTP_IMAGE3_intensity)]       
        sub.sort_values('Infection_intensity', ascending=False, inplace=True)
        if len(sub) >= 1: #1 or more cell!
            cell = sub.iloc[idx].squeeze() #select one!
            break
        if interval >= 1: #security
            break
        
    return cell

def make_cell_img(cell, img, img2, ax, bbox_x_cor, bbox_y_cor):
    pan = 5 #add some pixels on the edge
    bbox = [[max(0, np.min(cell.cyto_coords[:,0])-pan+bbox_x_cor), min(np.max(cell.cyto_coords[:,0])+pan+bbox_x_cor, img.shape[1])],
            [max(0, np.min(cell.cyto_coords[:,1])-pan+bbox_y_cor), min(np.max(cell.cyto_coords[:,1])+pan+bbox_y_cor, img.shape[2])]]
    sub_ICP4_mPTP = img[0, bbox[0][0]:bbox[0][1],bbox[1][0]:bbox[1][1]]
    sub_ICP4_mPTP = resize(sub_ICP4_mPTP, np.array(sub_ICP4_mPTP.shape)*resize_f, preserve_range=True)
    bbox2 = [[max(0, np.min(cell.cyto_coords[:,0])-pan), min(np.max(cell.cyto_coords[:,0])+pan, img.shape[1])],
            [max(0, np.min(cell.cyto_coords[:,1])-pan), min(np.max(cell.cyto_coords[:,1])+pan, img.shape[2])]]
    sub_ICP4 = img2[0, bbox2[0][0]:bbox2[0][1],bbox2[1][0]:bbox2[1][1]]
    sub_ICP4 = resize(sub_ICP4, np.array(sub_ICP4.shape)*resize_f, preserve_range=True)
    #cyto label
    cyto_b = np.zeros(img.shape[1:3], dtype=bool)
    for pxl in cell.cyto_coords:
        cyto_b[pxl[0]+bbox_x_cor,pxl[1]+bbox_y_cor] = True
    cyto_b = cyto_b[bbox[0][0]:bbox[0][1],bbox[1][0]:bbox[1][1]]
    cyto_b = resize(cyto_b, np.array(cyto_b.shape)*resize_f, preserve_range=True)
    cyto_b = cyto_b ^ ndimage.binary_erosion(cyto_b, iterations=1)
    #nucleus label
    nucleus_b = np.zeros(img.shape[1:3], dtype=bool)
    for pxl in cell.nucleus_coords:
        nucleus_b[pxl[0]+bbox_x_cor,pxl[1]+bbox_y_cor] = True
    nucleus_b = nucleus_b[bbox[0][0]:bbox[0][1],bbox[1][0]:bbox[1][1]]
    nucleus_b = resize(nucleus_b, np.array(nucleus_b.shape)*resize_f, preserve_range=True)
    nucleus_edge = nucleus_b ^ ndimage.binary_erosion(nucleus_b, iterations=1)
    #removing ICP4 signal
    if not np.array_equal(sub_ICP4, sub_ICP4_mPTP): #infected
        sub_mPTP = np.where(nucleus_b,
                            sub_ICP4_mPTP,
                            np.clip(sub_ICP4_mPTP.astype(float)-sub_ICP4.astype(float), 1, 100))
    else: #NI
        sub_mPTP = sub_ICP4_mPTP
    #remove the nuclei from mPTP
    sub_mPTP = np.where(nucleus_b, 0, sub_mPTP)
    #mitochondria
    mito = img[2, bbox[0][0]:bbox[0][1],bbox[1][0]:bbox[1][1]]
    mito = resize(mito, np.array(mito.shape)*resize_f, preserve_range=True)
    #color merger
    ax.imshow(u.img_color_merger(
                                 #b=nucleus_edge*255,
                                 gr=cyto_b*255,
                                 c=sub_mPTP*2.5,
                                 #y=mito
                                 ),
              interpolation=None
              )
    ax.add_artist(add_scalebar(pxl_size, ax))   
    
def extract_data(data, key):
    #extract the column
    to_df = []
    for x, d in enumerate(data):
        if x != 0: #NI
            sub = d[d.Infection_intensity > infection_threshold]
            sub = sub[key]
        else:
            sub = d['AVG_mPTP_IMAGE3_intensity']
        sub = np.array(sub)[np.abs(stats.zscore(sub)) < 3]
        to_df.append(sub)
    
    
    to_df = {condition_name_list[idx]:x for idx, x in enumerate(to_df)}
    to_df = pd.DataFrame(dict([(k, pd.Series(v)) for k,v in to_df.items()]))
    to_df = to_df[condition_name_list]
    
    return to_df 


div = 1.5

fig = plt.figure(figsize=(10/div, 15/div))

#images
resize_f = 2.0
pxl_size = 0.3/resize_f #um
#NI
ax = fig.add_subplot(3,2,1)
#load the data
img = tifffile.imread('NI WS'+os.sep+'noninfceted control working solution 12042024.lif - control noninfected with working solution IMAGE 4_Processed001_Merging_001.tif')
data = pd.read_pickle('NI WS'+os.sep+'noninfceted control working solution 12042024.lif - control noninfected with working solution Processed001_Merging_001.pckl')

cell = pick_cell(data, 2, 0.1)
make_cell_img(cell, img, img, ax, 0, 0)
plt.axis('off')


#4hpi
ax = fig.add_subplot(3,2,2)
#load the data
img = tifffile.imread('4hpi'+os.sep+'4 hpi second trial 12022024.lif - 4 hpi IMAGE 3_Processed001_Merging_001.tif')
img2 = tifffile.imread('4hpi'+os.sep+'4 hpi second trial 12022024.lif - 4 hpi IMAGE 2_Processed001_Merging_001.tif')
data = pd.read_pickle('4hpi'+os.sep+'4 hpi second trial 12022024.lif - 4 hpi Processed001_Merging_001.pckl')

cell = pick_cell(data, 6, 0.3)
make_cell_img(cell, img, img2, ax, -18, 0)
plt.axis('off')

#8hpi
ax = fig.add_subplot(3,2,3)
#load the data
img = tifffile.imread('8hpi'+os.sep+'8 hpi 14022024.lif - 8 hpi IMAGE 3_Processed001_Merging_001.tif')
img2 = tifffile.imread('8hpi'+os.sep+'8 hpi 14022024.lif - 8 hpi IMAGE 2_Processed001_Merging_001.tif')
data = pd.read_pickle('8hpi'+os.sep+'8 hpi 14022024.lif - 8 hpi Processed001_Merging_001.pckl')

cell = pick_cell(data, 2, 0.7)
make_cell_img(cell, img, img2, ax, -16, -10)
plt.axis('off')

#8hpi
ax = fig.add_subplot(3,2,4)
#load the data
img = tifffile.imread('12hpi'+os.sep+'12 hpi 14022024.lif - 12 hpi IMAGE 3_Processed001_Merging_001.tif')
img2 = tifffile.imread('12hpi'+os.sep+'12 hpi 14022024.lif - 12 hpi IMAGE 2 _Processed001_Merging_001.tif')
data = pd.read_pickle('12hpi'+os.sep+'12 hpi 14022024.lif - 12 hpi IMAGE 2 _Processed001_Merging_001.pckl')

cell = pick_cell(data, 2, 0.3)
make_cell_img(cell, img, img2, ax, -4, 10)
plt.axis('off')

#graph part
#load the data
data_4h = pd.read_pickle('4hpi'+os.sep+'4 hpi second trial 12022024.lif - 4 hpi Processed001_Merging_001.pckl')
data_8h = pd.read_pickle('8hpi'+os.sep+'8 hpi 14022024.lif - 8 hpi Processed001_Merging_001.pckl')
data_12h = pd.read_pickle('12hpi'+os.sep+'12 hpi 14022024.lif - 12 hpi IMAGE 2 _Processed001_Merging_001.pckl')
data_ni = pd.read_pickle('NI WS'+os.sep+'noninfceted control working solution 12042024.lif - control noninfected with working solution Processed001_Merging_001.pckl')

infection_threshold = data_ni.Infection_intensity.max()
data = [data_ni, data_4h, data_8h, data_12h]
condition_name_list = ['NI', '4 hpi', '8 hpi', '12 hpi']

pairs = [['NI', '4 hpi'],
         ['NI', '8 hpi'],
         ['NI', '12 hpi'],
         ['4 hpi', '8 hpi'],
         ['4 hpi', '12 hpi'],
         ['8 hpi', '12 hpi']
         ]
test = 't-test_ind'
dotsize = 5/div


ax = fig.add_subplot(3,2,5)
d = extract_data(data, 'AVG_mPTP_image_difference')
d.to_csv('Cytoplasmic Calcein.csv')
make_swarmplot(d, ax, dotsize)
ax.set_ylabel('Cytoplasmic Calcein (a.u.)')
ax.set_xticklabels(condition_name_list, rotation=45, ha='center')

ax = fig.add_subplot(3,2,6)
d = extract_data(data, 'Infection_intensity')
d.to_csv('Infection.csv')
make_swarmplot(d, ax, dotsize)
ax.set_ylabel('Infection (a.u.)')
ax.set_xticklabels(condition_name_list, rotation=45, ha='center')

# ax = fig.add_subplot(3,3,7)
# to_df = []
# for x, d in enumerate(data):
#     if x != 0: #NI
#         sub = d[d.Infection_intensity > infection_threshold]
#         sub = sub['cyto_mito_intensity']
#     else:
#         sub = d['cyto_mito_intensity']
#     sub = np.array(sub)[np.abs(stats.zscore(sub)) < 3]
#     to_df.append(sub)
# to_df = {condition_name_list[idx]:x for idx, x in enumerate(to_df)}
# to_df = pd.DataFrame(dict([(k, pd.Series(v)) for k,v in to_df.items()]))
# to_df = to_df[condition_name_list]

# make_swarmplot(to_df, ax, dotsize)
# ax.set_ylabel('Mitochondria intensity (a.u.)')
# ax.set_xticklabels(condition_name_list, rotation=45, ha='center')

plt.tight_layout()
plt.savefig('Figure.svg')
plt.show()