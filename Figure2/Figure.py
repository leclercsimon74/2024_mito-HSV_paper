# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 13:43:57 2023

@author: sleclerc
"""

import numpy as np
import tifffile
import pickle
from skimage import exposure, filters
import os
from scipy import stats
import seaborn as sns
import pandas as pd
from statannotations.Annotator import Annotator
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
mpl.rcParams['mathtext.default'] = 'regular'

import utils as u

color = u.get_color(4, 'original')
sns.set_palette(sns.color_palette(color, len(color)))

col = ['Parameters','Condition','nobs', 'minmax', 'mean', 'variance', 'skewness', 'kurtosis']
stats_df = []


with open(r"ALBA"+os.sep+"ALBA_data.pckl", 'rb') as f:
    alba_data = pickle.load(f)
    
with open(r"Berkeley"+os.sep+"Berkeley_data.pckl", 'rb') as f:
    berkeley_data = pickle.load(f)

def get_data_as_df(d, name):
    data = d[name]
    to_df = []
    i = 0
    for key in data.keys():
        print(key)
        if key == 'NI' and i == 0:
            NI = data[key]
            NI = np.array(NI)[~np.isnan(NI)]
            NI = np.array(NI)[np.abs(stats.zscore(NI)) < 3]
            to_df.append(NI)
            stats_df.append([name, 'Noninfected']+list(stats.describe(NI)))
            i = 1
        if key == 'hpi4' and i == 1:
            hpi4 = data[key]
            hpi4 = np.array(hpi4)[~np.isnan(hpi4)]
            hpi4 = np.array(hpi4)[np.abs(stats.zscore(hpi4)) < 3]
            to_df.append(hpi4)
            stats_df.append([name, '4 hpi']+list(stats.describe(hpi4)))
            i = 2
        if key == 'hpi8':
            if i == 1 or i == 2:
                hpi8 = data[key]
                hpi8 = np.array(hpi8)[~np.isnan(hpi8)]
                hpi8 = np.array(hpi8)[np.abs(stats.zscore(hpi8)) < 3]
                to_df.append(hpi8)
                stats_df.append([name, '8 hpi']+list(stats.describe(hpi8)))
                i = 2
        if key == 'hpi12' and i == 2:
            hpi12 = data[key]
            hpi12 = np.array(hpi12)[~np.isnan(hpi12)]
            hpi12 = np.array(hpi12)[np.abs(stats.zscore(hpi12)) < 3]
            to_df.append(hpi12)
            stats_df.append([name, '12 hpi']+list(stats.describe(hpi12)))
              
    to_df = {condition_name_list[idx]:x for idx, x in enumerate(to_df)}
    to_df = pd.DataFrame(dict([(k, pd.Series(v)) for k,v in to_df.items()]))
    to_df = to_df[condition_name_list]
    
    
    
    return to_df

def make_swarmplot(data, ax, dotsize):
    sns.boxplot(data=data, ax=ax, 
                      showcaps=False, boxprops={'facecolor':'None'},
                      showfliers=False, whiskerprops={'linewidth':0},
                      medianprops={'linewidth':0}, showmeans=True,
                      meanline=True, meanprops={'color':'black'})
    sns.swarmplot(data=data, size=dotsize, ax=ax, zorder=0)
    annotator = Annotator(ax, pairs, data=data)
    annotator.configure(test=test, verbose=True).apply_and_annotate()
  
    
  
    
  





#%% GRAPHS

test = 't-test_ind'    
test = 'Brunner-Munzel' #may be more adapted? Just change on the mitochondria volume, NI to 4 hpi, from **** to ** (making sense)

div = 2
with_title = False
fig = plt.figure(figsize=(32/div, 5/div))

#Berkeley
condition_name_list = ['Noninfected','4 hpi', '8 hpi']
pairs = [['Noninfected', '4 hpi'],
         ['Noninfected', '8 hpi'],
         ['4 hpi', '8 hpi']]

dotsize = 5
#mito number
ax = fig.add_subplot(1,6,1)
data = get_data_as_df(berkeley_data, 'mito_number')
data.to_csv('csv data'+os.sep+'Mito_number.csv')
make_swarmplot(data, ax, dotsize)
if with_title: ax.set_title('Mitochondria number by cell')
ax.set_ylabel('Mitochondria number')
ax.set_xticklabels(['NI', '4 hpi', '8 hpi'])
#LAC
ax = fig.add_subplot(1,6,2)
dotsize = 1
data = get_data_as_df(berkeley_data, 'LAC')
data.to_csv('csv data'+os.sep+'LAC.csv')
make_swarmplot(data, ax, dotsize)
if with_title:ax.set_title('Normalized LAC value of mitochondria')
ax.set_ylabel('Normalized LAC')
ax.set_xticklabels(['NI', '4 hpi', '8 hpi'])
#Volume
ax = fig.add_subplot(1,6,3)
ax.set_yscale('log')
data = get_data_as_df(berkeley_data, 'volume')
data.to_csv('csv data'+os.sep+'volume.csv')
make_swarmplot(data, ax, dotsize)
if with_title:ax.set_title('Mitochondria volume')
ax.set_ylabel('Volume (${\mu m^3}$)')
ax.set_xticklabels(['NI', '4 hpi', '8 hpi'])

dotsize = 1.5
#ALBA
condition_name_list = ['Noninfected', '8 hpi', '12 hpi']
pairs = [['Noninfected', '12 hpi'],
         ['Noninfected', '8 hpi'],
         ['8 hpi', '12 hpi']]
#Distance to nucleus
ax = fig.add_subplot(1,6,4)
data = get_data_as_df(alba_data, 'distance to nucleus')
data.to_csv('csv data'+os.sep+'distance to nucleus.csv')
make_swarmplot(data, ax, dotsize)
if with_title:ax.set_title('Mitochondria distance to nucleus')
ax.set_ylabel('Distance (${\mu m}$)')
ax.set_xticklabels(['NI', '8 hpi', '12 hpi'])
#Diameter
ax = fig.add_subplot(1,6,5)
data = get_data_as_df(alba_data, 'radius')
data = data * 2
data.to_csv('csv data'+os.sep+'diameter.csv')
make_swarmplot(data, ax, dotsize)
if with_title:ax.set_title('Diameter of mitochondria')
ax.set_ylabel('Average diameter (${\mu m}$)')
ax.set_xticklabels(['NI', '8 hpi', '12 hpi'])
#Roughness
ax = fig.add_subplot(1,6,6)
data = get_data_as_df(alba_data, 'SA/V')
data.to_csv('csv data'+os.sep+'roughness.csv')
make_swarmplot(data, ax, dotsize)
if with_title:ax.set_title('Mitochondria roughness')
ax.set_ylabel('Roughness (${\mu m}$)')
ax.set_xticklabels(['NI', '8 hpi', '12 hpi'])


#Finish
plt.tight_layout()
plt.savefig('mito analysis.png', dpi=600)
plt.savefig('mito analysis.svg')
plt.show()

#%% stats
col = ['Parameters','Condition','nobs', 'minmax', 'mean', 'variance', 'skewness', 'kurtosis']
stats_df = pd.DataFrame(stats_df, columns=col)
stats_df.to_pickle('stats.pckl')
stats_df.to_csv('stats.csv')
#%% Images
def add_scalebar(pxl_size, ax, scale_bar_size=2000):
    return AnchoredSizeBar(ax.transData,
                           scale_bar_size/pxl_size, '', 'lower right', 
                           borderpad = 5,
                           color='black',
                           frameon=False,
                           size_vertical=5,
                           fontproperties={'size':0})

#Berkeley
paths = {'NI':'Berkeley'+os.sep+'Berkely-NI_slice479-1.tif',
         '4 hpi':'Berkeley'+os.sep+'Berkely-4hpi_slice188-1.tif',
         '8 hpi':'Berkeley'+os.sep+'Berkely_8hpi_slice295-1.tif'}

resolution = {'NI':25.6,
              '4 hpi':32.7,
              '8 hpi':32.7}


fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1,3,1)
LAC_cor = 33 #LAC value correction factor
im = exposure.rescale_intensity(tifffile.imread(paths['NI'])*LAC_cor,
                                out_range=(0.0018*LAC_cor, 0.0139*LAC_cor))
im = ax.imshow(im, cmap='Greys', interpolation=None)
ax.add_artist(add_scalebar(resolution['NI'], ax))
plt.axis('off')

ax = fig.add_subplot(1,3,2)
im = exposure.rescale_intensity(tifffile.imread(paths['4 hpi'])*LAC_cor,
                                out_range=(0.0018*LAC_cor, 0.0139*LAC_cor))
im = ax.imshow(im, cmap='Greys', interpolation=None)
ax.add_artist(add_scalebar(resolution['4 hpi'], ax))
plt.axis('off')

ax = fig.add_subplot(1,3,3)
im = exposure.rescale_intensity(tifffile.imread(paths['8 hpi'])*LAC_cor,
                                out_range=(0.0018*LAC_cor, 0.0139*LAC_cor))
im2 = ax.imshow(im, cmap='Greys', interpolation=None)
ax.add_artist(add_scalebar(resolution['8 hpi'], ax))
plt.axis('off')


#plt.colorbar(im, ax=ax, location='top', ticks=[0.003, 0.006, 0.009, 0.012])
#plt.colorbar(im2, ax=ax, location='bottom', ticks=[0.1, 0.2, 0.3, 0.4])

plt.tight_layout()
plt.savefig('Berkeley_images.png', dpi=600)
plt.savefig('Berkeley_images.svg')
plt.show()

# draw a new figure and replot the colorbar there
fig,ax = plt.subplots(figsize=(3,5))
plt.colorbar(im2,ax=ax,location='bottom', ticks=[0.1, 0.2, 0.3, 0.4])
ax.remove()
plt.savefig('Berkeley_cbar.svg')

#%% ALBA

paths = {'NI':'ALBA'+os.sep+'ALBA NI-slice404-1.tif',
         '8 hpi':'ALBA'+os.sep+'ALBA 8hpi-slice340-1.tif',
         '12 hpi':'ALBA'+os.sep+'ALBA 12hpi-slice340-1.tif'}
resolution = 13

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1,3,1)

im = exposure.rescale_intensity(tifffile.imread(paths['NI']),
                                in_range=(0, 0.0025))
ax.imshow(im, cmap='Greys', interpolation=None)
ax.add_artist(add_scalebar(resolution, ax))
plt.axis('off')

ax = fig.add_subplot(1,3,2)
im = exposure.rescale_intensity(tifffile.imread(paths['8 hpi']),
                                in_range=(-2.14e-4, 1.25e-3))
ax.imshow(im, cmap='Greys_r', interpolation=None)
ax.add_artist(add_scalebar(resolution, ax))
plt.axis('off')

ax = fig.add_subplot(1,3,3)
im = exposure.rescale_intensity(tifffile.imread(paths['12 hpi']),
                                in_range=(-2.14e-4, 1.25e-3))
ax.imshow(im, cmap='Greys_r', interpolation=None)
ax.add_artist(add_scalebar(resolution, ax))
plt.axis('off')

plt.tight_layout()
plt.savefig('ALBA_images.png', dpi=600)
plt.savefig('ALBA_images.svg')
plt.show()