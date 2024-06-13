# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:06:50 2023

@author: sleclerc
"""

import pandas as pd
import seaborn as sns
from statannotations.Annotator import Annotator
import matplotlib.pyplot as plt
import matplotlib as mpl
import tifffile
import matplotlib
from matplotlib import cm
import numpy as np
from scipy import ndimage
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
mpl.rcParams['mathtext.default'] = 'regular'
from skimage.exposure import rescale_intensity


import utils as u

color = u.get_color(3, 'original')
sns.set_palette(sns.color_palette(color, len(color)))

pairs = [['Noninfected', '12 hpi'],
         ['Noninfected', '8 hpi'],
         ['8 hpi', '12 hpi']]
test = 't-test_ind'
#test = 'Brunner-Munzel'
spacing = [30,8,8] #in nm

nucleus_color = np.array([105, 172, 234])
mito_color = np.array([228, 159, 102])
cell_color = np.array([175, 175, 200])


thickness = 10



fig = plt.figure(figsize=(200, 100))

ax = fig.add_subplot(2,2,1)
scalebar = AnchoredSizeBar(ax.transData,
                           2000/spacing[1], '', 'lower right', 
                           borderpad = 5,
                           color='black',
                           frameon=False,
                           size_vertical=20,
                           fontproperties={'size':0})

mito_img = np.where(tifffile.imread("NI_slice_30_Mito_Seg.tif")>0, True, False)
nucl_img = np.where(tifffile.imread("NI_slice_30_Nucleus_Seg.tif")>0, True, False)
cell_img = np.where(tifffile.imread("NI_slice_30_Cell_Seg.tif")>0, True, False)
raw_img = tifffile.imread('NI_slice_30_Raw.tif')


filtered_img = np.where(cell_img, raw_img, False)
filtered_img = rescale_intensity(filtered_img, (0, 255))
ax.imshow(filtered_img, cmap='gray_r')


# nucl_c = nucl_img ^ ndimage.binary_erosion(nucl_img, iterations = thickness)
# nucl_c = np.where(nucl_c, nucl_c, np.nan)
# ax.imshow(nucl_c*255, cmap=matplotlib.colors.ListedColormap([nucleus_color/255, 'none']))


# if thickness <= 10:
#     mito_c = mito_img ^ ndimage.binary_erosion(mito_img, iterations = thickness)
#     mito_c = np.where(mito_c, 1, np.nan)
# else:
#     mito_c = np.where(mito_img, 1, np.nan)
# ax.imshow(mito_c*255, cmap=matplotlib.colors.ListedColormap([mito_color/255, 'none']))

ax.add_artist(scalebar)
plt.axis('off')

ax = fig.add_subplot(2,2,2)
scalebar = AnchoredSizeBar(ax.transData,
                           2000/spacing[1], '', 'lower right', 
                           borderpad = 5,
                           color='black',
                           frameon=False,
                           size_vertical=20,
                           fontproperties={'size':0})

mito_img = np.where(tifffile.imread("8hpi_slice_50_Mito_Seg.tif")>0, False, True)
nucl_img = np.where(tifffile.imread("8hpi_slice_50_Nucleus_Seg.tif")>0, True, False)
cell_img = np.where(tifffile.imread("8hpi_slice_50_Cell_Seg.tif")>0, True, False)
raw_img = tifffile.imread('8hpi_slice_50_Raw.tif')

#raw_img = rescale_intensity(raw_img, (0, 255))
filtered_img = np.where(cell_img, raw_img, False)
filtered_img = rescale_intensity(filtered_img, (0, 255))
ax.imshow(filtered_img, cmap='gray_r')


# nucl_c = nucl_img ^ ndimage.binary_erosion(nucl_img, iterations = thickness)
# nucl_c = np.where(nucl_c, nucl_c, np.nan)
# ax.imshow(nucl_c*255, cmap=matplotlib.colors.ListedColormap([nucleus_color/255, 'none']))


# if thickness <= 10:
#     mito_c = mito_img ^ ndimage.binary_erosion(mito_img, iterations = thickness)
#     mito_c = np.where(mito_c, 1, np.nan)
# else:
#     mito_c = np.where(mito_img, 1, np.nan)
# ax.imshow(mito_c*255, cmap=matplotlib.colors.ListedColormap([mito_color/255, 'none']))

ax.add_artist(scalebar)
plt.axis('off')

ax = fig.add_subplot(2,2,3)
scalebar = AnchoredSizeBar(ax.transData,
                           2000/spacing[1], '', 'lower right', 
                           borderpad = 5,
                           color='black',
                           frameon=False,
                           size_vertical=20,
                           fontproperties={'size':0})

mito_img = np.where(tifffile.imread("12hpi_slice_68_Mito_Seg.tif")>0, True, False)
nucl_img = np.where(tifffile.imread("12hpi_slice_68_Nucleus_Seg.tif")>0, True, False)
cell_img = np.where(tifffile.imread("12hpi_slice_68_Cell_Seg.tif")>0, True, False)
raw_img = tifffile.imread('12hpi_slice_68_Raw.tif')


filtered_img = np.where(cell_img, raw_img, False)
filtered_img = rescale_intensity(filtered_img, (0, 255))
ax.imshow(filtered_img, cmap='gray_r')


# nucl_c = nucl_img ^ ndimage.binary_erosion(nucl_img, iterations = thickness)
# nucl_c = np.where(nucl_c, nucl_c, np.nan)
# ax.imshow(nucl_c*255, cmap=matplotlib.colors.ListedColormap([nucleus_color/255, 'none']))


# if thickness <= 10:
#     mito_c = mito_img ^ ndimage.binary_erosion(mito_img, iterations = thickness)
#     mito_c = np.where(mito_c, 1, np.nan)
# else:
#     mito_c = np.where(mito_img, 1, np.nan)
# ax.imshow(mito_c*255, cmap=matplotlib.colors.ListedColormap([mito_color/255, 'none']))

ax.add_artist(scalebar)
plt.axis('off')



plt.tight_layout()
plt.savefig('images.svg')
plt.show()


#%%
fig = plt.figure(figsize=(3, 3))
ax1 = fig.add_subplot(1,2,1)
mito_l_mean = pd.read_csv('length.csv')
dotsize = 3
#ax1.set_title('')
sns.boxplot(data=mito_l_mean, ax=ax1, 
                  showcaps=False, boxprops={'facecolor':'None'},
                  showfliers=False, whiskerprops={'linewidth':0},
                  medianprops={'linewidth':0}, showmeans=True,
                  meanline=True, meanprops={'color':'black'})
plt.setp(ax1.lines, zorder=2000)
sns.swarmplot(data=mito_l_mean, size=dotsize, ax=ax1, zorder=0)
annotator = Annotator(ax1, pairs, data=mito_l_mean)
annotator.configure(test=test, verbose=True).apply_and_annotate()
ax1.set_ylabel('Length (${\mu m}$)')

ax1 = fig.add_subplot(1,2,2)
mito_l_mean = pd.read_csv('Volume.csv')
dotsize = 3
sns.boxplot(data=mito_l_mean, ax=ax1, 
                  showcaps=False, boxprops={'facecolor':'None'},
                  showfliers=False, whiskerprops={'linewidth':0},
                  medianprops={'linewidth':0}, showmeans=True,
                  meanline=True, meanprops={'color':'black'})
sns.swarmplot(data=mito_l_mean, size=dotsize, ax=ax1, zorder=0)
annotator = Annotator(ax1, pairs, data=mito_l_mean)
annotator.configure(test=test, verbose=True).apply_and_annotate()
ax1.set_ylabel('Volume (${\mu m^3}$)')
ax1.set_yscale('log')

plt.savefig('graph.svg')
plt.show()