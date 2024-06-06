# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 14:16:18 2023

@author: sleclerc
"""


import utils as u

import os
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
from skimage import measure
import pandas as pd

pixel_size = np.array([0.36, 0.16, 0.16]) #ZXY in um
mag_f = 12.0
spacing = pixel_size / mag_f

NI_path_list = u.found_img('NI')
INF_path_list = u.found_img('8 hpi')

final_df = pd.DataFrame()

for filename in NI_path_list:
    img = u.get_img(filename) #ZCXY format
    #split
    nucl_img = img[:,0]
    mito_img = img[:,1]
    vapb_img = img[:,2]
    del img
    
    #grab the full nucleus and substract it from mito and vapb signal
    nucl_img = sc.ndimage.gaussian_filter(nucl_img, (1,2,2))
    cleaned_data = u.boolean_remove_zeros(np.ndarray.flatten(nucl_img)) #w/o 0
    thr_min = u.thr_list(cleaned_data, 'otsu')
    nucl_binary = nucl_img > thr_min
    nucl_binary = sc.ndimage.binary_dilation(nucl_binary, iterations=5)
    nucl_binary = sc.ndimage.binary_fill_holes(nucl_binary)
    nucl_binary = sc.ndimage.binary_erosion(nucl_binary, iterations=5)
    u.ortho_view(nucl_binary, 'mean', os.path.basename(filename),ratio=[7,1])
    del nucl_img
    #remove the nucleus
    mito_img = np.where(nucl_binary == False, mito_img, 0)
    vapb_img = np.where(nucl_binary == False, vapb_img, 0)
    
    mito_img = sc.ndimage.gaussian_filter(mito_img, (.5, 1, 1))
    vapb_img = sc.ndimage.gaussian_filter(vapb_img, (.25,.5,.5))

    #distance map between mito and vapb
    #binary transform: all positive signal are positive!
    vapb_thr = u.thr_list(u.boolean_remove_zeros(vapb_img), 'otsu')
    label1 = np.array(np.where(vapb_img>=vapb_thr, True, False))
    #measure foci
    lbl_img = measure.label(label1)
    properties = ['label', 'area', 'centroid_weighted', 'intensity_mean']
    df = pd.DataFrame(measure.regionprops_table(lbl_img,
                                                intensity_image=vapb_img,
                                                properties=properties))
    #get the centroid position
    z_pos = np.array(df['centroid_weighted-0'])
    x_pos = np.array(df['centroid_weighted-1'])
    y_pos = np.array(df['centroid_weighted-2'])

    #assign the centroid as lbl 1
    label1 = np.zeros(vapb_img.shape, dtype=bool)
    centroid = (z_pos.astype(int), x_pos.astype(int), y_pos.astype(int))
    label1[centroid] = True
    
    #because no background, no thresholding
    label2 = np.array(np.where(mito_img>0, False, True))
    #calculate the distance
    dis = sc.ndimage.distance_transform_edt(label2, sampling=spacing)
    #grab the distance of the foci only
    df['Distance'] = dis[centroid]
    #add more details to the df
    df['Condition'] = 'NI'
    df['Name'] = os.path.basename(filename)
    #concatenate to final dataframe
    final_df = pd.concat([final_df, df], ignore_index=True)

    #visu
    fig, ax = plt.subplots(figsize=(15, 5) , nrows=1, ncols=3)
    ax[0].imshow(np.sum(lbl_img, axis=0)*10)
    ax[0].set_title('vapb labl')
    ax[0].set_yticks([])
    ax[0].set_xticks([])
    ax[1].imshow(np.sum(sc.ndimage.binary_dilation(label1, iterations=5), axis=0))
    ax[1].set_title('vapb foci')
    ax[1].set_yticks([])
    ax[1].set_xticks([])
    ax[2].hist(df['Distance'], bins=25)
    plt.title(os.path.basename(filename))
    plt.tight_layout()
    plt.show()
  
 
#%%
for filename in INF_path_list:
    img = u.get_img(filename) #ZCXY format
    #split
    nucl_img = img[:,0]
    mito_img = img[:,1]
    vapb_img = img[:,2]
    del img
    
    #grab the full nucleus and substract it from mito and vapb signal
    nucl_img = sc.ndimage.gaussian_filter(nucl_img, (1,2,2))
    cleaned_data = u.boolean_remove_zeros(np.ndarray.flatten(nucl_img)) #w/o 0
    thr_min = u.thr_list(cleaned_data, 'otsu')
    nucl_binary = nucl_img > thr_min
    nucl_binary = sc.ndimage.binary_dilation(nucl_binary, iterations=7)
    nucl_binary = sc.ndimage.binary_fill_holes(nucl_binary)
    nucl_binary = sc.ndimage.binary_erosion(nucl_binary, iterations=7)
    u.ortho_view(nucl_binary, 'mean', os.path.basename(filename),ratio=[7,1])
    del nucl_img
    #remove the nucleus
    mito_img = np.where(nucl_binary == False, mito_img, 0)
    vapb_img = np.where(nucl_binary == False, vapb_img, 0)
    
    mito_img = sc.ndimage.gaussian_filter(mito_img, (.5,1,1))
    vapb_img = sc.ndimage.gaussian_filter(vapb_img, (.25,.5,.5))

    #distance map between mito and vapb
    #binary transform: all positive signal are positive!
    vapb_thr = u.thr_list(u.boolean_remove_zeros(vapb_img), 'otsu')
    label1 = np.array(np.where(vapb_img>=vapb_thr, True, False))
    #measure foci
    lbl_img = measure.label(label1)
    properties = ['label', 'area', 'centroid_weighted', 'intensity_mean']
    df = pd.DataFrame(measure.regionprops_table(lbl_img,
                                                intensity_image=vapb_img,
                                                properties=properties))
    #get the centroid position
    z_pos = np.array(df['centroid_weighted-0'])
    x_pos = np.array(df['centroid_weighted-1'])
    y_pos = np.array(df['centroid_weighted-2'])
    #assign the centroid as lbl 1
    label1 = np.zeros(vapb_img.shape, dtype=bool)
    centroid = (z_pos.astype(int), x_pos.astype(int), y_pos.astype(int))
    label1[centroid] = True
    
    #because no background, no thresholding
    label2 = np.array(np.where(mito_img>0, False, True))
    #calculate the distance
    dis = sc.ndimage.distance_transform_edt(label2, sampling=spacing)
    #grab the distance of the foci only
    df['Distance'] = dis[centroid]
    #add more details to the df
    df['Condition'] = '8 hpi'
    df['Name'] = os.path.basename(filename)
    #concatenate to final dataframe
    final_df = pd.concat([final_df, df], ignore_index=True)

    #visu
    fig, ax = plt.subplots(figsize=(15, 5) , nrows=1, ncols=3)
    ax[0].imshow(np.sum(lbl_img, axis=0)*10)
    ax[0].set_title('vapb labl')
    ax[0].set_yticks([])
    ax[0].set_xticks([])
    ax[1].imshow(np.sum(sc.ndimage.binary_dilation(label1, iterations=5), axis=0))
    ax[1].set_title('vapb foci')
    ax[1].set_yticks([])
    ax[1].set_xticks([])
    ax[2].hist(df['Distance'], bins=25)
    plt.title(os.path.basename(filename))
    plt.tight_layout()
    plt.show()

#save the dataframe as csv and pickle
final_df.to_csv('Foci_distance_analysis.csv')
final_df.to_pickle('Foci_distance_analysis.pckl')
    
    #%%
# img = nucl_img
# fig, ax = plt.subplots(figsize=(20, 5) , nrows=1, ncols=4)
# ax[0].imshow(np.sum(sc.ndimage.gaussian_filter(img, (2,5,5)), axis=0))
# ax[0].set_title('Gaussian (2,5,5)')
# ax[0].set_yticks([])
# ax[0].set_xticks([])
# ax[1].imshow(np.sum(sc.ndimage.gaussian_filter(img, (.5,.5,.5)), axis=0))
# ax[1].set_title('Gaussian 0.5, 0.5, 0.5')
# ax[1].set_yticks([])
# ax[1].set_xticks([])
# ax[2].imshow(np.sum(sc.ndimage.gaussian_filter(img, (.5,1,1)), axis=0))
# ax[2].set_title('Gaussian 0.5, 1, 1')
# ax[2].set_yticks([])
# ax[2].set_xticks([])
# ax[3].imshow(np.sum(sc.ndimage.gaussian_filter(img, (1,2,2)), axis=0))
# ax[3].set_title('Gaussian 1, 2, 2')
# ax[3].set_yticks([])
# ax[3].set_xticks([])
# plt.tight_layout()
# plt.show()


