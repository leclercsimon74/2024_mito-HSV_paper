# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 08:34:25 2023

@author: sleclerc
"""

import utils as u

import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from skimage import measure
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
import os
import pandas as pd

sampling = (0.1797678,0.0632473,0.0632473) #resolution in um ZXY


def PLA_analysis(condition, result_df):
    """
    Analyse the images in the folder named 'condition'.
    Detect the nucleus (channel:0) by otsu segmentation
    Substract the nucleus seg from the pca segmentation (triangle)
    Watershed the pla seg, label it and measure each dot for its area and average signal intensity
    Segment the mitochondria with otsu
    Measure the distance between the pca seg foci to the closest mitochondria
    
    Save everything in the result_df
    

    Parameters
    ----------
    condition : str
        Folder name, no iterative
    result_df : pandas.DataFrame
        To save and merge the result

    Returns
    -------
    result_df : pandas.DataFrame
        To save and merge the result

    """
    paths = u.sorted_nicely(u.found_img(condition))
    for img_path in paths:
        print('Processing of '+os.path.basename(img_path))
        img = u.get_img(img_path) #ZCXY
        nucl_img = img[:,0]
        pca_img = img[:,2]
        mito_img = img[:,4]
        
        del img
        
        nucl_seg, s = u.nucleus_detection(nucl_img, '')
        fig = plt.figure(figsize=(15,5))
        ax = fig.add_subplot(151)
        ax.imshow(np.sum(nucl_seg, axis=0))
        ax.set_title('Nucleus SUM')
        plt.axis('off')
        
        del nucl_img, s
        #nucleus detection
        pca_img = sc.ndimage.gaussian_filter(pca_img, (0.5,1,1))
        pca_seg = pca_img >= u.thr_list(u.boolean_remove_zeros(pca_img), 'triangle')
        pca_seg = np.where(nucl_seg==False, pca_seg, False) #remove the nucleus
        pca_seg = sc.ndimage.binary_opening(pca_seg, iterations=1)
        
        ax = fig.add_subplot(152)
        ax.imshow(np.max(pca_seg, axis=0))
        ax.set_title('PLA MAX')
        plt.axis('off')
        
        #Watershed
        distance = sc.ndimage.distance_transform_edt(pca_seg, sampling=sampling)
        coords = peak_local_max(distance, footprint=np.ones((3, 3, 3)),
                                labels=pca_seg)
        mask = np.zeros(distance.shape, dtype=bool)
        mask[tuple(coords.T)] = True
        markers, _ = sc.ndimage.label(mask)
        labels = watershed(-distance, markers, mask=pca_seg)
        
        properties = ['label', 'area', 'centroid_weighted', 'intensity_mean']
        df = pd.DataFrame(measure.regionprops_table(labels,
                                                    intensity_image=pca_img,
                                                    properties=properties,
                                                    spacing=sampling
                                                    ))
        
        df['Name'] = os.path.basename(img_path)
        df['Condition'] = condition
    
        del distance, coords, mask, markers, labels
        
        #get the centroid position
        z_pos = np.array(df['centroid_weighted-0'])
        x_pos = np.array(df['centroid_weighted-1'])
        y_pos = np.array(df['centroid_weighted-2'])
    
        #assign the centroid in label = center of the foci
        label = np.zeros(pca_img.shape, dtype=bool)
        centroid = (z_pos.astype(int), x_pos.astype(int), y_pos.astype(int))
        label[centroid] = True
        #quick visu
        ax = fig.add_subplot(153)
        ax.imshow(np.sum(sc.ndimage.binary_dilation(label, iterations=5), axis=0))
        ax.set_title('PLA foci')
        plt.axis('off')
        
        #mito detection
        mito_img = sc.ndimage.gaussian_filter(mito_img, (0.5,1,1))
        mito_seg = mito_img >= u.thr_list(u.boolean_remove_zeros(mito_img), 'otsu')
        mito_seg = np.where(nucl_seg==False, mito_seg, False) #remove the nucleus
        mito_seg = sc.ndimage.binary_opening(mito_seg, iterations=1)
        ax = fig.add_subplot(154)
        ax.imshow(np.sum(mito_seg, axis=0))
        ax.set_title('Mito')
        plt.axis('off')    
        
        #distance PCA to mito
        dis = sc.ndimage.distance_transform_edt(mito_seg, sampling=sampling)
        df['Distance'] = dis[centroid]
        #save the result in the df
        result_df = pd.concat([result_df, df], ignore_index=True)
        
        ax = fig.add_subplot(155)
        l = np.zeros(dis.shape)
        l[centroid] = dis[centroid] +1 #to see even the 0 distance!
        ax.imshow(np.max(sc.ndimage.grey_dilation(l, footprint=np.ones((5,5,5))), axis=0))
        ax.set_title('Distance to Mito')
        plt.axis('off')
    
        plt.tight_layout()
        plt.savefig('Results_IMG'+os.sep+condition+'_'+os.path.basename(img_path))
        plt.show()
        
        print('Found '+str(len(df))+' PLA foci')
    
    return result_df
        


result_df = pd.DataFrame()
result_df = PLA_analysis('NI', result_df)
result_df = PLA_analysis('4 hpi', result_df)
result_df = PLA_analysis('8 hpi', result_df)
result_df = PLA_analysis('12 hpi', result_df)

result_df.to_csv('PLA_foci_distance.csv')
result_df.to_pickle('PLA_foci_distance.pckl')