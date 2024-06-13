# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 11:12:32 2024

@author: LECLERC


Code to extract data from the mPTP experiment

mPTP experiments notes

Acquisition protocol (Salla's email)
- Change media to imaging media 
- Take a 5 x 5 stack  (IMAGE 1.) NB each 5 x 5 stack took 10 minutes to acquire.. 
- Remove media and add a working solution containing the Calcein and CoCl2 as a quencher.
- Take a 5 x 5 stack (IMAGE 2.) --> this image equals to ~0 min of treatment 
- Incubate for 25 min 
- Take a 5 x 5 stack  (IMAGE 3.)--> this image equals to ~25 min of treatment
- Wash with dilution buffer 
- Take a 5 x 5 stack  (IMAGE 4.) in principle better DIC.
- Add ionomycin and incubate for 5 min. 
- Take a 5 x 5 stack  (IMAGE 5.)

Channel 1 is mPTP/infection
Channel 2 is DIC
Channel 3 is MitoTracker

It is possible to do the analysis from IMAGES 2. and 3. or IMAGES 1 and 4. (Salla's comment)

Will be difficult to do the montage since the Z-plane moved quite a bit (Z drift)

Channel 1 show the most changes: in image 1 and 2, only see the ICP4, while in image 3 and 4, can see the mitochondria.
Channel 2 (DIC) will be difficult to use. Image 1 is good, but the all the other presents some aggregate. Image 2 may be useable though
Channel 3 is presents only in image 2 to 4, and very strong. Should label the mitochondria, but label the mitochondria AND the cell. In some condition, perfect labeling of the mitochondria.

Analysis
Based on Image 2 and 3 couple, extract the infection status and level of the cell by the total fluorescence intensity in the nucleus, excluding NI cell. Substract the signal from image 2 to image 3 to grab only the mPTP fluorescence value. Grab the cytoplasm of each mitochondria based on the mitotracker (or the DIC?), then check that the mPTP signal is coming from the mitochondria. As for visualization, mPTP quantity only, infection status, normalization of the mPTP versus the infection (other?)
Protocol
-Max project Image 2 and 3
-On image 3, channel 3, threshold (normal/default)
	- found the number of nuclei based on the number of big holes
	- Bonus point if it is possible to segment them
-On image 2, channel 1, threshold (Otsu) -> infected nucleus
-Combine the two threshold image or on the nuclei segmentation
-Binary Open, fill holes then erode the image to remove noise
-Watershed based on the number of nuclei ->cell cytoplasm
-Label each cell
-On image 3, channel 3, segment the mitochondria on the MAX projection (Otsu)
-For each cell, measure the quantity of mPTP (channel 1) signal in the cytoplasm (exclude the nucleus) on both image 2 and 3
-Measure the infection status in (total signal quantity) in each cell in image 2 only, nucleus only
-Save and export the data as dataframe

Fuse the image to have only one big DONE
Change the algo to process just this one instead!

TODO process the NI (different!)

"""

#IMPORTS
import numpy as np
import tifffile
import os
import matplotlib.pyplot as plt
from scipy import ndimage

from skimage.filters import threshold_local
from skimage.measure import label, regionprops_table
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage.transform import resize

import pandas as pd

import utils as u
#%%

#FUNCTIONS

def pair_nuclei_cyto(nuc, data_cyto, strech):
    nuc0 = (nuc['centroid-0']-strech//2, nuc['centroid-0']+strech//2)
    nuc1 = (nuc['centroid-1']-strech//2, nuc['centroid-1']+strech//2)
    potential_label = []
    for idy, cyt in data_cyto.iterrows():
        if nuc0[0] < cyt['centroid-0'] and cyt['centroid-0']< nuc0[1]:
            if nuc1[0] < cyt['centroid-1'] and cyt['centroid-1']< nuc1[1]:
                potential_label.append(cyt.label)
    return potential_label

def channel_smooth(img, sigma):
    smooth_img = np.zeros(img.shape, dtype=float)
    for channel in range(len(img2)):
        smooth_img[channel] = ndimage.gaussian_filter(img[channel], sigma)
    return smooth_img

def nuclei_detection_from_mitotracker(img, block_size=120):
    """
    Roughly identify nuclei from the mitotracker label from the holes they are creating.
    Start by doing a local thresholding followed by a binary closing to clean the data
    Label the nuclei and remove too small and too big ones, as well as too loose one

    Parameters
    ----------
    img : np.array
        2d array
    block_size : int, optional
        The size of the block for the local segmentation. The default is 120.

    Returns
    -------
    nuclei_b : np.array
        2d array, the nuclei binary image
    local_img : np.array
        2d array, the locally thresholded image

    """
    block = 120
    if block%2 == 0: block += 1 #block size need to be an odd number
    local_img = img > threshold_local(img,
                                          block_size=block,
                                          offset=0,
                                          method='gaussian'
                                          )

    #binary manipulation to clean the data
    local_img = ndimage.binary_closing(local_img, iterations=1)
       
    #invert the data
    local_img_i = np.where(local_img, False, True)
    #label it
    local_img_l = label(local_img_i)
    #measure it
    props = regionprops_table(local_img_l, img,
                              properties=['label', 'area', 'solidity'])
    data = pd.DataFrame(props)
    data = data[(data.area>300) & (data.area<10000)] #size filter
    data = data[data.solidity > 0.7]
    
    nuclei_b = np.zeros(local_img_l.shape) #recipient
    #regenerate the img without the junk
    for n, row in data.iterrows():
        #grab the label only
        n_b = np.where(local_img_l==row.label, True, False)
        n_b = ndimage.binary_dilation(n_b)
        n_b = ndimage.binary_fill_holes(n_b)
        nuclei_b += n_b
    
    return nuclei_b, local_img


#%%
#MAIN
condition = ['4hpi', '8hpi', '12hpi']
for c in condition:
    print('Processing condition '+c)
    img_path = u.found_img(c+os.sep)
    img2_path = [x for x in img_path if 'IMAGE 2' in x][0]
    img3_path = [x for x in img_path if 'IMAGE 3' in x][0]
    #open the image
    img2 = tifffile.imread(img2_path)
    img3 = tifffile.imread(img3_path)
    #check if smooth img2 and 3 have the same size, if not resize to 2
    if img2.shape != img3.shape:
        img3 = resize(img3, img2.shape, preserve_range=True)
    
    #smooth the data
    s = 2.0
    smooth_img2 = channel_smooth(img2, s)
    smooth_img3 = channel_smooth(img3, s)
 
    #visualize both images
    fig = plt.figure(figsize=(20,20))
    fig.add_subplot(221)
    plt.imshow(smooth_img2[0])
    plt.title('Image 2 - ICP4')
    plt.axis('off')
    fig.add_subplot(222)
    plt.imshow(smooth_img2[2])
    plt.title('Image 2 - Mito')
    plt.axis('off')
    fig.add_subplot(223)
    plt.imshow(smooth_img3[0])
    plt.title('Image 3 - ICP4+mPTP')
    plt.axis('off')
    fig.add_subplot(224)
    plt.imshow(smooth_img3[2])
    plt.title('Image 3 - Mito')
    plt.axis('off')
    
    plt.tight_layout()
    plt.show()


    #%%
    #Objective: found the number and localization of nuclei   
    #local thresholding perform better than global

    nuclei_b, local_img2 = nuclei_detection_from_mitotracker(smooth_img2[2], block_size=120)
    #label the nuclei
    nuclei = label(nuclei_b)
    #measure it
    props = regionprops_table(nuclei, img2[2],
                              properties=['label', 'area', 'centroid',
                                          'eccentricity', 'intensity_mean',
                                          'solidity', 'coords'])
    data_nuclei = pd.DataFrame(props)
    
    #remove any small artifact, security
    for n, row in data_nuclei.iterrows():
        if row.area < 300: #minsize
            nuclei = np.where(nuclei==row.label, 0, nuclei) #remove the label
    data_nuclei = data_nuclei[data_nuclei.area>300] #remove from the listing as well
    
    #remove the nuclei label from the cytoplasm (due to fill holes)
    local_img2 = ndimage.binary_erosion(local_img2, iterations=2)
    cyto = local_img2 & np.where(nuclei>0, False, True)
    nuclei_cyto_b = cyto+ndimage.binary_dilation(np.where(nuclei>0, True, False))
    nuclei_cyto_b = ndimage.binary_fill_holes(nuclei_cyto_b)
    
    #visualization of detected nuclei    
    plt.imshow(u.img_color_merger(g=img2[2]*2, #'raw'
                                  r=np.where(nuclei>0, 205, 0), #nuclei
                                  b=np.where(cyto>0, 205, 0)), #cytoplasm
               interpolation='none')
    plt.axis('off')
    plt.show()
    #%%
    # Watershed + associate nucleus to cytplasm
    distance = ndimage.distance_transform_edt(nuclei_cyto_b)
    mask = np.zeros(distance.shape, dtype=bool)
    coo = np.array((data_nuclei['centroid-0'], data_nuclei['centroid-1']), dtype=int).T
    mask[coo[:,0], coo[:,1]] = True
    markers, _ = ndimage.label(mask)
    labels = watershed(-distance, markers,
                       connectivity=[[0,1,0],[1,1,1],[0,1,0]],
                       compactness=0.5,
                       mask=nuclei_cyto_b)

    #infinite loop to breakdown loose cell
    counter = 0
    len_result = 0
    while True:
        counter += 1
        #measure the different labeled cell
        props = regionprops_table(labels, img2[2],
                                  properties=['label', 'area', 'centroid',
                                              'eccentricity', 'intensity_mean',
                                              'solidity', 'coords'])
        data_cyto = pd.DataFrame(props)
        
        result = []
        for idx, row in data_cyto.iterrows():
            #only loose cell, eccentricity close to 1 and low solidity
            if row.eccentricity > 0.90 and row.solidity < 0.65:
                result.append(False) #to repeat the cycle
                n_b = np.where(labels == row.label, True, False) #grab the cell
                #watershed the cell to split it in two pieces
                distance = ndimage.distance_transform_edt(n_b)
                coords = peak_local_max(distance,
                                        footprint=np.ones((3, 3)),
                                        min_distance=65)
                mask = np.zeros(distance.shape, dtype=bool)
                mask[tuple(coords.T)] = True
                markers, _ = ndimage.label(mask)
                ls = watershed(-distance,
                               markers=markers,
                               compactness=0.1,
                               mask=n_b)
                
                #remove old label
                to_replace_idx = np.where(labels == row.label)
                labels[to_replace_idx] = 0
                
                uu, c = np.unique(ls, return_counts=True)
                for x, n in enumerate(uu):
                    if n == 0: continue #background
                    if c[x] < 10: continue #too small, ignore
                    to_replace_idx = np.where(ls == n)
                    #replace the old label with new label
                    labels[to_replace_idx] = np.max(labels)+1
            # compact cell
            else:
                result.append(True)
        #check the progression
        if len_result != len(np.unique(labels)):
            len_result = len(np.unique(labels))
        else:
            print('No change in watershed anymore, abort')
            break
        #if only compact cell are left, exit the loop
        if all(result) == True:
            print('Recurrent watershed completed')
            break
        if counter == 10: #security
            print('Abort watershed for security')
            break
        
    #%%
    
    #visualization
    name = os.path.basename(img2_path.replace('.tif', '.png'))
    name = name.replace('IMAGE 2_', '')
    
    fig, axes = plt.subplots(ncols=5, figsize=(40, 10), sharex=True, sharey=True)
    ax = axes.ravel()

    ax[0].imshow(smooth_img2[0])
    ax[0].set_title('ICP4')
    ax[1].imshow(np.clip(smooth_img3[0]-smooth_img2[0], 0, 255), interpolation=None)
    ax[1].set_title('mPTP')    
    ax[2].imshow(smooth_img2[2], interpolation=None)
    ax[2].set_title('Mitotracker')
    ax[3].imshow(nuclei, cmap=plt.cm.gist_ncar)
    ax[3].set_title('Nucleus Label')
    ax[4].imshow(labels, cmap=plt.cm.gist_ncar)
    ax[4].set_title('Cell Label')
    
    for a in ax:
        a.set_axis_off()
    
    fig.tight_layout()
    plt.savefig('result_img'+os.sep+name, dpi=150)
    plt.show()
    
    
    
    #%%Fused the cyto and nuclei dataframe
    #using the centroid of cytoplasm and nuclei, assign an label for both    
    cyto_idx = []
    strech = 10
    for idx, nuc in data_nuclei.iterrows():
        while True:
            potential_label = pair_nuclei_cyto(nuc, data_cyto, strech)
            #check and remove if potential label is already attributed
            potential_label = [x for x in potential_label if x not in cyto_idx]
            if len(potential_label) != 0:
                break
            else:
                strech += 5
            if strech>100: break #security
        
        if len(potential_label) > 1: #multiple choice, choose the biggest not attributed
            sub = data_cyto[data_cyto.label.isin(potential_label)]
            cyto_idx.append(int(sub[sub.area==sub.area.max()].label))                
        elif len(potential_label) == 1: #perfect, attribute the cyto to this nuclei
            cyto_idx.append(potential_label[0])
        else: #no cytoplasm associated to this nucleus
            cyto_idx.append(-1)
            print('Nuclei ID '+str(nuc.label)+' does not pair with a cytoplasm. Deleted')
            continue
            
        #print('Nuclei '+str(nuc.label)+' fit with cyto label '+str(potential_label))
    
    data_nuclei['Cyto_label'] = cyto_idx
    
    #remove abnormal entries
    data_nuclei = data_nuclei[data_nuclei.Cyto_label != -1]

    
    #%% full data fusion
    full_data = []
    for idx, nuc in data_nuclei.iterrows():
        d = {}
        d['nucleus_ID'] = nuc.label
        d['nucleus_area'] = int(nuc.area)
        d['nucleus_coords'] = nuc.coords
        
        d['cyto_area'] = int(data_cyto[data_cyto.label==nuc.Cyto_label].area)
        d['cyto_coords'] = data_cyto[data_cyto.label==nuc.Cyto_label].coords.squeeze()
        
        full_data.append(d)
    full_data = pd.DataFrame(full_data)
    
    #%%extracting data from image and saving in the dataframe
    #Check if nuclei is infected, looking at the 95 percentile value
    infected_int = []
    for idx, cell in full_data.iterrows():
        intensity = []
        for pxl in cell.nucleus_coords:
            intensity.append(img2[[0],pxl[0],pxl[1]])
        infected_int.append(np.percentile(intensity, 95))
    full_data['Infection_intensity'] = infected_int
    
    #check the cytoplasm intensity for the mitotracker on Image 3
    mitochondria_intensity = [] #average
    mPTP_before = []
    mPTP_after = []
    mPTP_diff = []
    for idx, cell in full_data.iterrows():
        nucleus_empty = np.zeros(img2[0].shape, dtype=bool)
        cyto_empty = np.zeros(img2[0].shape, dtype=bool)
        for pxl in cell.nucleus_coords:
            nucleus_empty[pxl[0],pxl[1]] = True
        for pxl in cell.cyto_coords:
            cyto_empty[pxl[0],pxl[1]] = True
        #extract the image of th cytoplasm for the mitochondria of the target cell while removing the nucleus
        cyto_mito = np.where(cyto_empty ^ nucleus_empty, img2[2], 0)
        cyto_idx = np.where(cyto_mito > 0)
        #mitochondria intensity
        mitochondria_intensity.append(np.mean(img2[2][cyto_idx])) #average intensity    
        #extract noisy mean mPTP/ICP4 in the mitochondria BEFORE treatment
        mPTP_before.append(np.mean(img2[0][cyto_idx]))
        #extract mean mPTP/ICP4 in the mitochondria AFTER treatment
        mPTP_after.append(np.mean(img3[0][cyto_idx]))   
        #the mean of the difference of the mPTP3-2
        mPTP_diff.append(np.mean(img3[0][cyto_idx]-img2[0][cyto_idx]))
    
    full_data['cyto_mito_intensity'] = mitochondria_intensity
    full_data['AVG_mPTP_IMAGE2_intensity'] = mPTP_before
    full_data['AVG_mPTP_IMAGE3_intensity'] = mPTP_after
    full_data['AVG_mPTP_image_difference'] = mPTP_diff
    #%%Quick visualization
    fig = plt.figure(figsize=(10,10))
    fig.add_subplot(231)
    plt.scatter(full_data['Infection_intensity'], full_data['AVG_mPTP_image_difference'])
    plt.ylabel('mPTP')
    plt.xlabel('Infection')
    
    fig.add_subplot(232)
    plt.scatter(full_data['cyto_mito_intensity'], full_data['AVG_mPTP_image_difference'])
    plt.ylabel('mPTP')
    plt.xlabel('Mitochondria')
    
    fig.add_subplot(233)
    plt.scatter(full_data['Infection_intensity'], full_data['cyto_mito_intensity'])
    plt.xlabel('Infection')
    plt.ylabel('Mitochondria')    
    
    fig.add_subplot(234)
    plt.scatter(full_data['AVG_mPTP_IMAGE2_intensity'], full_data['AVG_mPTP_IMAGE3_intensity'])
    plt.xlabel('Before')
    plt.ylabel('After')

    fig.add_subplot(235)
    plt.scatter(full_data['Infection_intensity'], full_data['AVG_mPTP_IMAGE3_intensity'])
    plt.xlabel('Infection')
    plt.ylabel('After')    
    
    plt.tight_layout()
    plt.savefig('result_img'+os.sep+name.replace('.png', '_graph.png'), dpi=150)
    plt.show()
    
    #%% Save the dataFrame for the field of view
    name = img2_path.replace('.tif', '.pckl')
    name = name.replace('IMAGE 2_', '')
    full_data.to_pickle(name)






#%%Process the NI as well



print('Processing NI WS')


NI_path = 'NI WS'+os.sep+'noninfceted control working solution 12042024.lif - control noninfected with working solution IMAGE 4_Processed001_Merging_001.tif'
ni = tifffile.imread(NI_path)
smooth_ni = channel_smooth(ni, s)
nuclei_b, local_img2 = nuclei_detection_from_mitotracker(smooth_ni[2], block_size=120)
nuclei = label(nuclei_b)
props = regionprops_table(nuclei, ni[2],
                          properties=['label', 'area', 'centroid',
                                      'eccentricity', 'intensity_mean',
                                      'solidity', 'coords'])
data_nuclei = pd.DataFrame(props)
#remove any small artifact, security
for n, row in data_nuclei.iterrows():
    if row.area < 300: #minsize
        nuclei = np.where(nuclei==row.label, 0, nuclei) #remove the label
data_nuclei = data_nuclei[data_nuclei.area>300] #remove from the listing as well

#remove the nuclei label from the cytoplasm (due to fill holes)
local_img2 = ndimage.binary_erosion(local_img2, iterations=2)
cyto = local_img2 & np.where(nuclei>0, False, True)
nuclei_cyto_b = cyto+ndimage.binary_dilation(np.where(nuclei>0, True, False))
nuclei_cyto_b = ndimage.binary_fill_holes(nuclei_cyto_b)
#visualization of detected nuclei    
plt.imshow(u.img_color_merger(g=ni[2]*2, #'raw'
                              r=np.where(nuclei>0, 205, 0), #nuclei
                              b=np.where(cyto>0, 205, 0)), #cytoplasm
           interpolation='none')
plt.axis('off')
plt.show()
#%% Watershed + associate nucleus to cytplasm
distance = ndimage.distance_transform_edt(nuclei_cyto_b)
mask = np.zeros(distance.shape, dtype=bool)
coo = np.array((data_nuclei['centroid-0'], data_nuclei['centroid-1']), dtype=int).T
mask[coo[:,0], coo[:,1]] = True
markers, _ = ndimage.label(mask)
labels = watershed(-distance, markers,
                   connectivity=[[0,1,0],[1,1,1],[0,1,0]],
                   compactness=0.5,
                   mask=nuclei_cyto_b)

#infinite loop to breakdown loose cell
counter = 0
len_result = 0
while True:
    counter += 1
    #measure the different labeled cell
    props = regionprops_table(labels, ni[2],
                              properties=['label', 'area', 'centroid',
                                          'eccentricity', 'intensity_mean',
                                          'solidity', 'coords'])
    data_cyto = pd.DataFrame(props)
    
    result = []
    for idx, row in data_cyto.iterrows():
        #only loose cell, eccentricity close to 1 and low solidity
        if row.eccentricity > 0.90 and row.solidity < 0.65:
            result.append(False) #to repeat the cycle
            n_b = np.where(labels == row.label, True, False) #grab the cell
            #watershed the cell to split it in two pieces
            distance = ndimage.distance_transform_edt(n_b)
            coords = peak_local_max(distance,
                                    footprint=np.ones((3, 3)),
                                    min_distance=65)
            mask = np.zeros(distance.shape, dtype=bool)
            mask[tuple(coords.T)] = True
            markers, _ = ndimage.label(mask)
            ls = watershed(-distance,
                           markers=markers,
                           compactness=0.1,
                           mask=n_b)
            
            #remove old label
            to_replace_idx = np.where(labels == row.label)
            labels[to_replace_idx] = 0
            
            uu, c = np.unique(ls, return_counts=True)
            for x, n in enumerate(uu):
                if n == 0: continue #background
                if c[x] < 10: continue #too small, ignore
                to_replace_idx = np.where(ls == n)
                #replace the old label with new label
                labels[to_replace_idx] = np.max(labels)+1
        # compact cell
        else:
            result.append(True)
    #check the progression
    if len_result != len(np.unique(labels)):
        len_result = len(np.unique(labels))
    else:
        print('No change in watershed anymore, abort')
        break
    #if only compact cell are left, exit the loop
    if all(result) == True:
        print('Recurrent watershed completed')
        break
    if counter == 10: #security
        print('Abort watershed for security')
        break
#visualization
name = os.path.basename(NI_path.replace('.tif', '.png'))
name = name.replace('IMAGE 4_', '')

fig, axes = plt.subplots(ncols=5, figsize=(40, 10), sharex=True, sharey=True)
ax = axes.ravel()

ax[0].imshow(smooth_ni[0])
ax[0].set_title('ICP4')
ax[1].imshow(np.clip(smooth_ni[0], 0, 255), interpolation=None)
ax[1].set_title('mPTP')    
ax[2].imshow(smooth_ni[2], interpolation=None)
ax[2].set_title('Mitotracker')
ax[3].imshow(nuclei, cmap=plt.cm.gist_ncar)
ax[3].set_title('Nucleus Label')
ax[4].imshow(labels, cmap=plt.cm.gist_ncar)
ax[4].set_title('Cell Label')

for a in ax:
    a.set_axis_off()

fig.tight_layout()
plt.savefig('result_img'+os.sep+name, dpi=150)
plt.show()


#%%Fused the cyto and nuclei dataframe
#using the centroid of cytoplasm and nuclei, assign an label for both    
cyto_idx = []
strech = 10
for idx, nuc in data_nuclei.iterrows():
    while True:
        potential_label = pair_nuclei_cyto(nuc, data_cyto, strech)
        #check and remove if potential label is already attributed
        potential_label = [x for x in potential_label if x not in cyto_idx]
        if len(potential_label) != 0:
            break
        else:
            strech += 5
        if strech>100: break #security
    
    if len(potential_label) > 1: #multiple choice, choose the biggest not attributed
        sub = data_cyto[data_cyto.label.isin(potential_label)]
        cyto_idx.append(int(sub[sub.area==sub.area.max()].label))                
    elif len(potential_label) == 1: #perfect, attribute the cyto to this nuclei
        cyto_idx.append(potential_label[0])
    else: #no cytoplasm associated to this nucleus
        cyto_idx.append(-1)
        print('Nuclei ID '+str(nuc.label)+' does not pair with a cytoplasm. Deleted')
        continue
        
    #print('Nuclei '+str(nuc.label)+' fit with cyto label '+str(potential_label))

data_nuclei['Cyto_label'] = cyto_idx

#remove abnormal entries
data_nuclei = data_nuclei[data_nuclei.Cyto_label != -1]


#%% full data fusion
full_data = []
for idx, nuc in data_nuclei.iterrows():
    d = {}
    d['nucleus_ID'] = nuc.label
    d['nucleus_area'] = int(nuc.area)
    d['nucleus_coords'] = nuc.coords
    
    d['cyto_area'] = int(data_cyto[data_cyto.label==nuc.Cyto_label].area)
    d['cyto_coords'] = data_cyto[data_cyto.label==nuc.Cyto_label].coords.squeeze()
    
    full_data.append(d)
full_data = pd.DataFrame(full_data)

#%%extracting data from image and saving in the dataframe
#Check if nuclei is infected, looking at the 95 percentile value
infected_int = []
for idx, cell in full_data.iterrows():
    intensity = []
    for pxl in cell.nucleus_coords:
        intensity.append(ni[[0],pxl[0],pxl[1]])
    infected_int.append(np.percentile(intensity, 95))
full_data['Infection_intensity'] = infected_int

#check the cytoplasm intensity for the mitotracker on Image 3
mitochondria_intensity = [] #average
mPTP_before = []
mPTP_after = []
mPTP_diff = []
for idx, cell in full_data.iterrows():
    nucleus_empty = np.zeros(ni[0].shape, dtype=bool)
    cyto_empty = np.zeros(ni[0].shape, dtype=bool)
    for pxl in cell.nucleus_coords:
        nucleus_empty[pxl[0],pxl[1]] = True
    for pxl in cell.cyto_coords:
        cyto_empty[pxl[0],pxl[1]] = True
    #extract the image of th cytoplasm for the mitochondria of the target cell while removing the nucleus
    cyto_mito = np.where(cyto_empty ^ nucleus_empty, ni[2], 0)
    cyto_idx = np.where(cyto_mito > 0)
    #mitochondria intensity
    mitochondria_intensity.append(np.mean(ni[2][cyto_idx])) #average intensity    
    #extract noisy mean mPTP/ICP4 in the mitochondria BEFORE treatment
    mPTP_before.append(np.mean(ni[0][cyto_idx]))
    #extract mean mPTP/ICP4 in the mitochondria AFTER treatment
    mPTP_after.append(np.mean(ni[0][cyto_idx]))   
    #the mean of the difference of the mPTP3-2
    mPTP_diff.append(np.mean(ni[0][cyto_idx]-ni[0][cyto_idx]))

full_data['cyto_mito_intensity'] = mitochondria_intensity
full_data['AVG_mPTP_IMAGE2_intensity'] = mPTP_before
full_data['AVG_mPTP_IMAGE3_intensity'] = mPTP_after
full_data['AVG_mPTP_image_difference'] = mPTP_diff
#%%Quick visualization
fig = plt.figure(figsize=(10,10))
fig.add_subplot(231)
plt.scatter(full_data['Infection_intensity'], full_data['AVG_mPTP_image_difference'])
plt.ylabel('mPTP')
plt.xlabel('Infection')

fig.add_subplot(232)
plt.scatter(full_data['cyto_mito_intensity'], full_data['AVG_mPTP_image_difference'])
plt.ylabel('mPTP')
plt.xlabel('Mitochondria')

fig.add_subplot(233)
plt.scatter(full_data['Infection_intensity'], full_data['cyto_mito_intensity'])
plt.xlabel('Infection')
plt.ylabel('Mitochondria')    

fig.add_subplot(234)
plt.scatter(full_data['AVG_mPTP_IMAGE2_intensity'], full_data['AVG_mPTP_IMAGE3_intensity'])
plt.xlabel('Before')
plt.ylabel('After')

fig.add_subplot(235)
plt.scatter(full_data['Infection_intensity'], full_data['AVG_mPTP_IMAGE3_intensity'])
plt.xlabel('Infection')
plt.ylabel('After')    

plt.tight_layout()
plt.savefig('result_img'+os.sep+name.replace('.png', '_graph.png'), dpi=150)
plt.show()

#%% Save the dataFrame for the field of view
name = NI_path.replace('.tif', '.pckl')
name = name.replace('IMAGE 4_', '')
full_data.to_pickle(name)










