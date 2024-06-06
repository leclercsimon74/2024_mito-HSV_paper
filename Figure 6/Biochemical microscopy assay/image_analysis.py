# -*- coding: utf-8 -*-
"""
Created on Fri May 19 09:57:59 2023

@author: sleclerc

Description:

The provided code is a simple and efficient tool for extracting data from TIFF 
images. It is designed to process a collection of images, each consisting of 
multiple channels. The code performs various image analysis and segmentation 
tasks to extract relevant information from the images.

The runtime of the code is approximately 5 minutes when processing 13 images,
each with a resolution of 3789x3789 pixels.

The code assumes that the images are composed of several channels, with each
channel representing a different aspect of the cellular structures or markers. 
The channel assignments are as follows:
- Channel 0 corresponds to the nuclear label, typically stained with DAPI.
- Channel 1 represents a marker for infection, specifically the nuclear 
    localization, and is visualized using YFP fluorescence.
- Channel 2 is associated with a mitochondrial function marker.
- Channel 3 contains information from the mitotracker.
- Channel 4, although not utilized in this code, is designated for DIC
    (differential interference contrast) imaging.

The code follows a step-by-step approach to extract data from the images:

1. Nucleus Segmentation: The code performs image segmentation to isolate the 
 nucleus in the image. It applies histogram equalization and smoothing techniques
 to enhance the nuclear region and then applies thresholding to separate it from
 the background. Additional cleaning operations, such as filling holes and opening,
 are performed to refine the segmentation. Finally, the labeled regions are
 obtained using the "measure.label" function.
2. Infected Nucleus Detection: The code identifies the infected nuclei by 
 analyzing the marker of infection (Channel 1). It applies histogram equalization
 and smoothing to improve the marker's visibility and uses thresholding to separate
 the infected regions from the background. Similar to the nucleus segmentation
 step, cleaning operations are performed to refine the segmentation.
3. Cytoplasm Segmentation: The code applies segmentation techniques to identify
 the cytoplasmic regions using the mitotracker channel (Channel 3). It performs
 histogram equalization, smoothing, and thresholding to extract the cytoplasmic
 regions. Additional opening operations are carried out for further refinement.
4. Nucleus-Cytoplasm Association: To link each nucleus to its corresponding
 cytoplasmic region, the code utilizes the centroid coordinates of the infected
 nuclei. It creates a mask using the centroid positions and applies the watershed
 algorithm using the cytoplasmic regions and the mask. This process ensures
 accurate association between nuclei and their respective cytoplasmic regions.
5. Measurement and Data Extraction: The code calculates various measurements and
 extracts relevant data from the original image. It calculates the areas of the
 nucleus, cytoplasm, and infected regions. Additionally, it computes the mean
 and standard deviation of the intensity values for different channels within
 each cell. The results are stored in a Pandas DataFrame for easy analysis and
 further processing.
6. Visualization: The code generates visualizations to aid in understanding and
 verifying the segmentation and association results. It creates a figure with
 subplots showing the segmented nucleus, segmented infection regions, and
 labeled cells after the watershed algorithm. These figures are saved as PNG
 files for later reference.

"""


from scipy import ndimage
from skimage import exposure, filters, measure, segmentation
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tifffile
import os

files = os.listdir(os.getcwd())
files = [x for x in files if x.endswith('.tif')]

for f in files:
    print('Processing of '+f)
    img = tifffile.imread(os.getcwd()+os.sep+f)
    #%% Nucleus detection
    #equalize the histogram
    nuc_o = exposure.equalize_adapthist(img[0], kernel_size=256, clip_limit=0.1)
    #smooth the histogram
    nuc = filters.gaussian(nuc_o, sigma=3)
    #threshold
    thr = filters.threshold_otsu(nuc)
    nuc = nuc > thr
    #Cleaning
    nuc = ndimage.binary_fill_holes(nuc)
    nuc = ndimage.binary_opening(nuc, iterations=10)
    #label the image!
    nuc_lbl = measure.label(nuc)
       
    #%% Detect if nucleus is infected
    inf_o = exposure.equalize_adapthist(img[1], clip_limit=0.01)
    #smooth the histogram
    inf = filters.gaussian(inf_o, sigma=4)
    #threshold
    thr = filters.threshold_triangle(inf)
    inf = inf > thr
    #cleaning
    inf = ndimage.binary_fill_holes(inf)
    inf = ndimage.binary_opening(inf, iterations=5)
    
    #measure the nucleus vs the infection marker
    props = ('label', 'area', 'intensity_mean', 'centroid')
    nuc_df = pd.DataFrame(measure.regionprops_table(nuc_lbl, intensity_image=inf, properties=props))
    #select only nuclei with yfp signal in it
    if 'NI' in f:
        pass
    else:
        nuc_df = nuc_df[nuc_df.intensity_mean > 0]
        nuc_df.reset_index(inplace=True)
    #%% Detect cytoplasm using mitotracker
    cyt_o = exposure.equalize_adapthist(img[3], clip_limit=0.05)
    #smooth the histogram
    cyt = filters.gaussian(cyt_o, sigma=2)
    #threshold
    thr = filters.threshold_triangle(cyt)
    cyt = cyt > thr
    #cleaning
    cyt = ndimage.binary_opening(cyt, iterations=10)
    # Watershed + associate nucleus to cytplasm
    distance = ndimage.distance_transform_edt(cyt)
    mask = np.zeros(distance.shape, dtype=bool)
    coo = np.array((nuc_df['centroid-0'], nuc_df['centroid-1']), dtype=int).T #using infected nuclei!
    mask[coo[:,0], coo[:,1]] = True
    markers, _ = ndimage.label(mask)
    labels = segmentation.watershed(-distance, markers, mask=cyt)
    
    #%%Visualization
    fig, axes = plt.subplots(ncols=3, figsize=(15, 5), sharex=True, sharey=True)
    ax = axes.ravel()
    
    ax[0].imshow(nuc)
    ax[0].set_title('All Seg. Nucleus')
    ax[1].imshow(inf)
    ax[1].set_title('Seg. Infection')
    ax[2].imshow(labels, cmap=plt.cm.nipy_spectral)
    ax[2].set_title('Cell Label')
    
    for a in ax:
        a.set_axis_off()
    
    fig.tight_layout()
    plt.savefig(f.replace('_XY1.ome.tif', 'fig.png'))
    plt.show()
    
    #%% Measure signal
    df = pd.DataFrame(columns=('label', 'area_nucleus', 'area_cytoplasm', 'area_infection',
                      'mean_int_nucleus', 'mean_int_cytoplasm', 'mean_int_infection',
                      'mean_int_mito_signal', 'std_int_nucleus', 'std_int_cytoplasm',
                      'std_int_infection', 'std_int_mito_signal', 'centroid-0', 'centroid-1'))
    
    for idx, row in nuc_df.iterrows():
        label = idx+1
        c = labels == label
        n = nuc_lbl == row.label   
        
        df.loc[idx] = [label, int(row.area), np.sum(c), np.sum(inf & c),
                       np.mean(img[0][n]), np.mean(img[3][c^n]), np.mean(img[1][n]),
                       np.mean(img[2][c^n]), np.std(img[0][n]), np.std(img[3][c^n]),
                       np.std(img[1][n]), np.std(img[2][c^n]), row['centroid-0'], row['centroid-1']
                       ]
    
    #save as pickle with a cleaner name than the image 
    df.to_pickle(f.replace('_XY1.ome.tif', '.pckl'))

