# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 16:15:03 2023

@author: sleclerc


Cristae comparison between NI and 8 hpi

Segmentation in tif with value:
- 0 is background
- 1 is mitochondria
- 2 is cristae
- 3 (only for the 1st image) is the mitochondria outer membrane

3D volume
Warning, mutually exculsive, but cristae is of course in the mitochondria!
Need to do a fill hole slice by slice (2D) - some cristae touch the top/bottom of the selection.

Analysis:
- Ratio Volume cristae / Volume Mito - swarmplot/violin
- Cristae 3D shape descriptor (thickness, shape) - swarmplot/violin
- Cristae orientation inside mitochondria - swarmplot/violin
- Surface area of cristae (full, no subseparation)
- lenght of cristae

Using the knowledge from the cristae average thickness, we can decide on a
 threshold to take in account only the inflated cristae

"""

import os
import tifffile
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy as sc
from scipy import stats
from scipy import ndimage as ndi
import localthickness as lcthk
import pandas as pd
from skimage.measure import label, regionprops_table, mesh_surface_area, marching_cubes
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from statannotations.Annotator import Annotator
from skimage import draw
from skimage.morphology import binary_erosion
import matplotlib as mpl
mpl.rcParams['mathtext.default'] = 'regular'

import utils as u

def calculate_local_thickness(array):
    edt = lcthk.distance_transform_edt(array)
    ridge = lcthk.distance_ridge(edt)
    local_thickness = lcthk.local_thickness(edt, ridge) #in pixel!
    local_thickness = np.max(local_thickness[local_thickness > 0])
    return local_thickness

def measure2D (binary):
    binary = label(binary, background = 0)
    to_measure = ['label', 'area', 'centroid', 'convex_area',
                  'major_axis_length', 'minor_axis_length', 'orientation']
    
    props = pd.DataFrame.from_dict(regionprops_table(binary, properties=to_measure))
    return props

def angle_between_v1_v2(vector1, vector2):
    # Calculate the dot product of the two vectors
    dot_product = np.dot(vector1, vector2)
    
    # Calculate the magnitudes (norms) of the two vectors
    magnitude1 = np.linalg.norm(vector1)
    magnitude2 = np.linalg.norm(vector2)
    
    # Calculate the cosine of the angle between the two vectors
    cosine_theta = dot_product / (magnitude1 * magnitude2)
    
    # Calculate the angle in radians
    angle_radians = np.arccos(cosine_theta)
    
    # Convert the angle from radians to degrees
    return np.degrees(angle_radians)

def found_intersection(mask, center_point, vector):
    # Define your starting position (x0, y0) and direction vector (dx, dy)
    x0, y0 = int(center_point[0]), int(center_point[1])  # Starting position (inside the mask)
    # Direction vector
    dx, dy = int(vector[0]), int(vector[1])
    # Initialize step size
    step_size = max(abs(dx), abs(dy))

    # Normalize the direction vector
    dx /= step_size
    dy /= step_size

    # Iterate along the vector to find when it exits the mask in BOTH direction
    exit_x, exit_y = 0, 0
    x1, y1 = x0, y0
    x2, y2 = x0, y0
    i, ii = False, False
    while True:
        if 0 <= y1 < mask.shape[0] and 0 <= x1 < mask.shape[1]:
            # Check if the current point in the mask is False (exiting)
            if not mask[int(y1), int(x1)]:
                exit_x, exit_y = x1, y1
                break  # Stop the loop when the vector exits the mask
            # Move along the vector
            x1 += dx
            y1 += dy
        else:
            i = True
            
        if 0 <= y2 < mask.shape[0] and 0 <= x2 < mask.shape[1]:
            # Check if the current point in the mask is False (exiting)
            if not mask[int(y2), int(x2)]:
                exit_x, exit_y = x2, y2
                break  # Stop the loop when the vector exits the mask
            # Move along the vector
            x2 -= dx
            y2 -= dy
        else:
            ii = True
        
        if i and ii:
            break
    
    return int(exit_x), int(exit_y)
    
def get_tang_mask(edge, intersection, radius):
    #circle
    rr, cc = draw.circle_perimeter(intersection[1], intersection[0], 5)
    zeros = np.zeros(edge.shape)
    zeros[rr,cc] = 1
    
    return np.where(zeros + edge == 2)


def extract_data(filename):
    thickness = []
    orientation = []
    lenght = []
    surface_area = []

    img = np.array(tifffile.imread(filename))
    mito = np.where(img>0, True, False)
    cristae = np.where(img==2, True, False)
    del img
    
    #quick visu
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1,2,1)
    ax.imshow(np.sum(mito, axis=0), interpolation=None)
    ax.set_title('Mitochondria')
    plt.axis('off')
    ax = fig.add_subplot(1,2,2)
    ax.imshow(np.sum(cristae, axis=0), interpolation=None)
    ax.set_title('Cristae')
    plt.axis('off')
    plt.tight_layout()
    plt.show()
    
    #volume ratio
    ratio = np.sum(cristae)/np.sum(mito) #in pixel!
    
    #watershed
    distance = ndi.distance_transform_edt(cristae)
    coords = peak_local_max(distance, footprint=np.ones((3, 3, 3)),
                            num_peaks=27) #TODO, check peak number! 30 better for the moment
    #peak num: 20 good, 27 very good,30 very good, 33 very good, 40 good
    mask = np.zeros(distance.shape, dtype=bool)
    mask[tuple(coords.T)] = True
    markers, _ = ndi.label(mask)
    labels = watershed(-distance, markers, mask=cristae)
    
    print('Number of labels is: '+str(np.max(labels)))
    
    #orientation only in 2D!
    mito_2D = np.where(np.sum(mito, axis=0)>0, True, False) #flatten mito
    mito_prop = measure2D(mito_2D)
    mito_prop = mito_prop[mito_prop.area == mito_prop.area.max()].squeeze()
   
    for i in range(np.max(labels)+1):
        if i == 0: continue #background
        new_sub = np.where(labels==i, True, False)
        #print('Cristae volume: '+str(np.sum(new_sub)))
        if np.sum(new_sub) < 500: continue #size limit
        #thickness
        thickness.append(calculate_local_thickness(new_sub))

        #surface area
        verts, faces, normals, values = marching_cubes(new_sub*255, 1,
                                                       step_size=1,
                                                       allow_degenerate=False)
        surface_area.append(mesh_surface_area(verts, faces))

        cristae_2D = np.where(np.sum(new_sub, axis=0)>0, True, False) #flatten mito
        cristae_prop = measure2D(cristae_2D)
        #length
        lenght.append(cristae_prop['major_axis_length'] * 2)
        
        #orientation OLD!
        # angle1 = abs(np.rad2deg(mito_prop['orientation']-cristae_prop['orientation']).values[0])
        # angle2 = abs(180 - abs(np.rad2deg(mito_prop['orientation']-cristae_prop['orientation']).values[0]))
        
        #orientation.append(min([angle1, angle2]))
        
        #orientation, NEW!
        center_point = [int(cristae_prop['centroid-1']), int(cristae_prop['centroid-0'])]
        # Direction vector 
        dx = int(np.sin(cristae_prop['orientation'])*cristae_prop['major_axis_length']/2+cristae_prop['centroid-1'])
        dy = int(np.cos(cristae_prop['orientation'])*cristae_prop['major_axis_length']/2+cristae_prop['centroid-0'])

        intersection = found_intersection(mito_2D, center_point, [dx, dy])

        #edge of the mitochondria
        edge = binary_erosion(mito_2D)
        edge = mito_2D ^ edge
        
        try:
            inter_x, inter_y = get_tang_mask(edge, intersection, 5)
        except IndexError:
            try:
                inter_x, inter_y = get_tang_mask(edge, intersection, 2)
            except IndexError:
                continue
        if len(inter_x) != 2 and len(inter_y) != 2: #bad vector
            continue
        #vector to 0
        v1 = [center_point[0]-dx, center_point[1]-dy]
        v2 = [inter_x[0]-inter_x[1], inter_y[0]-inter_y[1]]
        

        angle = angle_between_v1_v2(v2, v1)
        orientation.append(min(180-angle, angle))
    
    #quick visu last cristae
    print('Angle is '+str(min(180-angle, angle)))
    plt.figure()
    plt.imshow(cristae_2D + edge)
    plt.scatter(center_point[0], center_point[1], c='r')
    plt.plot([center_point[0],intersection[0]], [center_point[1],intersection[1]], c='r')
    plt.scatter(intersection[0], intersection[1], c='g')
    plt.plot(inter_y, inter_x, c='r')
    plt.title(os.path.basename(filename))
    plt.axis('off')
    plt.show()
    
    return ratio, surface_area, thickness, orientation, lenght


spacing = [5,5,5] #nm

path_NI = 'cristae NI'
path_INF8 = 'cristae 8hpi'

img_path_NI = u.recursive_found_img(path_NI)
img_path_INF8 = u.recursive_found_img(path_INF8)

#%% measure

NI_ratio_volume = []
NI_thickness = []
NI_orientation = []
NI_SA = []
NI_lenght = []

INF_ratio_volume = []
INF_thickness = []
INF_orientation = []
INF_SA = []
INF_lenght = []

for filename in img_path_NI:
    if 'raw' in filename:
        continue
    ratio, SA, thickness, orientation, lenght = extract_data(filename)
    NI_ratio_volume.append(ratio)
    NI_SA += SA
    NI_thickness += thickness
    NI_orientation += orientation
    NI_lenght += lenght
    
    
for filename in img_path_INF8:
    if 'raw' in filename:
        continue
    ratio, SA, thickness, orientation, lenght = extract_data(filename)
    INF_ratio_volume.append(ratio)
    INF_SA += SA
    INF_thickness += thickness
    INF_orientation += orientation
    INF_lenght += lenght


#%% Make Figure
def make_swarmplot(data, ax, dotsize):
    sns.boxplot(data=data, ax=ax, 
                      showcaps=False, boxprops={'facecolor':'None'},
                      showfliers=False, whiskerprops={'linewidth':0},
                      medianprops={'linewidth':0}, showmeans=True,
                      meanline=True, meanprops={'color':'black'})
    sns.swarmplot(data=data, size=dotsize, ax=ax, zorder=0)
    annotator = Annotator(ax, pairs, data=data)
    annotator.configure(test=test, verbose=True).apply_and_annotate()

def add_scalebar(pxl_size, ax, scale_bar_size=500):
    return AnchoredSizeBar(ax.transData,
                           scale_bar_size/pxl_size, '', 'lower right', 
                           borderpad = 5,
                           color='black',
                           frameon=False,
                           size_vertical=3,
                           fontproperties={'size':0})

def get_edge(arr):
    #get the edge of a binary image
    edge = binary_erosion(arr)
    return arr ^ edge

color = u.get_color(4, 'original')
sns.set_palette(sns.color_palette(color, len(color)))

pairs = [['NI', '8 hpi']]
test = 't-test_ind'
#test = 'Brunner-Munzel'
from PIL import Image
import tifffile
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
#%%
div = 2
fig = plt.figure(figsize=(15/div, 15/div))
#part A is raw data with segmentation (edge only)
#open tif image of the original


cristae_c = np.array((90,200,120))
outer_m_c = np.array((215,215,215))


ax = fig.add_subplot(3,2,1)
path_ni_raw = 'cristae NI\m1_NI_raw.tif'
path_ni_seg = 'cristae NI\m1_NI.tif'
z = 12
img = tifffile.imread(path_ni_raw)
seg = tifffile.imread(path_ni_seg)
cristae = get_edge(np.where(seg[z]==2, True, False))
outer_m = get_edge(np.where(seg[z]==3, True, False))
ax.imshow(u.img_color_merger(gr=img[z], g=cristae*255, y=outer_m*255))
ax.add_artist(add_scalebar(5, ax))
plt.axis('off')

ax = fig.add_subplot(3,2,3)
path_inf_raw = 'cristae 8hpi\m1_8hpi_raw.tif'
path_inf_seg = 'cristae 8hpi\m1_8hpi.tif'
z = 12
img = tifffile.imread(path_inf_raw)
seg = tifffile.imread(path_inf_seg)
cristae = get_edge(np.where(seg[z]==2, True, False))
outer_m = get_edge(np.where(seg[z]==3, True, False))
ax.imshow(u.img_color_merger(gr=img[z], g=cristae*255, y=outer_m*255))
ax.add_artist(add_scalebar(5, ax))
plt.axis('off')



#part B is 3D reconstruction
path_3D_ni = 'DragonFly and Blender\FIB-SEM_NI.png'
path_3D_inf = 'DragonFly and Blender\FIB-SEM_8hpi.png'

ax = fig.add_subplot(3,2,2)
ax.imshow(Image.open(path_3D_ni))
plt.axis('off')

ax = fig.add_subplot(3,2,4)
ax.imshow(Image.open(path_3D_inf))
plt.axis('off')

#Part C is analysis


# ax = fig.add_subplot(3,3,7)
# dotsize = 10
# data = {'NI':NI_ratio_volume,
#         '8 hpi':INF_ratio_volume}
# data = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in data.items()]))
# make_swarmplot(data, ax, dotsize)
# ax.set_title('Ratio cristae/mitochondria volume')
# ax.set_ylabel('Ratio cristae/mitochondria volume')


ax = fig.add_subplot(3,3,7)
dotsize = 3
data = {'NI':np.array(NI_thickness)[np.abs(stats.zscore(NI_thickness)) < 3],
        '8 hpi':np.array(INF_thickness)[np.abs(stats.zscore(INF_thickness)) < 3]}
data = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in data.items()]))
data = data*5
data.to_csv('csv'+os.sep+'Max thickness.csv')
make_swarmplot(data, ax, dotsize) #adjsut to the spacing
ax.set_title('Cristae thickness')
ax.set_ylabel('Max. Thickness (nm)')

# ax = fig.add_subplot(1,5,3)
# data = {'NI':np.array(NI_orientation)[np.abs(stats.zscore(NI_orientation)) < 3],
#         '8 hpi':np.array(INF_orientation)[np.abs(stats.zscore(INF_orientation)) < 3]}
# data = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in data.items()]))
# make_swarmplot(data, ax, dotsize)
# ax.set_title('Cristae orientation')
# ax.set_ylabel('Angle (degree)')


ax = fig.add_subplot(3,3,8)
data = {'NI':np.array(NI_SA)[np.abs(stats.zscore(NI_SA)) < 3],
        '8 hpi':np.array(INF_SA)[np.abs(stats.zscore(INF_SA)) < 3]}
data = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in data.items()]))
data = data*(5*5)/1000
data.to_csv('csv'+os.sep+'Surface area.csv')
make_swarmplot(data, ax, dotsize)
ax.set_title('Cristae Surface Area')
ax.set_ylabel('Surface Area (${\mu m^2}$)')

ax = fig.add_subplot(3,3,9)
data = {'NI':np.array(NI_lenght)[np.abs(stats.zscore(NI_lenght)) < 3],
        '8 hpi':np.array(INF_lenght)[np.abs(stats.zscore(INF_lenght)) < 3]}
data = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in data.items()]))
data = data*(5)/1000
data.to_csv('csv'+os.sep+'Length.csv')
make_swarmplot(data, ax, dotsize)
ax.set_title('Cristae length')
ax.set_ylabel('Cristae length (${\mu m}$)')


plt.tight_layout()
plt.savefig('Figure.png', dpi=600)
plt.savefig('Figure.svg')
plt.show()
