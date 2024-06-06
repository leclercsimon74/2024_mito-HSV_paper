# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 14:59:50 2023

@author: sleclerc
"""

import open3d as o3d
import tifffile
from skimage import measure
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import ndimage
import matplotlib


spacing = (30,8,8) #in nm

nucleus_color = np.array([105, 172, 234])
mito_color = np.array([228, 159, 102])
er_max_distance = 700 # in nm, need to check other!
contact_size_dist = 30 #in nm


mito_path = "Segmented_Mitochondria_Labels_MEF_WT_uninfected_220720_R3_8x8x30_85-sub.tif"
nucl_path = "Segmented_Nucleus_Labels_MEF_WT_uninfected_220720_R3_8x8x30_85-sub.tif"
endr_path = "Segmented_Reticulum_Labels_MEF_WT_uninfected_220720_R3_8x8x30_85-sub.tif"
ori_path = 'MEF_WT_uninfected_220720_R3_8x8x30_85-sub.tif' #<Too big for Github!
ori_path = 'MEF_WT_uninfected_220720_R3_8x8x30_85-sub-slice7.tif' #<Just one slice

def extend(img):
    zeros = np.zeros(np.array(img.shape)+2, dtype=bool)
    zeros[1:-1,1:-1,1:-1] = img
    return zeros


def add_mesh(img, color):
    verts, faces, normals, values = measure.marching_cubes(img*255, 1,
                                                           step_size=1,
                                                           allow_degenerate=False,
                                                           spacing=spacing)

    #create the mito mesh in open3D
    mesh = o3d.geometry.TriangleMesh()
    mesh.triangles = o3d.utility.Vector3iVector(faces)
    mesh.vertices = o3d.utility.Vector3dVector(verts)
    mesh = mesh.filter_smooth_simple(number_of_iterations=1)
    mesh.compute_vertex_normals()
    mesh.paint_uniform_color(color/255)
    
    return mesh




mito_img = tifffile.imread(mito_path)
nucl_img = tifffile.imread(nucl_path)
endr_img = tifffile.imread(endr_path)

mito_img = np.where(mito_img==1, True, False)
nucl_img = np.where(nucl_img==1, True, False)
endr_img = np.where(endr_img==1, True, False)

#extend to avoid edge effect
mito_img = extend(mito_img)
nucl_img = extend(nucl_img)
endr_img = extend(endr_img)

#%%
meshes = []

#meshes.append(add_mesh(mito_img, mito_color))
meshes.append(add_mesh(nucl_img, nucleus_color))
meshes.append(add_mesh(endr_img, np.array([175,175,175])))

#mito_dst = sc.ndimage.distance_transform_edt(np.invert(mito_img), sampling=spacing)
endr_dst = sc.ndimage.distance_transform_edt(np.invert(endr_img), sampling=spacing)

verts, faces, normals, values = measure.marching_cubes(mito_img*255, 1,
                                                       step_size=2,
                                                       allow_degenerate=False,
                                                       spacing=spacing)

#create the mesh in open3D
mesh = o3d.geometry.TriangleMesh()
mesh.triangles = o3d.utility.Vector3iVector(faces)
mesh.vertices = o3d.utility.Vector3dVector(verts)
mesh.compute_vertex_normals()


vert = np.asarray(mesh.vertices)
col = []
for idx, v in enumerate(vert):
    vv = (np.round(v)//np.array(spacing)).astype(int)
    c = endr_dst[vv[0], vv[1], vv[2]]
    if c < contact_size_dist: #contact site, distance in nm
        c = 0
    col.append(c)

col = np.array(col)
print(np.max(col))
#save the distance
np.save('NI_dist_map.pckl', col)
col = np.abs((col - np.min(col))/(er_max_distance - np.min(col)) - 1) #normalize and invert
i = np.where(col == 1) #found contact site
col = cm.viridis(col)
col = col[:,0:3] #remove transparency channel
col[i[0]] = np.array([1,0,0]) #assign color red to the contact site

mesh.vertex_colors = o3d.utility.Vector3dVector(col.astype(float))
meshes.append(mesh)

o3d.visualization.draw_geometries(meshes)


#%% Draw figure with contour

contour_thickness = 5
s = 7 #slice
enr_color = np.array([175, 175, 200]) #RGB


#ori_img = tifffile.imread(ori_path)
#zeros = np.zeros(np.array(ori_img.shape)+2, dtype='int8')
#zeros[1:-1,1:-1,1:-1] = ori_img
#ori_img = zeros
#del zeros
#ori_img = np.where(ori_img < 0, np.max(ori_img)+np.abs(ori_img), ori_img) #wrapping around...

ori_img = tifffile.imread(ori_path)

from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar



fig = plt.figure(figsize=(30,15))
ax = fig.add_subplot(1,2,1)
scalebar = AnchoredSizeBar(ax.transData,
                           500/spacing[1], '', 'lower right', 
                           borderpad = 5,
                           color='black',
                           frameon=False,
                           size_vertical=5,
                           fontproperties={'size':0})

ax.imshow(ori_img, cmap='gray_r')


ax.add_artist(scalebar)
plt.axis('off')

ax = fig.add_subplot(1,2,2)
ax.imshow(ori_img[s], cmap='gray_r')

cmap1 = plt.cm.viridis_r
mito_c = mito_img[s] ^ ndimage.binary_erosion(mito_img[s], iterations = contour_thickness)
#mito_c = mito_img[s]
mito_c = np.where(mito_c, endr_dst[s], np.nan)
im = ax.imshow(mito_c, cmap = 'viridis_r', vmin=contact_size_dist, vmax=er_max_distance)
im.cmap.set_under('r')

nucl_c = nucl_img[s] ^ ndimage.binary_erosion(nucl_img[s], iterations = contour_thickness)
nucl_c = np.where(nucl_c, nucl_img[s], np.nan)
ax.imshow(nucl_c, cmap=matplotlib.colors.ListedColormap([nucleus_color/255, 'black']))

endr_c = endr_img[s] ^ ndimage.binary_erosion(endr_img[s], iterations = contour_thickness)
endr_c = np.where(endr_c, endr_img[s], np.nan)
ax.imshow(endr_c, cmap=matplotlib.colors.ListedColormap([enr_color/255, 'black']))
scalebar = AnchoredSizeBar(ax.transData,
                           500/spacing[1], '', 'lower right', 
                           borderpad = 5,
                           color='black',
                           frameon=False,
                           size_vertical=5,
                           fontproperties={'size':0})
ax.add_artist(scalebar)

plt.axis('off')

plt.tight_layout()
plt.savefig('NI_3DView.svg')
plt.savefig('NI_3DView.png', dpi=300)
plt.show()

#%% colorbar

cmap = matplotlib.colormaps.get_cmap('viridis_r') #from 0 to 1!
norm = matplotlib.colors.Normalize(vmin=contact_size_dist, vmax=er_max_distance) #normalize

fig,ax = plt.subplots(figsize=(3,0.3))
cbar = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                                        orientation='horizontal')

plt.xticks([0, 30, 100, 300, 600])
plt.savefig("colorbar.svg", bbox_inches='tight')
plt.show



#%% Grab the number of contact site/mito and area of contact site
from skimage.morphology import binary_erosion

mito_l = measure.label(mito_img) #grab by mitochondria

n = 0
area_contact_site = []
contact_site_n = []
mito_SA = []
for i in range(1, np.max(mito_l)+1, 1):
    n += 1
    mito = np.where(mito_l == i, True, False)

    #Calculate the surface are of each mitochondria
    verts, faces, normals, values = measure.marching_cubes(mito*255, 1,
                                                           step_size=1,
                                                           allow_degenerate=False,
                                                           spacing=spacing)
    mito_SA.append(measure.mesh_surface_area(verts, faces)) #in nm
    
    #calculate the number of contact site
    mito_d = np.where(mito, endr_dst, np.nan)
    mito_d = np.where(mito_d < contact_size_dist, True, False)
    mito_lbl = measure.label(mito_d)
    if np.max(mito_lbl) == 0: #in case of NO contact site on this mitochondria
        contact_site_n.append(0)
        continue
    contact_site_n.append(np.max(mito_lbl))
    
    #calculate the surface area of the contact site
    #get the edge only
    edge = binary_erosion(mito)
    edge = mito ^ edge
    #grab the distance
    mito_d = np.where(edge, endr_dst, np.nan)
    mito_d = np.where(mito_d < contact_size_dist, True, False)
    mito_lbl = measure.label(mito_d)
    ids, count = np.unique(mito_lbl, return_counts=True)
    #append the volume in contact with ER
    for l in range(len(count)):
        if l == 0: #background
            continue
        if count[l] > 10: #minimum contact size in voxel
            #measure the contact site area
            contact_site = np.where(mito_lbl==l, True, False)
            verts, faces, normals, values = measure.marching_cubes(contact_site*255, 1,
                                                                   step_size=1,
                                                                   allow_degenerate=False,
                                                                   spacing=spacing)
            #divide by two since it takes both side
            area_contact_site.append(measure.mesh_surface_area(verts, faces)/2) #in nm
    
    

#save all lists
np.save('ni_area_contact_site', np.array(area_contact_site))
np.save('ni_contact_site_n', np.array(contact_site_n))
np.save('SA_mito', np.array(mito_SA))
