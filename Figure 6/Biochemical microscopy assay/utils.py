# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 09:44:14 2021

@author: sleclerc
Last updated one

"""

import numpy as np
import scipy as sc
from scipy.ndimage import morphology
from scipy.stats import gaussian_kde
from skimage import filters, measure
import re
import tifffile
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import itertools
import pickle


#image function
def get_img(filename):
    img = np.array(tifffile.imread(filename))
    if np.min(img) < 0:
        img = np.where(img<0, img+256, img)
    return img


def found_img(path, extension='.tif'):
    """
    Found all images with the file extension. Can remplace the file extension
    with a search string, such as 'mock cell 001' will found all image with
    this in its name /b in the folder b/. This is not recursive, meaning that 
    it will not search in subfolder.
    Parameters
    ----------
    path : str
        local or global string path
    extension : str, optional
        A search string option. The default is '.tif'. '' will deactivate it.

    Returns
    -------
    img_path : list
        a list of string indicate where to found each selected img

    """
    imgs = [f for f in os.listdir(path) if f.endswith(extension)]
    imgs_path = [path+"\\"+f for f in imgs]
    return imgs_path


def recursive_found_img(path, extension='.tif'):
    """
    Like found_img, but recursive.

    Parameters
    ----------
    path : str
        local or global string path
    extension : str, optional
        A search string option. The default is '.tif'. '' will deactivate it.

    Returns
    -------
    l : list
        a list of string indicate where to found each selected img

    """
    l = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(extension):
                l.append(root+os.sep+file)
    return l


def open_img(filename, img_list):
    """
    Using the filename, open all images sharing its name (so all channels).
    If the image is a grayscale image, open only this image
    Parameters
    ----------
    filename : str
        Name of the image to open (from get_img_path_name)
    img_list : list of string
        List of images path (from found_img)

    Returns
    -------
    data : numpy array
        Concatenate all channels in one image

    """
    imgs = [f for f in img_list if os.path.basename(filename) in os.path.basename(f)]
    C_str = ["C0", "C1", "C2", "C3", "C4", "C5"]
    data = []
    for c in C_str:
        for img in imgs:
            if c in img:
                data.append(get_img(img))
        
    data = np.array(data)
    
    if len(data) == 0: #empty, image is likely a grayscale img (alone)
        data = get_img(filename)
        
    if data.shape[-1] < 4: #order in axis not good
        data = np.fliplr(np.rot90(np.swapaxes(data,0,-1)[0], k=3))
            
    return data


def get_img_path_name(imgs_path, extension='.tif'):
    """
    Get to group the image by their channel.
    Rely on the name, like C0 to C5.

    Parameters
    ----------
    imgs_path : list of string
        gotten from found_img
    extension : str, optional
        Image extension. The default is '.tif'.

    Returns
    -------
    img_name : list of str
        str containing the global name for each image independant of the channel

    """
    C_str = ["C0", "C1", "C2", "C3", "C4", "C5"]
    img_name = []
    for f in imgs_path:
        i = False
        if i: break # found this image as a multichannel image
        for c in C_str:
            if c in f:
                #suppose that the C# is at the beggining or end of the name
                f = f.replace(c, '')
                f = f.replace(extension, '')
                if f not in img_name: #to avoid duplicate
                    img_name.append(f)
                    i = True
        if f not in img_name and i == False: #not channel in name, assume a grayscale img
            img_name.append(f)
            
    return img_name


def get_img_full_path(img_path, extension='.tif'):
    """do the inverse of get_img_path_name"""
    C_str = ["C0", "C1", "C2", "C3", "C4", "C5"]
    path = os.path.dirname(img_path)
    name = os.path.basename(img_path)
    name = name.replace(extension, '')
    if name.startswith(tuple(C_str)):
        name = name[2:]
    elif name.endswith(tuple(C_str)):
        name = name[:-2]
    img_channel = []
    for c in C_str: #multichannel
        if c in name:
            img_channel.append(path+os.sep+name+extension)
        elif os.path.isfile(path+os.sep+c+name+extension):
            img_channel.append(path+os.sep+c+name+extension)
        elif os.path.isfile(path+os.sep+name+c+extension):
            img_channel.append(path+os.sep+name+c+extension)
    if os.path.isfile(path+os.sep+name+extension): #grayscale
        img_channel.append(path+os.sep+name+extension)
    
    return img_channel


def rand_dist_mask(binary, mask, graph=False):
    """
    Function that randomize the position of positive pixel from a binary image
    inside a mask image. Somewhat quick. Work for 3D image, should work for 2D
    and more than 3D as well. The binary and mask should have the same dimension.

    Parameters
    ----------
    binary : np.array
        The array to randomize. 3D numpy array.
    mask : np.array
        The mask to apply. 3D numpy array, same shape as binary.
    graph : bool, optional
        If display a average z projection. The default is False.

    Returns
    -------
    random_array : np.array
        Same shape than mask. bool type.

    """
    
    #extract the working pixel (generally the nucleus)
    s = np.array(np.where(mask))
    s = s.T #transpose for easy selection
    # randomly choose some pixel from the working pixel, same size than the ROI
    new = np.random.choice(np.arange(len(s)), size=len(np.where(binary>0)[0]))
    # replace choosen pixel
    random_array = np.zeros(binary.shape, dtype=bool) #in a new array
    new = s[new] #select
    new = tuple(new.T) #transpose and transform in tuple for quick insertion
    random_array[new] = True #change the pixel value at the given position
    if graph:
        if len(random_array.shape) == 3:
            plt.imshow(np.mean(random_array, axis=0), interpolation='none')
        elif len(random_array.shape) == 2:
            plt.imshow(random_array, interpolation='none')
        plt.title('Random distribution of label in ROI')
        plt.axis('off')
        plt.tight_layout()
        plt.show()
    
    return random_array

def density(arr, variance=0.15): #1D
    data = np.ndarray.flatten(arr)
    data = data[data != -1]
    data_sorted = np.sort(data)
    
    density = gaussian_kde(data_sorted)
    density.covariance_factor = lambda : variance
    density._compute_covariance()
    
    return density


def density_estimation(m1, m2): #2D
    X, Y = np.mgrid[np.min(m1):np.max(m1):100j, np.min(m2):np.max(m2):100j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                               
    values = np.vstack([m1, m2])                                                                        
    kernel = gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z


def norm(arr):
    return (arr - np.min(arr))/np.ptp(arr)


def dst_map(label1, label2, ROI, res, graph=False):
    """
    The minimum distance between pixel of label1 from label2

    Parameters
    ----------
    label1 : np.array
        Bool array in 3D containing label1 localization.
    label2 : np.array
        Bool array in 3D containing label2 localization.
    ROI : np.array
        Bool array in 3D containing the ROI to focus as True.
        
    Returns
    -------
    dis : np.array
        The distance map of label1 to label2.
    rd_dis : np.array
        A random distance map if Label1 were dispersed inside the ROI.

    """
    #binary transform: all positive signal are positive!
    label1 = np.array(np.where(label1>0, 1, 0), dtype=bool)
    label2 = np.array(np.where(label2>0, 1, 0), dtype=bool)
    ROI = np.array(np.where(ROI>0, 1, 0), dtype=bool)
    
    #mask to take in account the distance to another object - negative
    mask = np.array(np.where(label2>0, 0, 1), dtype=bool)
    #calculate the distance
    dis = sc.ndimage.distance_transform_edt(mask, sampling=res)
    #keep distance value of the object of interest only
    dis = np.where(label1 == 1, dis, -1)
    
    #random distribution of label accross the nucleus
    rd = rand_dist_mask(label2, ROI)
    rd_neg = np.array(np.where(rd>0, 0, 1), dtype=bool) #invert it!
    rd_dis = sc.ndimage.distance_transform_edt(rd_neg, sampling=res)
    rd_dis = np.where(label1 == 1, rd_dis, -1)
    
    
    if graph:
        fig, ax = plt.subplots(figsize=(10, 5) , nrows=1, ncols=2)
        ax[0].imshow(dis, interpolation='none')
        ax[0].set_title('Label1 to Label2')
        ax[0].set_yticks([])
        ax[0].set_xticks([])
        ax[1].imshow(rd_dis, interpolation='none')
        ax[1].set_title('Label1 to random label2')
        ax[1].set_yticks([])
        ax[1].set_xticks([])
        plt.tight_layout()
        plt.show()

    return dis, rd_dis


def dist_correlation(label_seg, label_name, nucl_binary, img, res, channel, label_to_img, graph=True):
    label_combi = list(itertools.combinations(label_seg, 2))
    name_combi = list(itertools.combinations(label_name, 2))
    
    dist_data = {}
    intensity_data = {}
    
    #calculation
    for idx, label in enumerate(label_combi):
        name = name_combi[idx]
        
        #correspondance image channel from segmentation name for img intensity 
        if name[0] in channel:
            n0 = name[0]
            t0 = img[channel[name[0]]-1]
        else: #other name/treatment done
            n0 = label_to_img[label_to_img.index(label_to_img[label_name.index(name[0])])]
            t0 = img[channel[label_to_img[label_to_img.index(label_to_img[label_name.index(name[0])])]]-1]
        
        if name[1] in channel: 
            n1 = name[1]
            t1 = img[channel[name[1]]-1]
        else: #other name/treatment done
            n1 = label_to_img[label_to_img.index(label_to_img[label_name.index(name[1])])]
            t1 = img[channel[label_to_img[label_to_img.index(label_to_img[label_name.index(name[1])])]]-1]
        
        #Calculation for label1=>label2 direction
        dis, rd_dis = dst_map(label[0], label[1], nucl_binary, res)
        data = np.ndarray.flatten(dis)
        data = data[data != -1]
        rd_data = np.ndarray.flatten(rd_dis)
        rd_data = rd_data[rd_data != -1]
        dist_data[name[0]+'-'+name[1]] = data
        dist_data[name[0]+'-'+name[1]+'_rd'] = rd_data      
        
        intensity = np.where(dis>=0, t0, -1)
        intensity = np.ndarray.flatten(intensity)
        intensity = intensity[intensity != -1]
        rd_intensity = np.where(rd_dis>=0, t0, -1)
        rd_intensity = np.ndarray.flatten(rd_intensity)
        rd_intensity = rd_intensity[rd_intensity != -1]        
        intensity_data[name[0]+'-'+name[1]] = intensity
        intensity_data[name[0]+'-'+name[1]+'_rd'] = rd_intensity
        
        #Calculation for label2=>label1 direction
        dis, rd_dis = dst_map(label[1], label[0], nucl_binary, res)
        data = np.ndarray.flatten(dis)
        data = data[data != -1]
        rd_data = np.ndarray.flatten(rd_dis)
        rd_data = rd_data[rd_data != -1]
        dist_data[name[1]+'-'+name[0]] = data
        dist_data[name[1]+'-'+name[0]+'_rd'] = rd_data
        
        intensity = np.where(dis>=0, t1, -1)
        intensity = np.ndarray.flatten(intensity)
        intensity = intensity[intensity != -1]
        rd_intensity = np.where(rd_dis>=0, t1, -1)
        rd_intensity = np.ndarray.flatten(rd_intensity)
        rd_intensity = rd_intensity[rd_intensity != -1]        
        intensity_data[name[1]+'-'+name[0]] = intensity
        intensity_data[name[1]+'-'+name[0]+'_rd'] = rd_intensity        

    
    if graph:
        ncol = 6
        nrow = len(label_combi)
        fig, ax = plt.subplots(figsize=(5.0*ncol, 5.0*nrow) , nrows=nrow, ncols=ncol)

        for row in range(nrow):
            name = name_combi[row]
            label = label_combi[row]
            #correspondance image channel from segmentation name for img intensity 
            if name[0] in channel:
                n0 = name[0]
                t0 = img[channel[name[0]]-1]
            else: #other name/treatment done
                n0 = label_to_img[label_to_img.index(label_to_img[label_name.index(name[0])])]
                t0 = img[channel[label_to_img[label_to_img.index(label_to_img[label_name.index(name[0])])]]-1]
            if name[1] in channel: 
                n1 = name[1]
                t1 = img[channel[name[1]]-1]
            else: #other name/treatment done
                n1 = label_to_img[label_to_img.index(label_to_img[label_name.index(name[1])])]
                t1 = img[channel[label_to_img[label_to_img.index(label_to_img[label_name.index(name[1])])]]-1]
        
        
            for col in range(ncol):
                if col == 0: #one side
                    dis, rd_dis = dst_map(label[0], label[1], nucl_binary, res)
                    if len(dis.shape) == 3: 
                        RGB = np.dstack(np.array((norm(np.mean(dis + 1, axis=0)),
                                                  norm(np.mean(label[1], axis=0)),
                                                  np.zeros((np.mean(dis + 1, axis=0).shape), dtype=float))))
                    elif len(dis.shape) == 2: 
                        RGB = np.dstack(np.array((dis + 1,
                                                  label[1],
                                                  np.zeros((dis.shape), dtype=float))))
                        
                    ax[row, col].imshow(RGB, interpolation='none')
                    ax[row, col].set_yticks([])
                    ax[row, col].set_xticks([])
                
                elif col == 3: #other side
                    dis, rd_dis = dst_map(label[1], label[0], nucl_binary, res)       
                    if len(dis.shape) == 3: 
                        RGB = np.dstack(np.array((norm(np.mean(dis + 1, axis=0)),
                                                  norm(np.mean(label[0], axis=0)),
                                                  np.zeros((np.mean(dis + 1, axis=0).shape), dtype=float))))
                    elif len(dis.shape) == 2: 
                        RGB = np.dstack(np.array((dis + 1,
                                                  label[0],
                                                  np.zeros((dis.shape), dtype=float))))
                    
                    ax[row, col].imshow(RGB, interpolation='none')
                    ax[row, col].set_yticks([])
                    ax[row, col].set_xticks([])           
                
                if col == 1 or col == 4: #CDF
                    data = np.ndarray.flatten(dis)
                    data = data[data != -1]
                    # sort the data:
                    data_sorted = np.sort(data)
                    # calculate the proportional values of samples
                    p = 1. * np.arange(len(data)) / (len(data) - 1)
                    ax[row, col].plot(data_sorted, p, 'b', label='real') #plot
                    
                    rd_data = np.ndarray.flatten(rd_dis)
                    rd_data = rd_data[rd_data != -1]
                    # sort the data:
                    rd_data_sorted = np.sort(rd_data)
                    # calculate the proportional values of samples
                    rd_p = 1. * np.arange(len(rd_data)) / (len(rd_data) - 1)
                    ax[row, col].plot(rd_data_sorted, rd_p, 'c--', label='random') #plot
                    ax[row, col].set_xlabel(r"$\mu$m")
                    ax[row, col].set_ylabel(r"$\rho$")
                    ax[row, col].legend()
                    
                    #half distance calculation and display 'real'
                    if len(data) != 0:
                        x = [data_sorted[int(len(data_sorted)/2)]]
                        y = [0.5]
                        ax[row, col].vlines(x, 0, y, linestyle="dashed")
                        ax[row, col].hlines(y, 0, x, linestyle="dashed")
                        ax[row, col].scatter(x, y, c='b', zorder=20)
                        ax[row, col].annotate(round(x[0], 3), (x[0], y[0]))
                        
                        #half distance calculation and display 'random'
                        x = [rd_data_sorted[int(len(rd_data_sorted)/2)]]
                        y = [0.5]
                        ax[row, col].vlines(x, 0, y, linestyle="dashed")
                        ax[row, col].hlines(y, 0, x, linestyle="dashed")
                        ax[row, col].scatter(x, y, c='c', zorder=20)
                        ax[row, col].annotate(round(x[0], 3), (x[0], y[0]))
            
                        ax[row, col].set_xlim(0,None)
                        ax[row, col].set_ylim(0,None)
        
                    if col == 1:
                        ax[row, col].set_title('Distance '+name[0]+'-'+name[1])
                    elif col == 4:
                        ax[row, col].set_title('Distance '+name[1]+'-'+name[0])
                

                    
                if col == 2 or col ==5:#contour density plot
                    if col == 2: 
                        intensity = np.where(dis>=0, t0, -1)
                    elif col == 5:              
                        intensity = np.where(dis>=0, t1, -1)                    
                    
                    intensity = np.ndarray.flatten(intensity)
                    intensity = intensity[intensity != -1]       
                    
                    if len(intensity) > 20: #need a minimum amount of point!
                        #if too numerous data points, take a long time!
                        lenght = 1000 #Random sampling to reduce the calculation
                        if len(data) > lenght:
                            data = np.resize(data, lenght)
                            intensity = np.resize(intensity, lenght)
                        try:
                            X, Y, Z = density_estimation(data, intensity)
                        except:
                            continue
                    
                    ax[row, col].scatter(data, intensity, c='black', marker='.', s=2)
                    ax[row, col].contour(X, Y, Z)
                    ax[row, col].set_xlabel('Distance to '+n1+r" ($\mu$m)")
                    ax[row, col].set_ylabel('Intensity')    
                    
                    if col == 2: #contour density plot title
                        ax[row, col].set_title(n0)

                    elif col == 5: #contour density plot title          
                        ax[row, col].set_title(n1)
             

        fig.set_facecolor('w')
        plt.tight_layout()
        plt.show()

    return dist_data, intensity_data


#visualization
def meanZ_graph(R, G, B, img_name, multi=(1,1,1), label=['Nucleolus','NS2','Nucleus']):
    """
    False 3D visualization. Each channel should be a 3D image binarize.
    A function will remplace all value higher than 0 by the corresponding
    multiplicator value, then realize a mean calculation along the z axis. 
    This mean that the signal intensity is proportional to the thickness.
    By convention, the nucleus is Blue, and to avoid a too strong signal,
    its multiplicator should be smaller that 1.
    Parameters
    ----------
    R : 3D numpy array
        Binarized or labeled image.
    G : 3D numpy array
        Binarized or labeled image.
    B : 3D numpy array
        Binarized or labeled image.
    img_name : str
        Add to the title of the image for identification.
    multi : tuple, optional
        Corresponding multiplicator value for the different channel. The default is (5,2,0.33).
    label : list of str
        list if label for the legend. The default is ['Nucleolus','NS2','Nucleus'].
    
    Returns
    -------
    None.

    """
    #Show 2D results
    fig = plt.figure(figsize=(5,5))
    #Create a RGB image, where all value higher than 0 are considered
    RGB = np.dstack(np.array((np.mean(np.where(R > 0, multi[0], 0), axis=0),
                              np.mean(np.where(G > 0, multi[1], 0), axis=0),
                              np.mean(np.where(B > 0, multi[2], 0), axis=0))))
    # create a patch (proxy artist) for every color
    patches = [mpatches.Patch(color=[1,0,0], label=label[0]),
               mpatches.Patch(color=[0,1,0], label=label[1]),
               mpatches.Patch(color=[0,0,1], label=label[2])]
    # put those patched as legend-handles into the legend
    plt.legend(handles=patches, loc='best', borderaxespad=0., framealpha=1)
    plt.title('Mean z axis of '+img_name)
    plt.axis('off')
    plt.tight_layout()
    plt.imshow(RGB)
    fig.set_facecolor('w')
    plt.show()


def maxZ_graph(R, G, B, img_name, multi=(2,1,0.33), label=['Nucleolus','NS2','Nucleus']):
    """
    False 3D visualization. Each channel should be a 3D image binarize.
    A function will remplace all value higher than 0 by the corresponding
    multiplicator value, then realize a mean calculation along the z axis. 
    This mean that the signal intensity is proportional to the thickness.
    By convention, the nucleus is Blue, and to avoid a too strong signal,
    its multiplicator should be smaller that 1.
    Parameters
    ----------
    R : 3D numpy array
        Binarized or labeled image.
    G : 3D numpy array
        Binarized or labeled image.
    B : 3D numpy array
        Binarized or labeled image.
    img_name : str
        Add to the title of the image for identification.
    multi : tuple, optional
        Corresponding multiplicator value for the different channel. The default is (5,2,0.33).
    label : list of str
        list if label for the legend. The default is ['Nucleolus','NS2','Nucleus'].
    
    Returns
    -------
    None.

    """
    #Show 2D results
    fig = plt.figure(figsize=(5,5))
    #Create a RGB image, where all value higher than 0 are considered
    RGB = np.dstack(np.array((np.max(R, axis=0),
                              np.max(G, axis=0),
                              np.max(B, axis=0))))
    # create a patch (proxy artist) for every color
    patches = [mpatches.Patch(color=[1,0,0], label=label[0]),
               mpatches.Patch(color=[0,1,0], label=label[1]),
               mpatches.Patch(color=[0,0,1], label=label[2])]
    # put those patched as legend-handles into the legend
    plt.legend(handles=patches, loc='best', borderaxespad=0., framealpha=1)
    plt.title('Max z axis of '+img_name)
    plt.axis('off')
    plt.tight_layout()
    plt.imshow(to_8bits(RGB))
    fig.set_facecolor('w')
    plt.show()


def oneC_meanZ(arr, title):
    fig = plt.figure(figsize=(5,5))
    plt.title('Mean z axis of '+title)
    plt.axis('off')
    im = plt.imshow(to_8bits(np.mean(np.where(arr > 0, 1, 0), axis=0)), interpolation='none')
    fig.set_facecolor('w')
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.show()    
    
    
def ortho_view(img, kind, title):
    """
    Generate an ortho view of the 3D img. Can choose between a mean and max view.

    Parameters
    ----------
    img : 3D numpy array
        3 dimensions numpy array (ZXY)
    kind : str
        Kind of projection, only 'mean' and 'max'. mean by default/error
    title : str
        Graph title if any. '' will not add a title

    Returns
    -------
    None.

    """
    if len(img.shape) != 3:
        print('error, bad shape!')
        return
    
    if kind == 'max':
        top = np.max(img, axis=0)
        side = np.rot90(np.max(img, axis=2), k=-1)
        bottom = np.max(img, axis=1)
    else:
        top = np.mean(img, axis=0)
        side = np.rot90(np.mean(img, axis=2), k=-1)
        bottom = np.mean(img, axis=1)
    
    fig, ax = plt.subplots(figsize=(10, 10) , nrows=2, ncols=2,
                           gridspec_kw={'width_ratios': [10, 1],
                                        'height_ratios': [10, 1]})
    
    for a in ax: #remove all axis
        for b in a:
            b.set_axis_off()
    
    ax[0,0].imshow(top)
    ax[0,1].imshow(side)
    ax[1,0].imshow(bottom)
    
    plt.tight_layout() 
    if title != '':
        fig.suptitle(title, fontsize=20)
        fig.subplots_adjust(top=0.95)
    plt.show()
    
    
#segmentation
def thr_list(data, thr):
    """
    Return the threshold value of the threshold method applied on data
    Parameters
    ----------
    data : array
        array, whatever is the dimension (flatten).
    thr : string
        Threshold name. minimum, yen, otsu, mean, triangle, isodata and li for the moment.

    Returns
    -------
    thr_value : int or float
        Value of the thresholded data.

    """
    data = np.ndarray.flatten(data)
    if thr == 'minimum':
        thr_min = filters.threshold_minimum(data)
    if thr == 'yen':
        thr_min = filters.threshold_yen(data)
    if thr == 'otsu':
        thr_min = filters.threshold_otsu(data)
    if thr == 'mean':
        thr_min = filters.threshold_mean(data)
    if thr == 'triangle':
        thr_min = filters.threshold_triangle(data)
    if thr == 'isodata':
        thr_min = filters.threshold_isodata(data)
    if thr == 'li':
        thr_min = filters.threshold_li(data)
    
    return thr_min



def nucleus_detection(img, log, thr_method='otsu', bin_iter=1, nucleus_size=200000):
    """
    ?Misleading name?
    Process an image to segment it using one of the three threshold methods.
    In addition, an opening step is realized to remove any single pixel.
    It finish with a filter that keep the biggest structure (the first one).
    Parameters
    ----------
    img : 3D numpy array
        3D image to segment.
    log : str
        string containing the log information
    thr_method : str, optional
        Different method for thresholding. otsu, yen or minimum. The default is 'otsu'.
    bin_iter : int, optional
        Number of iteration of opening, in order to remove unwanted signal. The default is 1.
    nucleus_size: int, optionnal
        Minimum volume of the nucleus in pixel. Do not loop through all structure

    Returns
    -------
    binary : 3D numpy array
        Resulting segmented image.
    log : str
        Log data.

    """
    #calculate the global threshold accross the z stack while ignoring the 0
    cleaned_data = boolean_remove_zeros(np.ndarray.flatten(np.mean(img, axis=0))) #w/o 0
    if thr_method == 'yen':
        thr_min = filters.threshold_yen(cleaned_data)
    elif thr_method == 'minimum':
        thr_min = filters.threshold_minimum(cleaned_data)
    else:
        thr_method = 'otsu'
        thr_min = filters.threshold_otsu(cleaned_data)
    
    log += '- - - - - - - - - - - - - - - -\n'
    log += 'nucleus setting:\nThreshold method: '+thr_method+'\nbin_iter: '+str(bin_iter)
    log += '\nnucleus size: '+str(nucleus_size)+'\nCalculated threshold value: '+str(thr_min)+'\n'
    log += '- - - - - - - - - - - - - - - -\n'
    
    #binarise and simple operation
    binary = img > thr_min
    binary = morphology.binary_fill_holes(binary) #fill holes, if any
    binary = morphology.binary_opening(binary, iterations=bin_iter) #remove isolated pixel
    binary = morphology.binary_closing(binary, iterations=bin_iter)
    binary = morphology.binary_fill_holes(binary) #fill holes, if any
    binary = morphology.binary_erosion(binary, iterations=bin_iter)

    #filter by size, keeping the largest structure
    markers = measure.label(binary, background = 0)
    nucleus = np.zeros(img.shape, dtype=bool)
    
    if len(nucleus.shape) == 3: min_size = nucleus_size
    if len(nucleus.shape) == 2: min_size = nucleus_size/10
    
    for x in range(1, np.max(markers)+1):
        if len(markers[markers == x]) > min_size:
            nucleus[markers == x] = True
            break #work if there is only one nucleus

    return nucleus, log


def signal_detection(img, nucl_binary, log, thr_method='otsu', n_iter=1, exclude_zeros=True, focus_nucleus=True):

    if focus_nucleus: #select only nucleus pixel directly on the 3D image
        data = np.where(nucl_binary == 1, img, 0)
    else:
        # data = img[1]
        data = img
        
    if exclude_zeros: #remove all zeros pixel (area outside the nucleus)
        data_thr = boolean_remove_zeros(np.ndarray.flatten(data)) #w/o 0
    else:
        data_thr = np.ndarray.flatten(data)
    
    thr_min = thr_list(data_thr, thr_method)

    #simple threshold with opening/closing to remove single pixel and single hole
    binary = data > thr_min
    binary = sc.ndimage.binary_opening(binary, iterations=n_iter)
    binary = sc.ndimage.binary_closing(binary, iterations=n_iter)
    
    log += 'Threshold settings:\nThreshold method: '+thr_method+'\nn_iter: '+str(n_iter)
    log += '\nexclude zeros: '+str(exclude_zeros)+'\nfocus on nucleus: '+str(focus_nucleus)
    log += '\nCalculated threshold value: '+str(thr_min)+'\n'
    log += '- - - - - - - - - - - - - - - -\n'
    
    return binary, log


def measure_stats(binary, img):
    """
    Get some parameters from nucleolin measurment, using the regionprops from scikit

    Parameters
    ----------
    nucleolin_binary : 3D numpy array
        The segmented and labeled nucleolin 3D image.
    img : 3D numpy array
        The channel to measure intensity from
    img_name: str
        Path of the image from which the analysis is done

    Returns
    -------
    props : Dataframe
        Pandas dataframe, use the concatenate function to merge with other.
        Add the path of the file as 'path'.

    """
    binary = measure.label(binary, background = 0)
    
    to_measure = ['label', 'area', 'centroid', 'convex_area', 'euler_number',
                  'inertia_tensor', 'major_axis_length', 'minor_axis_length',
                  'max_intensity', 'mean_intensity','min_intensity', 'solidity',
                  'moments_central', 'bbox']
    # to_measure = ['label', 'area', 'centroid', 'major_axis_length', 'minor_axis_length',
    #               'max_intensity', 'mean_intensity','min_intensity','bbox']
    props = pd.DataFrame.from_dict(measure.regionprops_table(binary, img,
                                                            properties=to_measure))

    
    return props


def vertical_best_count(matrices, n=5, graph=False):
    """
    Mask a 3D array of bool with the continuous value in Z higher than n.
    original: https://stackoverflow.com/questions/44439703/python-find-consecutive-values-in-3d-numpy-array-without-using-groupby
    Parameters
    ----------
    matrices : 3D array of bool
        numpy array to get the best count
    n : int, optional
        Number of coninuous value. The default is 5.
    graph : bool, optional
        To display the result in false 3d. The default is False.

    Returns
    -------
    results : 3D array of bool
        Same shape than original, but with value not corresping removed

    """
    bests = np.zeros(matrices.shape[1:])
    counter = np.zeros(matrices.shape[1:])
    
    for depth in range(matrices.shape[0]):
        this_level = matrices[depth, :, :]
        counter = counter + this_level
        bests = (np.stack([bests, counter], axis=0)).max(axis=0)
    
    results = np.zeros(matrices.shape)
    for depth in range(matrices.shape[0]):
        this_level = matrices[depth, :, :]
        results[depth] = np.logical_and(bests > n, this_level)
    
    if graph:
        plt.imshow(np.sum(results, axis=0), interpolation='none')
        plt.title('results with n='+str(n))
        plt.axis('off')
        plt.colorbar(fraction=0.046, pad=0.04)
        plt.show()
    
    return results


def size_filter(img, limit=300):
    """
    remove structure smaller than the size limit

    Parameters
    ----------
    img : 2D or 3D array
        The image ( as numpy bool array) to remove from
    limit : int
        The minimum number of pixel composing a structure. The default is 300.

    Returns
    -------
    Numpy bool array
        Same shape than original img, without the smaller structure

    """
    img = measure.label(img)
    for x in range(1, np.max(img)+1):
        t = np.where(img == x)
        if len(t[0]) < limit:
            img[t] = 0
        
    return np.array(np.where(img > 0, 1, 0), dtype=bool)


#utils
def sorted_nicely(l): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


def boolean_remove_zeros(arr):
    """Mini function to remove all zeros from an array"""
    return arr[arr != 0]


#data management
def to_8bits(array):
    """Normalize an array between 0 and 1 before transforming to a 8 bit"""
    array = (array - np.min(array))/(np.max(array)-np.min(array))
    return np.array(array*255, dtype='uint8')


def Zto_RGB(array):
    """Simple function to make a max projection of a 3D multichannel to a 2D
    matplotlib readable format."""
    RGB = np.dstack(np.array((np.max(array[2], axis=2),
                              np.max(array[1], axis=2),
                              np.max(array[0], axis=2))))
    
    return np.array(RGB, dtype='uint8')


def to_RGB(array):
    """Simple function to transform an 2D multichannel image in a 2D matplotlib
    readable format"""
    if array.shape[0] == 3:
        RGB = np.dstack((array[2], array[1], array[0]))
    elif array.shape[0] == 1:
        RGB = np.dstack((array[0], array[0], array[0]))
    else:
        print(array.shape)
    return np.array(RGB, dtype='uint8')


def img_color_merger(r=[], g=[], b=[], gr=[], c=[], m=[], y=[]):
    """
    Merge color channel together to obtain a RGB image for plotting (x, y, c).
    Merging is done by maximum choices.
    Tested only with 2D array.
    Parameters
    ----------
    r : 2D array, optional
        DESCRIPTION. The default is [].
    g : 2D array, optional
        DESCRIPTION. The default is [].
    b : 2D array, optional
        DESCRIPTION. The default is [].
    gr : 2D array, optional
        DESCRIPTION. The default is [].
    c : 2D array, optional
        DESCRIPTION. The default is [].
    m : 2D array, optional
        DESCRIPTION. The default is [].
    y : 2D array, optional
        DESCRIPTION. The default is [].

    Returns
    -------
    3D array int8
        RGB merged color array, in the x,y,3 shape

    """
    a = [r,g,b,gr,c,m,y]
    for color in a: #get the shape of the image
        if len(color) != 0:
            blend = np.zeros((3,)+color.shape)
    #classical RGB
    if len(r) != 0:
        blend[0] = r
    if len(g) != 0:
        blend[1] = g
    if len(b) != 0:
        blend[2] = b
    #gray
    if len(gr) != 0:
        blend = np.maximum.reduce([blend, np.array([gr, gr, gr])])
    #less classical CMY color
    if len(c) != 0:
        c = np.array([np.zeros(c.shape), c, c])
        blend = np.maximum.reduce([blend, c])
    if len(m) != 0:
        m = np.array([m, np.zeros(m.shape), m])
        blend = np.maximum.reduce([blend, m])
    if len(y) != 0:
        y = np.array([y, y, np.zeros(y.shape)])
        blend = np.maximum.reduce([blend, y])
    #transform in matplotlib readable format
    blend = np.dstack((blend[0], blend[1], blend[2]))
    return np.array(blend, dtype='uint8')


# Saving/loading function.
class Data_Distance:
    def __init__(self, name, nucleus_r, nucleolus_r, dist_data, intensity_data, seg, log):
        self.name = name
        self.nucleus_r = nucleus_r
        self.nucleolus_r = nucleolus_r
        self.dist_data = dist_data
        self.intensity_data = intensity_data
        self.seg = seg
        self.log = log
        
        
def save_data(obj, dim):
    obj.log += 'Save data as: '+obj.name+'_'+str(dim)+'D_results.pckl'
    print(obj.log)
    with open(obj.name+'_'+str(dim)+'D_results.pckl', 'wb') as f:
        pickle.dump(obj, f)


def load_data(path_name):
    with open(path_name, 'rb') as f:
        obj = pickle.load(f)
    return obj

#plot function
def get_color(nb, kind, dic=False):
    """
    Get a list of colors (RGB 0-1) of lenght nb. Two style, bars and points.
    Actually have 8 differents colors, but cycle through them to generate 
    enough points. Based on https://healthdataviz.com/2012/02/02/optimal-colors-for-graphs/

    Parameters
    ----------
    nb : int
        The number of color required
    kind : str
        Choose between bars and points (if empty/invalid, return points)
    dic : bool, optional
        If True, return the base dictionnary with the color value (len(8)). The default is False.

    Returns
    -------
    colors : list
        List of 3-tuples (RGB) of lenght nb.

    """
    bars = {'blue':(114/255, 147/255, 203/255),
           'orange':(225/255, 151/255, 76/255),
           'green':(132/255, 186/255, 91/255),
           'red':(211/255, 94/255, 96/255),
           'grey':(128/255, 133/255, 133/255),
           'purple':(144/255, 103/255, 167/255),
           'bordeau':(171/255, 104/255, 87/255),
           'gold':(204/255, 194/255, 16/255)
        }

    points = {'blue':(57/255, 106/255, 177/255),
           'orange':(218/255, 124/255, 48/255),
           'green':(62/255, 150/255, 81/255),
           'red':(204/255, 37/255, 41/255),
           'grey':(83/255, 81/255, 84/255),
           'purple':(107/255, 76/255, 154/255),
           'bordeau':(146/255, 36/255, 40/255),
           'gold':(148/255, 139/255, 61/255)
        }
    
    if kind == 'bars':
        colors = bars
    else:
        colors = points
    
    if dic: return colors
    
    nb = int(nb) #security
    if nb <= len(colors):
        colors = list(colors.values())[:nb]
    else:
        colors = [list(colors.values()) for x in range(int(np.ceil((nb/len(colors)))))]
        colors = [item for sublist in colors for item in sublist]
    
    return colors




