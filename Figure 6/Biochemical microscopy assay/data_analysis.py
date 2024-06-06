# -*- coding: utf-8 -*-
"""
Created on Mon May 22 09:15:27 2023

@author: sleclerc


This script performs exploratory data analysis and generates graphs for different conditions.

The code performs the following tasks:
    Get a list of files in the current directory with the extension '.pckl'.
    Define variables for dot size, statistical test, and pairs of conditions for comparison.
    Define two functions: calculate_data and grab_and_show.
    Define a function exploratory_figure that generates exploratory figures.
    Group the data based on file names and perform calculations and corrections on each group.
    Generate exploratory figures for each group.
    Generate a final figure with three subplots.
    Save the figures in the 'Results' directory.
    
The code analyze and visualize data related to mitochondria, including nucleus area, MitoTracker signal, infection signal, and mitochondria signal intensity. It applies statistical tests and generates boxplots, swarmplots, and annotations to compare different conditions. The resulting figures are saved for further analysis and consultation.

"""

import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from statannotations.Annotator import Annotator
import matplotlib as mpl
import numpy as np
from scipy import stats
mpl.rcParams['mathtext.default'] = 'regular'

files = os.listdir(os.getcwd())
files = [x for x in files if x.endswith('.pckl')]

dotsize = 3
test = 't-test_ind' #statistical test to do
#condition to test
pairs = [['Noninfected', '4 hpi'],
         ['Noninfected', '8 hpi'],
         ['4 hpi', '8 hpi']]

#functions
def calculate_data(dic):
    """
    Calculae and correct some value for the DataFrames

    Parameters
    ----------
    dic : dict
        Dictionnary of DataFrame, one by condition

    Returns
    -------
    dic : dict
        Results of the cacluation

    """
    for key in dic: #correction by nucleus intensity
        df = dic[key]
        df['nucleus_total'] = df.area_nucleus * df.mean_int_nucleus
        df['nuc_norm'] = df.nucleus_total.mean() / df.nucleus_total
        df['norm_int_nucleus'] = df.mean_int_nucleus# * df.nuc_norm
        df['norm_int_cytoplasm'] = df.mean_int_cytoplasm# * df.nuc_norm
        df['norm_int_infection'] = df.mean_int_infection# * df.nuc_norm
        df['norm_int_mito'] = df.mean_int_mito_signal# * df.nuc_norm


        #Corrected value        
        df['cyto_total'] = df['norm_int_cytoplasm'] * df.area_cytoplasm
        df['inf_total'] = df['norm_int_infection'] * (df.area_infection+1)
        df['mito_total'] = df['norm_int_mito'] * df.area_cytoplasm
        
        df['mito_total_adj'] = df['mito_total'] / df['cyto_total']
        
        df = df[(np.abs(stats.zscore(df.nucleus_total)) < 3)]
        df = df[(np.abs(stats.zscore(df.cyto_total)) < 3)]
        
        dic[key] = df  # Assign the modified DataFrame back to the dictionary
        
    return dic

def grab_and_show(table, condition, title, ax):
    """
    extract the condition from the table and create a graph in ax with the given title
    Parameters
    ----------
    table : dict
        Dictionnary of DataFrame
    condition : str
        Name of the column to extract
    title : str
        Title to give to the graph
    ax : plt.axis
        Valid matplolib axis object

    Returns
    -------
    None.

    """
    result = [[],[],[]]
    for key in table:
        if key == 'NI': idx=0
        elif key == '4hpi': idx=1
        elif key == '8hpi': idx=2
        df = table[key]
        result[idx] += list(df[condition])
    
    condition_name_list = ['Noninfected', '4 hpi', '8 hpi']
    result = {condition_name_list[idx]:x for idx, x in enumerate(result)}
    result = pd.DataFrame(dict([(k, pd.Series(v)) for k,v in result.items()]))
    result = result[['Noninfected', '4 hpi','8 hpi']]
    
    ax.set_title(title)
    sns.boxplot(data=result, ax=ax, 
                     showcaps=False, boxprops={'facecolor':'None'},
                     showfliers=False, whiskerprops={'linewidth':0},
                     medianprops={'linewidth':0}, showmeans=True,
                     meanline=True, meanprops={'color':'black'})
    plt.setp(ax.lines, zorder=2000)
    sns.swarmplot(data=result, size=dotsize, ax=ax)
    #sns.violinplot(data=mito_l_mean, ax=ax, scale='count', inner=None, cut=0)
    annotator = Annotator(ax, pairs, data=result)
    annotator.configure(test=test, verbose=True).apply_and_annotate()
    ax.set_ylabel('Corrected intensity')


def exploratory_figure(dic, cond):
    """
    Generate 4 swarmplots to show the nucleus area, mitotracker signal, 
    infection signal and mitochondria signal. Save the resulting graph
    for latter consultation.

    Parameters
    ----------
    dic : dict
        of DataFrame, one entry by condition
    cond : str
        Name of the condition

    Returns
    -------
    None.

    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(2,2,1)
    grab_and_show(dic, 'area_nucleus', cond+', Nucleus Area', ax)
    ax = fig.add_subplot(2,2,2)
    grab_and_show(dic, 'cyto_total', cond+', Corrected MitoTracker signal', ax)
    ax = fig.add_subplot(2,2,3)
    grab_and_show(dic, 'inf_total', cond+', Corrected infection signal', ax)
    ax = fig.add_subplot(2,2,4)
    grab_and_show(dic, 'mito_total_adj', cond+', Corrected Mitochondria signal intensity', ax)
    
    plt.tight_layout()
    plt.savefig('Results'+os.sep+cond+' analysis.png', dpi=300)
    plt.show()
 


#group the data
m_only = {}
m_t = {}
m_s = {}
m_r = {}
for f in files:
    df = pd.read_pickle(f)
    if 'NI' in f: name = 'NI'
    elif '4hpi' in f: name = '4hpi'
    elif '8hpi' in f: name = '8hpi'
    if f.endswith('M.pckl'): m_only[name] = df
    elif f.endswith('M+T.pckl'): m_t[name] = df
    elif f.endswith('M+S.pckl'): m_s[name] = df
    elif f.endswith('M+R.pckl'): m_r[name] = df
        
# Add column to all data
m_only = calculate_data(m_only)
m_t = calculate_data(m_t)
m_s = calculate_data(m_s)
m_r = calculate_data(m_r)


# check all data
exploratory_figure(m_only, 'Mitotracker only')
exploratory_figure(m_t, 'TMRM')
exploratory_figure(m_s, 'MitoSox')
exploratory_figure(m_r, 'Rhod-2 (Ca)')

#make figure
fig = plt.figure(figsize=(9, 3))
dotsize = 1
grab_and_show(m_t, 'mito_total_adj', 'Membrane potential - TMRM', fig.add_subplot(1,3,1))
grab_and_show(m_s, 'mito_total_adj', 'SuperOxide - MitoSox', fig.add_subplot(1,3,2))
grab_and_show(m_r, 'mito_total_adj', 'Calcium - Rhod2', fig.add_subplot(1,3,3))
plt.tight_layout()
plt.savefig('Results'+os.sep+'Final analysis.png', dpi=300)
plt.savefig('Results'+os.sep+'Final analysis.svg')
plt.show()

#%%



