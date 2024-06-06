# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 13:08:53 2023

@author: sleclerc
"""

import numpy as np
import os
from statannotations.Annotator import Annotator
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import utils as u
from scipy import stats
import matplotlib.transforms as transforms
import matplotlib as mpl
mpl.rcParams['mathtext.default'] = 'regular'

def make_swarmplot(data, ax, dotsize):
    sns.boxplot(data=data, ax=ax, 
                     showcaps=False, boxprops={'facecolor':'None'},
                     showfliers=False, whiskerprops={'linewidth':0},
                     medianprops={'linewidth':0}, showmeans=True,
                     meanline=True, meanprops={'color':'black'})
    sns.swarmplot(data=data, size=dotsize, ax=ax, zorder=0)
    annotator = Annotator(ax, pairs, data=data)
    annotator.configure(test=test, verbose=True).apply_and_annotate()

def get_data(data, zscore=3):
    if zscore == 0:
        data = {x:data[x] for x in data}
    else:
        data = {x:np.array(data[x])[np.abs(stats.zscore(data[x])) < zscore] for x in data}
    data = pd.DataFrame(dict([(k, pd.Series(v)) for k,v in data.items()]))
    return data


dotsize = 1
test = 't-test_ind'

color = u.get_color(3, 'original')
sns.set_palette(sns.color_palette(color, len(color)))


path_ni = "NI"+os.sep+"NI_dist_map.pckl.npy"
path_8hpi = "8hpi"+os.sep+"8hpi_dist_map.pckl.npy"
path_12hpi = "12hpi"+os.sep+"12hpi_dist_map.pckl.npy"

data_ni = np.load(path_ni)
data_8hpi = np.load(path_8hpi)
data_12hpi = np.load(path_12hpi)

# #random selection dataset
# data_ni = np.random.choice(data_ni, size=len(data_ni)//50, replace=False)
# # data_8hpi = np.random.choice(data_8hpi, size=len(data_8hpi)//50, replace=False)
# data_12hpi = np.random.choice(data_12hpi, size=len(data_12hpi)//50, replace=False)

# data_ni_contact_site = data_ni[data_ni <= 30]
# #data_8hpi_contact_site = data_8hpi[data_8hpi <= 30]
# data_12hpi_contact_site = data_12hpi[data_12hpi <= 30]


condition_name_list = ['Noninfected',
                       '8 hpi',
                       '12 hpi']
pairs = [['Noninfected', '8 hpi'],
          ['Noninfected', '12 hpi'],
          ['8 hpi', '12 hpi']]


# key_measured = [data_ni_contact_site,
#                 #data_8hpi_contact_site,
#                 data_12hpi_contact_site]

#remove outliers
data_ni = np.array(data_ni)[np.abs(stats.zscore(data_ni)) < 3]
data_8hpi = np.array(data_8hpi)[np.abs(stats.zscore(data_8hpi)) < 3]
data_12hpi = np.array(data_12hpi)[np.abs(stats.zscore(data_12hpi)) < 3]

#full dataset
key_measured = [data_ni,
                data_8hpi,
                data_12hpi]

key_measured = {condition_name_list[idx]:x for idx, x in enumerate(key_measured)}
key_measured = pd.DataFrame(dict([(k, pd.Series(v)) for k,v in key_measured.items()]))
key_measured = key_measured[['Noninfected',
                             '8 hpi',
                             '12 hpi']]
#%%
test = 'Brunner-Munzel'
div = 2
fig = plt.figure(figsize=(15/div,5/div))
ax = fig.add_subplot(131)
key_measured.to_csv('Distance ER-Mito.csv')
sns.violinplot(data=key_measured, cut=0, ax=ax, scale='count', inner=None, zorder=0)
graph = sns.boxplot(data=key_measured, ax=ax, 
                    showcaps=False, boxprops={'facecolor':'None'},
                    showfliers=False, whiskerprops={'linewidth':0},
                    medianprops={'linewidth':0}, showmeans=True,
                    meanline=True, meanprops={'color':'black'}, zorder=500)
plt.setp(ax.lines, zorder=2000)

#add the 30 nm line for contact site
graph.axhline(30, c='r', linestyle='-.')
trans = transforms.blended_transform_factory(ax.get_yticklabels()[0].get_transform(), ax.transData)
ax.text(0,30, "{:.0f}".format(30), color="red", transform=trans, 
        ha="right", va="center")

annotator = Annotator(ax, pairs, data=key_measured)
annotator.configure(test=test, verbose=True).apply_and_annotate()

ax.set_title('ER-mitochondria distance')
ax.set_ylabel('Distance (nm)')


ax = fig.add_subplot(132)

n_ni = np.load("NI"+os.sep+'ni_contact_site_n.npy')
n_8hpi = np.load("8hpi"+os.sep+"8hpi_contact_site_n.npy")
n_12hpi = np.load("12hpi"+os.sep+"12hpi_contact_site_n.npy")

m_ni = np.load("NI"+os.sep+'SA_mito.npy')
m_8hpi = np.load("8hpi"+os.sep+"8hpi_SA_mito.npy")
m_12hpi = np.load("12hpi"+os.sep+"12hpi_SA_mito.npy")

data = {'Noninfected':n_ni/(m_ni/1000000),
        '8 hpi':n_8hpi/(m_8hpi/1000000),
        '12 hpi':n_12hpi/(m_12hpi/1000000)}

data = get_data(data, zscore=0)
data.to_csv('contact site per mito SA.csv')
make_swarmplot(data, ax, 4) #to get the unit in um2 (from nm2)

ax.set_title('Number of contact site/mitochondria')
ax.set_ylabel('Nb contact site/mitochondria ${\mu m^2}$')
#ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

ax = fig.add_subplot(133)

n_ni = np.load("NI"+os.sep+'ni_area_contact_site.npy')
n_8hpi = np.load("8hpi"+os.sep+"8hpi_area_contact_site.npy")
n_12hpi = np.load("12hpi"+os.sep+"12hpi_area_contact_site.npy")

data = {'Noninfected':n_ni,
        '8 hpi':n_8hpi,
        '12 hpi':n_12hpi}

data = get_data(data, zscore=1.5) 
data /= 1000000 #to get the unit in um2 (from nm2)
data.to_csv('contact site area.csv')
make_swarmplot(data, ax, 1.5)

ax.set_title('Area of contact site')
ax.set_ylabel('Contact site area (${\mu m^2}$)')
#ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.tight_layout()
plt.savefig('ER-mitochondria distance.png', dpi=300)
plt.savefig('ER-mitochondria distance.svg')
plt.show()
