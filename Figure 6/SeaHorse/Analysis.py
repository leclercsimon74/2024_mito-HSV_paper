# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 10:41:13 2023

@author: sleclerc
"""


import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
mpl.rcParams['mathtext.default'] = 'regular'
from statannotations.Annotator import Annotator

import os
import pandas as pd
import seaborn as sns
from scipy import stats
import numpy as np

import utils as u

color = u.get_color(4, 'original')
sns.set_palette(sns.color_palette(color, len(color)))


def make_swarmplot(data, ax, dotsize):
    sns.boxplot(data=data, ax=ax, 
                     showcaps=False, boxprops={'facecolor':'None'},
                     showfliers=False, whiskerprops={'linewidth':0},
                     medianprops={'linewidth':0}, showmeans=True,
                     meanline=True, meanprops={'color':'black'})
    sns.swarmplot(data=data, size=dotsize, ax=ax, zorder=0)
    annotator = Annotator(ax, pairs, data=data)
    annotator.configure(test=test, verbose=True).apply_and_annotate()

def get_data(cat, zscore=3):
    data = {x:np.array(df[df['Group Name'] == x][cat]) for x in condition_name_list}
    data = {x:np.array(data[x])[np.abs(stats.zscore(data[x])) < zscore] for x in data}
    data = pd.DataFrame(dict([(k, pd.Series(v)) for k,v in data.items()]))
    return data

df = pd.read_csv('Results.csv')


test = 't-test_ind' 
dotsize = 10
zscore = 1.5
with_title = False
condition_name_list = ['Noninfected', '4 hpi', '8 hpi']
pairs = [['Noninfected', '4 hpi'],
         ['Noninfected', '8 hpi'],
         ['4 hpi', '8 hpi']]




div = 3
fig = plt.figure(figsize=(20/div, 20/div))



#2
ax = fig.add_subplot(2,2,1)
cat = 'Basal Respiration'
data = get_data(cat, zscore=zscore)
make_swarmplot(data, ax, dotsize)
if with_title: ax.set_title(cat)
ax.set_ylabel('Basal respiration OCR (pmol/min)')

#4
ax = fig.add_subplot(2,2,2)
cat = 'Proton Leak'
data = get_data(cat, zscore=zscore)
make_swarmplot(data, ax, dotsize)
if with_title: ax.set_title(cat)
ax.set_ylabel('Proton leak OCR (pmol/min)')
#5
ax = fig.add_subplot(2,2,3)
cat = 'ATP Production'
data = get_data(cat, zscore=zscore)
data.to_csv('csv'+os.sep+'ATP production.csv')
make_swarmplot(data, ax, dotsize)
if with_title: ax.set_title(cat)
ax.set_ylabel('ATP production OCR (pmol/min)')

#8
ax = fig.add_subplot(2,2,4)
cat = 'Coupling Efficiency'
data = get_data(cat, zscore=zscore)
data = data * 100
data.to_csv('csv'+os.sep+'Coupling efficiency.csv')
make_swarmplot(data, ax, dotsize)
if with_title: ax.set_title(cat)
ax.set_ylabel('% Coupling efficiency')



plt.savefig('Main graphes.png', dpi=600)
plt.savefig('Main graphes.svg')
plt.tight_layout()
plt.show()






fig = plt.figure(figsize=(20/div, 20/div))
#1
ax = fig.add_subplot(2,2,1)
cat = 'Non-Mitochondrial Oxygen Consumption'
data = get_data(cat, zscore=zscore)
make_swarmplot(data, ax, dotsize)
if with_title: ax.set_title(cat)
ax.set_ylabel('OCR (pmol/min)')

#3
ax = fig.add_subplot(2,2,2)
cat = 'Maximal Respiration'
data = get_data(cat, zscore=zscore)
make_swarmplot(data, ax, dotsize)
if with_title: ax.set_title(cat)
ax.set_ylabel('OCR (pmol/min)')

#6
ax = fig.add_subplot(2,2,3)
cat = 'Spare Respiratory Capacity'
data = get_data(cat, zscore=zscore)
make_swarmplot(data, ax, dotsize)
if with_title: ax.set_title(cat)
ax.set_ylabel('OCR (pmol/min)')
#7
ax = fig.add_subplot(2,2,4)
cat = 'Spare Respiratory Capacity as a %'
data = get_data(cat, zscore=zscore)
make_swarmplot(data*100, ax, dotsize)
if with_title: ax.set_title(cat)
ax.set_ylabel('%')

plt.tight_layout()
plt.savefig('Supplementary graphes.png', dpi=600)
plt.savefig('Supplementary graphes.svg')
plt.show()


#%% Main figure!

df2 = pd.read_csv('Rate_data.csv')


X = np.unique(df2['Time'])


ni_data = df2[df2['Group'] == 'Noninfected']
ni_data = ni_data.drop(ni_data[ni_data['Well']=='D01'].index)

hpi_4_data = df2[df2['Group'] == '4 hpi']
hpi_4_data = hpi_4_data.drop(hpi_4_data[hpi_4_data['Well']=='B03'].index)

hpi_8_data = df2[df2['Group'] == '8 hpi']
hpi_8_data = hpi_8_data.drop(hpi_8_data[hpi_8_data['Well']=='A03'].index)

dotsize = 8
capsize = 3
div = 3.5
fig = plt.figure(figsize=(20/div, 40/div))

#Global graph
ax = fig.add_subplot(2,1,1)

ax.errorbar(X, 
            ni_data.groupby('Measurement').mean()['OCR'], 
            ni_data.groupby('Measurement').std()['OCR']/len(ni_data.groupby('Measurement'))**0.5,
            label='Noninfected',
            capsize=capsize)

ax.errorbar(X, 
            hpi_4_data.groupby('Measurement').mean()['OCR'], 
            hpi_4_data.groupby('Measurement').std()['OCR']/len(hpi_4_data.groupby('Measurement'))**0.5,
            label='4 hpi',
            capsize=capsize)

ax.errorbar(X, 
            hpi_8_data.groupby('Measurement').mean()['OCR'], 
            hpi_8_data.groupby('Measurement').std()['OCR']/len(hpi_8_data.groupby('Measurement'))**0.5,
            label='8 hpi',
            capsize=capsize)

ax.set_ylabel('OCR (pmol/min)')
ax.set_xlabel('Time (min)')
#line to indicate drug injection timepoint
ax.axvline(35.58, linestyle='--', color='k', zorder=0)
ax.axvline(61.36, linestyle='--', color='k', zorder=0)
ax.axvline(87.13, linestyle='--', color='k', zorder=0)

ax.legend(frameon=False)

with_title = False
#2
ax = fig.add_subplot(4,2,5)
cat = 'Basal Respiration'
data = get_data(cat, zscore=zscore)
data.to_csv('csv'+os.sep+'Basal respiration.csv')
make_swarmplot(data, ax, dotsize)
if with_title: ax.set_title(cat)
ax.set_ylabel('Basal respiration OCR (pmol/min)')

#4
ax = fig.add_subplot(4,2,6)
cat = 'Proton Leak'
data = get_data(cat, zscore=zscore)
data.to_csv('csv'+os.sep+'Proton leak.csv')
make_swarmplot(data, ax, dotsize)
if with_title: ax.set_title(cat)
ax.set_ylabel('Proton leak OCR (pmol/min)')
#5
ax = fig.add_subplot(4,2,7)
cat = 'ATP Production'
data = get_data(cat, zscore=zscore)
data.to_csv('csv'+os.sep+'ATP production.csv')
make_swarmplot(data, ax, dotsize)
if with_title: ax.set_title(cat)
ax.set_ylabel('ATP production OCR (pmol/min)')
#8
ax = fig.add_subplot(4,2,8)
cat = 'Coupling Efficiency'
data = get_data(cat, zscore=zscore)
data *= 100
data.to_csv('csv'+os.sep+'Coupling efficiency.csv')
make_swarmplot(data, ax, dotsize)
if with_title: ax.set_title(cat)
ax.set_ylabel('% Coupling efficiency')


plt.tight_layout()
plt.savefig('Main graphes2.png', dpi=600)
plt.savefig('Main graphes2.svg')
plt.show()





