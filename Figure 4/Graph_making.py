# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 14:55:47 2023

@author: sleclerc
"""

import utils as u

import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
import pandas as pd
from statannotations.Annotator import Annotator
import matplotlib.transforms as transforms
import matplotlib as mpl
mpl.rcParams['mathtext.default'] = 'regular'

color = u.get_color(4, 'original')
sns.set_palette(sns.color_palette(color, len(color)))

def outlier_removal(data, zscore):
    #data is a dict
    for key in data:
        data[key] = np.array(data[key])[np.abs(stats.zscore(data[key])) < zscore]
    return data


div = 2
fig = plt.figure(figsize=(6/div, 20/div))
#violin plot of the whole dataset of foci ExM


df = pd.read_pickle('ExpM'+os.sep+'Foci_distance_analysis.pckl')
size = 5
dist_thr = 0.060 #in um
df_filter = df[df.area>size]
pairs = [['NI', '8 hpi']]
test = 'Brunner-Munzel'
f = lambda b: np.log(b + 1)
g = lambda a: np.exp(a) - 1

data = {'NI':df_filter[df_filter.Condition=='NI']['Distance'],
        '8 hpi':df_filter[df_filter.Condition=='8 hpi']['Distance']}
data = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in data.items()]))
data.to_csv('csv'+os.sep+'VAPB distance.csv')
ax = fig.add_subplot(411)
sns.violinplot(data=data, cut=0, ax=ax, scale='count', inner=None, zorder=0)
graph = sns.boxplot(data=data, ax=ax, 
                    showcaps=False, boxprops={'facecolor':'None'},
                    showfliers=False, whiskerprops={'linewidth':0},
                    medianprops={'linewidth':0}, showmeans=True,
                    meanline=True, meanprops={'color':'black'}, zorder=500)
plt.setp(ax.lines, zorder=2000)
#add the line for contact site
graph.axhline(dist_thr, c='r', linestyle='-.')
trans = transforms.blended_transform_factory(ax.get_yticklabels()[0].get_transform(), ax.transData)
ax.text(0,dist_thr, "{:.3f}".format(dist_thr), color="red", transform=trans, 
        ha="right", va="center")
annotator = Annotator(ax, pairs, data=data)
annotator.configure(test=test, verbose=True).apply_and_annotate()
#ax.set_title('VAPB foci to mitochondria distance')
ax.set_ylabel('VAPB distance (um)')
ax.set_yscale(mpl.scale.FuncScale(ax, (f, g)))



#PLA
size = 0.02 #size in um3 of the PLA foci
zscore = 3
df = pd.read_pickle('PLA'+os.sep+'PLA_foci_distance.pckl')
df_filter = df[df.area>size]
cell_dotsize = 3
test = 't-test_ind'

pairs = [['NI', '4 hpi'],
         ['NI', '8 hpi'],
         ['NI', '12 hpi'],
         ['4 hpi', '8 hpi'],
         ['4 hpi', '12 hpi'],
         ['8 hpi', '12 hpi']]

def make_swarmplot(data, ax, dotsize, pairs):
    sns.boxplot(data=data, ax=ax, 
                      showcaps=False, boxprops={'facecolor':'None'},
                      showfliers=False, whiskerprops={'linewidth':0},
                      medianprops={'linewidth':0}, showmeans=True,
                      meanline=True, meanprops={'color':'black'})
    sns.swarmplot(data=data, size=dotsize, ax=ax, zorder=0)
    annotator = Annotator(ax, pairs, data=data)
    annotator.configure(test=test, verbose=True).apply_and_annotate()

# pairs = [['8 hpi', '12 hpi']]
ax = fig.add_subplot(412)
data = {
        'NI':df_filter[df_filter.Condition=='NI'].groupby('Name').count()['Distance'],
        '4 hpi':df_filter[df_filter.Condition=='4 hpi'].groupby('Name').count()['Distance'],
        '8 hpi':df_filter[df_filter.Condition=='8 hpi'].groupby('Name').count()['Distance'],
        '12 hpi':df_filter[df_filter.Condition=='12 hpi'].groupby('Name').count()['Distance']
        }
data = outlier_removal(data, zscore)
data = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in data.items()]))
data.to_csv('csv'+os.sep+'Foci number.csv')
make_swarmplot(data, ax, cell_dotsize, pairs)
#ax.set_title('Number of PLA foci by cell')
ax.set_ylabel('Foci number')

# pairs = [['NI', '4 hpi'],
#          ['NI', '8 hpi'],
#          ['NI', '12 hpi'],
#          ['8 hpi', '12 hpi']]

ax = fig.add_subplot(413)
data = {
        'NI':df_filter[df_filter.Condition=='NI'].groupby('Name').mean()['Distance'],
        '4 hpi':df_filter[df_filter.Condition=='4 hpi'].groupby('Name').mean()['Distance'],
        '8 hpi':df_filter[df_filter.Condition=='8 hpi'].groupby('Name').mean()['Distance'],
        '12 hpi':df_filter[df_filter.Condition=='12 hpi'].groupby('Name').mean()['Distance']
        }
data = outlier_removal(data, zscore)
data = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in data.items()]))
make_swarmplot(data, ax, cell_dotsize, pairs)
#ax.set_title('PLA foci distance to mito by cell')
ax.set_ylabel('PLA distance to mitochondria (${\mu m}$)')



#pairs = [['8 hpi', '12 hpi']]
ax = fig.add_subplot(414)
data = {
        'NI':df_filter[df_filter.Condition=='NI'].groupby('Name').mean()['area'],
        '4 hpi':df_filter[df_filter.Condition=='4 hpi'].groupby('Name').mean()['area'],
        '8 hpi':df_filter[df_filter.Condition=='8 hpi'].groupby('Name').mean()['area'],
        '12 hpi':df_filter[df_filter.Condition=='12 hpi'].groupby('Name').mean()['area']
        }
data = outlier_removal(data, zscore)
data = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in data.items()]))
data.to_csv('csv'+os.sep+'Volume PLA foci.csv')
make_swarmplot(data, ax, cell_dotsize, pairs)
#ax.set_title('PLA foci volume by cell')
ax.set_ylabel('Volume PLA foci (${\mu m^3}$)')



plt.tight_layout()
plt.savefig('Figure-graph.png', dpi=600)
plt.savefig('Figure-graph.svg')
plt.show()
