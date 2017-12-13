# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 10:46:01 2017

@author: hliu
"""
import os
import pandas as pd
from researchcode.plotting.plot_set import *

def xvg2Dat(filePath):
    specialLetters = ['#', '@']
    s = open(filePath, 'r').readlines()
    raw_data = [i for i in s if i[0] not in specialLetters]
    dim = len(raw_data[0].strip().split())
    data = np.fromstring(''.join(raw_data), sep='\n')
    data = data.reshape((len(data)/dim, dim))
    df = pd.DataFrame(data)
    return df


analysis = ['rmsf', 'pca', 'interdomain_motion']
conditions = ['apo', 'holo']
tops = ['A', 'B']
trajs = ['1', '2', '3', '1-3']

levels = (conditions, tops, trajs)

files = {'rmsf': 'rmsf.xvg',
         'pca': 'pc1-2.xvg',
         'hbond_dist': 'dist.dat',
         'hole_size': 'hole_rg.dat',
         'hole_contact': 'hole_contact.dat',
         'interdomain_motion': 'coords.dat'}
         
def buildLayout(folders, data):
    return 0

analysisFolderDefaultLayout = {'rmsf': levels.append('rmsf.xvg'),
                               'pca': levels.append('pc1-2.xvg'),
                               'hbond_dist': levels.append('dist.dat'),
                               'interdomain_motion': levels.append('coords.dat'),
                               'hole_contact': levels.append('hole_contact.dat'),
                               'hole_size': levels.append('hole_rg.dat')}

data = {}

AnaHOME = '/home/hliu/Projects/PHGDH/analysis'

def loadAnaData(ana, levels):
    data = []
    conds, tops, trajs = levels
    for c in conds:
        for t in tops:
            for traj in trajs:
                f = files[ana]
                f_path = os.sep.join([AnaHOME, ana, c, t, traj, f])
                if f.split('.')[-1] == 'xvg':
                    data.append(xvg2Dat(f_path))
                else:
                    data.append(pd.read_table(f_path, sep=' ', header=None))
    return data


def mergeAnaData(data, col):

    data = pd.concat(data,axis=1)
    data_cols = list(levels) + [col]
    data_cols = pd.MultiIndex.from_product(data_cols)
    data.columns = data_cols
    return data

rmsf = mergeAnaData(loadAnaData('rmsf', levels), ['resid', 'rmsf'])
pca = mergeAnaData(loadAnaData('pca', levels), ['pc1','pc2'])
IDA = mergeAnaData(loadAnaData('interdomain_motion', levels), ['open','twist'])
hdist = mergeAnaData(loadAnaData('hbond_dist', levels), ['h1', 'h2', 'h3'])

hsize = mergeAnaData(loadAnaData('hole_size', levels), ['G1', 'G2', 'G3', 'G4'])
hcontact = mergeAnaData(loadAnaData('hole_contact', levels), ['G1', 'G2', 'G3', 'G4'])

#################
# Hole PLOTTING #
################
colors = ['b', 'g', 'r', 'k']

def set_yaxis_majortick_interval(ax, interval):
    ticklocs = ax.yaxis.get_ticklocs()
    maxvalue = ticklocs[-1]
    current = ticklocs[0]
    newtickloc = []
    newticklabel = []
    newtickloc.append(current)
    newticklabel.append(str(current))
    while current < maxvalue:
        current += interval
        newtickloc.append(current)
        newticklabel.append(str(current))
    ax.set_yticks(newtickloc,newticklabel)


def aver(l):
        #global skip
        start = 0
        skip = 100
        smooth = []
        while 1:
                if (start + skip) < len(l):
                        a = float(sum(l[start:start+skip])/skip)
                        start = start + skip
                        smooth.append(a)
                else:
                        a = float(sum(l[start:len(l)])/(len(l)-start))
                        smooth.append(a)
                        break
        smooth = np.array(smooth).reshape((len(smooth),1))
        return smooth


for c in conditions:
    f, axarr = plt.subplots(4,2,sharex=True)
    for i, t in enumerate(tops):
        for j, n in enumerate(['G1', 'G2', 'G3', 'G4']):
            for k, traj in enumerate(trajs[:-1]):
                y = aver(hsize[c][t][traj][n].dropna())
                x = np.arange(len(y))
                axarr[j,i].plot(x, y, label=str(traj), c=colors[k])
                axarr[j,i].legend(frameon=False)
                set_yaxis_majortick_interval(axarr[j,i], 1)
            axarr[j,0].set_ylabel(n)
        axarr[3,i].set_xlabel(t)
    plt.subplots_adjust(top=0.95, bottom=0.15, left=0.10, right=0.95, hspace=0.1,
                    wspace=0.3)
    f.savefig('_'.join([c,'.tif']))


for c in conditions:
    f, axarr = plt.subplots(4,2,sharex=True)
    for i, t in enumerate(tops):
        for j, n in enumerate(['G1', 'G2', 'G3', 'G4']):
            for k, traj in enumerate(trajs[:-1]):
                y = aver(hcontact[c][t][traj][n].dropna())
                x = np.arange(len(y))
                axarr[j,i].plot(x, y, label=str(traj), c=colors[k])
                axarr[j,i].legend(frameon=False)
                #set_yaxis_majortick_interval(axarr[j,i], 1)
            axarr[j,0].set_ylabel(n)
        axarr[3,i].set_xlabel(t)
    plt.subplots_adjust(top=0.95, bottom=0.15, left=0.10, right=0.95, hspace=0.1,
                    wspace=0.3)
    f.savefig('_'.join([c,'.tif']))

for c in conditions:
    for t in tops:
        f, axarr = plt.subplots(4,4)
        for i, traj in enumerate(trajs):
            for j, n in enumerate(['G1', 'G2', 'G3', 'G4']):
                x = hsize[c][t][traj][n].dropna()
                y = hcontact[c][t][traj][n].dropna()
                fit = np.polyfit(x, y, deg=1)
                corrcoef = np.corrcoef(x,y)[0][1]
                axarr[i,j].plot(x, fit[0]*x+fit[1],c='r')
                axarr[i,j].scatter(x, y, c='k', label=str('%.2f'%corrcoef))
                axarr[i,j].yaxis.set_ticklabels([])
                axarr[i,j].xaxis.set_ticklabels([])
                axarr[i,j].legend(frameon=False, loc=1)
            axarr[i,0].set_ylabel(traj)
        f.savefig('_'.join([c,t,'.tif']))

ind = np.arange(4)
width=0.25

data = hsize
holob = data.holo.B['1-3'].mean()
holoberr = data.holo.B['1-3'].std()
holoa = data.holo.A['1-3'].mean()
holoaerr = data.holo.A['1-3'].std()
apob = data.apo.B['1-3'].mean()
apoberr = data.apo.B['1-3'].std()

plt.bar(ind, holob,width,color='r',yerr=holoberr,ecolor='k', label='Holo B')
plt.bar(ind+width,holoa,width,yerr=holoaerr,ecolor='k', label='Holo A')
plt.legend(frameon=False)
plt.xticks([0.25,1.25,2.25,3.25],['G1','G2','G3','G4'])



plt.bar(ind, holob,width,color='r',yerr=holoberr,ecolor='k',label='Holo B')
plt.bar(ind+width,apob,width,yerr=apoberr,ecolor='k',label='Apo B')
plt.legend(frameon=False)
plt.xticks([0.25,1.25,2.25,3.25],['G1','G2','G3','G4'])

#################
# RMSF PLOTTING #
################

# traj compares

for c in conditions:
    for t in tops:
        plt.figure()
        for traj in trajs[:-1]:
            d = rmsf[c][t][traj]
            plt.plot(d.resid, d.rmsf*10, label=traj)
        plt.xlabel('resid')
        plt.ylabel('RMSF (A)')
        plt.tight_layout()
        plt.legend(frameon=False)
        plt.savefig('_'.join([c, t, '.tif']))

for c in conditions:
    for traj in trajs:
        plt.figure()
        for t in tops:
            d = rmsf[c][t][traj]
            plt.plot(d.resid, d.rmsf*10, label=t)
        plt.xlabel('resid')
        plt.ylabel('RMSF (A)')
        plt.tight_layout()
        plt.legend(frameon=False)
        plt.savefig('_'.join([c, traj, '.tif']))

for t in tops:
    plt.figure()
    for c in conditions:
        d = rmsf[c][t]['1-3']
        plt.plot(d.resid, d.rmsf*10, label=c)
    plt.xlabel('resid')
    plt.ylabel('RMSF (A)')
    plt.tight_layout()
    plt.legend(frameon=False)
    plt.savefig('_'.join([t, '.tif']))

#################
# PCA PLOTTING #
###############
colors = ['b', 'g', 'r']

for c in conditions:
    for t in tops:
        plt.figure()
        for idx, traj in enumerate(trajs[:-1]):
            d = pca[c][t][traj]
            plt.scatter(d['pc1'], d['pc2'], label=traj, c=colors[idx])
        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.tight_layout()
        plt.legend(frameon=False)
        plt.savefig('_'.join([c, t, '.tif']))

#################
# IDA PLOTTING #
###############

scatter_settings = {'s': 60,
                    'alpha': 0.4}

for c in conditions:
    for t in tops:
        plt.figure()
        for idx, traj in enumerate(trajs[:-1]):
            d = IDA[c][t][traj]
            plt.scatter(d['open'], d['twist'], label=traj, c=colors[idx],**scatter_settings)
        plt.xlabel('open')
        plt.ylabel('twist')
        plt.tight_layout()
        plt.legend(frameon=False)
        plt.savefig('_'.join([c, t, '.tif']))

for t in tops:
    plt.figure()
    for idx, c in enumerate(conditions):
        d = IDA[c][t]['1-3']
        plt.scatter(d['open'], d['twist'], label=c, c=colors[idx],**scatter_settings)
    plt.xlabel('open')
    plt.ylabel('twist')
    plt.tight_layout()
    plt.legend(frameon=False)
    plt.savefig('_'.join([t, '.tif']))


for c in conditions:
    for t in tops:
        for traj in trajs:
            plt.figure()
            d = IDA[c][t][traj].dropna()
            x = d['open']
            y = d['twist']
            fit = np.polyfit(x, y, deg=1)
            corrcoef = np.corrcoef(x,y)[0][1]
            plt.plot(x, fit[0]*x+fit[1],c='r')
            plt.scatter(x, y, c='k', label=str('%.2f'%corrcoef),**scatter_settings)
            plt.xlabel('open')
            plt.ylabel('twist')
            plt.tight_layout()
            plt.legend(frameon=False)
            plt.savefig('_'.join([c, t, traj, '.tif']))

for c in conditions:
    for traj in trajs:
        plt.figure()
        for idx, t in enumerate(tops):
            d = IDA[c][t][traj]
            plt.scatter(d['open'], d['twist'], label=t, c=colors[idx],**scatter_settings)
        plt.xlabel('open')
        plt.ylabel('twist')
        plt.tight_layout()
        plt.legend(frameon=False)
        plt.savefig('_'.join([c, traj, '.tif']))

#####################
# IDA PCA PLOTTING #
###################

for c in conditions:
    for t in tops:
        for traj in trajs:
            f, axarr = plt.subplots(2,2)
            for i, a in enumerate(['pc1','pc2']):
                for j, d in enumerate(['open','twist']):
                    x = pca[c][t][traj][a].dropna()
                    y = IDA[c][t][traj][d].dropna()
                    fit = np.polyfit(x, y, deg=1)
                    corrcoef = np.corrcoef(x,y)[0][1]
                    axarr[i,j].plot(x, fit[0]*x+fit[1],c='r')
                    axarr[i,j].scatter(x, y, c='k', label=str('%.2f'%corrcoef))
                    axarr[1,j].set_xlabel(d)
                    axarr[i,0].set_ylabel(a)
                    axarr[i,j].legend(frameon=False)
            f.savefig('_'.join([c,t,traj,'.tif']))

#####################
# hDist IDA PLOTTING #
###################

for c in conditions:
    for t in tops:
        for traj in trajs:
            f, axarr = plt.subplots(2,3)
            for i, a in enumerate(['open','twist']):
                for j, d in enumerate(['h1','h2','h3']):
                    x = IDA[c][t][traj][a].dropna()
                    y = hdist[c][t][traj][d].dropna()
                    fit = np.polyfit(x, y, deg=1)
                    corrcoef = np.corrcoef(x,y)[0][1]
                    axarr[i,j].plot(x, fit[0]*x+fit[1],c='r')
                    axarr[i,j].scatter(x, y, c='k', label=str('%.2f'%corrcoef))
                    axarr[1,j].set_xlabel(d)
                    axarr[i,0].set_ylabel(a)
                    axarr[i,j].legend(frameon=False)
            f.savefig('_'.join([c,t,traj,'.tif']))

#####################
# hDist PCA PLOTTING #
###################

for c in conditions:
    for t in tops:
        for traj in trajs:
            f, axarr = plt.subplots(2,3)
            for i, a in enumerate(['pc1','pc2']):
                for j, d in enumerate(['h1','h2','h3']):
                    x = pca[c][t][traj][a].dropna()
                    y = hdist[c][t][traj][d].dropna()
                    fit = np.polyfit(x, y, deg=1)
                    corrcoef = np.corrcoef(x,y)[0][1]
                    axarr[i,j].plot(x, fit[0]*x+fit[1],c='r')
                    axarr[i,j].scatter(x, y, c='k', label=str('%.2f'%corrcoef))
                    axarr[1,j].set_xlabel(d)
                    axarr[i,0].set_ylabel(a)
                    axarr[i,j].legend(frameon=False)
            f.savefig('_'.join([c,t,traj,'.tif']))

### for mlt conds

conditions = ['apo_mlt']
tops = ['A', 'B']
trajs = ['1', '2', '3', '4']



levels = (conditions, tops, trajs)

rmsf2 = mergeAnaData(loadAnaData('rmsf', levels), ['resid', 'rmsf'])
IDA2 = mergeAnaData(loadAnaData('interdomain_motion', levels), ['open','twist'])
hsize2 = mergeAnaData(loadAnaData('hole_size', levels), ['G1', 'G2', 'G3', 'G4'])
hcontact2 = mergeAnaData(loadAnaData('hole_contact', levels), ['G1', 'G2', 'G3', 'G4'])



for c in conditions:
    f, axarr = plt.subplots(4,2,sharex=True)
    for i, t in enumerate(tops):
        for j, n in enumerate(['G1', 'G2', 'G3', 'G4']):
            for k, traj in enumerate(trajs):
                y = aver(hsize2[c][t][traj][n].dropna())
                x = np.arange(len(y))
                axarr[j,i].plot(x, y, label=str(traj), c=colors[k])
                axarr[j,i].legend(frameon=False)
                set_yaxis_majortick_interval(axarr[j,i], 1)
            axarr[j,0].set_ylabel(n)
        axarr[3,i].set_xlabel(t)
    plt.subplots_adjust(top=0.95, bottom=0.15, left=0.10, right=0.95, hspace=0.1,
                    wspace=0.3)
    f.savefig('_'.join([c,'.tif']))


for c in conditions:
    f, axarr = plt.subplots(4,2,sharex=True)
    for i, t in enumerate(tops):
        for j, n in enumerate(['G1', 'G2', 'G3', 'G4']):
            for k, traj in enumerate(trajs):
                y = aver(hcontact2[c][t][traj][n].dropna())
                x = np.arange(len(y))
                axarr[j,i].plot(x, y, label=str(traj), c=colors[k])
                axarr[j,i].legend(frameon=False)
                #set_yaxis_majortick_interval(axarr[j,i], 1)
            axarr[j,0].set_ylabel(n)
        axarr[3,i].set_xlabel(t)
    plt.subplots_adjust(top=0.95, bottom=0.15, left=0.10, right=0.95, hspace=0.1,
                    wspace=0.3)
    f.savefig('_'.join([c,'.tif']))

for c in conditions:
    for t in tops:
        f, axarr = plt.subplots(4,4)
        for i, traj in enumerate(trajs):
            for j, n in enumerate(['G1', 'G2', 'G3', 'G4']):
                x = hsize2[c][t][traj][n].dropna()
                y = hcontact2[c][t][traj][n].dropna()
                fit = np.polyfit(x, y, deg=1)
                corrcoef = np.corrcoef(x,y)[0][1]
                axarr[i,j].plot(x, fit[0]*x+fit[1],c='r')
                axarr[i,j].scatter(x, y, c='k', label=str('%.2f'%corrcoef))
                axarr[i,j].yaxis.set_ticklabels([])
                axarr[i,j].xaxis.set_ticklabels([])
                axarr[i,j].legend(frameon=False, loc=1)
            axarr[i,0].set_ylabel(traj)
        f.savefig('_'.join([c,t,'.tif']))

#################
# RMSF PLOTTING #
################

# traj compares

for c in conditions:
    for t in tops:
        plt.figure()
        for traj in trajs:
            d = rmsf2[c][t][traj]
            plt.plot(d.resid, d.rmsf*10, label=traj)
        plt.xlabel('resid')
        plt.ylabel('RMSF (A)')
        plt.tight_layout()
        plt.legend(frameon=False)
        plt.savefig('_'.join([c, t, '.tif']))

for c in conditions:
    for traj in trajs:
        plt.figure()
        for t in tops:
            d = rmsf2[c][t][traj]
            plt.plot(d.resid, d.rmsf*10, label=t)
        plt.xlabel('resid')
        plt.ylabel('RMSF (A)')
        plt.tight_layout()
        plt.legend(frameon=False)
        plt.savefig('_'.join([c, traj, '.tif']))

#################
# IDA PLOTTING #
###############

colors = ['b', 'g', 'r', 'k']

scatter_settings = {'s': 60,
                    'alpha': 0.4}

for c in conditions:
    for t in tops:
        plt.figure()
        for idx, traj in enumerate(trajs):
            d = IDA2[c][t][traj]
            plt.scatter(d['open'], d['twist'], label=traj, c=colors[idx],**scatter_settings)
        plt.xlabel('open')
        plt.ylabel('twist')
        plt.tight_layout()
        plt.legend(frameon=False)
        plt.savefig('_'.join([c, t, '.tif']))

for t in tops:
    plt.figure()
    for idx, c in enumerate(conditions):
        d = IDA[c][t]['1-3']
        plt.scatter(d['open'], d['twist'], label=c, c=colors[idx],**scatter_settings)
    plt.xlabel('open')
    plt.ylabel('twist')
    plt.tight_layout()
    plt.legend(frameon=False)
    plt.savefig('_'.join([t, '.tif']))


for c in conditions:
    for t in tops:
        for traj in trajs:
            plt.figure()
            d = IDA[c][t][traj].dropna()
            x = d['open']
            y = d['twist']
            fit = np.polyfit(x, y, deg=1)
            corrcoef = np.corrcoef(x,y)[0][1]
            plt.plot(x, fit[0]*x+fit[1],c='r')
            plt.scatter(x, y, c='k', label=str('%.2f'%corrcoef),**scatter_settings)
            plt.xlabel('open')
            plt.ylabel('twist')
            plt.tight_layout()
            plt.legend(frameon=False)
            plt.savefig('_'.join([c, t, traj, '.tif']))

for c in conditions:
    for traj in trajs:
        plt.figure()
        for idx, t in enumerate(tops):
            d = IDA[c][t][traj]
            plt.scatter(d['open'], d['twist'], label=t, c=colors[idx],**scatter_settings)
        plt.xlabel('open')
        plt.ylabel('twist')
        plt.tight_layout()
        plt.legend(frameon=False)
        plt.savefig('_'.join([c, traj, '.tif']))