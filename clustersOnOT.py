#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 16:45:11 2018

@author: hliu
"""

ida = sio.loadmat('/Volumes/Data/Projects/PHGDH/analysis/IDA_scaled.mat')
k = []
for c in ['apo', 'holo']:
    for t in ['A', 'B']:
        k.append(ida['{}_{}'.format(c, t)])
d = np.concatenate(k, axis=1)
ida_df = pd.DataFrame(d)
ida_df.columns=pd.MultiIndex.from_product([['apo', 'holo'], ['A', 'B'], ['open', 'twist']])

cond='apo'
chain = 'B'
num_macro = 6
for i in range(num_macro):
    centers = np.loadtxt('macro_%d_cluster_centers_index.dat'%(i+1))
    sns.kdeplot(ida_df[cond][chain].open, ida_df[cond][chain].twist,cmap='Blues',shade=True)
    plt.scatter(0,0,s=60,c='r')
    
    points = ida_df[cond][chain].ix[centers][:]
    plt.scatter(points.open, points.twist,s=60,c='g')
    plt.savefig('macro%d_onOT.png'%(i+1))

plt.figure()
colors = ['r','orange','y','w','c','b']
sns.kdeplot(ida_df[cond][chain].open, ida_df[cond][chain].twist,cmap='Blues',shade=True)
for i in range(num_macro):
    lines = open('macro_%d_cluster_centers_index.dat'%(i+1),'r').readlines()
    idx = int(lines[0].split(' ')[-1])
    points = ida_df[cond][chain].ix[idx][:]
    plt.scatter(points.open, points.twist,s=60, label='c%d'%(i+1), c=colors[i])
plt.legend()
plt.savefig('macroCenter_onOT.png')

    #ffff14149393