# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 13:32:43 2017

@author: hliul
"""

import pandas as pd
import numpy as np


topology = dict(
                E1 = 'resid 7:10',
                E2 = 'resid 28:31',
                E3 = 'resid 49:52',
                E4 = 'resid 71:74',
                E5 = 'resid 93:96',
                H1 = 'resid 16:25',
                H2 = 'resid 37:46',
                H3 = 'resid 60:65',
                H4 = 'resid 84:90',
                
                E6 = 'resid 146:150',
                E7 = 'resid 169:173',
                E8 = 'resid 201:204',
                E9 = 'resid 228:232',
                E10 = 'resid 253:258',
                E11 = 'resid 277:279',
                
                H5 = 'resid 101:118',
                H6 = 'resid 120:128',
                H7 = 'resid 154:164',
                H8 = 'resid 179:183',
                H9 = 'resid 192:198',
                H10 = 'resid 218:223',
                H11 = 'resid 241:250',
                
                H12 = 'resid 288:305',
                
                L1 = 'resid 11:15',
                L2 = 'resid 32:36',
                L3 = 'resid 53:59',
                L4 = 'resid 75:83',
                L5 = 'resid 96:101',

                L6 = 'resid 151:153',
                L7 = 'resid 174:178',
                L8 = 'resid 205:217',
                L9 = 'resid 233:240',
                L10 = 'resid 259:270',
                L11 = 'resid 280:287',
                L12 = 'resid 130:142',)

selString2List = lambda x: [int(i) for i in x[6:].split(':')]
list2Range = lambda x: np.arange(x[0], x[1]+1)
new_top = {k:list2Range(selString2List(topology[k])) for k in topology}


def assignRegion(resid):
    for k in new_top:
        if resid in new_top[k]:
            return k
    else:
        return np.nan

def cleanLine(l):
    '''ignore hydrogen type main or side chain'''
    l = l.replace('Seg{}-', '')
    l = l.replace('%', '')
#    l = [i for i in l.split() if i != sys.argv[2]]
    return l

def ignoreTypeResName(df):
    for hType in ['D', 'A']:
        df[hType] = df[hType].str.replace(r'[A-Z]{3}', '')
        for bType in ['-Side', '-Main']:
            df[hType] = df[hType].str.replace(bType, '')
        
    return df


def sort_group_byLVG(df, gbCols, valCol):
    '''LVG means largest values in group,
       sort groups by key of their largest value in valCol,
       gbCols are group by cols,
       groups are sorted ascendingly, values with group are also ascending'''

    def fillwith(a, valCol):
        a = a.sort_values(by=valCol, ascending=False)
        a[valCol] = a[valCol].max()
        return a
    
    cols = df.columns
    df = df.groupby(gbCols)[cols].apply(lambda x:x.sort_values(valCol,ascending=False))
    df.index = list(range(len(df)))

    df.index.name = 'ID'
    df = df.reset_index()
    g_reorderIDX = df.groupby(gbCols)[[valCol, 'ID']].apply(lambda x: fillwith(x,valCol))
    g_reorderIDX = g_reorderIDX.reset_index().sort_values([valCol,'ID'],ascending=[False,True])
    g_reorderIDX = g_reorderIDX['ID'].values
    df = df.reindex(index=g_reorderIDX)[cols]
    return df

def dealWith(lines):
    '''Read lines of a hbonds details file
       clean file and save raw; then clear close contact, bond type, adding topology
       output as by topology sorted by percentage'''
    lines = [cleanLine(l) for l in lines[2:]]
    lines = [l.strip().split() for l in lines]
    raw = pd.DataFrame(lines)
    p_index = raw.columns[-1]
    raw[p_index] = raw[p_index].astype(float)
    raw = raw.sort_values(by=p_index, ascending=False)

    data = raw.copy()
    data.columns = ['D', 'A', 'percent']
    data = ignoreTypeResName(data)
    data.D = data.D.astype(int)
    data.A = data.A.astype(int)

    data['dist'] = (data.D-data.A).apply(abs)
    data = data[data.dist>4]


    data['L'] = data[['D','A']].min(axis=1)
    data['U'] = data[['D','A']].max(axis=1)
    mask = (data.L.isin(new_top['L12'])) | (data.U.isin(new_top['H12']))
    data.loc[mask, ['L','U']] = data[['L','U']][mask].values[:,[1,0]]                

    data = data[['L','U','percent']]
    data = data.groupby(['L', 'U']).aggregate(sum).reset_index().sort_values('percent',ascending=False)
    data['LR'] = data.L.apply(assignRegion)
    data['UR'] = data.U.apply(assignRegion)
    sorted_df = sort_group_byLVG(data, ['LR'], 'percent').set_index(['LR','UR'])
    return raw, data, sorted_df


conds = ['apo', 'holo']
trajs = ['1', '2', '3']
chains = ['A', 'B']

import os

out = pd.ExcelWriter('hbonds_results.xls')

for c in conds:
    for t in trajs:
        for ch in chains:

            fNames = 'hbonds-details_DD_%s.dat' % ch
            fPath = os.sep.join([os.path.abspath(os.curdir), c, t, fNames])
            lines = open(fPath, 'r').readlines()

            r, d, s = dealWith(lines)
            r.to_excel(out, '{}_{}_{}_raw'.format(c,t,ch))
            d.to_excel(out, '{}_{}_{}_clean'.format(c,t,ch))
            s.to_excel(out, '{}_{}_{}_sorted'.format(c,t,ch))
out.save()

out = pd.ExcelWriter('hbonds_byRegion.xls')

for c in conds:
    for t in trajs:
        for ch in chains:

            fNames = 'hbonds-details_DD_%s.dat' % ch
            fPath = os.sep.join([os.path.abspath(os.curdir), c, t, fNames])
            lines = open(fPath, 'r').readlines()

            r, d, s = dealWith(lines)
            s = s.groupby(['LR','UR'],axis=0)[['percent']].aggregate(np.sum).reset_index()
            s = sort_group_byLVG(s,['LR'],'percent').set_index(['LR','UR'])
            s.to_excel(out, '{}_{}_{}'.format(c,t,ch))

out.save()

out = pd.ExcelWriter('hbonds_clear1.xls')

for c in conds:
    for t in trajs:
        for ch in chains:

            fNames = 'hbonds-details_DD_%s.dat' % ch
            fPath = os.sep.join([os.path.abspath(os.curdir), c, t, fNames])
            lines = open(fPath, 'r').readlines()

            r, d, s = dealWith(lines)
            s = s.reset_index()
            s = s[~(s.U==101)] # remove bad data
            s = s.set_index(['LR','UR'])
            a = s.groupby(['LR','UR'],axis=0)[['percent']].aggregate(np.sum).reset_index()
            a = a[a.percent>1]
            a = sort_group_byLVG(a,['LR'],'percent').set_index(['LR','UR'])
            a.to_excel(out, '{}_{}_{}_r'.format(c,t,ch))
            s = s[s.percent>1]
            s.to_excel(out, '{}_{}_{}'.format(c,t,ch))
out.save()

def hbondsSum(df):

    def fillwith(a):
        a = a.sort_values(by='mean', ascending=False)
        a['mean'] = a['mean'].max()
        return a
    
    df = df.sort_values(by=['percent','L'], ascending=[False, True])
    df = df.groupby(by=['L','U']).aggregate([np.sum, np.std])['percent']
    
    t = df.reset_index()
    t = t.groupby('L')[['U', 'mean', 'std']].apply(lambda x: x.sort_values(by='mean', ascending=False))
    t = t.reset_index()[['L','U','mean','std']]
    t.index.name = 'ID'
    reorder = t.groupby('L')[['mean']].apply(lambda x:fillwith(x)).reset_index().sort_values(['mean','ID'],ascending=[False,True])
    reorder = reorder.set_index(['mean','L']).ID.values
    t = t.reindex(reorder)
    t = t.reindex(reorder).set_index(['L','U'])
    return t


hbondsA = pd.read_excel('Apo_hbonds-details.xlsx', 'A')
hbondsB = pd.read_excel('Apo_hbonds-details.xlsx', 'B')

def do(data):
    data.columns = ['D', 'A', 'percent']
    data_clear = ignoreType(data)
    data_sum = hbondsSum(data)
    data_clear_sum = hbondsSum(data_clear)
    return data_sum, data_clear_sum



A_sum, A_clear = do(hbondsA)
B_sum, B_clear = do(hbondsB)

out = pd.ExcelWriter('out.xls')
A_sum.to_excel(out, 'A')
B_sum.to_excel(out, 'B')
A_clear.to_excel(out, 'A_c')
A_clear.to_excel(out, 'B_c')

out.save()