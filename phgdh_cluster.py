#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 21:36:05 2018

@author: hliu
"""

import os
import sys
from msmbuilder.featurizer import RawPositionsFeaturizer

code_pth = os.environ['phgdh']+os.sep+'code'
sys.path.append(code_pth)
from proj_config import *
from p53_cluster import *


xyz = dataset(os.sep.join([TRAJPATH, 'HLIU_Prot_CoF/', '[0-9]_dropF50.xtc']),
              topology=os.sep.join([STRUCTPATH,'protein_NAD.pdb']))
ref = md.load(os.sep.join([STRUCTPATH,'protein_NAD.pdb']),
              top=os.sep.join([STRUCTPATH,'protein_NAD.pdb']))
atom_idxs = ref.top.select('chainid 0 and (residue 77 or residue 206 or residue 174 or residue 211 or residue 154 or residue 233 or residue 259)')
featurizer = RawPositionsFeaturizer(atom_idxs, ref)
rawPs = xyz.fit_transform_with(featurizer, 'rawPs', fmt='dir-npy')
scaler = RobustScaler()
scaled_rawPs = rawPs.fit_transform_with(scaler, 'scaled_rawPs/', fmt='dir-npy')
tica_model = tICA(lag_time=1, n_components=100)
tica_model = scaled_rawPs.fit_with(tica_model)
tica_model.score(scaled_rawPs)
xyz_all = dataset(os.sep.join([TRAJPATH, 'Prot_Lig/', '1-3_*.xtc']),
                  topology=os.sep.join([STRUCTPATH,'protein_NAD.pdb']))


dt = 10
micro_num = 100
tica_dim = 2
tica_trajs = featurizeData(xyz, tica_dim)
txx = np.concatenate(tica_trajs)
clustered_trajs, clusterer = clusterData(tica_trajs, micro_num)



drawMicroCluster(txx, clusterer)
micro_index = np.concatenate(clustered_trajs)
micro_pop = getPop(micro_index)

lagtimes = list(range(1, 200, 50))
timescales, n_timescales = calc_ImpliedTimescale(lagtimes, clustered_trajs)
timescale = timescales*dt 
draw_ImpliedTimescale(timescales, n_timescales)
plt.ylim([5,10000])

msm = MarkovStateModel(lag_time=1, n_timescales=20)
msm.fit(clustered_trajs)
assignments = clusterer.partial_transform(txx)
assignments = msm.partial_transform(assignments)
draw_MicroCluster_FreeEnergy(txx, msm, clusterer, assignments)

msm = MarkovStateModel(lag_time=150, n_timescales=20)
msm.fit(clustered_trajs)

for i in range(2, 11):
    fig = plt.figure(i)
    n_timescales = i
    check_MacroNum(msm, n_timescales)
    fig.savefig('%d_timescales' % i)

msm_build_macro = MarkovStateModel(lag_time=1, n_timescales=20)
msm_build_macro.fit(clustered_trajs)
assignments = clusterer.partial_transform(txx)
assignments = msm_build_macro.partial_transform(assignments)
macro_clusters = {}
pcca_clusters = {}
for i in range(2, 11):
    fig = plt.figure(i)
    pcca = PCCAPlus.from_msm(msm_build_macro, i)
    pcca_clusters[i] = pcca
    macro_clusters[i] = MacroCluster(pcca, i, clustered_trajs)
    macro_clusters[i].draw_MacroCluster_FreeEnergy(txx, clusterer, pcca_clusters[i], assignments)
    fig.savefig('%d_macrostates.png' % i)

def stat_trj_source(macro_trajs, cluster_num):
    results = {}
    for i in range(cluster_num):
        num_of_each_traj = [sum(t==i) for t in macro_trajs]
        results[i+1] = np.array(num_of_each_traj)/sum(num_of_each_traj)*100

    return pd.DataFrame(results)

    
macro_num = 6
savePth = '{}_macro_trajs'.format(macro_num)
if os.path.exists(savePth):
    os.system('rm -rf {}'.format(savePth))
os.mkdir(savePth)
saveTrj(xyz_all, macro_clusters[macro_num].macro_index, savePth)