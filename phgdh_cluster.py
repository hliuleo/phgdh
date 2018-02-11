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
#from p53_cluster import *

def prep_cluster(ref=None, atom_idxs=None):

    featurizer = RawPositionsFeaturizer(atom_idxs, ref)
    scaler = RobustScaler()

    if os.path.exists('rawPs'):
        os.system('rm -rf rawPs')
    rawPs = xyz.fit_transform_with(featurizer, 'rawPs', fmt='dir-npy')

    if os.path.exists('scaled_rawPs'):
        os.system('rm -rf scaled_rawPs')
    scaled_rawPs = rawPs.fit_transform_with(scaler, 'scaled_rawPs/', fmt='dir-npy')

    print(scaled_rawPs[0].shape[1]) 

    ss = []
    for i in range(1, scaled_rawPs[0].shape[1]+1):
        tica_model = tICA(lag_time=1, n_components=i)
        tica_model = scaled_rawPs.fit_with(tica_model)
        ss.append(tica_model.score(scaled_rawPs))
    plt.plot(range(len(ss)), ss)
    plt.savefig('tica_n_comp.png')
    return scaled_rawPs


def run_micro_cluster(scaled_rawPs, tica_dim=None, micro_num=None, dt=None,
                      lagtimes=list(range(1, 1000, 50))):

    tica_model = tICA(lag_time=1, n_components=tica_dim)
    tica_model = scaled_rawPs.fit_with(tica_model)
    if os.path.exists('ticas'):
        os.system('rm -rf ticas')
    tica_trajs = scaled_rawPs.transform_with(tica_model, 'ticas/', fmt='dir-npy')

    txx = np.concatenate(tica_trajs)
    clustered_trajs, clusterer = clusterData(tica_trajs, micro_num)

    plt.figure()
    drawMicroCluster(txx, clusterer)
    plt.savefig('micro_cluster.png')
#    micro_index = np.concatenate(clustered_trajs)
#    micro_pop = getPop(micro_index)

    plt.figure()
    
    timescales, n_timescales = calc_ImpliedTimescale(lagtimes, clustered_trajs)
#    timescale = timescales*dt
    draw_ImpliedTimescale(timescales, n_timescales)
    plt.savefig('impliedT.png')
    #plt.ylim([5,10000])

    plt.figure()
    msm = MarkovStateModel(lag_time=1, n_timescales=20)
    msm.fit(clustered_trajs)
    assignments = clusterer.partial_transform(txx)
    assignments = msm.partial_transform(assignments)
    draw_MicroCluster_FreeEnergy(txx, msm, clusterer, assignments)
    plt.savefig('micro_cluster_eigen.png')
    return clustered_trajs, clusterer, txx


def run_macro_cluster(clustered_trajs, clusterer, txx):

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
    return macro_clusters


def save_macro_clusters(xyz_all, macro, clusterer):
    micro_centerIdx = np.array([find_clusterCenter_IDX(c, txx)[0] for c in clusterer.cluster_centers_])
    saveCenter(xyz_all, macro.microstate_mapping_, clusterer.counts_, micro_centerIdx, './')


xyz = dataset(os.sep.join([TRAJPATH, 'HLIU_Prot_CoF/', '[0-9]_dropF50.xtc']),
              topology=os.sep.join([STRUCTPATH,'protein_NAD.pdb']))
xyz_all = dataset(os.sep.join([TRAJPATH, 'HLIU_Prot_CoF/', '1-3_*.xtc']),
                  topology=os.sep.join([STRUCTPATH,'protein_NAD.pdb']))


ref_xyz = dataset(os.sep.join([STRUCTPATH,'protein_NAD.pdb']),
              topology=os.sep.join([STRUCTPATH,'protein_NAD.pdb']))

ref = ref_xyz[0]

nad_B_CA_atom_idxs = ref.top.select('name CA and chainid 1 and (residue 77 or residue 206 or residue 174 or residue 211 or residue 154 or residue 233 or residue 259)')
mlt_B_CA_atom_idxs = ref.top.select('name CA and (chainid 1 and (residue 77 or residue 53 or residue 235 or residue 282 or residue 54) or (chainid 0 and residue 134))')

nad_B_BB_atom_idxs = ref.top.select('backbone and chainid 1 and (residue 77 or residue 206 or residue 174 or residue 211 or residue 154 or residue 233 or residue 259)')
mlt_B_BB_atom_idxs = ref.top.select('backbone and (chainid 1 and (residue 77 or residue 53 or residue 235 or residue 282 or residue 54) or (chainid 0 and residue 134))')

nad_BB_atom_idxs = ref.top.select('backbone and (residue 77 or residue 206 or residue 174 or residue 211 or residue 154 or residue 233 or residue 259)')
mlt_BB_atom_idxs = ref.top.select('backbone and ((residue 77 or residue 53 or residue 235 or residue 282 or residue 54) or (residue 134))')

nad_A_BB_atom_idxs = ref.top.select('backbone and chainid 0 and (residue 77 or residue 206 or residue 174 or residue 211 or residue 154 or residue 233 or residue 259)')
mlt_A_BB_atom_idxs = ref.top.select('backbone and (chainid 0 and (residue 77 or residue 53 or residue 235 or residue 282 or residue 54) or (chainid 1 and residue 134))')

A_CA_atom_idxs = ref.top.select('name CA and chainid 0 and residue 7 to 305')
B_CA_atom_idxs = ref.top.select('name CA and chainid 1 and residue 7 to 305')

atom_idxs = A_CA_atom_idxs

prep_parms = dict(ref=ref,
                  atom_idxs=atom_idxs)

scaled_rawPs = prep_cluster(**prep_parms)

micro_parms = dict(tica_dim=50,
                   dt=100,
                   lagtimes=list(range(1,1000,50)),
                   micro_num=100)

clustered_trajs, clusterer, txx = run_micro_cluster(scaled_rawPs, **micro_parms)

macro_clusters = run_macro_cluster(clustered_trajs, clusterer, txx)

save_macro_clusters(xyz_all, macro_clusters[6], clusterer)

def stat_trj_source(macro_trajs, cluster_num):
    results = {}
    for i in range(cluster_num):
        num_of_each_traj = [sum(t==i) for t in macro_trajs]
        results[i+1] = np.array(num_of_each_traj)/sum(num_of_each_traj)*100

    return pd.DataFrame(results)

##################################
macro_num = 6
savePth = '{}_macro_trajs'.format(macro_num)
if os.path.exists(savePth):
    os.system('rm -rf {}'.format(savePth))
os.mkdir(savePth)
saveTrj(xyz_all, macro_clusters[macro_num].macro_index, savePth)