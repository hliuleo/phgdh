# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 15:09:03 2017

@author: hliu
""" 

import os
import sys
from researchcode.plotting.plot_set import *
from researchcode.plotting.plot_set import *
import MDAnalysis as md


code_pth = os.environ['phgdh']+os.sep+'code'
sys.path.append(code_pth)
from proj_config import *


##############################################
#          LOAD Top and Traj                #

protTop = os.sep.join([STRUCTPATH, 'protein.pdb'])
prot_NAD_Top = os.sep.join([STRUCTPATH, 'protein_NAD.pdb'])
prot_NAD_LIG_Top = os.sep.join([STRUCTPATH, 'protein_NAD_LIG2201.pdb'])

apoTraj_Info = {'main': 'HLIU_Prot_CoF',
               'num': 3,
               'top': prot_NAD_Top}
holoTraj_Info = {'main': 'HLIU_Prot_CoF_LIG2201',
                'num': 3,
                'top': prot_NAD_LIG_Top}

def loadTrajs(info):
    f_pth = os.sep.join([TRAJPATH, info['main']])
    num = info['num']
    trajs = {}
    top = info['top']
    traj_names = []
    for i in range(1, num+1):
        trj = os.sep.join([f_pth, '{}/{}.xtc'.format(i, i)])
        trajs[i] = md.Universe(top, trj)
        traj_names.append(trj)
    trajs['all'] = md.Universe(top, traj_names)
    return trajs


def writeTrajs(u, outName, top=None, beg=0, end=-1, skip=1):
    W = md.Writer(outName)
    if top is not None:
        selection = u.select_atoms(top)

        for ts in u.trajectory[beg:end:skip]:
            W.write(selection)
    else:
        for ts in u.trajectory[beg:end:skip]:
            W.write(u)
     
apoTrajs = loadTrajs(apoTraj_Info)
holoTrajs = loadTrajs(holoTraj_Info)

#########################################
#    RMSF                              #
from MDAnalysis.analysis import rms


apoRMSFs = {}
holoRMSFs = {}
for n in apoTrajs:
    u = apoTrajs[n]
    protein = u.select_atoms("protein")
    R = rms.RMSF(protein)
    R.run()
    apoRMSFs[n] = R.rmsf

#########################################
#    RMSD                              #

from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

chainA_btm = and_str.join([build_new_top(and_str, ['chA', 'protein', 'btmSecStruct']),
                           'backbone'])
chainA_top = and_str.join([build_new_top(and_str, ['chA', 'protein', 'topSecStruct']),
                           'backbone'])
chainB_btm = and_str.join([build_new_top(and_str, ['chB', 'protein', 'btmSecStruct']),
                           'backbone'])
chainB_top = and_str.join([build_new_top(and_str, ['chB', 'protein', 'topSecStruct']),
                           'backbone'])
                           


ref = md.Universe(prot_NAD_Top)
protein_ref = ref.select_atoms('protein and backbone')
chainA_ref = ref.select_atoms('segid A and protein and backbone')
chainB_ref = ref.select_atoms('segid B and protein and backbone')

chainA = 'segid A and protein and backbone'
chainB = 'segid B and protein and backbone'
chainA_btm = build_new_top(and_str, ['chA', 'protein', 'btm'])+' and backbone'
chainA_top = build_new_top(and_str, ['chA', 'protein', 'top'])+' and backbone'
chainB_btm = build_new_top(and_str, ['chB', 'protein', 'btm'])+' and backbone'
chainB_top = build_new_top(and_str, ['chB', 'protein', 'top'])+' and backbone'
chainA_btm_ref = ref.select_atoms(chainA_btm)
chainA_top_ref = ref.select_atoms(chainA_top)
chainB_btm_ref = ref.select_atoms(chainB_btm)
chainB_top_ref = ref.select_atoms(chainB_top)

apoRMSD = {}

for n in apoTrajs:
    u = apoTrajs[n]
#    md.analysis.align.A(u, ref, select='protein and backbone')

    R = md.analysis.rms.RMSD(u, ref, select='segid A and protein and backbone',
                             groupselections=[chainA,
                                              chainA_btm,
                                              chainA_top])
    R.run()
    apoRMSD[n] = R.rmsd

def aver(l):
        #global skip
        start = 0
        skip = 10
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


for n in apoRMSD:
    if n != 'all':
        d = apoRMSD[n][:,5]
        y = aver(d)
        x = range(len(y))
        plt.plot(x, y, label='traj %d'%n)
        plt.xlabel('Simulation time (ns)')
        plt.ylabel('RMSD')
        plt.legend(frameon=False)