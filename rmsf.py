# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 16:30:57 2017

@author: hliu
"""

import os
import MDAnalysis as md
import MDAnalysis.analysis.rms as rms
from researchcode.plotting.plot_set import *

WORKDIR = os.sep.join([os.environ['phgdh'], 'analysis', 'rmsf'])
TRJDIR = os.sep.join([os.environ['phgdh'], 'data', 'trj'])
STRUCTDIR = os.sep.join([os.environ['phgdh'], 'data', 'struct'])
os.chdir(WORKDIR)

'''

trajsPath = {'HL': 
                  {'traj':
                      {'apo': (os.sep.join([TRJDIR, 'HL_Prot_CoF']), 1),
                       'lig2201': (os.sep.join([TRJDIR, 'HL_Prot_CoF_LIG2201']), 1)
                      },
                   'top':
                      {'apo': os.sep.join([STRUCTDIR, 'protein_NAD.pdb']),
                       'lig2201': os.sep.join([STRUCTDIR, 'protein_NAD_LIG2201.pdb'])
                      }
                  },
             'PL':
                  {'traj':                  
                      {'apo': (os.sep.join([TRJDIR, 'Prot_Only']), 4),
                       'lig2201': (os.sep.join([TRJDIR, 'PL_Prot_LIG2201']), 3)
                      },
                  'top':
                      {'apo': os.sep.join([STRUCTDIR, 'protein.pdb']),
                       'lig2201': os.sep.join([STRUCTDIR, 'protein_LIG2201.pdb'])
                      }                    
                   }
            }


u = md.Universe(trajsPath['PL']['top']['apo'], trajsPath['PL']['traj']['apo'][0]+'/1.xtc')
chainA = u.select_atoms('segid A and not (resname NDP)')

for u in ['HL', 'PL']:
    for t in ['apo', 'lig2201']:
        top = trajsPath[u]['top'][t]
        trajPth, trajNum = trajsPath[u]['traj'][t]
        for i in range(trajNum):
            u = md.Universe(top, trajPth+'/%d.xtc' %(i+1))
            with md.Writer("%s", u.atoms.n_atoms) as W:
                for ts in u.trajectory[5001::5]:
                    W.write()
'''
    
num_skipRow = 16


class rmsf_Comparer():

    def __init__(self):
        self.trjs = {}

    def loadTrj(self, filePth, name):
        data = np.loadtxt(filePth, skiprows=num_skipRow)
        self.trjs[name] = rmsf_Comparer.splitChain(data)

    def getTrjNames(self):
        return self.trjs.keys()

    @staticmethod
    def splitChain(data):
        splitted = {}
        splitted = {'A': data[:302],
                    'B': data[302:]}
        return splitted

    def compareTrjsRMSF(self, plotting_trjs, chain='A'):
        for trj in plotting_trjs:
            plt.plot(self.trjs[trj][chain][:, 0], self.trjs[trj][chain][:, 1]*10, label=trj)
        plt.legend(loc=1, frameon=False)
        plt.xlabel('Residue')
        plt.ylabel(r'RMSF ($\AA$)')
        plt.title('chain %s' % chain)

    def compareChainRMSF(self, plotting_trj):
        for c in ['A', 'B']:
            plt.plot(self.trjs[plotting_trj][c][:, 0], self.trjs[plotting_trj][c][:, 1]*10, label=c)
        plt.title(plotting_trj)
        plt.legend(loc=1, frameon=False)
        plt.xlabel('Residue')
        plt.ylabel(r'RMSF ($\AA$)')
        

rmsf = rmsf_Comparer()
'''
trjInfo = {'HL_apo': 1,
           'HL_lig2201': 1,
           'PL_apo': 4,
           'PL_lig2201': 3}

filePths = []
for t in trjInfo:
    trj_num = trjInfo[t]
    for i in range(trj_num):
        filePths.append(os.sep.join([WORKDIR, 'rmsf_'+t+'_%d.xvg' % (i+1)]))
trjNames = []

for n in filePths:
    if 'HL' in n:
        d = 'new'
    elif 'PL' in n:
        d = 'old'
    t = n.split('_')[2]
    i = n.split('_')[-1][0]
    trjNames.append('_'.join([d, t, i]))

filePths = [os.sep.join([WORKDIR, 'rmsf', 'rmsf_apo.xvg']),
            os.sep.join([WORKDIR, 'rmsf', 'rmsf_holo.xvg']),
            os.sep.join([WORKDIR, 'PLIU_Traj', 'rmsf', 'rmsf_apo.xvg']),
            os.sep.join([WORKDIR, 'PLIU_Traj', 'rmsf', 'rmsf_holo.xvg'])]
'''
filePths = [os.sep.join([WORKDIR, 'holo', 'rmsf_1.xvg']),
            os.sep.join([WORKDIR, 'holo', 'rmsf_2.xvg']),
            os.sep.join([WORKDIR, 'holo', 'rmsf_3.xvg'])]
trjNames = ['holo 1', 'holo 2', 'holo 3']
for n, f in zip(trjNames, filePths):
    rmsf.loadTrj(f, n)

rmsf.compareTrjsRMSF(trjNames, chain='A')

genName = lambda u, t, i: '_'.join([u, t, str(i)])
rmsf.compareTrjsRMSF([genName('old', 'apo', 1), genName('old', 'apo', 2), genName('old', 'apo', 3), genName('old', 'apo', 4), genName('new', 'apo', 1)],chain='B')
rmsf