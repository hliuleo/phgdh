# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import numpy as np
import pandas as pd

homeDir = os.environ['phgdh']
workDir = homeDir + os.sep + 'data/system_prep/SYS'

os.chdir(workDir)

pdbFmt = '%-6s%5d  %-4s%s %s %d%12.3f%8.3f%8.3f%6.2f%6.2f           %-2s'
        

template = pd.read_csv('NAD_chainA.pdb', sep='\s+', header=None)
atom_name = template[2]
atom_type = template[11]

pdb = pd.read_csv('NAD_chainB.pdb', sep='\s+', header=None)
pdb.ix[pdb[3]=='NAD', 3] = 'NDP'
pdb[2] = atom_name
pdb[11] = atom_type

np.savetxt('test.pdb', pdb.values, fmt=pdbFmt)