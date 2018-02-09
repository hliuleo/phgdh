#!/home/hliu/anaconda2/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 09:54:59 2017

@author: hliu
"""

import re
import argparse
from researchcode.pdb.mol import *


PDBFMT = ['{}',
          '{:>7}',
          '{:>4}',
          '{:>5} ',
          '{:}',
          '{:>4}',
          '{:>12.3f}',
          '{:>8.3f}',
          '{:>8.3f}',
          '{:>6.2f}',
          '{:>6.2f}',
          '{:>12}']


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='file')
    parser.add_argument('--fix',
                        dest='fix',
                        action='store_true')
    parser.add_argument('--split',
                        dest='split',
                        action='store_true')
    args = parser.parse_args()

    f = open(args.file, 'r')
    f_str = ''.join(f.readlines())
    
    splitter = 'MODEL        1\n'
    
    models = re.split(splitter, f_str)

    '''Have to fix here to make it more general'''
    def split_atomInfo(model):
        lines = model.split('\n')
        atoms = []
        for idx, l in enumerate(lines):
            if l.startswith('ATOM'):
                atoms.append(l)
            else:
                other_idx = idx
                break
        atom_str = '\n'.join(atoms)
        other_str = '\n'.join(lines[other_idx:])
        return atom_str, other_str

    def fix_chainID(model):
    
        
        atom_str, other_str = split_atomInfo(model)
        atom_df = molStr2DF(atom_str)
        seq = fetchSeq(atom_df)
        index = sorted(seq.index, key=lambda x: seq.ix[x, 'ResID'])
        period_idx = (index[0], index[2])
        atom_df.ix[period_idx[0]:period_idx[1], 4] = 'A'
        atom_df.ix[period_idx[1]:, 4] = 'B'
        FMT = ''.join(PDBFMT[:len(atom_df.columns)])
        fixed_atom_str = '\n'.join([FMT.format(*atom_df.ix[i]) for i in atom_df.index])
        return '\n'.join([fixed_atom_str,other_str])

    if args.fix:

        new_m = [models[0]]
        for m in models[1:]:
            new_m.append(fix_chainID(m))
        new_str = splitter.join(new_m)
        
        f_out = open(args.file, 'w')
        f_out.write(new_str)

    if args.split:
        prefix = args.file.split('.')[0]
        for idx, m in enumerate(models[1:]):
            atom_str, other_str = split_atomInfo(m)
            f_out = open(prefix+'_{}.pdb'.format(idx+1), 'w')
            f_out.write(atom_str)
            f.close()
