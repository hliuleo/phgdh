# -*- coding: utf-8 -*-
# Doc path of projetct p53

import os

PROJPATH = os.environ['phgdh']
CODEPATH = os.sep.join([PROJPATH, 'code'])
STRUCTPATH = os.sep.join([PROJPATH, 'data', 'struct'])
TRAJPATH = os.sep.join([PROJPATH, 'data', 'traj'])

and_str = ' and '
or_str = ' or '

def build_new_top(ways, parts):
    parts_str = [topology[p] for p in parts]
    if isinstance(ways, str):
        return '('+ways.join(parts_str)+')'
    elif isinstance(ways, list):
        strs = parts[0]
        for idx, p in enumerate(parts[1:]):
            strs = strs+ways[idx]+p
        return '('+strs+')'

topology = dict(chA = "segid A",
                chB = "segid B",
                protein = "protein",
                NAD = "resname NDP",
                LIG = "resname LIG",
                
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
                H9 = 'resid 192:194',
                H10 = 'resid 218:223',
                H11 = 'resid 241:250',
                
                H12 = 'resid 288:305',
                
                L1 = 'resid 53:59',
                L2 = 'resid 75:82',
                L3 = 'resid 174:178',
                L4 = 'resid 205:217',
                L5 = 'resid 233:240',
                L6 = 'resid 259:276',
                L7 = 'resid 130:139',
                
                btm = 'resid 7:95',
                top = 'resid 146:279',
                spine = '(resid 101:138 or resid 138:287)')

_btmBeta = ['E1', 'E2', 'E3', 'E4', 'E5']
_topBeta = ['E6', 'E7', 'E8', 'E9', 'E10', 'E11']
_btmHelix = ['H1', 'H2', 'H3', 'H4']
_topHelix = ['H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11']
_btmSecStruct = _btmBeta+_btmHelix
_topSecStruct = _topBeta+_topHelix

multi_segs = dict(btmBeta = build_new_top(or_str, _btmBeta),
                  topBeta = build_new_top(or_str, _topBeta),
                  btmHelix = build_new_top(or_str, _btmHelix),
                  topHelix = build_new_top(or_str, _topHelix),
                  btmSecStruct = build_new_top(or_str, _btmSecStruct),
                  topSecStruct = build_new_top(or_str, _topSecStruct)
                 )
                  
                                          
topology.update(multi_segs)



            