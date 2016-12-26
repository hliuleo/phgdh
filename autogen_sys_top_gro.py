import re
import os
import itertools
chainID = ['protein', 'MLT_A', 'MLT_B', 'NDP_A', 'NDP_B', 'LIG2101_A', 'LIG2201_A']
moleculeType = ['protein', 'MLT', 'NDP', 'LIG2101', 'LIG2201']
packFiles = [['topol_Protein_chain_A.itp',
              'topol_Protein_chain_B.itp',
              'posre_Protein_chain_A.itp',
              'posre_Protein_chain_B.itp'],
             ['MLT.atp', 'MLT.itp', 'posre_MLT.itp'],
             ['NDP.atp', 'NDP.itp', 'posre_NDP.itp'],
             ['LIG2101.atp', 'LIG2101.itp', 'posre_LIG2101.itp'],
             ['LIG2201.atp', 'LIG2201.itp', 'posre_LIG2201.itp']
            ]

structFiles = ['protein_processed.gro',
               'MLT_chainA.gro',
               'MLT_chainB.gro',
               'NDP_chainA.gro',
               'NDP_chainB.gro',
               'LIG2101_chainA.gro',
               'LIG2201_chainA.gro']

resFiles = ['posre_MLT.inp',
            'posre_NDP.inp',
            'posre_LIG2101.inp',
            'posre_LIG2201.inp']

ffFiles = ['MLT.atp',
           'NDP.atp',
           'LIG2101.atp',
           'LIG2201.atp']

topFiles = ['MLT.itp',
            'NDP.itp',
            'LIG2101.itp',
            'LIG2201.itp']

structs = {x:y for x,y in zip(chainID, structFiles)}
tops = {x:y for x,y in zip(moleculeType[1:], topFiles)}
ffs = {x:y for x,y in zip(moleculeType[1:], ffFiles)}
restraints = {x:y for x,y in zip(moleculeType[1:], resFiles)}

class molecule:
    
    def __init__(self, chainID):
        self.name = chainID
        self.type = self.getMolType()
        self.loadMolInfo()
    
    def getMolType(self):
        for t in moleculeType:
            if t in self.name:
                return t

    def loadMolInfo(self):
        content = open(structs[self.name], 'r').readlines()
        self.header = content[0]
        self.atom_num = int(content[1].strip())
        self.size = content[-1]
        self.coor = ''.join(content[2:-1])

class topFile:

    def __init__(self, topPath):
        self.segNames = ['header', 'ff', 'top', 'restraint', 'chainTop']
        self.splitTop(topPath)

    def splitTop(self, topPath):       
        topTemplate = open('system_topol_template.top', 'r').read()
        splitter = re.compile('; Include')
        contents = splitter.split(topTemplate)
        for idx, n in enumerate(self.segNames):
            if n != 'header':
                c = contents[idx]
                c = '; Include'+c
                setattr(self, n, c[:-1])
            else:
                setattr(self, n, contents[idx][:-1])


def readGros(chainTop):
    gros = {}
    for t in chainTop:
        gros[t] = molecule(t)
    mols = set([gros[k].type for k in gros.keys()])
    gros['molTypes'] = mols
    ligNum, subNum, cofNum = (0, 0, 0)
    for m in mols:
        if 'LIG' in m:
            gros['ligName'] = m
            ligNum += 1
        elif 'MLT' in m:
            subNum += 1
        elif 'NDP' in m:
            cofNum += 1
    gros['ligNum'] = ligNum
    gros['cofNum'] = cofNum
    gros['subNum'] = subNum
    return gros

def buildSysGro(chainTop, gros):
    atomNum_fmt =  ' %d\n'
    header = gros['protein'].header
    coor = reduce(lambda x,y: x+y, [getattr(gros[t], 'coor') for t in chainTop])
    atomNum = atomNum_fmt % (sum([getattr(gros[t], 'atom_num') for t in chainTop]))
    size = gros['protein'].size
    sys_gro = header+atomNum+coor+size
    return sys_gro
    

def buildSysTop(chainTop, gros):
    top = topFile('system_topol_template.top')
    revised = {'top':'', 'ff':'', 'chainTop':''}
    file_fmt = '#include "%s"\n'
    sys_fmt  = '%-20s1\n'

    for t in gros['molTypes']:
        if t != 'protein':
            revised['ff'] += file_fmt % ffs[t]
            revised['top'] += file_fmt % tops[t]
    for t in chainTop:
        if t != 'protein':
            if 'LIG' in t:
                revised['chainTop'] += sys_fmt % 'LIG'
            else:
                revised['chainTop'] += sys_fmt % gros[t].type

    for r in revised:
        old = getattr(top, r)
        setattr(top, r, old+revised[r])
    
    newTop = '\n'.join([getattr(top, i) for i in top.segNames])
    return newTop


def packFile(gros, destDir):
    if os.path.exists(destDir):
        os.system('rm -rf %s' % destDir)
    os.mkdir(destDir)
    molFiles = {x:y for x,y in zip(moleculeType, packFiles)}
    forPacking = [molFiles[mol] for mol in gros['molTypes']]
    forPacking = list(itertools.chain.from_iterable(forPacking))
    for i in forPacking:
        os.system('cp %s %s' % (i, destDir))
    os.system('mv %s %s %s' % ('system.gro', 'topol.top', destDir))
    os.system('cp %s %s %s' % ('ions.mdp', 'system_setup.sh', destDir))
    if gros['subNum'] == 0 and gros['ligNum'] > 0:
        os.system('cp noMLT/%s %s' % (ffs[gros['ligName']], destDir))

def runGromax(destDir):
    os.chdir(destDir)
    os.system('bash system_setup.sh')


def prepSys(chainTop, destDir):
    groFile = open('system.gro', 'w')
    topFile = open('topol.top', 'w')
    
    gros = readGros(chainTop)
    sysGro = buildSysGro(chainTop, gros)
    sysTop = buildSysTop(chainTop, gros)
    groFile.write(sysGro)
    topFile.write(sysTop)
    packFile(gros, destDir)
    # runGromax(destDir)
     

chainTops = [['protein', 'MLT_A', 'NDP_A', 'MLT_B', 'NDP_B', 'LIG2101_A'],
             #['protein', 'MLT_A', 'NDP_A', 'LIG2101_A'],
             ['protein', 'MLT_A', 'NDP_A', 'MLT_B', 'NDP_B', 'LIG2201_A'],
             #['protein', 'MLT_A', 'NDP_A', 'LIG2201_A'],
             ['protein', 'NDP_A', 'NDP_B', 'LIG2101_A'],
             #['protein', 'NDP_A', 'LIG2101_A'],
             ['protein', 'NDP_A', 'NDP_B', 'LIG2201_A'],
             #['protein', 'NDP_A', 'LIG2201_A'],
             ['protein', 'MLT_A', 'NDP_A', 'MLT_B', 'NDP_B'],
             #['protein', 'MLT_A', 'NDP_A'],
             ['protein', 'NDP_A', 'NDP_B'],
             #['protein', 'NDP_A']
            ]
            

folders = ['Prot_CoF_SubS_dimer_LIG2101',
           #'Prot_CoF_SubS_monomer_LIG2101',
           'Prot_CoF_SubS_dimer_LIG2201',
           #'Prot_CoF_SubS_monomer_LIG2201',
           'Prot_CoF_dimer_LIG2101',
           #'Prot_CoF_monomer_LIG2101',
           'Prot_CoF_dimer_LIG2201',
           #'Prot_CoF_monomer_LIG2201',
           'Prot_CoF_SubS_dimer',
           #'Prot_CoF_SubS_monomer',
           'Prot_CoF_dimer',
           #'Prot_CoF_monomer',           
           ]



