export PATH=$PATH:/home/hliu/Projects/PHGDH/code

DATA='/home/hliu/Projects/PHGDH/data'
CodeHOME='/home/hliu/Projects/PHGDH/code'
SYS_SETUP='/home/hliu/Projects/PHGDH/system_prep/SYS'

TrjHOME='/home/hliu/Projects/PHGDH/data/traj/'
ApoTrj=$TrjHOME/'HLIU_Prot_CoF'
ApoLigTrj=$TrjHOME/'HLIU_Prot_CoF_LIG2201'
ApoMltTrj=$TrjHOME'HLIU_Prot_CoF_SubS_dimer'
ApoMltLigTrj=$TrjHOME'HLIU_Prot_CoF_SubS_dimer_LIG2201'

StructHOME='/home/hliu/Projects/PHGDH/data/struct/'
Prot=$StructHOME/'protein_NAD.pdb'
ProtLig=$StructHOME/'protein_NAD_LIG2201.pdb'
ProtMlt=$StructHOME/'protein_NAD_MLT2.pdb'
ProtMltLig=$StructHOME/'protein_NAD_MLT2_LIG2201.pdb'
Prot_index=$StructHOME/'protein_NAD.ndx'
ProtLig_index=$StructHOME/'protein_NAD_LIG2201.ndx'
ProtMlt_index=$StructHOME/'protein_NAD_MLT2.ndx'
ProtMltLig_index=$StructHOME/'protein_NAD_MLT2_LIG2201.ndx'

if [ $condition == 'protLig' ]
then
TrjDIR=$ApoLigTrj
struct=$ProtLig
index=$ProtLig_index
fi

if [ $condition == 'prot' ]
then
TrjDIR=$ApoTrj
struct=$Prot
index=$Prot_index
fi

if [ $condition == 'protMlt' ]
then
TrjDIR=$ApoMltTrj
struct=$ProtMlt
index=$ProtMlt_index
fi

if [ $condition == 'protMltLig' ]
then
TrjDIR=$ApoMltLigTrj
struct=$ProtMltLig
index=$ProtMltLig_index
fi

