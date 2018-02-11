source $phgdh/code/environ_set.sh
condition='prot'
set_parms

# atomIndex:
# Prot Prot_A_CA 27 Prot_B_CA 28
# ProtLig Prot_A_CA 29 Prot_B_CA 30
# ProtMlt Prot_A_CA 29 Prot_B_CA 30
# ProtMltLig Prot_A_CA 31 Prot_B_CA 32
atomIndex=28

covar=true
extreStruct=false
proj=true

for i in 1 2 3 1-3
do
cd ${i}

if $covar
then
gmx covar -s $struct -f ${TrjDIR}/${i}_dropF50.xtc -last 100 -n $index << EOF
$atomIndex
$atomIndex
EOF
fi

if $extreStruct
then
gmx anaeig -s $struct -f ${TrjDIR}/${i}_dropF50.xtc -extr extreme1.pdb -first 1 -last 1 -n $index << EOF
$atomIndex
$atomIndex
EOF
gmx anaeig -s $struct -f ${TrjDIR}/${i}_dropF50.xtc -extr extreme2.pdb -first 2 -last 2 -n $index<< EOF
$atomIndex
$atomIndex
EOF
fi

if $proj
then
gmx anaeig -s $struct -f ${TrjDIR}/${i}_dropF50.xtc -proj pc1.xvg -first 1 -last 1 -n $index << EOF
$atomIndex
$atomIndex
EOF
gmx anaeig -s $struct -f ${TrjDIR}/${i}_dropF50.xtc -proj pc2.xvg -first 2 -last 2 -n $index<< EOF
$atomIndex
$atomIndex
EOF
gmx anaeig -s $struct -f ${TrjDIR}/${i}_dropF50.xtc -proj -2d pc1-2.xvg -first 1 -last 2 -n $index<< EOF
$atomIndex
$atomIndex
EOF
fi

cd ..
done
