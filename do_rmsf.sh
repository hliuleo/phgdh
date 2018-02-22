source $phgdh/code/environ_set.sh
condition='prot'
set_parms

# atom index
# protMlt chain A 27 chian B 28

atomIndex=27


for i in 1 2 3 1-3
do
  if [ ! -d "${i}" ]
  then
      mkdir ${i}
  fi
cd ${i}
gmx rmsf -s $struct -f ${TrjDIR}/${i}_dropF50.xtc -n $index -res << EOF
$atomIndex
EOF
cd ..
done
