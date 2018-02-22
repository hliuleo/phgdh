source $phgdh/code/environ_set.sh
condition='prot'
set_parms

for i in 1 2 3
do
  if [ ! -d "${i}" ]
  then
      mkdir ${i}
  fi
  cd ${i}
  vmd -dispdev text -e $CodeHOME/hole_rg.tcl -args $struct ${TrjDIR}/${i}_dropF50.xtc A
  more hole_rg.dat >> ../hole_rg.dat
  cd ..
done

mkdir 1-3
mv hole_rg.dat 1-3/
