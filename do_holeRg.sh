condition='prot'
source /home/hliu/Projects/PHGDH/code/environ_set.sh

for i in 1 2 3
do
  cd ${i}
  vmd -dispdev text -e $CodeHOME/hole_rg.tcl -args $struct ${TrjDIR}/${i}_dropF50.xtc A
  more hole_rg.dat >> ../hole_rg.dat
  cd ..
done

mkdir 1-3
mv hole_rg.dat 1-3/
