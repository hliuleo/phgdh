condition='prot'
source /home/hliu/Projects/PHGDH/code/environ_set.sh

for i in 1 2 3
do
  cd ${i}
  vmd -dispdev text -e $CodeHOME/hbond_dist.tcl -args $struct ${TrjDIR}/${i}_dropF50.xtc A
  more dist.dat >> ../dist.dat
  cd ..
done

mkdir 1-3
mv dist 1-3/
