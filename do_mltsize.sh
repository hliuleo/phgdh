condition='protMlt'
source /home/hliu/Projects/PHGDH/code/environ_set.sh

for i in 1 2 3 4
do
  cd ${i}
  vmd -dispdev text -e $CodeHOME/mlt_size.tcl -args $struct ${TrjDIR}/${i}.xtc A
  more dist.dat >> ../dist.dat
  cd ..
done

mkdir 1-4
mv dist 1-4/
