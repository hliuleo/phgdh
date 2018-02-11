source $phgdh/code/environ_set.sh
condition='prot'
set_parms

for i in 1 2 3
do
  cd ${i}
  vmd -dispdev text -e $CodeHOME/com_plane.tcl -args $struct ${TrjDIR}/${i}_dropF50.xtc A
  python $CodeHOME/get_angles.py com.dat
  awk '{print $1}' coords.dat > o.dat
  awk '{print $2}' coords.dat > t.dat
  more coords.dat >> ../coords.dat
  more o.dat >> ../o.dat
  more t.dat >> ../t.dat
  cd ..
done

mkdir 1-3
mv coords.dat o.dat t.dat 1-3/
