source $phgdh/code/environ_set.sh
condition='protMlt'
set_parms

for i in 1 2 3 4
do
  if [ ! -d "${i}" ]
  then
      mkdir ${i}
  fi
  cd ${i}
  vmd -dispdev text -e $CodeHOME/mlt_size.tcl -args $struct ${TrjDIR}/${i}.xtc A
  more dist.dat >> ../dist.dat
  cd ..
done

mkdir 1-4
mv dist 1-4/
