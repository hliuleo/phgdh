# Use to extract cluster snapshot and save to pdb

## check if there are correct number of arguments
if { [llength $argv] != 3} {
puts "There are $argc arguments to this script"
puts "The name of this script is $argv0"
if { [llength $argv] > 0} {puts "The other arguments are: $argv" }

exit
}

set SysPdb [lindex $argv 0]
set TrjDCD [lindex $argv 1]
set chain  [lindex $argv 2]

mol new $SysPdb
animate read xtc $TrjDCD waitfor all 0

set l { 9214 9219 }

if { $chain == "A" } {
    set l { 9131 9136 }
}

set m "chain $chain and resname MLT"
set a1 "chain $chain and resname MLT and name O2"
set a2 "chain $chain and resname MLT and name O4"

set out [open dist.dat w]

# write out each frame as a pdb file with protein and the NME terminus.
set nf [molinfo 0 get numframes]

for {set i 1} {$i < $nf} {incr i} {

     set g1 [atomselect 0 $m frame $i]
     set g2 [atomselect 0 $a1 frame $i]
     set g3 [atomselect 0 $a2 frame $i]
     set rg [measure rgyr $g1]
#     set i1 [$g2 get index]
#     set i2 [$g3 get index]
#     puts $i1 $i2
#     set l  {$i2 $i2}
     set dist [measure bond $l]
     puts $out "$rg $dist"
}

exit
