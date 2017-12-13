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

proc calDist {com1 com2} {
     set sum 0
     foreach i $com1 j $com2 {
             set diff [expr $i-$j]
             set diff [expr $diff * $diff]
             set sum [expr $sum + $diff]
     }
     set dist [expr sqrt($sum)]
}

set ll "chain $chain and resid 57 and sidechain"
set lm "chain $chain and resid 80 and sidechain"
set lr "chain $chain and resid 99 and sidechain"
set ul "chain $chain and resid 264 and sidechain"
set um "chain $chain and resid 235 and sidechain"
set ur "chain $chain and resid 154 and sidechain"

set out [open dist.dat w]

# write out each frame as a pdb file with protein and the NME terminus.
set nf [molinfo 0 get numframes]
for {set i 1} {$i < $nf} {incr i} {

     set g1 [atomselect 0 $ll frame $i]
     set g2 [atomselect 0 $ul frame $i]
     set g3 [atomselect 0 $lm frame $i]
     set g4 [atomselect 0 $um frame $i]
     set g5 [atomselect 0 $lr frame $i]
     set g6 [atomselect 0 $ur frame $i]
     set com1 [measure center $g1]
     set com2 [measure center $g2]
     set com3 [measure center $g3]
     set com4 [measure center $g4]
     set com5 [measure center $g5]
     set com6 [measure center $g6]
     set dist1 [calDist $com1 $com2]
     set dist2 [calDist $com3 $com4]
     set dist3 [calDist $com5 $com6]
     puts $out "$dist1 $dist2 $dist3"
}

exit
