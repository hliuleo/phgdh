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

set chainS "A"
if { $chain eq "A" } {
   set chainS "B"
}

set s1 "chain $chain and resid 53 to 57"
set s2 "chain $chain and resid 259 to 264 282 285 286"
set s3 "chain $chain and resid 77 to 80"
set s4 "chain $chain and resid 233 to 237"
set s5 "chain $chain and resid 12 to 15 33 34 288"
set s6 "chain $chainS and resid 134 135 138"
set s7 "chain $chain and resid 53 to 57 259 to 264 282 285 286"
set s8 "chain $chain and resid 77 to 80 233 to 237"

set out [open hole_contact.dat w]

# write out each frame as a pdb file with protein and the NME terminus.
set nf [molinfo 0 get numframes]
for {set i 1} {$i < $nf} {incr i} {

     set g1 [atomselect 0 $s1 frame $i]
     set g2 [atomselect 0 $s2 frame $i]
     set g3 [atomselect 0 $s3 frame $i]
     set g4 [atomselect 0 $s4 frame $i]
     set g5 [atomselect 0 $s5 frame $i]
     set g6 [atomselect 0 $s6 frame $i]
     set g7 [atomselect 0 $s7 frame $i]
     set g8 [atomselect 0 $s8 frame $i]
     set c1 [llength [lindex [measure contacts 6 $g1 $g2] 0]]
     set c2 [llength [lindex [measure contacts 6 $g3 $g4] 0]]
     set c3 [llength [lindex [measure contacts 6 $g5 $g6] 0]]
     set c4 [llength [lindex [measure contacts 6 $g7 $g8] 0]]
# to be consistent with the order of Rg's results
     puts $out "$c4 $c1 $c3 $c2" 
}

exit
