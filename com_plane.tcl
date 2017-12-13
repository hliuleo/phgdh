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

set upper "chain $chain and resid 101 to 285 and name CA"
set lower "chain $chain and resid 6 to 96 288 to 304 and name CA"
set line "chain $chain and resid 241 to 250 and name CA"
set tp "chain $chain and resid 84 to 90 and name CA"
set op "chain $chain and resid 218 to 223 and name CA"

set out [open com.dat w]

# write out each frame as a pdb file with protein and the NME terminus.
set nf [molinfo 0 get numframes]
for {set i 1} {$i < $nf} {incr i} {

     set g1 [atomselect 0 $upper frame $i]
     set g2 [atomselect 0 $lower frame $i]
     set g3 [atomselect 0 $line frame $i]
     set g4 [atomselect 0 $tp frame $i]
     set g5 [atomselect 0 $op frame $i]
     set com1 [measure center $g1]
     set com2 [measure center $g2]
     set com3 [measure center $g3]
     set com4 [measure center $g4]
     set com5 [measure center $g5]
     puts $out "$com1 $com2 $com3 $com4 $com5"
}

exit
