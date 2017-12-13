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

set h1 "chain $chain and resid 53 to 57 259 to 264 282 285 286 77 to 80 233 to 237"
set h2 "chain $chain and resid 53 to 57 259 to 264 282 285 286"
set h3 "chain $chain and resid 12 to 15 33 34 288 or (chain $chainS and resid 134 135 138)"
set h4 "chain $chain and resid 77 to 80 233 to 237"

set out [open hole_rg.dat w]

# write out each frame as a pdb file with protein and the NME terminus.
set nf [molinfo 0 get numframes]
for {set i 1} {$i < $nf} {incr i} {

     set g1 [atomselect 0 $h1 frame $i]
     set g2 [atomselect 0 $h2 frame $i]
     set g3 [atomselect 0 $h3 frame $i]
     set g4 [atomselect 0 $h4 frame $i]
     set rg1 [measure rgyr $g1]
     set rg2 [measure rgyr $g2]
     set rg3 [measure rgyr $g3]
     set rg4 [measure rgyr $g4]
     puts $out "$rg1 $rg2 $rg3 $rg4"
}

exit
