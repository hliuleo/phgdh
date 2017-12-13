# Use to extract cluster snapshot and save to pdb

## check if there are correct number of arguments
if { [llength $argv] != 2} {
puts "There are $argc arguments to this script"
puts "The name of this script is $argv0"
if { [llength $argv] > 0} {puts "The other arguments are: $argv" }

exit
}

set SysPdb [lindex $argv 0]
set TrjDCD [lindex $argv 1]

mol new $SysPdb
animate read xtc $TrjDCD waitfor all 0

set upperA "chain A and resid 101 to 285 and name CA"
set upperB "chain B and resid 101 to 285 and name CA"
set lowerA "chain A and resid 6 to 96 288 to 304 and name CA"
set lowerB "chain B and resid 6 to 96 288 to 304 and name CA"
set lineA "chain A and resid 241 to 250 and name CA"
set lineB "chain B and resid 241 to 250 and name CA"
set tpA "chain A and resid 84 to 90 and name CA"
set tpB "chain B and resid 84 to 90 and name CA"
set opA "chain A and resid 218 to 223 and name CA"
set opB "chain B and resid 218 to 223 and name CA"

set out1 [open A_com.dat w]
set out2 [open B_com.dat w]


# write out each frame as a pdb file with protein and the NME terminus.
set nf [molinfo 0 get numframes]
for {set i 1} {$i < $nf} {incr i} {

     set g1A [atomselect 0 $upperA frame $i]
     set g2A [atomselect 0 $lowerA frame $i]
     set g3A [atomselect 0 $lineA frame $i]
     set g4A [atomselect 0 $tpA frame $i]
     set g5A [atomselect 0 $opA frame $i]
     set com1A [measure center $g1A]
     set com2A [measure center $g2A]
     set com3A [measure center $g3A]
     set com4A [measure center $g4A]
     set com5A [measure center $g5A]
     puts $out1 "$com1A $com2A $com3A $com4A $com5A"
     set g1B [atomselect 0 $upperB frame $i]
     set g2B [atomselect 0 $lowerB frame $i]
     set g3B [atomselect 0 $lineB frame $i]
     set g4B [atomselect 0 $tpB frame $i]
     set g5B [atomselect 0 $opB frame $i]
     set com1B [measure center $g1B]
     set com2B [measure center $g2B]
     set com3B [measure center $g3B]
     set com4B [measure center $g4B]
     set com5B [measure center $g5B]
     puts $out2 "$com1B $com2B $com3B $com4B $com5B"
}

exit
