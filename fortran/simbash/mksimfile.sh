#!/bin/bash
# parse arguments to variables
FILE=$1 # file location and name
NUM=$2 # number of times the simulation has been executed, typically one
SIMID=$3 # unique id associated with simulation
EVENT1=$4 # number of events before equilibrium period
EVENT2=$5 # number of events during eqilibrium period
CELL=$6 # number of cells in one dimension
FRAC=$7 # volume of area fraction of simulation
TEMP=$8 # temperature of simulation

# write variables to simfile
printf " %s\n" $FILE
printf " %s\n %s\n %s\n %s\n %s\n %s\n %s\n" $NUM $SIMID $EVENT1 $EVENT2 $CELL $FRAC $TEMP > $FILE

exit 0
