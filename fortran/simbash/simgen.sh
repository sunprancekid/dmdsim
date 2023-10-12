#!/bin/bash

# AUTHOR: Matthew Dorsey
# DATE: 2022-07-21

# This program generates nested directories used for molecular simulation jobs on HPC systems.
# The directories that are generated are hirearchical according to simulation parameters that are used during execution.
# This script utilizes a series of other scripts that are stored in the simbin directory.
# Each final directory contains a file which stores simulation parameters.

# ARGUMENTS
# 1 : location of / path to bin containing bash scripts (simbin) and molecular simulation software (molsim)
# 2 : location of home directory to generate and execute molecular simulations from

# PATHWAY AND FILES
SIMBIN=$1"simbin/" # location of executables and molecular simulation software
HOME="~/share/hall2/madorse2/" # location where directory generation and simulation exection occurs
MOD="hardspheremodule.mod" # title of fortran module to use
# SIM=

# SIMULATION SETTINGS
EVENT1=700000000 # number of events before equilibrium period
EVENT2=300000000 # number of events during equilibrium period
TEMPINIT=1.5 # initial simulation temperature

# establish DOE parameters
MODULE=("HS") # defines the MD module to use
CELL=("05" "10" "15" "20" "25" "30" "35" "40") # array containing cell sizes used in DOE
ETA=("05" "10" "15") # array containing volume / area fractions used in DOE
TEMP=("1.5") # array containing temperatures used in DOE

# loop through arrays, generate directories and simid
DIRECTORY="./"
SIMID=""

for D1 in ${MODULE[@]}
do

SIMID=$SIMID$D1

    for D2 in ${CELL[@]}
    do

    SIMID=$SIMID"c"$D2

        for D3 in ${ETA[@]}
        do

        SIMID=$SIMID"e"$D3

            for D4 in ${TEMP[@]}
            do

                # generate directories
                DIRECTORY="./"$D1"/c"$D2"/e"$D3"/"$D4"/"
                mkdir -p $DIRECTORY
                # bash "$SIMBIN"mkdirectories.sh $D1 "c"$D2 "e"$D3 $D4
                
                # generate the simulation file
                bash "$SIMBIN"mksimfile.sh $DIRECTORY$D1".sim" 0 $SIMID $EVENT1 $EVENT2 $D2 $D3 $D4
                
                # HENRY2: generate bsub file and execute
                bash "$SIMBIN"mkbsub.sh $DIRECTORY $SIMID "molecularsimulation.f90" "mrsec" "5760"

            done
        done
    done
done

exit 0
