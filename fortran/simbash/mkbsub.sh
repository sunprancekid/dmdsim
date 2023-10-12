#!/bin/bash

# AUTHOR: Matthew Dorsey
# DATE: 2022-07-21

# PURPOSE
# program that generates henry2 executable files

# ARGUMENTS
# 1 : directory to store bsub file in
# 2 : unqiue simulation ID (needs to be 10 characters or less in order to be displayed correctly in bjobs output)
# 3 : title of main molecular simulation program to compile and execut
# 4 (OPTIONAL) : queue for job submission. if not specified, serial (max time is 4 days or 5760 minutes).
# 5 (OPTIONAL) : wall clock time for job. if not specified, 1 day (1440 minutes). keep in mind queue specific clock limits.

# CONSTANTS
BSUB=$1"test.bsub" # title of bsub file containing executable instructions
SIMID=$2 # unique id used to identify simulation during in HPC jobs
MOLSIM=$3 # title of main molecular simulation program that each job should compile and execute
QUEUE=""
if [ $# -ge 4 ]
then
    QUEUE="$4"
else
    QUEUE="serial"
fi
WALLCLOCK=""
if [ $# -ge 5 ]
then
    WALLCLOCK="$5"
else
    WALLCLOCK="1440"
fi

# BSUB PARAMETERS
echo "#!/bin/tcsh" > $BSUB # headbang tag for HENRY2 cluster and executables, overwrites existing bsub files
echo "#BSUB -n 1" >> $BSUB # number of cores
echo "#BSUB -W $WALLCLOCK" >> $BSUB # max wall clock time of molecular simulation
echo "#BSUB -o out.$SIMID.%J" >> $BSUB # name of standard out file from HPC job (%J corresponds to HENRY2 job id)
echo "#BSUB -e err.$SIMID.%J" >> $BSUB # name of standard err file from HPC job (%J corresponds to HENRY2 job id)
echo "#BSUB -q $QUEUE" >> $BSUB # queue to submit job to on HENRY2 cluster
echo "#BSUB -J $SIMID" >> $BSUB # unique id associated with job submissiob
echo "" >> $BSUB

# BSUB SCRIPT
echo "module load PrgEnv-intel" >> $BSUB
echo "ifort -no-wrap-margin -O3 -fp-model=precise $MOLSIM" >> $BSUB
echo "./a.out" >> $BSUB
