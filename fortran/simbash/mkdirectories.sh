#!/bin/bash
# loop through aguments passed to executable
DIRECTORY="./"
for i in "$@"
do
    # add the argument to the directory
    DIRECTORY=$DIRECTORY"/"$i
    # establish that the directory has been established
    if test ! -d $DIRECTORY
    then
        mkdir $DIRECTORY
    fi
done
exit 0
