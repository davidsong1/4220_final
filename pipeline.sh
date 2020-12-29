#!/bin/bash

settings=$1
# checks if the job directory is provided
if [[ ! -z $2 ]]
then
    #creates directory from argument
    mkdir -p $2
    cd $2
else
    #creates a tmp directory
    mkdir -p tmp
    cd tmp
fi
#parse settings for get_seq.sh
echo "Running get_seq.sh"
../parse_settings.sh ../$settings get_seq.sh

