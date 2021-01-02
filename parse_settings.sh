#!/bin/bash

settings=$1
step=$2
#gets the extension
extension="${step##*.}"
#iterate through the settings file
cat $settings | while read line
do
    #get the pipeline step for each line
    pipeline=`echo $line | cut -d ',' -f 1`
    #check if the pipeline step matched the input step
    if [[ $pipeline == $step ]]
    then
        IFS=$'\n'
        #initializes command array
        command=()
        #gets variable between = and ;
        for var in `echo $line | grep -oP '(?<==).*?(?=;)'`
        do
            #concatenates variables into the command
            command+=("$var")
        done
        #runs the script with the arguments
        if [[ $extension == "sh" ]]
        then
            . ../$step "${command[@]}"
            echo "./$step ${command[@]}" >> commands.txt
        elif [[ $extension == "py" ]]
        then
            python ../$step "${command[@]}"
            echo "./$step ${command[@]}" >> commands.txt
        fi
    fi
done




