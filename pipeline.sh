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
#echo "Running get_seq.sh"
#../parse_settings.sh ../$settings get_seq.sh

#parse settings for make_align.sh
#echo "Running make_align.sh"
#../parse_settings.sh ../$settings make_align.sh

#parse settings for make_phylo.sh
#echo "Running make_phylo.sh"
#../parse_settings.sh ../$settings make_phylo.sh

#parse settings for make_mol_stats.py
#echo "Running make_mol_stats.py"
#../parse_settings.sh ../$settings make_mol_stats.py

#parse settings for make_dnds.py
#echo "Running make_dnds.py"
#../parse_settings.sh ../$settings make_dnds.py

#parse settings for get_hydrophobic_image.py
#echo "Running get_hydrophobic_image.py"
#../parse_settings.sh ../$settings get_hydrophobic_image.py

#parse settings for generate_hydrophilicity_curve.py
#echo "Running generate_hydrophilicity_curve.py"
#../parse_settings.sh ../$settings generate_hydrophilicity_curve.py

#parse settings for make_results.py
echo "Runngin make_results.py"
../parse_settings.sh ../$settings make_results.py
