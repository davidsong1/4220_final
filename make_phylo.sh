#!/bin/bash

#get the alignment file
ALIGN=$1
#get the phylo tool
TOOL=$2
#get the basename of alignment file
BASE=`basename $ALIGN .fasta`
#get the prefix for phylo
PRE=$BASE.phylo_$TOOL
#create directories based on the tool
mkdir -p phylo
mkdir -p phylo/$BASE

#function to generate display and log file
generate_log_and_display() {
    #generate the nw_display
    nw_display ./phylo/$BASE/$TOOL/$PRE.tre > ./phylo/$BASE/$TOOL/$PRE.nw_display.txt
    #get name of file with phylogenetic estimate
    echo "Name of file: ./phylo/$BASE/$TOOL/$PRE.tre" > ./phylo/$BASE/$TOOL/$PRE.log
    #get the command used
    echo "Command: $1" >> ./phylo/$BASE/$TOOL/$PRE.log
    #get when phylogeny was created
    echo "Date: $(date)" >> ./phylo/$BASE/$TOOL/$PRE.log
}

#checks if the tool is fasttree
if [[ $TOOL == fasttree ]]
then
    mkdir -p phylo/$BASE/$TOOL
    #checks if the option is -gtr
    if [[ $3 == -gtr ]]
    then
        #runs fasttree with gtr
        fasttree -nt -noml -gtr $ALIGN > ./phylo/$BASE/$TOOL/$PRE.tre
        CMD="fasttree -nt -noml -gtr $ALIGN > ./phylo/$BASE/$TOOL/$PRE.tre"
        generate_log_and_display "$CMD"
    else
        #runs default fasstree
        fasttree -nt -noml $ALIGN > ./phylo/$BASE/$TOOL/$PRE.tre
        CMD="fasttree -nt -noml $ALIGN > ./phylo/$BASE/$TOOL/$PRE.tre"
        generate_log_and_display "$CMD"
    fi
fi
