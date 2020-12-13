#!/bin/bash
#get the alignment file
ALIGN=$1
#get the phylo tool
TOOL=$2
#get the basename of alignment file
BASE=`basename $ALIGN .fasta`
#get the prefix for phylo
PRE=$BASE.phylo_$TOOL
#get the path
path=$(dirname "$ALIGN")
#create a tmp file with no spaces
cat $ALIGN | tr -d " " > tmp
#replace the alignment with the tmp
mv tmp $ALIGN
rm tmp
#create directories based on the tool
mkdir -p ./$path/phylo
mkdir -p ./$path/phylo/$BASE

#function to generate display and log file
generate_log_and_display() {
    #generate the nw_display
    nw_display $path/phylo/$BASE/$TOOL/$PRE.tre > $path/phylo/$BASE/$TOOL/$PRE.nw_display.txt
    #get name of file with phylogenetic estimate
    echo "Name of file: $path/phylo/$BASE/$TOOL/$PRE.tre" > $path/phylo/$BASE/$TOOL/$PRE.log
    #get the command used
    echo "Command: $1" >> $path/phylo/$BASE/$TOOL/$PRE.log
    #get when phylogeny was created
    echo "Date: $(date)" >> $path/phylo/$BASE/$TOOL/$PRE.log
}

#checks if the tool is fasttree
if [[ $TOOL == fasttree ]]
then
    mkdir -p $path/phylo/$BASE/$TOOL
    #checks if the option is -gtr
    if [[ $3 == -gtr ]]
    then
        #runs fasttree with gtr
        fasttree -nt -noml -gtr $ALIGN > $path/phylo/$BASE/$TOOL/$PRE.tre
        CMD="fasttree -nt -noml -gtr $ALIGN > $path/phylo/$BASE/$TOOL/$PRE.tre"
        generate_log_and_display "$CMD"
    else
        #runs default fasstree
        fasttree -nt -noml $ALIGN > $path/phylo/$BASE/$TOOL/$PRE.tre
        CMD="fasttree -nt -noml $ALIGN > $path/phylo/$BASE/$TOOL/$PRE.tre"
        generate_log_and_display "$CMD"
    fi
elif [[ $TOOL == iqtree ]]
then
    mkdir -p $path/phylo/$BASE/$TOOL
    #create a tmp directory
    mkdir -p tmp
    #check if options GTR and JC are provided
    if [[ $3 == -m ]] && [[ $4 == GTR || $4 == JC ]]
    then
        #run iqtree and store output in tmp
        iqtree -s $ALIGN -m $4 -pre tmp/$PRE
        #move the treefile to the appropriate folder
        mv tmp/$PRE.treefile $path/phylo/$BASE/$TOOL/$PRE.tre
        #delete tmp folder
        rm -r tmp
        CMD="iqtree -s $ALIGN -m $4 -pre $path/phylo/$BASE/$TOOL/$PRE"
        generate_log_and_display "$CMD"
    else
        #run iqtree on default
        iqtree -s $ALIGN -pre tmp/$PRE
        #move the treefile to the appropriate folder
        mv tmp/$PRE.treefile $path/phylo/$BASE/$TOOL/$PRE.tre
        #delete tmp folder
        rm -r tmp
        CMD="iqtree -s $ALIGN -pre $path/phylo/$BASE/$TOOL/$PRE"
        generate_log_and_display "$CMD"
    fi
elif [[ $TOOL == mpboot ]]
then
    mkdir -p tmp
    mkdir -p $path/phylo/$BASE/$TOOL
    #run mpboot default settings
    mpboot -s $ALIGN -pre tmp/$PRE
    mv tmp/$PRE.treefile $path/phylo/$BASE/$TOOL/$PRE.tre
    #delete tmp folder
    rm -r tmp
    CMD="mpboot -s $ALIGN -pre $path/phylo/$BASE/$TOOL/$PRE"
    generate_log_and_display "$CMD"
fi
