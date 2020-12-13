#!/bin/bash

#get the working directory
CURRENTDIR=$(pwd)
#get the input tool
TOOL=$2
#change to input directory
cd $1
mkdir -p alignment
#concatenate all fasta files in the sequence directory
cat ./*.fasta > ./alignment/seq.fasta
cd ./alignment
#gets the name of the directory
NAME=$(echo $1 | rev | cut -d "/" -f2)

#get the filename
FILENAME="$NAME.align_$TOOL.fasta"

#function to generate log files
generate_log() {
    SEQUENCE=`cat $FILENAME | grep -o ">" | wc -l`
    SITE=`cat $FILENAME | grep -v ">" | wc -c`
    echo "Name of alignment file: $FILENAME" > $NAME.align_$TOOL.log
    echo "Command: $1" >> $NAME.align_$TOOL.log
    echo "Date: $(date)" >> $NAME.align_$TOOL.log
    echo "Number of sequences: $SEQUENCE" >> $NAME.align_$TOOL.log
    echo "Number of sites: $SITE" >> $NAME.align_$TOOL.log
}

#check if the align tool is muscle
if [[ $2 == muscle ]]
then
  #checks if gapopen flag is provided
  if [[ $3 == "-gapopen" ]]
  then
    #run the muscle alignment and all the generate_log function
    muscle -in seq.fasta -out $FILENAME -gapopen $4
    CMD="muscle -in seq.fasta -out $FILENAME -gapopen $4"
    generate_log "$CMD"
  else
    #run default muscle
    muscle -in seq.fasta -out $FILENAME
    CMD="muscle -in seq.fasta -out $FILENAME"
    generate_log "$CMD"
  fi

#check if the align tool is mafft
elif [[ $2 == mafft ]]
then
  #defaul parameters for op and ep if no options are provided
  op=1.53
  ep=0.123
  #updates op and ep if there are values given
  if [[ $3 == "--op" ]]
  then
    op=$4
  elif [[ $5 == "--op" ]]
  then
    op=$6
  fi
  if [[ $3 == "--ep" ]]
  then
    ep=$4
  elif [[ $5 == "--ep" ]]
  then
    ep=$6
  fi
  #run the mafft alignment on input
  mafft --op $op --ep $ep seq.fasta > $FILENAME
  CMD="mafft --op $op --ep $ep seq.fasta > $FILENAME"
  generate_log "$CMD"

#check if the align tool is prank
elif [[ $2 == prank ]]
then
  #default parameters for gaprate and gapext
  gaprate=0.025
  gapext=0.75
  #updates gaprate and gapext if there are values given
  if [[ $3 == "-gaprate" ]]
  then
    gaprate=$4
  elif [[ $5 == "-gaprate" ]]
  then
    gaprate=$6
  fi
  if [[ $3 == "-gapext" ]]
  then
    gapext=$4
  elif [[ $5 == "-gapext" ]]
  then
    gapext=$6
  fi
  prank -gaprate=$gaprate -gapext=$gapext -d=seq.fasta -o=$FILENAME
  mv $FILENAME.best.fas $FILENAME
  CMD="prank -gaprate=$gaprate -gapext=$gapext -d=seq.fasta -o=$FILENAME"
  generate_log "$CMD"
else
  #tells the user if the alignment tool is not supported and adds to warnings.log
  echo "$2 is not an alignment tool that is supported"
  echo "$2 is not an alignment tool that is supported" >> $CURRENTDIR/warnings.log
fi
rm seq.fasta
cd $CURRENTDIR
