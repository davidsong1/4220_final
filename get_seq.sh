#!/bin/bash
# read in accession list
ACCESSIONS=`cat $1 | cut -d"," -f2`
SEQ_DIR=$2

# enter sequence directory
mkdir -p $SEQ_DIR
cd $SEQ_DIR

# get each accession
for ACC in $ACCESSIONS;
do
    # check if the fasta file already exists and overwrite argument is not true
    if [[ -f "$ACC.fasta" ]] && [[ $3 != true ]]
    then
        # check if first line starts with the file description and line 2 contains nucleotides
        if [[ `cat $ACC.fasta | sed -n '1p'` =~ ^'>' ]] && [[ `cat $ACC.fasta | sed -n '2p'` =~ ^[ACTG]+$ ]]
        then
            echo "$ACC.fasta exists and is valid"
        else
            #deletes the invalid accession file and downloads it
            echo "$ACC.fasta exists but is not valid"
            rm $ACC.fasta
            echo "$ACC.fasta was deleted" >> ../warnings.log
            efetch -db nuccore -format fasta_cds_na -id $ACC > $ACC.fasta
        fi
    else
        # fetch the accession, separating each coding gene into its own fasta entry
        echo "Downloading $ACC"
        efetch -db nuccore -format fasta_cds_na -id $ACC > $ACC.fasta
    fi
done

# loop through all fasta files in seq_dir
for FILE in `ls *.fasta`;
do
    # get accession from filename
    ACC=`echo $FILE | cut -d '.' -f 1`
    # transform text, using `!` as special delimiter, to be cut later
    TEXT=`cat $FILE | tr -d '\n' | sed 's/>/\n>/g' | rev | sed 's/]/!]/' | rev`
    # loop through each entry
    IFS=$'\n'
    for LINE in $TEXT;
    do
        # get each gene name (field 1)
        GENE=`echo $LINE | cut -d '!' -f 1 | grep -o -e "\[gene=[A-Za-z0-9 ]*\]" | tr -d '[]' | cut -d '=' -f 2`
        # get the sequence info (field 2)
        SEQ=`echo $LINE | cut -d '!' -f 2`
        # if gene name is nonempty
        if [[ -n $GENE ]];
        then
            # make dir, if it doesn't exist
            mkdir -p $GENE
            # create file with fasta entry for gene
            NAME=$ACC'_'$GENE
            echo -e '> '$NAME'\n'$SEQ > $GENE/$NAME.fasta
        fi
    done
done

# loop through all gene_dirs in seq_dir, then make concat files
for GENE_DIR in `ls -d */`;
do
    GENE_DIR=`echo $GENE_DIR | tr -d '/'`
    cat $GENE_DIR/*.fasta > "all_"$GENE_DIR".fasta"
done

# return to original directory
cd ..
