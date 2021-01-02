# 4220_final

pipeline.sh
Usage: ./pipeline.sh SETTINGS_FILE [JOB_DIR]
The script pipeline.sh itself runs the other pipeline step

parse_settings.sh
Usage: ./parse_settings.sh SETTINGS_FILE PIPELINE_STEP
This script will parse analysis settings from a setting file. Users will provide two arguments: 
(1) the file path to the pipeline settings file, and (2) the name of the pipeline step to parse.

get_seq.sh
Usage: ./get_seq ACCESSION_FILE SEQUENCE_DIR [OVERWRITE]
The get_seq.sh manages and downloads fasta-formatted accessions from GenBank. As input, the script accepts two arguments: 
(1) a list of accessions, and (2) a directory where the sequences are managed. 
The script will then check whether each accession has already been downloaded into the managed directory, download any missing sequences, and 
append any issues to the file warnings.log

make_align.sh
Usage: ./make_align SEQUENCE_DIR ALIGN_TOOL [ALIGN_TOOL_OPTIONS]
This script will align a set of fasta sequences located in a target directory.

make_phylo.sh
Usage: ./make_phylo ALIGN_FILE PHYLO_TOOL [PHYLO_TOOL_OPTIONS]
This script will estimate a phylogeny from a multiple sequence alignment.

make_mol_stats.py
Usage: ./make_mol_stats.py ALIGN_FILE
The make_mol_stats.py script generates a report of various summary statistics and transformations for a multiple sequence alignment.

make_dnds.py
Usage: ./make_dnds.py ALIGN_FILE PHYLO_FILE
The make_dnds.sh script will test for the molecular signature of positive selection using the modeling software, PAML.

make_results.py
Usage: ./make_results.py SEQUENCE_DIR
This script should generate a README.md file in SEQUENCE_DIR that lists the analysis settings and the output files for each step

get_hydrophobic_image.py
Usage: ./get_hydrophobic_image.py [Protein]
This script uses PyMOL to generate a hydrophobic surface of the protein specified

generate_hydrophilicity_curve.py
Usage: ./generate_hydrophilicity_curve.py [Fasta1] [Fasta2] [Start] [End]
This script uses matplotlib to generate a hydrophilicty curve comparison between two fasta files within a window
