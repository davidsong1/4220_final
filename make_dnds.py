#!/usr/bin/python
import fasta_to_paml
import sys
import os
from pathlib import Path

#get the alignment file
alignment = sys.argv[1]
#get the phylo file
phylo = sys.argv[2]
#get prefix
prefix = Path(alignment).stem
#get path
path = os.path.dirname(alignment)
output_file = path + "/dnds/" + prefix + "/" + prefix + ".paml.tre"
#create dnds directory
os.system("mkdir -p " + path + "/dnds")
os.system("mkdir -p " + path + "/dnds/" + prefix)
#copy the codeml file
#os.system("cp codeml.ctl " + path + "/dnds/" + prefix)
fasta_to_paml.convert_paml(alignment, path+"/dnds/" + prefix+"/")

#opens the codeml.ctl file
with open("codeml.ctl", 'r') as f:
    lines = f.readlines()
    lines[0] = "      seqfile = " + alignment + "    * sequence data file name\n"
    lines[1] = "     treefile = " + phylo + "    * tree structure file name\n"
    lines[3] = "      outfile = " + output_file + "    * main result file name\n"
with open("codeml.ctl", 'w') as f:
    f.writelines(lines)
os.system("codeml")
os.system("mv rst " + path + "/dnds/" + prefix + "/")
#remove extra output files
os.system("rm 2NG.dN")
os.system("rm 2NG.dS")
os.system("rm 2NG.t")
os.system("rm lnf")
os.system("rm rst1")
os.system("rm rub")
