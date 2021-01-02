#!/usr/bin/python
import fasta_to_paml
import sys
import os
import re
import csv
from pathlib import Path

#get the alignment file
alignment = sys.argv[1]
#get the phylo file
phylo = sys.argv[2]
#get prefix
prefix = Path(alignment).stem
#get path
path = os.path.dirname(alignment)
output_file = "dnds/" + prefix + "/" + prefix + ".paml.tre"
#create dnds directory
os.system("mkdir -p dnds")
os.system("mkdir -p dnds/" + prefix)
#copy the codeml file
os.system("cp ../codeml.ctl ./")
fasta_to_paml.convert_paml(alignment, "dnds/" + prefix+"/")

#opens the codeml.ctl file
with open("codeml.ctl", 'r') as f:
    lines = f.readlines()
    lines[0] = "      seqfile = " + alignment + "    * sequence data file name\n"
    lines[1] = "     treefile = " + phylo + "    * tree structure file name\n"
    lines[3] = "      outfile = " + output_file + "    * main result file name\n"
with open("codeml.ctl", 'w') as f:
    f.writelines(lines)
os.system("codeml")
os.system("mv rst dnds/" + prefix + "/")
#remove extra output files
os.system("rm 2NG.dN")
os.system("rm 2NG.dS")
os.system("rm 2NG.t")
os.system("rm lnf")
os.system("rm rst1")
os.system("rm rub")

#read dn/ds info
with open("dnds/" + prefix + "/rst") as f:
    lines = f.readlines()
    #initializes model 1 list
    model_1 = []
    model_2 = []
    model_1_bool = True
    for line in lines:
        if "Model 2" in line:
            model_1_bool = False
            continue
        if re.match(r'\s*\d', line) and model_1_bool:
            rows = line.split()
            del rows[4]
            rows[4] = rows[4][:-1]
            model_1.append(rows)
        if "Positively" in line:
            break
        if re.match(r'\s*\d', line) and not model_1_bool:
            rows = line.split()
            del rows[5]
            rows[5] = rows[5][:-1]
            model_2.append(rows)
    model_1_fields = ["site_pos", "ref_AA", "class_1_prob", "class_2_prob", "most_prob_class", "mean_dN/dS"]
    model_2_fields = ["site_pos", "ref_AA", "class_1_prob", "class_2_prob", "class_3_prob", "most_prob_class", "mean_dN/dS", "prob_dN/dS>1"]
    with open("dnds/" + prefix + "/" + prefix + ".site_dnds.csv", "w") as f:
            write = csv.writer(f)
            write.writerow(model_1_fields)
            write.writerows(model_1)
