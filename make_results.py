#!/usr/bin/python
import os
import subprocess
with open("README.md", "w") as f:
    with open("commands.txt", "r") as c:
        #get_seq.sh command
        seq = c.readline()
        f.write(seq)
        #list output files from get_seq
        seq_directory = subprocess.run(["find", seq.split()[-2], "-type", "f"], stdout=subprocess.PIPE).stdout.decode('utf-8')
        f.write(seq_directory)
        #make_align command
        f.write(c.readline())
        #list output files from make_align
        align_directory = subprocess.run(["find", "alignment", "-type", "f"], stdout=subprocess.PIPE).stdout.decode('utf-8')
        f.write(align_directory)
        #make_phylo command
        f.write(c.readline())
        #list output files from make_phylo
        phylo_directory = subprocess.run(["find", "phylo", "-type", "f"], stdout=subprocess.PIPE).stdout.decode('utf-8')
        f.write(phylo_directory)
        #make_mol_stats.py command
        f.write(c.readline())
        #list output files from make_mol_stats
        phylo_directory = subprocess.run(["find", "mol_stats", "-type", "f"], stdout=subprocess.PIPE).stdout.decode('utf-8')
        f.write(phylo_directory)
        #make_dnds command
        f.write(c.readline())
        #list output files from make_dnds
        dnds_directory = subprocess.run(["find", "dnds", "-type", "f"], stdout=subprocess.PIPE).stdout.decode('utf-8')
        f.write(dnds_directory)
        #get_hydrophobic_image.py command
        f.write(c.readline())
        #generate_hydrophilicity_curve.py
        f.write(c.readline())
        images_dir = subprocess.run(["find", "images", "-type", "f"], stdout=subprocess.PIPE).stdout.decode('utf-8')
        f.write(images_dir)
        c.close()
    f.close()
