#!/usr/bin/python
import sys
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from collections import Counter

#get the alignment file
align = sys.argv[1]
#get prefix
prefix = Path(align).stem
#create iterable container of SeqRecords
f = SeqIO.parse(align, 'fasta')
#create mol_stats directory
os.system("mkdir -p mol_stats")
#initializes sequence length
seq_length = 0
#codon table
codon_table = {
        'A': ('GCT', 'GCC', 'GCA', 'GCG'),
        'C': ('TGT', 'TGC'),
        'D': ('GAT', 'GAC'),
        'E': ('GAA', 'GAG'),
        'F': ('TTT', 'TTC'),
        'G': ('GGT', 'GGC', 'GGA', 'GGG'),
        'H': ('CAT', 'CAC'),
        'I': ('ATT', 'ATC', 'ATA'),
        'K': ('AAA', 'AAG'),
        'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
        'M': ('ATG',),
        'N': ('AAT', 'AAC'),
        'P': ('CCT', 'CCC', 'CCA', 'CCG'),
        'Q': ('CAA', 'CAG'),
        'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
        'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
        'T': ('ACT', 'ACC', 'ACA', 'ACG'),
        'V': ('GTT', 'GTC', 'GTA', 'GTG'),
        'W': ('TGG',),
        'Y': ('TAT', 'TAC')
    }
#compute the GC-richness for each sequence
with open("./mol_stats/"+prefix+".seq_GC.csv", 'w') as file:
    for row in f:
        file.write(row.id+","+str(GC(row.seq))+'\n')
        seq_length = len(row.seq)
print("Generated " +prefix+ ".seq_GC.csv")

#compute the GC-richness for each site
with open("./mol_stats/"+prefix+".site_GC.csv", 'w') as file:
    #iterates over each site
    for x in range(seq_length):
        site = ""
        f = SeqIO.parse(align, 'fasta')
        #iterates over each sequence
        for row in f:
            #adds the base from each site of the sequences
            site += row.seq[x]
        file.write(str(x)+","+str(GC(Seq(site)))+'\n')
print("Generated " +prefix+ ".site_GC.csv")

#determine if each site is phylogenetically informative
with open("./mol_stats/"+prefix+".site_phylo_inf.csv", 'w') as file:
    #iterates over each site
    for x in range(seq_length):
        site = ""
        counter = 0
        informative = False
        f = SeqIO.parse(align, 'fasta')
        #iterates over each sequence
        for row in f:
            site += row.seq[x]
        #checks if any site has more than 2 occurrences
        if Seq(site).count('A') >= 2:
            counter += 1
        if Seq(site).count('C') >= 2:
            counter += 1
        if Seq(site).count('T') >= 2:
            counter += 1
        if Seq(site).count('G') >= 2:
            counter += 1
        #checks if there are more than 2 bases that have more than 2 occurrences
        if counter >= 2:
            #sets informative boolean to true
            informative = True
        file.write(str(x)+","+str(informative)+'\n')
print("Generated " +prefix+ ".site_phylo_inf.csv")

#determine the codon frequencies for each sequence
with open("./mol_stats/"+prefix+".seq_codon.csv", 'w') as file:
    f = SeqIO.parse(align, 'fasta')
    #iterates over ever sequence
    for row in f:
        #gets a list of codons from the sequence
        codons = [str(row.seq[i:i+3]) for i in range(0, len(row.seq), 3)]
        #counts the number of occurences of each codon
        c = Counter(codons)
        #iterates through the Counter and adds the frequency to the csv
        for k,v in c.items():
            #only add valid codons
            if "-" not in k:
                file.write(row.id + "," + k + "," + str(v/len(codons)) + '\n')
print("Generated " +prefix+ ".seq_codon.csv")

#compute codon frequencies for each site
with open("./mol_stats/"+prefix+".site_codon.csv", 'w') as file:
    #iterates through codon sites
    for x in range(0, seq_length, 3):
        site_codon = []
        f = SeqIO.parse(align, 'fasta')
        #iterates through every seq
        for row in f:
            #add the codon to the site_codon list
            site_codon.append(str(row.seq[x:x+3]))
        #counts the number of occurences of each codon
        c = Counter(site_codon)
        #iterates through the Counter and adds the frequency to the csv
        for k,v in c.items():
            #only add valid codons
            if "-" not in k:
                file.write(str(int((x+3)/3)) + "," + k + "," + str(v/len(site_codon)) + '\n')
print("Generated " +prefix+ ".site_codon.csv")

#Compute biased codon usage proportions across amino acids, sites, and sequences
with open("./mol_stats/"+prefix+".codon_usage.csv", 'w') as file:
    f = SeqIO.parse(align, 'fasta')
    #iterates through every sequence
    for row in f:
        #gets every codon from the sequence
        codons = [str(row.seq[i:i+3]) for i in range(0, len(row.seq), 3)]
        keys = []
        #iterates through every codon
        for c in codons:
            #translates the codons to amino acids
            keys += [key for key, value in codon_table.items() if c in value]
        #creates a list of amino acids without duplicates as keys
        dkeys = list(set(keys))
        d = {}
        #iterates through every amino acid
        for aa in dkeys:
            freq = {}
            #iterates through codon table
            for codon in codon_table[aa]:
                #gets the frequency for each codon
                freq[codon] = codons.count(codon)/keys.count(aa)
            d[aa] = freq
        #itreates through dictionary
        for k,v in d.items():
            file.write(row.id+"," + k+","+str(v)+'\n')
print("Generated " +prefix+ ".codon_usage.csv")

