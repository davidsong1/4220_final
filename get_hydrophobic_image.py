#!/usr/bin/python
import pymol
from pymol import cmd
import sys
import os
import __main__
__main__.pymol_argv = [ 'pymol', '-qc']
pymol.finish_launching()

# creates images directory to store pngs
os.system('mkdir -p images')

# function to highlight both hydrophobicty and charge in protein structures
def yrb(selection='all'):
    # sets colors
	cmd.remove("hydro")
	cmd.set_color('yellow',[0.950,0.78,0.0])
	cmd.set_color('grey',[0.95,0.95,0.95])
	cmd.set_color('red',[1.0,0.4,0.4])
	cmd.set_color('blue',[0.2,0.5,0.8])
    # charged oxygens of glutamate and aspartate are highlighted red and the charged nitrogens of arginine and lysine are highlighted blue
    # carbon atoms not bound to nitrogen or oxygen atoms are colored yellow, remaining atoms are white
	mapping = {}
	mapping['arg'] = [ ('NE,NH2,NH1', 'blue'), ('CD,CZ', 'grey'), ('CG', 'yellow') ]
	mapping['asn'] = [ ('CG,OD1,ND2', 'grey') ]
	mapping['asp'] = [ ('CG', 'grey'), ('OD2,OD1', 'red')  ]
	mapping['cys'] = [ ('SG', 'grey') ]
	mapping['gln'] = [ ('CG', 'yellow'), ('CD,OE1,NE2', 'grey') ]
	mapping['glu'] = [ ('CG', 'yellow'), ('CD', 'grey'), ('OE1,OE2', 'red') ]
	mapping['his'] = [ ('CG,CD2,ND1,NE2,CE1', 'grey') ]	
	mapping['ile'] = [ ('CG1,CG2,CD1', 'yellow') ]
	mapping['leu'] = [ ('CG,CD1,CD2', 'yellow') ]
	mapping['lys'] = [ ('CG,CD', 'yellow'), ('CE', 'grey'), ('NZ', 'blue') ]
	mapping['met'] = [ ('CG,CE', 'yellow'), ('SD', 'grey') ]
	mapping['phe'] = [ ('CG,CD1,CE1,CZ,CE2,CD2', 'yellow') ]
	mapping['pro'] = [ ('CG', 'yellow'), ('CD', 'grey') ]
	mapping['ser'] = [ ('CB,OG', 'grey') ]
	mapping['thr'] = [ ('CB,OG1', 'grey'), ('CG2', 'yellow') ]
	mapping['trp'] = [ ('CG,CD2,CZ2,CH2,CZ3,CE3', 'yellow'), ('CD1,NE1,CE2', 'grey') ]
	mapping['tyr'] = [ ('CG,CE1,CD1,CE2,CD2', 'yellow'), ('CZ,OH', 'grey') ]
	mapping['val'] = [ ('CG1,CG2', 'yellow') ]

	obj_list = cmd.get_names('objects')
	for obj in obj_list:
		if (obj == selection or selection == 'all'):
			cmd.color('grey','(n. N,C,CA,O and ' + obj + ')')
			cmd.color('yellow','(n. CB and ' + obj + ')')
			for key in mapping:
				for (atom, color) in mapping[key]:
					cmd.color(color, '( n. ' + atom + ' and r. ' + key + ' and ' + obj + ' )')

# checks if the argument is the envelope protein
if sys.argv[1] == "E":
    #loads in the pdb file
    cmd.load("./pdb/7k3g.pdb")
    #calls the coloring function
    yrb()
    #shows the surface
    cmd.hide("all")
    cmd.show("surface")
    cmd.set_view ("\
         0.545077741,    0.838331640,   -0.009660738,\
         0.245126382,   -0.170380458,   -0.954402864,\
        -0.801753283,    0.517852783,   -0.298365384,\
         0.000000000,    0.000000000, -220.324935913,\
         25.228977203,    0.014108658,    0.018569946,\
         184.920150757,  255.729721069,  -20.000000000")
    os.system('mkdir -p ./images/E')
    cmd.png("./images/E/envelope_protein_1.png")
    cmd.load("./pdb/7k3g.pdb")
    yrb()
    cmd.hide("all")
    cmd.show("surface")
    cmd.set_view ("\
         0.001396583,    0.035106514,    0.999381900,\
         0.997273207,    0.073681392,   -0.003983754,\
        -0.073777318,    0.996663272,   -0.034911018,\
         0.000000000,    0.000000000, -220.324935913,\
       25.228977203,    0.014108658,    0.018569946,\
      -497.505432129,  938.155151367,  -20.000000000")
    cmd.png("./images/E/envelope_protein_2.png")
    cmd.quit()

#checks if the input is the spike protein
elif sys.argv[1] == "S":
    cmd.load("./pdb/6vxx.pdb")
    yrb()
    cmd.hide("all")
    cmd.show("surface")
    cmd.set_view ("\
         0.065237045,   -0.076570503,   -0.994927406,\
        -0.246505573,   -0.967386425,    0.058287665,\
        -0.966942012,    0.241452962,   -0.081985921,\
         0.000000000,    0.000000000, -515.517028809,\
       209.999954224,  210.000427246,  206.595550537,\
       406.437561035,  624.596496582,  -20.000000000")
    os.system('mkdir -p ./images/S')
    cmd.png("./images/S/spike_protein.png")
    cmd.quit()

#checks if the input is the nucleocapsid protein
elif sys.argv[1] == "N":
    cmd.load("./pdb/6wzo.pdb")
    yrb()
    cmd.hide("all")
    cmd.show("surface")
    cmd.set_view ("\
        -0.533524752,   -0.667001605,   -0.520050287,\
        -0.639863133,   -0.083785713,    0.763901412,\
        -0.553101182,    0.740324259,   -0.382091016,\
         0.000000000,    0.000000000, -238.985290527,\
        -1.367336273,   27.149894714,    7.354507446,\
       188.417816162,  289.552764893,  -20.000000000")
    os.system('mkdir -p ./images/N')
    cmd.png("./images/N/nucleocapsid_protein.png")
    cmd.quit()
