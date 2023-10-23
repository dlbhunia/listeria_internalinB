# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 11:44:45 2022

@author: Harrison Helmick
"""

# import the necessary libraries from Chimera
from chimera import openModels, Molecule
from chimera import runCommand as rc
from chimera import selection 

# define the location and output file name
outf = open(r'PATH\TO\example.txt','w')

#outf = open(r'C:\Users\hhelmick\Desktop\gromacs\dongqi\STAR_METHODS\output\example.txt','w')

# define the open model as the model of interest
m = openModels.list(modelTypes=[Molecule])[0]

# this is a for loop that operates on each from of the MD simulation
for fn, cs in m.coordSets.items():    
    # define the active coordset you want
    m.activeCoordSet = cs
    # select the residues of the internalin
    rc('sel :2550-2788')
    # name that selection for the zone command
    rc('namesel int')
    # unselect everything
    rc('~sel')
    # select atoms within 5A of the internalin B
    rc('zone int 5')
    # unselect the ions and waters, they are not the relevant atoms
    rc('~sel :ion :W')
    # unselect anything that is part of the internalin chain
    rc('~sel :2550-2788')
    # generate labels for exports
    rc('rlabel sel')
    # define this selection to output the atom names
    temp = selection.currentResidues()

    # write all the atoms to a list that will write to the text file
    l1 = []
    for ele in temp:
        l1.append(ele.label)
    
    outf.write('{}'.format(l1) + '\n')    

outf.close()
