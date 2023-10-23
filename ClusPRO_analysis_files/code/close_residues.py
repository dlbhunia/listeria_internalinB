# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 09:11:14 2021

@author: Harrison Helmick
"""

import os

import chimera
from chimera import runCommand as rc
from chimera import selection 

# This is the name of the file that this program will generate
outf = open(PATH\TO\close.csv','w')

# Change the work directory to the folder that contains the ClusPRO results. This should be a folder with 
# .pdb files that were generated
os.chdir(PATH\TO\example_data')

#find all files that end in .pdb in the working directory
files = [fn for fn in os.listdir('.') if fn.endswith('.pdb')]

#write the fist row in the file to keep track of things.  
outf.write('file,' + '\n')

# open a for loop to analyze the .pdn files in the folder
# the code in the for loop selects the residues that are within 5A of the LAP protein. 
# it then writes the names of these atoms to a .csv file. This code can be easily 
# reversed to select the LAP residues within 5A of the internalin as well. 
for fn in files:

    rc('open ' + fn)
    temp_list = []
    rc('split')
    rc('sel #0.2 & :.A')
    #rc('sel invert') # if reversing this code to select INT run this line
    rc('namesel int') # if reversing this code to select INT change to ('namesel lap')
    rc('zone int 5') # if reversing this code to select LAP change to ('zone int 5')
    rc('rlabel')
    temp = selection.currentResidues()
    
    for ele in temp:
        temp_list.append(ele.label)
    
    outf.write(fn + ',')
    outf.write(str(temp_list) + '\n')
    
    rc('close all')

rc('close all')
outf.close()

