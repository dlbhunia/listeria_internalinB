# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 09:49:05 2022

@author: Harrison Helmick 
"""

'''
This analysis file is to be used with the outputs of Chain_distances.py
'''

# In[IMPORT LIBRARIES]
import pandas as pd
import matplotlib.pyplot as plt

# In[READ AND CLEAN THE DATA]

# known simulation time from the original
simulation_total_time = 208 #ns
# number of frames that is opened in Chimera
frames = 6403

# Read in the File that was generated from the Chimera Analysis
df = pd.read_csv(r'PATH\TO\example.txt', header = None, sep = '\n')

# the \ symbol before a symbol, as in \[ indicates to look for that symbol and drop it. Otherwise the symbol is ignored
df[0] = df[0].str.replace(r'[\[\]]', '')
df[0] = df[0].str.replace('u', '')
df1 = df[0].str.split(',', expand = True)

# get a time column based on how long your simulation was and the number of frames loaded
df1['time'] = df1.index * (simulation_total_time/frames)

# In[MAKING A DF OF A PARTICULAR CHAIN IF NEEDED]

# redefine the df so you don't lose your original df one
df2 = df1
# drop the time axis, since this breaks the extraction line
df2 = df2.drop(['time'], axis = 1)
# mask everything that isn't the E chain. This is a Boolean DF that you can use to drop everything from the other df
bChain = df2.apply(lambda row: row.astype(str).str.contains('\.B'), axis=1)
cChain = df2.apply(lambda row: row.astype(str).str.contains('\.C'), axis=1)
# drop everything that is masked
bChainDF = df2[bChain]
cChainDF = df2[cChain]

# put the unmasked data into new DFs
bChainDF = pd.DataFrame(bChainDF)
cChainDF = pd.DataFrame(cChainDF)

# In[MASKING RESIDUES TO PLOT SPECIFIC CATEGORIES]

######### NEGATIVE COUNTS ##########
BmaskASP = bChainDF.apply(lambda row: row.astype(str).str.contains('ASP'), axis=1)
CmaskASP = cChainDF.apply(lambda row: row.astype(str).str.contains('ASP'), axis=1)

BcountASP = BmaskASP.sum(True)
CcountASP = CmaskASP.sum(True)

countASP = BcountASP + CcountASP

BmaskGLU = bChainDF.apply(lambda row: row.astype(str).str.contains('GLU'), axis=1)
CmaskGLU = cChainDF.apply(lambda row: row.astype(str).str.contains('GLU'), axis=1)

BcountGLU = BmaskGLU.sum(True)
CcountGLU = CmaskGLU.sum(True)

countASP = BcountASP + CcountASP
countGLU = BcountGLU + CcountGLU

neg = countASP + countGLU

######### POSITIVE COUNTS ##########
### LYS ###
BmaskLYS = bChainDF.apply(lambda row: row.astype(str).str.contains('LYS'), axis=1)
CmaskLYS = cChainDF.apply(lambda row: row.astype(str).str.contains('LYS'), axis=1)

BcountLYS = BmaskLYS.sum(True)
CcountLYS = CmaskLYS.sum(True)

countLYS = BcountLYS + CcountLYS

### ARG ###
BmaskARG = bChainDF.apply(lambda row: row.astype(str).str.contains('ARG'), axis=1)
CmaskARG = cChainDF.apply(lambda row: row.astype(str).str.contains('ARG'), axis=1)

BcountARG = BmaskARG.sum(True)
CcountARG = CmaskARG.sum(True)

countARG = BcountASP + CcountASP

pos = countLYS + countARG

# electostatically charged interacting residues
electro = pos + neg

######### HYDROPHOBIC COUNTS ##########
BmaskVAL = bChainDF.apply(lambda row: row.astype(str).str.contains('VAL'), axis=1)
CmaskVAL = cChainDF.apply(lambda row: row.astype(str).str.contains('VAL'), axis=1)

BmaskLEU = bChainDF.apply(lambda row: row.astype(str).str.contains('LEU'), axis=1)
CmaskLEU = cChainDF.apply(lambda row: row.astype(str).str.contains('LEU'), axis=1)

BmaskILE = bChainDF.apply(lambda row: row.astype(str).str.contains('ILE'), axis=1)
CmaskILE = cChainDF.apply(lambda row: row.astype(str).str.contains('ILE'), axis=1)

BmaskPHE = bChainDF.apply(lambda row: row.astype(str).str.contains('PHE'), axis=1)
CmaskPHE = cChainDF.apply(lambda row: row.astype(str).str.contains('PHE'), axis=1)

BmaskCYS = bChainDF.apply(lambda row: row.astype(str).str.contains('CYS'), axis=1)
CmaskCYS = cChainDF.apply(lambda row: row.astype(str).str.contains('CYS'), axis=1)

BmaskMET = bChainDF.apply(lambda row: row.astype(str).str.contains('MET'), axis=1)
CmaskMET = cChainDF.apply(lambda row: row.astype(str).str.contains('MET'), axis=1)

BmaskALA = bChainDF.apply(lambda row: row.astype(str).str.contains('ALA'), axis=1)
CmaskALA = cChainDF.apply(lambda row: row.astype(str).str.contains('ALA'), axis=1)

BcountVAL = BmaskVAL.sum(True)
CcountVAL = CmaskVAL.sum(True)

BcountLEU = BmaskLEU.sum(True)
CcountLEU = CmaskLEU.sum(True)

BcountILE = BmaskILE.sum(True)
CcountILE = CmaskILE.sum(True)

BcountPHE = BmaskPHE.sum(True)
CcountPHE = CmaskPHE.sum(True)

BcountCYS = BmaskCYS.sum(True)
CcountCYS = CmaskCYS.sum(True)

BcountMET = BmaskMET.sum(True)
CcountMET = CmaskMET.sum(True)

BcountALA = BmaskALA.sum(True)
CcountALA = CmaskALA.sum(True)

hydrophobic = BcountVAL + CcountVAL + BcountLEU + CcountLEU + BcountILE + CcountILE + BcountPHE + CcountPHE + BcountCYS + CcountCYS + BcountMET + CcountMET + BcountALA + CcountALA

######### OTHER COUNTS ##########
BmaskASN = bChainDF.apply(lambda row: row.astype(str).str.contains('ASN'), axis=1)
CmaskASN = cChainDF.apply(lambda row: row.astype(str).str.contains('ASN'), axis=1)

BmaskGLN = bChainDF.apply(lambda row: row.astype(str).str.contains('GLN'), axis=1)
CmaskGLN = cChainDF.apply(lambda row: row.astype(str).str.contains('GLN'), axis=1)

BmaskGLY = bChainDF.apply(lambda row: row.astype(str).str.contains('GLY'), axis=1)
CmaskGLY = cChainDF.apply(lambda row: row.astype(str).str.contains('GLY'), axis=1)

BmaskHIS = bChainDF.apply(lambda row: row.astype(str).str.contains('HIS'), axis=1)
CmaskHIS = cChainDF.apply(lambda row: row.astype(str).str.contains('HIS'), axis=1)

BmaskPRO = bChainDF.apply(lambda row: row.astype(str).str.contains('PRO'), axis=1)
CmaskPRO = cChainDF.apply(lambda row: row.astype(str).str.contains('PRO'), axis=1)

BmaskSER = bChainDF.apply(lambda row: row.astype(str).str.contains('SER'), axis=1)
CmaskSER = cChainDF.apply(lambda row: row.astype(str).str.contains('SER'), axis=1)

BmaskTHR = bChainDF.apply(lambda row: row.astype(str).str.contains('THR'), axis=1)
CmaskTHR = cChainDF.apply(lambda row: row.astype(str).str.contains('THR'), axis=1)

BmaskTYR = bChainDF.apply(lambda row: row.astype(str).str.contains('TYR'), axis=1)
CmaskTYR = cChainDF.apply(lambda row: row.astype(str).str.contains('TYR'), axis=1)

BmaskTRP = bChainDF.apply(lambda row: row.astype(str).str.contains('TRP'), axis=1)
CmaskTRP = cChainDF.apply(lambda row: row.astype(str).str.contains('TRP'), axis=1)

BcountASN = BmaskASN.sum(True)
CcountASN = CmaskASN.sum(True)

BcountGLN = BmaskGLN.sum(True)
CcountGLN = CmaskGLN.sum(True)

BcountGLY = BmaskGLY.sum(True)
CcountGLY = CmaskGLY.sum(True)

BcountHIS = BmaskHIS.sum(True)
CcountHIS = CmaskHIS.sum(True)

BcountPRO = BmaskPRO.sum(True)
CcountPRO = CmaskPRO.sum(True)

BcountSER = BmaskSER.sum(True)
CcountSER = CmaskSER.sum(True)

BcountTHR = BmaskTHR.sum(True)
CcountTHR = CmaskTHR.sum(True)

BcountTYR = BmaskTYR.sum(True)
CcountTYR = CmaskTYR.sum(True)

BcountTRP = BmaskTRP.sum(True)
CcountTRP = CmaskTRP.sum(True)

other = BcountASN + CcountASN + BcountGLN + CcountGLN + BcountGLY + CcountGLY + BcountHIS + CcountHIS + BcountPRO + CcountPRO + BcountSER + CcountSER + BcountTHR + CcountTHR + BcountTYR + CcountTYR + BcountTRP + CcountTRP 

# In[APPLY A ROLLING MEAN]

# apply the rolling mean for data smoothing
df1['posSMA'] = pos.rolling(10).mean()
df1['negSMA'] = neg.rolling(10).mean()
df1['phobicSMA'] = hydrophobic.rolling(10).mean()
df1['otherSMA'] = other.rolling(10).mean()

# In[GENERATE THE NORMALIZED VALUES]

normDF = pd.DataFrame()

# lenght of chain being analyzed and count of amino acid / amino acid types
# length of the B and C chains together
chain_length = 818
# count of Arg in B and C chain. Other counts are simliarly indicated. Phobic is sum of hydrophobic. Other is from the "other" category above
ARG = 90
LYS = 68
ASP = 58
GLU = 50
PHOBIC = 332
OTHER = chain_length - ARG - LYS - ASP - GLU - PHOBIC

# normalize the count of the negative and positive amino acids to their frequencies in the B and C chain
# and apply a rolling mean
negNorm = (countGLU + countASP) / (GLU + ASP) *100
posNorm = (countARG + countLYS)/(ARG + LYS) *100
normDF['neg'] = negNorm
normDF['pos'] = posNorm
normDF['negMSA'] = negNorm.rolling(10).mean()
normDF['posMSA'] = posNorm.rolling(10).mean()

# normalize the count of the hydrophobic amino acids to their frequencies in the B and C chain
# and apply a rolling mean
hydrophobicNorm = (hydrophobic / PHOBIC) * 100
normDF['phobic'] = hydrophobicNorm
normDF['phobicMSA'] = hydrophobicNorm.rolling(10).mean()

# normalize the count of the other amino acids to their frequencies in the B and C chain
# and apply a rolling mean
otherNorm = (other / chain_length) * 100
normDF['other'] = otherNorm
normDF['otherMSA'] = otherNorm.rolling(10).mean()

# add in the time, since that does not change with normalization 
normDF['time'] = df1['time']

# In[PLOTTING INDIVIDUAL INTEACTIONS- COUNTS]

plt.figure(0)
plt.style.use('seaborn')
plt.plot(df1['time'], pos, alpha = 0.5)
plt.plot(df1['time'], df1['posSMA'], color = 'blue', label = 'Positive Interacting Residues', linestyle = '-')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('Count of LAP Residues within 5Å of InIB')
plt.rc('font', size = 20)
plt.legend(bbox_to_anchor = (1,1))

plt.figure(1)
plt.style.use('seaborn')
plt.plot(df1['time'], neg, alpha = 0.5, color = 'red')
plt.plot(df1['time'], df1['negSMA'], color = 'red', label = 'Negative Interacting Residues', linestyle = '-')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('Count of LAP Residues within 5Å of InIB')
plt.rc('font', size = 20)
plt.legend(bbox_to_anchor = (1,1))

plt.figure(2)
plt.plot(df1['time'], hydrophobic, alpha = 0.5, color = 'green')
plt.plot(df1['time'], df1['phobicSMA'], color = 'green', label = 'Hydrophobic Interacting Residues', linestyle  = '-')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('Count of LAP Residues within 5Å of InIB')
plt.rc('font', size = 20)
plt.legend(bbox_to_anchor = (1,1))

plt.figure(3)
plt.plot(df1['time'], other, alpha = 0.5, color = 'grey')
plt.plot(df1['time'], df1['otherSMA'], color = 'grey', label = 'Other Interacting Residues', linestyle = '-')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('Count of LAP Residues within 5Å of InIB')
plt.rc('font', size = 20)
plt.legend(bbox_to_anchor = (1,1))

# In[PLOTTING ALL INTERACTIONS- COUNTS]

plt.figure(4)
plt.style.use('seaborn')
plt.plot(df1['time'], pos, alpha = 0.5)
plt.plot(df1['time'], df1['posSMA'], color = 'blue', label = 'Positive Interacting Residues', linestyle = '-')

plt.style.use('seaborn')
plt.plot(df1['time'], neg, alpha = 0.5, color = 'red')
plt.plot(df1['time'], df1['negSMA'], color = 'red', label = 'Negative Interacting Residues', linestyle = '-')

plt.plot(df1['time'], hydrophobic, alpha = 0.5, color = 'green')
plt.plot(df1['time'], df1['phobicSMA'], color = 'green', label = 'Hydrophobic Interacting Residues', linestyle  = '-')

plt.plot(df1['time'], other, alpha = 0.5, color = 'grey')
plt.plot(df1['time'], df1['otherSMA'], color = 'grey', label = 'Other Interacting Residues', linestyle = '-')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('Count of LAP Residues within 5Å of InIB')
plt.rc('font', size = 20)
plt.legend(bbox_to_anchor = (1,1))

# In[PLOTTING ALL INTERACTIONS- NORMALIZED]

plt.figure(5)
plt.style.use('seaborn')
plt.plot(normDF['time'], normDF['pos'], alpha = 0.5)
plt.plot(normDF['time'], normDF['posMSA'], color = 'blue', label = 'Positive Interacting Residues', linestyle = '-')

plt.style.use('seaborn')
plt.plot(normDF['time'], normDF['neg'], alpha = 0.5, color = 'red')
plt.plot(normDF['time'], normDF['negMSA'], color = 'red', label = 'Negative Interacting Residues', linestyle = '-')

plt.plot(normDF['time'], normDF['phobic'], alpha = 0.5, color = 'green')
plt.plot(normDF['time'], normDF['phobicMSA'], color = 'green', label = 'Hydrophobic Interacting Residues', linestyle  = '-')

plt.plot(normDF['time'], normDF['other'], alpha = 0.5, color = 'grey')
plt.plot(normDF['time'], normDF['otherMSA'], color = 'grey', label = 'Other Interacting Residues', linestyle = '-')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('Percent of LAP Residues within 5Å of InIB')
plt.rc('font', size = 20)
plt.legend(bbox_to_anchor = (1,1))

plt.savefig(r'PATH\TO\OUTPUT\norm_example.png', bbox_inches = 'tight', dpi = 200)
