# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 10:59:30 2021

@author: Harrison Helmick
"""

# In[IMPORT LIBRARIES, DEFINE FUNCTIONS]
import pandas as pd

def readNclean (path):
    df = pd.read_csv(path, sep = '\n', header = None)
    df1 = df[0].str.split(',', expand = True)
    df2 = df1.transpose()
    return df2
    

# gets the counts from the dataframe you put in
def get_counts(df):
    cols = df.columns

    l1 = []
    for col in cols:
        l1.append(df[col])
    
    l2 = []
    for ele in l1:
        for i in ele:
            l2.append(i)
    
    df1 = pd.DataFrame(l2, columns = ['acids'])
    df1['value_counts'] = df1['acids'].str.strip('[]')
    df2 = df1.dropna()
    df3 = df2['value_counts'].value_counts().sort_values(ascending = False)

    return df3

# gets the counts from the dataframe you put in
def list_columns(df):
    cols = df.columns

    l1 = []
    for col in cols:
        l1.append(df[col])
    
    l2 = []
    for ele in l1:
        for i in ele:
            l2.append(i)
    return l2

# In[RUN THE CODE]

# this should be the file that was exported from chimera that you generated from the close_to_internalin file
path_in = PATH\TO\close.csv
clean_data = readNclean(balance_in)
frequency = get_counts(clean_data) 
# this should be the place you want to export the cleaned file
frequency.to_csv(PATH\TO\close_processed.csv)















































































