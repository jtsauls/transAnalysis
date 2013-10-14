# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import cobra as cp
import numpy as np
import pandas as pd
from pickle import load # Can't get SBML to work, using pickle

%cd /Users/jt/Desktop/palsson_rotation/data/

filename_fpkm = 'anaerobe1_isoforms.fpkm_tracking'
filename_iJO1366 = '/Users/jt/Desktop/palsson_rotation/data/iJO1366.pickle'

# <codecell>

# Get fpkm values per gene and convert into dict. 
fpkms = pd.read_csv(filename_fpkm, delimiter='\t')
#print(fpkms.loc[:10,['tracking_id','FPKM']], type(fpkms.loc[:10,['tracking_id','FPKM']])) # Looking at what matters. 
# Convert to dict
fpkm_gene_dict = dict(zip(fpkms.tracking_id, fpkms.FPKM))

fpkm_gene_dict['b2215']
#fpkm_gene_dict['s0001']

# <codecell>

# Get minspan. This is pretty hard coded... weak. 
minspan_df = pd.read_excel('MinSpanPathways.xlsx', 'iJO1366', index_col=None, na_values=['NA'])

# This processing should be automated
# Remove any rows that aren't gene names (aka the last line with pathway lenght)
minspan_df = minspan_df.drop(minspan_df.index[-1])
# Create a list of list of all the minspans. 
minspan_list = []

for col in minspan_df.columns[2:]: # One minspan per column, starts on 3 column
    #             Entire col of names    Index of rows with values (not NaN)         No labels
    now_minspan = minspan_df.loc[:,'Rxn'][np.logical_not(np.isnan(minspan_df.loc[:,col]))].values
    now_minspan = now_minspan.tolist() # No need to be numpy arrays. 
    # Change from unicode... not sure if this is wise
    now_minspan = [now_minspan[i].encode('ascii','ignore') for i in range(len(now_minspan))]
    minspan_list.append(now_minspan)
    
# Format and type check.     
#print minspan_list[:5], type(minspan_list[:5])

# <codecell>

# Get iJO
f = open(filename_iJO1366, "rb")
iJO1366 = load(f)
f.close()

# <codecell>

# Create a dict with reactions as keys and fpkms as values.
fpkm_reaction_dict = {} # Will hold fpkm values per reaction based on combination of genes.

for reaction in iJO1366.reactions:
    gene_list = reaction.get_gene()
    #print(gene_list)
    # Combine fpkm data for involved genes. 
    fpkm_combined = 0.0
    for gene in gene_list:
        #print("gene:", gene, "fpkm:", fpkm_gene_dict[str(gene)])
        try:
            fpkm_combined += fpkm_gene_dict[str(gene)] # Just add them. This makes since for or, but not and...
        except:
            pass
            #print "%s (from %s) is not in fpkm_gene_dict." % (str(gene), str(reaction.id))
    fpkm_reaction_dict[reaction.id] = fpkm_combined
    #print fpkm_combined

#print fpkm_reaction_dict

# <codecell>

# Now to iterate over minspans, averaging fpkm values from involved reactions. 
fpkm_minspan_list = []
for minspan in minspan_list:
    minspan_fpkm = 0.0
    for reaction in minspan:
        try: 
            minspan_fpkm += fpkm_reaction_dict[reaction]
        except: # A lot don't match becaues of SBML naming vs names from the minspan file
            pass
            #print "%s is not in fpkm_reaction_dict." % str(reaction)
    minspan_fpkm = minspan_fpkm / len(minspan) # Just aveage for now. 
    fpkm_minspan_list.append(minspan_fpkm)

#print fpkm_minspan_list

# Create list of tuples with fpkm values and then the minspan list
ranked_minspans = zip(fpkm_minspan_list, minspan_list)
ranked_minspans = sorted(ranked_minspans, reverse=True)
print ranked_minspans[:5]

# <codecell>


# <codecell>


# <codecell>

# scrap 

#fpkm_reaction_list.append(fpkm_combined)
#print type(fpkm_reaction_list), fpkm_reaction_list

# Marry reaction names and their values to make dictionary 
#fpkm_reaction_dict = dict(zip(iJO1366.reactions, fpkm_reaction_list))

