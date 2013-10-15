# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import cobra as cp
import numpy as np
import pandas as pd
from pickle import load # Can't get SBML to work, using pickle

def run_from_ipython():
    '''Just checking if we can use ipython things.'''
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

if run_from_ipython():
    # I wanted to be clever but python doesn't like this even if it is in an if or try.
    #%cd /Users/jt/Desktop/palsson_rotation/transAnalysis/data
    print 'Run from ipython'
else:
    print 'Not run from ipython.'

minspan_filename = 'MinSpanPathways.xlsx'
sheetname = 'iJO1366'
fpkm_filename = 'bop27isoforms.fpkm_tracking' 
#fpkm_filename = 'anaerobe1_isoforms.fpkm_tracking'
#fpkm_filename = 'anaerobe2_isoforms.fpkm_tracking'
#fpkm_filename = 'rpoBisoforms.fpkm_tracking'
model_filename = 'iJO1366.pickle'

# <codecell>

# Get minspan.
minspan_df = pd.read_excel(minspan_filename, sheetname, index_col=None, na_values=['NA'])

# This processing should be automated.
# Remove any rows that aren't gene names (aka the last line with pathway length)
minspan_df = minspan_df.drop(minspan_df.index[-1])

# Create a list of list of all the minspans. 
minspan_list = []

for col in minspan_df.columns[2:]: # One minspan per column, starts on 3rd column.
    #             Entire col of names    Index of rows with values (not NaN)         No labels
    now_minspan = minspan_df.loc[:,'Rxn'][np.logical_not(np.isnan(minspan_df.loc[:,col]))].values
    now_minspan = now_minspan.tolist() # No need to be numpy arrays. 
    # Change from unicode... not sure if this is wise.
    now_minspan = [now_minspan[i].encode('ascii','ignore') for i in range(len(now_minspan))]
    minspan_list.append(now_minspan)
    

# <codecell>

# Get model.
f = open(model_filename, "rb")
model = load(f)
f.close()

# <codecell>

# Get fpkm values per gene and convert into dict. 
fpkm_df = pd.read_csv(fpkm_filename, delimiter='\t') 
# Convert to dict of ids and FPKM values.
gene_fpkm_dict = dict(zip(fpkm_df.tracking_id, fpkm_df.FPKM))

# <codecell>

# Create a dict with reactions as keys and fpkms as values.
reaction_fpkm_dict = {} # Will hold fpkm values per reaction based on combination of genes.

for reaction in model.reactions:
    gene_list = reaction.get_gene()
    # Combine fpkm data for involved genes. 
    fpkm_combined = 0.0
    for gene in gene_list:
        try:
            fpkm_combined += gene_fpkm_dict[str(gene)] # Just add them. This makes since for or, but not and...
        except:
            pass
            #print "%s (from %s) is not in fpkm_gene_dict." % (str(gene), str(reaction.id))
    reaction_fpkm_dict[reaction.id] = fpkm_combined

# <codecell>

# Now to iterate over minspans, averaging fpkm values from involved reactions. 
minspan_fpkm_list = []
for minspan in minspan_list:
    minspan_fpkm = 0.0
    no_data_count = 0 # Don't average in 0's where data does not exist. 
    for reaction in minspan:
        try: 
            minspan_fpkm += fpkm_reaction_dict[reaction]
        except: # A lot don't match becaues of SBML naming vs names from the minspan file
            no_data_count += 1
            #print "%s is not in fpkm_reaction_dict." % str(reaction)
    try:
        minspan_fpkm = minspan_fpkm / (len(minspan) - no_data_count) # Just aveage for now.
    except:
        # There is no data for this minspan
        pass
    minspan_fpkm_list.append(minspan_fpkm)

# Create list of tuples with fpkm values and then the minspan list
ranked_minspans = zip(minspan_fpkm_list, minspan_list)
ranked_minspans = sorted(ranked_minspans, reverse=True)
ranked_minspans = zip(range(1,len(ranked_minspans)+1), ranked_minspans)

# <codecell>

# Print to file
output_filename = 'output_minspan_rank_%s.txt' % fpkm_filename
f = open(output_filename, 'w')
for item in ranked_minspans:
    f.write("%s\n" % str(item))
f.close()

# <codecell>

# Just for fun find the fpkm ranking per reacion. 
ranked_reactions = [(b, a) for a, b in reaction_fpkm_dict.items()]
ranked_reactions = sorted(ranked_reactions, reverse=True)

output_filename = 'output_gene_rank_%s.txt' % fpkm_filename
f = open(output_filename, 'w')
for item in ranked_reactions:
    f.write("%s\n" % str(item))
f.close()

# <codecell>

# scrap 

#rxn_test = model.reactions.get_by_id('GLCptspp')
#print rxn_test
#print rxn_test.get_gene(), type(rxn_test.get_gene())
#print rxn_test.annotation
#print rxn_test.gene_reaction_rule, type(rxn_test.gene_reaction_rule)
#print rxn_test.notes
#print rxn_test.parse_gene_association(rxn_test.gene_reaction_rule)

#fpkm_reaction_list.append(fpkm_combined)
#print type(fpkm_reaction_list), fpkm_reaction_list

# Marry reaction names and their values to make dictionary 
#fpkm_reaction_dict = dict(zip(iJO1366.reactions, fpkm_reaction_list))

# Format and type check.     
#print minspan_list[:5], type(minspan_list[:5])

#print ranked_minspans[:5]
#type(str(ranked_minspans[1]))
#print str(ranked_minspans[1])

#print fpkm_minspan_list

    #print fpkm_combined

#print fpkm_reaction_dict

#print("gene:", gene, "fpkm:", fpkm_gene_dict[str(gene)])

#print(gene_list)

# Format and type check.     
#print minspan_list[:5], type(minspan_list[:5])

#print(fpkms.loc[:10,['tracking_id','FPKM']], type(fpkms.loc[:10,['tracking_id','FPKM']])) # Looking at what matters.

#print fpkm_gene_dict['b2215']
#fpkm_gene_dict['s0001']

