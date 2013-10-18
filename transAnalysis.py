# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import re
import numpy as np
import pandas as pd
import cobra as cp

# <codecell>

def make_minspan_list(minspan_filename):
    '''Reads the minspan .xls file returns minspans
    INPUT
        minspan_filename    .xls with minspans
    OUTPUT
        minspan_list        A list of list of minspans
    '''
    minspan_list = []
    sheetname = 'iJO1366' # Could make this an arg. 
    
    # Get minspan.
    minspan_df = pd.read_excel(minspan_filename, sheetname, 
                               index_col=None, na_values=['NA'])
    # Remove any rows that aren't gene names 
    # (aka the last line with pathway length).
    minspan_df = minspan_df.drop(minspan_df.index[-1])
    
    # One minspan per column, starts on 3rd column.
    for col in minspan_df.columns[2:]: 
        # Index of rows with values (not NaN)  
        value_index = np.isnan(minspan_df.loc[:,col])
        value_index = np.logical_not(value_index)
        now_minspan = minspan_df.loc[:,'Rxn'][value_index]
        now_minspan = now_minspan.values # No labels
        now_minspan = now_minspan.tolist() # No numpy arrays. 
        # Change from unicode... not sure if this is wise.
        now_minspan = [now_minspan[i].encode('ascii','ignore') \
                       for i in range(len(now_minspan))]
        minspan_list.append(now_minspan)
            
    return minspan_list

# <codecell>

def make_gene_fpkm_dict(fpkm_filename):
    '''Reads a fpkm tracking file and returns
    a dictionary with genes as keys and fpkms 
    as values.
    '''
    # Get fpkm values per gene and convert into dict. 
    fpkm_df = pd.read_csv(fpkm_filename, delimiter='\t') 
    # Convert to dict of ids and FPKM values.
    gene_fpkm_dict = dict(zip(fpkm_df.tracking_id, fpkm_df.FPKM))
    
    return gene_fpkm_dict

# <codecell>

def make_reaction_fpkm_dict(model, gene_fpkm_dict):
    '''Create a dict with reactions as keys and fpkms as values.
    INPUT
        model              Cobra model
        gene_fpkm_dict     Dictionary mapping fpkm values to 
                           genes. The output of make_gene_fpkm_dict
    OUTPUT
        reaction_fpkm_dict
    '''
    reaction_fpkm_dict = {}
    
    for reaction in model.reactions:
        gene_bool = reaction.gene_reaction_rule
        gene_bool = re.sub('\(|\)', ' ', gene_bool)
        gene_bool = gene_bool.split('or')
        fpkm_combined = 0
        for orgroup in gene_bool:
            andgroup = orgroup.split('and')
            andgroup_fpkm = 0
            avgby = len(andgroup)
            for genename in andgroup:
                genename = genename.strip()
                try:
                    andgroup_fpkm += gene_fpkm_dict[genename]
                except:
                    avgby -= 1
                    #print "%s is not in reaction_fpkm_dict." \
                    #% str(genename)
            try:
                andgroup_fpkm = andgroup_fpkm / avgby
            except:
                pass
            fpkm_combined += andgroup_fpkm
        reaction_fpkm_dict[reaction.id] = fpkm_combined
    return reaction_fpkm_dict

# <codecell>

#THIS IS THE OLD ONE THAT JUST ADDS THEM UP
def make_reaction_fpkm_dict1(model, gene_fpkm_dict):
    '''Create a dict with reactions as keys and fpkms as values.
    INPUT
        model              Cobra model
        gene_fpkm_dict     Dictionary mapping fpkm values to 
                           genes. The output of make_gene_fpkm_dict
    OUTPUT
        reaction_fpkm_dict
    '''
    reaction_fpkm_dict = {}
    
    for reaction in model.reactions:
        gene_list = reaction.get_gene()
        # Combine fpkm data for involved genes. 
        fpkm_combined = 0.0
        for gene in gene_list:
            try:
                # Just add them. This makes since for or, not and
                fpkm_combined += gene_fpkm_dict[str(gene)]
            except:
                pass
                #print "%s (from %s) is not in fpkm_gene_dict." % \
                #(str(gene), str(reaction.id))
        reaction_fpkm_dict[reaction.id] = fpkm_combined
    
    return reaction_fpkm_dict

# <codecell>

def make_minspan_fpkm_list(minspan_list, reaction_fpkm_dict, 
                           sort=False):
    '''Rank minspans based on the average fpkm values of their
    constituite reactions.
    INPUT
        minspan_list
        reaction_fpkm_dict
        sort                 Indcates if the minspans should be ranked.
    OUTPUT
        ranked_minspans
    '''
    # Iterate over minspans, avg fpkm values from involved reactions. 
    minspan_fpkm_list = []
    for minspan in minspan_list:
        minspan_fpkm = 0.0
        no_data_count = 0 # Don't average in 0's. 
        for reaction in minspan:
            try: 
                minspan_fpkm += reaction_fpkm_dict[reaction]
            except: 
                # A lot don't match becaues of SBML naming vs names 
                # from the minspan file.
                no_data_count += 1
                #print "%s is not in reaction_fpkm_dict." \
                #% str(reaction)
        try:
             # Just average for now.
            minspan_fpkm = minspan_fpkm / (len(minspan) - 
                                           no_data_count)
        except:
            # There is no data for this minspan
            pass
        minspan_fpkm_list.append(minspan_fpkm)
    
    # Create list of tuples with fpkm values and then the minspan list
    minspan_fpkm_list = zip(minspan_fpkm_list, minspan_list)
    if sort:
        minspan_fpkm_list = sorted(minspan_fpkm_list, reverse=True)
        
    return minspan_fpkm_list

# <codecell>

def rank_reactions(reaction_fpkm_dict):
    '''Ranks reactions by their fpkm value.
    INPUT
        reaction_fpkm_dict
    OUTPUT
        ranked_reactions
    '''
    ranked_reactions = [(b, a) for a, b in reaction_fpkm_dict.items()]
    ranked_reactions = sorted(ranked_reactions, reverse=True)
    
    return ranked_reactions

# <codecell>


