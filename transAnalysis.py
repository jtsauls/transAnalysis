# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import re
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
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

def make_gene_fpkm_dict(fpkm_filename, lowcut = 0, highcut = None):
    '''Reads a fpkm tracking file and returns
    a dictionary with genes as keys and fpkms 
    as values.
    INPUT
        fpkm_filename    An fpkm_tracking file.
        lowcut           Reduce fpkm values belof this value to 0.
        highcut          Set all fpkm values above this value to x.
    OUTPUT
        gene_fpkm_dict   Dictionary as per above descritpion. 
    '''
    # Get fpkm values per gene and convert into dict. 
    fpkm_df = pd.read_csv(fpkm_filename, delimiter='\t') 
    # Process values, though not sure this is wise here
    
    
    # Convert to dict of ids and FPKM values.
    gene_fpkm_dict = dict(zip(fpkm_df.tracking_id, fpkm_df.FPKM))
    
    return gene_fpkm_dict

# <codecell>

def make_reaction_fpkm_dict(model, gene_fpkm_dict, mode=0):
    '''Create a dict with reactions as keys and fpkms as values.
    INPUT
        model              Cobra model
        gene_fpkm_dict     Dictionary mapping fpkm values to 
                           genes. The output of make_gene_fpkm_dict
    OUTPUT
        reaction_fpkm_dict
    '''
    reaction_fpkm_dict = {}
    
    if mode == 0:
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
    elif mode == 1:
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

def fpkm_comparison(fpkm_filename1, fpkm_filename2, 
                    minspan_list, model, mode=0):
    '''Graphs a comparison betweet to fpkm data set on the same 
    minspan ranking
    '''
    
    filename_list = [fpkm_filename1, fpkm_filename2]
    minspan_fpkm_vals = [0, 0]
    minspan_fpkm_ranks = [0, 0]
    x = range(1,len(minspan_list)+1)
    
    for i in range(2):
        gene_fpkm_dict = make_gene_fpkm_dict(
                         filename_list[i])
        reaction_fpkm_dict = make_reaction_fpkm_dict(
                             model, gene_fpkm_dict, mode)
        minspan_fpkm_list = make_minspan_fpkm_list(
                            minspan_list, reaction_fpkm_dict)
        # Change (value, minspan) tuple structure to two lists and
        # pull out values.
        minspan_fpkm_vals[i] = zip(*minspan_fpkm_list)[0] 
        
        # Find the rank of each minspan. 
        rank = zip(minspan_fpkm_list, x)
        rank = sorted(rank, reverse=True)
        minspan_fpkm_ranks[i] = zip(*rank)[1] # Just pull out ranks
    
    # Correlation
    m, b, r, p_val, std_err = stats.linregress(
                              minspan_fpkm_vals[0], 
                              minspan_fpkm_vals[1])
    print('R^2 for fpkm values is: ',r)
    m, b, r, p_val, std_err = stats.linregress(
                              minspan_fpkm_ranks[0], 
                              minspan_fpkm_ranks[1])
    print('R^2 for rank order is: ',r)
    
    # Plotting
    plt.subplots(2, 2, figsize=(10,10), dpi=100)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('title')
    
    # Value comparison.
    plt.subplot(2,2,1)
    plt.title('Minspan number (unranked) vs. fpkm score.')
    plt.xlabel('rank'); plt.ylabel('fpkm')
    #plt.legend(["data1", "data2"])
    plt.scatter(x, minspan_fpkm_vals[0], c='b')
    plt.scatter(x, minspan_fpkm_vals[1], c='r')
    plt.subplot(2,2,2)
    plt.title('fpkm vs. fpkm')
    plt.xlabel('fpkm 1'); plt.ylabel('fpkm 2') 
    plt.scatter(minspan_fpkm_vals[0], minspan_fpkm_vals[1], c='g')
    
    # Rank comparison
    plt.subplot(2,2,3)
    plt.title('Minspan number (unranked) vs. rank.');
    plt.xlabel('rank'); plt.ylabel('rank')
    #plt.legend(["data1", "data2"])
    plt.scatter(x, minspan_fpkm_ranks[0], c='b')
    plt.scatter(x, minspan_fpkm_ranks[1], c='r')
    plt.subplot(2,2,4)
    plt.title('rank vs rank')
    plt.xlabel('minspan rank'); plt.ylabel('minspan rank')
    plt.scatter(minspan_fpkm_ranks[0], minspan_fpkm_ranks[1], c='y')
    plt.show() # I'm not sure why fig is not needen as an arg here.

