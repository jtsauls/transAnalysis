# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>
from __future__ import division
import re
import math
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import cobra as cp  # Used in Eflux
#import fluxAnalysis
import pfba_example


def make_minspan_list(minspan_filename, sheetname='iJO1366'):
    '''Reads the minspan .xls file returns minspans
    INPUT
        minspan_filename    .xls with minspans
    OUTPUT
        minspan_list        A list of list of minspans
    '''
    minspan_list = []

    # Get minspan.
    minspan_df = pd.read_excel(minspan_filename, sheetname,
                               index_col=None, na_values=['NA'])
    # Remove any rows that aren't gene names
    # (aka the last line with pathway length).
    minspan_df = minspan_df.drop(minspan_df.index[-1])

    # Change (e) at end of reaction to _e
    #minspan_df.replace('\(e\)$', '_e')
    #for reacton in minspan_df[0]:
    #    reaction = re.sub('\(e\)$', '_e', reaction)

    # One minspan per column, starts on 3rd column.
    for col in minspan_df.columns[2:]:
        # Index of rows with values (not NaN)
        value_index = np.isnan(minspan_df.loc[:, col])
        value_index = np.logical_not(value_index)
        now_minspan = minspan_df.loc[:, 'Rxn'][value_index]
        now_minspan = now_minspan.values  # No labels
        now_minspan = now_minspan.tolist()  # No numpy arrays.
        # Change from unicode... not sure if this is wise.
        now_minspan = [now_minspan[i].encode('ascii', 'ignore')
                       for i in range(len(now_minspan))]
        # Change (e) to _e
        now_minspan = [reaction.replace('(e)', '_e') for reaction in
                       now_minspan]
        now_minspan = [reaction.replace('-', '__') for reaction in
                       now_minspan]

        minspan_list.append(now_minspan)

    return minspan_list


def make_gene_fpkm_dict(fpkm_filename, norm=False, lowcut=0, highcut=None):
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
    if norm:  # Normalize values
        fpkm_df.FPKM = fpkm_df.FPKM / max(fpkm_df.FPKM)

    # Convert to dict of ids and FPKM values.
    gene_fpkm_dict = dict(zip(fpkm_df.tracking_id, fpkm_df.FPKM))
    return gene_fpkm_dict


def make_reaction_fpkm_dict(model, gene_fpkm_dict, mode=0):
    '''Create a dict with reactions as keys and fpkms as values.
    INPUT
        model              Cobra model
        gene_fpkm_dict     Dictionary mapping fpkm values to
                           genes. The output of make_gene_fpkm_dict
        mode               Indictates how to combine fpkms according to
                           gene reaction association.
                           0 averages and's and sums or's
                           1 sums everything
    OUTPUT
        reaction_fpkm_dict
    '''
    reaction_fpkm_dict = {}

    # Averages and's and sums or's
    if mode == 0:
        for reaction in model.reactions:
            gene_bool = reaction.gene_reaction_rule
            gene_bool = re.sub('\(|\)', ' ', gene_bool)
            gene_bool = gene_bool.split('or')
            fpkm_combined = 0
            data_count = 0
            for orgroup in gene_bool:
                andgroup = orgroup.split('and')
                andgroup_fpkm = 0
                avgby = len(andgroup)
                for genename in andgroup:
                    genename = genename.strip()
                    try:
                        andgroup_fpkm += gene_fpkm_dict[genename]
                        data_count += 1
                    except:
                        avgby -= 1
                        #print "%s is not in reaction_fpkm_dict." \
                        #    % str(genename)
                try:
                    andgroup_fpkm = andgroup_fpkm / avgby
                    fpkm_combined += andgroup_fpkm
                except:
                    #print('No data for this andgroup')
                    pass
            if data_count == 0:
                reaction_fpkm_dict[reaction.id] = float('NaN')
                #print('no data for ', reaction.id)
            else:
                reaction_fpkm_dict[reaction.id] = fpkm_combined
    # This will just add them up
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

    # Should report how many reactions from model not in fpkm file
    # and vice versa
    return reaction_fpkm_dict


def make_minspan_fpkm_list(minspan_list, reaction_fpkm_dict,
                           sort=False):
    '''Rank minspans based on the average fpkm values of their
    constituite reactionsk.
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
        no_data_count = 0  # Don't average in no data.
        for reaction in minspan:
            try:
                reaction_val = reaction_fpkm_dict[reaction]
                if math.isnan(reaction_val):
                    no_data_count += 1
                else:
                    minspan_fpkm += reaction_val
            except:
                # A lot don't match becaues of SBML naming vs names
                # from the minspan file.
                no_data_count += 1
                #print "%s is not in reaction_fpkm_dict." \
                #% str(reaction)
        try:
            # Just average for now.
            minspan_fpkm = minspan_fpkm / (len(minspan) - no_data_count)
        except:
            # There is no data for this minspan. Nan fucks up sorted()
            minspan_fpkm = -1
        minspan_fpkm_list.append(minspan_fpkm)

    # Create list of tuples with fpkm values and then the minspan list
    minspan_fpkm_list = zip(minspan_fpkm_list, minspan_list)
    if sort:
        minspan_fpkm_list = sorted(minspan_fpkm_list, reverse=True)

    return minspan_fpkm_list


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


def fpkm_comparison(fpkm_filename1, fpkm_filename2,
                    minspan_list, model, mode=0):
    '''Graphs a comparison betweet to fpkm data set on the same
    minspan ranking
    '''

    filename_list = [fpkm_filename1, fpkm_filename2]
    minspan_fpkm_vals = [0, 0]
    minspan_fpkm_ranks = [0, 0]
    x = range(1, len(minspan_list) + 1)  # Just a range for rank #'s.

    for i in range(2):
        gene_fpkm_dict = make_gene_fpkm_dict(filename_list[i])
        reaction_fpkm_dict = make_reaction_fpkm_dict(
                             model, gene_fpkm_dict, mode)
        minspan_fpkm_list = make_minspan_fpkm_list(
                            minspan_list, reaction_fpkm_dict, False)
        # Change (value, minspan) tuple structure to two lists and
        # pull out values.
        minspan_fpkm_vals[i] = zip(*minspan_fpkm_list)[0]

        # Find the rank of each minspan.
        rank = zip(minspan_fpkm_list, x)
        rank = sorted(rank, reverse=True)
        minspan_fpkm_ranks[i] = zip(*rank)[1]  # Just pull out ranks

    # Correlation
    m, b, r, p_val, std_err = stats.linregress(
                              minspan_fpkm_vals[0],
                              minspan_fpkm_vals[1])
    print('R^2 for fpkm values is: ', r)
    m, b, r, p_val, std_err = stats.linregress(
                              minspan_fpkm_ranks[0],
                              minspan_fpkm_ranks[1])
    print('R^2 for rank order is: ', r)

    # Plotting
    plt.subplots(2, 2, figsize=(10, 10), dpi=100)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('title')

    # Value comparison.
    plt.subplot(2, 2, 1)
    plt.title('Minspan number (unranked) vs. fpkm score.')
    plt.xlabel('rank')
    plt.ylabel('fpkm')
    #plt.legend(["data1", "data2"])
    plt.scatter(x, minspan_fpkm_vals[0], c='b')
    plt.scatter(x, minspan_fpkm_vals[1], c='r')
    plt.subplot(2, 2, 2)
    plt.title('fpkm vs. fpkm')
    plt.xlabel('fpkm 1')
    plt.ylabel('fpkm 2')
    plt.scatter(minspan_fpkm_vals[0], minspan_fpkm_vals[1], c='g')

    # Rank comparison
    plt.subplot(2, 2, 3)
    plt.title('Minspan number (unranked) vs. rank.')
    plt.xlabel('rank')
    plt.ylabel('rank')
    #plt.legend(["data1", "data2"])
    plt.scatter(x, minspan_fpkm_ranks[0], c='b')
    plt.scatter(x, minspan_fpkm_ranks[1], c='r')
    plt.subplot(2, 2, 4)
    plt.title('rank vs rank')
    plt.xlabel('minspan rank')
    plt.ylabel('minspan rank')
    plt.scatter(minspan_fpkm_ranks[0], minspan_fpkm_ranks[1], c='y')
    plt.show()  # I'm not sure why fig is not needen as an arg here.


def Eflux(model, reaction_fpkm_dict, pFBA=False,
          scale=True, scale_reaction_id='EX_glc_e', scale_value=-10):
    '''Optimizes the model with bounds adjusted by fpkm data. If a
    scale reaction is provided, bounds will only be dictated by fpkm
    values and afterwards the fluxes will be weighted by the flux ratio
    of the scale reaction without and with Eflux. e.g. Provide 'EX_glc_e'
    as scale reaction with flux under normal condidions of -10. Eflux
    will be performed with bounds as dictatedby fpkm levels.
    Eflux_sol.x will then be adjusted such that the flux of 'EX_glc_e'
    is again -10. This makes it comparible to the normal solution, but
    you lose specific control over the bounds.
    INPUT
        model
        reaction_fpkm_dict  From make_reaction_fpkm_dict
        scale               Boolean indicates to scale fluxes or not
        scale_reaction_id   Name of reaction from which fluxes are scaled
    OUTPUT
        sol                 Normal Cobra solution with pFBA
        Eflux_sol           Cobra solution under Eflux and pFBA
    '''

    # Adjust bounds
    k = list(reaction_fpkm_dict.keys())
    v = list(reaction_fpkm_dict.values())
    v = abs(v / max(v))  # abs should be superfulous
    fpkm_array = zip(k, v)

    # Eflux ignores human set uptake bounds
    for reaction in model.reactions:
        if reaction.lower_bound < 0:
            reaction.lower_bound = -1000
        else:
            reaction.lower_bound = 0
        if reaction.upper_bound > 0:
            reaction.upper_bound = 1000
        else:
            reaction.upper_bound = 0

    # Here we scale each bound by the relative fpkm values
    for rxn in fpkm_array:
        if not np.isnan(rxn[1]):
            rxn_name = rxn[0]
            cobrarxn = model.reactions.get_by_id(rxn_name)
            cobrarxn.lower_bound = cobrarxn.lower_bound * rxn[1]
            cobrarxn.upper_bound = cobrarxn.upper_bound * rxn[1]
        else:
            pass

    # Compute flux
    biomass_rxn = 'Ec_biomass_iJO1366_core_53p95M'
    model.reactions.get_by_id(biomass_rxn).lower_bound = 0
    model.reactions.get_by_id(biomass_rxn).upper_bound = 1000

    # pFBA does not work
    if pFBA:
        # This should find the objective function automatically
        pfba_sol = pfba_example.run_pfba(model,
                                         biomass_rxn)
    else:
        model.optimize()

    # Weight solution and bounds by scaling reaction
    if scale:
        if not pFBA:
            scale_flux = model.solution.x_dict[scale_reaction_id]
        else:
            scale_flux = pfba_sol[scale_reaction_id]

        #print scale_value, scale_flux

        if scale_flux == 0:
            print 'Scale flux is 0'
            return pfba_sol

        scale_by = abs(scale_value / scale_flux)

    model.solution.x = list(np.array(model.solution.x) * scale_by)
    model.solution.f = model.solution.f * scale_by

    for rxn in model.reactions:
        rxn.lower_bound = rxn.lower_bound * scale_by
        rxn.upper_bound = rxn.upper_bound * scale_by

        model.solution.x_dict[rxn.id] = \
                model.solution.x_dict[rxn.id] * scale_by

        #for k, b in pfba_sol.iteritems():
        #    pfba_sol[k] = v * scale_by

    if pFBA:
        Eflux_sol = pfba_sol
    else:
        Eflux_sol = model.solution

    return Eflux_sol


def relax_bounds(model):
    # Test the impact of Eflux bounds
    unbound_reactions = []
    count = 0

    relaxed_value = \
        model.reactions.Ec_biomass_iJO1366_core_53p95M.upper_bound

    # Eflux orginal solution
    Eflux_sol = model.solution

    for rxn in model.reactions:
        # Save current bounds.
        ub_Eflux = model.reactions.get_by_id(rxn.id).upper_bound
        lb_Eflux = model.reactions.get_by_id(rxn.id).lower_bound

        # Change bounds back to that of original model
        model.reactions.get_by_id(rxn.id).upper_bound = relaxed_value
        model.reactions.get_by_id(rxn.id).lower_bound = -relaxed_value

        # Calculate new flux
        model.optimize()

        # Save data for when f was different
        if abs(model.solution.f - Eflux_sol.f) > 0.001:
            unbound_reactions.append((rxn.id,
                                      model.solution.f - Eflux_sol.f))

        # Revert bounds to normal for next iteration
        model.reactions.get_by_id(rxn.id).upper_bound = ub_Eflux
        model.reactions.get_by_id(rxn.id).lower_bound = lb_Eflux

        count += 1
        if count % 100 == 0:
            print count

    return


def make_minspan_k_dict(minspan_list):
    '''
    INPUT
        minspan_list
    OUTPUT
        minspan_k_dict    minspan id as key, list of reactions as value
    '''
    minspan_ids = range(1, len(minspan_list) + 1)
    minspan_k_dict = {}
    for k in minspan_ids:
    # The keys will one day be ints, so I will start it that way.
        minspan_k_dict[str(k)] = minspan_list[k - 1]

    return minspan_k_dict


def make_reaction_k_dict(model, minspan_k_dict, report_wo=False):
    '''
    INPUTS
        model
        minspan_k_dict
        report_wo    Return reaction_wo_minspans or not.
    OUTPUTS
        reaction_k_dict    reaction name as k, with minspans involved
                           as values.
        reaction_wo_minspans    Reactions in the model that are not in
                                a minspan.
    '''
    reaction_k_dict = {}
    reaction_wo_minspans = []
    for reaction in model.reactions:
        temp_reaction_list = []
        for k, v in minspan_k_dict.iteritems():
            if reaction.id in v:
                temp_reaction_list.append(k)
        if len(temp_reaction_list) == 0:
            reaction_wo_minspans.append(reaction.id)
        reaction_k_dict[reaction.id] = temp_reaction_list

    if report_wo:
        return reaction_k_dict, reaction_wo_minspans
    else:
        return reaction_k_dict


def make_reaction_k_dict1(minspan_k_dict):
    ''' This makes a reaction dict but uses reactions that are only 
    in the minspans
    INPUTS
        minspan_k_dict
    OUTPUTS
        reaction_k_dict    reaction name as k, with minspans involved
                           as values.
    '''
    reaction_k_dict = {}

reaction_list = []
for minspan in minspan_list:
    for rxn in minspan:
        if rxn not in reaction_list:
            reaction_list.append(rxn)
        else:
            pass
    for reaction in model.reactions:
        temp_reaction_list = []
        for k, v in minspan_k_dict.iteritems():
            if reaction.id in v:
                temp_reaction_list.append(k)
        if len(temp_reaction_list) == 0:
            reaction_wo_minspans.append(reaction.id)
        reaction_k_dict[reaction.id] = temp_reaction_list

    if report_wo:
        return reaction_k_dict, reaction_wo_minspans
    else:
        return reaction_k_dict
