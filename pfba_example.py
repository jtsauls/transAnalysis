from cobra import Reaction

def make_model_irrev(model):
    irrev_model=model.copy()
    rxns_to_add=[]
    reaction_to_metabolite={}
    old_to_new_rxn_map={}
    for rxn in irrev_model.reactions:
        if rxn.lower_bound<0:
            new_rxn=Reaction(name=rxn.id+'_reverse')
            new_rxn.id=new_rxn.name
            new_rxn._metabolites=dict([(met,-val) for met, val in rxn._metabolites.iteritems()])
            new_rxn.lower_bound=0
            new_rxn.upper_bound=-1*(rxn.lower_bound)
            rxns_to_add.append(new_rxn)

            rxn.lower_bound=0
            if rxn.upper_bound<0:
                new_rxn.lower_bound=-1*(rxn.upper_bound)
                rxn.upper_bound=0
            old_to_new_rxn_map[rxn.id]={rxn.id:1,rxn.id+'_reverse':-1}
        else:
            old_to_new_rxn_map[rxn.id]={rxn.id:1}
    irrev_model.add_reactions(rxns_to_add)
    return irrev_model,old_to_new_rxn_map

def run_pfba(model,biomass_function_name):
    model.optimize(new_objective={biomass_function_name:1})

    model.reactions.get_by_id(biomass_function_name).upper_bound=model.solution.x_dict[biomass_function_name]
    model.reactions.get_by_id(biomass_function_name).lower_bound=model.solution.x_dict[biomass_function_name]
    
    irrev_model,old_to_new_rxn_map=make_model_irrev(model)
    obj=dict([(k.id,1) for k in irrev_model.reactions])
    irrev_model.optimize(new_objective=obj,objective_sense='minimize')
    
    sol={}
    for rxn in model.reactions:
        sol[rxn.id]=sum([irrev_model.solution.x_dict[nrxn]*stoich for nrxn,stoich in old_to_new_rxn_map[rxn.id].items()])
    
    return sol
