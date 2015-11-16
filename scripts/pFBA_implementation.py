from copy import deepcopy
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from cobra.manipulation.modify import convert_to_irreversible
from cobra.io.sbml import create_cobra_model_from_sbml_file
import pandas as pd

def perform_pFBA(model, cs, gr, ur):

    model = deepcopy(model)
    convert_to_irreversible(model)            

    rxns = dict([(r.id, r) for r in model.reactions])
    rxns['EX_glc_e'].lower_bound = 0 # uptake of carbon source reaction is initialized    
    try:
        rxns['EX_' + cs + '_e'].lower_bound = -ur # redefine sole carbon source uptake reaction in mmol/gr/h
    except:
        print cs, ur
        rxns['EX_glc_e'].lower_bound = -ur
    rxns['Ec_biomass_iJO1366_core_53p95M'].upper_bound = gr            
    print "solving pFBA",
    optimize_minimal_flux(model, already_irreversible=True)
    
    flux_dist = pd.DataFrame(model.solution.x_dict.items()).set_index(0)
    
    return flux_dist   
    
if __name__ == "__main__":
    model_fname = "../data/iJO1366.xml"
    model = create_cobra_model_from_sbml_file(model_fname)
    convert_to_irreversible(model)
    reactions = map(lambda x: x.id, model.reactions)
    fluxes = perform_pFBA(model, 'glc', 0.5, 18.5)