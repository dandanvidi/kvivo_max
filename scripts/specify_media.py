from cPickle import load
from cobra.manipulation import initialize_growth_medium
from cobra.test import salmonella_pickle
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
import pickle
import pandas as pd

def perform_pFBA(model, media, gr):
    
    initialize_growth_medium(model, media)
    rxns = dict([(r.id, r) for r in model.reactions])
    rxns['Ec_biomass_iJO1366_core_53p95M'].upper_bound = gr            
    print "solving pFBA",
    optimize_minimal_flux(model, already_irreversible=True)
    print "-> %.2f" %model.solution.f
    flux_dist = pd.DataFrame(model.solution.x_dict.items()).set_index(0)
    
    return flux_dist    


media_dict = pickle.load( open( "../data/media.p", "rb" ) )
gc = pd.DataFrame.from_csv('../data/growth_conditions.csv')
model = create_cobra_model_from_sbml_file('../data/iJO1366.xml')
convert_to_irreversible(model)

fluxes = pd.DataFrame(index=map(lambda x: x.id, model.reactions), columns=gc.index)
for c in gc.iterrows():
    gr = c[1]['growth rate [h-1]']
    try:
        media = media_dict[c[0]]
    except KeyError:
        raise Exception("media composition for %s not specified" %c[1]['media'])
    fluxes[c[0]] = perform_pFBA(model, media, gr)
    
fluxes.index.name = 'reaction'
fluxes.to_csv('../data/flux[mmol_gCDW_h].csv')
