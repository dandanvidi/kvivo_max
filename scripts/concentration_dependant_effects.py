import ast
import uncertainties.unumpy as unumpy  
import pandas as pd
from cobra.io.sbml import create_cobra_model_from_sbml_file
import numpy as np
from thermodynamics_for_cobra import reaction_thermodynamics

class MM_kinetics(reaction_thermodynamics):

    def __init__(self, reactions, condition):
        
        reaction_thermodynamics.__init__(self, reactions)
        
        self.condition = condition        
        self.reactions = reactions        
        self.km_data = pd.DataFrame.from_csv('../data/km_values.csv', sep='\t')
        self.km_data.km_sparse = map(ast.literal_eval, self.km_data.km_sparse)
        self.metab_conc = pd.DataFrame.from_csv('../data/bennett_metabolite_concentrations[mM].csv') # concentrations in mM
        self.metab_conc.dropna(how='all', inplace=True)
        self.measured_cids = set(self.metab_conc.index) & set(self.Kmodel.cids)
        self.conc = self._set_metabolite_concentrations()
        self.cid_to_conc = dict(zip(self.Kmodel.cids, self.conc[0]))
        self.get_km_values()
        self.dG_prime = self.get_dG_prime()        
        self.backward_flux = self.backward_flux_effects()
        self.under_saturation = self.substrate_saturation_effect()
#        self.substrates = self.get_reaction_substrates()
#        self.Ks = self.get_known_Ks_values()

    def _set_metabolite_concentrations(self):        
        #initialize concentrations of metabolites as 100uM
        conc = np.ones((1, len(self.Kmodel.cids))) * 1e-4 # concentrations in M
        for cid  in self.measured_cids: # set the conc of known metabolites
            if cid in self.Kmodel.cids:
                i = self.Kmodel.cids.index(cid)
                # set the conc of water as 1M
                if cid == 'C00001':
                    conc[0, i] = 1
                else:
                    c = self.metab_conc[self.condition][cid]
                    if not np.isnan(c):
                        conc[0, i] = c  * 1e-3  # concentrations in M
        return conc

    def get_dG_prime(self):
        dG_prime = self.dG0_prime + self.R * self.T * np.dot(np.log(self.conc), 
                                                                     self.Kmodel.S)                
        return dG_prime

    def backward_flux_effects(self):
        tmp = self.dG_prime
#        tmp[tmp>0] = 0
        return -unumpy.expm1(tmp / (self.R * self.T))

    def get_km_values(self):
        for r in self.reactions:
            if r.id in self.km_data.index:
                sparse = self.km_data.loc[r.id].values[0]
                sparse = {k:unumpy.uarray(np.mean(km), np.std(km)) / 1000
                            for k, km in sparse.iteritems()} 
            else:
                sparse = None
            r.km_sparse = sparse # KM values are in mM

    def _get_prod_s_by_ks(self, use_default_conc=False):
        prod = []
        for r in self.reactions:
            try:
                substrates = [k for k,v in r.kegg_reaction.iteritems() if v<0]
                ks = {k:v for k,v in r.km_sparse.iteritems() if k in substrates}
                if use_default_conc:
                    conc = {k:v for k,v in self.cid_to_conc.iteritems() 
                                            if k in ks}                
                else:
                    conc = {k:v for k,v in self.cid_to_conc.iteritems() 
                                            if k in ks and k in self.measured_cids}
                s_by_ks = {k:conc[k]/ks[k] for k in ks if k in conc}
                prod.append(np.prod(s_by_ks.values()))
            except AttributeError:
                prod.append(unumpy.uarray(np.nan, np.nan))
        return np.matrix(prod) 

    def substrate_saturation_effect(self):
        prod = self._get_prod_s_by_ks()
        return prod / (1+prod)
    
if __name__ == "__main__":
    
    from model_addons import add_to_model
    model_fname = "../data/iJO1366.xml"
    model = create_cobra_model_from_sbml_file(model_fname)
    add_to_model(model)
    
    reactions = ['MDH','AIRC2','TPI','FBP','PGM','SERAT','TMDS','DBTS','DM_4CRSOL']
    reactions = map(model.reactions.get_by_id, reactions)    
    mm = MM_kinetics(reactions, 'glc')