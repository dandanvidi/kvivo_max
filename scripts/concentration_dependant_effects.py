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
        self.km_data = pd.DataFrame.from_csv('../data/km_values.csv', sep='\t')
        self.km_data.km_sparse = map(ast.literal_eval, self.km_data.km_sparse)
        self.metab_conc = pd.DataFrame.from_csv('../data/bennett_metabolite_concentrations[mM].csv') # concentrations in mM
        self.metab_conc.dropna(how='all', inplace=True)
        self.measured_cids = set(self.metab_conc.index) & set(self.Kmodel.cids)
        self.conc = self._set_metabolite_concentrations() # in Molar
        self.cid_to_conc = dict(zip(self.Kmodel.cids, self.conc[0]))
        self.get_km_values()
        self.dG_prime = self.get_dG_prime()        
        self.backward_flux = self.backward_flux_effects().T
        self.under_saturation = self.substrate_saturation_effect().T
#        self.substrates = self.get_reaction_substrates()
#        self.Ks = self.get_known_Ks_values()

    def _set_metabolite_concentrations(self):        
        #initialize concentrations of metabolites as 100uM
        conc = np.ones((1, len(self.Kmodel.cids))) * 1e-4 # concentrations in M
        for cid  in self.measured_cids: # set the conc of known metabolites
            i = self.Kmodel.cids.index(cid)
            # set the conc of water as 1M
            if cid == 'C00001':
                conc[0, i] = 1
            else:
                c = self.metab_conc[self.condition][cid]
                if not np.isnan(c):
                    conc[0, i] = c * 1e-3  # concentrations in M
        return conc

    def get_dG_prime(self):
        dG_prime = self.dG0_prime + self.R * self.T * np.dot(np.log(self.conc), 
                                                                     self.Kmodel.S)                
        return dG_prime

    def backward_flux_effects(self):
        tmp = self.dG_prime
        tmp = np.clip(tmp, -200, 200, out=tmp)
        t = -unumpy.expm1(tmp / (self.R * self.T))
        for i,r in enumerate(self.reactions):
            if not unumpy.isnan(t)[0,i]:
                r.backward = t[0,i]
            else:
                r.backward = np.nan
        return t    

    def get_km_values(self):
        for r in self.reactions:
            if r.id in self.km_data.index:
                sparse = self.km_data.loc[r.id].values[0]
                sparse = {k:unumpy.uarray(np.mean(km), np.std(km)) / 1e6
                            for k, km in sparse.iteritems()} 
            else:
                sparse = None
            r.km_sparse = sparse # KM values are in units of Molar (same as metabolite concentrations)

    def _get_prod_s_by_ks(self):
        prod = []
        for r in self.reactions:
            if r.km_sparse:
                substrates = [k for k,v in r.kegg_reaction.iteritems() if v<0]
                ks = {k:v for k,v in r.km_sparse.iteritems() if k in substrates}
                                
                #check if all ks values are known, else assign nan to the product
                tmp = np.array([v.n for v in ks.itervalues()])
                if np.isnan(tmp).any():
                    prod.append(unumpy.uarray(np.nan, np.nan))                    
                    continue
                #check if the concentration of all substrates is known
                elif not set(substrates).issubset(self.measured_cids):
                    prod.append(unumpy.uarray(np.nan, np.nan))                    
                    continue
                
                else:
                    conc = {k:v for k,v in self.cid_to_conc.iteritems() if k in substrates}
                    s_by_ks = {k:conc[k]/ks[k] for k in ks if k in conc}
                    prod.append(np.prod(s_by_ks.values()))
            else:
                prod.append(unumpy.uarray(np.nan, np.nan))                    
        return np.matrix(prod)

    def substrate_saturation_effect(self):
        prod = self._get_prod_s_by_ks()
        s = prod / (1+prod)
        for i,r in enumerate(self.reactions):
            if not unumpy.isnan(s)[0,i]:
                r.saturation = s[0,i]
            else:
                r.saturation = np.nan
        return s
        
    def get_reactions_with_all_known_KM(self):
        reaction_list = []        
        for r in self.reactions:
            if r.km_sparse:
                if not any(unumpy.isnan(r.km_sparse.values())):
                    reaction_list.append(r)
        return reaction_list

    def get_reactions_with_all_known_S(self):
        reaction_list = []        
        for r in self.reactions:
            cids = set(r.kegg_reaction.keys())
            if cids.issubset(set(self.measured_cids)):
                reaction_list.append(r)
        return reaction_list
        
if __name__ == "__main__":
    
    from catalytic_rates import rates
    R = rates()    
    
    index = R.kcat.index & R.kmax.index
    reactions = [R.rxns[r] for r in index]
    for gc in ['glc']:#, 'ac', 'glyc']:                
        out = pd.DataFrame(index=index, columns=['under saturation', 'backward flux'])
        mm = MM_kinetics(reactions, gc)
        for r in mm.reactions:
            if not unumpy.isnan(r.saturation) and not unumpy.isnan(r.backward):
                out['under saturation'][r.id] = r.saturation.n
                out['backward flux'][r.id] = r.backward.n
        out.dropna(inplace=True)
#        out.to_csv("../res/conc_dependant_effects_on_%s.csv" %gc, sep='\t')
