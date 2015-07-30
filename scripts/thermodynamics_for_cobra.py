import warnings
import pandas as pd
import numpy as np
import uncertainties.unumpy as unumpy  
from cobra.io.sbml import create_cobra_model_from_sbml_file
from copy import deepcopy
from component_contribution.kegg_reaction import KeggReaction
from component_contribution.kegg_model import KeggModel
from component_contribution.component_contribution import ComponentContribution
from component_contribution.thermodynamic_constants import R, default_T

class THERMODYNAMICS_FOR_COBRA(object):

    def __init__(self, reactions):

        self.cc = ComponentContribution.init()
        
        self.pH = 7.5
        self.I = 0.2
        self.T = default_T

        self._not_balanced = []        
        self.Kmodel = self.generate_kegg_model()
        self.get_udG0_prime()
        self.get_log_Keq()
        self.get_log_RI()
        
    def generate_kegg_model(self):
        
        rstrings = []
        for r in reactions:
            k = r.kegg_reaction
            if k:
                if k.is_balanced() and not k.is_empty():
                    rstrings.append(k.write_formula())
            else:
                self._not_balanced.append(r)
        return KeggModel.from_formulas(rstrings)
        
    def get_udG0_prime(self):
        '''
            Calculates the dG0 of a list of a reaction.
            Uses the component-contribution package (Noor et al) to estimate
            the standard Gibbs Free Energy of reactions based on 
            component contribution  approach and measured values (NIST and Alberty)
        '''
        
        self.Kmodel.add_thermo(self.cc)
        dG0_prime, dG0_cov = self.Kmodel.get_transformed_dG0(pH=7.5, I=0.2, T=298.15)
        dG0_std = 1.96*np.diag(dG0_cov.round(1))
        arr = unumpy.uarray(dG0_prime.flat, dG0_std.flat)
        for i, r in enumerate(reactions):
            if r in self._not_balanced:                
                r.dG0_prime = np.NaN
            else:
                r.dG0_prime = arr[i]

    def get_log_Keq(self):
        '''
            Calculates the equilibrium constants of a reaction, using dG0_prime.
        '''
        for r in reactions:
            # cap the maximal dG0 at 200
            tmp = r.dG0_prime
            if tmp > 200: tmp = 200
            if tmp < -200: tmp = -200
            r.logKeq = -tmp / (R*default_T)
            
    def get_log_RI(self, fixed_conc=0.1):
        '''
            Calculates the reversibility index (RI) of a reaction.
            The RI represent the change in concentrations of metabolites
            (from equal reaction reactants) that will make the reaction reversible.
            That is, the higher RI is, the more irreversible the reaction.
            A convenient threshold for reversibility is RI>=1000, that is a change of
            1000% in metabolite concentrations is required in ordeer to flip the
            reaction direction. 
        '''
        for r in reactions:   
            N_P = sum([v  for v in r.metabolites.itervalues() if v > 0])
            N_S = sum([-v for v in r.metabolites.itervalues() if v < 0])
            N = N_P + N_S
            r.logRI = (2/N) * (r.logKeq + (N_P - N_S)*np.log(fixed_conc))
            
if __name__ == "__main__":
    
    from model_addons import add_to_model
    model_fname = "../data/iJO1366.xml"
    model = create_cobra_model_from_sbml_file(model_fname)
    add_to_model(model)
    
    reactions = ['MDH','FBA','TPI','FBP','PGM','SERAT','TMDS','DBTS','DM_4CRSOL']
    reactions = map(model.reactions.get_by_id, reactions)    
    
    TFC = THERMODYNAMICS_FOR_COBRA(reactions)

