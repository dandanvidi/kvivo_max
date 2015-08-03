import numpy as np
import uncertainties.unumpy as unumpy  
from cobra.io.sbml import create_cobra_model_from_sbml_file
from component_contribution.kegg_model import KeggModel
from component_contribution.component_contribution import ComponentContribution
from component_contribution.thermodynamic_constants import R, default_T

class reaction_thermodynamics(object):

    def __init__(self, reactions):

        self.cc = ComponentContribution.init()
        
        self.pH = 7.5
        self.I = 0.2
        self.T = default_T
        self.R = R
        self._not_balanced = []        
        self.Kmodel = self.generate_kegg_model(reactions)
        self.dG0_prime = self.add_thermodynamics(reactions)
        
    def generate_kegg_model(self, reactions):
        
        rstrings = []
        for r in reactions:
            k = r.kegg_reaction
            if k:
                if k.is_balanced() and not k.is_empty():
                    rstrings.append(k.write_formula())
            else:
                self._not_balanced.append(r)
        return KeggModel.from_formulas(rstrings)
        
    def add_thermodynamics(self, reactions):
        '''
            Calculates the dG0 of a list of a reaction.
            Uses the component-contribution package (Noor et al) to estimate
            the standard Gibbs Free Energy of reactions based on 
            component contribution  approach and measured values (NIST and Alberty)

            Calculates the reversibility index (RI) of a reaction.
            The RI represent the change in concentrations of metabolites
            (from equal reaction reactants) that will make the reaction reversible.
            That is, the higher RI is, the more irreversible the reaction.
            A convenient threshold for reversibility is RI>=1000, that is a change of
            1000% in metabolite concentrations is required in ordeer to flip the
            reaction direction. 
        '''
        
        self.Kmodel.add_thermo(self.cc)
        dG0_prime, dG0_cov = self.Kmodel.get_transformed_dG0(pH=7.5, I=0.2, T=298.15)
        dG0_std = 1.96*np.diag(dG0_cov.round(1))
        dG0_prime = unumpy.uarray(dG0_prime.flat, dG0_std.flat)

        conc = np.ones((1, len(self.Kmodel.cids))) * 1e-3 # concentrations in M
        dGm_prime = dG0_prime + self.R * self.T * np.dot(np.log(conc), self.Kmodel.S)
        
        for i, r in enumerate(reactions):
            if r in self._not_balanced:                
                r.dG0_prime = np.NaN
            else:
                try:
                    # d0_prime
                    r.dG0_prime = dG0_prime[i]
                    # dGm_prime
                    r.dGm_prime = dGm_prime[0,i]
                    
                except IndexError:
                    r.dG0_prime = np.NaN
                    r.dGm_prime = np.NaN
                # equilibrium constant
                tmp = r.dG0_prime
                if tmp > 200: tmp = 200
                if tmp < -200: tmp = -200
                r.logKeq = -tmp / (self.R * self.T)
                # reversibility index
                N_P = sum([v  for v in r.metabolites.itervalues() if v > 0])
                N_S = sum([-v for v in r.metabolites.itervalues() if v < 0])
                N = N_P + N_S
                fixed_conc=0.1                
                r.logRI = (2/N) * (r.logKeq + (N_P - N_S)*np.log(fixed_conc))        

        return dG0_prime
            
if __name__ == "__main__":
    
    from model_addons import add_to_model
    model_fname = "../data/iJO1366.xml"
    model = create_cobra_model_from_sbml_file(model_fname)
    add_to_model(model)
    
    reactions = ['MDH','FBA','TPI','FBP','PGM','SERAT','TMDS','DBTS','DM_4CRSOL']
    reactions = map(model.reactions.get_by_id, reactions)    
#    reactions = model.reactions
    Th = reaction_thermodynamics(reactions)

