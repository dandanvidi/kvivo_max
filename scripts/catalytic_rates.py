import pandas as pd
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import numpy as np
from model_addons import add_to_model

class rates(object):
    
    def __init__(self):
        
        self.model = create_cobra_model_from_sbml_file("../data/iJO1366.xml")

        # Modify model
        convert_to_irreversible(self.model)  
        self.rxns = dict([(r.id, r) for r in self.model.reactions])
        self.genes = dict([(g.id, g) for g in self.model.genes])
        add_to_model(self.model)
        self.include_specific_isozmyes()

        self.gc = pd.DataFrame.from_csv("../data/growth_conditions.csv")

        self.flux_data = pd.DataFrame.from_csv('../data/flux[mmol_gCDW_s].csv')
        # PPKr_reverse reaction is used for ATP generation from ADP 
        # in the FBA model. Nevertheless, acording to EcoCyc, it is used to 
        # to generate polyP (inorganic phosphate) chains from ATP and it is not
        # part of the oxidative phosphorilation, thus removed from rate calculations
        if 'PPKr_reverse' in self.flux_data.index:            
            self.flux_data.drop('PPKr_reverse', axis=0, inplace=True)
            
        self.expression_data = pd.DataFrame.from_csv('../data/abundance[mmol_gCDW].csv')

        self.enzymatic_reactions = self.enzymatic_reactions()       
        self.homomeric_reactions = self.reactions_by_homomeric_enzymes() 
        self.kapp = self.get_kapp() # per subunit
        self.SA = self.get_specific_activity()

        self.kcat = pd.DataFrame.from_csv("../data/kcat_data.csv") 
        self.p_per_as = (self.kcat['polypeptides per complex'] 
                                    / self.kcat['catalytic sites per complex'])
        self.kmax = self.get_kmax(self.kapp)
        self.SAmax = self.get_maximum_specific_activity(self.SA)             

    def enzymatic_reactions(self):
        '''
            Returns a list of cobra Reaction objects catalyzed by enzymes.
        '''
        return filter(lambda r:len(r.genes)>=1, self.model.reactions)

    def include_specific_isozmyes(self):
        '''
            Possible add-ons to the list of unique homomeric enzymes
            obtained by the function "reactions_to_unique_enzyme". 
            These isoenzymes are known to have only one active isoenzyme
            across all tested conditions and therefore were manually added
        '''
        pairs = [
                ('METS','b3829'),# metE - cobalamin-independent homocysteine transmethylase
                ('HCO3E','b0126'),# can - carbonic anhydrase
                ('PFK','b3916'), #  6-phosphofructokinase
                ('RPI','b2914') #  ribose-5-phosphate isomerase A
                ] 
        for (r,g) in pairs:
            
            self.rxns[r]._genes = [self.model.genes.get_by_id(g)]
            self.rxns[r].gene_reaction_rule = '('+g+')'
            
    def reactions_by_unique_enzyme(self):
        '''
            Returns a list of reactions (as cobra REACTION objects)
            in the model catalyzed by unique enzymes. Enzymes can either
            be homomeric or hetorometic complexes.
        '''        
        one_enzyme_reac = filter(lambda r: 'or' not in r.gene_reaction_rule, 
                                 self.enzymatic_reactions)
        return one_enzyme_reac
#        


    def reactions_by_homomeric_enzymes(self):
        '''
            Returns a list of reactions (as cobra REACTION objects)
            in the model catalyzed by unique enzymes which are composed
            of a single polypeptide chain, i.e., unique homomeric enzymes.
        '''
        homomers = filter(lambda r: len(r.genes)==1, 
                                 self.enzymatic_reactions)
        return homomers

    def _convert_copies_fL_to_mmol_gCDW(self, expression_data):
        '''
            Convertes the units of proteomics data (usually reported in 
            copies per fL of cytoplasm) to units of mmol per gCDW.
            This unit conversion is performed to match flux units from
            metabolic models (usually given in mmol/gCDW/h)
        '''
        rho = 1100 # average cell density gr/liter
        DW_fraction = 0.3 # fraction of DW of cells
        Avogadro = 6.02214129 # Avogadro's number "exponent-less"
        expression_data[expression_data<10] = np.nan
        expression_data /= (Avogadro*1e5)
        expression_data /= (rho * DW_fraction)
        expression_data.to_csv('../cache/abundance[mmol_gCDW].csv')
        return expression_data

    def _convert_mmol_gCDW_h_to_mmol_gCDW_s(self, flux_data):
        '''
            Convertes the units of flux data (usually reported in 
            mmol/gCDW/h) to units of mmol/gCDW per second.
            This unit conversion is performed to allow calculation of
            turnover rates in units of s^-1, as traditioanlly excepted. 
        '''
        flux_data /= 3600        
        flux_data.to_csv('../cache/flux[mmol_gCDW_s].csv')
        return flux_data
        
    def _convert_mmol_gCDW_to_mg_gCDW(self, expression_data):
        
        genes = set(self.genes.keys()) & (set(expression_data.index))
        mass = [self.genes[g].MW for g in genes]        
        MW = pd.Series(index=genes, data=mass)
        return expression_data.loc[MW.index].mul(MW, axis=0) 
        
    def get_kapp(self):
        '''
            Calculates the catalytic rate of a single subunit of a homomeric
            enzyme for a given reaction, by dividing the flux through the 
            reaction by the abundance of the polypeptide chain that comprises 
            the enzyme.
            
            Arguments:
                flux [mmol/gCDW/s]
                proteomics [mmol/gCDW]
            
            Returns:
                pandas dataframe with catalytic rates per polypeptide chain
                in units of s^-1. Rows are reactions, columns are conditions 
        '''
        gc = self.flux_data.columns & self.expression_data.columns
        
        index = map(lambda x: x.id, self.homomeric_reactions)
        rate = pd.DataFrame(index=index, 
                            columns=gc)
        for r in self.reactions_by_homomeric_enzymes():
            try:
                rate.loc[r.id] = (self.flux_data.loc[r.id] / 
                                  self.expression_data.loc[list(r.genes)[0].id])
            except KeyError:
                continue
            
        rate.replace([0, np.inf, -np.inf], np.nan, inplace=True)
        rate.dropna(how='all', inplace=True)
        
        return rate 

    def get_kmax(self, kapp, minimal_conditions=5):
        '''
            Take the maximum rate of a given enzyme-reaction pair 
            across all conditions.
            
            Arguments:
                catalytic rate of enzyme-reaction pairs across conditions
                as a pandas dataframe. Rows are reactions, columns are conditions
            
            Returns:
                Maximal rate for each enzyme-reaction pair, the condition in 
                which it was found, the metabolic pathway associated with 
                the reaction and the carbon source on which the cells were grown.
                Notice that maximal rates are given per polypeptide chain and
                per active site in two seperate columns.
                Rate units are s^-1.
        '''
        kapp.dropna(thresh=minimal_conditions, inplace=True)
        kmax = pd.DataFrame(index=kapp.index)
        
        subsystems =[self.rxns[r].subsystem for r in kmax.index]
        genes = [list(self.rxns[r].genes)[0].id for r in kmax.index]
        names = [list(self.rxns[r].genes)[0].name for r in kmax.index]
        kmax.index.name = 'reaction'
        kmax['bnumber'] = genes
        kmax['primary gene name (uniprot)'] = names
        kmax['kmax per chain [s^-1]'] = kapp.max(axis=1)

        tmp = self.kapp.loc[kmax.index].mul(self.p_per_as[kmax.index], axis=0)
        kmax['kmax per active site [s^-1]'] = tmp.max(axis=1)
        
        kmax['subsystem'] = subsystems
        kmax['condition'] = kapp.idxmax(axis=1)
        
        return kmax
        
    def get_specific_activity(self):
        '''
            Calculates the specific activity in units of umol/mg/min
            for all reactions in the model. The sum of all associated
            polypeptide chains is used as the molecular weight of the enzyme
            and the flux through the reaction is divided by this weight.
            
            Notice that if a reaction can be carried by several different enzymes,
            i.e., isoenzymes, the returned values are a weighted average of the 
            rate of the enzymes by their mass.
            
            Arguments:
                flux [mmol/gCDW/s]
                proteomics [mmol/gCDW]

            Returns:
                pandas dataframe with specific activeites of enzymes
                in units of umol/mg/min. Rows are reactions, columns are conditions 
        '''
        weighted_mass = self._convert_mmol_gCDW_to_mg_gCDW(self.expression_data)
        reactions = map(lambda x: x.id, self.enzymatic_reactions)
        
        SA = pd.DataFrame(index=reactions, columns=self.gc.index)

        for r in self.enzymatic_reactions:
            genes = map(lambda x: x.id, r.genes)
            try:
                SA.loc[r.id] = self.flux_data.loc[r.id] / weighted_mass.loc[genes].sum()
            except KeyError:
                continue
        
        SA.replace([0, np.inf, -np.inf], np.nan, inplace=True)
        SA.dropna(how='all', inplace=True)
        
        return SA * 1000 * 60
        
    def get_maximum_specific_activity(self, specific_activity, minimal_conditions=5):

        '''
            Take the maximum rate of a given enzyme-reaction pair 
            across all conditions.
            
            Arguments:
                specific activities of enzyme-reaction pairs across conditions
                as a pandas dataframe. Rows are reactions, columns are conditions
            
            Returns:
                Maximal specific activity for each enzyme-reaction pair, 
                the condition in which it was found, the metabolic pathway 
                associated with the reaction and the carbon source on which 
                the cells were grown.

                Notice that maximal specific activities are given for the sum 
                of all associated enzymes, thus represent the weighted average
                of the specific activites of the isoenzymes. Being a weighted
                average, it means that the values underestimate the maximal 
                potential rate.
        '''

        specific_activity.dropna(thresh=minimal_conditions, inplace=True)

        SAmax = pd.DataFrame(index=specific_activity.index)
        reactions = map(self.model.reactions.get_by_id, SAmax.index)
        subsystems = map(lambda r: r.subsystem, reactions)

        SAmax['max specific activity [umol/mg/min]'] = specific_activity.max(axis=1)        
        SAmax['subsystem'] = subsystems
        SAmax['condition'] = specific_activity.idxmax(axis=1)
        SAmax['carbon source'] = map(lambda x: self.gc['carbon source'][x], SAmax.condition)

        return SAmax

    def get_second_max(self):
        '''
            Finds the second maximal kapp value by reaction
            
            Arguments:
                self
            Returns:
                Pandas Series with reactions as index and snd max as values
                
        '''

        rate = self.kapp.mul(self.p_per_as, axis=0)
        rate.dropna(how='all', inplace=True)
        second = pd.Series(index=rate.index)
        for r in rate.index:
            array = sorted(rate.loc[r])
            second[r] = array[-2]

        return second
        
if __name__ == "__main__":
    R = rates()
