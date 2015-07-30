from cobra.io.sbml import create_cobra_model_from_sbml_file
import pandas as pd
from component_contribution.kegg_reaction import KeggReaction

class add_to_model(object):
    
    def __init__(self, cobra_model):
    
        self.gene_info = pd.DataFrame.from_csv('../data/model_genes.csv',sep='\t')
        self.metab_info = pd.DataFrame.from_csv("../data/model_metabolites.csv",sep='\t')        
        self.metab_info.dropna(inplace=True)
        self.add_to_genes(cobra_model)
        self.add_to_metab(cobra_model)
        self.add_to_reac(cobra_model)
        
    def add_to_genes(self ,cobra_model):

        for g in cobra_model.genes:
            g.name = self.gene_info['uniprot_primary_name'][g.id]
            g.MW = self.gene_info['MW_Da'][g.id]
            g.length = self.gene_info['length_aa'][g.id]

    def add_to_metab(self, cobra_model):

        for m in cobra_model.metabolites:
            try:
                m.CID = self.metab_info.kegg_id[m[:-2]]
            except KeyError:
                m.CID = None
        
    def add_to_reac(self, cobra_model):
    
        for r in cobra_model.reactions:
            CIDS = dict(zip(r.metabolites.keys(), map(lambda x: x.CID, r.metabolites.keys())))
            if None in CIDS.values():
                r.kegg_reaction = None
            else:
                sparse = {CIDS[m]:v for m,v in r.metabolites.iteritems()
                                            if CIDS[m]!='C00080'}        
                r.kegg_reaction = KeggReaction(sparse)
    
if __name__ == "__main__":
    model = create_cobra_model_from_sbml_file('../data/iJO1366.xml')
    add_to_model(model)
