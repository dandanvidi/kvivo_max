import pandas as pd
import pickle
from cobra.test import salmonella_pickle
import numpy as np
from cPickle import load
from cobra.manipulation import initialize_growth_medium
from cobra.test import salmonella_pickle

with open(salmonella_pickle) as in_file:
    cobra_model = load(in_file)
    
initialize_growth_medium(cobra_model, 'LB')
LB = cobra_model.media_compositions['LB']

gc = pd.DataFrame.from_csv('../data/growth_conditions.csv')

media = {}
for c in gc.iterrows():
    condition = c[0]    
    c = c[1]
    if c['media_key'] not in ['lb', 'glyc_aa']:
        k = 'EX_'+c['media_key']+'_e'
        v = -c['uptake rate [mmol gCDW-1 h-1]']
        if np.isnan(v):
            v = -1000
        media[condition] = {'EX_o2_e':-1000, k:v}
lb = pickle.load( open( "../data/LB_media.p", "rb" ) )
lb = {k:-1000 for k in lb.keys()}
lb['EX_o2_e'] = -1000
glyc_aa =['glyc', 'ala_L', 'arg_L', 'asp_L','cys_L','glu_L','gly','his_L',
         'ile_L','leu_L','lys_L','met_L','phe_L','pro_L','ser_L','thr_L',
         'trp_L','tyr_L','val_L','asn_L','gln_L','ade','ura']
glyc_aa = {'EX_'+k+'_e':-1000 for k in glyc_aa}
glyc_aa['EX_o2_e'] = -1000
media['LB_BATCH_mu=1.9_S'] = lb
media['GLYC+AA_BATCH_mu=1.27_S'] = glyc_aa

pickle.dump( media, open( "../data/media.p", "wb" ) )