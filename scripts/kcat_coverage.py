# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 11:59:28 2015

@author: dan
"""
from catalytic_rates import rates
import pandas as pd
from xml.dom.minidom import parse
import numpy as np
from cobra.io.sbml import create_cobra_model_from_sbml_file
from model_addons import add_to_model

def map_model_reactions_to_EC(model_fname):
    
    ''' parse the sbml model file to extract EC numbers of reactions '''
    
    document = parse(model_fname)
    ec_list = []
    for r_elem in document.getElementsByTagName('reaction'):
        for p_elem in r_elem.getElementsByTagName('p'):
            val = p_elem.childNodes[0].nodeValue
            if val.find('EC Number:') != -1:
                ec = val[11:].encode('utf-8') 
                ec_list.append(ec)
            
    return ec_list
    
a = map_model_reactions_to_EC("../data/iJO1366.xml")
m = create_cobra_model_from_sbml_file("../data/iJO1366.xml")
add_to_model(m)

#b = set(filter(None, a))
brenda = pd.DataFrame.from_csv("../data/BRENDA_data.csv")
brenda = brenda[brenda["Organism ID"]==6] # ecoli reactions
brenda = brenda[brenda["kcat"]>0]
measured_ecs = set(brenda.index)

out = pd.DataFrame(index=a)
out['reaction'] = [r.id for r in m.reactions]
out['subsystem'] = [r.subsystem for r in m.reactions]
out['kcat'] = [1 if ec in brenda.index else 'NaN' for ec in out.index]
out[out['kcat']!=1]
out.drop('', inplace=True)
out.drop('', inplace=True)
out['EC'] = out.index
out.set_index('reaction', inplace=True)

out = out[['EC', 'subsystem', 'kcat']]
out.to_csv('../res/missing_kcat.csv', sep='\t')