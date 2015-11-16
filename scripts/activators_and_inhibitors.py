# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 11:17:15 2015

@author: dan
"""

from catalytic_rates import rates
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
from scipy import stats

R = rates()
index = R.kcat.index & R.kmax.index
res = (R.kcat['kcat per active site [s-1]']/R.kmax['kmax per active site [s-1]'])[index]

reactions = R.kcat['EC'][index]
effectors = pd.DataFrame.from_csv("../data/ecoli_effectors_by_ec.csv", sep='\t')

activators = effectors['activators'].dropna()
inhibitors = effectors['inhibitors'].dropna()

count_activators = Counter(activators.index)
count_inhibitors = Counter(inhibitors.index)

number_of_inhibitors = pd.DataFrame.from_dict(count_inhibitors.items()).set_index(0)
number_of_inhibitors.replace({0:{v:k for k,v in reactions.to_dict().iteritems()}})

y = pd.DataFrame.from_dict(count_inhibitors.items()).set_index(0)[1]
y = y[y.index & reactions.index]
x = res[y.index]

a = len(x[x<1]) # above the line and have inhibitors
b = len(x[x>1]) # below the line and have inhibitors
c = len(res[res<1]) - len(x[x<1]) #above the line but dont have inhibitors
d = len(res[res>1]) - len(x[x>1]) #below the line but dont have inhibitors

oddsratio, pvalue = stats.fisher_exact([[a,b],[c,d]])
print pvalue

