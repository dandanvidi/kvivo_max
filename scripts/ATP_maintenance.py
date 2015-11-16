# bugs not fixed

from catalytic_rates import rates
import pandas as pd
import matplotlib.pyplot as plt
from figure_correlation import generate_figure
from pFBA_implementation import perform_pFBA
import numpy as np
from scipy import stats, odr

R = rates()

def customize_ATPM(x):

    R.rxns['ATPM'].lower_bound = 3.15*x #mmol/gCDW/h
    fluxes = pd.DataFrame(index=R.rxns.keys(), columns=R.gc.index)
    for c in R.gc.iterrows():
        gr = c[1]['growth rate [h-1]']
        cs = c[1]['media_key']
        ur = c[1]['uptake rate [mmol gCDW-1 h-1]']
        if np.isnan(ur):
            ur = 18.5
        fluxes[c[0]] = perform_pFBA(R.model, cs, gr, ur)
        print "- %s" %c[0]
    
    fluxes.index.name = 'reaction'
    fluxes.to_csv('../res/flux[mmol_h_gCDW]_maintenance_zero.csv') #export results
    R.rxns['ATPM'].lower_bound = 3.15 #reset bounds for the ATPM reaction

def calculate_new_kmax(fluxes):
    fluxes = fluxes / 3600
#    kapp = R.get_rate_per_subunit(fluxes, R.p)
    kmax = R.get_kmax(kapp)
    return kmax

def fit(x,y):
    #Fit the data using scipy.odr
    fit_func = lambda B,x: B[0]*x + B[1]
    Model = odr.Model(fit_func)
    Data = odr.RealData(x, y)
    Odr = odr.ODR(Data, Model, beta0=[1,1])
    return Odr.run()
    
if __name__ == "__main__":
    
#    ATP_maintenance()    
    old_kmax = R.kmax['kmax per active site [s-1]']
    R.flux_data = pd.DataFrame.from_csv('../res/flux[mmol_h_gCDW]_maintenance_zero.csv') / 3600
    R.kapp = R.get_kapp()
    R.kmax = R.get_kmax(R.kapp)
    zero_kmax = R.kmax['kmax per active site [s-1]']

    tmp = (old_kmax / zero_kmax).dropna()
    f = fit(np.log(old_kmax[tmp.index]), np.log(zero_kmax[tmp.index]))
#    main_zero = pd.DataFrame.from_csv('../res/flux[mmol_h_gCDW]_maintenance_zero.csv')       
#    main_high = pd.DataFrame.from_csv('../res/flux[mmol_h_gCDW]_maintenance_high.csv')   
    
#    kmax = R.kmax['kmax per chain [s-1]']
#    kmax_mz = calculate_new_kmax(main_zero)['kmax per chain [s^-1]'][kmax.index]
#    kmax_mh = calculate_new_kmax(main_high)['kmax per chain [s^-1]'][kmax.index]
    
