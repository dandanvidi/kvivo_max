from catalytic_rates import rates
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import sys, os
from collections import defaultdict
sys.path.append(os.path.expanduser('~/git/across-projects'))
from color import ColorMap

R = rates()
#gc = R.gc
gc = R.gc[R.gc['growth mode']=='batch']
remove = R.gc['comments'].dropna()
gc = gc.drop(remove.index)

index = gc.index & R.kapp.columns
gr = gc['growth rate [h-1]'][index]
enzymes = R._convert_mmol_gCDW_to_mg_gCDW(R.expression_data)[gr.index]
V = R.flux_data[gr.index]
kmax = R.kmax['kmax per active site [s-1]'].dropna()
kcat = R.kcat['kcat per active site [s-1]'][kmax.index]

kapp = R.kapp[gr.index].loc[kcat.index]

# kapp / kcat - a measure of the catalytic efficiency (CA) of enzymes 
tmp = kapp.mul(R.p_per_as, axis=0)
CA = tmp.div(kcat, axis=0).dropna(how='all')


WCE = pd.DataFrame(index=CA.index, 
                 columns=CA.columns) # mg/gCDW
E = pd.DataFrame(index=CA.index, 
                 columns=CA.columns) # mg/gCDW
for r in CA.index:
    b = list(R.rxns[r].genes)[0].id
    E.loc[r] = enzymes.loc[b]
    WCE.loc[r] = E.loc[r] * CA.loc[r]

cond = gr.index[-1]
a = WCE[cond].dropna()
a.sort()
x = a[-1:].index

#WCA.drop(x, axis=0, inplace=True)
#E.drop(x, axis=0, inplace=True)

matric = WCE.div(E.sum(), axis=1)

fig = plt.figure(figsize=(10,6))
ax = plt.axes(axisbg='0.95')

z = {R.rxns[r].id:R.rxns[r].subsystem for r in WCE.index}
subsystems = defaultdict(list)
for key, value in sorted(z.iteritems()):
    subsystems[value].append(key)

colors = ColorMap(subsystems.keys())
for k,v in subsystems.iteritems():
    array = matric.loc[v]
    narray = array.div(array.mean(axis=1), axis=0)
#    narray.dropna(how='any', inplace=True)
    g = gr[narray.columns]
    print len(g), len(narray.columns)
    
    (a,b), cov = curve_fit(lambda a,b,x: a*x+b, g, narray.sum())
#    print k, b
    ax.plot(g, a*g+b, c=colors[k], marker='o', label=k)
    
#plt.scatter(gr, matric.sum(), c='#4DB8FF', edgecolor='none', 
#            s=50, label='relative to $k_{max}^{vivo}$')
#            
#labels = [gc['media'][c] for c in gr.index]
#for i, txt in enumerate(labels):
#    ax.annotate(txt, (gr[i],matric.sum()[i]))
#
#ax.set_xlim(0.1,0.8)
#ax.set_ylim(0.1,0.45)
plt.grid()
ax.set_xlabel('growth rate [h$^{-1}$]', size=15)
ax.set_ylabel('enzyome efficiency', size=15)
ax.legend(scatterpoints=1)
plt.tight_layout()
#x = a[-10:].index
#colors = ColorMap(x)
##plt.figure()
##for r in x:
##    plt.plot(gr, CA.loc[r]/CA.loc[r][0], label=r, c = colors[r], marker='o')
##    plt.legend()
#
#plt.figure()
#for r in x:
#    plt.plot(gr, E.loc[r], label=r, c = colors[r], marker='o')
#    plt.legend()
#    
#plt.figure()
#for r in x:
#    plt.plot(gr, V.loc[r], label=r, c = colors[r], marker='o')
#    plt.legend()
#    plt.figure()
#    plt.hist(a, 40)
#    print a['PSERT'] / gr[c]
#
#[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
#[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]
#
#plt.tight_layout()
#
##fig = plt.figure(figsize=(6,6))
##ax = plt.axes()
##
##for i, r in enumerate(efficiency.index):
##    x = weighted_concentration.loc[r].astype('float').dropna()
##    y = efficiency.loc[r].astype('float').dropna()
##
##    z = gr[x.index]
##    z.sort
##    x = x[z.index]    
##    y = y[z.index]
##    plt.scatter(x, y, c='b')
##    if i >= 0:
##        break
###ax.set_xscale('log')
##
##
###                
##fig = plt.figure(figsize=(6,6))
##ax = plt.axes()
##conditions = list(efficiency.columns)
###conditions = [conditions[0:12]]# + [conditions[-1]]
##cm = plt.cm.get_cmap('Greens')
###cm = ['r']#, 'b']
##for i, j in enumerate(conditions):
##    x = weighted_concentration[j].astype('float').dropna()
##    y = efficiency[j].astype('float').dropna()
##    
##    index = x.index & y.index
##    x = x.loc[index]
##    y = y.loc[index]
##    
##    plt.scatter(x, y, c=cm(1.0*i/len(conditions)), 
##                edgecolor='none')
##
##[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
##[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]
##ax.set_xlabel('log E [mg/gCDW]', size=15)
##ax.set_ylabel('catalytic efficiency', size=15)
##ax.set_xscale('log')
##ax.set_xlim(1e-4,1e1)
##ax.set_ylim(0,1.1)
##
##
#
##plt.scatter(gr, matric_theoretical, c='#8500AD', edgecolor='none', 
##            s=50, label='relative to $k_{cat}$')
##
##a, p = curve_fit(lambda x, a: a*x, gr, matric_theoretical)
##
###gr = np.append(gr,1/a)
###plt.plot(gr, a* gr)
##
###for i, c in enumerate(R.gc.index):
###    if R.gc['growth mode'][c] == 'batch':
####        plt.scatter(gr[i], matric[c], c='k', edgecolor='none', 
####                    s=50)
###        plt.scatter(gr[i], matric_theoretical[c], c='k', edgecolor='none', 
###                s=50)
##ax.set_xlim(0,1)
##ax.set_ylim(0,0.4)
###plt.legend(scatterpoints=1, loc=3, fontsize=15)#ax.plot([0,0.8],[0,0.8], '#993333')
##plt.grid()
##ax.set_xlabel('growth rate [h$^{-1}$]', size=15)
##ax.set_ylabel('enzyome saturation', size=15)
##
##[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
##[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]
##
##plt.tight_layout()
##
##plt.savefig('../res/mass_efficiency.pdf')
##
##
###fig = plt.figure(figsize=(6,6))
###ax = plt.axes()
###conditions = list(efficiency.columns)
###conditions = [conditions[-10:-1]]# + [conditions[-1]]
###cm = plt.cm.get_cmap('Greens')
###cm = ['b']#, 'b']
###for i, j in enumerate(conditions):
###    x = weighted_concentration[j].astype('float').dropna()
###    y = efficiency[j].astype('float').dropna()
###    
###    index = x.index & y.index
###    x = x.loc[index]
###    y = y.loc[index]
###    
###    plt.scatter(x, y, c=cm[i], 
###                edgecolor='none')
###
###[tick.label.set_fontsize(15) for tick in ax.xaxis.get_major_ticks()]
###[tick.label.set_fontsize(15) for tick in ax.yaxis.get_major_ticks()]
###ax.set_xlabel('log E [mg/gCDW]', size=15)
###ax.set_ylabel('catalytic efficiency', size=15)
###ax.set_xscale('log')
###ax.set_xlim(1e-4,1e1)
###ax.set_ylim(0,1.1)
