from catalytic_rates import rates
import pandas as pd
import matplotlib.pyplot as plt
from figure_correlation import generate_figure

R = rates()


fva_0 = pd.DataFrame.from_csv('../data/flux_variability[mmol_gCDW_s]_relaxation=0.csv')

kcat = R.kcat['kcat per active site [s-1]']
kmax = R.kmax['kmax per active site [s-1]']

ind = kcat.index & kmax.index

x = kcat[ind]
y = kmax[ind]

v = pd.Series(index=ind)
for r in y.index:
    c = R.kmax['condition'][r]
    v[r] = R.flux_data[c][r]
E = v[ind]/y

min_kmax = y - fva_0.minimum[ind]/E 
max_kmax = fva_0.maximum[ind]/E - y

fig, ax = plt.subplots(1, 1, figsize=(6,6), sharey=True)
fontsize=20

labels = {'PDX5POi':'pdxH'}
report = generate_figure(x, y, fig, ax,labels=labels, yerr = [min_kmax, max_kmax],fit=True)

ax.set_ylabel(r'$k_{\mathrm{max}}^{\mathrm{vivo}}\,\left[s^{-1}\right]$ (in vivo)', 
              size=fontsize, style='italic')
ax.set_xlabel(r'$k_{\mathrm{cat}}\,\left[s^{-1}\right]$ (in vitro)', 
              size=fontsize, style='italic')
ax.tick_params(axis='both', which='both', top='off', right='off')

[tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]
ax.set_xlim(1e-3/4,4*1e3)
ax.set_ylim(1e-3/4,4*1e3)
plt.tight_layout()
