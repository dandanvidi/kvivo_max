from catalytic_rates import rates
from concentration_dependant_effects import MM_kinetics
from figure_correlation import generate_figure
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from figure_correlation import generate_figure
from scipy import stats
from uncertainties import ufloat_fromstr

fontsize = 30
R = rates()

#reactions = ['CYTK1','DAPE','DHORTS_reverse','G3PD2_reverse',
#             'GLUDy_reverse','GLUR_reverse','HSST','MDH','PGI',
#             'PRAGSr','PTPATi','SERAT']
index = R.kcat.index & R.kmax.index
reactions = [R.rxns[x] for x in index]
kcat = R.kcat['kcat per active site [s-1]'][index]
kmax = R.kmax['kmax per active site [s-1]'][index]

saturation = pd.DataFrame(index=index, columns=['glc','ac','glyc'])
backwrdflx = pd.DataFrame(index=index, columns=['glc','ac','glyc'])

for c in ['glc','ac','glyc']:
    f = pd.DataFrame.from_csv('../res/conc_dependant_effects_on_%s.csv' %c, sep='\t')
    saturation[c] = (f['under saturation'][index])
    backwrdflx[c] = (f['backward flux'][index])

s = pd.Series(index=index)
t = pd.Series(index=index)
for r in index:        
    c = R.gc['media_key'][R.kmax['condition'][r]] 
    if c in ['glc','ac','glyc','xyl_D']:
        if c == 'xyl_D': c = 'glc'
        if saturation[c][r]>0 and backwrdflx[c][r]>0:
            s[r] = saturation[c][r]
            t[r] = backwrdflx[c][r]


invivo_kcat = kmax / (s * t)
new_index = invivo_kcat.dropna().index
names = {k:list(R.rxns[k].genes)[0].name for k in new_index}
fig = plt.figure(figsize=(12.5,12.5))
ax= plt.axes()
report_before = generate_figure(kcat[new_index], kmax[new_index], fig, ax, 
                                labels=names, hide_overlap=False)

r1, pval1 = stats.pearsonr(np.log10(kcat[new_index]), np.log10(kmax[new_index]))

report_after = generate_figure(kcat[new_index], invivo_kcat[new_index], fig, ax, fit=True, color='k')

r2, pval2 = stats.pearsonr(np.log10(kcat[new_index]), np.log10(invivo_kcat[new_index]))

def stacked_residual(x, y, saturation_effect, thermodynamic_effect, ax):
    for (x0, y0, s, t) in zip(x, y, saturation_effect, thermodynamic_effect):
        ax.vlines(x0, y0, y0/s, lw=5, colors='#8383FF')
        ax.vlines(x0, y0/s, y0/(s*t), lw=5, colors='#FFDA47')

stacked_residual(kcat[new_index], kmax[new_index], s[new_index], t[new_index], ax)

ax.set_xlim(1e-1*2,1e3*4)
ax.set_ylim(1e-1*2,1e3*4)

ax.set_ylabel(r'$k_{\mathrm{max}}^{\mathrm{vivo}}\,\left[s^{-1}\right]$ (in vivo)', 
              size=fontsize, style='italic')
ax.set_xlabel(r'$k_{\mathrm{cat}}\,\left[s^{-1}\right]$ (in vitro)', 
              size=fontsize, style='italic')
ax.tick_params(axis='both', which='both', top='off', right='off')


[tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]
plt.tight_layout()

