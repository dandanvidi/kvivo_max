from catalytic_rates import rates
from scipy import stats, odr
import numpy as np
import matplotlib.pyplot as plt
import math

def generate_figure(x, y, fig, ax, color='#9E7E5E', edge='none', 
                               yerr='none', labels={}, scatter_size=40, hide_overlap=True,
                               fontsize=18, fit=False, zorder=3):
    
    logx = np.log10(x)
    logy = np.log10(y)
    
    ax.scatter(x, y,s=scatter_size, c=color, marker='o', edgecolor=edge, zorder=zorder)
    
    if yerr != 'none':
        ax.errorbar(x, y, 
                    yerr=yerr, barsabove=False, 
                    fmt=None, ecolor='0.3', zorder=10)
                
    ax.plot([1e-4, 1e4], [1e-4,1e4], 'k', ls='-', lw=2, zorder=5)
     
    #Define function for scipy.odr
    fit_func = lambda B,x: B[0]*x + B[1]
    
    #Fit the data using scipy.odr
    Model = odr.Model(fit_func)
    Data = odr.RealData(logx, logy)
    Odr = odr.ODR(Data, Model, beta0=[1,1])
    output = Odr.run()
    #output.pprint()
    beta = output.beta                

    if fit:  
        edge = np.array([-4, 4])
        ax.plot([1e-4, 1e4], 10**fit_func(beta, edge), color='#9E7E5E', ls=':', lw=3, zorder=1)
        
        
    ax.set_xscale('log', nonposx='clip')
    ax.set_yscale('log', nonposy='clip')
                    
    if labels!={}:
        add_labels(x, y, labels, ax, fig, fontsize, hide_overlap)
                                
    return output

def add_labels(x, y, labels, ax, fig, fontsize=18, hide_overlap=True):
    ann = []
    for r, name in labels.iteritems():
        if x[r]>y[r]:
            ann.append(ax.text(x[r], y[r]/1.1, name, 
                                ha='center', va='top', zorder=5, size=fontsize))
        if x[r]<y[r]:
            ann.append(ax.text(x[r], y[r]*1.1, name,
                                ha='center', va='bottom', zorder=5, size=fontsize))
                                    
        mask = np.zeros(fig.canvas.get_width_height(), bool)
        fig.canvas.draw()
        for i, a in enumerate(ann):
            bbox = a.get_window_extent()
            x0 = int(bbox.x0)
            x1 = int(math.ceil(bbox.x1))
            y0 = int(bbox.y0)
            y1 = int(math.ceil(bbox.y1))
        
            s = np.s_[x0:x1, y0:y1]
            if hide_overlap:
                if np.any(mask[s]):
                    a.set_visible(False)
                else:
                    mask[s] = True
            else:
                mask[s] = True

if __name__ == "__main__":
    
    R = rates()

    fontsize = 30

    fig = plt.figure(figsize=(12.5,12.5))
    ax = plt.axes()
    x = R.kcat['kcat per active site [s-1]'].dropna()
    y = R.kmax['kmax per active site [s-1]'].dropna()
    index = x.index & y.index
    x = x[index]
    y = y[index]
    res = np.abs(np.log10(x) - np.log10(y))
    names = {k:list(R.rxns[k].genes)[0].name for k in index}
    report = generate_figure(x, y, fig, ax, fit=True, labels=names, hide_overlap=False)
    
    rmse = np.sqrt( report.sum_square / len(x) )
    r, pval = stats.pearsonr(np.log10(x), np.log10(y))
    
    ax.set_ylabel(r'$k_{\mathrm{max}}^{\mathrm{vivo}}\,\left[s^{-1}\right]$ - in vivo', 
                  size=fontsize, style='italic')
    ax.set_xlabel(r'$k_{\mathrm{cat}}\,\left[s^{-1}\right]$ - in vitro', 
                  size=fontsize, style='italic')
    ax.tick_params(axis='both', which='both', top='off', right='off')
    
    [tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]

    ax.set_xlim(1e-3/4,4*1e3)
    ax.set_ylim(1e-3/4,4*1e3)
    
    plt.tight_layout()
