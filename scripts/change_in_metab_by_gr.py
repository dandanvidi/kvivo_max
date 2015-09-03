import pandas as pd
import matplotlib.pylab as plt
import plotly.plotly as py
import numpy as np

py.sign_in('DanDavidi', '4noole85cy')

metab = pd.DataFrame.from_csv('../data/bennett_metabolite_concentrations[mM].csv')
metab.drop(['comments'], axis=1, inplace=True)
metab.dropna(how='any', inplace=True)
metab_names = pd.DataFrame.from_csv('../data/model_metabolites.csv', sep='\t')

fig = plt.figure(figsize=(9,6))
ax = plt.axes()
gr = [0.58,0.47,0.3]
for row in metab.iterrows():
    cid = row[0]
    values = (row[1]/row[1].mean()).values
    name = str(cid)
    try:
        name=metab_names[metab_names.kegg_id == cid]['name'].values[0]
    except:
        if name == "C00077":
            name = 'L-Ornithine'
        else:
            print "%s not found" %name
    
    increase = np.sort(values)
    decrease = increase[::-1]

    if (values==increase).all():
        print name        
        ax.plot(gr,values,'o-', lw=1.5, color='r', label=name, zorder=3)
    elif (values==decrease).all():
        ax.plot(gr,values,'o-', lw=1.5, color='b', label=name, zorder=3)
    else:
        ax.plot(gr,values,'o-', lw=1.5, color='0.5', alpha=0.5, label=name, zorder=0)
ax.set_xlabel('growth rate [1/h]', size=15)
ax.set_ylabel('normalized metabolite concentration', size=15)
ax.tick_params(axis='x', labelsize=15)
ax.tick_params(axis='y', labelsize=15)
ax.set_xlim(0.25,0.63)
ax.set_xticks(np.arange(0.3, 0.61,0.05))
ax.set_ylim(0,3.2)
ax.text(0.30,3,'acetate', ha='center', size=15)
ax.text(0.47,3,'glycerol', ha='center', size=15)
ax.text(0.58,3,'glucose', ha='center', size=15)
#plot_url = py.plot_mpl(fig, file_name="metabolite concentrations across conditions")


plt.savefig('../res/metab_change.pdf')