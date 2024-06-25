# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 09:58:01 2023

@author: arch9809
"""

# In[]

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px


# %%
dfGH2O = -237.2  # kJ/mol
dfGO2aq = 16.5  # kJ/mol
# DeltaG of complete oxidation of organic matter
gamma_O2 = 4
dGred_O2 = 2 * dfGH2O - dfGO2aq


def dGox_func(gamma):
    return 60.3 - 28.5 * (4 - gamma)


def dCG_O2(gamma):
    return dGox_func(gamma) + (gamma / gamma_O2) * dGred_O2


def efficiency(gamma):
    # Check if gamma is a scalar or a vector
    if np.isscalar(gamma):
        gamma = np.array([gamma])

    if np.any(gamma > 8):
        print("Error. Wrong values of DR of Substrate or Product")
        return

    # Ts = 298
    # pre-assign std G of formation from CHNOSZ
    dfGH2O = -237.2  # kJ/mol
    dfGO2aq = 16.5  # kJ/mol
    gamma_B = 4.2  # e- mol/Cmol
    # define electron acceptor and some other DR
    gamma_eA = 4
    dGred_eA = 2 * dfGH2O - dfGO2aq

    # growth yield calculations
    dGrX = np.where(gamma < 4.67, -(666.7 / gamma + 243.1), -(157 * gamma - 339))  # kJ/Cmol biomass

    # Anabolic reaction
    dCGX = dCG_O2(gamma_B)
    dGana = (gamma_B / gamma) * dCG_O2(gamma) - dCGX
    dGana1 = dGana

    dG_ox = dGox_func(gamma)
    dGcat = dG_ox + gamma * dGred_eA / gamma_eA
    Y = dGcat / (dGrX - dGana1 + gamma_B / gamma * dGcat)

    return Y


efficiency(4)
# %% plantC_geolocations
plt.close('all')
plant_data = pd.read_excel('../collated data/NMR_MMM.xlsx')
corrdata = pd.read_excel('../collated data/corrdata.xlsx')
len(plant_data['Study'].unique())
len(corrdata['Study'].unique())

# %% Figue S1
# Filter non-NaN values for 'Latitude' and 'Longitude' columns in 'srcdata'
lat = corrdata.loc[~np.isnan(corrdata['Latitude']), 'Latitude']
long = corrdata.loc[~np.isnan(corrdata['Longitude']), 'Longitude']

loc = pd.DataFrame({'Latitude': lat, 'Longitude': long})
fig = px.scatter_geo(loc, lat='Latitude', lon='Longitude', color_discrete_sequence=['red'])
fig.update_geos(projection_type="natural earth")  # You can choose a different projection type
fig.update_layout(title_text='Your Title', margin=dict(l=0, r=0, b=0, t=0))
fig.write_html('first_figure.html', auto_open=True)
fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
export_params = dict(format='png', width=1200, height=800, scale=3)
# fig.write_image("../results/plantC_geolocations.svg", **export_params)


# %%

replacements = {'needle':'needles','leaf': 'broadleaves', 'stem': 'wood',
                'lichen': 'lichen+moss', 'moss': 'lichen+moss', 'grass+herbs': 'grasses+herbs',
                'grass+herbs':'grasses+herbs'}
corrdata['Csource'] = corrdata['Csource'].replace(replacements)

sns.set(rc={'axes.labelsize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'grid.linewidth': 0.25})
sns.set_style("whitegrid")
# Extract columns for C fractions
col_c = ['Carbohydrates', 'Proteins', 'Lignins', 'Lipids', 'Carbonyls']
T2_c = pd.melt(corrdata, id_vars=['Csource'], value_vars=col_c,
               var_name='NMR', value_name='C fractions')

T2_c['Csource'].unique()
# Extract columns for PA fractions
col_pa = ['labile', 'AS', 'AIS']
T2_pa = pd.melt(corrdata, id_vars=['Csource'], value_vars=col_pa,
                var_name='NMR', value_name='PA fractions')

hue_counts = corrdata['Csource'].value_counts().sort_index()
hueorder = list(hue_counts.index)


# Set up subplots
fig, axes = plt.subplots(2, 1, figsize=(11, 8))

# Plot for C fractions
ax1 = axes[1]
grped_bplot_c = sns.boxplot(x='NMR', y='C fractions', hue="Csource", hue_order=hueorder,data=T2_c, palette="Set1", ax=ax1)
sns.stripplot(x='NMR', y='C fractions', hue="Csource", hue_order=hueorder, data=T2_c, jitter=True, dodge=True, size=4, marker='o', palette="Pastel1", alpha=0.5, ax=ax1)
ax1.set_xlabel('')
ax1.set_ylabel('NMR fractions (g/g litter)')

handles, labels = grped_bplot_c.get_legend_handles_labels()
legend_labels = labels[0:7]
legend_labels_with_count = [f"{label} (n={count})" for label, count in zip(legend_labels, hue_counts)]
l = ax1.legend(handles[0:7], legend_labels_with_count, fontsize='12')
# l.set_frame_on(False)
# specify just one legend


# Plot for PA fractions
ax2 = axes[0]
grped_bplot_pa = sns.boxplot(x='NMR', y='PA fractions', hue="Csource", hue_order=hueorder,data=T2_pa, palette="Set1", ax=ax2)
sns.stripplot(x='NMR', y='PA fractions', hue="Csource", hue_order=hueorder, data=T2_pa, jitter=True, dodge=True, size=4, marker='o',palette="Pastel1", alpha=0.5, ax=ax2)

ax2.legend('')
ax2.set_xlabel('')
ax2.set_xticklabels(['Labile (nonpolar and water soluble)', 'Acid soluble', 'Acid insoluble'])
ax2.set_ylabel('PA fractions (g/g litter)')
ax2.get_legend().set_frame_on(False)
plt.tight_layout()
plt.show()
plt.savefig("../results/Figure2.png", dpi=300)
plt.savefig("../results/Figure2.svg", dpi=300)

#%% Figure S3
gamma = 4-corrdata['Cox']
corrdata['cue'] = efficiency(gamma.values)
hue_counts = corrdata['Csource'].value_counts().sort_index()

sns.set(rc={'axes.labelsize': 14, 'xtick.labelsize': 12, 'ytick.labelsize': 12})
sns.set_style("whitegrid")

plt.figure(figsize=(9, 4.5))
ax1 = plt.subplot(2, 1, 1)
sns.violinplot(x="Csource", y="Cox", data=corrdata, palette="Set1", ax=ax1, order=hueorder)
sns.stripplot(x="Csource", y="Cox", data=corrdata,order=hueorder, jitter=True,
              color="black", size=3, alpha=0.25, ax=ax1)

plt.xlabel('')
plt.ylabel('NOSC')
plt.tight_layout()
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.ylim((-1, 0.25))
ax1.grid(True, linewidth=0.25)  # Retain grid lines with specified line width


gamma = 4-plant_data['Cox']
plant_data['cue'] = efficiency(gamma.values)
sns.set_style("whitegrid")
ax2 = plt.subplot(2, 1, 2)
sns.violinplot(x="Csource", y="cue", data=corrdata, palette="Set1", order=hueorder)
sns.stripplot(x="Csource", y="cue", data=corrdata, order=hueorder, jitter=True,
              color="black", size=3, alpha=0.25)
plt.xlabel('')
plt.ylabel('CUE')
plt.tight_layout()
ax2.grid(True, linewidth=0.25)
# sns.despine()
plt.ylim((0.45, 0.65))
# Add counts below each x-label
for i, count in enumerate(hue_counts):
    # Calculate the position for annotation
    x_pos = i
    y_pos = ax2.get_ylim()[0] - (ax2.get_ylim()[1] - ax2.get_ylim()[0]) * 0.2  # Adjust this factor for the position of the annotation
    # Annotate the plot
    ax2.text(x_pos, y_pos, f'(n={count})', ha='center', va='top', fontsize=10)

plt.show()
plt.tight_layout()
plt.savefig("../results/FigureA2.png", dpi=300)
plt.savefig("../results/FigureA2.svg", dpi=300)
