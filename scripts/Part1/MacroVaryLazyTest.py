#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import sys

from numpy.random import normal

import cobra

sys.path.append('../')
from uncBiomFuncs import file_name, make_pfba_csv, make_fluxes_dataframe, get_coeff_without_gam


# In[2]:


coeff_std = 0.05 # Set this to 0.05, 0.2, and 0.3 and rerun
size = 10000


# In[3]:


dir_model = 'EcoliGecko.json'
model = cobra.io.load_json_model(dir_model)
model.solver = 'gurobi'

#revBetaOx = ['ACACT1r', 'ACACT2r', 'ACACT3r', 'ACACT4r', 'ACACT5r',
#             'ACACT6r', 'ACACT7r', 'ACACT8r']

revBetaOx = ['arm_ACACT1r', 'arm_ACACT1r_REV', 'arm_ACACT2r', 'arm_ACACT2r_REV', 'arm_ACACT3r', 'arm_ACACT3r_REV',
            'arm_ACACT4r', 'arm_ACACT4r_REV', 'arm_ACACT5r', 'arm_ACACT5r_REV', 'arm_ACACT6r', 'arm_ACACT6r_REV',
            'arm_ACACT7r', 'arm_ACACT7r_REV', 'arm_ACACT8r', 'arm_ACACT8r_REV']
for rxnid in revBetaOx:
    rxn = model.reactions.get_by_id(rxnid)
    rxn.upper_bound = 0
    
#model.reactions.POR5.lower_bound = 0
model.reactions.arm_POR5.lower_bound = 0

config = cobra.Configuration()
config.tolerance = 1e-9
zero_tol = 1e-9

biomId = 'BIOMASS_Ec_iML1515_core_75p37M'
gam_val = 75.37723
dir_case = './output/'
import os
if not os.path.exists(dir_case):
    os.makedirs(dir_case)
    
bMets = get_coeff_without_gam(model, biomId, gam_val)

aas = ['ala__L','arg__L','asn__L','asp__L','cys__L','gln__L','glu__L','gly','his__L','ile__L',
       'leu__L','lys__L','met__L','phe__L','pro__L','ser__L','thr__L','trp__L','tyr__L','val__L']
aas = {k + '_c':1 for k in aas}

nus = ['datp','dctp','dgtp','dttp','ctp','gtp','utp','atp']
nus = {k + '_c':1 for k in nus}


# In[4]:


df_biom = pd.read_excel('Biom_frac_annotated.xlsx', sheet_name='biom')
macro_frac = dict()
for macro in set(df_biom.type):
    macro_frac[macro] = round(df_biom[df_biom.type == macro].biom_frac.sum(), 9)

cols = ['Sum_frac', 'Amino acids', 'Cell wall components', 'DNA', 'RNA',
        'Lipids', 'Cofactors and prosthetic groups', 'Inorganic ions']
df_macro = pd.DataFrame(columns=cols, index=range(0, size+1))
i = 0
for k,v in macro_frac.items():
    df_macro.loc[i, k] = v

seed = -1
for col in df_macro.columns[1:]:
    c0 = df_macro.loc[0, col]
    seed += 1
    np.random.seed(seed)
    coeff_norm = c0 + coeff_std*c0*normal(size=size)
    # Only allow positive value for biom fraction
    coeff_norm = [coeff if coeff > zero_tol else zero_tol for coeff in coeff_norm]
    df_macro.loc[range(1,size+1), col] = coeff_norm
    
df_macro['Sum_frac'] = df_macro.iloc[:, 1:].sum(axis=1)

# Normalize
df_macro.iloc[:, 1:] = df_macro.iloc[:, 1:].div(df_macro.Sum_frac, axis=0)


# In[5]:


met_in_macro = {}
for macro in set(df_biom.type):
    met_in_macro[macro] = df_biom[df_biom.type == macro].id.to_list()
    
round_digit = -int(np.log10(zero_tol))
bMets = get_coeff_without_gam(model, biomId, gam_val)

dfCoeff = pd.DataFrame(index = range(0,size+1), columns=['biomMW'] + list(bMets.keys()))

i = 0
for col in dfCoeff.columns[1:]:
    dfCoeff.loc[i, col] = bMets[col]


# In[6]:


for i in range(1, size+1):
    for macro,mets in met_in_macro.items():
        dfCoeff.loc[i,mets] = dfCoeff.loc[0,mets] * df_macro.loc[i, macro] / df_macro.loc[0, macro]
        
for i in range(1, size+1):
    coeff = 0
    for k,v in aas.items():
        coeff += v*dfCoeff.loc[i, k]
    dfCoeff.loc[i, 'h2o_c'] = -coeff

    coeff = 0
    dfCoeff.loc[i, 'ppi_c'] = 0
    for k,v in nus.items():
        coeff += v*dfCoeff.loc[i, k]
    dfCoeff.loc[i, 'ppi_c'] = -coeff
    
mets_mw = [model.metabolites.get_by_id(i).formula_weight for i in dfCoeff.columns[1:]]
dfCoeff['biomMW'] = -dfCoeff.iloc[:, 1:].multiply(mets_mw, axis=1).sum(axis=1)

# Add GAM back
for col in ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']:
    if col not in dfCoeff.columns:
        dfCoeff[col] = [0.]*dfCoeff.shape[0] # Create array of zeros
for metId in ['atp_c', 'h2o_c']:
    dfCoeff.loc[:, metId] -= gam_val
for metId in ['adp_c', 'pi_c', 'h_c']:
    dfCoeff.loc[:, metId] += gam_val
    
# Rearrange
atpm = ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']
biomrxn = model.reactions.get_by_id(biomId)
mets = atpm + sorted([met.id for met in biomrxn.metabolites.keys() if met.id not in atpm])
dfCoeff = dfCoeff.reindex(['biomMW'] + mets, axis=1)

# Save
fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfCoeff')
dfCoeff.to_csv(dir_case + fname + '.csv', sep='\t')


# #### Step 2 - run pFBA

# In[7]:


fname = file_name(coeff_std=coeff_std, gam_std=None)
dir_pfba = dir_case + 'pFBA/' + fname + '/'
make_pfba_csv(model, dfCoeff, dir_pfba, biomId)


# ####  Step 3 - compile flux dataframe from pFBA csv

# In[8]:


fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfCoeff')
dfCoeff = pd.read_csv(dir_case + fname + '.csv', sep='\t', index_col=0, header=0)

fname = file_name(coeff_std=coeff_std, gam_std=None)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfFlux')
index = dfCoeff.index.to_list()
dfFlux = make_fluxes_dataframe(index, dir_pfba)
dfFlux = dfFlux.loc[:, (dfFlux.abs() > zero_tol).any(axis=0)]
dfFlux.to_csv(dir_case + fname + '.csv', sep='\t')

# Remove pFBA files
import shutil
fname = file_name(coeff_std=coeff_std, gam_std=None)
shutil.rmtree(dir_pfba)


# In[ ]:


coeff_std = 0.1 # Set this to 0.05, 0.2, and 0.3 and rerun
size = 10000


# In[ ]:


df_biom = pd.read_excel('Biom_frac_annotated.xlsx', sheet_name='biom')
macro_frac = dict()
for macro in set(df_biom.type):
    macro_frac[macro] = round(df_biom[df_biom.type == macro].biom_frac.sum(), 9)

cols = ['Sum_frac', 'Amino acids', 'Cell wall components', 'DNA', 'RNA',
        'Lipids', 'Cofactors and prosthetic groups', 'Inorganic ions']
df_macro = pd.DataFrame(columns=cols, index=range(0, size+1))
i = 0
for k,v in macro_frac.items():
    df_macro.loc[i, k] = v

seed = -1
for col in df_macro.columns[1:]:
    c0 = df_macro.loc[0, col]
    seed += 1
    np.random.seed(seed)
    coeff_norm = c0 + coeff_std*c0*normal(size=size)
    # Only allow positive value for biom fraction
    coeff_norm = [coeff if coeff > zero_tol else zero_tol for coeff in coeff_norm]
    df_macro.loc[range(1,size+1), col] = coeff_norm
    
df_macro['Sum_frac'] = df_macro.iloc[:, 1:].sum(axis=1)

# Normalize
df_macro.iloc[:, 1:] = df_macro.iloc[:, 1:].div(df_macro.Sum_frac, axis=0)


# In[ ]:


for i in range(1, size+1):
    for macro,mets in met_in_macro.items():
        dfCoeff.loc[i,mets] = dfCoeff.loc[0,mets] * df_macro.loc[i, macro] / df_macro.loc[0, macro]
        
for i in range(1, size+1):
    coeff = 0
    for k,v in aas.items():
        coeff += v*dfCoeff.loc[i, k]
    dfCoeff.loc[i, 'h2o_c'] = -coeff

    coeff = 0
    dfCoeff.loc[i, 'ppi_c'] = 0
    for k,v in nus.items():
        coeff += v*dfCoeff.loc[i, k]
    dfCoeff.loc[i, 'ppi_c'] = -coeff
    
mets_mw = [model.metabolites.get_by_id(i).formula_weight for i in dfCoeff.columns[1:]]
dfCoeff['biomMW'] = -dfCoeff.iloc[:, 1:].multiply(mets_mw, axis=1).sum(axis=1)

# Add GAM back
for col in ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']:
    if col not in dfCoeff.columns:
        dfCoeff[col] = [0.]*dfCoeff.shape[0] # Create array of zeros
for metId in ['atp_c', 'h2o_c']:
    dfCoeff.loc[:, metId] -= gam_val
for metId in ['adp_c', 'pi_c', 'h_c']:
    dfCoeff.loc[:, metId] += gam_val
    
# Rearrange
atpm = ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']
biomrxn = model.reactions.get_by_id(biomId)
mets = atpm + sorted([met.id for met in biomrxn.metabolites.keys() if met.id not in atpm])
dfCoeff = dfCoeff.reindex(['biomMW'] + mets, axis=1)

# Save
fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfCoeff')
dfCoeff.to_csv(dir_case + fname + '.csv', sep='\t')


# In[ ]:


fname = file_name(coeff_std=coeff_std, gam_std=None)
dir_pfba = dir_case + 'pFBA/' + fname + '/'
make_pfba_csv(model, dfCoeff, dir_pfba, biomId)


# In[ ]:


fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfCoeff')
dfCoeff = pd.read_csv(dir_case + fname + '.csv', sep='\t', index_col=0, header=0)

fname = file_name(coeff_std=coeff_std, gam_std=None)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfFlux')
index = dfCoeff.index.to_list()
dfFlux = make_fluxes_dataframe(index, dir_pfba)
dfFlux = dfFlux.loc[:, (dfFlux.abs() > zero_tol).any(axis=0)]
dfFlux.to_csv(dir_case + fname + '.csv', sep='\t')

# Remove pFBA files
import shutil
fname = file_name(coeff_std=coeff_std, gam_std=None)
shutil.rmtree(dir_pfba)


# In[ ]:


coeff_std = 0.2 # Set this to 0.05, 0.2, and 0.3 and rerun
size = 10000


# In[ ]:


df_biom = pd.read_excel('Biom_frac_annotated.xlsx', sheet_name='biom')
macro_frac = dict()
for macro in set(df_biom.type):
    macro_frac[macro] = round(df_biom[df_biom.type == macro].biom_frac.sum(), 9)

cols = ['Sum_frac', 'Amino acids', 'Cell wall components', 'DNA', 'RNA',
        'Lipids', 'Cofactors and prosthetic groups', 'Inorganic ions']
df_macro = pd.DataFrame(columns=cols, index=range(0, size+1))
i = 0
for k,v in macro_frac.items():
    df_macro.loc[i, k] = v

seed = -1
for col in df_macro.columns[1:]:
    c0 = df_macro.loc[0, col]
    seed += 1
    np.random.seed(seed)
    coeff_norm = c0 + coeff_std*c0*normal(size=size)
    # Only allow positive value for biom fraction
    coeff_norm = [coeff if coeff > zero_tol else zero_tol for coeff in coeff_norm]
    df_macro.loc[range(1,size+1), col] = coeff_norm
    
df_macro['Sum_frac'] = df_macro.iloc[:, 1:].sum(axis=1)

# Normalize
df_macro.iloc[:, 1:] = df_macro.iloc[:, 1:].div(df_macro.Sum_frac, axis=0)


# In[ ]:


for i in range(1, size+1):
    for macro,mets in met_in_macro.items():
        dfCoeff.loc[i,mets] = dfCoeff.loc[0,mets] * df_macro.loc[i, macro] / df_macro.loc[0, macro]
        
for i in range(1, size+1):
    coeff = 0
    for k,v in aas.items():
        coeff += v*dfCoeff.loc[i, k]
    dfCoeff.loc[i, 'h2o_c'] = -coeff

    coeff = 0
    dfCoeff.loc[i, 'ppi_c'] = 0
    for k,v in nus.items():
        coeff += v*dfCoeff.loc[i, k]
    dfCoeff.loc[i, 'ppi_c'] = -coeff
    
mets_mw = [model.metabolites.get_by_id(i).formula_weight for i in dfCoeff.columns[1:]]
dfCoeff['biomMW'] = -dfCoeff.iloc[:, 1:].multiply(mets_mw, axis=1).sum(axis=1)

# Add GAM back
for col in ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']:
    if col not in dfCoeff.columns:
        dfCoeff[col] = [0.]*dfCoeff.shape[0] # Create array of zeros
for metId in ['atp_c', 'h2o_c']:
    dfCoeff.loc[:, metId] -= gam_val
for metId in ['adp_c', 'pi_c', 'h_c']:
    dfCoeff.loc[:, metId] += gam_val
    
# Rearrange
atpm = ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']
biomrxn = model.reactions.get_by_id(biomId)
mets = atpm + sorted([met.id for met in biomrxn.metabolites.keys() if met.id not in atpm])
dfCoeff = dfCoeff.reindex(['biomMW'] + mets, axis=1)

# Save
fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfCoeff')
dfCoeff.to_csv(dir_case + fname + '.csv', sep='\t')


# In[ ]:


fname = file_name(coeff_std=coeff_std, gam_std=None)
dir_pfba = dir_case + 'pFBA/' + fname + '/'
make_pfba_csv(model, dfCoeff, dir_pfba, biomId)


# In[ ]:


fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfCoeff')
dfCoeff = pd.read_csv(dir_case + fname + '.csv', sep='\t', index_col=0, header=0)

fname = file_name(coeff_std=coeff_std, gam_std=None)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfFlux')
index = dfCoeff.index.to_list()
dfFlux = make_fluxes_dataframe(index, dir_pfba)
dfFlux = dfFlux.loc[:, (dfFlux.abs() > zero_tol).any(axis=0)]
dfFlux.to_csv(dir_case + fname + '.csv', sep='\t')

# Remove pFBA files
import shutil
fname = file_name(coeff_std=coeff_std, gam_std=None)
shutil.rmtree(dir_pfba)


# In[ ]:


coeff_std = 0.3 # Set this to 0.05, 0.2, and 0.3 and rerun
size = 10000


# In[ ]:


df_biom = pd.read_excel('Biom_frac_annotated.xlsx', sheet_name='biom')
macro_frac = dict()
for macro in set(df_biom.type):
    macro_frac[macro] = round(df_biom[df_biom.type == macro].biom_frac.sum(), 9)

cols = ['Sum_frac', 'Amino acids', 'Cell wall components', 'DNA', 'RNA',
        'Lipids', 'Cofactors and prosthetic groups', 'Inorganic ions']
df_macro = pd.DataFrame(columns=cols, index=range(0, size+1))
i = 0
for k,v in macro_frac.items():
    df_macro.loc[i, k] = v

seed = -1
for col in df_macro.columns[1:]:
    c0 = df_macro.loc[0, col]
    seed += 1
    np.random.seed(seed)
    coeff_norm = c0 + coeff_std*c0*normal(size=size)
    # Only allow positive value for biom fraction
    coeff_norm = [coeff if coeff > zero_tol else zero_tol for coeff in coeff_norm]
    df_macro.loc[range(1,size+1), col] = coeff_norm
    
df_macro['Sum_frac'] = df_macro.iloc[:, 1:].sum(axis=1)

# Normalize
df_macro.iloc[:, 1:] = df_macro.iloc[:, 1:].div(df_macro.Sum_frac, axis=0)


# In[ ]:


for i in range(1, size+1):
    for macro,mets in met_in_macro.items():
        dfCoeff.loc[i,mets] = dfCoeff.loc[0,mets] * df_macro.loc[i, macro] / df_macro.loc[0, macro]
        
for i in range(1, size+1):
    coeff = 0
    for k,v in aas.items():
        coeff += v*dfCoeff.loc[i, k]
    dfCoeff.loc[i, 'h2o_c'] = -coeff

    coeff = 0
    dfCoeff.loc[i, 'ppi_c'] = 0
    for k,v in nus.items():
        coeff += v*dfCoeff.loc[i, k]
    dfCoeff.loc[i, 'ppi_c'] = -coeff
    
mets_mw = [model.metabolites.get_by_id(i).formula_weight for i in dfCoeff.columns[1:]]
dfCoeff['biomMW'] = -dfCoeff.iloc[:, 1:].multiply(mets_mw, axis=1).sum(axis=1)

# Add GAM back
for col in ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']:
    if col not in dfCoeff.columns:
        dfCoeff[col] = [0.]*dfCoeff.shape[0] # Create array of zeros
for metId in ['atp_c', 'h2o_c']:
    dfCoeff.loc[:, metId] -= gam_val
for metId in ['adp_c', 'pi_c', 'h_c']:
    dfCoeff.loc[:, metId] += gam_val
    
# Rearrange
atpm = ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']
biomrxn = model.reactions.get_by_id(biomId)
mets = atpm + sorted([met.id for met in biomrxn.metabolites.keys() if met.id not in atpm])
dfCoeff = dfCoeff.reindex(['biomMW'] + mets, axis=1)

# Save
fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfCoeff')
dfCoeff.to_csv(dir_case + fname + '.csv', sep='\t')


# In[ ]:


fname = file_name(coeff_std=coeff_std, gam_std=None)
dir_pfba = dir_case + 'pFBA/' + fname + '/'
make_pfba_csv(model, dfCoeff, dir_pfba, biomId)


# In[ ]:


fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfCoeff')
dfCoeff = pd.read_csv(dir_case + fname + '.csv', sep='\t', index_col=0, header=0)

fname = file_name(coeff_std=coeff_std, gam_std=None)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfFlux')
index = dfCoeff.index.to_list()
dfFlux = make_fluxes_dataframe(index, dir_pfba)
dfFlux = dfFlux.loc[:, (dfFlux.abs() > zero_tol).any(axis=0)]
dfFlux.to_csv(dir_case + fname + '.csv', sep='\t')

# Remove pFBA files
import shutil
fname = file_name(coeff_std=coeff_std, gam_std=None)
shutil.rmtree(dir_pfba)

