#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import sys

from numpy.random import normal

import cobra

sys.path.append('../')
from uncBiomFuncs import file_name, make_normal_distributed_gam, make_pfba_csv, make_fluxes_dataframe 


# In[2]:


gam_std = 0.05 # Set this to 0.05, 0.1, 0.2, and 0.3 and rerun
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
model.reactions.arm_POR5.lower_bound = 0
#model.reactions.POR5.lower_bound = 0

config = cobra.Configuration()
config.tolerance = 1e-9
zero_tol = 1e-9

biomId = 'BIOMASS_Ec_iML1515_core_75p37M'
#biomId = 'BIOMASS_Ec_iML1515_WT_75p37M'
gam_val = 75.37723
dir_case = './output/'
import os
if not os.path.exists(dir_case):
    os.makedirs(dir_case)


# In[4]:


# Step 1 - vary GAM
fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfCoeff')
dfCoeff = make_normal_distributed_gam(model, gam_val, gam_std=gam_std,
                biomId=biomId, size=size)
dfCoeff.to_csv(dir_case + fname + '.csv', sep='\t')


# In[5]:


# Step 2 - run pFBA
fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfCoeff')
dfCoeff = pd.read_csv(dir_case + fname + '.csv', sep='\t', index_col=0, header=0)

fname = file_name(coeff_std=None, gam_std=gam_std)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

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
model.reactions.arm_POR5.lower_bound = 0
#model.reactions.POR5.lower_bound = 0

config = cobra.Configuration()
config.tolerance = 1e-9
zero_tol = 1e-9

make_pfba_csv(model, dfCoeff, dir_pfba, biomId)


# In[6]:


# Step 3 - compile flux dataframe from pFBA csv
fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfCoeff')
dfCoeff = pd.read_csv(dir_case + fname + '.csv', sep='\t', index_col=0, header=0)

fname = file_name(coeff_std=None, gam_std=gam_std)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfFlux')
index = dfCoeff.index.to_list()
dfFlux = make_fluxes_dataframe(index, dir_pfba)
dfFlux = dfFlux.loc[:, (dfFlux.abs() > zero_tol).any(axis=0)]
dfFlux.to_csv(dir_case + fname + '.csv', sep='\t')

# Remove pFBA files
import shutil
fname = file_name(gam_std=gam_std, coeff_std=None)
shutil.rmtree(dir_pfba)


# In[ ]:


gam_std = 0.1 # Set this to 0.05, 0.1, 0.2, and 0.3 and rerun
size = 10000


# In[ ]:


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
model.reactions.arm_POR5.lower_bound = 0
#model.reactions.POR5.lower_bound = 0

config = cobra.Configuration()
config.tolerance = 1e-9
zero_tol = 1e-9

biomId = 'BIOMASS_Ec_iML1515_core_75p37M'
#biomId = 'BIOMASS_Ec_iML1515_WT_75p37M'
gam_val = 75.37723
dir_case = './output/'
import os
if not os.path.exists(dir_case):
    os.makedirs(dir_case)


# In[ ]:


# Step 1 - vary GAM
fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfCoeff')
dfCoeff = make_normal_distributed_gam(model, gam_val, gam_std=gam_std,
                biomId=biomId, size=size)
dfCoeff.to_csv(dir_case + fname + '.csv', sep='\t')


# In[ ]:


# Step 2 - run pFBA
fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfCoeff')
dfCoeff = pd.read_csv(dir_case + fname + '.csv', sep='\t', index_col=0, header=0)

fname = file_name(coeff_std=None, gam_std=gam_std)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

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
model.reactions.arm_POR5.lower_bound = 0
#model.reactions.POR5.lower_bound = 0

config = cobra.Configuration()
config.tolerance = 1e-9
zero_tol = 1e-9

make_pfba_csv(model, dfCoeff, dir_pfba, biomId)


# In[ ]:


# Step 3 - compile flux dataframe from pFBA csv
fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfCoeff')
dfCoeff = pd.read_csv(dir_case + fname + '.csv', sep='\t', index_col=0, header=0)

fname = file_name(coeff_std=None, gam_std=gam_std)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfFlux')
index = dfCoeff.index.to_list()
dfFlux = make_fluxes_dataframe(index, dir_pfba)
dfFlux = dfFlux.loc[:, (dfFlux.abs() > zero_tol).any(axis=0)]
dfFlux.to_csv(dir_case + fname + '.csv', sep='\t')

# Remove pFBA files
import shutil
fname = file_name(gam_std=gam_std, coeff_std=None)
shutil.rmtree(dir_pfba)


# In[ ]:


gam_std = 0.2 # Set this to 0.05, 0.1, 0.2, and 0.3 and rerun
size = 10000


# In[ ]:


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
model.reactions.arm_POR5.lower_bound = 0
#model.reactions.POR5.lower_bound = 0

config = cobra.Configuration()
config.tolerance = 1e-9
zero_tol = 1e-9

biomId = 'BIOMASS_Ec_iML1515_core_75p37M'
#biomId = 'BIOMASS_Ec_iML1515_WT_75p37M'
gam_val = 75.37723
dir_case = './output/'
import os
if not os.path.exists(dir_case):
    os.makedirs(dir_case)


# In[ ]:


# Step 1 - vary GAM
fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfCoeff')
dfCoeff = make_normal_distributed_gam(model, gam_val, gam_std=gam_std,
                biomId=biomId, size=size)
dfCoeff.to_csv(dir_case + fname + '.csv', sep='\t')


# In[ ]:


# Step 2 - run pFBA
fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfCoeff')
dfCoeff = pd.read_csv(dir_case + fname + '.csv', sep='\t', index_col=0, header=0)

fname = file_name(coeff_std=None, gam_std=gam_std)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

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
model.reactions.arm_POR5.lower_bound = 0
#model.reactions.POR5.lower_bound = 0

config = cobra.Configuration()
config.tolerance = 1e-9
zero_tol = 1e-9

make_pfba_csv(model, dfCoeff, dir_pfba, biomId)


# In[ ]:


# Step 3 - compile flux dataframe from pFBA csv
fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfCoeff')
dfCoeff = pd.read_csv(dir_case + fname + '.csv', sep='\t', index_col=0, header=0)

fname = file_name(coeff_std=None, gam_std=gam_std)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfFlux')
index = dfCoeff.index.to_list()
dfFlux = make_fluxes_dataframe(index, dir_pfba)
dfFlux = dfFlux.loc[:, (dfFlux.abs() > zero_tol).any(axis=0)]
dfFlux.to_csv(dir_case + fname + '.csv', sep='\t')

# Remove pFBA files
import shutil
fname = file_name(gam_std=gam_std, coeff_std=None)
shutil.rmtree(dir_pfba)


# In[ ]:


gam_std = 0.3 # Set this to 0.05, 0.1, 0.2, and 0.3 and rerun
size = 10000


# In[ ]:


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
model.reactions.arm_POR5.lower_bound = 0
#model.reactions.POR5.lower_bound = 0

config = cobra.Configuration()
config.tolerance = 1e-9
zero_tol = 1e-9

biomId = 'BIOMASS_Ec_iML1515_core_75p37M'
#biomId = 'BIOMASS_Ec_iML1515_WT_75p37M'
gam_val = 75.37723
dir_case = './output/'
import os
if not os.path.exists(dir_case):
    os.makedirs(dir_case)


# In[ ]:


# Step 1 - vary GAM
fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfCoeff')
dfCoeff = make_normal_distributed_gam(model, gam_val, gam_std=gam_std,
                biomId=biomId, size=size)
dfCoeff.to_csv(dir_case + fname + '.csv', sep='\t')


# In[ ]:


# Step 2 - run pFBA
fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfCoeff')
dfCoeff = pd.read_csv(dir_case + fname + '.csv', sep='\t', index_col=0, header=0)

fname = file_name(coeff_std=None, gam_std=gam_std)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

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
model.reactions.arm_POR5.lower_bound = 0
#model.reactions.POR5.lower_bound = 0

config = cobra.Configuration()
config.tolerance = 1e-9
zero_tol = 1e-9

make_pfba_csv(model, dfCoeff, dir_pfba, biomId)


# In[ ]:


# Step 3 - compile flux dataframe from pFBA csv
fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfCoeff')
dfCoeff = pd.read_csv(dir_case + fname + '.csv', sep='\t', index_col=0, header=0)

fname = file_name(coeff_std=None, gam_std=gam_std)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

fname = file_name(coeff_std=None, gam_std=gam_std, ftype='dfFlux')
index = dfCoeff.index.to_list()
dfFlux = make_fluxes_dataframe(index, dir_pfba)
dfFlux = dfFlux.loc[:, (dfFlux.abs() > zero_tol).any(axis=0)]
dfFlux.to_csv(dir_case + fname + '.csv', sep='\t')

# Remove pFBA files
import shutil
fname = file_name(gam_std=gam_std, coeff_std=None)
shutil.rmtree(dir_pfba)

