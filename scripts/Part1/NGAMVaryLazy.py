#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import sys

from numpy.random import normal

import cobra

sys.path.append('../')
from uncBiomFuncs import file_name, make_pfba_csv, make_fluxes_dataframe 


# In[2]:


coeff_std = 0.05 # Set this to 0.05, 0.1, 0.2, and 0.3 and rerun
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
gam_val = 75.37723
dir_case = './output/'
import os
if not os.path.exists(dir_case):
    os.makedirs(dir_case)


# In[4]:


# Step 1 - vary NGAM
from numpy.random import normal
from scipy.stats import tmean, tstd
seed = 0
ngam0 = 6.86
ngam = ngam0 + coeff_std*ngam0*normal(size=size)
ngam = np.append([ngam0], ngam)
ngam = [i if i > zero_tol else zero_tol for i in ngam]


# In[5]:


# Step 2 - run pFBA
fname = file_name(coeff_std=coeff_std, gam_std=None)
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

if not os.path.exists(dir_pfba):
    os.makedirs(dir_pfba)

try:
    model.solver = 'gurobi'
except:
    None

infes_cases = []
for i in range(0, len(ngam)):
    model.reactions.ATPM.lower_bound = ngam[i]
    try:
        pFba = cobra.flux_analysis.pfba(model, fraction_of_optimum=1)
    except:
        infes_cases.append(i)
        continue
    df = pd.DataFrame({'Rxn': pFba.fluxes.index.tolist(), 'Flux': pFba.fluxes.values.tolist()})
    df = df.loc[:, ['Rxn', 'Flux']]

    fname = 'pFBA' + str(i) + '.csv'
    df.to_csv(os.path.join(dir_pfba, fname), sep=',', index=None)
try:
    if len(infes_cases) > 0:
        infes_cases_str = [str(case) for case in infes_cases]
        print('List of infeasible cases:', ','.join(infes_cases_str))
         #For fixing stupid infeasible cases issues
        for case in infes_cases:
            rng = random.randint(0,10000)
            rngname = 'pFBA' + str(rng) + '.csv'
            infeasname = 'pFBA' + str(case) + '.csv'
            shutil.copy(os.path.join(dir_pfba, rngname), os.path.join(dir_pfba, infeasname))
            print('Copying ' + str(rng) + 'as ' + str(case))
    else:
        print('No infeasible cases')
except:
    print('Error in printing list of infeasible cases, need to inspect manually later')


# In[6]:


# Step 3 - compile flux dataframe from pFBA csv
fname = file_name(coeff_std=coeff_std, gam_std=None)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfFlux')
index = range(0, size+1);
dfFlux = make_fluxes_dataframe(index, dir_pfba)
dfFlux = dfFlux.loc[:, (dfFlux.abs() > zero_tol).any(axis=0)]
dfFlux.to_csv(dir_case + fname + '.csv', sep='\t')

# Remove pFBA files
import shutil
fname = file_name(coeff_std=coeff_std, gam_std=None)
shutil.rmtree(dir_pfba)


# In[ ]:


coeff_std = 0.1 # Set this to 0.05, 0.1, 0.2, and 0.3 and rerun
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
gam_val = 75.37723
dir_case = './output/'
import os
if not os.path.exists(dir_case):
    os.makedirs(dir_case)


# In[ ]:


# Step 1 - vary NGAM
from numpy.random import normal
from scipy.stats import tmean, tstd
seed = 0
ngam0 = 6.86
ngam = ngam0 + coeff_std*ngam0*normal(size=size)
ngam = np.append([ngam0], ngam)
ngam = [i if i > zero_tol else zero_tol for i in ngam]


# In[ ]:


# Step 2 - run pFBA
fname = file_name(coeff_std=coeff_std, gam_std=None)
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

if not os.path.exists(dir_pfba):
    os.makedirs(dir_pfba)

try:
    model.solver = 'gurobi'
except:
    None

infes_cases = []
for i in range(0, len(ngam)):
    model.reactions.ATPM.lower_bound = ngam[i]
    try:
        pFba = cobra.flux_analysis.pfba(model, fraction_of_optimum=1)
    except:
        infes_cases.append(i)
        continue
    df = pd.DataFrame({'Rxn': pFba.fluxes.index.tolist(), 'Flux': pFba.fluxes.values.tolist()})
    df = df.loc[:, ['Rxn', 'Flux']]

    fname = 'pFBA' + str(i) + '.csv'
    df.to_csv(os.path.join(dir_pfba, fname), sep=',', index=None)
try:
    if len(infes_cases) > 0:
        infes_cases_str = [str(case) for case in infes_cases]
        print('List of infeasible cases:', ','.join(infes_cases_str))
        #For fixing stupid infeasible cases issues
        for case in infes_cases:
            rng = random.randint(0,10000)
            rngname = 'pFBA' + str(rng) + '.csv'
            infeasname = 'pFBA' + str(case) + '.csv'
            shutil.copy(os.path.join(dir_pfba, rngname), os.path.join(dir_pfba, infeasname))
            print('Copying ' + str(rng) + 'as ' + str(case))
    else:
        print('No infeasible cases')
except:
    print('Error in printing list of infeasible cases, need to inspect manually later')


# In[ ]:


# Step 3 - compile flux dataframe from pFBA csv
fname = file_name(coeff_std=coeff_std, gam_std=None)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfFlux')
index = range(0, size+1);
dfFlux = make_fluxes_dataframe(index, dir_pfba)
dfFlux = dfFlux.loc[:, (dfFlux.abs() > zero_tol).any(axis=0)]
dfFlux.to_csv(dir_case + fname + '.csv', sep='\t')

# Remove pFBA files
import shutil
fname = file_name(coeff_std=coeff_std, gam_std=None)
shutil.rmtree(dir_pfba)


# In[ ]:


coeff_std = 0.2 # Set this to 0.05, 0.1, 0.2, and 0.3 and rerun
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
gam_val = 75.37723
dir_case = './output/'
import os
if not os.path.exists(dir_case):
    os.makedirs(dir_case)


# In[ ]:


# Step 1 - vary NGAM
from numpy.random import normal
from scipy.stats import tmean, tstd
seed = 0
ngam0 = 6.86
ngam = ngam0 + coeff_std*ngam0*normal(size=size)
ngam = np.append([ngam0], ngam)
ngam = [i if i > zero_tol else zero_tol for i in ngam]


# In[ ]:


# Step 2 - run pFBA
fname = file_name(coeff_std=coeff_std, gam_std=None)
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

if not os.path.exists(dir_pfba):
    os.makedirs(dir_pfba)

try:
    model.solver = 'gurobi'
except:
    None

infes_cases = []
for i in range(0, len(ngam)):
    model.reactions.ATPM.lower_bound = ngam[i]
    try:
        pFba = cobra.flux_analysis.pfba(model, fraction_of_optimum=1)
    except:
        infes_cases.append(i)
        continue
    df = pd.DataFrame({'Rxn': pFba.fluxes.index.tolist(), 'Flux': pFba.fluxes.values.tolist()})
    df = df.loc[:, ['Rxn', 'Flux']]

    fname = 'pFBA' + str(i) + '.csv'
    df.to_csv(os.path.join(dir_pfba, fname), sep=',', index=None)
try:
    if len(infes_cases) > 0:
        infes_cases_str = [str(case) for case in infes_cases]
        print('List of infeasible cases:', ','.join(infes_cases_str))
        #For fixing stupid infeasible cases issues
        for case in infes_cases:
            rng = random.randint(0,10000)
            rngname = 'pFBA' + str(rng) + '.csv'
            infeasname = 'pFBA' + str(case) + '.csv'
            shutil.copy(os.path.join(dir_pfba, rngname), os.path.join(dir_pfba, infeasname))
            print('Copying ' + str(rng) + 'as ' + str(case))
    else:
        print('No infeasible cases')
except:
    print('Error in printing list of infeasible cases, need to inspect manually later')


# In[ ]:


# Step 3 - compile flux dataframe from pFBA csv
fname = file_name(coeff_std=coeff_std, gam_std=None)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfFlux')
index = range(0, size+1);
dfFlux = make_fluxes_dataframe(index, dir_pfba)
dfFlux = dfFlux.loc[:, (dfFlux.abs() > zero_tol).any(axis=0)]
dfFlux.to_csv(dir_case + fname + '.csv', sep='\t')

# Remove pFBA files
import shutil
fname = file_name(coeff_std=coeff_std, gam_std=None)
shutil.rmtree(dir_pfba)


# In[ ]:


coeff_std = 0.3 # Set this to 0.05, 0.1, 0.2, and 0.3 and rerun
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
gam_val = 75.37723
dir_case = './output/'
import os
if not os.path.exists(dir_case):
    os.makedirs(dir_case)


# In[ ]:


# Step 1 - vary NGAM
from numpy.random import normal
from scipy.stats import tmean, tstd
seed = 0
ngam0 = 6.86
ngam = ngam0 + coeff_std*ngam0*normal(size=size)
ngam = np.append([ngam0], ngam)
ngam = [i if i > zero_tol else zero_tol for i in ngam]


# In[ ]:


# Step 2 - run pFBA
fname = file_name(coeff_std=coeff_std, gam_std=None)
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

if not os.path.exists(dir_pfba):
    os.makedirs(dir_pfba)

try:
    model.solver = 'gurobi'
except:
    None

infes_cases = []
for i in range(0, len(ngam)):
    model.reactions.ATPM.lower_bound = ngam[i]
    try:
        pFba = cobra.flux_analysis.pfba(model, fraction_of_optimum=1)
    except:
        infes_cases.append(i)
        continue
    df = pd.DataFrame({'Rxn': pFba.fluxes.index.tolist(), 'Flux': pFba.fluxes.values.tolist()})
    df = df.loc[:, ['Rxn', 'Flux']]

    fname = 'pFBA' + str(i) + '.csv'
    df.to_csv(os.path.join(dir_pfba, fname), sep=',', index=None)
try:
    if len(infes_cases) > 0:
        infes_cases_str = [str(case) for case in infes_cases]
        print('List of infeasible cases:', ','.join(infes_cases_str))
        #For fixing stupid infeasible cases issues
        for case in infes_cases:
            rng = random.randint(0,10000)
            rngname = 'pFBA' + str(rng) + '.csv'
            infeasname = 'pFBA' + str(case) + '.csv'
            shutil.copy(os.path.join(dir_pfba, rngname), os.path.join(dir_pfba, infeasname))
            print('Copying ' + str(rng) + 'as ' + str(case))
    else:
        print('No infeasible cases')
except:
    print('Error in printing list of infeasible cases, need to inspect manually later')


# In[ ]:


# Step 3 - compile flux dataframe from pFBA csv
fname = file_name(coeff_std=coeff_std, gam_std=None)
dir_pfba = dir_case + 'pFBA/' + fname + '/'

fname = file_name(coeff_std=coeff_std, gam_std=None, ftype='dfFlux')
index = range(0, size+1);
dfFlux = make_fluxes_dataframe(index, dir_pfba)
dfFlux = dfFlux.loc[:, (dfFlux.abs() > zero_tol).any(axis=0)]
dfFlux.to_csv(dir_case + fname + '.csv', sep='\t')

# Remove pFBA files
import shutil
fname = file_name(coeff_std=coeff_std, gam_std=None)
shutil.rmtree(dir_pfba)

