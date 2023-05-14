#!/usr/bin/env python
# coding: utf-8

# In[8]:


import cobra, os
import pandas as pd
import numpy as np


# #### Other coeff vary

# In[9]:


path = './output/'

stds_coeff = ['0.05', '0.1', '0.2', '0.3']

dir_model = 'EcoliGecko.json'
model = cobra.io.load_json_model(dir_model)

#Create the stociomatrixes needed for ATP, NADH and NADHP production
ATPStochioDict = {}
for x in model.metabolites.get_by_id("atp_c").reactions:
    looking = model.reactions.get_by_id(x.id)
    for y in looking.metabolites:
        if str(y) == "atp_c":
            ATPStochioDict[x.id] = looking.metabolites[y]

NADHStochioDict = {}
for x in model.metabolites.get_by_id("nadh_c").reactions:
    looking = model.reactions.get_by_id(x.id)
    for y in looking.metabolites:
        if str(y) == "nadh_c":
            NADHStochioDict[x.id] = looking.metabolites[y]            

NADPHStochioDict = {}
for x in model.metabolites.get_by_id("nadph_c").reactions:
    looking = model.reactions.get_by_id(x.id)
    for y in looking.metabolites:
        if str(y) == "nadph_c":
            NADPHStochioDict[x.id] = looking.metabolites[y]            

for std in stds_coeff:

    # Load dataframe
    fname = 'ftype-dfFlux_CoeffStd-' + std + '.csv'
    fpath = os.path.join(path, fname)
    df = pd.read_csv(fpath, sep='\t', index_col=0)
    df = df.fillna(0)

    # Fix
    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['QMO2No1', 'QMO3No1']:
        if rxn in df.columns:
            sumval += df[rxn]
    df['sumQMOx'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['PYKNo1', 'PYKNo2', 'PYK2No1', 'PYK3No1', 'PYK4No1', 'PYK6No1']:
        if rxn in df.columns:
            sumval += df[rxn]
    df['sumPYKx'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['NADH16ppNo1', 'NADH17ppNo1']:
        if rxn in df.columns:
            sumval += df[rxn]
    df['sumNADHxpp'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['arm_FDH4pp', 'arm_FDH5pp']:
        if rxn in df.columns:
            sumval += df[rxn]
    df['sumFDHxpp'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['DHORD2No1', 'DHORD5No1', 'DHORDfum']:
        if rxn in df.columns:
            sumval += df[rxn]
    df['sumDHORDx'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['GCALDDNo1', '4HTHRANo1', '4HTHRA_RevNo1']:
    # Downstream involvement: ['GCALDD', '4HTHRA', 'GLYCTO2',  'MALS', '4HTHRK',
    # 'OHPBAT', 'PERD', 'E4PD', 'GLYCLTDx', 'GLYCLTDy', 'GLYCTO2', 'GLYCTO3', 'GLYCTO4']
        if rxn in df.columns:
            if rxn in ['4HTHRANo1']:
                sumval -= df[rxn]
            else:
                sumval += df[rxn]
    df['sumGCALDdeg'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['GARTNo1', 'GARFTNo1', 'GARFT_REVNo1']:
        if rxn in df.columns:
            sumval += df[rxn]
    df['sumFGAMsyn'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['THRt4pp', 'GLYCLTt4pp', 'PROt4pp', 'GLUt4pp', 'ACt4pp',
                'ALAt4pp', 'GLYt4pp', 'SERt4pp']:
        if rxn in df.columns:
            sumval += df[rxn]
    df['sumSodiumImport'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['THRt2pp', 'GLYCLTt2rpp', 'PROt2rpp', 'GLUt2rpp', 'ACt2rpp',
                'ALAtpp', 'GLYtpp', 'SERtpp']:
        if rxn in df.columns:
            if rxn in ['GLYCLTt2rpp', 'PROt2rpp', 'GLUt2rpp', 'ACt2rpp']:
                sumval -= df[rxn]
            else:
                sumval += df[rxn]
    df['sumSodiumImportBalance'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['ASPO3No1', 'ASPO4No1', 'ASPO5No1', 'ASPO6No1']:
        if rxn in df.columns:
            sumval += df[rxn]
    df['sumASPOx'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['arm_NDPK1', 'PYK3No1']:
        if rxn in df.columns:
            sumval += df[rxn]
    df['sumGTPsyn'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['arm_NDPK2', 'PYK2No1']:
        if rxn in df.columns:
            sumval += df[rxn]
    df['sumUTPsyn'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['arm_NDPK3', 'PYK4No1']:
        if rxn in df.columns:
            sumval += df[rxn]
    df['sumCTPsyn'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['arm_NDPK4', 'PYK6No1']:
        if rxn in df.columns:
            sumval += df[rxn]
    df['sumDTTPsyn'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ['ADK1No1', 'ADK3No1']:
        if rxn in df.columns:
            sumval += df[rxn]
    df['sumAMPkinase'] = sumval
    
    #Add together the fluxes of all reactions capable of producing ATP, multiplied by the stochiomatrix of ATP aslong as this value is positive
    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in ATPStochioDict:
        if rxn in df.columns:
            fluxes = df[rxn]
            scaledfluxes = [flux * ATPStochioDict[rxn] for flux in fluxes]
            if scaledfluxes[1] > 0:
                sumval += scaledfluxes
    df['TotalATPFlux'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in NADHStochioDict:
        if rxn in df.columns:
            fluxes = df[rxn]
            scaledfluxes = [flux * NADHStochioDict[rxn] for flux in fluxes]
            if scaledfluxes[1] > 0:
                sumval += scaledfluxes
    df['TotalNADHFlux'] = sumval

    sumval = pd.Series(0, index=np.arange(0,10001))
    for rxn in NADPHStochioDict:
        if rxn in df.columns:
            fluxes = df[rxn]
            scaledfluxes = [flux * NADPHStochioDict[rxn] for flux in fluxes]
            if scaledfluxes[1] > 0:
                sumval += scaledfluxes
    df['TotalNADPHFlux'] = sumval

    # Trim columns with all zeros (tolerance 10^-8)
    rxns_drop = []
    for rxn in df.columns:
        if all(df[rxn].abs() < 1e-8):
            rxns_drop.append(rxn)
    df = df.drop(rxns_drop, axis=1)

    # Save
    fname = 'ftype-dfFlux_CoeffStd-' + std + '_added.csv'
    fpath = os.path.join(path, fname)
    df.to_csv(fpath, sep='\t', index=True)


# In[ ]:




