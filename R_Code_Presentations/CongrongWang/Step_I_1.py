# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 15:24:52 2021

@author: lucp10491
"""

from mofapy2.run.entry_point import entry_point
import pyreadr
import numpy as np
import pandas as pd
 
# generates list of 4 float values
seeds = np.random.random_integers(1, 1000, size = 50)
ELBO_list = []
Nfactor = []

# Use a list of numpy array
data_dt = pyreadr.read_r("OmicsMatrices_Mval_0.05filtered.RData") # also works for Rds
keys = list(data_dt.keys())


for s in seeds:

    # initialise the entry point
    ent = entry_point()
    
    data_mat = [[None for g in range(1)] for m in range(4)]
    for m in range(4):
        for g in range(1):
            data_mat[m][g] = data_dt[keys[m]].transpose().to_numpy()

    # Set data options
    ent.set_data_options(
        scale_groups = False, 
        scale_views = True
        )

    # Add the data
    ent.set_data_matrix(data_mat, likelihoods = ["gaussian"]*4)# Set model options

    # Set model options
    ent.set_model_options(
        factors = 10, 
        spikeslab_weights = True, 
        ard_factors = True,
        ard_weights = True
        )

    # Set training options
    ent.set_train_options(
        iter = 1000, 
        convergence_mode = "fast", # "slow", "medium" or "fast", corresponding to 1e-5%, 1e-4% or 1e-3% deltaELBO change
        startELBO = 1, 
        freqELBO = 1, 
        dropR2 = 0.02, 
        gpu_mode = True, 
        verbose = False, 
        seed = s  
        )

# We do not want to use stochastic inference for this data set (n < the order of 1e4)
# ent.set_stochastic_options(
# batch_size = 0.5,
# learning_rate = 0.5, 
# forgetting_rate = 0.25
#)

    ent.build()
    ent.run()
    ELBO_list.append(ent.model.calculateELBO()["total"])
    Nfactor.append(ent.model.dim['K'])


res = pd.DataFrame({'Nfactor': Nfactor, 'Seed': seeds, 'ELBO':ELBO_list})
res[res['ELBO']==min(res['ELBO'])]
s = res['Seed'][res['ELBO']==min(res['ELBO'])]

res.to_csv("omics_MOFA_results.csv", index = False)

# Save the output
outfile = "omics_DNAm_0.05filtered.hdf5"
ent.save(outfile)
