# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 11:43:22 2021

@author: lucp10491
"""

from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np


data_dt = pd.read_csv("exposome.long.df.txt", sep ="\t")
data_dt.head()

# generates list of 4 float values
seeds = np.random.random_integers(1, 1000, size = 100)
ELBO_list = []
Nfactor = []

for s in seeds:
    
    # initialise the entry point
    ent = entry_point()

    # Set data options
    ent.set_data_options(
        scale_groups = False, 
        scale_views = True
        )
    # Add the data
    ent.set_data_df(data_dt, likelihoods = ["bernoulli", "gaussian"] * 2)
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
        convergence_mode = "fast", 
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

res.to_csv("exposome_MOFA_results.csv", index = False)

# Save the output
outfile = "exposome_4mat.hdf5"
ent.save(outfile)
