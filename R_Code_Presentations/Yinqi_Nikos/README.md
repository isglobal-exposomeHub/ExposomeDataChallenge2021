#  Latent Unknown Clustering with Integrated Multi-Omics Data (LUCID) 
  
Yinqi Zhao (yinqiz@usc.edu), Nikolaos Stratakis (nstratak@usc.edu)

For this data challenge, we addressed the research question: Does early-life organochlorine (OC) exposure increase obesity risk, and if yes, what are the underlying mechanistic pathways?

We conducted an integrated analysis by implementing a statistical model called LUCID. By integrating exposure profiles to OC, proteomics, serum and urine metabolites, we identified 4 sbugroups characterized by distinguished omics variables and estimated their associations to BMI (measured in z-score). The R script `data_challenge_LUCID.R` contains codes for this analysis. The LUCID model is implemented in an R package called `LUCIDus`. It is availabe on `CRAN` and an updated version can be founc on [Github](https://github.com/USCbiostats/LUCIDus/tree/stable-v1). A copy of the package source file for LUCIDus is also attached in this repo.

For LUCID method paper, please refer to 
+ Peng, Cheng, et al. "A latent unknown clustering integrating multi-omics data (LUCID) with phenotypic traits." Bioinformatics 36.3 (2020): 842-850.

For latest LUCID application, please refer to 
+ Stratakis, Nikos, et al. "Prenatal exposure to perfluoroalkyl substances associated with increased susceptibility to liver injury in children." Hepatology 72.5 (2020): 1758-1770.