**NOTE for selected talks**: Read [this file](https://github.com/isglobal-exposomeHub/ExposomeDataChallenge2021/blob/main/R_Code_Presentations/README.md) to see how to upload your R code.  

# Exposome Data Challenge 2021

The exposome, described as "the totality of human environmental exposures from conception onwards", recognizes that individuals are exposed simultaneously to a multitude of different environmental factors and takes a holistic approach to the discovery of etiological factors for disease. The exposome’s main advantage over traditional ‘one-exposure-one-disease’ study approaches is that it provides an unprecedented conceptual framework for the study of multiple environmental hazards (urban, chemical, lifestyle, social) and their combined effects.

The objective of this event (described [here](https://www.isglobal.org/-/exposome-data-analysis-challenge)) is to promote innovative statistical, data science, or other quantitative approaches to studying the health effects of complex high-throughput measurement of exposure indicators (exposome). Detailed challenge examples are given on [this link](https://docs.google.com/document/d/1ul3v-sIniLuTjFB1F1CrFQIX8mrEXVnvSzOF7BCOnpQ/edit). 

These are the availalbe datasets to propose data analyses to address any challenge:

- **Exposome data (n=1301)**:  Rdata file [without missings](https://github.com/isglobal-brge/brgedata/blob/master/data/ExposomeDataChallenge2021/exposome.RData) and [with missings](https://github.com/isglobal-brge/brgedata/blob/master/data/ExposomeDataChallenge2021/exposome_NA.RData) containing three objects:
     - 1 object for exposures: `exposome`
     - 1 object for covariates: `covariates`
     - 1 object for outcomes: `phenotype`

The three tables can be linked using **ID** variable. See the [codebook](https://github.com/isglobal-brge/brgedata/blob/master/data/ExposomeDataChallenge2021/codebook.xlsx) for variable description (variable name, domain, type of variable, transformation, ...)


- **omic data**: Exposome and omic data can be linked using **ID** variable. 
     - [Proteome](https://github.com/isglobal-brge/brgedata/blob/master/data/ExposomeDataChallenge2021/proteome.Rdata): ExpressionSet called `metabol_serum` of **1170 individuals** and **39 proteins** (log-transformed) that are annotated in the `ExpressionSet` object (use `fData(proteome)` after loading `Biobase` Bioconductor package).
     - [Serum Metabolome](https://github.com/isglobal-brge/brgedata/blob/master/data/ExposomeDataChallenge2021/metabol_serum.Rdata): ExpressionSet called `metabol_serum` of **1198 individuals** and **177 metabolites** (log-transformed) (see [here](https://github.com/isglobal-brge/brgedata/blob/master/data/ExposomeDataChallenge2021/HELIX_serum_metabol_report_IC_v4_APS_2017_04_06.pdf) for a descripton).
     - [Urine Metabolome](https://github.com/isglobal-brge/brgedata/blob/master/data/ExposomeDataChallenge2021/metabol_urine.Rdata): ExpressionSet called `metabol_urine` of **1192 individuals** and **44 metabolites** (see [here](https://github.com/isglobal-brge/brgedata/blob/master/data/ExposomeDataChallenge2021/HELIX_urine_metabol_report_IC_v3_CHL_2017_01_26.pdf) for a descripton). 
     - [Gene expression](https://figshare.com/s/571c8cff7acf5167f343): ExpressionSet called `genexpr`  (see [here](https://isglobal-brge.github.io/Master_Bioinformatics/bioconductor.html#expressionset) what an ExpressionSet is) of **1007 individuals** and **28,738 transcripts** with annotated gene symbols. 
     - [Methylation](https://figshare.com/s/46e6a1d66ff135bb15c8): GenomicRatioSet called `methy` (see [here](https://www.rdocumentation.org/packages/minfi/versions/1.18.4/topics/GenomicRatioSet-class) what a GenomicRatioSet is) of **918 individuals** and **386,518 CpGs**

The variables that are available in the metadata are:

> 1. **ID**: identification number
> 2. **e3_sex**: gender (male, female)
> 3. **age_sample_years**: age (in years)
> 4. **h_ethnicity_cauc**: caucasic? (yes, no)
> 5. **ethn_PC1**: first PCA to address population stratification
> 6. **ethn_PC2**: second PCA to address population stratification
> 7. **Cell-type estimates** (only for methylation): NK_6, Bcell_6, CD4T_6, CD8T_6, Gran_6, Mono_6
 



