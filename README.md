# NOTE for publications using data from this challenge

Overall, users of data are strongly encouraged to publish their results in peer-reviewed journals and to present research findings at scientific meetings, etc. Investigators planning to conduct analyses similar to those described at https://www.projecthelix.eu/ may contact Consortium members at helixdata@isglobal.org  to discuss collaborations, if so desired. The raw data supporting the current study are available on request subject to ethical and legislative review, see details here: https://www.projecthelix.eu/index.php/es/data-inventory 
These data may also be used for educational purposes.

The following text should be added to any publications based on this data:
This data were created as part of the ISGlobal Exposome data challenge 2021, presented in this publication (preprint: https://arxiv.org/abs/2202.01680 - under review in Env. Int.). The HELIX study [Vrijheid, Slama, et al. EHP 2014; Maitre et al. 2018 BMJ Open] represents a collaborative project across six established and ongoing longitudinal population-based birth cohort studies in six European countries (France, Greece, Lithuania, Norway, Spain, and the United Kingdom).  The research leading to these results has received funding from the European Community’s Seventh Framework Programme (FP7/2007-2013) under grant agreement no 308333 – the HELIX project and the H2020-EU.3.1.2. - Preventing Disease Programme under grant agreement no 874583 (ATHLETE project). The data used for the analyses described in this manuscript were obtained from: Figshare https://figshare.com/account/home#/projects/98813 (project number 98813 accesed on MM/DD/YYYY) and github https://github.com/isglobal-exposomeHub/ExposomeDataChallenge2021/.

# Exposome Data Challenge 2021

The exposome, described as "the totality of human environmental exposures from conception onwards", recognizes that individuals are exposed simultaneously to a multitude of different environmental factors and takes a holistic approach to the discovery of etiological factors for disease. The exposome’s main advantage over traditional ‘one-exposure-one-disease’ study approaches is that it provides an unprecedented conceptual framework for the study of multiple environmental hazards (urban, chemical, lifestyle, social) and their combined effects.

The objective of this event (described [here](https://www.isglobal.org/-/exposome-data-analysis-challenge)) is to promote innovative statistical, data science, or other quantitative approaches to studying the health effects of complex high-throughput measurement of exposure indicators (exposome). Detailed challenge examples are given on [this link](https://docs.google.com/document/d/1ul3v-sIniLuTjFB1F1CrFQIX8mrEXVnvSzOF7BCOnpQ/edit). 

These are the availalbe datasets to propose data analyses to address any challenge:

- **Exposome data (n=1301)**:  Rdata file [without missings](https://github.com/isglobal-brge/brge_data_large/blob/master/data/ExposomeDataChallenge2021/exposome.RData) and [with missings](https://github.com/isglobal-brge/brge_data_large/blob/master/data/ExposomeDataChallenge2021/exposome_NA.RData) containing three objects:
     - 1 object for exposures: `exposome`
     - 1 object for covariates: `covariates`
     - 1 object for outcomes: `phenotype`

The three tables can be linked using **ID** variable. See the [codebook](https://github.com/isglobal-brge/brge_data_large/blob/master/data/ExposomeDataChallenge2021/codebook.xlsx) for variable description (variable name, domain, type of variable, transformation, ...)


- **omic data**: Exposome and omic data can be linked using **ID** variable. 
     - [Proteome](https://github.com/isglobal-brge/brge_data_large/blob/master/data/ExposomeDataChallenge2021/proteome.Rdata): ExpressionSet called `metabol_serum` of **1170 individuals** and **39 proteins** (log-transformed) that are annotated in the `ExpressionSet` object (use `fData(proteome)` after loading `Biobase` Bioconductor package).
     - [Serum Metabolome](https://github.com/isglobal-brge/brge_data_large/blob/master/data/ExposomeDataChallenge2021/metabol_serum.Rdata): ExpressionSet called `metabol_serum` of **1198 individuals** and **177 metabolites** (log-transformed) (see [here](https://github.com/isglobal-brge/brge_data_large/blob/master/data/ExposomeDataChallenge2021/HELIX_serum_metabol_report_IC_v4_APS_2017_04_06.pdf) for a descripton).
     - [Urine Metabolome](https://github.com/isglobal-brge/brge_data_large/blob/master/data/ExposomeDataChallenge2021/metabol_urine.Rdata): ExpressionSet called `metabol_urine` of **1192 individuals** and **44 metabolites** (see [here](https://github.com/isglobal-brge/brgedata/blob/master/data/ExposomeDataChallenge2021/HELIX_urine_metabol_report_IC_v3_CHL_2017_01_26.pdf) for a descripton). 
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
 



