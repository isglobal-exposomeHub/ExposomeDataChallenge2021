# Exposome Data Challenge code submission

## Authors / Institutions

Columbia University Mailman School of Public Health

Department of Environmental Health Sciences

All code contained in this directory was written by [Jaime Benavides](https://github.com/jaime-benavides), [Yanelli Nunez](https://github.com/yanellinunez), and [Lawrence Chillrud](https://github.com/lawrence-chillrud), with advice from [Marianthi-Anna Kioumourtzoglou](https://github.com/marianthi) and [Elizabeth Gibson](https://github.com/lizzyagibson).

## Project Overview

Title: Pre- and postnatal urban exposure patterns and childhood neurobehavior

Abstract: Several studies have reported detrimental neurobehavioral effects of urban air pollution in children and protective effects of access to green space. However, the influence of the urban exposome on children’s neurobehavior remains largely unexplored. We aim to identify pre- and postnatal exposure patterns that may be associated with adverse neurobehavior in children. We propose to use an unsupervised dimensionality reduction and pattern recognition algorithm, Principal Component Pursuit (PCP), a robust computer vision algorithm. PCP decomposes the exposure matrix into a low-rank matrix to identify consistent exposure patterns and a sparse matrix to isolate unique exposure events. This feature allows the separation of rare or extreme events from shared patterns in the study population. Thus, patterns are not influenced by outlying values, and unique exposures may still provide useful information for health analyses. The main PCP advantage over traditionally-used factor analytic approaches is that it is not influenced by outlying values and that it can still recover the low-rank matrix even in presence of missingness in the data. Although PCP performs well in computer vision applications, it has not been applied to environmental and urban exposome data that are likely noisier and may lack a prominent low-rank structure. We will subsequently include the PCP-identified patterns in health models and adjust for potential confounders to evaluate their association with neurobehavior.

## Contents

1. [code](code) is a directory containing all `.R` and `.Rmd` scripts used for our submission to the [2021 Exposome Data Challenge](https://github.com/isglobal-exposomeHub/ExposomeDataChallenge2021). Below we provide a brief overview of each file:
    * [0\_00\_run\_analysis.R](code/0_00_run_analysis.R) is the master script of the project. This script calls all other scripts to conduct our analysis, in the sequence they should be called in. In other words, this should be the only file you need to run to conduct our analysis.
    * [0\_01\_init\_directory\_structure.R](code/0_01_init_directory_structure.R) creates all the directories that will be needed for our code to run. This includes an `output` directory and a `generated_data` sub-directory inside the [data](data) directory.
    * Scripts with the prefix `A_0` relate to data preparation, cleaning, wrangling, formatting, etc. along with exploratory data analysis that inform later decisions in our analysis.
    * The script with the prefix `B_0` applies PCP to & subsequently performs pattern recognition on the generated data from the `A_0` family.
    * Scripts with the prefix `C_0` run the health models for this analysis, using the patterns and data outputted from the `B_0` script.
    * Scripts with the prefix `F` serve as helper functions that are called by scripts in the `A_0` and `B_0` families. 
    * Each script is documented to give anyone interested the general idea of the workflow. It is assumed you have watched / heard / read through our presentation when reading the code & its documentation.

2. [data](data) is a directory containing the data we used for our analysis. This directory contains the following sub-directories:
    * [raw\_challenge\_data](data/raw_challenge_data) contains the raw data provided by the 2021 Exposome Data Challenge.
    * `generated_data` is a directory created when running the analysis. It is not included in this repo.

3. [Jaime\_Benavides.Rproj](Jaime_Benavides.Rproj) is the `.Rproj` file for this directory.

4. [R\_session\_info.txt](R_session_info.txt) contains the packages used for this analysis, along with their versions.

## Package Dependencies

Our code makes use of many packages that can be installed from CRAN directly. See [R\_session\_info.txt](R_session_info.txt) for package versions. Here we note two packages not yet available on CRAN that are integral to our analysis:

1. [pcpr](https://github.com/Columbia-PRIME/pcpr) is a package that contains the source code for all PCP functions used throughout our analysis. We are adapting PCP to environmental health data, and have made our repository public as we continue to develop this exciting method. Please note that as this package is still under constant development, any aspects of the package may change in the future, such as documentation, usage, scope, etc. We would appreciate any feedback if you use this package on your own data! We will also be publishing a paper on this package very soon, so please be on the look out for that as well.

2. [PCPhelpers](https://github.com/Columbia-PRIME/PCPhelpers) is a complimentary package to `pcpr` that provides useful functions for testing / working with PCP, such as gridsearches, etc. This package is also under heavy development. Any aspects of the package may change in the future, such as documentation, usage, scope, etc. And again, we would appreciate any feedback if you use this package when working with PCP!

These above highlighted packages should be installed before running our analysis.
