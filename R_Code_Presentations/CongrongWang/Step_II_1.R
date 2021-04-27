########## Data processing for  ##########
########## MOFA on the Exposome ##########
# Aim: divide the exposome matrix into 4 sub-matrices

library(tidyverse)
library(fastDummies)
library(MOFA2)
load('../Data/exposome_NA.RData')

# ----- Obtain variable names based on codebook -----

# Check codebook:
table(codebook$domain)
table(codebook$var_type)
table(codebook$period)
table(codebook$transformation)

expo.domain <- names(table(codebook$domain))[c(1,3:5)]

# Names of exposure variables that are categorical:
fac.expo <- codebook %>% 
  filter((var_type == "factor" | transformation == "Dichotomous") # some variables are "numerica" in var_type but "Dichotomous" in transformation
          & domain %in% expo.domain ) %>%
  pull(variable_name)

# Names of exposure variables that are numeric:
num.expo <- codebook %>% 
  filter(!variable_name %in% fac.expo & domain %in% expo.domain) %>%
  pull(variable_name)

# Names of exposure variables that are measured during pregnancy:
pre.expo <- codebook %>%
  filter(period == "Pregnancy" & domain %in% expo.domain) %>%
  pull(variable_name)

# Names of exposure variables that are measured in postnatal period:
post.expo <- codebook %>%
  filter(period == "Postnatal" & domain %in% expo.domain) %>%
  pull(variable_name)


# ----- Create the matrics -----

##### Prenatal numeric exposures #####
exposome_pre_num <- exposomeNA %>% 
  remove_rownames() %>%
  select("ID", intersect(pre.expo, num.expo)) %>%
  column_to_rownames(var = "ID") %>% t()

##### Postnatal numeric exposures #####
exposome_post_num <- exposomeNA %>% 
  remove_rownames() %>%
  select("ID", intersect(post.expo, num.expo)) %>%
  column_to_rownames(var = "ID") %>% t()

##### Prenatal categorical exposures #####
exposome_pre_cat <- exposomeNA %>% 
  select("ID", intersect(pre.expo, fac.expo))
# Check the factor levels, change the reference level if needed (easier for interpretation)
exposome_pre_cat[2] <- as.factor(exposome_pre_cat[2])
lapply(exposome_pre_cat[-1], function(x) levels(as.factor(x)))
exposome_pre_cat$h_pavig_t3_None <- relevel(exposome_pre_cat$h_pavig_t3_None, 
                                            ref = "Low")
exposome_pre_cat$hs_tl_mdich_None <- relevel(exposome_pre_cat$hs_tl_mdich_None, 
                                             ref = "Undetected")
# Create dummy variables:
exposome_pre_bin <- exposome_pre_cat %>%
  dummy_cols(remove_first_dummy = TRUE, remove_selected_columns = TRUE) %>%
  remove_rownames() %>%
  column_to_rownames(var = "ID") %>% t()
  
##### Postnatal categorical exposures #####
exposome_post_cat <- exposomeNA %>% 
  select("ID", intersect(post.expo, fac.expo))
# Check the factor levels, change the reference level if needed (easier for interpretation)
exposome_post_cat[2:3] <- lapply(exposome_post_cat[2:3], as.factor)
lapply(exposome_post_cat[-1], levels)
exposome_post_cat[c("hs_tl_cdich_None", "hs_dmdtp_cdich_None", "hs_cotinine_cdich_None")] <- 
  lapply(exposome_post_cat[c("hs_tl_cdich_None", "hs_dmdtp_cdich_None", "hs_cotinine_cdich_None")],
  function(x) relevel(x, ref = "Undetected"))
exposome_post_cat$hs_contactfam_3cat_num_None <- relevel(exposome_post_cat$hs_contactfam_3cat_num_None, 
                                                 ref = "Less than once a week")
exposome_post_cat$hs_globalexp2_None <- relevel(exposome_post_cat$hs_globalexp2_None,
                                                ref = "no exposure")
exposome_post_cat$hs_smk_parents_None <- relevel(exposome_post_cat$hs_smk_parents_None,
                                                ref = "neither")
# Create dummy variables:
exposome_post_bin <- exposome_post_cat %>%
  dummy_cols(remove_first_dummy = TRUE, remove_selected_columns = TRUE) %>%
  remove_rownames() %>%
  column_to_rownames(var = "ID") %>% t()


# ----- Make a long data frame of the matrices -----

matrix.list <- list(exposome_pre_num, exposome_pre_bin,
                    exposome_post_num, exposome_post_bin)
names(matrix.list) <- c('exposome_pre_num', 'exposome_pre_bin',
                        'exposome_post_num', 'exposome_post_bin') 
for (i in 1:length(matrix.list)) {
  dimnames.m <- dimnames(matrix.list[[i]])
  matrix.list[[i]] <- apply(matrix.list[[i]], 2, as.numeric)
  dimnames(matrix.list[[i]]) <- dimnames.m
  matrix.list[[i]] <- data.frame(reshape2::melt(matrix.list[[i]], 
                                                varnames = c("feature", "sample")),  
                                 view = names(matrix.list)[i],
                                 group = "single_group")
}

matrix.long.df <- do.call(rbind, matrix.list)
rownames(matrix.long.df) <- NULL
write.table(matrix.long.df, "exposome.long.df.txt", row.names = F, sep = "\t")

# ---> go to python and use "mofapy2" with the long df there. 
# ---> python program: Step_II_1.py
# ---> save the model as "exposome.hdf5"
