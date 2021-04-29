#####################################
##### MultiOmics preprocessing ######
#####################################
rm(list = ls())
library(limma)
library(minfi)

################################
###read the MultiOmics
for (file in list.files(pattern = "*.RData|*.Rdata")) {load(file)}

################################
## Gene expression data
## In this step we filter the transcripts associated (via robust linear regression) with zBMI with a p-value<0.10

genes_expr <- exprs(genexpr) # retrieve expression data from eSet

de.design <- model.matrix(~ 0 + phenotypeNA$hs_zbmi_who[match(genexpr@phenoData@data$ID, phenotypeNA$ID)])
genexpr.fit <- lmFit(genes_expr, design = de.design, method = 'robust') # fit robust regression between transcritome and zBMI
genexpr.fit <- eBayes(genexpr.fit)
genexpr.fit.df <- topTable(genexpr.fit, number = nrow(genes_expr))
genexpr.fit.df <- genexpr.fit.df[genexpr@featureData@data$CallRate > 95,] #select features with call rate >95
tx <- rownames(genexpr.fit.df)[genexpr.fit.df$P.Value < 0.05]
genes_expr_filtered <- genes_expr[tx,]


################################
## DNA methylation data 
## In this step we filter the CpGs associated (via robust linear regression) with zBMI with a p-value<0.10 and we get rid of the the noise of cell type (via linear model)

DNAm <- getM(methy) # retrieve M values data from MethylSet

# Filtering
pheno_methyl <- as.data.frame(pData(methy))
dm.design <- model.matrix(~ 0 + phenotypeNA$hs_zbmi_who[match(colnames(DNAm), phenotypeNA$ID)] +
                            NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6,
                          data = pheno_methyl[match(colnames(DNAm), pheno_methyl$ID),])
methy.fit <- lmFit(DNAm, design = dm.design, method = 'robust')
methy.fit <- eBayes(methy.fit)
methy.fit.df <- topTable(methy.fit, number = nrow(DNAm))
dnam <- rownames(methy.fit.df)[methy.fit.df$P.Value < 0.05] # select features with p val<0.10 and save them
DNAm_filterd <- betas[dnam,] # filter the betas based on the former CpG list

# Denoise each probe for cell type
DNAm_filtered_denoised <- matrix(NA, nrow(DNAm_filterd), ncol(DNAm_filterd),
                                 dimnames = c(list(rownames(DNAm_filterd)), list(colnames(DNAm_filterd))) )
all(rownames(pheno_methyl)==colnames(DNAm_filterd))
for (i in 1:nrow(DNAm_filterd)) {
  DNAm_filtered_denoised[i, ] <- 
    resid(lm(DNAm_filterd[i,] ~ NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 + Mono_6, data = pheno_methyl))
  print(rownames(DNAm_filterd)[i])
}


################################
## Tidy the MultiOmics data 
## In this step we complete each Omic dataset adding a row for each individual present at least in one Omic dataset

# retrieve urine metabolites data
metabolomics_urine<-exprs(metabol_urine) 
metabolomics_serum<-exprs(metabol_serum)

# list of individuals present at least in one omic dataset
ids <- unique(c(colnames(DNAm_filtered_denoised),
                colnames(genes_expr_filtered),
                colnames(metabolomics_serum),
                colnames(metabolomics_urine)))

# list of individuals missing in each omic data compared to "ids"
ids_DNAm <- ids[!(ids%in%colnames(DNAm_filtered_denoised))]
ids_metab_urine <- ids[!(ids%in%colnames(metabolomics_urine))]
ids_metab_serum <- ids[!(ids%in%colnames(metabolomics_serum))]
ids_genes <- ids[!(ids%in%colnames(genes_expr_filtered))]

# generate a matrix filled with NAs for the missing individuals in each omic data
DNAm_add <- matrix(NA, nrow(DNAm_filtered_denoised), length(ids_DNAm), 
                   dimnames = c(list(rownames(DNAm_filtered_denoised)), 
                                list(ids_DNAm)) )
metabolomics_urine_add <- matrix(NA,nrow(metabolomics_urine),length(ids_metab_urine), 
                                 dimnames = c(list(rownames(metabolomics_urine)), 
                                              list(ids_metab_urine)) )
metabolomics_serum_add <- matrix(NA,nrow(metabolomics_serum),length(ids_metab_serum), 
                                 dimnames = c(list(rownames(metabolomics_serum)), 
                                              list(ids_metab_serum)) )
genes_expr_add <- matrix(NA, nrow(genes_expr_filtered), length(ids_genes), 
                         dimnames = c(list(rownames(genes_expr_filtered)), 
                                      list(ids_genes)) )

# add the NA matrix to the original omic data
DNAm_filtered_denoised <- cbind(DNAm_filtered_denoised, DNAm_add)
metabolomics_urine <- cbind(metabolomics_urine, metabolomics_urine_add)
metabolomics_serum <- cbind(metabolomics_serum,metabolomics_serum_add)
genes_expr_filtered <- cbind(genes_expr_filtered, genes_expr_add)

# re-order the datasets
metabolomics_urine <- metabolomics_urine[,ids]
metabolomics_serum <- metabolomics_serum[,ids]
genes_expr_filtered <- genes_expr_filtered[,ids]
DNAm_filtered_denoised <- DNAm_filtered_denoised[,ids]


################################
##Save the tidy MultiOmics data
#In this step we save the data
save(metabolomics_serum, metabolomics_urine, 
     genes_expr_filtered, DNAm_filtered_denoised,
     file = "Data/OmicsMatrices_filtered.RData")

################################
##Run MultiOmics MOFA using python ---> go to Step_I_1.py
#Please, go to Step_I_1.py to run MultiOmics MOFA