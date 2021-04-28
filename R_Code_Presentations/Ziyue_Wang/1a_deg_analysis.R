### ISGlobal Exposome Challenge 2021 - differential expressed gene analysis ###
### This script identify DE genes associated with Asthma
### Author: Dillon Lloyd


library(tidyverse)
library(limma)
library(edgeR)
library(plotly)
library(RColorBrewer)

DATA = "~/NIEHS/Exposome_Data_Challenge/DATA/"
CODE = "~/NIEHS/Exposome_Data_Challenge/CODE/"
RESULT = "~/NIEHS/Exposome_Data_Challenge/RESULT/"


## Input Data ##
load(paste0(DATA, "exposome.rdata"))
load(paste0(DATA, "genexpr.rdata"))

# Create Model Data
pheno <- genexpr@phenoData@data %>% 
  mutate(ID = as.numeric(ID))

outcome <- phenotype %>% 
  dplyr::select(ID,hs_asthma) 

model_data <- left_join(pheno,outcome, by = 'ID') %>% mutate(hs_asthma = as.factor(hs_asthma))


## Preprocess Counts ##
counts <- as.data.frame(genexpr@assayData[["exprs"]])

# Drop multiple gene maps 
gene_names <- genexpr@featureData@data %>% 
  dplyr::filter(CallRate >= 50) %>% 
  dplyr::select(transcript_cluster_id,GeneSymbol_Affy) %>% 
  separate(GeneSymbol_Affy,c('Gene_Name','Second_Name'), sep = ';') %>% 
  dplyr::filter(is.na(Second_Name)) %>% 
  dplyr::select(-Second_Name)

length(unique(gene_names$Gene_Name))

gene_counts <- as.data.frame(genexpr@assayData[["exprs"]]) %>%
  rownames_to_column('transcript_id') %>%
  left_join(.,gene_names, by = c('transcript_id' = 'transcript_cluster_id')) %>%
  dplyr::select(transcript_id,Gene_Name,everything()) %>%
  column_to_rownames('transcript_id') %>%
  dplyr::group_by(Gene_Name) %>%
  dplyr::summarise_all('mean')

gene_counts <- gene_counts %>% dplyr::filter(!is.na(Gene_Name))
save(gene_counts, file = paste0(RESULT,"gene_counts.Rdata"))

# Remove columns with all zeros (None are removed)
load(file = paste0(RESULT,"gene_counts.Rdata"))

allcounts <- gene_counts %>% column_to_rownames('Gene_Name')
qc_counts <-  allcounts[,apply(allcounts,2,function(x) !all(x==0))]


## Fit Limma Model ##
design <- model.matrix(~ model_data$hs_asthma + model_data$e3_sex + model_data$age_sample_years + model_data$h_ethnicity_cauc)
fit <- lmFit(qc_counts,design)
fit <- eBayes(fit)
topTable(fit)

all_data <- topTable(fit, coef=1, number=9999999) %>% 
  rownames_to_column('Gene_ID') %>% 
  mutate(fdr_p = p.adjust(P.Value, method = 'BH'),
         sig = as.factor(ifelse(fdr_p <= 0.05,1,0)))

degs <- all_data %>% dplyr::filter(sig == 1) 

save(degs,file = '../DEGS_No_PCS.Rdata')


## Visualizaton ##
# violin plot of DE genes
p <- ggplot(data=all_data, aes(x=logFC, y=-log10(adj.P.Val), color = sig,id = Gene_ID)) + 
  ggtitle('DEGs for Asthma, FDR = 0.05') + 
  labs(subtitle = paste0('Number of DEGs = ',degs)) +
  xlab('Log 2 Fold Change') +
  ylab('-log10(P Value)') +
  geom_point() +
  scale_color_brewer(palette = 'Dark2') +
  theme_classic() 

ggplotly(p)


# PCA on count file and label color by sex and ashtma 
c <- allcounts
c <- c[,apply(c, 2, var, na.rm=TRUE) != 0]

X <- t(scale(t( c ))) 
Y <- X[complete.cases(X),]
p <- prcomp(Y, center = T)

pc <- unclass(p$rotation)[, 1:5]
p_var <- (round(100 * (p$sdev ^ 2 / sum(p$sdev ^ 2)), 1))[1:5]

trt_data <- model_data %>% dplyr::select(ID,hs_asthma,e3_sex) %>% 
  mutate(ID = as.character(ID))

pc %>% as.data.frame() %>% rownames_to_column(. , var = "Sample") %>% rowwise() %>%
  arrange(Sample) %>% as.data.frame() %>% 
  left_join(.,trt_data, by = c('Sample' = 'ID')) %>% 
  mutate(hs_asthma = as.factor(hs_asthma))-> test

test %>% ggplot(aes(x = PC1, y = PC2,   color = hs_asthma, shape = e3_sex)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = 'Set2') +
  
  labs(title = "PCA of all Samples",
       subtitle = 'Color Indicates Asthma Status and Shape Indicates Sex',
       x = paste("PC1 (", p_var[1], "%)"),
       y = paste("PC2 (", p_var[2], "%)")) + theme_classic() -> p
ggplotly(p)  
