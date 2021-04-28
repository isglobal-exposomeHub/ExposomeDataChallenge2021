### ISGlobal Exposome Challenge 2021 - pathway analysis and visualization ###
### This script conduct gene set enrichment analysis for DE genes and mediator genes
### Author: Dillon Lloyd


library(tidyverse)
library(XGR)

DATA = "~/NIEHS/Exposome_Data_Challenge/DATA/"
CODE = "~/NIEHS/Exposome_Data_Challenge/CODE/"
RESULT = "~/NIEHS/Exposome_Data_Challenge/RESULT/"


## Input Data ##
load(paste0(RESULT,"DEGS_No_PCS.rdata"))
load(paste0(DATA, "genexpr.rdata"))

med_genes <- readRDS(paste0(RESULT, "gene_med.rds"))
med_genes_ers <- readRDS(paste0(RESULT, "gene_med_ERS.rds"))

## For DE Genes ##
gene_names <- genexpr@featureData@data %>% 
  dplyr::filter(CallRate >= 50) %>% 
  dplyr::select(transcript_cluster_id,GeneSymbol_Affy) %>% 
  separate(GeneSymbol_Affy,c('Gene_Name','Second_Name'), sep = ';') %>% 
  dplyr::filter(is.na(Second_Name)) %>% 
  dplyr::select(-Second_Name)

contrasts <- degs %>% mutate(Gene = toupper(sapply(Gene_ID, function(x) unlist(str_split(x, "_"))[1])),
                                 contrast = 'ast_vs_control', trt = 'ast_vs_control') %>% rename(padj = adj.P.Val,log2FoldChange = logFC)
background <- unique(gene_names$Gene_Name)

contrasts %>% dplyr::filter(padj < .1) %>%
  dplyr::select(Gene, padj) %>% mutate(symbol = toupper(Gene)) %>%
  pull(symbol) -> mysymbol


# REACTOME pathway
myxgr <- xEnricherGenes(data = mysymbol,
                        ontology = 'MsigdbC2REACTOME', 
                        background = background)

# GO BP 
myxgr <- xEnricherGenes(data = mysymbol,
                        ontology = 'GOBP', 
                        background = background)

# GO CC
myxgr <- xEnricherGenes(data = mysymbol,
                        ontology = 'GOCC', 
                        background = background)
                        
# GO MF
myxgr <- xEnricherGenes(data = mysymbol,
                        ontology = 'GOMF', 
                        background = background)

# output results
pw_data <- xEnrichViewer(myxgr, top_num = 250, sortBy = "adjp", details = T) %>%
  mutate(Contrast = names(myxgr)[1]) %>% 
  dplyr::select(Contrast, PW = name, PW.genes = nAnno, 
                Enr.Genes = nOverlap, FDR = adjp, 
                Overlap_Genes = members_Overlap) %>% 
  mutate(GeneRatio = Enr.Genes/PW.genes)


## Visualization ##
# bar plot
p <- xEnrichBarplot(myxgr, displayBy="fdr", FDR.cutoff = as.numeric(.1), wrap.width = 45, top = 50,bar.color = 'green-purple') 
print(p)

# dot plot
p <- pw_data %>% 
  dplyr::mutate(sig = as.factor(ifelse(FDR <= .10,1,0))) %>% 
  dplyr::filter(sig ==1) %>% 
  ggplot(.) +
  geom_point(aes(x = reorder(PW,GeneRatio), y = GeneRatio,color = FDR,size = Enr.Genes,label = PW)) +
  scale_color_gradientn(colours = terrain.colors(10)) +
  ggtitle('DE Pathways') +
  xlab('Pathways') +
  theme_classic() +
  coord_flip()
print(p)

# network plot
network_input <- contrasts %>% dplyr::filter(padj < .1) %>%
  mutate(symbol = toupper(Gene)) %>% 
  dplyr::select(symbol, padj) 

network <- xSubneterGenes(data = network_input)

subg <- network
pattern <- -log10(as.numeric(V(network)$significance))
xVisNet(g=subg, pattern=pattern, vertex.shape="sphere", vertex.label.font=2, newpage=F)


## For Mediator Genes ##
med_gene_list <- paste0(med_genes$Gene_ID,collapse = '|')
string_vec <- str_detect(pw_data$Overlap_Genes, med_gene_list)
med_pws <- pw_data[string_vec,]

med_gene_list_ers <- paste0(med_genes_ers$Gene_ID,collapse = '|')
string_vec <- str_detect(pw_data$Overlap_Genes, med_gene_list_ers)
med_pws_ers <- pw_data[string_vec,]

med_pws %>% 
  dplyr::mutate(sig = as.factor(ifelse(FDR <= .10,1,0))) %>% 
  ggplot(.) +
  geom_point(aes(x = reorder(PW,GeneRatio), y = GeneRatio,color = FDR,size = Enr.Genes,label = PW)) +
  scale_color_gradientn(colours = terrain.colors(10)) +
  ggtitle('Pathways Containing proposed Mediator Genes') +
  xlab('Pathways') +
  theme_classic() +
  coord_flip()

med_pws_ers %>% 
  dplyr::mutate(sig = as.factor(ifelse(FDR <= .10,1,0))) %>% 
  ggplot(.) +
  geom_point(aes(x = reorder(PW,GeneRatio), y = GeneRatio,color = FDR,size = Enr.Genes,label = PW)) +
  scale_color_gradientn(colours = terrain.colors(10)) +
  ggtitle('Pathways Containing proposed Mediator Genes w/ Risk Score') +
  xlab('Pathways') +
  theme_classic() +
  coord_flip()
