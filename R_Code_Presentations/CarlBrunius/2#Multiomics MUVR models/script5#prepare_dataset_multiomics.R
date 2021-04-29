############################################################
# Exposome Challenge 
# Harmonizing and combining all omics datasets
############################################################
#Combining the results from the sPLS analysis of methylation and gene expression data into
#data frames specific for each exposure and printing them as .csv-files


# multiomics
rm(list=ls())
library(GEOquery)
library(minfi)

load(file="methy.RData")
b <- "getM"(methy)
x<-t(b)
methy<-x[order(as.integer(rownames(x))),]

load(file="metabol_serum.RData")
metabol_serum <- t(exprs(metabol_serum)) 
metabol_serum<-metabol_serum[order(as.integer(rownames(metabol_serum))),]

load(file="metabol_urine.RData")
metabol_urine <- t(exprs(metabol_urine))
metabol_urine<-metabol_urine[order(as.integer(rownames(metabol_urine))),]

load(file="proteome.RData")
proteome <- t(exprs(proteome))
proteome <- proteome[order(as.integer(rownames(proteome))),]


load(file="genexpr.RData")
x<-t(exprs(genexpr))
gene<-x[order(as.integer(rownames(x))),]

all_id <- Reduce(intersect, list(rownames(metabol_serum), rownames(metabol_urine),rownames(proteome),rownames(gene), rownames(methy)))
#877 IDs incommon for all datasets

save(all_id,file="incommon_id.rda")

metabol_serum<-metabol_serum[which (rownames(metabol_serum) %in% all_id),]
metabol_urine <-metabol_urine[which (rownames(metabol_urine) %in% all_id),]
proteome <-proteome[which (rownames(proteome) %in% all_id),]
methy<-methy[which (rownames(methy) %in% all_id),]
gene<-gene[which (rownames(gene) %in% all_id),]


load("selected_exposures_partcor.Rdata") # From the R file partcor_uniquevariables.R

# selected variables from spls
load("gene_whichElected.rda") 
whichElected_gene <- whichElected
load("methy_whichElected_M.rda")
whichElected_methyl <- whichElected


multiOmicsXVars<-list()
for(i in 1:26){
  gene_selected <- gene[,whichElected_gene[[i]]] # 1:26, for each selected exposure
  methyl_selected <- methy[,whichElected_methyl[[i]]] # 1:26, for each selected exposure
  multiOmicsXVars[[i]]<- data.frame(metabol_serum, metabol_urine, proteome,gene_selected,methyl_selected)
  write.csv2(multiOmicsXVars[[i]], file=paste0("data_merged",allPartCorr[i],".csv"))
}

save(multiOmicsXVars,file="multiOmicsXVars.rda")

