################################################################
##### Downstream analysis on the MultiOmics MOFA model #########
################################################################
library(MOFA2)
library(ggplot2)
library(dplyr)
library(forcats)
library(lme4)

################################
#### Load MultiOmics MOFA trained model, add metadata and give an overview of the data
# In this step we load the MultiOmics MOFA trained model, shape and add metadata to the model and plot an overview of the data

model <- load_model("omics_DNAm_0.05filtered.hdf5") #load the MultiOmics MOFA trained model

load("OmicsMatrices_Mval_0.05filtered.RData") #load the tidy MultiOmics data
views_names(model) <- c("metabolomics_serum", "metabolomics_urine", "genexpr_filtered", "DNAm_denoised_filtered")
features_names(model)[[1]] <- rownames(metabolomics_serum)
features_names(model)[[2]] <- rownames(metabolomics_urine)
features_names(model)[[3]] <- rownames(genes_expr_filtered)
features_names(model)[[4]] <- rownames(DNAm_filtered_denoised)
all.equal(colnames(metabolomics_serum), colnames(metabolomics_urine), 
          colnames(genes_expr_filtered), colnames(DNAm_filtered_denoised))
samples_names(model)[[1]] <- colnames(metabolomics_serum)

# add metadata
load("../Data/exposome_NA.RData")

Nsamples <- sum(model@dimensions$N)
sample_metadata <- data.frame(
  sample = colnames(metabolomics_serum),
  covariatesNA[match(colnames(metabolomics_serum), covariatesNA$ID), -1], 
  phenotypeNA[match(colnames(metabolomics_serum), phenotypeNA$ID), -1]
)

samples_metadata(model) <- sample_metadata
head(model@samples_metadata, n=3)

plot_data_overview(model)#plot data overview


################################
#### Variance decomposition
# In this step we plot the variance explainde by Omic type and by every factor in each Omic

head(model@cache$variance_explained$r2_total[[1]])#total variance explained per Omic
plot_variance_explained(model,  y="factor", plot_total = T)[[2]]#Plot total variance explained per Omic

head(model@cache$variance_explained$r2_per_factor[[1]]) #Variance explained for every factor in each Omic
plot_variance_explained(model, x="view", y="factor")#Plot variance explained for every factor in each Omic


################################
#### Association with zBMI
# In this step we select the MultiOmics factors generated via MOFA and test their association with zBMI

factors <- as.data.frame(get_factors(model, factors = "all")[[1]])
dim(factors)
all(rownames(factors) == sample_metadata$sample)
factors$sample <- rownames(factors)
dataset <- merge(factors, sample_metadata, by="sample")

vars <- colnames(factors)[1:6] #generate a vector including the factor names
results <- matrix(NA, (ncol(factors)-1), 4, #generate a matrix to stor the results from lmm
                  dimnames=list(colnames(factors)[1:6],
                                c("nobs", "coef", "coef.se", "pval")))
for (i in (1:(ncol(factors)-1))) {#run a loop in which for each factor we fit a linear mixed model to test the association with zBMI
  var<-vars[[i]]
  covars <- c("(1|h_cohort)",var ,"e3_sex_None","e3_gac_None", "h_native_None", "h_parity_None","h_mbmi_None","h_edumc_None", "h_age_None","hs_child_age_None" )
  f <- sprintf("hs_zbmi_who ~ %s", paste0(covars, collapse=" + "))
  model_lm <- try(lmer(as.formula(f), data=(dataset))
                  , silent=TRUE)
  if (!inherits(model_lm, "try-error")) {
    coefs<-coef(summary(model_lm))
    if (var %in% rownames(coefs)) {
      results[i,"nobs"] <- as.numeric(nobs(model_lm))
      results[i,"coef"] <- coefs[var,"Estimate"]
      results[i,"coef.se"] <- coefs[var,"Std. Error"]
      results[i,"pval"] <- 
        2*pt(abs(coefs[var,"t value"]),
             df=df.residual(model_lm),lower.tail=FALSE)
               }
  }
}
results <- as.data.frame(results)


################################
#### Visualization of features weights
#In this step we visualize the weights of the features per each factor (overall and only focousing on the top based on z weight >0.9)

weights <- get_weights(model, views = "all", factors = "all", scale = T)#get the scaled weights per factor
W <- reshape2::melt(weights)#reshape the weights in order to plot them
colnames(W) <- c("feature", "factor", "value", "view")
W$feature <- as.character(W$feature)

Wscaled <- W %>% # use absolute weights and add a sign column (weight direction)
    mutate(value.abs = abs(value),
           sign = ifelse(value > 0, '+', '-'))

W.plot <- Wscaled %>%
  group_by(factor) %>% 
  top_n(5, value.abs)  # plot the top 5 for each factor
ggplot(W.plot, 
       aes(x = value.abs,
           y = fct_reorder(feature, feature),
           col = view)) +
  geom_point(size=2) +
  geom_segment(aes_string(xend="value.abs", x=0.00, 
                          yend="feature"), 
               size=0.75) +
  facet_grid(. ~ factor) +
  scale_color_manual(
    name = "View", 
    values = c("#E94F37",
               "#F2CD5D",
               "#587B7F",
               "#99d98c"),
    labels = c("DNA methylome",
               "gene expression",
               "serum metabolome",
               "urine metabolome")) +
  scale_x_continuous(
    name="relative weights", 
    limits=c(0, 1.2),
    breaks = c(0, .5, 1)) +
  geom_text(label=fct_reorder(W.plot$sign, W.plot$feature), 
            x=max(W.plot$value.abs)+0.1, size=4, col = "black") +
  theme_bw() + 
  theme(
    axis.title.x = element_text(color='black'),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=rel(1.1), hjust=1, color='black'),
    axis.text.x = element_text(color='black'),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(),
    legend.position = 'top',
    legend.title = element_blank(),
    legend.text = element_text(color="black"),
    legend.key = element_rect(fill='transparent'),
    
    # facets
    strip.background = element_rect(fill = NA),
    
    # gridlines
    panel.grid.major.y = element_blank(),
  ) 

################################################################################
##### Circular visualization of top features of  MultiOmics MOFA model #########
################################################################################

################################
#### Generate annotation data
#In this step we generate metadata to annotate to the features of the MultiOmics MOFA model
rm(list = ls())
readRDS("weights_omics_scaled_DNAm_0.05filtered.RDS")->weights#load the scaled weights
subset(weights, weights$factor=="Factor3"|weights$factor=="Factor2")->weights#subset only the factors related to BMI
library(data.table)
weights$match<-paste0(weights$view,"_", weights$factor)#tidy and select only top 5 per each omic per factor
weights <- data.table(weights, key="match")
weights.out <- weights[, .SD[value.abs %in% tail(sort(unique(value.abs)), 5)], by=match]
library(dplyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)#load the annotation for methylation and subset only features that are in our ideogram (the axis of the plot)
anno450k<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno450k<-as.data.frame(anno450k[c(anno450k$Name%in%weights.out$feature), c("Name", "UCSC_RefGene_Group","UCSC_RefGene_Name") ])#note that we only select Gene group and Gene name's columns
anno450k$UCSC_RefGene_Group<-sapply(strsplit(anno450k$UCSC_RefGene_Group, ";", fixed=TRUE), `[`, 1)
anno450k$UCSC_RefGene_Name<-sapply(strsplit(anno450k$UCSC_RefGene_Name, ";", fixed=TRUE), `[`, 1)
table(weights$view)
anno450k$chr<-"DNAm_denoised_filtered"
colnames(anno450k)<-c("Name", "band", "label", "chr")
anno450k[order(anno450k$band),]->anno450k
anno450k$coord<-c(1:nrow(anno450k))
load("../Data/genexpr.Rdata")#here we do the same also for gene expression
anno <- fData(genexpr)
anno<-as.data.frame(anno[(anno$probeset_id%in%weights.out$feature),c("probeset_id", "seqname", "GeneSymbol_Affy")])
anno$GeneSymbol_Affy<-sapply(strsplit(anno$GeneSymbol_Affy, ";", fixed=TRUE), `[`, 1)
anno$chr<-"genexpr_filtered"
colnames(anno)<-c("Name", "band", "label", "chr")
anno$band1<-gsub("chr","",anno$band)
anno[order(anno$band1),]->anno
anno$band1<-NULL
anno$coord<-c(1:nrow(anno))
load("../Data/metabol_serum.Rdata")#and also for metabolomics in serum
metab_anno <- fData(metabol_serum)
metab_anno<-as.data.frame(metab_anno[(metab_anno$Rvar.log%in%weights.out$feature),c("Rvar.log", "Class", "Rvar")])
metab_anno$chr<-"metabolomics_serum"
colnames(metab_anno)<-c("Name", "band", "label", "chr")
metab_anno[order(metab_anno$band),]->metab_anno
metab_anno$coord<-c(1:nrow(metab_anno))
metab_anno$band[metab_anno$band=="biogenicamines"]<-"BA"#we change a bit the labels to make them more readable
metab_anno$band[metab_anno$band=="glycerophospholipids"]<-"GPL"
load("../Data/metabol_urine.Rdata")#and finally we do the same also for metabolomics in urine
metab_urine <- fData(metabol_urine)
metab_urine<-as.data.frame(metab_urine[(metab_urine$Rvar%in%weights.out$feature),c("Rvar",  "var")])
metab_urine$band<-"metabolomics_urine"
metab_urine$chr<-"metabolomics_urine"
metab_urine[, c(1,3,2,4)]->metab_urine
colnames(metab_urine)<-c("Name", "band", "label", "chr")
metab_urine[order(metab_urine$band),]->metab_urine
metab_urine$coord<-c(1:nrow(metab_urine))
factors <-as.data.frame(unique(weights.out$factor))
factors$Name<-factors$`unique(weights.out$factor)`
factors$band<-factors$Name
factors$label<-factors$`unique(weights.out$factor)`
factors$chr<-"Factors"
factors<-factors[,c("Name", "band", "label", "chr")]
factors$coord<-c(1:nrow(factors))
rbind(anno, anno450k, metab_anno,metab_urine, factors)->annot
annot$label[annot$label=="3-aminoisobutyrate"]<-"BAIBA"#here we just make the labels more readable
annot$label[annot$label=="Trimethylamine oxide"]<-"TMAO"
annot$label[annot$label=="N-methyl-2-pyridone-5-carboxamide"]<-"NMPC"
annot$label[annot$label=="3-hydroxybutyrate/3-aminoisobutyrate"]<-"3-OHbutyrate"
annot$label[annot$label=="Proline betaine"]<-"Proline_betaine"
annot$label[annot$label=="4-deoxyerythronic acid"]<-"4-DEA"
annot$label[annot$label=="3-hydroxyisovalerate"]<-"3-OHvalerate"
saveRDS(annot,'Data/circos_files/ideogram_annotation0.05.rds')#save the annotation

################################
#### Circular plot of the data 
#In this step we visualize the top 5 features per each omic per each factor in a circular plot using Circos
rm(list = ls())
readRDS("Data/circos_files/ideogram_annotation0.05.rds")->annot##read the annotation
writeKaryotype <- function(annot, fileout="Data/circos_files/karyotype.txt") {#function to generate kariotype txt file
  cat("",file=fileout)
  chr <- unique(annot$chr)
  colors <- c("lblue","dyellow","vvlporange","lred", "lgreen")#define colors
  for(ii in 1:length(chr)) {
    N <- nrow(annot[annot$chr == chr[ii],])
    txt <- paste('chr -', chr[ii] , chr[ii], 0, N, colors[ii])
    cat(txt, file=fileout, append=T)
    cat("\n", file=fileout, append=T)
  }
  for(ii in 1:length(chr)) {
    band <- annot[annot$chr == chr[ii],]
    min_coord<-tapply(band$coord,band$band,min)
    max_coord<-tapply(band$coord,band$band,max)
    min_coord<-min_coord[match(unique(band$band), names(min_coord))]
    max_coord<-max_coord[match(unique(band$band), names(max_coord))]
    band <- band[!duplicated(band$band),]
    band$min_coord<-min_coord
    band$max_coord<-max_coord
    colors <- c("white","grey")#define colors
    for(jj in 1:length(band$band)) {
      if(length(band$band)==1) {
        bands <- band$band
        col <-"white"
      } else      { 
        bands <- gsub(" ", "_", band$band[jj])
        col <- ifelse((jj%%2),colors[[1]],colors[[2]])
      }
      txt <- paste('band', chr[ii], bands,bands, band$min_coord[jj]-1, band$max_coord[jj], col)
      cat(txt, file=fileout,append=T)
      cat("\n", file=fileout, append=T)
    }	
  }
}
writeKaryotype(annot)#run the function to create the Karyotype txt file
writeLabels <- function(annot, fileout="Data/circos_files/labels.txt") {#function to generate labels txt file
  cat("",file=fileout)
  for(ii in 1:nrow(annot)) {
    label <- annot$label[ii]
    label <- gsub(" ", "_", label)
    cat(sprintf('%s\n', paste(annot$chr[ii], annot$coord[ii]-1, annot$coord[ii], label)),file=fileout,append=T)
  }
}
writeLabels(annot)#run the function to create the ideogram labels
writeBandLabels <- function(annot, fileout="Data/circos_files/band_labels.txt") {#generate the function to create the band_labels txt file
  cat("",file=fileout)
  chr <- unique(annot$chr)
  for(ii in 1:length(chr)) {
    band <- annot[annot$chr == chr[ii],]
    min_coord<-tapply(band$coord, band$band, FUN=min)
    max_coord<-tapply(band$coord, band$band, FUN=max)
    min_coord<-min_coord[match(unique(band$band), names(min_coord))]
    max_coord<-max_coord[match(unique(band$band), names(max_coord))]
    band <- band[!duplicated(band$band),]
    band$min_coord<-min_coord
    band$max_coord<-max_coord
    for(jj in 1:length(band$band)) {
      bands <- gsub(" ", "_", band$band[jj])
      txt <- paste( chr[ii], band$min_coord[jj]-1, band$max_coord[jj], bands)
      cat(txt, file=fileout,append=T)
      cat("\n", file=fileout, append=T)
    }	
  }
}
writeBandLabels(annot)#run the function to generate the bands labels
readRDS("weights_omics_scaled_DNAm_0.05filtered.RDS")->weights#load the scaled weights
subset(weights, weights$factor=="Factor3"|weights$factor=="Factor2")->weights#subset only the factor related to BMI
library(data.table)
weights$match<-paste0(weights$view,"_", weights$factor)#tidy and select only top 5 per each omic per factor
weights <- data.table(weights, key="match")
weights.out <- weights[, .SD[value.abs %in% tail(sort(unique(value.abs)), 5)], by=match]
readRDS("Data/circos_files/ideogram_annotation0.05.rds")->annot#read the annotation
weights.out$Factor<-"Factors"
as.data.frame(weights.out)->dataset
dataset$value.rel<-dataset$value
writeLinks <- function(dataset, annot, fileout="Data/circos_files/links.txt") {#write the fuction to draw the link
  cat("",file=fileout)
  for(i in 1:nrow(dataset)) {
    coord1 <-  subset(annot, annot$chr == dataset[i, 'view'] & annot$Name == dataset[i, 'feature'])$coord[1]
    coord2 <- subset(annot, annot$chr == dataset[i, 'Factor'] & annot$Name == dataset[i, 'factor'])$coord[1]
    txt1 <- paste(dataset[i, 'view'], coord1-1, coord1)
    txt2 <- paste(dataset[i, 'Factor'], coord2-1, coord2)
    optionString <- sprintf(" corr=%f", dataset[i, 'value.rel'])
    cat(sprintf("%s %s%s\n", txt1, txt2, optionString),file=fileout,append=T)
  }
}
writeLinks(dataset,annot)#run the function

shell(cmd = "/circos-0.69-9/bin/circos -conf /circos-0.69-9/etc/links_omics.conf -outputfile links.png -outputdir /circos-0.69-9/etc/")#if you save all the file generate to the ect folder of circos in your computer and add the configuration file we provided then you will be able to generate the plot
