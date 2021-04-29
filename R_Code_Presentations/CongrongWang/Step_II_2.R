########## Downstream Analysis ##########
########## Exposome MOFA model ##########

library(ggplot2)
library(MOFA2)
library(dplyr)
library(tidyverse)
library(psych)
library(lme4)
library(lmerTest)


# ----- Load the MOFA model -----

model <- load_model("exposome_4mat.hdf5")


# ----- Add meta data -----

Nsamples <- sum(model@dimensions$N)

omics_factors <- readRDS("factors_DNAm_0.05filtered.RDS")
colnames(omics_factors)[2:7] <- paste0("Omics.Factor", 1:6)

sample_metadata <- data.frame(
  sample = samples_names(model)[[1]],
  covariatesNA[match(samples_names(model)[[1]], covariatesNA$ID), -1], 
  phenotypeNA[match(samples_names(model)[[1]], phenotypeNA$ID), -1],
  omics_factors[match(samples_names(model)[[1]], omics_factors$sample), 2:4] # use only omics factor1:3
)

samples_metadata(model) <- sample_metadata
head(model@samples_metadata, n=3)


########## Downstream analysis ##########

# ----- Plot the data -----

plot_data_overview(model)

# ----- Variance decomposition -----

# Total variance explained per view
model@cache$variance_explained$r2_total$single_group
plot_variance_explained(model, x="view", y="factor", plot_total = T)[[2]]

# Variance explained for every factor in per view
model@cache$variance_explained$r2_per_factor$single_group 
plot_variance_explained(model, x="view", y="factor") 


# ----- Extract factors -----

factors <- get_factors(model, factors = "all", groups = 1)[[1]]
dim(factors)
colnames(factors) <- paste0("Expo.Factor", 1:3)


# ----- Regression of BMI/reduced multi-omics on the 3 factors -----

data <- data.frame(sample_metadata, 
                   factors[match(sample_metadata$sample, rownames(factors)), 1:3])

# with BMI?
results_1 <- matrix(NA, (ncol(factors)), 4,
                  dimnames=list(colnames(factors),
                                c("nobs", "coef", "coef.se", "pval")))
for (i in (1:(ncol(factors)))) {
  factor <- colnames(factors)[i]
  covars <- c("(1|h_cohort)", factor,"e3_sex_None","e3_gac_None", "h_native_None", "h_parity_None","h_mbmi_None","h_edumc_None", "h_age_None","hs_child_age_None" )
  f <- sprintf("hs_zbmi_who ~ %s", paste0(covars, collapse=" + "))
  
  model_lm <- lmer(as.formula(f), data=data)#, silent=TRUE)
  if (!inherits(model_lm, "try-error")) {
    coefs <- coef(summary(model_lm))
    if (factor %in% rownames(coefs)) {
      results_1[i,"nobs"] <- as.numeric(nobs(model_lm))
      results_1[i,"coef"] <- coefs[factor, "Estimate"]
      results_1[i,"coef.se"] <- coefs[factor, "Std. Error"]
      results_1[i,"pval"] <- coefs[factor, "Pr(>|t|)"]
    }
  }
}
write.csv(results_1, "results_expo_bmi.csv")

# with Omics factors 1:3?
results_2 <- data.frame()
for (i in (1:(ncol(factors)))) {
  for (j in 1:3) {
    omic.factor <- paste0("Omics.Factor", j)
    factor <- colnames(factors)[i]
    
    covars <- c("(1|h_cohort)", factor,"e3_sex_None","e3_gac_None", "h_native_None", "h_parity_None","h_mbmi_None","h_edumc_None", "h_age_None","hs_child_age_None" )
    f <- sprintf("%s ~ %s", omic.factor, paste0(covars, collapse=" + "))
  
    model_lm <- lmer(as.formula(f), data=data)
    if (!inherits(model_lm, "try-error")) {
      coefs <- coef(summary(model_lm))
      if (factor %in% rownames(coefs)) {
        results_2 <- rbind(
          results_2,
          data.frame(
            Expo.factor = colnames(factors)[i], 
            Omics.factor = omic.factor,
            nobs = as.numeric(nobs(model_lm)), 
            coef = coefs[factor, "Estimate"],
            coef.se = coefs[factor, "Std. Error"],
            pval = coefs[factor, "Pr(>|t|)"]
            )
          )
      }
    }
  }
}


# ----- Weights and features -----

weights <- get_weights(model, views = "all", factors = "all", scale = T) # Scaled weights
W <- reshape2::melt(weights)
colnames(W) <- c("feature", "factor", "value", "view")

# Correct the feature names in data:
W$feature <- as.character(W$feature)
for (i in 1:nrow(W)) {
  W$feature[i] <- gsub(x = W$feature[i], pattern = W$view[i],
                       replacement = "", fixed = TRUE)
}
for (f in 1:ncol(factors)) {
  
  factor.name <- paste0("Factor", f)
  
  for (i in 1:length(matrix.list)) {
    
    view.name <- names(matrix.list)[i]
    
    for (j in 1:nrow(matrix.list[[i]])) {
      
      pat <- rownames(matrix.list[[i]])[j]
      x <- W$feature[W$view==view.name & W$factor==factor.name][j]
      ind <- grep(pat, x, fixed = T)
      print(ind)
      
      W$feature[W$view==view.name & W$factor==factor.name][j] <- 
        ifelse(ind, pat, NA)
    }
  }
}

Wscaled <- W %>%              # use absolute weights and add a sign column (weight direction)
    mutate(value.abs = abs(value),
           sign = ifelse(value > 0, '+', '-'))

# Plot the top 5 features for each factor
W.plot <- Wscaled %>%
  group_by(factor) %>% 
  top_n(5, value.abs)  
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
    labels = c("postnatal binary exposures",
               "postnatal numerical exposures",
               "prenatal binary exposures",
               "prenatal numerical exposures")) +
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
##### Circular visualization of top features of  Exposome MOFA model #########
################################################################################

################################
#### Generate annotation data
#In this step we generate metadata to annotate to the features of the Exposome MOFA model
rm(list = ls())
readRDS("./weights_exposome_scaled.rds")->weights#read the scaled weights
library(data.table)
weights<-subset(weights,weights$factor%in%c("Factor1", "Factor2"))#subset only the factors related to the reduced multiomics 
weights$match<-paste0(weights$view,"_", weights$factor)#tidy and select only top 5 per each omic per factor
weights <- data.table(weights, key="match")
weights.out <- weights[, .SD[value.abs %in% tail(sort(unique(value.abs)), 5)], by=match]
library(dplyr)
gsub("exposome_", "",weights.out$view)->weights.out$chr
weights.out$Name<-weights.out$feature
weights.out$band<-weights.out$view
weights.out$label<-weights.out$feature
exp<-weights.out[,c("Name", "band", "label", "chr")]
readxl::read_xlsx("../Data/codebook.xlsx", sheet = 1)->codebook#load the annotation for methylation and subset only features that are in our ideogram
codebook[,c("variable_name", "domain", "family")]->codebook
gsub('[(].*',"",exp$label)->exp$variable_name
exp$variable_name <- sub("\\_$", "", exp$variable_name)
merge(exp, codebook, by="variable_name", all.x=T)->exp
exp$band<-exp$family
gsub("h_","",exp$label)->exp$label
gsub("hs_","",exp$label)->exp$label
gsub("_.*","",exp$label)->exp$label
exp<-exp[,c("Name", "band", "label", "chr")]
exp$match<-paste0(exp$Name, exp$chr)
exp[!duplicated(exp$match),]->exp
exp$match<-NULL
as.character(exp$label)->exp$label
as.character(exp$band)->exp$band
exp$label[duplicated(exp$label)]
exp$label[exp$Name=="hs_total_bread_Ter_(17.5,Inf]"]<-"bread"#we change a bit the labels to make them more readable
exp$label[exp$Name=="hs_total_fish_Ter_(3,Inf]"]<-"fish"
exp$label[exp$Name=="h_accesslines300_preg_dic0"]<-"Haccesslines300"
exp$label[exp$Name=="hs_accesslines300_s_dic0_1"]<-"Saccesslines300"
exp$label[exp$Name=="hs_sd_wk_None"]<-"physActiv"
exp$label[exp$Name=="h_cereal_preg_Ter_(9,27.3]"]<-"cereal2"
exp$label[exp$Name=="h_cereal_preg_Ter_(27.3,Inf]"]<-"cereal3"
exp$label[exp$Name=="h_legume_preg_Ter_(0.5,2]"]<-"legume2"
exp$label[exp$Name=="h_legume_preg_Ter_(2,Inf]"]<-"legume3"
exp$label[exp$Name=="hs_builtdens300_h_Sqrt"]<-"builtdens30H"
exp$label[exp$Name=="hs_builtdens300_s_Sqrt"]<-"builtdens30S"
exp$label[exp$Name=="hs_popdens_h_Sqrt"]<-"popdensH"
exp$label[exp$Name=="hs_popdens_s_Sqrt"]<-"popdensS"
exp$band[is.na(exp$band)]<-c("Social","Lifestyle", "Builtenv","Builtenv","Smoke","Nature")
exp$band[exp$band=="Air Pollution"]<-"Pollution"
exp$band[exp$band=="Built environment"]<-"Builtenv"
exp$band[exp$band=="Tobacco Smoke"]<-"Smoke"
exp$band[exp$band=="Meteorological"]<-"Meteo"
exp[order(exp$band),]->exp
exp_pre_bin<-subset(exp,exp$chr=="pre_bin")
exp_post_bin<-subset(exp,exp$chr=="post_bin")
exp_pre_num<-subset(exp,exp$chr=="pre_num")
exp_post_num<-subset(exp,exp$chr=="post_num")
exp_pre_bin$coord<-c(1:nrow(exp_pre_bin))
exp_post_bin$coord<-c(1:nrow(exp_post_bin))
exp_pre_num$coord<-c(1:nrow(exp_pre_num))
exp_post_num$coord<-c(1:nrow(exp_post_num))
factors <-as.data.frame(unique(weights.out$factor))
factors$Name<-factors$`unique(weights.out$factor)`
factors$band<-factors$Name
factors$label<-factors$`unique(weights.out$factor)`
factors$chr<-"Factors"
factors<-factors[,c("Name", "band", "label", "chr")]
factors$coord<-c(1:nrow(factors))
rbind(exp_pre_bin,exp_pre_num,exp_post_bin, exp_post_num, factors)->annotation
as.data.frame(annotation)->annot
saveRDS(annot,'Data/circos_files/exposome/ideogram_annotation_exp005.rds')

################################
#### Circular plot of the data 
#In this step we visualize the top 5 features per each exposome view per each factor in a circular plot using Circos
rm(list = ls())
readRDS("Data/circos_files/exposome/ideogram_annotation_exp005.rds")->annot##read the annotation
writeKaryotype <- function(annot, fileout="Data/circos_files/exposome/karyotype.txt") {##generate the function to create the Karyotype txt file
  cat("",file=fileout)
  chr <- unique(annot$chr)
  colors <- c("lblue","dyellow","vvlporange","lred", "lgreen")#define the colors
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
    colors <- c("white","grey")#define the colors
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
writeKaryotype(annot)##run the function to create the Karyotype txt file
writeLabels <- function(annot, fileout="Data/circos_files/exposome/labels.txt") {##generate the function to create the ideogram labels
  cat("",file=fileout)
  for(ii in 1:nrow(annot)) {
    label <- annot$label[ii]
    label <- gsub(" ", "_", label)
    cat(sprintf('%s\n', paste(annot$chr[ii], annot$coord[ii]-1, annot$coord[ii], label)),file=fileout,append=T)
  }
}
writeLabels(annot)##run the function to create the ideogram labels
writeBandLabels <- function(annot, fileout="Data/circos_files/exposome/band_labels.txt") {##generate the function to create the band_labels txt file
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
readRDS("./weights_exposome_scaled.rds")->weights#read the scaled weights
library(data.table)
weights<-subset(weights,weights$factor%in%c("Factor1", "Factor2"))#subset only the factors related to the reduced multiomics 
weights$match<-paste0(weights$view,"_", weights$factor)#tidy and select only top 5 per each omic per factor
weights <- data.table(weights, key="match")
weights.out <- weights[, .SD[value.abs %in% tail(sort(unique(value.abs)), 5)], by=match]
gsub("exposome_", "",weights.out$view)->weights.out$chr
readRDS("Data/circos_files/exposome/ideogram_annotation_exp005.rds")->annot##read the annotation
weights.out$Factor<-"Factors"
as.data.frame(weights.out)->dataset
dataset$value.rel<-dataset$value
annot$Name<-as.character(annot$Name)
writeLinks <- function(dataset, annot, fileout="Data/circos_files/exposome/links.txt") {#write the fuction to draw the link
  cat("",file=fileout)
  for(i in 1:nrow(dataset)) {
    coord1 <-  subset(annot, annot$chr == dataset[i, 'chr'] & annot$Name == dataset[i, 'feature'])$coord[1]
    coord2 <- subset(annot, annot$chr == dataset[i, 'Factor'] & annot$Name == dataset[i, 'factor'])$coord[1]
    txt1 <- paste(dataset[i, 'chr'], coord1-1, coord1)
    txt2 <- paste(dataset[i, 'Factor'], coord2-1, coord2)
    optionString <- sprintf(" corr=%f", dataset[i, 'value.rel'])
    cat(sprintf("%s %s%s\n", txt1, txt2, optionString),file=fileout,append=T)
  }
}
writeLinks(dataset,annot)#run the function

shell(cmd = "/circos-0.69-9/bin/circos -conf /circos-0.69-9/etc/circos_files_omics/exp/links_exp.conf -outputfile links.png -outputdir C:/circos-0.69-9/etc/")#if you save all the file generate to the ect folder of circos in your computer and add the configuration file we provided then you will be able to generate the plot
