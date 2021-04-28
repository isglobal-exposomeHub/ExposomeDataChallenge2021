setwd("/Users/qwu/Dropbox/AwardCompetitionFellow/isglobal_2021/Submit_files")
library(R.matlab)
library(pracma)

load("metabol_serum.Rdata")
metabol_serum$ID

load("metabol_urine.Rdata")
metabol_urine$ID

overlapID = intersect(metabol_serum$ID,metabol_urine$ID)
overlapIDX = match(overlapID,metabol_serum$ID)
size(overlapIDX)
serum = exprs(metabol_serum)
serum = serum[,overlapIDX]
urine = exprs(metabol_urine)
meta = rbind(serum,urine)

serum2 = as(metabol_serum,"data.frame")
urine2 = as(metabol_urine,"data.frame")
meta_names = c(names(serum2)[1:177],names(urine2)[1:44])
writeMat('meta_merge.mat',meta=meta,meta_names=meta_names)

load("exposome_NA.Rdata")
num_cols = unlist(lapply(exposomeNA, is.numeric))         # Identify numeric columns
num_cols
exposomeNA_num = exposomeNA[,num_cols]
overlapID2 = intersect(metabol_urine$ID,exposomeNA_num$ID)
overlapIDX2 = match(overlapID2,exposomeNA_num$ID)
size(overlapIDX2)
exposomeNA_num = data.matrix(exposomeNA_num, rownames.force = NA)
exposomeNA_num2 = t(exposomeNA_num[overlapIDX2,])
exposomeNA_num2 = exposomeNA_num2[-1,]

expos_names = names(exposomeNA)[num_cols]
expos_names = expos_names[-1]

writeMat('exposomeNA_num.mat',expos=exposomeNA_num[,-1])
writeMat('expos_merge.mat',expos=exposomeNA_num2,expos_names=expos_names)
