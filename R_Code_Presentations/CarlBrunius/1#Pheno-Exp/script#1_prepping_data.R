load("exposome.Rdata")

### Converting strings to integer factors and same for intervals
exposome$h_pamod_t3_None<-as.character(exposome$h_pamod_t3_None)
exposome$h_pamod_t3_None[which(as.character(exposome$h_pamod_t3_None)=="None")]<-0
exposome$h_pamod_t3_None[which(as.character(exposome$h_pamod_t3_None)=="Sometimes")]<-1
exposome$h_pamod_t3_None[which(as.character(exposome$h_pamod_t3_None)=="Often")]<-2
exposome$h_pamod_t3_None[which(as.character(exposome$h_pamod_t3_None)=="Very Often")]<-3
exposome$h_pamod_t3_None<-as.integer(exposome$h_pamod_t3_None)
exposome$h_pamod_t3_None<-as.factor(exposome$h_pamod_t3_None)

exposome$h_pavig_t3_None<-as.character(exposome$h_pavig_t3_None)
exposome$h_pavig_t3_None[which(as.character(exposome$h_pavig_t3_None)=="Low")]<-0
exposome$h_pavig_t3_None[which(as.character(exposome$h_pavig_t3_None)=="Medium")]<-1
exposome$h_pavig_t3_None[which(as.character(exposome$h_pavig_t3_None)=="High")]<-2
exposome$h_pavig_t3_None<-as.integer(exposome$h_pavig_t3_None)
exposome$h_pavig_t3_None<-as.factor(exposome$h_pavig_t3_None)

exposome$hs_pet_None<-as.character(exposome$hs_pet_None)
exposome$hs_pet_None[which(as.character(exposome$hs_pet_None)=="No")]<-0
exposome$hs_pet_None[which(as.character(exposome$hs_pet_None)=="Yes")]<-1
exposome$hs_pet_None<-as.integer(exposome$hs_pet_None)
exposome$hs_pet_None<-as.factor(exposome$hs_pet_None)

exposome$hs_tl_cdich_None<-as.character(exposome$hs_tl_cdich_None)
exposome$hs_tl_cdich_None[which(as.character(exposome$hs_tl_cdich_None)=="Undetected")]<-0
exposome$hs_tl_cdich_None[which(as.character(exposome$hs_tl_cdich_None)=="Detected")]<-1
exposome$hs_tl_cdich_None<-as.integer(exposome$hs_tl_cdich_None)
exposome$hs_tl_cdich_None<-as.factor(exposome$hs_tl_cdich_None)

exposome$hs_tl_mdich_None<-as.character(exposome$hs_tl_mdich_None)
exposome$hs_tl_mdich_None[which(as.character(exposome$hs_tl_mdich_None)=="Undetected")]<-0
exposome$hs_tl_mdich_None[which(as.character(exposome$hs_tl_mdich_None)=="Detected")]<-1
exposome$hs_tl_mdich_None<-as.integer(exposome$hs_tl_mdich_None)
exposome$hs_tl_mdich_None<-as.factor(exposome$hs_tl_mdich_None)

exposome$hs_dmdtp_cdich_None<-as.character(exposome$hs_dmdtp_cdich_None)
exposome$hs_dmdtp_cdich_None[which(as.character(exposome$hs_dmdtp_cdich_None)=="Undetected")]<-0
exposome$hs_dmdtp_cdich_None[which(as.character(exposome$hs_dmdtp_cdich_None)=="Detected")]<-1
exposome$hs_dmdtp_cdich_None<-as.integer(exposome$hs_dmdtp_cdich_None)
exposome$hs_dmdtp_cdich_None<-as.factor(exposome$hs_dmdtp_cdich_None)

exposome$FAS_cat_None<-as.character(exposome$FAS_cat_None)
exposome$FAS_cat_None[which(as.character(exposome$FAS_cat_None)=="Low")]<-0
exposome$FAS_cat_None[which(as.character(exposome$FAS_cat_None)=="Middle")]<-1
exposome$FAS_cat_None[which(as.character(exposome$FAS_cat_None)=="High")]<-2
exposome$FAS_cat_None<-as.integer(exposome$FAS_cat_None)
exposome$FAS_cat_None<-as.factor(exposome$FAS_cat_None)

exposome$hs_contactfam_3cat_num_None<-as.character(exposome$hs_contactfam_3cat_num_None)
exposome$hs_contactfam_3cat_num_None[which(as.character(exposome$hs_contactfam_3cat_num_None)=="Less than once a week")]<-0
exposome$hs_contactfam_3cat_num_None[which(as.character(exposome$hs_contactfam_3cat_num_None)=="Once a week")]<-1
exposome$hs_contactfam_3cat_num_None[which(as.character(exposome$hs_contactfam_3cat_num_None)=="(almost) Daily")]<-2
exposome$hs_contactfam_3cat_num_None<-as.integer(exposome$hs_contactfam_3cat_num_None)
exposome$hs_contactfam_3cat_num_None<-as.factor(exposome$hs_contactfam_3cat_num_None)

exposome$hs_participation_3cat_None<-as.character(exposome$hs_participation_3cat_None)
exposome$hs_participation_3cat_None[which(as.character(exposome$hs_participation_3cat_None)=="None")]<-0
exposome$hs_participation_3cat_None[which(as.character(exposome$hs_participation_3cat_None)=="1 organisation")]<-1
exposome$hs_participation_3cat_None[which(as.character(exposome$hs_participation_3cat_None)=="2 or more organisations")]<-2
exposome$hs_participation_3cat_None<-as.integer(exposome$hs_participation_3cat_None)
exposome$hs_participation_3cat_None<-as.factor(exposome$hs_participation_3cat_None)

exposome$hs_cotinine_cdich_None<-as.character(exposome$hs_cotinine_cdich_None)
exposome$hs_cotinine_cdich_None[which(as.character(exposome$hs_cotinine_cdich_None)=="Undetected")]<-0
exposome$hs_cotinine_cdich_None[which(as.character(exposome$hs_cotinine_cdich_None)=="Detected")]<-1
exposome$hs_cotinine_cdich_None<-as.integer(exposome$hs_cotinine_cdich_None)
exposome$hs_cotinine_cdich_None<-as.factor(exposome$hs_cotinine_cdich_None)

exposome$hs_cotinine_mcat_None<-as.character(exposome$hs_cotinine_mcat_None)
exposome$hs_cotinine_mcat_None[which(as.character(exposome$hs_cotinine_mcat_None)=="Non-smokers")]<-0
exposome$hs_cotinine_mcat_None[which(as.character(exposome$hs_cotinine_mcat_None)=="SHS smokers")]<-1
exposome$hs_cotinine_mcat_None[which(as.character(exposome$hs_cotinine_mcat_None)=="Smokers")]<-2
exposome$hs_cotinine_mcat_None<-as.integer(exposome$hs_cotinine_mcat_None)
exposome$hs_cotinine_mcat_None<-as.factor(exposome$hs_cotinine_mcat_None)

exposome$hs_globalexp2_None<-as.character(exposome$hs_globalexp2_None)
exposome$hs_globalexp2_None[which(as.character(exposome$hs_globalexp2_None)=="no exposure")]<-0
exposome$hs_globalexp2_None[which(as.character(exposome$hs_globalexp2_None)=="exposure")]<-1
exposome$hs_globalexp2_None<-as.integer(exposome$hs_globalexp2_None)
exposome$hs_globalexp2_None<-as.factor(exposome$hs_globalexp2_None)

exposome$hs_smk_parents_None<-as.character(exposome$hs_smk_parents_None)
exposome$hs_smk_parents_None[which(as.character(exposome$hs_smk_parents_None)=="neither")]<-0
exposome$hs_smk_parents_None[which(as.character(exposome$hs_smk_parents_None)=="one")]<-1
exposome$hs_smk_parents_None[which(as.character(exposome$hs_smk_parents_None)=="both")]<-2
exposome$hs_smk_parents_None<-as.integer(exposome$hs_smk_parents_None)
exposome$hs_smk_parents_None<-as.factor(exposome$hs_smk_parents_None)

exposome[,48]
exposome[,48] <- as.character(exposome[,48])
exposome[which(exposome[,48]=="(0,10.8]"),48] <- 0
exposome[which(exposome[,48]=="(10.8,34.9]"),48] <- 1
exposome[which(exposome[,48]=="(34.9,Inf]"),48] <- 2
exposome[,48] <- as.integer(exposome[,48])
exposome[,48] <- as.factor(exposome[,48])


exposome[,49]
exposome[,49] <- as.character(exposome[,49])
exposome[which(exposome[,49]=="(0,9]"),49] <- 0
exposome[which(exposome[,49]=="(9,27.3]"),49] <- 1
exposome[which(exposome[,49]=="(27.3,Inf]"),49] <- 2
exposome[,49] <- as.integer(exposome[,49])
exposome[,49] <- as.factor(exposome[,49])

exposome[,50]
exposome[,50] <- as.character(exposome[,50])
exposome[which(exposome[,50]=="(0,17.1]"),50] <- 0
exposome[which(exposome[,50]=="(17.1,27.1]"),50] <- 1
exposome[which(exposome[,50]=="(27.1,Inf]"),50] <- 2
exposome[,50] <- as.integer(exposome[,50])
exposome[,50] <- as.factor(exposome[,50])

exposome[,51]
exposome[,51] <- as.character(exposome[,51])
exposome[which(exposome[,51]=="(0,0.25]"),51] <- 0
exposome[which(exposome[,51]=="(0.25,0.83]"),51] <- 1
exposome[which(exposome[,51]=="(0.83,Inf]"),51] <- 2
exposome[,51] <- as.integer(exposome[,51])
exposome[,51] <- as.factor(exposome[,51])

exposome[,52]
exposome[,52] <- as.character(exposome[,52])
exposome[which(exposome[,52]=="(0,1.9]"),52] <- 0
exposome[which(exposome[,52]=="(1.9,4.1]"),52] <- 1
exposome[which(exposome[,52]=="(4.1,Inf]"),52] <- 2
exposome[,52] <- as.integer(exposome[,52])
exposome[,52] <- as.factor(exposome[,52])

exposome[,54]
exposome[,54] <- as.character(exposome[,54])
exposome[which(exposome[,54]=="(0,0.6]"),54] <- 0
exposome[which(exposome[,54]=="(0.6,18.2]"),54] <- 1
exposome[which(exposome[,54]=="(18.2,Inf]"),54] <- 2
exposome[,54] <- as.integer(exposome[,54])
exposome[,54] <- as.factor(exposome[,54])

exposome[,55]
exposome[,55] <- as.character(exposome[,55])
exposome[which(exposome[,55]=="(0,0.5]"),55] <- 0
exposome[which(exposome[,55]=="(0.5,2]"),55] <- 1
exposome[which(exposome[,55]=="(2,Inf]"),55] <- 2
exposome[,55] <- as.integer(exposome[,55])
exposome[,55] <- as.factor(exposome[,55])

levels(exposome[,56])
exposome[,56] <- as.character(exposome[,56])
exposome[which(exposome[,56]=="(0,6.5]"),56] <- 0
exposome[which(exposome[,56]=="(6.5,10]"),56] <- 1
exposome[which(exposome[,56]=="(10,Inf]"),56] <- 2
exposome[,56] <- as.integer(exposome[,56])
exposome[,56] <- as.factor(exposome[,56])

levels(exposome[,59])
exposome[,59] <- as.character(exposome[,59])
exposome[which(exposome[,59]=="(0,8.8]"),59] <- 0
exposome[which(exposome[,59]=="(8.8,16.5]"),59] <- 1
exposome[which(exposome[,59]=="(16.5,Inf]"),59] <- 2
exposome[,59] <- as.integer(exposome[,59])
exposome[,59] <- as.factor(exposome[,59])

levels(exposome[,60])
exposome[,60] <- as.character(exposome[,60])
exposome[which(exposome[,60]=="(0,2]"),60] <- 0
exposome[which(exposome[,60]=="(2,6]"),60] <- 1
exposome[which(exposome[,60]=="(6,Inf]"),60] <- 2
exposome[,60] <- as.integer(exposome[,60])
exposome[,60] <- as.factor(exposome[,60])

levels(exposome[,61])
exposome[,61] <- as.character(exposome[,61])
exposome[which(exposome[,61]=="(0,0.132]"),61] <- 0
exposome[which(exposome[,61]=="(0.132,1]"),61] <- 1
exposome[which(exposome[,61]=="(1,Inf]"),61] <- 2
exposome[,61] <- as.integer(exposome[,61])
exposome[,61] <- as.factor(exposome[,61])

levels(exposome[,62])
exposome[,62] <- as.character(exposome[,62])
exposome[which(exposome[,62]=="(0,1.1]"),62] <- 0
exposome[which(exposome[,62]=="(1.1,5.5]"),62] <- 1
exposome[which(exposome[,62]=="(5.5,Inf]"),62] <- 2
exposome[,62] <- as.integer(exposome[,62])
exposome[,62] <- as.factor(exposome[,62])

levels(exposome[,63])
exposome[,63] <- as.character(exposome[,63])
exposome[which(exposome[,63]=="(0,0.132]"),63] <- 0
exposome[which(exposome[,63]=="(0.132,Inf]"),63] <- 1
exposome[,63] <- as.integer(exposome[,63])
exposome[,63] <- as.factor(exposome[,63])

levels(exposome[,64])
exposome[,64] <- as.character(exposome[,64])
exposome[which(exposome[,64]=="(0,14.6]"),64] <- 0
exposome[which(exposome[,64]=="(14.6,25.6]"),64] <- 1
exposome[which(exposome[,64]=="(25.6,Inf]"),64] <- 2
exposome[,64] <- as.integer(exposome[,64])
exposome[,64] <- as.factor(exposome[,64])

levels(exposome[,65])
exposome[,65] <- as.character(exposome[,65])
exposome[which(exposome[,65]=="(0,0.132]"),65] <- 0
exposome[which(exposome[,65]=="(0.132,0.5]"),65] <- 1
exposome[which(exposome[,65]=="(0.5,Inf]"),65] <- 2
exposome[,65] <- as.integer(exposome[,65])
exposome[,65] <- as.factor(exposome[,65])

levels(exposome[,68])
exposome[,68] <- as.character(exposome[,68])
exposome[which(exposome[,68]=="(0,0.132]"),68] <- 0
exposome[which(exposome[,68]=="(0.132,1]"),68] <- 1
exposome[which(exposome[,68]=="(1,Inf]"),68] <- 2
exposome[,68] <- as.integer(exposome[,68])
exposome[,68] <- as.factor(exposome[,68])

levels(exposome[,72])
exposome[,72] <- as.character(exposome[,72])
exposome[which(exposome[,72]=="(0,1.5]"),72] <- 0
exposome[which(exposome[,72]=="(1.5,4]"),72] <- 1
exposome[which(exposome[,72]=="(4,Inf]"),72] <- 2
exposome[,72] <- as.integer(exposome[,72])
exposome[,72] <- as.factor(exposome[,72])

levels(exposome[,73])
exposome[,73] <- as.character(exposome[,73])
exposome[which(exposome[,73]=="(0,0.132]"),73] <- 0
exposome[which(exposome[,73]=="(0.132,0.5]"),73] <- 1
exposome[which(exposome[,73]=="(0.5,Inf]"),73] <- 2
exposome[,73] <- as.integer(exposome[,73])
exposome[,73] <- as.factor(exposome[,73])

levels(exposome[,75])
exposome[,75] <- as.character(exposome[,75])
exposome[which(exposome[,75]=="(0,7]"),75] <- 0
exposome[which(exposome[,75]=="(7,17.5]"),75] <- 1
exposome[which(exposome[,75]=="(17.5,Inf]"),75] <- 2
exposome[,75] <- as.integer(exposome[,75])
exposome[,75] <- as.factor(exposome[,75])

levels(exposome[,76])
exposome[,76] <- as.character(exposome[,76])
exposome[which(exposome[,76]=="(0,14.1]"),76] <- 0
exposome[which(exposome[,76]=="(14.1,23.6]"),76] <- 1
exposome[which(exposome[,76]=="(23.6,Inf]"),76] <- 2
exposome[,76] <- as.integer(exposome[,76])
exposome[,76] <- as.factor(exposome[,76])

levels(exposome[,77])
exposome[,77] <- as.character(exposome[,77])
exposome[which(exposome[,77]=="(0,1.5]"),77] <- 0
exposome[which(exposome[,77]=="(1.5,3]"),77] <- 1
exposome[which(exposome[,77]=="(3,Inf]"),77] <- 2
exposome[,77] <- as.integer(exposome[,77])
exposome[,77] <- as.factor(exposome[,77])

levels(exposome[,78])
exposome[,78] <- as.character(exposome[,78])
exposome[which(exposome[,78]=="(0,7]"),78] <- 0
exposome[which(exposome[,78]=="(7,14.1]"),78] <- 1
exposome[which(exposome[,78]=="(14.1,Inf]"),78] <- 2
exposome[,78] <- as.integer(exposome[,78])
exposome[,78] <- as.factor(exposome[,78])

levels(exposome[,79])
exposome[,79] <- as.character(exposome[,79])
exposome[which(exposome[,79]=="(0,3]"),79] <- 0
exposome[which(exposome[,79]=="(3,7]"),79] <- 1
exposome[which(exposome[,79]=="(7,Inf]"),79] <- 2
exposome[,79] <- as.integer(exposome[,79])
exposome[,79] <- as.factor(exposome[,79])

levels(exposome[,80])
exposome[,80] <- as.character(exposome[,80])
exposome[which(exposome[,80]=="(0,6]"),80] <- 0
exposome[which(exposome[,80]=="(6,9]"),80] <- 1
exposome[which(exposome[,80]=="(9,Inf]"),80] <- 2
exposome[,80] <- as.integer(exposome[,80])
exposome[,80] <- as.factor(exposome[,80])

levels(exposome[,81])
exposome[,81] <- as.character(exposome[,81])
exposome[which(exposome[,81]=="(0,3]"),81] <- 0
exposome[which(exposome[,81]=="(3,4]"),81] <- 1
exposome[which(exposome[,81]=="(4,Inf]"),81] <- 2
exposome[,81] <- as.integer(exposome[,81])
exposome[,81] <- as.factor(exposome[,81])

levels(exposome[,82])
exposome[,82] <- as.character(exposome[,82])
exposome[which(exposome[,82]=="(0,4.1]"),82] <- 0
exposome[which(exposome[,82]=="(4.1,8.5]"),82] <- 1
exposome[which(exposome[,82]=="(8.5,Inf]"),82] <- 2
exposome[,82] <- as.integer(exposome[,82])
exposome[,82] <- as.factor(exposome[,82])

levels(exposome[,83])
exposome[,83] <- as.character(exposome[,83])
exposome[which(exposome[,83]=="(0,6]"),83] <- 0
exposome[which(exposome[,83]=="(6,8.5]"),83] <- 1
exposome[which(exposome[,83]=="(8.5,Inf]"),83] <- 2
exposome[,83] <- as.integer(exposome[,83])
exposome[,83] <- as.factor(exposome[,83])

levels(exposome[,84])
exposome[,84] <- as.character(exposome[,84])
exposome[which(exposome[,84]=="(0,6]"),84] <- 0
exposome[which(exposome[,84]=="(6,8.5]"),84] <- 1
exposome[which(exposome[,84]=="(8.5,Inf]"),84] <- 2
exposome[,84] <- as.integer(exposome[,84])
exposome[,84] <- as.factor(exposome[,84])

save.image("exposomeEdited.Rdata")