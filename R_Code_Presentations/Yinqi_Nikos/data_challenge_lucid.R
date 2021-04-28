#### 1. load data and packages ####
# This script is constructed based on an updated version of LUCID
# This version corresponds to the stable-v1 branch on Github
# https://github.com/USCbiostats/LUCIDus/tree/stable-v1
# install.packages("LUCIDus_2.1.1.tar.gz", repos = NULL, type = "source")
library(Biobase)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mclust)
library(reshape2)
library(LUCIDus)
library(fmsb) # for radar plot
library(R.utils) # keep track the time
library(networkD3)

load("~/Documents/Research/exposome_data_challenge/data/exposome.rdata") # exposure without missingness
load("~/Documents/Research/exposome_data_challenge/data/metabol_serum.rdata")
load("~/Documents/Research/exposome_data_challenge/data/metabol_urine.rdata")
load("~/Documents/Research/exposome_data_challenge/data/proteome.rdata")

#### 2. data clean ####
# 1. individual level data
expo_names = as.character(codebook$variable_name[codebook$family == "Organochlorines"])
cova_names = c("h_mbmi_None", "e3_sex_None", "h_age_None", "h_cohort", "h_edumc_None")
expo = as_tibble(exposome[, c("ID", expo_names)])
cova = as_tibble(covariates[, c("ID", cova_names)])
phen = as_tibble(phenotype[, c("ID", "hs_bmi_c_cat", "hs_zbmi_who")])
# 2. serum data (1198 x 177 features)
expr_serum_raw = exprs(metabol_serum)
expr_serum_sample_id = colnames(expr_serum_raw)
expr_serum_feature_id = paste0("serum_", rownames(expr_serum_raw))
expr_serum = as_tibble(t(expr_serum_raw))
colnames(expr_serum) = expr_serum_feature_id
expr_serum$ID = as.integer(expr_serum_sample_id)
expr_serum = expr_serum[, c(ncol(expr_serum), 1:(ncol(expr_serum) - 1))]
# 3. urine data (1192 x 44 features): 
expr_urine_raw = exprs(metabol_urine)
expr_urine_sample_id = colnames(expr_urine_raw)
expr_urine_feature_id = paste0("urine_", rownames(expr_urine_raw))
expr_urine = as_tibble(t(expr_urine_raw)) 
colnames(expr_urine) = expr_urine_feature_id
expr_urine$ID = as.integer(expr_urine_sample_id)
expr_urine = expr_urine[, c(ncol(expr_urine), 1:(ncol(expr_urine) - 1))]
# 4. proteome data (1170 x 36 features)
expr_proteome_raw = exprs(proteome)
expr_proteome_sample_id = colnames(expr_proteome_raw)
expr_proteome_feature_id = rownames(expr_proteome_raw)
expr_proteome_annot = fData(proteome)
expr_proteome = as_tibble(t(expr_proteome_raw))
colnames(expr_proteome) = expr_proteome_feature_id
expr_proteome$ID = as.integer(expr_proteome_sample_id)
expr_proteome = expr_proteome[, c(ncol(expr_proteome), 1:(ncol(expr_proteome) - 1))]

full_dat = inner_join(expo, cova, by = "ID") %>%
  inner_join(., phen, by = "ID") %>% 
  inner_join(., expr_serum, by = "ID") %>% 
  inner_join(., expr_urine, by = "ID") %>% 
  inner_join(., expr_proteome, by = "ID")
dim(full_dat) # 1152 obs x 283 variables (ID, 36 + 44 + 177 = 257 metabolites, 2 pheno, 5 cova, 18 exposures)
sum(is.na(full_dat)) # 0 NA value

write.csv(full_dat, file = "data/dat_organochlorines.csv", quote = FALSE, row.names = FALSE)

# This part of the code is extracted from
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
cor_dat = as.data.frame(full_dat)
for (i in 1:ncol(cor_dat)) {
  if(is.factor(cor_dat[, i])) {
    cor_dat[, i] = as.numeric(cor_dat[, i])
  }
}
cormat = round(cor(cor_dat[, -1]), 2)
dim(cormat)

# Get lower triangle of the correlation matrix
get_lower_tri = function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
cormat_lower = get_lower_tri(cormat = cormat)
cormat_melt = melt(cormat, na.rm = TRUE)
head(cormat_melt)

ggplot(data = cormat_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(text = element_text(size=3),
        axis.text.x = element_text(angle = 90, vjust = 1, 
                                   hjust = 1))+
  coord_fixed()



#### For continuous outcome ####
# all the associations are adjusted for covariates
# 1.1 E and Y
expo_names = colnames(expo)[-1]
res_E_Y = NULL
for (i in 1:length(expo_names)) {
  temp_fit = lm(as.formula(paste("hs_zbmi_who ~", expo_names[i], " + ", paste(cova_names, collapse = " + "))), 
                data = full_dat)
  temp_sum = summary(temp_fit)
  res_E_Y = rbind(res_E_Y, temp_sum$coefficients[2, ])
}
rownames(res_E_Y) = expo_names
write.csv(round(res_E_Y, 3), "result/association_cont_Y_E.csv", quote = FALSE)

# 1.2 M and Y
M_names = c(colnames(expr_serum)[-1], colnames(expr_urine)[-1], colnames(expr_proteome)[-1])
res_M_Y = NULL
for (i in 1:length(M_names)) {
  temp_fit = lm(as.formula(paste("hs_zbmi_who ~", M_names[i], " + ", paste(cova_names, collapse = " + "))), 
                data = full_dat)
  temp_sum = summary(temp_fit)
  res_M_Y = rbind(res_M_Y, temp_sum$coefficients[2, ])
}
rownames(res_M_Y) = M_names
p_adj1 = p.adjust(res_M_Y[, 4], method = "fdr")
res_M_Y = cbind(res_M_Y, p_adj1)
colnames(res_M_Y)[5] = "p_fdr"
write.csv(round(res_M_Y, 3), "result/association_cont_Y_M.csv", quote = FALSE)
sum(res_M_Y[, 5] < 0.1)
sum(res_M_Y[, 5] < 0.05)

# 1.3 E and M
res_M_E = NULL
for (i in 1:length(M_names)) {
  for (j in 1:length(expo_names)) {
    temp_fit = lm(as.formula(paste(M_names[i], "~", expo_names[j], "+", paste(cova_names, collapse = " +"))),
                  data = full_dat)
    temp_sum = data.frame(est = summary(temp_fit)$coefficients[2, 1],
                          p = summary(temp_fit)$coefficients[2, 4])
    temp_sum$M = M_names[i]
    temp_sum$E = expo_names[j]
    res_M_E = rbind(res_M_E, temp_sum)
  }
}
res_M_E$p.adj = p.adjust(res_M_E$p, method = "fdr")
write.csv(res_M_E, file = "result/asociation_cont_M_E.csv", quote = FALSE, row.names = FALSE)





#### 2. preliminary screening ####
# 1. M should be significantly associated with at least 1 exposure
# 2. M should be significantly associated with the outcome
# cutoff: fdr = 0.1
name1 = unique(rownames(res_M_Y)[res_M_Y[, 5] < 0.1])
name2 = unique(res_M_E$M[res_M_E$p.adj < 0.1])
common_name = intersect(name1, name2) 
length(common_name) # 45 metabolites passed preliminary screening





#### 3. LUCID modeling ####
exposure = as.matrix(full_dat[, colnames(expo)[-1]])
metabolite = scale(as.matrix(full_dat[, common_name])) # scale the metabolite
outcome = as.matrix(full_dat[, "hs_zbmi_who"])
covariate = model.matrix(~., full_dat[, colnames(cova)[-1]])[, -1]

# correlation heatmap for exposure, 45 metabolites and the outcome
cor_dat2 = as.data.frame(cbind(exposure, metabolite))
for (i in 1:ncol(cor_dat2)) {
  if(is.factor(cor_dat2[, i])) {
    cor_dat2[, i] = as.numeric(cor_dat2[, i])
  }
}
cormat2 = round(cor(cor_dat2[, -1]), 2)
cormat_lower2 = get_lower_tri(cormat = cormat2)
cormat_melt2 = melt(cormat2, na.rm = TRUE)
head(cormat_melt2)

ggplot(data = cormat_melt2, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1))+
  coord_fixed() +
  xlab("") + 
  ylab("")
ggsave("plot/heatmapt_45_metabolites.pdf", height = 8, width = 8)

# test lucid model
set.seed(123)
fit_try1 = est.lucid(G = exposure, Z = metabolite, Y = outcome, 
                     useY = FALSE, K = 2, family = "normal")
summary(fit_try1)
hist(fit_try1$post.p[, 1])


# tune the model
set.seed(123)
res_K = 2:4
res_Rho_Z_InvCov = c(0.1, 0.2, 0.3)
res_Rho_Z_CovMu = c(30, 40, 50)
res_Rho_G = c(0.005, 0.01)
res_BIC = expand.grid(K = res_K, Z_InvCov = res_Rho_Z_InvCov, Z_Mu = res_Rho_Z_CovMu, E = res_Rho_G)
res_BIC = as.data.frame(res_BIC)
res_BIC$BIC = rep(0, nrow(res_BIC))
res_BIC$likelihood = rep(0, nrow(res_BIC))
for (i in 1:nrow(res_BIC)) {
  temp_fit = try(est.lucid(G = exposure, Z = metabolite, Y = outcome, CoY = covariate,
                           useY = FALSE, K = res_BIC$K[i], family = "normal", 
                           tune = def.tune(Select_Z = TRUE, Rho_Z_InvCov = res_BIC$Z_InvCov[i], Rho_Z_CovMu = res_BIC$Z_Mu[i],
                                           Select_G = TRUE, Rho_G = res_BIC$E[i])))
  if("try-error" %in% class(temp_fit)) {
    next
  } else {
    res_BIC$BIC[i] = summary(temp_fit)$BIC
    res_BIC$likelihood[i] = temp_fit$likelihood
  }
}
write.csv(res_BIC, file = "result/tune1.csv", quote = FALSE, row.names = FALSE)
res_BIC = res_BIC[res_BIC$BIC != 0, ]
res_BIC[res_BIC$BIC == min(res_BIC$BIC), ] # 4, 0.1, 50, 3 is the best
res_BIC = as_tibble(res_BIC)
res_BIC$Z_InvCov = as.factor(res_BIC$Z_InvCov)
res_BIC$Z_Mu = as.factor(res_BIC$Z_Mu)

# plot for tuning process
ggplot(data = res_BIC, aes(x = K, y = BIC, group = interaction(Z_InvCov, Z_Mu, E))) +
  geom_line(aes(color = Z_InvCov, linetype = as.factor(E))) +
  geom_point(aes(shape = Z_Mu)) +
  scale_x_continuous(name = "K (number of latent cluster)",
                     breaks = c(2, 3, 4)) +
  labs(color = "GLASSO penalty(Z)", shape = "LASSO penalty(Z)", linetype = "LASSO penalty(E)")
ggsave("plot/tune_contY.pdf", width = 6, height = 4)


# refit the model with the optimal combination
set.seed(123)
fit_try2 = est.lucid(G = exposure, Z = metabolite, Y = outcome, CoY = covariate,
                     useY = FALSE, K = 4, family = "normal", 
                     tune = def.tune(Select_Z = TRUE, Rho_Z_InvCov = 0.1, Rho_Z_CovMu = 50,
                                     Select_G = TRUE, Rho_G = 0.01))
summary(fit_try2)
hist(fit_try2$post.p[, 1])
hist(fit_try2$post.p[, 2])
hist(fit_try2$post.p[, 3])
hist(fit_try2$post.p[, 4])

# refit the model with the selected features
set.seed(123)
fit_try3 = est.lucid(G = exposure[, fit_try2$select$selectG], 
                     Z = metabolite[, fit_try2$select$selectZ], 
                     Y = outcome, CoY = covariate,
                     useY = FALSE, K = 4, family = "normal")
summary(fit_try3)
pred_fit_try3 = predict(fit_try3, newG = exposure[, fit_try2$select$selectG],
                        newZ = metabolite[, fit_try2$select$selectZ], 
                        CoY = covariate)
table(pred_fit_try3$pred.x)
# 1   2   3   4 
# 593  37 386 136 



# plot the sankey diagram
plot(fit_try3)
# personalize color by metabolites group
x = fit_try3
K <- x$K
var.names <- x$var.names
var.names$Ynames = "BMI z-score"
pars <- x$pars
dimG <- length(var.names$Gnames)
dimZ <- length(var.names$Znames)
valueGtoX <- as.vector(t(x$pars$beta[, -1]))
valueXtoZ <- as.vector(t(x$pars$mu))
valueXtoY <- as.vector(x$pars$gamma$beta)[1:K]
GtoX <- data.frame(source = rep(x$var.names$Gnames, K),
                   target = paste0("Latent Cluster", as.vector(sapply(1:K, function(x) rep(x, dimG)))),
                   value = abs(valueGtoX),
                   group = as.factor(valueGtoX > 0))
XtoZ <- data.frame(source = paste0("Latent Cluster", as.vector(sapply(1:K, function(x) rep(x, dimZ)))),
                   target = rep(var.names$Znames, K),
                   value = abs(valueXtoZ),
                   group = as.factor(valueXtoZ > 0))
XtoY <- data.frame(source = paste0("Latent Cluster", 1:K),
                   target = rep(var.names$Ynames, K),
                   value = abs(valueXtoY),
                   group = as.factor(valueXtoY > 0))
if(x$family == "binary"){
  XtoY$value <- exp(valueXtoY)
}
links <- rbind(GtoX, XtoZ, XtoY)
nodes <- data.frame(name = unique(c(as.character(links$source), as.character(links$target))),
                    group = as.factor(c(rep("exposure", dimG), 
                                        rep("lc", K), 
                                        rep("serum", 6),
                                        "urine",
                                        rep("protein", 7),
                                        "outcome")))
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
links$group2 = rep("a", nrow(links))
links$group2[links$group == FALSE] = "negative"
links$group2[links$group == TRUE] = "positive"
my_color <- 'd3.scaleOrdinal() .domain(["exposure", "lc", "serum", "urine", "protein", "outcome", "positive", "negative"]) .range(["dimgray", "#eb8c30", "red", "#2fa4da", "#2ECC71", "#afa58e", "#67928b", "#d1e5eb"])'
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name",
                   colourScale = my_color, LinkGroup ="group2", NodeGroup ="group",
                   sinksRight = FALSE, fontSize = 12)
p


# plot for posterior distribution of outcome
Y_fit_try3 = as.data.frame(cbind(cluster = as.factor(pred_fit_try3$pred.x), full_dat[, "hs_zbmi_who"]))
Y_fit_try3 = melt(Y_fit_try3, id.vars = "cluster")
ggplot(Y_fit_try3, aes(factor(variable), value, goup = factor(cluster))) + 
  geom_boxplot() + 
  xlab("")
ggsave("plot/posterior_boxplot_zbmi.pdf", width = 10, height = 7)



#### 4. omics profiles ####
# use stacked barplot for mean estimate
M_mean = as.data.frame(fit_try3$pars$mu)
M_mean$cluster = as.factor(1:4)

# layout 4 by 1
ggplot(M_mean_melt, aes(fill = variable, y = value, x = variable)) + 
  geom_bar(position="dodge", stat="identity") + 
  ggtitle("Omics profiles for 4 latent clusters") +
  facet_wrap(~cluster) +
  facet_grid(rows = vars(cluster)) + 
  theme(legend.position="none") +
  xlab("") +
  theme(text = element_text(size=5),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1))
ggsave("plot/omics_profile3.pdf", width = 10, height = 7)

# this plot is visually misleading and not recommended for presenting the omics profiles
# # 2. try the circular barplot 
# # (these circular plot codes are copied from https://www.r-graph-gallery.com/297-circular-barplot-with-groups.html)
# # add metabolites type
# M_mean_melt$type = as.factor(c(rep("serum-metabolite", 24), rep("urine-metabolite", 4), rep("proteomic", 28)))
# # add some NA rows in order to create separation
# empty_bar <- 4
# to_add <- data.frame( matrix(NA, nlevels(M_mean_melt$cluster) * empty_bar*nlevels(M_mean_melt$type), ncol(M_mean_melt)) )
# colnames(to_add) <- colnames(M_mean_melt)
# to_add$type <- rep(rep(levels(M_mean_melt$type), each = empty_bar), nlevels(M_mean_melt$cluster))
# to_add$cluster = rep(levels(M_mean$cluster), each = nlevels(M_mean_melt$cluster) * nlevels(M_mean_melt$type))
# M_mean_melt <- rbind(M_mean_melt, to_add)
# M_mean_melt <- M_mean_melt %>% arrange(cluster, type)
# M_mean_melt$id <- rep(seq(1, nrow(M_mean_melt) / 4), 4)
# for(i in 1:4) {
#   #Get the name and the y position of each label
#   label_data <- M_mean_melt[M_mean_melt$cluster == i, ]
#   number_of_bar <- nrow(label_data)
#   angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
#   label_data$hjust <- ifelse( angle < -90, 1, 0)
#   label_data$angle <- ifelse(angle < -90, angle+180, angle)
#   # plot for the 4 clusters separately
#   ggplot(M_mean_melt[M_mean_melt$cluster == i, ], aes(x=as.factor(id), y=value * 100, fill=type)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
#     geom_bar(stat="identity", alpha=0.5) +
#     ylim(-100,120) +
#     theme_minimal() +
#     theme(
#       legend.position = "none",
#       axis.text = element_blank(),
#       axis.title = element_blank(),
#       panel.grid = element_blank(),
#       plot.margin = unit(rep(-1,4), "cm") 
#     ) +
#     coord_polar() + 
#     geom_text(data=label_data, 
#               aes(x=id, y=(value > 0) * value * 100 +10, label=variable, hjust=hjust), 
#               color="black", fontface="bold",alpha=0.6, size=2.5, 
#               angle= label_data$angle, inherit.aes = FALSE ) + 
#     ggtitle(paste0("Metabolites profiles for cluster ", i))
#   
#   ggsave(paste0("plot/omics_profile4-", i, ".pdf"), width = 5, height = 5)
# }


#### 5. exposure profiles ####
E_mean = as.data.frame(fit_try3$pars$beta[, -1])
E_mean$cluster = as.factor(1:4)
E_mean_melt = melt(E_mean, id.vars = "cluster")
ggplot(E_mean_melt, aes(fill = variable, y = value, x = variable)) + 
  geom_bar(position="dodge", stat="identity") + 
  ggtitle("Exposure profiles for 4 latent clusters") +
  facet_wrap(~cluster) +
  facet_grid(rows = vars(cluster)) + 
  theme(legend.position="none") +
  xlab("") +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1))
ggsave("plot/exposure_profile.pdf", width = 10, height = 7)



#### 6. bootstrap CI ####
n_boot = 500 # not many of the iterations don't converge. Such iteration is immediately terminated
n_obs = nrow(exposure)
res = vector(mode = "list", length = n_boot)
for (i in 1:500) {
  set.seed(i)
  index = sample(1:n_obs, n_obs, replace = TRUE)
  fit_temp = withTimeout({
    try(est.lucid(G = exposure[index, fit_try2$select$selectG], 
                  Z = metabolite[index, fit_try2$select$selectZ], 
                  Y = as.matrix(outcome[index, ]), CoY = covariate[index, ],
                  useY = FALSE, K = 4, family = "normal", 
                  control = def.control(max_itr = 400, max_tot.itr = 400)))
  }, timeout = 60, onTimeout = "warning")
  if("try-error" %in% class(fit_temp)) {
    next
  } else if(is.null(fit_temp)){
    next
  } else {
    res[[i]] = fit_temp$pars
  }
  print(i)
}

# summarize the result
# 1. order the parameters based on the increasing order of BMI z-score
result_sum = NULL
for(i in 1:n_boot) {
  temp_sum = res[[i]]
  if(is.null(temp_sum)) {
    next
  } else{
    temp_beta = temp_sum$beta[, -1]
    temp_mu = temp_sum$mu
    temp_gamma = temp_sum$gamma$beta
    index_par = order(temp_gamma[1:4])
    temp_beta = temp_beta[index_par, ]
    for(i in 1:4) {
      temp_beta[i, ] = temp_beta[i, ] - temp_beta[1, ]
    }
    temp_mu = temp_mu[index_par, ]
    temp_res = c(as.vector(temp_beta[-1, ]), as.vector(temp_mu), temp_gamma)
    names(temp_res) = c(as.vector(sapply(colnames(temp_beta), function(x) {paste0(x,"_cluster", 2:4)})),
                        as.vector(sapply(colnames(temp_mu), function(x) {paste0(x,"_cluster", 1:4)})),
                        "LC1", names(temp_gamma)[-1])
    result_sum = rbind(result_sum, temp_res)
  }
}

# 3. include 300 iterations
sd = apply(result_sum, 2, sd) # some exposure have very large effect estimate and is not reliable
pars = c(as.vector(fit_try3$pars$beta[-1, -1]), as.vector(fit_try3$pars$mu), fit_try3$pars$gamma$beta)
dat_pars = data.frame(pars = pars,
                      sd = sd)
dat_pars$cluster = c(rep(paste0("cluster", 2:4), ncol(temp_beta)), rep(paste0("cluster", 1:4), ncol(temp_mu)),  rep(NA, 14))
dat_pars$name = c(rep(colnames(temp_beta), each = 3), rep(colnames(temp_mu), each = 4), "LC1", names(temp_gamma)[-1])
dat_pars$type = c(rep("exposure", (ncol(temp_beta) * 3)), rep("serum", 4 * 6), rep("urine", 4), rep("protein", 28), rep(NA, 14))
write.csv(dat_pars, file = "result/boot_se.csv", quote = FALSE, row.names = FALSE)
dat_pars$name = factor(dat_pars$name, levels = c(colnames(temp_beta), "HGF", "IL1beta", "IL6", "IL8", "INSULIN", "Leptin",
                                                 "TNFalfa", "urine_metab_5", "serum_metab_63", "serum_metab_79", "serum_metab_81", "serum_metab_95", "serum_metab_122", "serum_metab_153"))

# point estimates for exposures
ggplot(data = dat_pars[1:(ncol(temp_beta) * 3), ], aes(x = name, y = pars)) + 
  geom_point() + 
  facet_wrap(~cluster) + 
  facet_grid(rows = vars(cluster)) + 
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1),
        legend.position = "none") + 
  geom_errorbar(aes(ymin = pars - 1.96 * sd, ymax = pars + 1.96 * sd), width=.2,
                position=position_dodge(0.05)) + 
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") + 
  xlab("") + 
  ylab("") + 
  ggtitle("Point estimate for expsoure effect with 95% CI") + 
  geom_point(data = na.omit(dat_pars[dat_pars$name == "hs_hcb_cadj_Log2" | dat_pars$name == "hs_hcb_madj_Log2", ]),
             aes(x = name, y = pars, color = "red"))
ggsave("plot/exposure_point_estimate.pdf", width = 10, height = 7)

# barplot for exposures
ggplot(rbind(dat_pars[1:(ncol(temp_beta) * 3), c(1, 3, 4)],
             data.frame(pars = rep(NA, ncol(temp_beta)),
                        cluster = rep("cluster1", ncol(temp_beta)),
                        name = colnames(temp_beta))), 
       aes(y = pars, x = name)) + 
  geom_bar(position="dodge", stat="identity") + 
  ggtitle("Exposure profiles for 4 latent clusters") +
  facet_wrap(~cluster) +
  facet_grid(rows = vars(cluster)) + 
  xlab("") +
  ylab("Estimated Mean") + 
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1),
        legend.title = element_blank())
ggsave("plot/omics_profile3.pdf", width = 10, height = 7)


# point estimate for metabolites
ggplot(data = dat_pars[(ncol(temp_beta) * 3 + 1):(ncol(temp_beta) * 3 + ncol(temp_mu) * 4), ],
       aes(x = name, y = pars, color = type, shape = type)) + 
  geom_point() + 
  facet_wrap(~cluster) + 
  facet_grid(rows = vars(cluster)) +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1)) +
  geom_errorbar(aes(ymin = pars - 1.96 * sd, ymax = pars + 1.96 * sd), width=.2,
                position=position_dodge(0.05)) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  xlab("") + 
  ylab("") + 
  ggtitle("Point estimate for metabolite effect with 95% CI") 
ggsave("plot/metabolite_point_estimate.pdf", width = 10, height = 7)


# barplot for omics data for each cluster
ggplot(dat_pars[(ncol(temp_beta) * 3 + 1):(ncol(temp_beta) * 3 + ncol(temp_mu) * 4), ], 
       aes(fill = type, y = pars, x = name)) + 
  geom_bar(position="dodge", stat="identity") + 
  ggtitle("Omics profiles for 4 latent clusters") +
  facet_wrap(~cluster) +
  facet_grid(rows = vars(cluster)) + 
  xlab("") +
  ylab("Estimated Mean") + 
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1),
        legend.title = element_blank())
ggsave("plot/omics_profile3.pdf", width = 10, height = 7)
