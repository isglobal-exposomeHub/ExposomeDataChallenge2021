### ISGlobal Exposome Challenge 2021 - Risk score computation ###
### This script computes environmental risk scores for all participants using
### the Elastic-net selected environmental exposure variables.
### Author: Farida S. Akhtari
### Date: 4/26/2021


# Setup

library(librarian)
librarian::shelf(
  car, corrplot, data.table, glmnet, naniar, rms, RColorBrewer, rstudioapi,
  tidyverse, wrapr
)

cwd <- dirname(rstudioapi::getActiveDocumentContext()$path)

DATA = "~/NIEHS/Exposome_Data_Challenge/DATA/"
CODE = "~/NIEHS/Exposome_Data_Challenge/CODE/"
RESULT = "~/NIEHS/Exposome_Data_Challenge/RESULT/"

load(paste0(RESULT,"exposure_select.RData"))


# Functions

## Compute risk scores
compute_risk_scores <- function(exp_data, betas) {
  risk_scores <- as.matrix(exp_data) %*% betas
  return(as.numeric(risk_scores))
}


# Total/Single risk score
risk_scores <- compute_risk_scores(exp_data = exposure, betas = coef_fit)
summary(risk_scores)
hist(risk_scores)
boxplot(risk_scores)

## Scale risk scores
risk_scores_scaled <- as.numeric(scale(risk_scores))
summary(risk_scores_scaled)
hist(risk_scores_scaled)
boxplot(risk_scores_scaled)

## Model the risk score
rs_df <- covariate %>%
  add_column(outcome = outcome, .before = 1) %>%
  add_column(risk_score = as.numeric(risk_scores_scaled))

glm_model <- glm(
  outcome ~ .,
  family = "binomial",
  data = rs_df,
  singular.ok = FALSE
)
summary(glm_model)
anova_out <- Anova(glm_model)
print(anova_out)


# Risk score per family/category
exp_families <- unique(codebook_expo_sele$family)
family_risk_score <- data.frame(matrix(vector(mode = 'numeric'), nrow = 1007))

for (ef in exp_families) {
  cat("Exposure family:", ef, "\n")
  
  sel_exps <- codebook_expo_sele %>% filter(family == ef) %>% pull(variable_name)
  # Since we have dummy variables in the exposure data, select all
  # variables/columns starting with selected exposure variable string
  exp_indx <- str_detect(colnames(exposure), paste(paste0("^", sel_exps), collapse = "|"))
  exp_data <- exposure[, exp_indx, drop = FALSE]
  print("Selected exposures:")
  print(colnames(exp_data))
  
  # Extract betas for selected exposures
  exp_betas <- coef_fit[exp_indx]
  print(length(exp_betas))
  family_risk_score[, paste0(ef, "_rs")] <-
    as.numeric(scale(compute_risk_scores(exp_data = exp_data, betas = exp_betas)))
}
family_risk_score <- family_risk_score %>% rename_all(make.names)

## Model the risk scores
all_rs_df <- cbind(covariate, family_risk_score) %>%
  add_column(outcome = outcome, .before = 1) %>%
  add_column(Total_rs = risk_scores_scaled)
colnames(all_rs_df)
write.table(all_rs_df, paste0(cwd, "/risk_score_data.txt"))

glm_model <- glm(
  outcome ~ .,
  family = "binomial",
  data = all_rs_df %>% select(-Total_rs),
  singular.ok = FALSE
)
summary(glm_model)
anova_out <- Anova(glm_model)
print(anova_out)

glm_op <- coef(summary(glm_model))
glm_op <- cbind(variable = rownames(glm_op), glm_op)
rownames(glm_op) <- NULL
glm_op <- as.data.frame(glm_op) %>%
  rename_all((make.names)) %>%
  setnames(old = c('Std..Error','Pr...z..'), new = c('Std.Error','pvalue')) %>%
  mutate_at(qc(Estimate, Std.Error, z.value, pvalue), as.numeric) %>%
  arrange(pvalue)
write.table(glm_op, paste0(cwd, "/risk_score_glmop.txt"))

## Significant results
glm_op %>% filter(pvalue < 0.05) %>% select(variable, pvalue)


# Make plots
plot_df <- all_rs_df %>%
  select(outcome, ends_with("_rs"))

pdf(paste0(cwd, "/rs_plots.pdf"))

rs_cols <- colnames(plot_df)[str_detect(colnames(plot_df), "_rs$")]

for (rs in rs_cols) {
  ## Risk score boxplots by outcome
  print(
    ggplot(data = plot_df, aes_string(x = "outcome", y = rs, fill = "outcome")) +
      geom_boxplot() +
      scale_fill_brewer(palette = "Paired") +
      ylab(rs) +
      xlab("Asthma Outcome") +
      theme_classic() +
      theme(legend.position = "none") +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14)
      )
  )
  
  outcome_df <- plot_df %>% select(all_of(c("outcome", rs)))
  
  ## Risk score percentile plots
  prev_ptile <- plot_df %>%
    select(all_of(c("outcome", rs))) %>%
    mutate(percentile = ntile(as.numeric(.data[[rs]]), 100)) %>%
    group_by(percentile) %>%
    summarise(prev = mean(as.numeric(levels(outcome))[outcome]))
  
  print(
    ggplot(
      data = prev_ptile,
      aes(x = percentile, y = prev * 100, colour = percentile)
    ) +
      geom_point() +
      scale_color_gradientn(colours = brewer.pal(n = 9, name = "Blues")) +
      ylab("Asthma prevalence %") +
      xlab(paste(gsub("_rs$", "", rs), "risk score percentile")) +
      theme_classic()
  )
  
  ## Risk score decile plots
  prev_decile <- plot_df %>%
    select(all_of(c("outcome", rs))) %>%
    mutate(decile = ntile(as.numeric(.data[[rs]]), 10)) %>%
    group_by(decile) %>%
    summarise(prev = mean(as.numeric(levels(outcome))[outcome]))
  
  print(
    ggplot(
      data = prev_decile,
      aes(x = as.factor(decile), y = prev * 100, fill = decile)) +
      geom_col() +
      scale_fill_distiller(direction = 1) +
      ylab("Asthma prevalence %") +
      xlab(paste(gsub("_rs$", "", rs), "risk score decile")) +
      theme_classic() +
      theme(
        #legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14)
      )
  )
}

dev.off()
