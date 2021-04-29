##avoid source
no_function()

sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())
source("code/tools.R")

##load data
###load exposome with phenotye
load("data_analysis/correlation_network/within_child/exposome_vs_phenotype/exposome_phenotype_glm")
exposome_phenotype_glm = 
  exposome_phenotype_glm %>% 
  dplyr::filter(p.adjust < 0.05)

###load exposome with internal-ome
load("data_analysis/correlation_network/within_child/exposome_vs_internal_omics/exposome_internal_omics_cor")
exposome_internal_omics_cor = 
  exposome_internal_omics_cor %>% 
  dplyr::filter(p.adjust < 0.05)

###load internal-ome with phenotype
load("data_analysis/correlation_network/within_child/internal_ome_phenotype/internal_ome_phenotype_glm")
internal_ome_phenotype_glm = 
  internal_ome_phenotype_glm %>% 
  dplyr::filter(p.adjust < 0.05)

###load internal-ome data
load("data_analysis/transcriptome_data_analysiss/data_preparation/variable_info")
load("data_analysis/transcriptome_data_analysiss/data_preparation/expression_data")
load("data_analysis/transcriptome_data_analysiss/data_preparation/sample_info")
transcriptome_variable_info = variable_info
transcriptome_expression_data = expression_data
transcriptome_sample_info = sample_info

load("data_analysis/proteome_data_analysiss/data_preparation/variable_info")
load("data_analysis/proteome_data_analysiss/data_preparation/expression_data")
load("data_analysis/proteome_data_analysiss/data_preparation/sample_info")
proteome_variable_info = variable_info
proteome_expression_data = expression_data
proteome_sample_info = sample_info

load("data_analysis/serum_metabolome_data_analysiss/data_preparation/variable_info")
load("data_analysis/serum_metabolome_data_analysiss/data_preparation/expression_data")
load("data_analysis/serum_metabolome_data_analysiss/data_preparation/sample_info")
serum_metabolome_variable_info = variable_info
serum_metabolome_expression_data = expression_data
serum_metabolome_sample_info = sample_info

load("data_analysis/urine_metabolome_data_analysiss/data_preparation/variable_info")
load("data_analysis/urine_metabolome_data_analysiss/data_preparation/expression_data")
load("data_analysis/urine_metabolome_data_analysiss/data_preparation/sample_info")
urine_metabolome_variable_info = variable_info
urine_metabolome_expression_data = expression_data
urine_metabolome_sample_info = sample_info

###combine internal_ome data
intersect_name =
  Reduce(f = intersect, 
         x = list(colnames(transcriptome_expression_data),
                  colnames(proteome_expression_data),
                  colnames(serum_metabolome_expression_data),
                  colnames(urine_metabolome_expression_data)))

internal_ome_expression_data = 
  rbind(transcriptome_expression_data[,intersect_name],
        proteome_expression_data[,intersect_name],
        serum_metabolome_expression_data[,intersect_name],
        urine_metabolome_expression_data[,intersect_name])

internal_ome_variable_info = 
  rbind(transcriptome_variable_info[,c("variable_id", "GeneSymbol_Affy")] %>% 
          dplyr::rename(true_name = GeneSymbol_Affy) %>% 
          dplyr::mutate(class = "Transcriptome"),
        proteome_variable_info[,c("variable_id", "Gene_Symbol")] %>% 
          dplyr::rename(true_name = Gene_Symbol) %>% 
          dplyr::mutate(class = "Proteome"),
        serum_metabolome_variable_info[,c("variable_id", "var")] %>% 
          dplyr::rename(true_name = var) %>% 
          dplyr::mutate(class = "Serum_metabolome"),
        urine_metabolome_variable_info[,c("variable_id", "var")] %>% 
          dplyr::rename(true_name = var) %>% 
          dplyr::mutate(class = "Urine_metabolome"))

internal_ome_sample_info = 
  rbind(
    transcriptome_sample_info[,c("sample_id", "subject_id")]
  ) %>% 
  dplyr::filter(sample_id %in% colnames(internal_ome_expression_data))

dim(internal_ome_expression_data)

sum(colnames(internal_ome_expression_data) == internal_ome_sample_info$sample_id)
sum(rownames(internal_ome_expression_data) == internal_ome_variable_info$variable_id)


###load exposome data
load("data_analysis/exposome_air_data_analysis/data_preparation/variable_info")
load("data_analysis/exposome_air_data_analysis/data_preparation/expression_data")
load("data_analysis/exposome_air_data_analysis/data_preparation/sample_info")
exposome_air_variable_info = variable_info
exposome_air_sample_info = sample_info
exposome_air_expression_data = expression_data

load("data_analysis/exposome_outdoor_data_analysis/data_preparation/variable_info")
load("data_analysis/exposome_outdoor_data_analysis/data_preparation/expression_data")
load("data_analysis/exposome_outdoor_data_analysis/data_preparation/sample_info")
exposome_outdoor_variable_info = variable_info
exposome_outdoor_sample_info = sample_info
exposome_outdoor_expression_data = expression_data

load("data_analysis/exposome_chemical_data_analysis/data_preparation/variable_info")
load("data_analysis/exposome_chemical_data_analysis/data_preparation/expression_data")
load("data_analysis/exposome_chemical_data_analysis/data_preparation/sample_info")
exposome_chemical_variable_info = variable_info
exposome_chemical_sample_info = sample_info
exposome_chemical_expression_data = expression_data

###combine exposome data
colnames(exposome_air_expression_data)
colnames(exposome_outdoor_expression_data)
colnames(exposome_chemical_expression_data)

intersect_name = Reduce(f = intersect, x = list(
  colnames(exposome_air_expression_data),
  colnames(exposome_outdoor_expression_data),
  colnames(exposome_chemical_expression_data)
))

exposome_expression_data =
  rbind(
    exposome_air_expression_data[, intersect_name],
    exposome_outdoor_expression_data[, intersect_name],
    exposome_chemical_expression_data[, intersect_name]
  )

exposome_sample_info = exposome_air_sample_info %>% 
  dplyr::filter(sample_id %in% intersect_name)

exposome_variable_info =
  rbind(
    exposome_air_variable_info,
    exposome_outdoor_variable_info,
    exposome_chemical_variable_info
  )

sum(colnames(exposome_expression_data) == exposome_sample_info$sample_id)
sum(rownames(exposome_expression_data) == exposome_variable_info$variable_id)

dim(exposome_expression_data)



######only remain the overlapped samples in exposome and internal omics data
intersect_name = 
  intersect(colnames(exposome_expression_data),
            colnames(internal_ome_expression_data))

exposome_expression_data = 
  exposome_expression_data[,intersect_name]

exposome_sample_info = 
exposome_sample_info %>% 
  dplyr::filter(sample_id %in% intersect_name)

internal_ome_expression_data = 
internal_ome_expression_data[,intersect_name]

internal_ome_sample_info = 
  internal_ome_sample_info %>% 
  dplyr::filter(sample_id %in% intersect_name)

#####only remain the exposome which are associated with phenotype and internal-ome
interact_name = 
intersect(unique(exposome_phenotype_glm$variable_id),
          unique(exposome_internal_omics_cor$from))

exposome_phenotype_glm = 
exposome_phenotype_glm %>% 
  dplyr::filter(variable_id %in% interact_name)

exposome_internal_omics_cor = 
  exposome_internal_omics_cor %>% 
  dplyr::filter(from %in% interact_name)


dim(exposome_phenotype_glm)
dim(exposome_internal_omics_cor)

###only remain the internal_ome which are also associated with phenotye
exposome_internal_omics_cor = 
exposome_internal_omics_cor %>% 
  dplyr::filter(to %in% unique(internal_ome_phenotype_glm$variable_id))

#####mediantion analysis
####internla-ome ~ exposome (79) + covariates
###phenotype ~ exposome + internal_ome + covariates

setwd("data_analysis/correlation_network/within_child/mediation_analysis")

library(mediation)

# mediation_result = NULL
# 
# for(i in 4091:nrow(exposome_internal_omics_cor)) {
#   cat(i, "\n")
#   exposome_variable_id = exposome_internal_omics_cor$from[i]
#   internal_ome_variable_id = exposome_internal_omics_cor$to[i]
# 
#   temp_data =   
#   data.frame(exposome_value = as.numeric(exposome_expression_data[exposome_variable_id,]),
#              internal_ome_value = as.numeric(internal_ome_expression_data[internal_ome_variable_id,]),
#              exposome_sample_info
#              )
#   
#   glm_reg1 =
#     glm(
#       internal_ome_value ~ exposome_value + Child.sex + Year.of.birth + Maternal.BMI + Gestational.age.at.birth +
#         Maternal.age + Child.height + Child.weight + Birthweight,
#       family = gaussian,
#       temp_data
#     )
#   
#   temp_exposome_phenotype_glm = 
#     exposome_phenotype_glm %>% 
#     dplyr::filter(variable_id %in% exposome_variable_id)
#   
#   if(nrow(temp_exposome_phenotype_glm) == 0){
#     next()
#   }
#   
#   for (j in 1:nrow(temp_exposome_phenotype_glm)) {
#     cat(j, " ")
#     phenotype_id = temp_exposome_phenotype_glm$phenotype[j]
#   
#     temp_data2 = 
#       data.frame(phenotype_value = temp_data[,phenotype_id],
#                  temp_data)
#     
#     glm_reg2 =
#       glm(
#         phenotype_value ~ exposome_value + internal_ome_value + Child.sex + Year.of.birth + 
#           Maternal.BMI + Gestational.age.at.birth +
#           Maternal.age + Child.height + Child.weight + Birthweight,
#         family = gaussian,
#         temp_data2
#       )
#     
#     result = mediate(
#       glm_reg1 ,
#       glm_reg2,
#       treat = "exposome_value",
#       mediator = "internal_ome_value",
#       boot = F,
#       sims = 1000
#     )
#   
#   return_result = list(exposome_variable_id = exposome_variable_id,
#                        internal_ome_variable_id = internal_ome_variable_id,
#                        phenotype_id, 
#                        result = result)              
#   mediation_result = c(mediation_result, return_result)
#   }
# }
# 
# save(mediation_result, file = "mediation_result")
# load("mediation_result")
# length(mediation_result)
# 
# mediate =
#   purrr::map(4 * (1:(length(mediation_result) / 4)), function(idx) {
#     phenotype = mediation_result[[idx - 1]]
#     mediator = mediation_result[[idx - 2]]
#     treat = mediation_result[[idx - 3]]
#     if(phenotype == "Behavior" & mediator == "PAI1" & treat == "NDVI_home"){
#       cat(idx)
#     }
#     
#     results = summary(mediation_result[[idx]])
#     acme = results$d1
#     acme_ci_lower = results$d1.ci[1]
#     acme_ci_upper = results$d1.ci[2]
#     acme_p = results$d1.p
#     c(phenotype = phenotype, 
#       mediator = mediator, 
#       treat = treat, 
#       acme = acme, 
#       acme_ci_lower = acme_ci_lower, 
#       acme_ci_upper = acme_ci_upper, 
#       acme_p = acme_p)
#   }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame()
# 
# mediate$acme = as.numeric(mediate$acme)
# mediate$`acme_ci_lower.2.5%` = as.numeric(mediate$`acme_ci_lower.2.5%`)
# mediate$`acme_ci_upper.97.5%` = as.numeric(mediate$`acme_ci_upper.97.5%`)
# mediate$acme_p = as.numeric(mediate$acme_p)
# 
# mediate = 
# mediate %>% 
#   dplyr::filter(acme_p < 0.05 & acme > 0 & acme < 1 & 
#                   `acme_ci_upper.97.5%` < 1 & `acme_ci_lower.2.5%` > 0)
# 
# save(mediate, file = "mediate")
# 
# load("mediate")
# 
# ###get the correlation between treat with mediator, mediator and phenotype and treat with phenoty
# cor_value = 
# t(mediate) %>% 
#   as.data.frame() %>% 
# purrr::map(function(x){
#  ##treat vs mediator
#    cor_treat_mediator_test = 
#   cor.test(x = as.numeric(exposome_expression_data[x[3],]),
#            y = as.numeric(internal_ome_expression_data[x[2],]), 
#            method = "spearman")
#   cor_treat_mediator = cor_treat_mediator_test$estimate
#   cor_treat_mediator_p = cor_treat_mediator_test$p.value
#   
#   ##mediator vs phenotype
#   cor_mediator_phenotype_test = 
#     cor.test(x = as.numeric(exposome_sample_info[,x[1]]),
#              y = as.numeric(internal_ome_expression_data[x[2],]), 
#              method = "spearman")
#   cor_mediator_phenotype = cor_mediator_phenotype_test$estimate
#   cor_mediator_phenotype_p = cor_mediator_phenotype_test$p.value
#   
#   ##treat vs phenotype
#   cor_treat_phenotype_test = 
#     cor.test(x = as.numeric(exposome_sample_info[,x[1]]),
#              y = as.numeric(exposome_expression_data[x[3],]), 
#              method = "spearman")
#   cor_treat_phenotype = cor_treat_phenotype_test$estimate
#   cor_treat_phenotype_p = cor_treat_phenotype_test$p.value
#   
#   c(cor_treat_mediator = cor_treat_mediator,
#     cor_treat_mediator_p = cor_treat_mediator_p,
#     cor_mediator_phenotype = cor_mediator_phenotype,
#     cor_mediator_phenotype_p = cor_mediator_phenotype_p,
#     cor_treat_phenotype = cor_treat_phenotype,
#     cor_treat_phenotype_p = cor_treat_phenotype_p
#     )
# }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame()
# 
# mediate_result =
#   data.frame(mediate, cor_value)
# 
# ##only remain the relations with significant cor
# mediate_result = 
# mediate_result %>% 
#   dplyr::filter(cor_treat_mediator_p < 0.05 & 
#                   cor_mediator_phenotype_p < 0.05 &
#                   cor_treat_phenotype_p < 0.05) %>% 
#   dplyr::filter(((cor_treat_mediator.rho * cor_mediator_phenotype.rho) > 0 & cor_treat_phenotype.rho > 0) |
#                   ((cor_treat_mediator.rho * cor_mediator_phenotype.rho) < 0 & cor_treat_phenotype.rho < 0)) %>% 
#   dplyr::arrange(desc(acme))
# 
# save(mediate_result, file = "mediate_result")

load("mediate_result")
length(unique(mediate_result$phenotype))
length(unique(mediate_result$mediator))
length(unique(mediate_result$treat))


library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Arial Narrow")
addWorksheet(wb, sheetName = "BMI", gridLines = TRUE)
addWorksheet(wb, sheetName = "Behavior", gridLines = TRUE)
addWorksheet(wb, sheetName = "IQ", gridLines = TRUE)

bmi = 
  mediate_result %>% 
  dplyr::filter(phenotype == "Body.mass.index.z.score") %>% 
  dplyr::arrange(desc(acme)) %>% 
  dplyr::select(phenotype, mediator, treat, acme) %>% 
  dplyr::left_join(internal_ome_variable_info, by = c("mediator" = "variable_id")) %>% 
  dplyr::rename(mediator_name = true_name, mediator_class = class) %>% 
  dplyr::left_join(exposome_variable_info[,c("variable_id", "description", "domain")], by = c("treat" = "variable_id")) %>% 
  dplyr::rename(treat_name = description, treat_class = domain)

behavior = 
  mediate_result %>% 
  dplyr::filter(phenotype == "Behavior") %>% 
  dplyr::arrange(desc(acme)) %>% 
  dplyr::select(phenotype, mediator, treat, acme) %>% 
  dplyr::left_join(internal_ome_variable_info, by = c("mediator" = "variable_id")) %>% 
  dplyr::rename(mediator_name = true_name, mediator_class = class) %>% 
  dplyr::left_join(exposome_variable_info[,c("variable_id", "description", "domain")], by = c("treat" = "variable_id")) %>% 
  dplyr::rename(treat_name = description, treat_class = domain)

iq = 
  mediate_result %>% 
  dplyr::filter(phenotype == "Intelligence.quotient") %>% 
  dplyr::arrange(desc(acme)) %>% 
  dplyr::select(phenotype, mediator, treat, acme) %>% 
  dplyr::left_join(internal_ome_variable_info, by = c("mediator" = "variable_id")) %>% 
  dplyr::rename(mediator_name = true_name, mediator_class = class) %>% 
  dplyr::left_join(exposome_variable_info[,c("variable_id", "description", "domain")], by = c("treat" = "variable_id")) %>% 
  dplyr::rename(treat_name = description, treat_class = domain)
  

freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) ## freeze first row and column
freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE) ## freeze first row and column
freezePane(wb, sheet = 3, firstRow = TRUE, firstCol = TRUE) ## freeze first row and column
writeDataTable(wb, sheet = 1, x = bmi,colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 2, x = behavior, colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 3, x = iq, colNames = TRUE, rowNames = FALSE)
saveWorkbook(wb, "mediation_analysis_result.xlsx", overwrite = TRUE) ## save to working directory

mediate_result = 
mediate_result %>% 
  dplyr::arrange(desc(acme)) %>% 
  head(250)

#### network to show how exposome affect phenotypes via omics
edge_data = 
  rbind(mediate_result[,c("treat", "mediator", "cor_treat_mediator.rho", "cor_treat_mediator_p")] %>% 
          dplyr::rename(from = treat, to = mediator, cor = cor_treat_mediator.rho, p = cor_treat_mediator_p),
        mediate_result[,c("mediator", "phenotype", "cor_mediator_phenotype.rho", "cor_mediator_phenotype_p")] %>% 
          dplyr::rename(from = mediator, to = phenotype, cor = cor_mediator_phenotype.rho, p = cor_mediator_phenotype_p)
        )

node_data = 
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>% 
  dplyr::mutate(
    class1 = case_when(
      node %in% exposome_air_variable_info$variable_id ~ "Exposome_air",
      node %in% exposome_outdoor_variable_info$variable_id ~ "Exposome_outdoor",
      node %in% exposome_chemical_variable_info$variable_id ~ "Exposome_chemical",
      node %in% transcriptome_variable_info$variable_id ~ "Transcriptome",
      node %in% proteome_variable_info$variable_id ~ "Proteome",
      node %in% serum_metabolome_variable_info$variable_id ~ "Serum_metabolome",
      node %in% transcriptome_variable_info$variable_id ~ "Transcriptome",
      node %in% urine_metabolome_variable_info$variable_id ~ "Urine_metabolome",
      node %in% serum_metabolome_variable_info$variable_id ~ "Serum_metabolome",
      node %in% proteome_variable_info$variable_id ~ "Proteome",
      node %in% "Behavior" ~ "Behavior",
      node %in% "Intelligence.quotient" ~ "IQ",
      node %in% "Body.mass.index.z.score" ~ "BMI"
    )
  ) %>% 
  dplyr::mutate(class2 = 
                  case_when(
                    class1 %in% c("Exposome_air", "Exposome_outdoor", "Exposome_chemical") ~ "Exposome",
                    class1 %in% c("Transcriptome", "Proteome", "Urine_metabolome", "Serum_metabolome") ~ "Internal_ome",
                    class1 %in% c("Behavior", "IQ", "BMI") ~ "Phenotype",
                  )) %>% 
  dplyr::select(node, class1, class2) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::left_join(transcriptome_variable_info[, c("variable_id", "GeneSymbol_Affy")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(GeneSymbol_Affy) ~ GeneSymbol_Affy,
                                      TRUE ~ node)) %>% 
  dplyr::select(-GeneSymbol_Affy) %>% 
  dplyr::left_join(proteome_variable_info[, c("variable_id", "Gene_Symbol")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(Gene_Symbol) ~ Gene_Symbol,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-Gene_Symbol) %>% 
  dplyr::left_join(serum_metabolome_variable_info[, c("variable_id", "var")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(var) ~ var,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-var) %>% 
  dplyr::left_join(urine_metabolome_variable_info[, c("variable_id", "var")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(var) ~ var,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-var)  %>% 
  dplyr::left_join(exposome_air_variable_info[,c("variable_id", "description")],
                   by = c("node" = "variable_id")) %>% 
  dplyr::mutate(true_name = case_when(!is.na(description) ~ description,
                                      TRUE ~ true_name)) %>%
  dplyr::select(-description) %>% 
  dplyr::left_join(exposome_outdoor_variable_info[, c("variable_id", "labels")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(labels) ~ labels,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-labels) %>% 
  dplyr::left_join(exposome_chemical_variable_info[, c("variable_id", "labels")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(labels) ~ labels,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-labels)

library(tidygraph)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))


pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

library(ggraph)

my_layout <- create_layout(temp_data, 
                           layout = 'linear')

my_layout$y[my_layout$class2 == "Internal_ome"] <- 5
my_layout$y[my_layout$class2 == "Phenotype"] <- 10

my_layout1 <-
  my_layout

my_layout1$x <-
  my_layout$y

my_layout1$y <-
  my_layout$x

my_layout1$y[my_layout1$class2 == "Exposome"] <-
  seq(from = 1, to = 100, length.out = sum(my_layout1$class2 == "Exposome"))

my_layout1$y[my_layout1$class2 == "Internal_ome"] <-
  seq(from = 1, to = 100, length.out = sum(my_layout1$class2 == "Internal_ome"))

my_layout1$y[my_layout1$class2 == "Phenotype"] <-
  my_layout1$y[my_layout1$class2 == "Phenotype"] <-
  seq(from = 30, to = 70, length.out = sum(my_layout1$class2 == "Phenotype"))

plot <-
  ggraph(my_layout1) +
  geom_edge_link(aes(color = cor),
                 show.legend = TRUE, alpha = 0.5) +
  geom_node_point(aes(fill = class1,
                      size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  scale_fill_manual(
    values = c(omics_color, phenotype_color)
  ) +
  scale_color_manual(
    values = c(omics_color, phenotype_color)
  ) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1,
      label = true_name,
      hjust = ifelse(class2 == "Internal_ome", 1, 0),
      size = 3,
      colour = class1
    ),
    size = 3,
    alpha = 1, 
    show.legend = FALSE
  ) +
  guides(edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
         fill = guide_legend(title = "Class", 
                             override.aes = list(size = 4, linetype = "blank")),
         size = guide_legend(title = "Degree", override.aes = list(linetype = 0))) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
    # legend.position = c(1,0), legend.justification = c(1,0)
  )

plot

# ggsave(plot, filename = "mediation_analysis.pdf", width = 15, height = 7)


###some examples
####only for phenotype IQ
load("mediate_result")
length(unique(mediate_result$phenotype))
length(unique(mediate_result$mediator))
length(unique(mediate_result$treat))

mediate_result = 
  mediate_result %>% 
  dplyr::filter(phenotype == "Intelligence.quotient") %>% 
  dplyr::arrange(desc(acme)) %>% 
  head(100)

#### network to show how exposome affect phenotypes via omics
edge_data = 
  rbind(mediate_result[,c("treat", "mediator", "cor_treat_mediator.rho", "cor_treat_mediator_p")] %>% 
          dplyr::rename(from = treat, to = mediator, cor = cor_treat_mediator.rho, p = cor_treat_mediator_p),
        mediate_result[,c("mediator", "phenotype", "cor_mediator_phenotype.rho", "cor_mediator_phenotype_p")] %>% 
          dplyr::rename(from = mediator, to = phenotype, cor = cor_mediator_phenotype.rho, p = cor_mediator_phenotype_p)
  )

node_data = 
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>% 
  dplyr::mutate(
    class1 = case_when(
      node %in% exposome_air_variable_info$variable_id ~ "Exposome_air",
      node %in% exposome_outdoor_variable_info$variable_id ~ "Exposome_outdoor",
      node %in% exposome_chemical_variable_info$variable_id ~ "Exposome_chemical",
      node %in% transcriptome_variable_info$variable_id ~ "Transcriptome",
      node %in% proteome_variable_info$variable_id ~ "Proteome",
      node %in% serum_metabolome_variable_info$variable_id ~ "Serum_metabolome",
      node %in% transcriptome_variable_info$variable_id ~ "Transcriptome",
      node %in% urine_metabolome_variable_info$variable_id ~ "Urine_metabolome",
      node %in% serum_metabolome_variable_info$variable_id ~ "Serum_metabolome",
      node %in% proteome_variable_info$variable_id ~ "Proteome",
      node %in% "Behavior" ~ "Behavior",
      node %in% "Intelligence.quotient" ~ "IQ",
      node %in% "Body.mass.index.z.score" ~ "BMI"
    )
  ) %>% 
  dplyr::mutate(class2 = 
                  case_when(
                    class1 %in% c("Exposome_air", "Exposome_outdoor", "Exposome_chemical") ~ "Exposome",
                    class1 %in% c("Transcriptome", "Proteome", "Urine_metabolome", "Serum_metabolome") ~ "Internal_ome",
                    class1 %in% c("Behavior", "IQ", "BMI") ~ "Phenotype",
                  )) %>% 
  dplyr::select(node, class1, class2) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::left_join(transcriptome_variable_info[, c("variable_id", "GeneSymbol_Affy")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(GeneSymbol_Affy) ~ GeneSymbol_Affy,
                                      TRUE ~ node)) %>% 
  dplyr::select(-GeneSymbol_Affy) %>% 
  dplyr::left_join(proteome_variable_info[, c("variable_id", "Gene_Symbol")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(Gene_Symbol) ~ Gene_Symbol,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-Gene_Symbol) %>% 
  dplyr::left_join(serum_metabolome_variable_info[, c("variable_id", "var")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(var) ~ var,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-var) %>% 
  dplyr::left_join(urine_metabolome_variable_info[, c("variable_id", "var")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(var) ~ var,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-var)  %>% 
  dplyr::left_join(exposome_air_variable_info[,c("variable_id", "description")],
                   by = c("node" = "variable_id")) %>% 
  dplyr::mutate(true_name = case_when(!is.na(description) ~ description,
                                      TRUE ~ true_name)) %>%
  dplyr::select(-description) %>% 
  dplyr::left_join(exposome_outdoor_variable_info[, c("variable_id", "labels")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(labels) ~ labels,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-labels) %>% 
  dplyr::left_join(exposome_chemical_variable_info[, c("variable_id", "labels")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(labels) ~ labels,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-labels)

library(tidygraph)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))


pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

library(ggraph)

my_layout <- create_layout(temp_data, 
                           layout = 'linear')

my_layout$y[my_layout$class2 == "Internal_ome"] <- 5
my_layout$y[my_layout$class2 == "Phenotype"] <- 10

my_layout1 <-
  my_layout

my_layout1$x <-
  my_layout$y

my_layout1$y <-
  my_layout$x

my_layout1$y[my_layout1$class2 == "Exposome"] <-
  seq(from = 1, to = 100, length.out = sum(my_layout1$class2 == "Exposome"))

my_layout1$y[my_layout1$class2 == "Internal_ome"] <-
  seq(from = 1, to = 100, length.out = sum(my_layout1$class2 == "Internal_ome"))

my_layout1$y[my_layout1$class2 == "Phenotype"] <- 50

plot <-
  ggraph(my_layout1) +
  geom_edge_link(aes(color = cor),
                 show.legend = TRUE, alpha = 0.5) +
  geom_node_point(aes(fill = class1,
                      size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  scale_fill_manual(
    values = c(omics_color, phenotype_color)
  ) +
  scale_color_manual(
    values = c(omics_color, phenotype_color)
  ) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1,
      label = true_name,
      hjust = ifelse(class2 == "Internal_ome", 1, 0),
      size = 3,
      colour = class1
    ),
    size = 3,
    alpha = 1, 
    show.legend = FALSE
  ) +
  guides(edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
         fill = guide_legend(title = "Class", 
                             override.aes = list(size = 4, linetype = "blank")),
         size = guide_legend(title = "Degree", override.aes = list(linetype = 0))) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
    # legend.position = c(1,0), legend.justification = c(1,0)
  ) 

plot

# ggsave(plot, filename = "IQ_mediation_analysis.pdf", width = 7, height = 7)















####only for phenotype Behavior
load("mediate_result")
length(unique(mediate_result$phenotype))
length(unique(mediate_result$mediator))
length(unique(mediate_result$treat))

mediate_result = 
  mediate_result %>% 
  dplyr::filter(phenotype == "Behavior") %>% 
  dplyr::arrange(desc(acme)) %>% 
  head(100)

#### network to show how exposome affect phenotypes via omics
edge_data = 
  rbind(mediate_result[,c("treat", "mediator", "cor_treat_mediator.rho", "cor_treat_mediator_p")] %>% 
          dplyr::rename(from = treat, to = mediator, cor = cor_treat_mediator.rho, p = cor_treat_mediator_p),
        mediate_result[,c("mediator", "phenotype", "cor_mediator_phenotype.rho", "cor_mediator_phenotype_p")] %>% 
          dplyr::rename(from = mediator, to = phenotype, cor = cor_mediator_phenotype.rho, p = cor_mediator_phenotype_p)
  )

node_data = 
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>% 
  dplyr::mutate(
    class1 = case_when(
      node %in% exposome_air_variable_info$variable_id ~ "Exposome_air",
      node %in% exposome_outdoor_variable_info$variable_id ~ "Exposome_outdoor",
      node %in% exposome_chemical_variable_info$variable_id ~ "Exposome_chemical",
      node %in% transcriptome_variable_info$variable_id ~ "Transcriptome",
      node %in% proteome_variable_info$variable_id ~ "Proteome",
      node %in% serum_metabolome_variable_info$variable_id ~ "Serum_metabolome",
      node %in% transcriptome_variable_info$variable_id ~ "Transcriptome",
      node %in% urine_metabolome_variable_info$variable_id ~ "Urine_metabolome",
      node %in% serum_metabolome_variable_info$variable_id ~ "Serum_metabolome",
      node %in% proteome_variable_info$variable_id ~ "Proteome",
      node %in% "Behavior" ~ "Behavior",
      node %in% "Intelligence.quotient" ~ "IQ",
      node %in% "Body.mass.index.z.score" ~ "BMI"
    )
  ) %>% 
  dplyr::mutate(class2 = 
                  case_when(
                    class1 %in% c("Exposome_air", "Exposome_outdoor", "Exposome_chemical") ~ "Exposome",
                    class1 %in% c("Transcriptome", "Proteome", "Urine_metabolome", "Serum_metabolome") ~ "Internal_ome",
                    class1 %in% c("Behavior", "IQ", "BMI") ~ "Phenotype",
                  )) %>% 
  dplyr::select(node, class1, class2) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::left_join(transcriptome_variable_info[, c("variable_id", "GeneSymbol_Affy")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(GeneSymbol_Affy) ~ GeneSymbol_Affy,
                                      TRUE ~ node)) %>% 
  dplyr::select(-GeneSymbol_Affy) %>% 
  dplyr::left_join(proteome_variable_info[, c("variable_id", "Gene_Symbol")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(Gene_Symbol) ~ Gene_Symbol,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-Gene_Symbol) %>% 
  dplyr::left_join(serum_metabolome_variable_info[, c("variable_id", "var")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(var) ~ var,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-var) %>% 
  dplyr::left_join(urine_metabolome_variable_info[, c("variable_id", "var")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(var) ~ var,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-var)  %>% 
  dplyr::left_join(exposome_air_variable_info[,c("variable_id", "description")],
                   by = c("node" = "variable_id")) %>% 
  dplyr::mutate(true_name = case_when(!is.na(description) ~ description,
                                      TRUE ~ true_name)) %>%
  dplyr::select(-description) %>% 
  dplyr::left_join(exposome_outdoor_variable_info[, c("variable_id", "labels")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(labels) ~ labels,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-labels) %>% 
  dplyr::left_join(exposome_chemical_variable_info[, c("variable_id", "labels")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(labels) ~ labels,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-labels)

library(tidygraph)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))


pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

library(ggraph)

my_layout <- create_layout(temp_data, 
                           layout = 'linear')

my_layout$y[my_layout$class2 == "Internal_ome"] <- 5
my_layout$y[my_layout$class2 == "Phenotype"] <- 10

my_layout1 <-
  my_layout

my_layout1$x <-
  my_layout$y

my_layout1$y <-
  my_layout$x

my_layout1$y[my_layout1$class2 == "Exposome"] <-
  seq(from = 1, to = 100, length.out = sum(my_layout1$class2 == "Exposome"))

my_layout1$y[my_layout1$class2 == "Internal_ome"] <-
  seq(from = 1, to = 100, length.out = sum(my_layout1$class2 == "Internal_ome"))

my_layout1$y[my_layout1$class2 == "Phenotype"] <- 50

plot <-
  ggraph(my_layout1) +
  geom_edge_link(aes(color = cor),
                 show.legend = TRUE, alpha = 0.5) +
  geom_node_point(aes(fill = class1,
                      size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  scale_fill_manual(
    values = c(omics_color, phenotype_color)
  ) +
  scale_color_manual(
    values = c(omics_color, phenotype_color)
  ) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1,
      label = true_name,
      hjust = ifelse(class2 == "Internal_ome", 1, 0),
      size = 3,
      colour = class1
    ),
    size = 3,
    alpha = 1, 
    show.legend = FALSE
  ) +
  guides(edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
         fill = guide_legend(title = "Class", 
                             override.aes = list(size = 4, linetype = "blank")),
         size = guide_legend(title = "Degree", override.aes = list(linetype = 0))) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
    # legend.position = c(1,0), legend.justification = c(1,0)
  ) 

plot

# ggsave(plot, filename = "Behavior_mediation_analysis.pdf", width = 7, height = 7)
# 
# 




####only for phenotype BMI
load("mediate_result")
length(unique(mediate_result$phenotype))
length(unique(mediate_result$mediator))
length(unique(mediate_result$treat))

mediate_result = 
  mediate_result %>% 
  dplyr::filter(phenotype == "Body.mass.index.z.score") %>% 
  dplyr::arrange(desc(acme)) %>% 
  head(100)

#### network to show how exposome affect phenotypes via omics
edge_data = 
  rbind(mediate_result[,c("treat", "mediator", "cor_treat_mediator.rho", "cor_treat_mediator_p")] %>% 
          dplyr::rename(from = treat, to = mediator, cor = cor_treat_mediator.rho, p = cor_treat_mediator_p),
        mediate_result[,c("mediator", "phenotype", "cor_mediator_phenotype.rho", "cor_mediator_phenotype_p")] %>% 
          dplyr::rename(from = mediator, to = phenotype, cor = cor_mediator_phenotype.rho, p = cor_mediator_phenotype_p)
  )

node_data = 
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>% 
  dplyr::mutate(
    class1 = case_when(
      node %in% exposome_air_variable_info$variable_id ~ "Exposome_air",
      node %in% exposome_outdoor_variable_info$variable_id ~ "Exposome_outdoor",
      node %in% exposome_chemical_variable_info$variable_id ~ "Exposome_chemical",
      node %in% transcriptome_variable_info$variable_id ~ "Transcriptome",
      node %in% proteome_variable_info$variable_id ~ "Proteome",
      node %in% serum_metabolome_variable_info$variable_id ~ "Serum_metabolome",
      node %in% transcriptome_variable_info$variable_id ~ "Transcriptome",
      node %in% urine_metabolome_variable_info$variable_id ~ "Urine_metabolome",
      node %in% serum_metabolome_variable_info$variable_id ~ "Serum_metabolome",
      node %in% proteome_variable_info$variable_id ~ "Proteome",
      node %in% "Behavior" ~ "Behavior",
      node %in% "Intelligence.quotient" ~ "IQ",
      node %in% "Body.mass.index.z.score" ~ "BMI"
    )
  ) %>% 
  dplyr::mutate(class2 = 
                  case_when(
                    class1 %in% c("Exposome_air", "Exposome_outdoor", "Exposome_chemical") ~ "Exposome",
                    class1 %in% c("Transcriptome", "Proteome", "Urine_metabolome", "Serum_metabolome") ~ "Internal_ome",
                    class1 %in% c("Behavior", "IQ", "BMI") ~ "Phenotype",
                  )) %>% 
  dplyr::select(node, class1, class2) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::left_join(transcriptome_variable_info[, c("variable_id", "GeneSymbol_Affy")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(GeneSymbol_Affy) ~ GeneSymbol_Affy,
                                      TRUE ~ node)) %>% 
  dplyr::select(-GeneSymbol_Affy) %>% 
  dplyr::left_join(proteome_variable_info[, c("variable_id", "Gene_Symbol")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(Gene_Symbol) ~ Gene_Symbol,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-Gene_Symbol) %>% 
  dplyr::left_join(serum_metabolome_variable_info[, c("variable_id", "var")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(var) ~ var,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-var) %>% 
  dplyr::left_join(urine_metabolome_variable_info[, c("variable_id", "var")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(var) ~ var,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-var)  %>% 
  dplyr::left_join(exposome_air_variable_info[,c("variable_id", "description")],
                   by = c("node" = "variable_id")) %>% 
  dplyr::mutate(true_name = case_when(!is.na(description) ~ description,
                                      TRUE ~ true_name)) %>%
  dplyr::select(-description) %>% 
  dplyr::left_join(exposome_outdoor_variable_info[, c("variable_id", "labels")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(labels) ~ labels,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-labels) %>% 
  dplyr::left_join(exposome_chemical_variable_info[, c("variable_id", "labels")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(labels) ~ labels,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-labels)

library(tidygraph)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))


pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

library(ggraph)

my_layout <- create_layout(temp_data, 
                           layout = 'linear')

my_layout$y[my_layout$class2 == "Internal_ome"] <- 5
my_layout$y[my_layout$class2 == "Phenotype"] <- 10

my_layout1 <-
  my_layout

my_layout1$x <-
  my_layout$y

my_layout1$y <-
  my_layout$x

my_layout1$y[my_layout1$class2 == "Exposome"] <-
  seq(from = 1, to = 100, length.out = sum(my_layout1$class2 == "Exposome"))

my_layout1$y[my_layout1$class2 == "Internal_ome"] <-
  seq(from = 1, to = 100, length.out = sum(my_layout1$class2 == "Internal_ome"))

my_layout1$y[my_layout1$class2 == "Phenotype"] <- 50

plot <-
  ggraph(my_layout1) +
  geom_edge_link(aes(color = cor),
                 show.legend = TRUE, alpha = 0.5) +
  geom_node_point(aes(fill = class1,
                      size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  scale_fill_manual(
    values = c(omics_color, phenotype_color)
  ) +
  scale_color_manual(
    values = c(omics_color, phenotype_color)
  ) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1,
      label = true_name,
      hjust = ifelse(class2 == "Internal_ome", 1, 0),
      size = 3,
      colour = class1
    ),
    size = 3,
    alpha = 1, 
    show.legend = FALSE
  ) +
  guides(edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
         fill = guide_legend(title = "Class", 
                             override.aes = list(size = 4, linetype = "blank")),
         size = guide_legend(title = "Degree", override.aes = list(linetype = 0))) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
    # legend.position = c(1,0), legend.justification = c(1,0)
  ) 

plot

# ggsave(plot, filename = "BMI_mediation_analysis.pdf", width = 7, height = 7)









