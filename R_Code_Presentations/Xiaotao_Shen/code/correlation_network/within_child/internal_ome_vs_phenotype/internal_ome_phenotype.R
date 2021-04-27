##avoid source
no_function()

sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())
source("code/tools.R")

##load data
###load the variable information of all omics data
load("data_analysis/transcriptome_data_analysiss/data_preparation/variable_info")
transcriptome_variable_info = variable_info

load("data_analysis/proteome_data_analysiss/data_preparation/variable_info")
proteome_variable_info = variable_info

load("data_analysis/serum_metabolome_data_analysiss/data_preparation/variable_info")
serum_metabolome_variable_info = variable_info

load("data_analysis/urine_metabolome_data_analysiss/data_preparation/variable_info")
urine_metabolome_variable_info = variable_info

dir.create("data_analysis/correlation_network/within_child/internal_ome_phenotype")
setwd("data_analysis/correlation_network/within_child/internal_ome_phenotype")

load("../transcriptome_vs_phenotype/transcriptome_phenotype_glm")
load("../proteome_vs_phenotype/proteome_phenotype_glm")
load("../serum_metabolome_vs_phenotype/serum_metabolome_phenotype_glm")
load("../urine_metabolome_vs_phenotype/urine_metabolome_phenotype_glm")

transcriptome_phenotype_glm = 
  transcriptome_phenotype_glm %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  dplyr::mutate(class = "Transcriptome")

proteome_phenotype_glm = 
  proteome_phenotype_glm %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  dplyr::mutate(class = "Proteome")

serum_metabolome_phenotype_glm = 
  serum_metabolome_phenotype_glm %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  dplyr::mutate(class = "Serum_metabolome")

urine_metabolome_phenotype_glm = 
  urine_metabolome_phenotype_glm %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  dplyr::mutate(class = "Urine_metabolome")

internal_ome_phenotype_glm =
  rbind(
    transcriptome_phenotype_glm,
    proteome_phenotype_glm,
    serum_metabolome_phenotype_glm,
    urine_metabolome_phenotype_glm
  ) %>%
  tidyr::pivot_longer(
    cols = -c(variable_id, class),
    names_to = "phenotype",
    values_to = "p.adjust"
  )

save(internal_ome_phenotype_glm, file = "internal_ome_phenotype_glm")

load("internal_ome_phenotype_glm")

library(openxlsx)

internal_ome_phenotype_glm

temp_data = 
  internal_ome_phenotype_glm %>% 
  dplyr::mutate(significant = case_when(
    p.adjust < 0.05 ~ "YES",
    TRUE ~ "NO"
  )) %>% 
  dplyr::group_by(phenotype, class, significant) %>% 
  dplyr::summarise(n = n())

internal_ome_phenotype_glm %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::distinct(variable_id, .keep_all = TRUE) %>% 
  pull(class) %>% 
  table()



#####
library(ggalluvial)
temp_data = 
  internal_ome_phenotype_glm %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::mutate(variable_id = paste(variable_id, phenotype, sep = "_")) %>% 
  dplyr::mutate(freq = 1) 

temp_data$phenotype[temp_data$phenotype == "Body.mass.index.z.score"] = "BMI"
temp_data$phenotype[temp_data$phenotype == "Intelligence.quotient"] = "IQ"

plot = 
  temp_data %>%
  ggplot(aes(y = freq, axis1 = phenotype, axis2 = class)) +
  geom_alluvium(aes(fill = phenotype), width = 1 / 12, show.legend = FALSE) +
  geom_stratum(width = 1 / 12,
               # fill = "grey",
               aes(fill = phenotype),
               color = "black",
               show.legend = FALSE) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Phenotype", "Exposome"),
                   expand = c(.05, .05)) +
  scale_fill_manual(values = phenotype_color) +
  coord_flip() +
  theme_void()

plot

ggsave(plot, filename = "alluvium.pdf", width = 10, height = 7)

internal_ome_phenotype_glm %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dim()

internal_ome_phenotype_glm %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::distinct(.keep_all = TRUE) %>% 
  pull(phenotype) %>% 
  table()

internal_ome_phenotype_glm %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  pull(class) %>% 
  table()

library(ggvenn)

temp = list("IQ" = 
              internal_ome_phenotype_glm %>% 
              dplyr::filter(p.adjust < 0.05) %>% 
              dplyr::filter(phenotype == "Intelligence.quotient") %>% 
              dplyr::pull(variable_id),
            "BMI" = 
              internal_ome_phenotype_glm %>% 
              dplyr::filter(p.adjust < 0.05) %>% 
              dplyr::filter(phenotype == "Body.mass.index.z.score") %>% 
              dplyr::pull(variable_id),
            "Behavior" = 
              internal_ome_phenotype_glm %>% 
              dplyr::filter(p.adjust < 0.05) %>% 
              dplyr::filter(phenotype == "Behavior") %>% 
              dplyr::pull(variable_id))

plot = 
  ggvenn(
    data = temp,
    fill_color = unname(phenotype_color[c("IQ", "BMI", 'Behavior')]),
    stroke_color = "white", set_name_size = 10, text_size = 6, fill_alpha = 0.9, 
    text_color = "white"
  )
plot
# ggsave(plot, filename = "venn.pdf", width = 7, height = 7)



####network plot to show the correlations.
edge_data = 
  internal_ome_phenotype_glm %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::rename(from = variable_id, to = phenotype) %>% 
  dplyr::select(from, to) %>% 
  dplyr::mutate(Class = to) %>% 
  dplyr::mutate(Class = case_when(
    Class %in% "Behavior" ~ "Behavior",
    Class %in% "Intelligence.quotient" ~ "IQ",
    Class %in% "Body.mass.index.z.score" ~ "BMI"
  ))

node_data <-
  edge_data %>%
  dplyr::select(from, to) %>%
  tidyr::pivot_longer(cols = c(from, to),
                      names_to = "class",
                      values_to = "node") %>%
  dplyr::mutate(
    class1 = case_when(
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
  dplyr::select(node, class1) %>%
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::arrange(Class) %>%
  # dplyr::mutate(Class = factor(
  #   Class,
  #   levels = c(
  #     "Behavior",
  #     "BMI",
  #     "IQ",
  #     "Transcriptome",
  #     "Proteome",
  #     "Serum_metabolome"
  #   )
  # )) %>%
dplyr::left_join(transcriptome_variable_info[,c("variable_id", "description")],
                 by = c("node" = "variable_id")) %>% 
  dplyr::mutate(true_name = case_when(!is.na(description) ~ description,
                                      TRUE ~ node)) %>%
  dplyr::select(-description) %>% 
  dplyr::left_join(proteome_variable_info[, c("variable_id", "labels")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(labels) ~ labels,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-labels) %>% 
  dplyr::left_join(serum_metabolome_variable_info[, c("variable_id", "labels")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(labels) ~ labels,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-labels) %>% 
  dplyr::left_join(transcriptome_variable_info[, c("variable_id", "GeneSymbol_Affy")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(GeneSymbol_Affy) ~ GeneSymbol_Affy,
                                      TRUE ~ true_name)) %>% 
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
  dplyr::select(-var)  %>% 
  dplyr::left_join(urine_metabolome_variable_info[, c("variable_id", "var")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(var) ~ var,
                                      TRUE ~ true_name)) %>% 
  dplyr::select(-var) 



####output node data and edge data
library(openxlsx)
# wb = createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Arial Narrow")
# addWorksheet(wb, sheetName = "Node information", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Edge information", gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE)
# writeDataTable(wb, sheet = 1, x = node_data,
#                colNames = TRUE, rowNames = FALSE)
# writeDataTable(wb, sheet = 2, x = edge_data %>% dplyr::select(from, to, everything()),
#                colNames = TRUE, rowNames = FALSE)
# 
# saveWorkbook(wb, "internal_ome_internal_ome_network.xlsx", overwrite = TRUE)

total_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

plot <-
  ggraph(total_graph,
         layout = 'linear', 
         circular = TRUE) +
  geom_edge_arc(show.legend = FALSE,
                alpha = 1,
                aes(color = Class)) +
  geom_node_point(
    aes(fill = Class,
        size = Class),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  scale_fill_manual(values = c(omics_color, phenotype_color)) +
  scale_color_manual(values = c(omics_color, phenotype_color)) +
  scale_edge_color_manual(values = phenotype_color) +
  scale_size_manual(values = c(
    "Behavior" = 5,
    "BMI" = 5,
    "IQ" = 5,
    "Transcriptome" = 4,
    "Proteome" = 4,
    "Serum_metabolome" = 4,
    "Transcriptome" = 4,
    "Proteome" = 4,
    "Urine_metabolome" = 4,
    "Serum_metabolome" = 4
  )) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1.05,
      label = true_name,
      hjust = 'outward',
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      size = 3,
      colour = Class
    ),
    size = 3,
    alpha = 1,
    show.legend = FALSE
  ) +
  guides(
    edge_width = guide_legend(title = "-log10(FDR adjusted P value)",
                              override.aes = list(shape = NA)),
    edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
    fill = guide_legend(
      title = "Class",
      override.aes = list(size = 7, linetype = "blank")
    )
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

plot

save(total_graph, file = "total_graph")
load("total_graph")

# ggsave(
#   plot,
#   filename = "internal_ome_phenotype_correlation_network.pdf",
#   width = 8.5,
#   height = 7,
#   bg = "transparent"
# )








