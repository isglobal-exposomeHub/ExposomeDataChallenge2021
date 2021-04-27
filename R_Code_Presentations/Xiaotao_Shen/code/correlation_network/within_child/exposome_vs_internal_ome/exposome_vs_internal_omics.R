##avoid source
no_function()

sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())
source("code/tools.R")

###load data
load("data_analysis/correlation_network/within_child/exposome_air_vs_transcriptome/cor_value")
exposome_air_vs_transcriptome_cor = cor_value

load("data_analysis/correlation_network/within_child/exposome_air_vs_proteome/cor_value")
exposome_air_vs_proteome_cor = cor_value

load("data_analysis/correlation_network/within_child/exposome_air_vs_serum_metabolome/cor_value")
exposome_air_vs_serum_metabolome_cor = cor_value

load("data_analysis/correlation_network/within_child/exposome_air_vs_urine_metabolome/cor_value")
exposome_air_vs_urine_metabolome_cor = cor_value

load("data_analysis/correlation_network/within_child/exposome_outdoor_vs_transcriptome/cor_value")
exposome_outdoor_vs_transcriptome_cor = cor_value

load("data_analysis/correlation_network/within_child/exposome_outdoor_vs_proteome/cor_value")
exposome_outdoor_vs_proteome_cor = cor_value

load("data_analysis/correlation_network/within_child/exposome_outdoor_vs_serum_metabolome/cor_value")
exposome_outdoor_vs_serum_metabolome_cor = cor_value

load("data_analysis/correlation_network/within_child/exposome_outdoor_vs_urine_metabolome/cor_value")
exposome_outdoor_vs_urine_metabolome_cor = cor_value

load("data_analysis/correlation_network/within_child/exposome_chemical_vs_transcriptome/cor_value")
exposome_chemical_vs_transcriptome_cor = cor_value

load("data_analysis/correlation_network/within_child/exposome_chemical_vs_proteome/cor_value")
exposome_chemical_vs_proteome_cor = cor_value

load("data_analysis/correlation_network/within_child/exposome_chemical_vs_serum_metabolome/cor_value")
exposome_chemical_vs_serum_metabolome_cor = cor_value

load("data_analysis/correlation_network/within_child/exposome_chemical_vs_urine_metabolome/cor_value")
exposome_chemical_vs_urine_metabolome_cor = cor_value

exposome_internal_omics_cor = 
  rbind(exposome_air_vs_transcriptome_cor,
        exposome_air_vs_proteome_cor,
        exposome_air_vs_serum_metabolome_cor,
        exposome_air_vs_urine_metabolome_cor,
        exposome_outdoor_vs_transcriptome_cor,
        exposome_outdoor_vs_proteome_cor,
        exposome_outdoor_vs_serum_metabolome_cor,
        exposome_outdoor_vs_urine_metabolome_cor,
        exposome_chemical_vs_transcriptome_cor,
        exposome_chemical_vs_proteome_cor,
        exposome_chemical_vs_serum_metabolome_cor,
        exposome_chemical_vs_urine_metabolome_cor
        )

dim(exposome_internal_omics_cor)

###load the variable information of all omics data
load("data_analysis/transcriptome_data_analysiss/data_preparation/variable_info")
transcriptome_variable_info = variable_info

load("data_analysis/proteome_data_analysiss/data_preparation/variable_info")
proteome_variable_info = variable_info

load("data_analysis/serum_metabolome_data_analysiss/data_preparation/variable_info")
serum_metabolome_variable_info = variable_info

load("data_analysis/urine_metabolome_data_analysiss/data_preparation/variable_info")
urine_metabolome_variable_info = variable_info

load("data_analysis/exposome_air_data_analysis/data_preparation/variable_info")
exposome_air_variable_info = variable_info

load("data_analysis/exposome_outdoor_data_analysis/data_preparation/variable_info")
exposome_outdoor_variable_info = variable_info

load("data_analysis/exposome_chemical_data_analysis/data_preparation/variable_info")
exposome_chemical_variable_info = variable_info


sxtTools::setwd_project()
dir.create("data_analysis/correlation_network/within_child/exposome_vs_internal_omics", showWarnings = FALSE)
setwd("data_analysis/correlation_network/within_child/exposome_vs_internal_omics")

# save(exposome_internal_omics_cor, file = "exposome_internal_omics_cor")

library(igraph)
library(ggraph)
library(tidygraph)


edge_data <-
  exposome_internal_omics_cor %>%
  dplyr::rename(from = from,
                to = to,
                Correlation = cor) %>%
  dplyr::mutate(p.adjust = -log(p.adjust, 10))



node_data <-
  exposome_internal_omics_cor %>%
  dplyr::select(from, to) %>%
  tidyr::pivot_longer(cols = c(from, to),
                      names_to = "class",
                      values_to = "node") %>%
  dplyr::mutate(
    class1 = case_when(
      node %in% exposome_air_variable_info$variable_id ~ "Exposome_air",
      node %in% exposome_outdoor_variable_info$variable_id ~ "Exposome_outdoor",
      node %in% exposome_chemical_variable_info$variable_id ~ "Exposome_chemical",
      node %in% transcriptome_variable_info$variable_id ~ "Transcriptome",
      node %in% urine_metabolome_variable_info$variable_id ~ "Urine_metabolome",
      node %in% serum_metabolome_variable_info$variable_id ~ "Serum_metabolome",
      node %in% proteome_variable_info$variable_id ~ "Proteome"
    )
  ) %>%
  dplyr::select(node, class1) %>%
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::arrange(Class) %>%
  dplyr::mutate(Class = factor(
    Class,
    levels = c(
      "Exposome_air",
      "Exposome_outdoor",
      "Exposome_chemical",
      "Transcriptome",
      "Proteome",
      "Urine_metabolome",
      "Serum_metabolome"
    )
  )) %>% 
  dplyr::left_join(exposome_air_variable_info[,c("variable_id", "description")],
                   by = c("node" = "variable_id")) %>% 
  dplyr::mutate(true_name = case_when(!is.na(description) ~ description,
                                      TRUE ~ node)) %>%
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

node_data %>% 
  dplyr::filter(Class == "Exposome_air")

node_data %>% 
  dplyr::filter(Class == "Exposome_chemical")

node_data %>% 
  dplyr::filter(Class == "Exposome_outdoor")

node_data %>% 
  dplyr::filter(Class == "Proteome")

node_data %>% 
  dplyr::filter(Class == "Serum_metabolome")

node_data %>% 
  dplyr::filter(Class == "Transcriptome")

node_data %>% 
  dplyr::filter(Class == "Urine_metabolome")


####output node data and edge data
library(openxlsx)
wb = createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Arial Narrow")
addWorksheet(wb, sheetName = "Node information", gridLines = TRUE)
addWorksheet(wb, sheetName = "Edge information", gridLines = TRUE)
freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) 
freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE) 
writeDataTable(wb, sheet = 1, x = node_data,
               colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 2, x = edge_data %>% dplyr::select(from, to, everything()),
               colNames = TRUE, rowNames = FALSE)

saveWorkbook(wb, "exposome_internal_ome_network.xlsx", overwrite = TRUE)


alpha_value <-
  c(
    "Exposome_air" = 1,
    "Exposome_outdoor" = 1,
    "Exposome_chemical" = 1,
    "Proteome" = 1,
    "Serum_metabolome" = 1,
    "Urine_metabolome" = 1,
    "Transcriptome" = 1
  )


total_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

plot <-
  ggraph(total_graph,
         layout = 'stress') +
  geom_edge_link(aes(color = Correlation),
                 show.legend = TRUE,
                 alpha = 0.5) +
  geom_node_point(
    aes(fill = Class,
        size = Degree,
        alpha = Class),
    shape = 21,
    # alpha = 0.6,
    show.legend = TRUE
  ) +
  scale_fill_manual(values = omics_color) +
  scale_color_manual(values = omics_color) +
  scale_alpha_manual(values = alpha_value) +
  # geom_node_text(
  #   aes(
  #     x = x * 1.05,
  #     y = y * 1.05,
  #     label = true_name,
  #     hjust = 'outward',
  #     angle = -((-node_angle(x, y) + 90) %% 180) + 90,
  #     size = 3,
  #     colour = Class
  #   ),
  #   size = 3,
#   alpha = 1,
#   show.legend = FALSE
# ) +
guides(
  edge_width = guide_legend(title = "-log10(FDR adjusted P value)",
                            override.aes = list(shape = NA)),
  edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
  fill = guide_legend(
    title = "Class",
    override.aes = list(size = 7, linetype = "blank")
  ),
  size = guide_legend(title = "Degree", override.aes = list(linetype = 0))
) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(0.5, 5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

plot


# save(total_graph, file = "total_graph")
load("total_graph")

ggsave(
  plot,
  filename = "inter_omics_correlation_network.pdf",
  width = 10,
  height = 7,
  bg = "transparent"
)










sum(edge_data$from %in% node_data$node[node_data$Class == "Exposome_air"])
sum(edge_data$from %in% node_data$node[node_data$Class == "Exposome_chemical"])
sum(edge_data$from %in% node_data$node[node_data$Class == "Exposome_outdoor"])


sum(edge_data$to %in% node_data$node[node_data$Class == "Transcriptome"])
sum(edge_data$to %in% node_data$node[node_data$Class == "Proteome"])
sum(edge_data$to %in% node_data$node[node_data$Class == "Serum_metabolome"])
sum(edge_data$to %in% node_data$node[node_data$Class == "Urine_metabolome"])

###edge distributation
temp_data <-
  exposome_internal_omics_cor %>%
  dplyr::select(from, to) %>%
  dplyr::mutate(
    var1_class = case_when(
      from %in% exposome_air_variable_info$variable_id ~ "Exposome_air",
      from %in% exposome_outdoor_variable_info$variable_id ~ "Exposome_outdoor",
      from %in% exposome_chemical_variable_info$variable_id ~ "Exposome_chemical",
      from %in% transcriptome_variable_info$variable_id ~ "Transcriptome",
      from %in% urine_metabolome_variable_info$variable_id ~ "Urine_metabolome",
      from %in% serum_metabolome_variable_info$variable_id ~ "Serum_metabolome",
      from %in% proteome_variable_info$variable_id ~ "Proteome"
    )
  ) %>%
  dplyr::mutate(
    var2_class =
      case_when(
        to %in% exposome_air_variable_info$variable_id ~ "Exposome_air",
        to %in% exposome_outdoor_variable_info$variable_id ~ "Exposome_outdoor",
        to %in% exposome_chemical_variable_info$variable_id ~ "Exposome_chemical",
        to %in% transcriptome_variable_info$variable_id ~ "Transcriptome",
        to %in% urine_metabolome_variable_info$variable_id ~ "Urine_metabolome",
        to %in% serum_metabolome_variable_info$variable_id ~ "Serum_metabolome",
        to %in% proteome_variable_info$variable_id ~ "Proteome"
      )
  )

library(ggalluvial)

temp_data1 <- temp_data %>%
  dplyr::mutate(id = paste(from, to, var1_class, var2_class, sep = "_"))

temp_data2 <-
  temp_data1 %>%
  dplyr::select(id, var1_class, var2_class) %>%
  tidyr::pivot_longer(cols = -id,
                      names_to = "class",
                      values_to = "data")  %>%
  dplyr::mutate(freq = 1) %>%
  dplyr::mutate(data = factor(
    data,
    levels = c(
      "Exposome_air",
      "Exposome_outdoor",
      "Exposome_chemical",
      "Transcriptome",
      "Proteome",
      "Urine_metabolome",
      "Serum_metabolome"
    )
  ))

plot <-
  ggplot(temp_data2,
         aes(
           x = class,
           y = freq,
           stratum = data,
           alluvium = id,
           fill = data,
           label = data
         )) +
  scale_x_discrete(expand = c(.1, .1)) +
  ggalluvial::geom_stratum(alpha = 1,
                           color = "black",
                           show.legend = FALSE) +
  ggalluvial::geom_flow(show.legend = FALSE) +
  geom_text(stat = "stratum", size = 5, color = "white") +
  labs(x = "", y = "") +
  scale_fill_manual(values = omics_color) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(
  plot,
  file = file.path("edge_information.pdf"),
  width = 8,
  height = 7,
  bg = "transparent"
)


####node distributation
plot <-
  node_data %>%
  dplyr::group_by(Class) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Class = factor(
    Class,
    levels = c(
      "Exposome_air",
      "Exposome_outdoor",
      "Exposome_chemical",
      "Transcriptome",
      "Proteome",
      "Urine_metabolome",
      "Serum_metabolome"
    ) %>%
      rev()
  )) %>%
  ggplot(aes(x = 2, y = Class)) +
  geom_point(shape = 21,
             aes(fill = Class, size = n),
             show.legend = FALSE) +
  theme_void() +
  scale_size_continuous(range = c(5, 30)) +
  scale_fill_manual(values = omics_color) +
  geom_text(aes(x = 2, y = Class, label = n), color = "white")

plot

# ggsave(
#   plot,
#   file = file.path("node_information.pdf"),
#   width = 5,
#   height = 7,
#   bg = "transparent"
# )

total_graph


