omics_color <- 
  c("Exposome_air" = ggsci::pal_d3()(10)[1],
    "Exposome_outdoor" = ggsci::pal_d3()(10)[5],
    "Exposome_chemical" = ggsci::pal_d3()(10)[2],
    "Proteome" = ggsci::pal_d3()(10)[4],
    "Serum_metabolome" = ggsci::pal_d3()(10)[3],
    "Urine_metabolome" = ggsci::pal_d3()(10)[6],
    "Transcriptome" = ggsci::pal_d3()(10)[7]
  )




phenotype_color = 
  c("Behavior" = RColorBrewer::brewer.pal(n = 3,name = "Accent")[1],
    "BMI" = RColorBrewer::brewer.pal(n = 3,name = "Accent")[2],
    "IQ" = RColorBrewer::brewer.pal(n = 3,name = "Accent")[3])
