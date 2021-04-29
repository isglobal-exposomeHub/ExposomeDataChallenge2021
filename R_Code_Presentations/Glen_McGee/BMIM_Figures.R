library(bsmim2)
library(tidyverse)
library(patchwork)

# load data
load("DataPrep/index_exposure_list.rda")

# load results from full model
load("~/Projects/ISGlobalExposomeChallenge/BMIM/Output/bmim_bmi_base_select.Rdata")
ls()

#------------------------------
# Plot PIPs
#------------------------------

# get PIPs and weight estimates
ThetaTable <- bind_rows(summarize_thetas(fit))

# add index and component labels
ThetaTable$Index <- index_identifiers$index_name
ThetaTable$Component <- index_identifiers$labelsshort

# make marginal PIPs
ThetaTable$PIPmarginal <- ThetaTable$PIP * ThetaTable$PIP_RHO

# plot component PIPs
Index_PIPs <- ThetaTable %>% 
  select(PIP_RHO,Index) %>% 
  distinct()

p_pip_index <- ggplot(Index_PIPs) + 
  geom_point(aes(x=reorder(Index,PIP_RHO), y=PIP_RHO)) + 
  ylim(0,1) + coord_flip() +
  ylab("Component PIP") + xlab(NULL) +
  theme_light()

pdf("BMIM/Figures/BMIM_pip_index.pdf", height=5, width=5)
print(p_pip_index)
dev.off()

# plot component PIPs for the selected indices
Component_PIPs <- ThetaTable %>% 
  filter(PIP_RHO>0.5)  %>% 
  mutate(Index2 = str_replace(Index, "_", "\n"))

p_pip_components <- ggplot(Component_PIPs) + 
  geom_point(aes(x=reorder(Component,PIPmarginal), y=PIPmarginal)) + 
  ylim(0,1) +  xlab(NULL) + ylab("Marginal component PIP") + 
  facet_grid(Index2~., scales="free_y", space="free_y") + 
  coord_flip() +
  theme_light() + 
  theme(strip.text.y = element_text(angle=0, colour="black"),
        strip.background = element_blank())


pdf("BMIM/Figures/BMIM_pip_components.pdf", height=5, width=5)
print(p_pip_components)
dev.off()
 


# plot component weights for the selected indices

p_pip_weights <- ggplot(Component_PIPs, aes(x=reorder(Component,PIPmarginal), y=est, ymin=lci, ymax=uci)) + 
  geom_point() + 
  geom_errorbar(width=.5) + 
  ylim(-1,1) +  xlab(NULL) + ylab("Component weight") + 
  facet_grid(Index2~., scales="free_y", space="free_y") + 
  coord_flip() +
  theme_light() + 
  theme(strip.text.y = element_text(angle=0, colour="black"),
        strip.background = element_blank())


pdf("BMIM/Figures/BMIM_pip_weights.pdf", height=5, width=6)
print(p_pip_weights)
dev.off()




#------------------------------
# Plot component wise exposure response functions
#------------------------------

# find selected curves
selected_components <- which(names(index_exposure_list)%in%Component_PIPs$Index)
selected_component_names <- names(index_exposure_list)[selected_components]
selected_component_names <- str_replace(selected_component_names,"_"," ")
# plot selected curves
plot_component_ER <- plot_univar_hnew_indexwise2(pred_ind,1)
plot_component_ER <- plot_component_ER[selected_components]

# add title
for(i in 1:length(plot_component_ER)){
  plot_component_ER[[i]] <- plot_component_ER[[i]] + ggtitle(selected_component_names[i])
}


pdf("BMIM/Figures/BMIM_index_exposure_response.pdf", height=5, width=8)
(plot_component_ER[[1]] + plot_component_ER[[2]]) /
  (plot_component_ER[[3]] + plot_component_ER[[4]]) 
dev.off()

#------------------------------
# Look at interaction plots
#------------------------------

pred_inter_selected <- pred_inter %>% 
  filter(var1%in%paste("Index",selected_components) & var2%in%paste("Index",selected_components))

selected_component_names <- str_replace(selected_component_names," ","\n")


pred_inter_selected <- left_join(pred_inter_selected,data.frame(var1=paste("Index",selected_components), name1=selected_component_names, order1=selected_components))
pred_inter_selected <- left_join(pred_inter_selected,data.frame(var2=paste("Index",selected_components), name2=selected_component_names, order2=selected_components))


pp_inter <- ggplot(pred_inter_selected, aes(grid, est)) + 
  geom_smooth(aes(col = as.factor(quantile)), stat = "identity",fill="white") + 
  facet_grid(reorder(name2,order2) ~ reorder(name1,order1), scales = "free_x") +
  labs(x="",y="",col="Quantile")+
  theme_light() + 
  theme(strip.text = element_text(colour = 'black'),
        strip.text.y = element_text(angle=0),
        strip.background = element_blank(), 
        legend.position = "bottom")


pdf("BMIM/Figures/BMIM_interactions.pdf", height=5, width=6)
print(pp_inter)
dev.off()

