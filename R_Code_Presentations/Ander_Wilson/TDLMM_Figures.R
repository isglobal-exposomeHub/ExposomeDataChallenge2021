rm(list=ls())
library(dlmtree)
library(ggplot2)
library(tidyverse)
library(viridis)

load("~/Projects/ISGlobalExposomeChallenge/TDLMM/Output/TDLMM_bmi_noself_50.Rdata")
sfit




#------------------------------
# Plot bayes factors and PIPs
#------------------------------


# main effect bayes factor

BFmain <- data.frame(exp=names(sfit$expBF), BF=sfit$expBF)
# change infinite BF to really high for plotting
BFmain$BF[which(BFmain$BF==Inf)] <- 10*max(BFmain$BF[which(BFmain$BF!=Inf)])
BFmain$selected <- c("No","Yes")[(BFmain$BF>10^.5)+1]

scale_fill_discrete <- function(...) scale_fill_manual(values=c("black","blue"))
scale_colour_discrete <- function(...) scale_colour_manual(values=c("black","blue"))

p_main_BF <- ggplot(BFmain, aes(x=reorder(exp,-BF), y=BF, colour=selected)) + 
  geom_point(size=4) + 
  theme_classic(base_size = 16) + 
  ylab("Bayes Factor") + xlab(NULL) +
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1), legend.position = "none") + 
  geom_hline(yintercept=10^.5, linetype=2, size=2) 
  

pdf("TDLMM/Figures/TDLMM_BF_main.pdf", height=5, width=10)
print(p_main_BF) 
dev.off()


# main effect PIPs

PIPmain <- data.frame(exp=names(sfit$expInc), PIP=sfit$expInc)
PIPmain$selected <- c("No","Yes")[(PIPmain$PIP>.5)+1]

p_main_PIP <- ggplot(PIPmain, aes(x=reorder(exp,-PIP), y=PIP, colour=selected)) + 
  geom_point(size=4) + 
  theme_classic(base_size = 16) + 
  ylab("Posterior inclusion probability (PIP)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1), legend.position = "none") + 
  geom_hline(yintercept=.5, linetype=2, size=2) +ylim(0,1) 


pdf("TDLMM/Figures/TDLMM_PIP_main.pdf", height=5, width=10)
print(p_main_PIP) 
dev.off()


# interaction  PIPs

PIPinteraction <- data.frame(exp=names(sfit$mixInc), PIP=sfit$mixInc)
PIPinteraction$selected <- c("No","Yes")[(PIPinteraction$PIP>.5)+1]

p_ineraction_PIP <- ggplot(PIPinteraction, aes(x=reorder(exp,-PIP), y=PIP, colour=selected)) + 
  geom_point(size=1) + 
  theme_classic(base_size = 16) + 
  ylab("Posterior inclusion probability (PIP)") + xlab("Pairwise interaction") +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) + 
  geom_hline(yintercept=.5, linetype=2) + ylim(0,1) + 
  theme(legend.position = "none", 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


pdf("TDLMM/Figures/TDLMM_PIP_interaction.pdf", height=5, width=10)
print(p_ineraction_PIP) 
dev.off()


PIPinteraction <- data.frame(exp=names(sfit$mixInc), PIP=sfit$mixInc)
PIPinteraction$selected <- c("No","Yes")[(PIPinteraction$PIP>.5)+1]
PIPinteraction_top1pct <- PIPinteraction %>% filter(PIP>quantile(PIPinteraction$PIP,.99))

p_ineraction_PIP_top1pct <- ggplot(PIPinteraction_top1pct, aes(x=reorder(exp,-PIP), y=PIP, colour=selected)) + 
  geom_point(size=4) + 
  theme_classic(base_size = 14) + 
  ylab("PIP") + xlab(NULL) +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1), legend.position = "none") + 
  geom_hline(yintercept=.5, linetype=2) 


pdf("TDLMM/Figures/TDLMM_PIP_interaction_top1pct.pdf", height=4, width=3)
print(p_ineraction_PIP_top1pct) 
dev.off()




#------------------------------
# Main effect exposure-response
#------------------------------

PIPmain_selected <- PIPmain %>% filter(PIP>.5)

main_effects <- data.frame()
for(i in PIPmain_selected$exp){
  main_effects <- rbind(main_effects,
                        data.frame(exposure=i,
                                   Est = sfit$DLM[[i]]$marg.matfit,
                    CIMin = sfit$DLM[[i]]$marg.cilower,
                    CIMax = sfit$DLM[[i]]$marg.ciupper,
                    X = c("Pre","Post"),
                    Xorder = c(1,2),
                    sig = sfit$DLM[[i]]$marg.cilower>0 | sfit$DLM[[i]]$marg.ciupper<0) 
  )
}


p_main_effects <- ggplot(main_effects, aes(x = reorder(X,Xorder), y = Est, ymin = CIMin, ymax = CIMax, fill=sig, color=sig)) +
  geom_hline(yintercept = 0, linetype=2, size=.5) +
  geom_errorbar(size=1.5, width=.5) +
  geom_point(size=4) +
  theme_bw() +
  facet_grid(.~exposure) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.text = element_text(size=12)) + 
  ylab("Parameter estimate") + xlab(NULL)



pdf("TDLMM/Figures/TDLMM_main_effects.pdf", height=5, width=10)
print(p_main_effects)
dev.off()





#------------------------------
# Interaction effects 
#------------------------------

PIPinteractionSelected <- PIPinteraction %>% filter(selected=="Yes")


plotDatTemp <- data.frame(X = rep(c("Pre","Post"),2),
                      Y = rep(c("Pre","Post"), each=2),
                      Xorder = c(1,2),
                      Yorder = c(1,2))
plotDat <- data.frame()

for(i in PIPinteractionSelected$exp){
  plotDat <- rbind(plotDat,
                   cbind.data.frame(plotDatTemp,
                                    combination=i,
                              Effect = c(sfit$MIX[[i]]$matfit),
                              CW = c(sfit$MIX[[i]]$cw.plot)))  
}


p_interactions <- ggplot(plotDat, aes(x = reorder(X,Xorder), y = reorder(Y,Yorder), 
                                      z = Effect, fill = Effect)) +
  geom_tile() +
  scale_fill_viridis() + 
  geom_point(data = plotDat[which(plotDat$CW != 0),],
             aes(x = reorder(X,Xorder), y = reorder(Y,Yorder), size = CW),
             color = "red", show.legend = FALSE) +
  theme_bw(base_size=16) +
  coord_equal()  +
  facet_wrap(~combination, ncol=3) +
  xlab("Period") + ylab("Period")



pdf("TDLMM/Figures/TDLMM_interaction_effects.pdf", height=5, width=10)
print(p_interactions)
dev.off()

