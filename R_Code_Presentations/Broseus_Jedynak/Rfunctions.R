################################################################################

packNames <- c('MASS', 'dplyr', 'magrittr', 'FactoMineR',
               'ggplot2', 'ordinalClust', 'ordinalForest')
packNames <- setdiff(packNames, installed.packages())

if(length(packNames)>0) for( p in packNames ) install.packages(p,dependencies = TRUE)

suppressPackageStartupMessages({
  
  library(magrittr)
  
  library(ggplot2)
  
  library(FactoMineR)
  
}
  
)

################################################################################
# Display

myTheme <- theme_minimal() + 
  theme(panel.background = element_rect(fill = "snow1", 
                                        colour = "snow1", 
                                        size =1,  
                                        linetype = "solid")) 

################################################################################
# For ordinal covariates: transform tertiles into ordered classes: 1<2<3

tertile2class <- function(string){
  
  string <- stringr::str_split(string = string, pattern = ",")[[1]]
  x <- grep(pattern = "Inf", x = string)
  
  if(length(x)>0) return( 3 )
  else if(string[1]=="(0") return( 1 )
  else return( 2 )
  
}

################################################################################
# Plot MCA of all categorical covariates in the input matrix
# return MCA plot as ggplot object

Lu_plotMCA <- function(catMatrix, ind.sup = NULL, quali.sup = NULL, title = "MCA plot"){
  
  # columns must be factors
  catMatrix <- apply(catMatrix,2, as.factor)
  cats <- apply(catMatrix, 2, function(x) nlevels(as.factor(x)))
  # MCA
  mca <- FactoMineR::MCA(catMatrix, 
                         ind.sup = ind.sup, 
                         quali.sup, graph = FALSE)
  
  # data frames for ggplot
  mca_vars_df <- data.frame(mca$var$coord, Variable = rep(names(cats), cats))
  mca_obs_df <- data.frame(mca$ind$coord)
  
  # plot of variable categories
  p <- ggplot(data = mca_vars_df, aes(x = Dim.1, y = Dim.2, label = rownames(mca_vars_df))) + 
       geom_density2d(colour = "gray80") +
       geom_hline(yintercept = 0, colour = "gray70") + 
       geom_vline(xintercept = 0, colour = "gray70") + 
       geom_text(aes(colour = Variable)) + 
       ggtitle(title)
  
  return( p )
}

################################################################################
# Ordinal logistic regression
# uses MASS::polr

Lu_ordinalLogistic <- function(data, 
                               mainEffect,
                               isOrdered = F,
                               refLevel = NULL,
                               confounders, 
                               outcome, 
                               level = 0.95,
                               verbose = FALSE){
  
  if(isOrdered) data[,mainEffect] <- ordered(data[,mainEffect])
  
  data[,outcome] <- ordered(data[,outcome])
  
  formula <- stats::as.formula(paste(outcome, "~",
                                   mainEffect, "+",
                                   paste(confounders, collapse = "+")))
  
  if(!is.null(refLevel)) data[,mainEffect] <- relevel(factor(data[,mainEffect]), ref = refLevel)

  res <- MASS::polr(formula = formula, data= data, Hess=TRUE)
  
  if(verbose) print(summary(res))
  
  suppressMessages( res <- confint(res, level = level) )
  
  if(verbose) print(res)
  
  sfit <- matrix(res[grep(pattern = mainEffect, x = rownames(res)),], ncol = 2)

  sfit <- round(sfit,3)
  sfit <- cbind.data.frame(variable_name = mainEffect, sfit)
  
  colnames(sfit) <- c("variable_name","conf_inf", "conf_sup")
  
  if(isOrdered) sfit <- cbind.data.frame(sfit, polynomial = c("L", "Q"))
  else sfit <- cbind.data.frame(sfit, levels = rownames(res)[grep(pattern = mainEffect, x = rownames(res))])
  
  return(sfit)
  
}

################################################################################
# Ordinal random forest
# return  a data.frame variable_name vs var_importance

Lu_ordRF.perform <- function(data, 
                             outcome.name, 
                             perffunction = "equal",
                             classweights = NA,
                             classimp = NA){
  
  data[, outcome.name] <- as.factor(data[, outcome.name])
  
  ordforres <- ordinalForest::ordfor(depvar = outcome.name, 
                                     data = data, 
                                     nsets=1000, 
                                     ntreeperdiv=100,
                                     ntreefinal=5000, 
                                     perffunction = perffunction,
                                     classimp = classimp, classweights = NA)
  
  # Study variable importance values
  x <- sort(ordforres$varimp, decreasing=TRUE)
  
  res.ordRF <- data.frame(variable_name = names(x), var_importance = x) %>% 
                    dplyr::arrange(desc(var_importance))
  
  res.ordRF$variable_name <- ordered(res.ordRF$variable_name, levels = res.ordRF$variable_name)
  
  return( res.ordRF )

}

Lu_ordRF.plot <- function(res.ordRF, 
                          selectedVar.names = NULL, 
                          title = "Ordinal Random Forest"){
  
  if(!is.null(selectedVar.names)){
    selectedVar.names <- as.vector(selectedVar.names)
    res.ordRF <- res.ordRF  %>% 
                 dplyr::filter(variable_name %in% selectedVar.names) 
  } 
  
  res.ordRF$variable_name <- ordered(res.ordRF$variable_name, 
                                     levels = res.ordRF$variable_name)
  
  p <- res.ordRF %>% ggplot(aes(x = variable_name, y = var_importance)) + 
    geom_point(color="black", size=3, alpha=0.6) +
    geom_segment(aes(x=variable_name, xend=variable_name, 
                     y=0, yend=var_importance), color="black") +
    ggtitle(title) +
    xlab("Predictor") + ylab("Variable importance") +
    coord_flip() + myTheme
  
  return( p )

}

################################################################################
# BOS classification
# import ordinalClust::bosclust
# initialized using kmeans

Lu_bosclust <- function(M,  #quantile matrix
                        kr, #number of clusters
                        m,  #input data is ordinal data with m levels
                        nbindmini = 2, 
                        nbSEM = 100, nbSEMburn = 90){
  
  init <- "kmeans"
  M <- as.matrix(M)
  M <- apply(M, 2, as.numeric)
  
  object <- ordinalClust::bosclust(x = as.matrix(M),
                                   kr = kr, 
                                   m = m, 
                                   nbindmini = nbindmini, 
                                   init = init, 
                                   nbSEM = nbSEM, 
                                   nbSEMburn = nbSEMburn)
  
  bosclusters <- apply(object@V,1, function(l) which(l==1))
  
  return(  bosclusters )
  
}

