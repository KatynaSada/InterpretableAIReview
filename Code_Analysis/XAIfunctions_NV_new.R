library(CVST) # Kernel Ridge regression
library(DiagrammeRsvg) # save plot niceTree
library(rsvg) # save plot niceTree

###########################################
## Functions for XAI
###########################################

source("Code_Analysis/BOSO.r")

### Optimizing trees
findsplitExp <- function(C, Expression, minimum = 1, weights=NULL, verbose = F) {
  # C: IC50
  # X: mutation matrix
  if (verbose) print("Split!")
  if (!is.null(weights)) {
    C <- C[weights == 1,]
    Expression <- Expression[weights == 1,]
  }
  if (nrow(C) <= minimum) return(NULL)
  
  Gene <- which.min(apply(Expression, 2, getsumic50v3, C)) # Find optimal split
  names(Gene) <- colnames(Expression[,Gene])
  Output <- getSplit(Gene, Expression, C)
  return(partysplit(varid = as.integer(Gene), 
                    breaks = Output$expressionSplit,
                    index = c(1L,2L),
                    info = list(treatments = c(names(Output$T1), names(Output$T2)))))
}

getTreatment <- function(C, weights) {
  # C: IC50
  C <- C[weights == 1,,drop=F]
  sumtreat <- colSums2(C)
  biomk <- which.min(sumtreat)
  names(biomk) <- colnames(C)[biomk]
  return(biomk)
}


growtreeExp <- function(id = 1L, response, data, minbucket = 10, weights = NULL) {
  if (is.null(weights))  weights <- rep(1L, nrow(data))
  
  ## for less than minbucket observations stop here
  if (sum(weights) < minbucket) return(partynode(id = id, info = names(getTreatment(response, weights))))
  
  ## find best split
  sp <- findsplitExp(response, data, minimum = minbucket, weights)
  datadf <- as.data.frame(data)
  ## no split found, stop here
  if (is.null(sp)) return(partynode(id = id, info = names(getTreatment(response, weights))))
  
  ## actually split the data
  kidids <- kidids_split(sp, data = datadf)
  if(length(unique(kidids[weights==1]))==1) 
    return(partynode(id = id, info = names(getTreatment(response, weights))))
  if(min(table(kidids[weights==1])) < minbucket) 
    return(partynode(id = id, info = names(getTreatment(response, weights))))
  
  ## set up all daugther nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))
  for (kidid in 1:length(kids)) {
    ## select observations for current node
    w <- weights
    w[kidids != kidid] <- 0
    ## get next node id
    if (kidid > 1) {
      myid <- max(nodeids(kids[[kidid - 1]]))
    } else {
      myid <- id
    }
    ## start recursion on this daugther node
    kids[[kidid]] <- growtreeExp(id = as.integer(myid + 1), response, 
                                 data, weights =w, minbucket = minbucket)
  }
  
  ## return nodes
  return(partynode(id = as.integer(id), split = sp, kids = kids,
                   info = ""))
}

getTreatment <- function(C, weights) {
  # C: IC50
  C <- C[weights == 1,,drop=F]
  sumtreat <- colSums2(C)
  biomk <- which.min(sumtreat)
  names(biomk) <- colnames(C)[biomk]
  return(biomk)
}

getSplit <- function(gene, Expression, C) {
  I <- order(Expression[,gene])
  Id <- order(Expression[,gene], decreasing = T)
  
  M1 <- rbind(0,colCumsums(C[I,]))
  M2 <- rbind(0,colCumsums(C[Id,]))
  M2 <- M2[c(nrow(M1):1),,drop=F]
  
  
  # plot(rowMins(M1) + rowMins(M2))
  th <- which.min(rowMins(M1) + rowMins(M2))
  sumic50 <- min(rowMins(M1) + rowMins(M2))
  T1 = max.col(-M1)[th]
  names(T1) = colnames(C)[T1]
  T2 = max.col(-M2)[th]
  names(T2) = colnames(C)[T2]
  
  sumic50
  return(list(sumic50 = sumic50, 
              T1 = T1,
              T2 = T2,
              split = th,
              expressionSplit = Expression[I[th],gene]))
}

getsumic50 <- function(gene, Expression, C) {
  I <- order(Expression[,gene])
  Id <- order(Expression[,gene], decreasing = T)
  
  M1 <- rbind(0,colCumsums(C[I,]))
  M2 <- rbind(0,colCumsums(C[Id,]))
  M2 <- M2[c(nrow(M1):1),,drop=F]
  
  th <- which.min(rowMins(M1) + rowMins(M2))
  sumic50 <- min(rowMins(M1) + rowMins(M2))
  # T1 = max.col(-M1)[th]
  # names(T1) = colnames(C)[T1]
  # T2 = max.col(-M2)[th]
  # names(T2) = colnames(C)[T2]
  return(sumic50)
}

getsumic50v2 <- function(gene, Expression, C) {
  I <- order(Expression[,gene])
  Id <- order(Expression[,gene], decreasing = T)
  
  M1 <- colCumsums(C[I,])
  M2 <- colCumsums(C[Id,])
  M2 <- M2[c(nrow(M2):1),,drop=F]
  dummy <- c(0,rowMins(M1)) + c(rowMins(M2),0)
  sumic50 <- min(dummy)
  # T1 = max.col(-M1)[th]
  # names(T1) = colnames(C)[T1]
  # T2 = max.col(-M2)[th]
  # names(T2) = colnames(C)[T2]
  return(sumic50)
}

getsumic50v3 <- function(geneExpression, C) {
  I <- order(geneExpression)
  # Id <- order(geneExpression, decreasing = T) # Equivalent to
  Id <- I[length(I):1]
  
  M1 <- colCumsums(C[I,])
  M2 <- colCumsums(C[Id,])
  M2 <- M2[c(nrow(M2):1),,drop=F]
  dummy <- c(0,rowMins(M1)) + c(rowMins(M2),0)
  sumic50 <- min(dummy)
  return(sumic50)
}


#############################################################
# ExpressionTrees
#############################################################

library(partykit)

trainTree <- function(Expression, C, minbucket = 20) {
  nodes <- growtreeExp(id = 1L, C, Expression, minbucket = minbucket)
  tree <- party(nodes, data = as.data.frame(Expression))
  tree$node$info <- NULL
  return(tree)
}

library(ggparty)
library(data.tree)

niceTree <- function(tree, folder = NULL) {
  treeNode <- as.Node(tree)
  colors <- c("", "#367592", "#39A7AE", "#96D6B6", "#FDE5B0", "#F3908B", "#E36192", "#8E4884", "#A83333")

  levels <- max(treeNode$Get(function(x) c(level = x$level)))
  # bucle
  for (i in 1:levels) {
    color_node <- colors[(i %% length(colors)) + 1]
    nodes_level <- Traverse(treeNode, filterFun = function(x) x$level == i)
    Do(nodes_level, function(node) {
      SetNodeStyle(node,
        label = function(node) paste0(node$splitname),
        tooltip = function(node) paste0(nrow(node$data), " observations"),
        fontname = "Roboto",
        shape = "diamond",
        style = "filled",
        color = color_node,
        fillcolor = paste(color_node, "88", sep = ""),
        fontcolor = "black"
      )
    })
  }

  SetEdgeStyle(treeNode,
    arrowhead = "none",
    label = function(node) node$splitLevel,
    fontname = "Roboto",
    penwidth = function(node) 12 * nrow(node$data) / nrow(node$root$data),
  )

  Do(treeNode$leaves, function(node) SetNodeStyle(node, shape = "box"))

  if (!is.null(folder)) {
    # Save plot
    tmp <- DiagrammeRsvg::export_svg(plot(treeNode))
    tmp <- charToRaw(tmp) # flatten
    rsvg::rsvg_png(tmp, paste(folder, "/tree_plot.png", sep = ""),height = 2000)
  }

  return(treeNode)
}

predictTree <- function(tree,C, Expression) {
  treatments <- unlist(nodeapply(tree, 
                                 predict.party(tree, as.data.frame(Expression))
                                 , info_node))
  TratamientoTree <- match(treatments, colnames(C))
  TratamientoTree <- factor(TratamientoTree, levels = 1:ncol(C))
}

#############################################################
# MutationTrees
#############################################################

findsplitMut <- function(C, X, minimum = 1, weights) {
  # C: IC50
  # X: mutation matrix
  C <- C[weights == 1,]
  X <- X[weights == 1,]
  if (nrow(C) <= minimum) return(NULL)
  
  mut <- max(X) # Mutation is the largest value of X
  WT <- min(X) # Wildtype is the smalles value of X
  
  tA <- t(C) %*% (X==mut) # treatment with mutation
  tB <- t(C) %*% (X==WT)  # Treatment for WT
  
  #   tA Sum of the IC50 of the mutated samples treated with each treatment
  #   tB Sum of the IC50 of the wild type samples treated with each treatment
  
  
  # The biomarker mutation
  optimal <- colMins(tA) + colMins(tB) # Optimal selection of treatments for each mutation
  biomk <- which.min(optimal) # Optimal mutation
  names(biomk) <- colnames(tA)[biomk]
  
  # If both mutated and WT samples return the same treatment, don't split
  if (length(unique(X[,biomk])) == 1) return(NULL) 
  
  # The treatments
  T1Mut <- which.min(tA[,biomk])
  T1WT <- which.min(tB[,biomk])
  
  # If both mutated and WT samples return the same treatment, don't split. Probably already done.
  if (T1Mut == T1WT) return(NULL)
  
  index = c(1L,2L)
  names(index) <- c(paste0(names(biomk),"mut"), paste0(names(biomk),"WT"))
  return(partysplit(varid = as.integer(biomk), 
                    index = index,
                    info = list(treatments = c(names(T1Mut), names(T1WT)))))
}

growtreeMut <- function(id = 1L, response, data, minbucket = 10, weights = NULL, findsplit = findsplitMut) { # nolint: line_length_linter.
  if (is.null(weights))  weights <- rep(1L, nrow(data))
  
  # Parameters:
  # id: unique Id of the node
  # response: IC50, or in general, the funciton to minminze.
  # data: Expression, mutations, or in general, the data to do the split
  # minbucket: minimum number of samples in child node.
  # weights: don't set it. Variable to state which are the samples under study.
  
  ## for less than "minbucket"  observations stop here
  if (sum(weights) < minbucket) return(partynode(id = id, info = names(getTreatment(response, weights)))) # nolint: line_length_linter.
  
  ## find best split
  sp <- findsplitMut(response, data, minimum = minbucket, weights)
  
  ## no split found, stop here
  if (is.null(sp)) return(partynode(id = id, info = names(getTreatment(response, weights))))
  
  
  ## actually split the data
  datadf <- as.data.frame(data)
  kidids <- kidids_split(sp, data = datadf)
  
  ## If only one kid, return
  if(length(unique(kidids))==1) 
    return(partynode(id = id, info = names(getTreatment(response, weights))))
  
  ## If any of the splits smaller than minbucket, return
  if(min(table(kidids[weights==1,drop=F])) < minbucket) 
    return(partynode(id = id, info = names(getTreatment(response, weights))))
  
  ## set up all daugther nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))
  for (kidid in 1:length(kids)) {
    ## select observations for current node
    w <- weights
    w[kidids != kidid] <- 0
    ## get next node id
    if (kidid > 1) {
      myid <- max(nodeids(kids[[kidid - 1]]))
    } else {
      myid <- id
    }
    ## start recursion on this daugther node
    kids[[kidid]] <- growtreeMut(id = as.integer(myid + 1), response, 
                                 data, weights =w, minbucket = minbucket, findsplit = findsplitMut)
  }
  
  ## return nodes
  return(partynode(id = as.integer(id), split = sp, kids = kids,
                   info = ""))
}


trainTreeMut <- function(Expression, C, minbucket = 20) {
  Expression <- Expression - min(Expression) + 1L
  mode(Expression) <- "integer"
  nodes <- growtreeMut(id = 1L, C, Expression, minbucket = minbucket)
  tree <- party(nodes, data = as.data.frame(Expression))
  tree$node$info <- NULL
  return(tree)
}

predictTreeMut <- function(tree,C, Expression) {
  Expression <- Expression - min(Expression) + 1L
  mode(Expression) <- "integer"
  predictTree(tree,C, Expression)
}

#############################################################
# Kernel Ridge Regression
#############################################################

trainKRR <- function(Expression, C, sigma=quantile(colSds(Expression),.1), lambda=.1/nrow(Expression)) {
  krr = constructKRRLearner()
  p = list(kernel="rbfdot", sigma=quantile(colSds(Expression),.1), lambda=.1/nrow(Expression))
  modellist <- list(krr = krr, p = p)
  for (n in 1:ncol(C)) {
    ns$y <- C[,n]
    ns$x <- Expression
    m <- krr$learn(ns, p)
    modellist[[n+2]] <- m
  }
  return(modellist)
}

predictKRR <- function(Salida, response, Expression) {
  Prediccion <- matrix(NA, nrow = nrow(Expression), ncol = (length(Salida)-2))
  krr <- constructKRRLearner()
  ns <- constructData(Expression, y = NA)
  for (n in 1:(length(Salida)-2)) {
    Prediccion[,n] <- krr$predict(Salida[[n+2]],ns)
  }
  TratamientoKRR <- apply(Prediccion, 1, which.min)
  TratamientoKRR <- factor(TratamientoKRR, levels = 1:ncol(response))
  return(TratamientoKRR)
  
}




#############################################################
# MultinomialLasso
#############################################################

trainMnLasso <- function(Expression, Y) {
  Salida <- cv.glmnet(Expression,Y, family = "multinomial", 
                      alpha = 1, type.multinomial = "grouped")
  return(Salida)
}

predictMnLasso <- function(Salida, response, Expression) {
  Prediccion <- drop(predict(Salida,Expression, s="lambda.min"))
  TratamientoMnLasso <- apply(Prediccion, 1, which.max)
  TratamientoMnLasso <- factor(TratamientoMnLasso, levels = 1:ncol(response))
  return(TratamientoMnLasso)
}

#############################################################
# Lasso
#############################################################

trainLasso <- function(Expression, C) {
  Salida <- cv.glmnet(Expression, C, family = "mgaussian")
  return(Salida)
}

predictLasso <- function(Salida, response, Expression) {
  Prediccion <- drop(predict(Salida,Expression, s="lambda.min"))
  TratamientoLasso <- apply(Prediccion, 1, which.min)
  TratamientoLasso <- factor(TratamientoLasso, levels = 1:ncol(response))
  return(TratamientoLasso)
}

#############################################################
# BOSO
#############################################################

trainBOSO <- function(Expression, Y,...) {
  x <- Expression[c(T,T,F),]
  xval <- Expression[c(F,F,T),]
  
  coeffs <- matrix(NA, ncol = ncol(Expression)+1, nrow = ncol(Y))
  start_time = Sys.time()
  for (n in 1:ncol(Y)) {
    cat(colnames(Y)[n], "\n")
    y <- Y[c(T,T,F),n]
    yval <- Y[c(F,F,T),n]
    obj <- myBOSO(x,y,xval,yval,
                  nlambda=50,
                  intercept= T,
                  Threads=4, verbose = 1, seed = 2021,...)
    coeffs[n,] <- coef(obj)
    tiempo <-Sys.time() - start_time 
    cat("Elapsed time:", tiempo, attr(tiempo, "units"),"\n")
  }
  return(coeffs)
}

predictBOSO <- function(coeffs, response, Expression) {
  Prediccion <- (t(t(Expression %*% t(coeffs[,-1])) + coeffs[,1]))
  TratamientoBOSO <- apply(Prediccion, 1, which.min)
  TratamientoBOSO <- factor(TratamientoBOSO, levels = 1:ncol(response))
  return(TratamientoBOSO)  
}

plotIntragroup <- function(chosen_method, color, Validation_Beat, drug_response, dataset) {
  # Function for plotting intragroup validation boxplots
  Plot_aux_method <- Validation_Beat[Validation_Beat$Method == chosen_method, ] # matriz de mutaciones que sea distinta tiene que tener 64
  Plot_aux_oracle <- Validation_Beat[Validation_Beat$Method == "ORACLE", ]
  Plot_aux <- merge(Plot_aux_method, Plot_aux_oracle, by = "patients")
  # Plot_aux$delta <- Plot_aux$IC50.x - Plot_aux$IC50.y
  
  if (length(which(is.na(Plot_aux$drug.x))) > 0) {
    Plot_aux <- Plot_aux[-which(is.na(Plot_aux$drug.x)), ] # delete patient if no drug is assigned
  }
  
  drug_plots <- NULL
  
  # Loop through each unique drug in the data set
  for (drug_idx in 1:length(unique(Plot_aux$drug.x))) {
    current_drug <- unique(Plot_aux$drug.x)[drug_idx] # drug name
    is_current_drug <- Plot_aux$drug.x == current_drug # indexes of cells with drug
    
    # Create an empty data frame for the plot
    current_drug_data <- data.frame(
      DrugIndex = drug_idx,
      DrugName = current_drug,
      Patients = Plot_aux$patients,
      Recommended = "No",
      IC50 = NA
    )
    
    drug_response <- drug_response[Plot_aux$patients, ]
    if (current_drug %in% colnames(drug_response)) {
      current_drug_data$Recommended[is_current_drug] <- "Yes"
      current_drug_data$Recommended[!is_current_drug] <- "No"
      current_drug_data$IC50[is_current_drug] <- drug_response[is_current_drug, current_drug] - Plot_aux$IC50.y[is_current_drug]
      current_drug_data$IC50[!is_current_drug] <- drug_response[!is_current_drug, current_drug] - Plot_aux$IC50.y[!is_current_drug]
      
    }
    
    # Add the current drug data to the results data frame
    drug_plots <- rbind(drug_plots, current_drug_data)
  }
  
  # Create the plot using ggplot2
  # drug_plots1 <- drug_plots %>%
  #   group_by(DrugName, Recommended) %>%
  #   mutate(N = n()) %>%
  #   mutate(N = ifelse(IC50 == min(IC50, na.rm = T), paste0("n=", N), NA))
  
  # Default 12 character target width.
  swr = function(string, nwrap=12) {
    paste(strwrap(string, width=nwrap), collapse="\n")
  }
  swr = Vectorize(swr)
  
  # Create line breaks in DrugName
  drug_plots$DrugName <- swr(drug_plots$DrugName)
  
  gg_drug_plots <- ggplot(drug_plots, aes(x = Recommended, y = IC50, fill = Recommended)) +
    geom_boxplot(outlier.shape = NA, lwd = 0.8) +
    # geom_text(nudge_y = -0.6, size = 6) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, alpha=0.8)+
    facet_wrap(~DrugName) +
    # geom_jitter(alpha = 0.9, colour = "black", size = 0.9) +
    theme_bw() +
    scale_fill_manual(values = c("gray", color)) + # Change colors to gray and blue
    stat_compare_means(vjust = 1, hjust = 0, hide.ns = FALSE, label = "p", size = 5) +
    ylab(expression(Delta * IC50)) +
    xlab("") +
    ggtitle(chosen_method) +
    theme(
      text = element_text(size = 22, face = "bold", family = "Roboto"), # Increase font size and make drug names bold
      panel.background = element_rect(fill = "transparent"), # Set background to transparent
      plot.background = element_rect(fill = "transparent"),
      legend.background = element_rect(fill = "transparent"),
      title = element_text(size = 22),
      strip.text = element_text(size=15),
      legend.position = "delete"
      
    )
  
  num_plots <- length(table(drug_plots$DrugName))
  if (num_plots > 12) {
    num_plots <- 12
  }
  # if (num_plots < 5) {
  #   num_plots <- 7
  # }
  # ggsave(paste(folder_dir, "/images/intragroup_", chosen_method, "_", dataset, ".png", sep = ""), gg_drug_plots, width = num_plots * 1.6, height = num_plots * 1.3, dpi = 1000)
  ggsave(paste(folder_dir, "/images/intragroup_", chosen_method, "_", dataset, ".png", sep = ""), gg_drug_plots, width = 12, height =12, dpi = 1000) # width = 5.5, height = 7.5 for beataml and h=6 for GDSC
  return(gg_drug_plots)
}
