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
  if (sum(weights) < minbucket) return(partynode(id = id, info = names(getTreatment(C, weights))))
  
  ## find best split
  sp <- findsplitExp(response, data, minimum = minbucket, weights)
  datadf <- as.data.frame(data)
  ## no split found, stop here
  if (is.null(sp)) return(partynode(id = id, info = names(getTreatment(C, weights))))
  
  ## actually split the data
  kidids <- kidids_split(sp, data = datadf)
  if(length(unique(kidids[weights==1]))==1) 
    return(partynode(id = id, info = names(getTreatment(C, weights))))
  if(min(table(kidids[weights==1])) < minbucket) 
    return(partynode(id = id, info = names(getTreatment(C, weights))))
  
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

niceTree <- function(tree) {
  treeNode <- as.Node(tree)
  SetNodeStyle(treeNode,
               label = function(node) paste0(node$splitname), 
               tooltip = function(node) paste0(nrow(node$data), " observations"),
               fontname = "Roboto",
               shape = "diamond")
  
  SetEdgeStyle(treeNode,
               arrowhead = "none",
               label = function(node) node$splitLevel,
               fontname = "Roboto",
               penwidth = function(node) 12 * nrow(node$data)/nrow(node$root$data),  )
  Do(treeNode$leaves, function(node) SetNodeStyle(node, shape = "box"))
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

growtreeMut <- function(id = 1L, response, data, minbucket = 10, weights = NULL, findsplit = findsplitMut) {
  if (is.null(weights))  weights <- rep(1L, nrow(data))
  
  # Parameters:
  # id: unique Id of the node
  # response: IC50, or in general, the funciton to minminze.
  # data: Expression, mutations, or in general, the data to do the split
  # minbucket: minimum number of samples in child node.
  # weights: don't set it. Variable to state which are the samples under study.
  
  ## for less than "minbucket"  observations stop here
  if (sum(weights) < minbucket) return(partynode(id = id, info = names(getTreatment(drug_response_matrix, weights))))
  
  ## find best split
  sp <- findsplitMut(response, data, minimum = minbucket, weights)
  
  ## no split found, stop here
  if (is.null(sp)) return(partynode(id = id, info = names(getTreatment(drug_response_matrix, weights))))
  
  
  ## actually split the data
  datadf <- as.data.frame(data)
  kidids <- kidids_split(sp, data = datadf)
  
  ## If only one kid, return
  if(length(unique(kidids))==1) 
    return(partynode(id = id, info = names(getTreatment(drug_response_matrix, weights))))
  
  ## If any of the splits smaller than minbucket, return
  if(min(table(kidids[weights==1,drop=F])) < minbucket) 
    return(partynode(id = id, info = names(getTreatment(drug_response_matrix, weights))))
  
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
