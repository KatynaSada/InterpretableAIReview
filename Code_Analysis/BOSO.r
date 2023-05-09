myBOSO <- function (x, y, xval, yval, IC = "eBIC", IC.blocks = NULL, nlambda = 100, 
          nlambda.blocks = 10, lambda.min.ratio = ifelse(nrow(x) < 
                                                           ncol(x), 0.01, 1e-04), lambda = NULL, intercept = TRUE, 
          standardize = TRUE, dfmax = NULL, maxVarsBlock = 10, costErrorVal = 1, 
          costErrorTrain = 0, costVars = 0, Threads = 0, timeLimit = 1e+75, 
          verbose = F, seed = NULL, warmstart = F, TH_IC = 0.001, 
          indexSelected = NULL) 
{
  if (!requireNamespace("cplexAPI", quietly = TRUE)) {
    stop("Package cplexAPI not installed (required here)!", 
         call. = FALSE)
  }
  if (!is(x, "matrix") & !is(x, "Matrix")) {
    stop("input x must be a matrix or a Matrix class")
  }
  if (!is(y, "numeric") & !is(y, "matrix") & !is(y, "array")) {
    stop("input y must be numeric")
  }
  if (!is(xval, "matrix") & !is(xval, "Matrix")) {
    stop("input xval must be a matrix or a Matrix class")
  }
  if (!is(yval, "numeric") & !is(yval, "matrix") & !is(yval, 
                                                       "array")) {
    stop("input yval must be numeric")
  }
  if (!is(IC, "character")) {
    stop("information criterion metric must be character")
  }
  if (!is(nlambda, "numeric")) {
    stop("nlambda must be numeric")
  }
  if (!is(nlambda.blocks, "numeric")) {
    stop("nlambda.blocks must be numeric")
  }
  if (!is(lambda.min.ratio, "numeric")) {
    stop("lambda.min.ratio must be numeric")
  }
  if (!is(nlambda, "numeric")) {
    stop("nlambda must be numeric")
  }
  if (!is(maxVarsBlock, "numeric")) {
    stop("maxVarsBlock must be numeric")
  }
  if (!is(TH_IC, "numeric")) {
    stop("TH_IC must be numeric")
  }
  x = as.matrix(x)
  y = as.numeric(y)
  xval = as.matrix(xval)
  yval = as.numeric(yval)
  n = nrow(x)
  nval = nrow(xval)
  p = ncol(x)
  if (standardize) {
    obj = standardize(x, y, intercept = T, normalize = T)
    x = obj$x
    y = obj$y
    mx = obj$mx
    my = obj$my
    sx = obj$sx
    obj = standardize(xval, yval, mx = mx, my = my, sx = sx)
    xval = obj$x
    yval = obj$y
    intercept = F
    mx = c(0, mx)
    sx = c(1, sx)
  }
  else {
    mx = rep(0, p + 1)
    my = 0
    sx = rep(1, p + 1)
  }
  if (is.null(dfmax)) {
    dfmax = p
  }
  if (is.null(lambda)) {
    lambda_max <- max(abs(t(y - mean(y) * (1 - mean(y))) %*% 
                            x))/n
    lambda_min <- lambda_max * lambda.min.ratio
    lambda <- exp(seq(log(lambda_max * 1000), log(lambda_min), 
                      length.out = nlambda))
    lambda.blocks <- exp(seq(log(lambda_max * 1000), log(lambda_min), 
                             length.out = nlambda.blocks))
  }
  else {
    nlambda <- length(lambda)
    nlambda.blocks <- length(lambda.blocks)
    lambda.blocks <- lambda.blocks
  }
  if (is.null(IC.blocks)) {
    IC.blocks <- ifelse(IC == "eBIC", "BIC", IC)
  }
  data.raw <- list(x = x, y = y, xval = xval, yval = yval)
  if (is.null(indexSelected)) {
    indexSelected <- 1:p
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  FinishProblem <- F
  ContinueBlocks <- T
  numIter <- 0
  indexSelectedIter <- list()
  while (ContinueBlocks & length(indexSelected) > maxVarsBlock * 
         1.5) {
    numIter <- numIter + 1
    idx <- sample(indexSelected)
    k <- ceiling(length(idx)/maxVarsBlock) * maxVarsBlock
    if (k != length(idx)) {
      idx[(length(idx) + 1):k] <- NA
    }
    idx <- matrix(idx, nrow = maxVarsBlock, byrow = T)
    idx <- lapply(lapply(apply(idx, 2, as.list), unlist), 
                  function(x) {
                    x[!is.na(x)]
                  })
    indexSelectedIter[[numIter]] <- list()
    time_block <- rep(NA, length(idx))
    if (verbose) {
      cat(paste0("Iteration ", numIter, "\t number of variables: ", 
                 length(indexSelected), " \t block size: ", maxVarsBlock, 
                 "\n"))
    }
    for (i in 1:length(idx)) {
      p.block <- length(idx[[i]])
      numVarArray <- seq(0, p.block)
      x = data.raw$x[, idx[[i]]]
      xval = data.raw$xval[, idx[[i]]]
      time_block[i] <- Sys.time()
      if (verbose > 1) {
        cat(paste0("Iteration ", numIter, "\t number of variables: ", 
                   length(indexSelected), "\t subproblem ", round(i * 
                                                                    100/length(idx)), "%\n"))
      }
      if (warmstart) {
        result.block <- BOSO.multiple.warmstart(x = x, 
                                                y = y, xval = xval, yval = yval, lambda = lambda.blocks, 
                                                intercept = intercept, standardize = F, dfmin = 0, 
                                                dfmax = p.block, costErrorVal = costErrorVal, 
                                                costErrorTrain = costErrorVal, costVars = costVars, 
                                                Threads = Threads, timeLimit = timeLimit, 
                                                verbose = max(verbose - 2, 0), IC = IC.blocks, 
                                                n.IC = n + nval, p.IC = p, TH_IC = TH_IC)
      }
      else {
        result.block <- BOSO.multiple.coldstart(x = x, 
                                                y = y, xval = xval, yval = yval, lambda = lambda.blocks, 
                                                intercept = intercept, standardize = F, dfmin = 0, 
                                                dfmax = p.block, costErrorVal = costErrorVal, 
                                                costErrorTrain = costErrorVal, costVars = costVars, 
                                                Threads = Threads, timeLimit = timeLimit, 
                                                verbose = max(verbose - 2, 0), IC = IC.blocks, 
                                                n.IC = n + nval, p.IC = p, TH_IC = TH_IC)
      }
      if (nrow(result.block$betas) > length(idx[[i]])) {
        result.block$betas <- result.block$betas[-1, 
        ]
      }
      indexSelectedIter[[numIter]][[i]] <- idx[[i]][result.block$betas[, 
                                                                       which.min(result.block$score)] != 0]
      time_block[i] <- as.numeric(Sys.time() - time_block[i])
      if (verbose > 1) {
        cat(paste0("Iteration ", numIter, "\t number of variables: ", 
                   length(idx[i]), "\t block size: ", maxVarsBlock, 
                   "\tSelected = ", length(indexSelectedIter[[numIter]][[i]]), 
                   "\tElapsed time= ", round(time_block[i], 3), 
                   "\n"))
      }
    }
    indexSelectedIter[[numIter]] <- sort(unlist(indexSelectedIter[[numIter]]))
    indexSelected <- indexSelectedIter[[numIter]]
    if (verbose >= 3) {
      print(indexSelected)
    }
    if (length(indexSelected) <= maxVarsBlock * 1.5) {
      ContinueBlocks = F
    }
    else if (numIter > 2) {
      if (length(indexSelectedIter[[numIter]]) < length(indexSelectedIter[[numIter - 
                                                                           1]])) {
        ContinueBlocks = T
      }
      else if (length(union(indexSelected, indexSelectedIter[[numIter - 
                                                              1]])) == length(indexSelected) & length(union(indexSelected, 
                                                                                                            indexSelectedIter[[numIter - 2]])) == length(indexSelected)) {
        ContinueBlocks = F
      }
    }
  }
  if (length(indexSelected) == 0) {
    result <- list(x = x, y = y, xval = xval, yval = yval, 
                   IC = IC, nlambda = nlambda, lambda = lambda, intercept = intercept, 
                   standardize = standardize, mx = mx, my = my, sx = sx, 
                   dfmax = dfmax, lambda.selected = 0, p = p, n = n, 
                   nval = nval)
    result$betas <- c(mean(y), rep(0, p))
    if (!intercept) {
      result$betas <- result$betas[-1]
    }
    result$betas <- matrix(result$betas, ncol = 1)
    if (intercept) {
      result$errorTrain = y - mean(y)
      result$errorVal = yval - mean(y)
    }
    else {
      result$errorTrain = y
      result$errorVal = yval
    }
    class(result) = "BOSO"
    object <- result
    return(result)
  }
  dfmax.raw <- dfmax
  p.final <- length(indexSelected)
  dfmax <- ifelse(p.final > dfmax, dfmax, p.final)
  numVarArray <- seq(0, dfmax)
  x = data.raw$x[, indexSelected, drop=F]
  xval = data.raw$xval[, indexSelected, drop =F]
  if (verbose) {
    cat(paste0("Final problem:\t number of variables: ", 
               dim(x)[2], " \t \n"))
  }
  if (warmstart) {
    result.final <- BOSO.multiple.warmstart(x = x, y = y, 
                                            xval = xval, yval = yval, lambda = lambda, intercept = intercept, 
                                            standardize = F, dfmin = 0, dfmax = dfmax, costErrorVal = costErrorVal, 
                                            costErrorTrain = costErrorVal, costVars = costVars, 
                                            Threads = Threads, timeLimit = timeLimit, verbose = max(verbose - 
                                                                                                      1, 0), IC = IC, n.IC = n + nval, p.IC = p, TH_IC = TH_IC)
  }
  else {
    result.final <- BOSO.multiple.coldstart(x = x, y = y, 
                                            xval = xval, yval = yval, lambda = lambda, intercept = intercept, 
                                            standardize = F, dfmin = 0, dfmax = dfmax, costErrorVal = costErrorVal, 
                                            costErrorTrain = costErrorVal, costVars = costVars, 
                                            Threads = Threads, timeLimit = timeLimit, verbose = max(verbose - 
                                                                                                      1, 0), IC = IC, n.IC = n + nval, p.IC = p, TH_IC = TH_IC)
  }
  idx = which.min(result.final$score)
  result <- list(x = x, y = y, xval = xval, yval = yval, IC = IC, 
                 nlambda = nlambda, lambda = lambda, intercept = intercept, 
                 standardize = standardize, mx = mx, my = my, sx = sx, 
                 dfmax = dfmax, result.final = result.final, errorTrain = result.final$errorTrain[, 
                                                                                                  idx], errorVal = result.final$errorVal[, idx], lambda.selected = result.final$lambda.selected[idx], 
                 p = p, n = n, nval = nval, blockStrategy = indexSelectedIter)
  result$betas <- rep(0, p + 1)
  result$betas[c(1, indexSelected + 1)] <- result.final$betas[, 
                                                              idx]
  result$betas <- matrix(result$betas, ncol = 1)
  class(result) = "BOSO"
  return(result)
}
