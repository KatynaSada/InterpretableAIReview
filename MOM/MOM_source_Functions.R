require(IHW, quietly = T)
require(impute, quietly = T)
require(dplyr, quietly = T)
require(tidyverse, quietly = T)

getBMeffect_F <- function(Drug, Mutations) {
  # Obtaining Wilcoxon P-Values-----------------------------------------------
  PVal <- data.frame(
    Biomarker = rep(colnames(Mutations), each = ncol(Drug)),
    Drug = rep(colnames(Drug), ncol(Mutations)),
    P.Value = NA
  )

  print("Getting Biomarker-drugs deppendencies")
  for (k in 1:nrow(PVal)) {
    mutsel <- Mutations[, as.character(PVal$Biomarker[k])]
    dsel_mut <- Drug[which(rownames(Drug) %in% names(mutsel[mutsel == 1])), as.character(PVal$Drug[k])]
    dsel_WT <- Drug[which(rownames(Drug) %in% names(mutsel[mutsel == 0])), as.character(PVal$Drug[k])]
    p <- wilcox.test(as.numeric(dsel_mut), as.numeric(dsel_WT))$p.value
    PVal$P.Value[k] <- p

    PVal$DeltaIC50[k] <- median(as.matrix(dsel_WT)) - median(as.matrix(dsel_mut))
    if (which.min(c(median(as.matrix(dsel_WT)), median(as.matrix(dsel_mut)))) == 2) {
      PVal$MutStatus_TE[k] <- "Mut"
    } else {
      PVal$MutStatus_TE[k] <- "WT"
    }
  }

  PVal <- PVal[order(PVal$P.Value, decreasing = F), ]

  print("Correcting Hypothesis")
  PVal$Biomarker <- as.factor(PVal$Biomarker)

  ihwRes <- suppressMessages({
    ihw(
      pvalues = PVal$P.Value,
      covariates = as.factor(PVal$Biomarker),
      alpha = 0.2,
      covariate_type = "nominal",
      # nbins = length(unique(PVal$Biomarker)),
      m_groups = NULL,
      quiet = TRUE,
      nfolds = 1L,
      nfolds_internal = 1L,
      nsplits_internal = 1L,
      lambdas = "auto",
      seed = 3L,
      distrib_estimator = "grenander",
      lp_solver = "lpsymphony",
      adjustment_type = "BH",
      return_internal = FALSE
    )
  })

  PVal$IHW_pvalue <- adj_pvalues(ihwRes)

  ihwResDf <- as.data.frame(ihwRes)

  aux <- ihwResDf %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(median(weight))
  colnames(aux)[c(1, 2)] <- c("Tumor", "Med.Weight")
  aux$Tumor <- as.character(aux$Tumor)
  aux <- aux %>% filter(Med.Weight > 0)

  aux$index <- 1:nrow(aux)
  rownames(aux) <- NULL
  Rnk_mutations <- aux

  rm(aux)

  print("Dephicering Biomarker Effect")

  PVal$histogram <- PVal$P.Value * .5 * (PVal$DeltaIC50 > 0) + (1 - PVal$P.Value) * .5 * (PVal$DeltaIC50 <= 0)
  PVal$Effect <- "Unknown"
  PVal$Score <- NA

  for (i in as.character(unique(PVal$Biomarker[which(PVal$Biomarker %in% Rnk_mutations$Tumor)]))) {
    if (2 * sum(2 * (PVal$histogram[PVal$Biomarker == i]) < 0.2) > sum(2 * (PVal$histogram[PVal$Biomarker == i]) == 0.2) & 2 * sum(2 * (PVal$histogram[PVal$Biomarker == i]) < 0.2) > sum(2 * (PVal$histogram[PVal$Biomarker == i]) > 0.2)) {
      PVal$Effect[PVal$Biomarker == i] <- "Sensitive"
      PVal$Score[PVal$Biomarker == i] <- 100 * (2 * sum(2 * (PVal$histogram[PVal$Biomarker == i]) < 0.05)) / nrow(PVal)
    }

    if (2 * sum(2 * (PVal$histogram[PVal$Biomarker == i]) > 0.8) > sum(2 * (PVal$histogram[PVal$Biomarker == i]) < 0.8) & sum(2 * (PVal$histogram[PVal$Biomarker == i]) > 0.8) > sum(2 * (PVal$histogram[PVal$Biomarker == i]) == 0.8)) {
      PVal$Effect[PVal$Biomarker == i] <- "Resistant"
      PVal$Score[PVal$Biomarker == i] <- 100 * (2 * sum(2 * (PVal$histogram[PVal$Biomarker == i]) > 0.95)) / nrow(PVal)
    }
  }


  PVal$Effect <- factor(PVal$Effect, levels = c("Unknown", "Resistant", "Sensitive"))

  return(PVal)
}

Get_Treatment_OptimizationMatrix <- function(BM_Result, Drug, Mutations) {
  Optmization_matrix <- vector(mode = "list", length = 2)

  ET <- BM_Result[(!(BM_Result$Effect %in% "Unknown") & BM_Result$IHW_pvalue < 0.05 & abs(BM_Result$DeltaIC50) > 0.2), ]
  ET$Code <- paste0(ET$Drug, "_", ET$Biomarker, "_", ET$MutStatus_TE)
  TxP <- matrix(0, ncol = nrow(Drug), nrow = nrow(ET))
  rownames(TxP) <- ET$Code
  colnames(TxP) <- rownames(Drug)

  for (l in 1:dim(ET)[1]) {
    g <- as.character(ET$Biomarker[l])
    for (m in 1:dim(Drug)[1]) {
      p <- rownames(Drug)[m]
      sol <- Mutations[p, g]
      if (sol > 0 & as.character(ET$MutStatus_TE[l]) == "Mut") {
        TxP[l, p] <- 1
      }
      if (sol < 1 & as.character(ET$MutStatus_TE[l]) == "WT") {
        TxP[l, p] <- 1
      }
    }
  }

  Response <- TxP


  d <- as.character(ET$Drug[ET$Code %in% rownames(TxP)])

  aux <- which(Response > 0, arr.ind = T)
  aux2 <- d[aux[, 1]]
  pat <- colnames(Response)[aux[, 2]]

  for (j in 1:dim(aux)[1]) {
    # If the patient has the biomarker then we select the effectiveness of the treatment
    if (Response[aux[j, 1], aux[j, 2]] > 0) {
      Response[aux[j, 1], aux[j, 2]] <- Drug[pat[j], aux2[j]]
    }
  }

  Optmization_matrix[[1]] <- t(Response)
  Optmization_matrix[[2]] <- t(TxP)
  names(Optmization_matrix) <- c("IC50*", "PT")
  return(Optmization_matrix)
}

MOM_MILP <- function(IC50 = IC50_log, PT = PT, S = S) {
  # Optimization model------------------------------------------
  # repl_python()
  py_run_string("import numpy as np")
  py_run_string("from docplex.mp.model import Model")
  # Total number of Steps to be considered
  py_run_string("S=r.S")
  py_run_string("S=int(S)")
  py_run_string("s=[i for i in range(1, S)]")
  py_run_string("s=[0]+s")
  py_run_string("print('s=',s)")
  # Total number of patients in the study
  py_run_string("P=r.P ")
  py_run_string("p=[i for i in range(1, P)]")
  py_run_string("p=[0]+p")
  py_run_string("print('P, number of patients equals ',P)")
  # Total number of possible treatments in the example
  py_run_string("T=r.Treat ")
  py_run_string("t=[i for i in range(1, T)]")
  py_run_string("t=[0]+t")
  py_run_string("print('T number of treatments equals ',T)")
  py_run_string("IC50_log=r.IC50")

  py_run_string("PT=r.PT")

  # MODEL DEFINITION
  py_run_string("mdl2=Model('AML_Real_Case')")
  # We define touples to create
  py_run_string("A=[(i,j,k) for i in s for j in t for k in p]")
  # the dictionary variable
  py_run_string("B=[(i,j) for i in s for j in t]")

  # Cplex internal binary variable definition
  py_run_string("Xstp=mdl2.binary_var_dict(A, name='Xstp')")
  py_run_string("Yst=mdl2.binary_var_dict(B, name='Yst')")

  # Objective Function
  py_run_string("mdl2.minimize(mdl2.sum(IC50_log[k][j]*Xstp[i,j,k]  for k in p for j in t for i in s)) ")

  # Contraints of the model
  py_run_string("mdl2.add_constraints_(mdl2.sum(Yst[i,j] for j in t)== 1 for i in s)")
  py_run_string("mdl2.add_constraints_(mdl2.sum(Xstp[i,j,k] for i in s for j in t)<=1 for k in p)")
  py_run_string("mdl2.add_constraints_(Xstp[i,j,k]<= PT[k][j]*Yst[i,j] for i in s  for j in t for k in p)")
  py_run_string("mdl2.add_constraints_((Xstp[i,j,k] + mdl2.sum(Xstp[m,n,k] for m in range(i) for n in t) >=Yst[i,j]*PT[k][j]) for i in s for j in t for k in p)")
  # mdl2.solve()
  py_run_string("print('Loaded')")

  # Solve the model
  py_run_string("solucion=mdl2.solve(log_output=True)")
  # Prints the computing details
  py_run_string("print(mdl2.get_solve_details())")
  # Prints the solution status: Is it feasible?
  py_run_string("print(mdl2.get_solve_status())")
  # Shows the solution

  py_run_string("Xstp_active=[i for i in A if Xstp[i].solution_value>0]")
  py_run_string("Yst_active=[j for j in B if Yst[j].solution_value>0]")

  # quit

  sol_pat <- py$Xstp_active

  return(sol_pat)
}

Extract_results <- function(Solution_MILP, PT) {
  Results_MILP <- vector(mode = "list", length = 2)
  nst <- NULL
  pos <- NULL
  patient <- matrix(0, nrow = length(Solution_MILP), ncol = 2)

  colnames(patient) <- c("TreatmentCode", "PatientCode")
  patient <- as.data.frame(patient)

  print("Extracting Subindexes")
  for (i in 1:length(Solution_MILP)) {
    nst <- c(nst, Solution_MILP[[i]][[1]])
    pos <- c(pos, Solution_MILP[[i]][[2]])
    patient$TreatmentCode[i] <- Solution_MILP[[i]][[2]]
    patient$PatientCode[i] <- Solution_MILP[[i]][[3]]
  }

  aux <- as.matrix(table(nst))
  patient <- patient + 1
  nst <- unique(nst)
  pos <- unique(pos) + 1


  patient$PatientCode <- rownames(PT)[patient$PatientCode]

  Results <- cbind(1:length(nst), colnames(PT)[pos], pos, aux)
  colnames(Results) <- c("Step", "Treatment Code", "Position in PT matrix", "Number of Patients Treated")
  Results <- as.data.frame(Results)
  patient$TreatmentName <- as.character(Results$`Treatment Code`)[match(
    patient$TreatmentCode,
    Results$`Position in PT matrix`
  )]
  Results_MILP[[1]] <- Results
  Results_MILP[[2]] <- patient
  names(Results_MILP) <- c("TreatmentGuideline", "Patient-Distribution")
  return(Results_MILP)
}

CV_Prediction2 <- function(MILP_classification, kfold_validation_drug, kfold_validation_mut) {
  # Corrected
  GL <- MILP_classification$TreatmentGuideline
  selectedTreat <- kfold_validation_mut[, 1, drop = FALSE] * NA
  colnames(selectedTreat) <- "Treatment"
  for (s in nrow(GL):1) {
    selectedDrug <- unlist(strsplit(GL[s, 2], split = "_"))[1]
    mut_step <- unlist(strsplit(GL[s, 2], split = "_"))[2]
    status <- factor(unlist(strsplit(GL[s, 2], split = "_"))[3], levels = c("WT", "Mut"))
    Poner <- kfold_validation_mut[, mut_step] == as.numeric(status) - 1
    selectedTreat[Poner, ] <- selectedDrug
  }
  if (!is.null(kfold_validation_drug)) {
    IC50 <- kfold_validation_drug[cbind(1:length(selectedTreat), match(selectedTreat, colnames(kfold_validation_drug)))]
    selectedTreat <- data.frame(selectedTreat, IC50 = IC50)
  }
  return(selectedTreat)
}

CV_Prediction <- function(MILP_classification, kfold_validation_drug, kfold_validation_mut) {
  Results <- MILP_classification$TreatmentGuideline

  validation_fold <- data.frame(
    Drug = sapply(as.character(Results$`Treatment Code`), FUN = function(X) {
      unlist(strsplit(X, split = "_"))[1]
    }),
    Mut = as.character(sapply(as.character(Results$`Treatment Code`), FUN = function(X) {
      unlist(strsplit(X, split = "_"))[2]
    })),
    MStatus = sapply(as.character(Results$`Treatment Code`), FUN = function(X) {
      unlist(strsplit(X, split = "_"))[3]
    }),
    Treatment = as.character(Results$`Treatment Code`)
  )

  validation_fold$Patient <- NULL

  validation_ic50 <- NULL

  for (i in 1:nrow(validation_fold)) {
    if (validation_fold$MStatus[i] == "Mut") {
      p <- rownames(kfold_validation_mut)[which(kfold_validation_mut[, as.character(validation_fold$Mut[i])] == 1)]
      ic50 <- kfold_validation_drug[p, validation_fold$Drug[i]]
    }
    if (validation_fold$MStatus[i] == "WT") {
      p <- rownames(kfold_validation_mut)[which(kfold_validation_mut[, as.character(validation_fold$Mut[i])] == 0)]
      ic50 <- kfold_validation_drug[p, validation_fold$Drug[i]]
    }
    aux <- data.frame(
      Drug = rep(validation_fold$Treatment[i], length(ic50)),
      Mut = rep(validation_fold$Mut[i], length(ic50)),
      IC50 = as.numeric(as.character(ic50)),
      Pat_name = as.character(p)
    )

    validation_ic50 <- rbind(validation_ic50, aux)

    kfold_validation_drug <- kfold_validation_drug[-c(which(rownames(kfold_validation_drug) %in% p)), , drop = FALSE]
    kfold_validation_mut <- kfold_validation_mut[-c(which(rownames(kfold_validation_mut) %in% p)), , drop = FALSE]
  }
  return(validation_ic50)
}