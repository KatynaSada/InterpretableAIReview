# Loading Libraries and Packages --------------------------------------------
#The first part is to load all the libraries and external functions that are required for the analysis. # nolint # nolint

library(readxl)
library(RColorBrewer)
library(matrixStats)
library(partykit)
library(glmnet)
library(BOSO)


source("Code_Analysis/XAIfunctions_NV_new.R")

library(impute)
library(knitr)
library(reticulate)
library(readxl)
library(readr)
library(limma)
library(ggplot2)
library(ggpattern)
library(ggpubr)
library(ggsci)
library(matrixStats)
library(maftools)
library(qvalue)
library(tibble)
library(Rtsne)
library(IHW)
library(matrixStats)
library(pheatmap)
library(DT)
library("dplyr")
library(tictoc)

source("./MOM/2022-05-24_MOM_source_Functions.R")
source("./MOM/predictMOM_ALL.R")

# Loading Data -------------------------------------------- 

#Then, we will load all the data for the analysis.
#Data selected was the BeatAML cohort for which we had mutational data, drug sensitivity data and gene expression data.

# Load mutational data
gene_variants <- read_excel("./data/input/41586_2018_623_MOESM3_ESM.xlsx",
  sheet = "Table S7-Variants for Analysis"
)

# Load clinical data
clinical <- read_excel("./data/input/41586_2018_623_MOESM3_ESM.xlsx",
  sheet = "Tabe S5-Clinical Summary"
)

# Load drug sensitivity data
drug_response <- read_excel("./data/input/41586_2018_623_MOESM3_ESM.xlsx",
  sheet = "Table S10-Drug Responses"
)

# Load gene expression data
expression <- read_excel("./data/input/41586_2018_623_MOESM3_ESM.xlsx",
  sheet = "Table S8-Gene Counts RPKM"
)


# Pre-processing drug and mutational information  --------------------------------------------

# Now it is time to pre-process the data in order to obtain two different matrices, a matrix containing all patients and their mutations available and a second matrix containing all patients and their sensitivity to the different drugs.
# For doing so, we will impute the missing values in the drug matrix.

# Build Drug Matrix: poner los 1s cuando hay mutacion
drug_matrix <- matrix(
  data = "NA", nrow = length(unique(drug_response$inhibitor)),
  ncol = length(unique(drug_response$lab_id))
) # empty matrix

colnames(drug_matrix) <- as.character(unique(drug_response$lab_id))
rownames(drug_matrix) <- unique(drug_response$inhibitor)

for (i in as.character(unique(drug_response$lab_id))) {
  ind <- which(drug_response$lab_id == i)
  D <- as.character(drug_response$inhibitor[ind])
  drug_matrix[D, i] <- drug_response$ic50[ind]
}

# Identifying missing values
is.na(drug_matrix) <- drug_matrix == "NA"

# change number formatting
for (i in 2:(dim(drug_matrix)[2] - 1)) {
  drug_matrix[, i] <- as.numeric(as.character(drug_matrix[, i]))
}
 
# Building the mutations matrix
mutations <- matrix(0,
  nrow = length(unique(gene_variants$labId)),
  ncol = length(unique(gene_variants$symbol))
)

rownames(mutations) <- as.character(unique(gene_variants$labId))
colnames(mutations) <- as.character(unique(gene_variants$symbol))

for (i in as.character(unique(gene_variants$labId))) {
  jj <- unique(as.character(gene_variants$symbol[which(gene_variants$labId == i)]))
  mutations[i, jj] <- 1
}

# Add Translocations 
# look for translocations:  filtering the rows to only those that contain "AML with inv" or "AML with t" 
translocations_clinical<-clinical[,c("LabId","specificDxAtAcquisition")][grep("AML with inv", clinical$specificDxAtAcquisition),]
translocations_clinical<-rbind(translocations_clinical,clinical[,c("LabId","specificDxAtAcquisition")][grep("AML with t", clinical$specificDxAtAcquisition),])

trans<-as.character(unique(translocations_clinical$specificDxAtAcquisition))
add<-matrix(0, nrow=nrow(mutations), ncol=length(trans)) # add columns for translocations
colnames(add)<-trans
rownames(add)<-rownames(mutations)

for(j in colnames(add)){
  p<-translocations_clinical$LabId[which(translocations_clinical$specificDxAtAcquisition ==j)]
  p<-p[which(p %in% rownames(add))]
  add[p,j]<-1
}

# Change column name and also add t: extract the second element of each translocation name (i.e., the specific type of translocation) and use it as the new column name in add. 
trans<-sapply(trans, function(X){unlist(strsplit(X, split = "; "))[2]})
colnames(add)<-trans
mutations<-cbind(mutations,add)

# Unifying imputed Drug and Mutation matrices with common patients
drug_matrix_2<-as.matrix(t(drug_matrix)) # transpose drug 
mutations<-as.matrix(mutations)

# impute missing values DUDA, PORQUE LO HACES ASI: imputar los valores, de todas las parejas drug-cell line, no todas tienen IC50, imputar las que no tienen
# imput knn
if (sum(is.na(drug_matrix_2)) > 0) {
  drug_matrix_3 <- drug_matrix_2[
    -which((rowSums(is.na(drug_matrix_2)) / dim(drug_matrix_2)[2]) > 0.8),
    -which((colSums(is.na(drug_matrix_2)) / dim(drug_matrix_2)[1]) > 0.7)
  ] # el porcentaje de nas por filas tiene que ser menor de un 20 y en col del 30
  Out <- impute.knn(as.matrix(drug_matrix_3),
    k = 10, rowmax = 0.8,
    colmax = 0.75, maxp = 1500, rng.seed = 362436069
  )
  drug_matrix_4 <- Out$data[, -which(colnames(Out$data) %in% c("Elesclomol", "JNJ-7706621"))] # Los quitamos por ser toxicos, no los queremos, cuando hicimos la normalizacion se deberian de haber ido solitos
} else {
  drug_matrix_4 <- drug_matrix_2[, -which(colnames(drug_matrix_2) %in% c("Elesclomol", "JNJ-7706621"))]
}


# Combining Drug and mutational Data: quedarme con los pacientes que tienen respuesta a farmaco
patients<-rownames(mutations)[which(rownames(mutations) %in% rownames(drug_matrix_4))]

# filtrar mutaciones
mutations<-mutations[patients, colSums(mutations[patients, ])>0.01*length(patients)] # EL 0.01, al menos un 1% de los pacientes tienen que tener la mutacion
sum(colSums(mutations)<1)

# calculo del IC50* 
drug<-drug_matrix_4[patients,]
drug<-log10(drug)-1
drug <- t(t(drug) - colMeans(drug, na.rm=T))  # unificar los farmacos por columnas, compensamos con la dosis

drug_response_matrix_MOM<-drug-max(drug,na.rm = T)  # MOM minimiza en vez de maximizar 
drug_response_matrix <- drug - min(drug) # por como funciona MOM los valores tiene que ser todos 


# Tidying up the variables
# borrar variables
rm(list = c("drug_matrix_2", "drug_matrix_3", "drug_matrix_4", "drug_response", "gene_variants"))
rm(list = c("Out", "translocations_clinical", "clinical", "add", "drug_matrix")) 

pseudop <- exp(-.5*drug_response_matrix) # Convert weights into pseudoprobabilities for Multinomial
pseudop <- pseudop /rowSums(pseudop) # Probabilities sum up one

#############################################################
# Mutation models
#############################################################

# Comparison of methods with Cross-Validation --------------------------------------------
"
We performed a 5-fold cross-validation using the BeatAML dataset.
We trained all models with genetic variants data from 319 patients, dividing the cohort between the training samples 4-folds and testing samples the selected 1-fold.
Each of the folds were tested, and the predicted IC50* for the 5-fold testing was compared for all the methods and compared against the Oracle -the drug with the optimum IC50\*.

We calculated the Oracle as the minimum IC50* value for each patient.
"

#Results from MOM and KRL were directly loaded into the workspace due to the computing limitations in R.
set.seed(2022)
folds <- 5
groups <- sample(1:folds, nrow(drug_response_matrix), replace = T)
treatmentMnLassoMut <- treatmentLassoMut <- treatmentODTMut <- treatmentODTMutSqrt <- rep(NA, length(groups))

# predict functions require the drug_response_matrix just to extract the names of the drugs
# train and test
for (group in 1:folds) {
  cat("Fold: ", group, "\n")
  dejar <- groups != group
  quitar <- !dejar
  
  # ODTMut
  cat("ODTMut\n")
  ODTMut <- trainTreeMut(mutations[dejar,], drug_response_matrix[dejar,], minbucket = 10)
  treatmentODTMut[quitar] <- predictTreeMut(ODTMut,drug_response_matrix[quitar,], mutations[quitar,])

  # ODTMut sqrt
  cat("ODTMut sqrt\n")
  ODTMutSqrt <- trainTreeMut(mutations[dejar,], (drug_response_matrix[dejar,])^.5, minbucket = 10)
  treatmentODTMutSqrt[quitar] <- predictTreeMut(ODTMutSqrt,drug_response_matrix[quitar,], mutations[quitar,])

  # Using MnLasso
  cat("Multinomial Lasso Mut\n")
  MnLassoMut <- trainMnLasso(mutations[dejar,], pseudop[dejar,])
  treatmentMnLassoMut[quitar] <- predictMnLasso(MnLassoMut,drug_response_matrix[quitar,],mutations[quitar,])

  # Using Lasso
  cat("Lasso Mut\n")
  LassoMut <- trainLasso(mutations[dejar,], drug_response_matrix[dejar,])
  treatmentLassoMut[quitar] <- predictLasso(LassoMut,drug_response_matrix[quitar,],mutations[quitar,])

}
# Use same groups as previously voy aqui
treatmentBOSOMut <- rep(NA, length(groups))
for (group in 1:folds) {
  cat("Fold: ", group, "\n")
  dejar <- groups != group
  quitar <- !dejar
  BOSOMut <- trainBOSO(mutations[dejar,], drug_response_matrix[dejar,],maxVarsBlock = 10, standardize=F, IC = "eBIC")
  treatmentBOSOMut[quitar] <- predictBOSO(BOSOMut, drug_response_matrix[quitar,], mutations[quitar,])
}

save(treatmentODTMut, treatmentODTMutSqrt, treatmentMnLassoMut, treatmentLassoMut, treatmentBOSOMut, file = "C:/Users/ksada/OneDrive - Tecnun/Paper XAI Methods/Rdata/mutations1.RData")
load("C:/Users/ksada/OneDrive - Tecnun/Paper XAI Methods/Rdata/mutations1.RData")

load("./data/output/mom_COMPARISON.RData")

# KRL
KRL_results_cv <- read_csv("data/output/KRL_results_5cv_samples_prediction.csv")
KRL_all<-read_csv("data/output/KRL_results_all_samples_prediction.csv")

KRL_results_cv<-KRL_results_cv[which(KRL_results_cv$Patients %in% rownames(mutations)),]
KRL_all<-KRL_all[which(KRL_all$Patients %in% rownames(mutations)),]

KRL_results_cv<-KRL_results_cv[which(KRL_results_cv$DrugName %in% colnames(drug)),] # nolint: line_length_linter.
KRL_all<-KRL_all[which(KRL_all$DrugName %in% colnames(drug)),]

KRL_all$`Best Drug`<-KRL_all$`Best Drug`+1
KRL_results_cv$Drug<-KRL_results_cv$Drug+1

KRL_all$IC50<-drug_response_matrix[cbind(KRL_all$Patients, KRL_all$DrugName)]
KRL_results_cv$IC50<-drug_response_matrix[cbind(KRL_results_cv$Patients, KRL_results_cv$DrugName)]

treatmentOracle <- apply(drug_response_matrix, 1, which.min)
 
### Plotting BeatAML Cross-validation Results and comparing means difference
treatment_plot<-data.frame(Method="ORACLE", IC50=drug_response_matrix[cbind(1:nrow(drug_response_matrix),treatmentOracle)])
treatment_plot<-rbind(treatment_plot,data.frame(Method="ODTMut", IC50=drug_response_matrix[cbind(1:nrow(drug_response_matrix), treatmentODTMut)]))
treatment_plot<-rbind(treatment_plot,data.frame(Method="ODTMutSqrt", IC50=drug_response_matrix[cbind(1:nrow(drug_response_matrix), treatmentODTMutSqrt)]))
treatment_plot<-rbind(treatment_plot,data.frame(Method="MnLassoMut", IC50=drug_response_matrix[cbind(1:nrow(drug_response_matrix), treatmentMnLassoMut)]))
treatment_plot<-rbind(treatment_plot,data.frame(Method="LassoMut", IC50=drug_response_matrix[cbind(1:nrow(drug_response_matrix), treatmentLassoMut)]))
treatment_plot<-rbind(treatment_plot,data.frame(Method="BOSOMut", IC50=drug_response_matrix[cbind(1:nrow(drug_response_matrix), treatmentBOSOMut)]))
treatment_plot<-rbind(treatment_plot,data.frame(Method="MOM", IC50=MOM$IC50))
treatment_plot<-rbind(treatment_plot,data.frame(Method="KRL", IC50=KRL_results_cv$IC50))

treatment_plot$Method<-factor(treatment_plot$Method, levels=c("ORACLE", "MOM", "BOSOMut", "ODTMutSqrt", "LassoMut", "ODTMut", "MnLassoMut","KRL"))

ggplot(treatment_plot, aes(x=Method, y=IC50, fill=Method)) +
geom_violin() +
theme_bw()+scale_fill_npg()+ylab("IC50*") +
ggtitle("BeatAML 5-fold Cross Validation")+
theme(text = element_text(size = 30, family = "Helvetica"))

treatment_plot_gg <- ggplot(treatment_plot, aes(x=Method, y=IC50, fill=Method)) +
                        geom_boxplot(lwd=0.8,position=position_dodge(0.8)) +
                        theme_bw()+scale_fill_npg()+ylab("IC50*") +
                        theme(text = element_text(size = 20, family = "Roboto"),
                        axis.text.x = element_text(size = 11),
                        plot.background = element_rect(fill = "transparent", size = 0, color = "black"),
                        panel.background = element_rect(fill = "transparent"))
treatment_plot_gg
ggsave("treatment_plot_mutations.png", treatment_plot_gg, width = 10, height = 6, dpi = 1000)

# compute wilcoxon tests with the Oracle
wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="MOM"], alternative="less")
wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="BOSOMut"], alternative="less")
wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="ODTMutSqrt"], alternative="less")
wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="LassoMut"], alternative="less")
wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="ODTMut"], alternative="less")
wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="MnLassoMut"], alternative="less")
wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="KRL"], alternative="less")

# Independent Cohort Validation --------------------------------------------
"
One of the main challenges of Machine Learning, including Precision Medicine, is generalization, i.e. the ability to adapt to new, previously unseen data.
All the methods were tested on the GDSC AML dataset to check their generalization ability.
The models were trained using the BeatAML dataset and were used to predict the optimal drug for AML cell lines from GDSC using its mutation files.
Each of the cell lines was recommended a drug, we compared the all-samples IC50 for all the models and against the Oracle (the drug with the minimum IC50 for each cell line).
"

# Training models in complete BeatAML data and timing
"
We trained the models in BeatAML cohort and measured the training time for each model except KRL and MOM, that were meassured independently.
Technical implementation refers to the computational burden and software that the method requires.Despite it could be considered less important, some of the algorithms require hours of computing time for the BeatAML of the subset of AML samples in GDSC -that be considered to be small/medium size.
" 
# Train with ALL the data
# ODTMut
cat("ODTMut\n")
tic()
ODTMut_All <- trainTreeMut(mutations, drug_response_matrix, minbucket = 10)
treatmentODTMut_ALL <- predictTreeMut(ODTMut_All,drug_response_matrix, mutations)
toc()

# ODTMut sqrt
cat("ODTMut sqrt\n")
tic()
ODTMutSqrt_All <- trainTreeMut(mutations, (drug_response_matrix)^.5, minbucket = 10)
treatmentODTMutSqrt_All <- predictTreeMut(ODTMutSqrt_All,drug_response_matrix, mutations)
toc()

# Using MnLasso
cat("Multinomial Lasso Mut\n")
tic()
MnLassoMut_All <- trainMnLasso(mutations, pseudop)
treatmentMnLassoMut_All <- predictMnLasso(MnLassoMut_All,drug_response_matrix,mutations)
toc()

# Using Lasso
cat("Lasso Mut\n")
tic()
LassoMut_All <- trainLasso(mutations, drug_response_matrix)
treatmentLassoMut_All <- predictLasso(LassoMut_All,drug_response_matrix,mutations)
toc()

# Using Bosso
cat("Bosso Mut\n")
tic()
BOSOMut_All <- trainBOSO(mutations, drug_response_matrix,maxVarsBlock = 10, standardize=F, IC = "eBIC")
treatmentBOSOMut_All <- predictBOSO(BOSOMut_All,drug_response_matrix,mutations)
toc()

### Validation in GDSC
"
The models were trained using the BeatAML dataset and were used to predict the optimal drug for AML cell lines from GDSC using its mutation files. 
Each of the cell lines was recommended a drug, we compared the all-samples IC50 for all the models and against the Oracle (the drug with the minimum IC50 for each cell line).
The function `predictMOM_ALL`.
Applies the decision tree generated by MOM in AML.

Expresion medida en GSDC puede que no se compare con la expresion de BeatAML, microArrays 
"
#### Load GDSC and CCLE data
GDSC1_ESS<-read.csv("./data/input/LAML_IC_Wed Nov 24 10_31_42 2021.csv")
GDSC1_MUT<-read.csv("./data/input/LAML_Genetic_feature_variant_Wed Nov 24 GDSC1.csv")

load("./data/input/DEMETER2/Processed/new_CCLE_DEMETER2_20Q3_SYNLETHDB_data_v1.Rdata")
CCLE_mutation_data <- demeter2_mut
# cambiar manualmente porque en el cima hicieron validacion y si que estaban mutadas 
CCLE_mutation_data[c("MV411_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","NOMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","OCIAML2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"),"TP53"]   <- 1
CCLE_mutation_data["HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","KDM6A"]  <- 1
CCLE_mutation_data[c("MOLM13_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","MV411_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"),"FLT3"]   <- 1
CCLE_mutation_data["NOMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","KRAS"]   <- 1

rm(demeter2_mut, demeter2_scores, demeter2_fus)

# Biomarker Pre-processing
# We build the mutations matrix and complete that information with CCLE mutational data.

# Biomarker pre-processing and adding CCLE information
aux_gdsc<-data.frame(cell_line_name="ME-1",
                     cosmic_sample_id="1330942",
                     gdsc_desc1="blood",
                     gdsc_desc2="acute_myeloid_leukaemia",
                     tcga_desc="LAML",
                     genetic_feature="CBFB-MYH11_mut",
                     is_mutated="1",
                     recurrent_gain_loss=NA,
                     genes_in_segment=NA)

BMVal<-rbind(GDSC1_MUT,aux_gdsc)

BMVal$cell_line_name<-sapply(BMVal$cell_line_name, FUN=function(X){return(paste(unlist(strsplit(X, split="-")), collapse=""))})

BMVal$cell_line_name<-paste0(BMVal$cell_line_name, "_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE")


BM_matrix<-matrix(0,
                  nrow=length(unique(BMVal$cell_line_name)), 
                  ncol=ncol(mutations))

colnames(BM_matrix)<-colnames(mutations)
rownames(BM_matrix)<-as.character(unique(BMVal$cell_line_name))

for(k in 1:nrow(BM_matrix)){
  c<-rownames(BM_matrix)[k]
  for(m in 1:ncol(BM_matrix)){
    bm_ccle<-colnames(BM_matrix)[m]
    bm<-paste0(bm_ccle, "_mut")
    aux<-BMVal$is_mutated[which(BMVal$cell_line_name==c & BMVal$genetic_feature==bm)]
    
    if((bm_ccle %in% colnames(CCLE_mutation_data))&& (c %in% rownames(CCLE_mutation_data))){
      aux_ccle<-CCLE_mutation_data[c, bm_ccle]
    }else{
      aux_ccle<-NULL
    }
    if(length(aux_ccle)<1 && length(aux)<1){next}
    if(max(aux,aux_ccle)==1){
      BM_matrix[k,m]<-as.numeric(max(aux,aux_ccle))
    }
    aux<-NULL
    aux_ccle<-NULL
  }
}

BM_matrix<-BM_matrix[,match( colnames(BM_matrix),colnames(mutations))]
# colnames(BM_matrix)<-sapply(colnames(BM_matrix), FUN=function(X){return(unlist(strsplit(X, split="_"))[1])})

mutations_AML_GDSC<-BM_matrix

#### Drug information pre-processing
GDSC1_ESS$Cell.line.name<-sapply(GDSC1_ESS$Cell.line.name, FUN=function(X){return(paste(unlist(strsplit(X, split="-")), collapse=""))})
GDSC1_ESS$Cell.line.name<-paste0(GDSC1_ESS$Cell.line.name, "_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE")

Ess_matrix<-matrix(50000,
                   nrow=length(unique(GDSC1_ESS$Cell.line.name)), 
                   ncol=length(unique(GDSC1_ESS$Drug.name)))


colnames(Ess_matrix)<-as.character(unique(GDSC1_ESS$Drug.name))
rownames(Ess_matrix)<-as.character(unique(GDSC1_ESS$Cell.line.name))

for(k in 1:nrow(Ess_matrix)){
  c<-rownames(Ess_matrix)[k]
  for(m in 1:ncol(Ess_matrix)){
    d<-colnames(Ess_matrix)[m]
    aux<-GDSC1_ESS$IC50[which(GDSC1_ESS$Cell.line.name==c & GDSC1_ESS$Drug.name==d)]
    aux<-mean(as.numeric(aux))
    if(length(aux)<1){next}
    Ess_matrix[k,m]<-aux
    aux<-NULL
  }
}

Ess_matrix[Ess_matrix=="NaN"]<-NA
# Impute Missing Data for Validation
Ess_matrix<-Ess_matrix[-which((rowSums(is.na(Ess_matrix))/dim(Ess_matrix)[2])>0.8),
                       -which((colSums2(is.na(Ess_matrix))/dim(Ess_matrix)[1])>0.70)]


Out_gdsc<-impute.knn(as.matrix(Ess_matrix), k = 10, rowmax = 0.75, 
                     colmax = 0.8, maxp = 1500, rng.seed=362436069)

colnames(Out_gdsc$data)[which(colnames(Out_gdsc$data)=="Venotoclax")]<-"Venetoclax"

GDSC_AML<-Out_gdsc$data
GDSC_AML<-t(t(GDSC_AML)-colMeans(GDSC_AML))

drugsGDSC<-sapply(colnames(drug_response_matrix), FUN=function(X){return(unlist(strsplit(X, split=" "))[1])})

GDSC_AML<-GDSC_AML[,na.omit(match(drugsGDSC, colnames(GDSC_AML)))]

colnamesGDSC<-colnames(drug_response_matrix)[match(colnames(GDSC_AML),drugsGDSC)]

colnames(GDSC_AML)<-colnamesGDSC

GDSC_AML<-GDSC_AML[which(rownames(GDSC_AML) %in% rownames(mutations_AML_GDSC)),]
GDSC_AML<-GDSC_AML[match(rownames(mutations_AML_GDSC), rownames(GDSC_AML)),]

identical(rownames(GDSC_AML), rownames(mutations_AML_GDSC))

# Make Predictions with the models

# We used the models trained in BeatAML to predict over GDSC.
 
identical(colnames(mutations_AML_GDSC), colnames(mutations)) #TRUE

Guideline_ODTMutSqrt<-predictTreeMut(ODTMutSqrt_All,drug_response_matrix, mutations_AML_GDSC)
Guideline_ODTMutSqrt<-colnames(drug_response_matrix)[Guideline_ODTMutSqrt]

Guideline_BOSO<-predictBOSO(BOSOMut_All,drug_response_matrix,mutations_AML_GDSC)
Guideline_BOSO<-colnames(drug_response_matrix)[Guideline_BOSO]

Guideline_ODTMut<-predictTreeMut(ODTMut_All,drug_response_matrix, mutations_AML_GDSC)
Guideline_ODTMut<-colnames(drug_response_matrix)[Guideline_ODTMut]

Guideline_Multinomial<-predictMnLasso(MnLassoMut_All,drug_response_matrix, mutations_AML_GDSC)
Guideline_Multinomial<-colnames(drug_response_matrix)[Guideline_Multinomial]

Guideline_Lasso<-predictLasso(LassoMut_All,drug_response_matrix, mutations_AML_GDSC)
Guideline_Lasso<-colnames(drug_response_matrix)[Guideline_Lasso]

Guideline_MOM<-predictMOM_ALL(mutations_AML_GDSC)
Guideline_MOM<-colnames(drug_response_matrix)[Guideline_MOM$Treatment]

Guideline_KRL<-read_csv("data/output/KRL_GDSC_pred.csv")
Guideline_KRL$DrugName[which(!(Guideline_KRL$DrugName %in% colnames(drug_response_matrix)))]<-NA
Guideline_KRL<-Guideline_KRL$DrugName

#### Obtain Prediction IC50 values

#We will now obtain the GDSC IC50 values corresponding to the predicted treatments
 
# Validate the Guidelines..................................................

identical(rownames(GDSC_AML), rownames(mutations_AML_GDSC)) # TRUE

GDSC_IC50<-data.frame(CL=rownames(mutations_AML_GDSC),
                      ODTMut=NA,
                      ODTMutSqrt=NA,
                      MnLassoMut=NA,
                      BOSOMut=NA,
                      LassoMut=NA,
                      MOM=NA,
                      KRL=NA,
                      ORACLE=NA)

for(i in 1:nrow(GDSC_IC50)){
  CL<-GDSC_IC50$CL[i]
  
  # ODTMut
  dODTMut<-Guideline_ODTMut[i]
  if(dODTMut %in% colnames(GDSC_AML)){
    GDSC_IC50$ODTMut[i]<-GDSC_AML[CL,dODTMut]
  }
  
  # ODTMutSqrt
  dODTMutSQrt<-Guideline_ODTMutSqrt[i]
  GDSC_IC50$ODTMutSqrt[i]<-GDSC_AML[CL,dODTMutSQrt]
  
  # Multinomial
  dMultinomial<-Guideline_Multinomial[i]
  if(dMultinomial %in% colnames(GDSC_AML)){
    GDSC_IC50$MnLassoMut[i]<-GDSC_AML[CL,dMultinomial]
  }
  
  # BOSO
  dBOSO<-Guideline_BOSO[i]
  if(dBOSO %in% colnames(GDSC_AML)){
    GDSC_IC50$BOSOMut[i]<-GDSC_AML[CL,dBOSO]
  }
  
  # Lasso
  dLasso<-Guideline_Lasso[i]
  if(dLasso %in% colnames(GDSC_AML)){
    GDSC_IC50$LassoMut[i]<-GDSC_AML[CL,dLasso]
  }
  
  #MOM
  dMOM<-Guideline_MOM[i]
  GDSC_IC50$MOM[i]<-GDSC_AML[CL,dMOM]

  #KRL
  dKRL<-Guideline_KRL[i]
  if(dKRL %in% colnames(GDSC_AML)){
    GDSC_IC50$KRL[i]<-GDSC_AML[CL,dKRL]
  }
  
  #ORACLE
  GDSC_IC50$ORACLE[i]<-min(GDSC_AML[CL,])
}

rownames(GDSC_IC50)<-GDSC_IC50$CL

#### Plotting GDSC Validation MALLLLLLLL
Validation_plot<-data.frame(Dataset=rep("GDSC", nrow(GDSC_IC50)*8),  #8 es el numero de metodos ncol(GDSC_IC50)-1
                                                   CellLine=rep(rownames(GDSC_IC50),8),
                                                   Method=rep(colnames(GDSC_IC50[,-1]), each=nrow(GDSC_IC50)), 
                                                   IC50=unlist(GDSC_IC50[,-1]))

rownames(Validation_plot)<-1:nrow(Validation_plot)

mycomparisons<-combn(unique(Validation_plot$Method),2)
mycomparisons<-as.list(data.frame(mycomparisons))
mycomparisons<-mycomparisons[-grep("ORACLE",mycomparisons)]

Plotting<-Validation_plot[Validation_plot$Dataset=="GDSC",]

library(ggsci)

Plotting$Method<-factor(Plotting$Method, levels=c("ORACLE", "MOM", "BOSOMut", "ODTMutSqrt", "LassoMut", "ODTMut", "MnLassoMut", "KRL"))

g1_gdsc<-ggplot(Plotting, aes(x=Method, y=IC50, fill=Method))+geom_violin()+
                        theme_bw()+scale_fill_npg()+ylab("GDSC IC50")+ggtitle("GDSC") +
                        theme(text = element_text(size = 20, family = "Roboto"),
                        axis.text.x = element_text(size = 11),
                        plot.background = element_rect(fill = "transparent", size = 0, color = "black"),
                        panel.background = element_rect(fill = "transparent"))

g2_gdsc<-ggplot(Plotting, aes(x=Method, y=IC50, fill=Method))+
                        geom_boxplot(lwd=0.8,position=position_dodge(0.8))+
                        theme_bw()+scale_fill_npg()+ylab("GDSC IC50")+ggtitle("GDSC")+
                        theme(text = element_text(size = 20, family = "Roboto"),
                        axis.text.x = element_text(size = 11),
                        plot.background = element_rect(fill = "transparent", size = 0, color = "black"),
                        panel.background = element_rect(fill = "transparent"))

g1_gdsc
g2_gdsc

ggsave("treatment_plot_GDSC.png", g2_gdsc, width = 10, height = 6, dpi = 1000)

# Intragroup Validation  -------------------------------------------- -- EL NOMBRE DE INITRAGROUP ME LIA, NO ME GUSTA NADA
"
We compared whether the IC50* of a drug in patients in whom it was recommended was lower than the IC50* in patients in whom it was not recommended.
Using this information we compared the sensitivity to a drug for a specific group against the sensitivity to that drug for the rest of the samples by using a 2-tailed Wilcoxon test.
This analysis was performed both for the BeatAML dataset (training dataset) and the GDSC AML cell lines cohort (predicted dataset).
"
### Obtaining MOM and KRL predicitions in all samples 
treatmentMOM_All<-predictMOM_ALL(mutations)
identical(treatmentMOM_All$CLs, rownames(drug_response_matrix)) # TRUE

treatmentKRL_All<- read_csv("data/output/KRL_results_all_samples_prediction.csv", show_col_types = FALSE)
treatmentKRL_All<-treatmentKRL_All[which(treatmentKRL_All$Patients %in% rownames(mutations)),]
treatmentKRL_All$DrugName[which(!(treatmentKRL_All$DrugName %in% colnames(drug_response_matrix)))]<-NA

### Intragroup Validation in BeatAML
 
# 1. BeatAML--------------------------------------------------------------------
 # MAGIC NUMBERS AGAIN!!! 7? (number of methods)-1 NO HAY ORACLE
Validation_Beat<-data.frame(patients=rep(rownames(drug_response_matrix),nlevels(treatment_plot$Method)-1), 
                            drug=c(colnames(drug_response_matrix)[treatmentODTMut_ALL],
                                   colnames(drug_response_matrix)[treatmentODTMutSqrt_All],
                                   colnames(drug_response_matrix)[treatmentMnLassoMut_All],
                                   colnames(drug_response_matrix)[treatmentLassoMut_All],
                                   colnames(drug_response_matrix)[treatmentBOSOMut_All],
                                   colnames(drug_response_matrix)[treatmentMOM_All$Treatment],
                                   treatmentKRL_All$DrugName), 
                            IC50=c(drug_response_matrix[cbind(1:nrow(drug_response_matrix), treatmentODTMut_ALL)],
                                   drug_response_matrix[cbind(1:nrow(drug_response_matrix), treatmentODTMutSqrt_All)],
                                   drug_response_matrix[cbind(1:nrow(drug_response_matrix), treatmentMnLassoMut_All)],
                                   drug_response_matrix[cbind(1:nrow(drug_response_matrix), treatmentLassoMut_All)],
                                   drug_response_matrix[cbind(1:nrow(drug_response_matrix), treatmentBOSOMut_All)],
                                   drug_response_matrix[cbind(1:nrow(drug_response_matrix), treatmentMOM_All$treatment)],
                                   drug_response_matrix[cbind(treatmentKRL_All$Patients,treatmentKRL_All$DrugName )]),
                            Method=rep(unique(treatment_plot$Method[!(treatment_plot$Method=="ORACLE")]), each=nrow(drug_response_matrix)))

#### MOM ME SALE PEOR
Plot_aux<-Validation_Beat[Validation_Beat$Method=="MOM",] # matriz de mutaciones que sea distinta tiene que tener 64

MOM_plots<-NULL

for(j in 1:length(unique(Plot_aux$drug))){
  d<-unique(Plot_aux$drug)[j]
  ind<-Plot_aux$drug==d
  plotDrug<-data.frame(Drug=j,
                       drugname=d,
                       Patients=Plot_aux$patients,
                       IC50=NA,
                       Recommended="No")
  
  plotDrug$Recommended[ind]<-"Yes"
  plotDrug$IC50[ind]<-Plot_aux$IC50[ind]
  plotDrug$IC50[!ind]<-drug_response_matrix[!ind, d]
  plotDrug$Recommended[!ind]<-"NO"
  
  MOM_plots<-rbind(MOM_plots,plotDrug)
  
}

gg1_validation_beat<-ggplot(MOM_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+
  scale_fill_npg()+stat_compare_means(vjust = 1)+ylab("IC50*")+ggtitle("MOM")

gg1_validation_beat
 

#### ODTMut

# TREE MUT
Plot_aux<-Validation_Beat[Validation_Beat$Method=="ODTMut",]

ODT_plots<-NULL

for(j in 1:length(unique(Plot_aux$drug))){
  d<-unique(Plot_aux$drug)[j]
  ind<-Plot_aux$drug==d
  plotDrug<-data.frame(Drug=j,
                       drugname=d,
                       Patients=Plot_aux$patients,
                       IC50=NA,
                       Recommended="No")
  
  plotDrug$Recommended[ind]<-"Yes"
  plotDrug$IC50[ind]<-Plot_aux$IC50[ind]
  plotDrug$IC50[!ind]<-drug_response_matrix[!ind, d]
  plotDrug$Recommended[!ind]<-"NO"
  
  ODT_plots<-rbind(ODT_plots,plotDrug)
  
}
gg3_validation_beat<-ggplot(ODT_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_npg()+
  stat_compare_means()+ylab("IC50*")+ggtitle("ODTMut")
gg3_validation_beat

#### ODTMutSqrt

# TREE MUT SQRT

Plot_aux<-Validation_Beat[Validation_Beat$Method=="ODTMutSqrt",]

ODTQR_plots<-NULL

for(j in 1:length(unique(Plot_aux$drug))){
  d<-unique(Plot_aux$drug)[j]
  ind<-Plot_aux$drug==d
  plotDrug<-data.frame(Drug=j,
                       drugname=d,
                       Patients=Plot_aux$patients,
                       IC50=NA,
                       Recommended="No")
  
  plotDrug$Recommended[ind]<-"Yes"
  plotDrug$IC50[ind]<-Plot_aux$IC50[ind]
  plotDrug$IC50[!ind]<-drug_response_matrix[!ind, d]
  plotDrug$Recommended[!ind]<-"NO"
  
  ODTQR_plots<-rbind(ODTQR_plots,plotDrug)
  
}

gg4_validation_beat<-ggplot(ODTQR_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_npg()+
  stat_compare_means()+ylab("IC50*")+ggtitle("ODTMutSqrt")

gg4_validation_beat
 
#### BOSO

#BOSO

Plot_aux<-Validation_Beat[Validation_Beat$Method=="BOSOMut",]

BOSO_plots<-NULL

for(j in 1:length(unique(Plot_aux$drug))){
  d<-unique(Plot_aux$drug)[j]
  ind<-Plot_aux$drug==d
  plotDrug<-data.frame(Drug=j,
                       drugname=d,
                       Patients=Plot_aux$patients,
                       IC50=NA,
                       Recommended="No")
  
  plotDrug$Recommended[ind]<-"Yes"
  plotDrug$IC50[ind]<-Plot_aux$IC50[ind]
  plotDrug$IC50[!ind]<-drug_response_matrix[!ind, d]
  plotDrug$Recommended[!ind]<-"NO"
  
  BOSO_plots<-rbind(BOSO_plots,plotDrug)
  
}

gg2_validation_beat<-ggplot(BOSO_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_npg()+
  stat_compare_means(vjust = 1)+ylab("IC50*")+ggtitle("BOSO")

gg2_validation_beat

#### Multinomial
 
# MULTINOMIAL

Plot_aux<-Validation_Beat[Validation_Beat$Method=="MnLassoMut",]

MnLasso_plots<-NULL

for(j in 1:length(unique(Plot_aux$drug))){
  d<-unique(Plot_aux$drug)[j]
  ind<-Plot_aux$drug==d
  plotDrug<-data.frame(Drug=j,
                       drugname=d,
                       Patients=Plot_aux$patients,
                       IC50=NA,
                       Recommended="No")
  
  plotDrug$Recommended[ind]<-"Yes"
  plotDrug$IC50[ind]<-Plot_aux$IC50[ind]
  plotDrug$IC50[!ind]<-drug_response_matrix[!ind, d]
  plotDrug$Recommended[!ind]<-"NO"
  
  MnLasso_plots<-rbind(MnLasso_plots,plotDrug)
  
}

gg5_validation_beat<-ggplot(MnLasso_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_npg()+
  stat_compare_means(vjust=1)+ylab("IC50*")+ggtitle("Multinomial")
gg5_validation_beat
 
#### KRL

# KRL

Plot_aux<-Validation_Beat[Validation_Beat$Method=="KRL",]

KRL_plots<-NULL

for(j in 1:length(unique(Plot_aux$drug))){
  d<-unique(Plot_aux$drug)[j]
  if(is.na(d)){next}
  ind<-na.omit(Plot_aux$drug==d)
  plotDrug<-data.frame(Drug=j,
                       drugname=d,
                       Patients=Plot_aux$patients,
                       IC50=NA,
                       Recommended="No")
  
  plotDrug$Recommended[ind]<-"Yes"
  plotDrug$IC50[ind]<-Plot_aux$IC50[ind]
  plotDrug$IC50[!ind]<-drug_response_matrix[!ind, d]
  plotDrug$Recommended[!ind]<-"NO"
  
 KRL_plots<-rbind(KRL_plots,plotDrug)
  
}

gg7_validation_beat<-ggplot(KRL_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_npg()+
  stat_compare_means(vjust=1)+ylab("IC50*")+ggtitle("KRL")

gg7_validation_beat

# Intragroup Validation in GDSC

identical(colnames(mutations_AML_GDSC), colnames(mutations)) #TRUE


Guideline_ODTMutSqrt<-predictTreeMut(ODTMutSqrt_All,drug_response_matrix, mutations_AML_GDSC)
Guideline_ODTMutSqrt<-colnames(drug_response_matrix)[Guideline_ODTMutSqrt]


Guideline_BOSO<-predictBOSO(BOSOMut_All,drug_response_matrix,mutations_AML_GDSC)
Guideline_BOSO<-colnames(drug_response_matrix)[Guideline_BOSO]

Guideline_ODTMut<-predictTreeMut(ODTMut_All,drug_response_matrix, mutations_AML_GDSC)
Guideline_ODTMut<-colnames(drug_response_matrix)[Guideline_ODTMut]

Guideline_Multinomial<-predictMnLasso(MnLassoMut_All,drug_response_matrix, mutations_AML_GDSC)
Guideline_Multinomial<-colnames(drug_response_matrix)[Guideline_Multinomial]

Guideline_Lasso<-predictLasso(LassoMut_All,drug_response_matrix, mutations_AML_GDSC)
Guideline_Lasso<-colnames(drug_response_matrix)[Guideline_Lasso]

Guideline_MOM<-predictMOM_ALL(mutations_AML_GDSC)
Guideline_MOM<-colnames(drug_response_matrix)[Guideline_MOM$Treatment]

Guideline_KRL<-read_csv("data/output/KRL_GDSC_pred.csv")
Guideline_KRL$DrugName[which(!(Guideline_KRL$DrugName %in% colnames(drug_response_matrix)))]<-NA
Guideline_KRL<-Guideline_KRL$DrugName

Validation_GDSC<-data.frame(patients=rep(rownames(GDSC_AML),7),
                            drug=c(Guideline_ODTMut,
                                   Guideline_ODTMutSqrt,
                                   Guideline_Multinomial,
                                   Guideline_Lasso,
                                   Guideline_BOSO,
                                   Guideline_MOM, 
                                   Guideline_KRL), 
                            IC50=NA,
                            Method=rep(unique(treatment_plot$Method[!(treatment_plot$Method=="ORACLE")]), each=nrow(GDSC_AML)))

#### MOM

# MOM

Plot_aux<-Validation_GDSC[Validation_GDSC$Method=="MOM",]

MOM_plots<-NULL

for(j in 1:length(unique(Plot_aux$drug))){
  d<-unique(Plot_aux$drug)[j]
  ind<-Plot_aux$drug==d
  plotDrug<-data.frame(Drug=j,
                       drugname=d,
                       Patients=Plot_aux$patients,
                       IC50=NA,
                       Recommended="No")
  
  plotDrug$Recommended[ind]<-"Yes"
  plotDrug$IC50[ind]<-GDSC_AML[ind,d]
  plotDrug$IC50[!ind]<-GDSC_AML[!ind, d]
  plotDrug$Recommended[!ind]<-"No"
  
  MOM_plots<-rbind(MOM_plots,plotDrug)
  
}

MOM_plots$Recommended<-factor(MOM_plots$Recommended, levels=c("No", "Yes"))

gg1_validation_gdsc<-ggplot(MOM_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_manual(values=c("#999999", "#E69F00"))+
  stat_compare_means()+ylab("IC50*")+ggtitle("MOM")

gg1_validation_gdsc

#### ODTMut
 
# TREE MUT

Plot_aux<-Validation_GDSC[Validation_GDSC$Method=="ODTMut",]

ODT_plots<-NULL

for(j in 1:length(unique(Plot_aux$drug))){
  d<-unique(Plot_aux$drug)[j]
  ind<-Plot_aux$drug==d
  plotDrug<-data.frame(Drug=j,
                       drugname=d,
                       Patients=Plot_aux$patients,
                       IC50=NA,
                       Recommended="No")
  
  plotDrug$Recommended[ind]<-"Yes"
  if(d %in% colnames(GDSC_AML)){
    plotDrug$IC50[ind]<-GDSC_AML[ind,d]
    plotDrug$IC50[!ind]<-GDSC_AML[!ind, d]
    plotDrug$Recommended[!ind]<-"No"
  }
  
  ODT_plots<-rbind(ODT_plots,plotDrug)
  
}

gg3_validation_gdsc<-ggplot(ODT_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_manual(values=c("#999999", "#E69F00"))+
  stat_compare_means()+ylab("IC50*")+ggtitle("ODTMut")

gg3_validation_gdsc
 

#### ODTMutSqrt

# TREE MUT SQRT

Plot_aux<-Validation_GDSC[Validation_GDSC$Method=="ODTMutSqrt",]

ODTQR_plots<-NULL

for(j in 1:length(unique(Plot_aux$drug))){
  d<-unique(Plot_aux$drug)[j]
  ind<-Plot_aux$drug==d
  plotDrug<-data.frame(Drug=j,
                       drugname=d,
                       Patients=Plot_aux$patients,
                       IC50=NA,
                       Recommended="No")
  if(d %in% colnames(GDSC_AML)){
    plotDrug$Recommended[ind]<-"Yes"
    plotDrug$IC50[ind]<-GDSC_AML[ind,d]
    plotDrug$IC50[!ind]<-GDSC_AML[!ind, d]
    plotDrug$Recommended[!ind]<-"No"
  }

  ODTQR_plots<-rbind(ODTQR_plots,plotDrug)
  
}

gg4_validation_gdsc<-ggplot(ODTQR_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_manual(values=c("#999999", "#E69F00"))+
  stat_compare_means()+ylab("IC50*")+ggtitle("ODTMutSqrt")

gg4_validation_gdsc

#############################################################
# Multiomics in BeatAML: Expression models
#############################################################
"
Some of the methods only accept as input binary variables.
Although, genomic variants can be transformed into binary variables, gene expression, methylation, or openness of the chromatin are intrinsically continuous variables.
We have included a table showing whether the algorithm accepts only binary inputs (genomic variants only) or whether it also accepts continuous data (gene expression, methylation, etc.).
For the methods that accept continuous variables, we assessed the performance of the predictions (5-fold cross-validation) in the BeatAML dataset using both data sources.
We state the statistical significance using a 2-tail Wilcoxon's test comparing the IC50 using as input either genetic variants or gene expression.
"

# Pre-process gene expression information  

# Keep genes with the largest variance
ngenes <- 1000
expression <- as.data.frame(expression)
rownames(expression) <- expression[,2] # set symbols as the name
expression <- expression[,-c(1,2)] # remove ensembl name and extra symbol col
sdevs <- rowSds(as.matrix(expression))
keep <- which(order(sdevs, decreasing = T) <= ngenes)
expression <- expression[keep,]

# Keep patients that have drug response, mutations and expression
commonPatients <- intersect(rownames(mutations), colnames(expression))
expression <- expression[,commonPatients]
expression <- t(expression)
# create new matrix for patients that have both: mutations and expression
drug_response_matrix_both <- drug_response_matrix[commonPatients, ] 
mutations<-mutations[commonPatients, ]

pseudop <- exp(-.5*drug_response_matrix_both) # Convert weights into pseudoprobabilities
pseudop <- pseudop /rowSums(pseudop) # Probabilities sum up one
 
### Train Models with expression Data
library(robustbase)
set.seed(2022)
folds <- 5
groups <- sample(1:folds, nrow(drug_response_matrix_both), replace = T)

treatmentMnLassoExp <- treatmentLassoExp <- treatmentODTExp <- treatmentODTExpSqrt <- rep(NA, length(groups))
treatmentMnLassoMut2 <- treatmentLassoMut2 <- treatmentODTMut2 <- treatmentODTMutSqrt2 <- rep(NA, length(groups))


for (group in 1:folds) {
  cat("Fold: ", group, "\n")
  dejar <- groups != group
  quitar <- !dejar

  # ODT
  cat("ODT\n")
  ODT <- trainTree(expression[dejar,], drug_response_matrix_both[dejar,], minbucket = 20)
  treatmentODTExp[quitar] <- predictTree(ODT,drug_response_matrix_both[quitar,], expression[quitar,])

  # ODT sqrt
  cat("ODT sqrt\n")
  ODTSqrt <- trainTree(expression[dejar,], (drug_response_matrix_both[dejar,])^.5, minbucket = 20)
  treatmentODTExpSqrt[quitar] <- predictTree(ODTSqrt,drug_response_matrix_both[quitar,], expression[quitar,])

  # Using MnLasso
  cat("Multinomial Lasso\n")
  MnLasso <- trainMnLasso(expression[dejar,], pseudop[dejar,])
  treatmentMnLassoExp[quitar] <- predictMnLasso(MnLasso,drug_response_matrix_both[quitar,],expression[quitar,])

  # Using Lasso
  cat("Lasso\n")
  Lasso <- trainLasso(expression[dejar,], drug_response_matrix_both[dejar,])
  treatmentLassoExp[quitar] <- predictMnLasso(Lasso,drug_response_matrix_both[quitar,],expression[quitar,])

  # ODTMut
  cat("ODTMut\n")
  ODTMut <- trainTreeMut(mutations[dejar,], drug_response_matrix_both[dejar,], minbucket = 20)
  treatmentODTMut2[quitar] <- predictTreeMut(ODTMut,drug_response_matrix_both[quitar,], mutations[quitar,])

  # ODTMut sqrt
  cat("ODTMut sqrt\n")
  ODTMutSqrt <- trainTreeMut(mutations[dejar,], (drug_response_matrix_both[dejar,])^.5, minbucket = 20)
  treatmentODTMutSqrt2[quitar] <- predictTreeMut(ODTMutSqrt,drug_response_matrix_both[quitar,], mutations[quitar,])

  # Using MnLasso
  cat("Multinomial Lasso Mut\n")
  MnLassoMut <- trainMnLasso(mutations[dejar,], pseudop[dejar,])
  treatmentMnLassoMut2[quitar] <- predictMnLasso(MnLassoMut,drug_response_matrix_both[quitar,],mutations[quitar,])

  # Using Lasso
  cat("Lasso Mut\n")
  LassoMut <- trainLasso(mutations[dejar,], drug_response_matrix_both[dejar,])
  treatmentLassoMut2[quitar] <- predictLasso(LassoMut,drug_response_matrix_both[quitar,],mutations[quitar,])

}

set.seed(2021)
folds <- 5
groups <- sample(1:5, 257, replace = T) # MAGIC NUMBER!!

# Use same groups as previously
treatmentBOSOExp <- treatmentBOSOMut2 <- rep(NA, length(groups))

for (group in 1:folds) {
  dejar <- groups != group
  quitar <- !dejar

  BOSOMut <- trainBOSO(mutations[dejar,], drug_response_matrix_both[dejar,],maxVarsBlock = 10, standardize=F, IC = "eBIC")
  treatmentBOSOMut2[quitar] <- predictBOSO(BOSOMut, drug_response_matrix_both[quitar,], mutations[quitar,])
  BOSO <- trainBOSO(expression[dejar,], drug_response_matrix_both[dejar,],maxVarsBlock = 10, standardize=F, IC = "eBIC")
  treatmentBOSOExp[quitar] <- predictBOSO(BOSO, drug_response_matrix_both[quitar,], expression[quitar,])
}

save(treatmentODTMut2, treatmentODTMutSqrt2, treatmentMnLassoMut2, treatmentLassoMut2, treatmentBOSOMut2, file = "C:/Users/ksada/OneDrive - Tecnun/Paper XAI Methods/Rdata/mutations2.RData")
save(treatmentODTExp, treatmentODTExpSqrt, treatmentMnLassoExp, treatmentLassoExp, treatmentBOSOExp, file = "C:/Users/ksada/OneDrive - Tecnun/Paper XAI Methods/Rdata/expression.RData")

load("C:/Users/ksada/OneDrive - Tecnun/Paper XAI Methods/Rdata/mutations2.RData")
load("C:/Users/ksada/OneDrive - Tecnun/Paper XAI Methods/Rdata/expression.RData")


treatmentOracle <- apply(drug_response_matrix_both, 1, which.min)

treatments <- cbind(drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentOracle)],
                      drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentODTExp)],
                      drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentODTExpSqrt)],
                      drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentMnLassoExp)],
                      drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentLassoExp)],
                      drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentBOSOExp)],
                      drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentODTMut2)],
                      drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentODTMutSqrt2)],
                      drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentMnLassoMut2)],
                      drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentLassoMut2)],
                      drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentBOSOMut2)]
)
#
colnames(treatments) <- c("Oracle","ODTExp","ODTSqrtExp","MnLassoExp","LassoExp","BOSOExp",
                            "ODTMut","ODTSqrtMut","MnLassoMut","LassoMut","BOSOMut")

# Plotting Results.... 
boxplot(treatments)
plot_Ex_Mut<-data.frame(Method=rep(colnames(treatments), each=nrow(treatments)))
plot_Ex_Mut$IC50<-unlist(as.data.frame(treatments))
plot_Ex_Mut$Type<-"GE"
plot_Ex_Mut$Type[grep("Mut", plot_Ex_Mut$Method)]<-"Mut" # add type "Mut" to mutation values


plot_Ex_Mut<-plot_Ex_Mut[order(plot_Ex_Mut$Method),]
plot_Ex_Mut$General_Method<-"unknown"
plot_Ex_Mut$General_Method[grep("Oracle", plot_Ex_Mut$Method)]<-"Oracle"
plot_Ex_Mut$General_Method[grep("BOSO", plot_Ex_Mut$Method)]<-"BOSO"
plot_Ex_Mut$General_Method[grep("ODT", plot_Ex_Mut$Method)]<-"ODT"
plot_Ex_Mut$General_Method[grep("ODTSqrt", plot_Ex_Mut$Method)]<-"ODTSqrt"
plot_Ex_Mut$General_Method[grep("Lasso", plot_Ex_Mut$Method)]<-"Lasso"
plot_Ex_Mut$General_Method[grep("MnLasso", plot_Ex_Mut$Method)]<-"MnLasso"

# plot_Ex_Mut<-rbind(plot_Ex_Mut, data.frame(Method=NA, IC50=5, Type="Mut", General_Method="MOM"))

plot_Ex_Mut$Method<-factor(plot_Ex_Mut$Method, levels=c("Oracle","BOSOExp","BOSOMut",
                                                        "ODTSqrtExp","ODTSqrtMut",
                                                        "LassoExp","LassoMut",
                                                        "ODTExp","ODTMut",
                                                        "MnLassoExp","MnLassoMut"))
plot_Ex_Mut$General_Method<-factor(plot_Ex_Mut$General_Method, levels=c("Oracle", "BOSO", "ODTSqrt", "Lasso",
                                                                        "ODT","MnLasso"))


mycomparisons<-list(c("ODTExp", "ODTMut"), c("ODTSqrtExp", "ODTSqrtMut"), c("LassoExp", "LassoMut"), c("BOSOExp", "BOSOMut"),
                    c("MnLassoExp", "MnLassoMut"))
mypal = pal_npg("nrc")(9)
mypal<-mypal[c(1,3:7)]
gg_ex_mut<-ggplot(plot_Ex_Mut, aes(x=Method, y=IC50, fill=General_Method, alpha=Type))+
  geom_boxplot(lwd=0.8,position=position_dodge(0.8))+
  scale_fill_manual(values=mypal)+theme_bw()+ggtitle("5-fold CV in BeatAML")+ylab("IC50*")+
  labs(fill="General Method")+
  stat_compare_means(aes(Method),vjust=0, hjust=-1, comparisons=mycomparisons,label.y = 5)+
  scale_alpha_manual(values=c( 0.5, 1))+
  theme(text = element_text(size = 20, family = "Roboto"),
  axis.text.x = element_text(size = 11),
  plot.background = element_rect(fill = "transparent", size = 0, color = "black"),
  panel.background = element_rect(fill = "transparent"))
gg_ex_mut

ggsave("treatment_plot_both.png", gg_ex_mut, width = 14, height = 6, dpi = 1000)

#############################################################
# Multiomics in BeatAML: Expression and Mutation models
#############################################################

# Join expression with mutations 
exp_mut <- merge(expression, mutations, by = "row.names", by.row = TRUE)
rownames(exp_mut) <- exp_mut$Row.names
exp_mut$Row.names <- NULL
exp_mut <- data.matrix(exp_mut, rownames.force = TRUE)

folds <- 5
groups <- sample(1:folds, nrow(drug_response_matrix_both), replace = T)
treatmentMnLassoBoth <- treatmentLassoBoth <- treatmentODTBoth <- treatmentODTSqrtBoth <- rep(NA, length(groups))

# Train Models with EXPRESSION AND MUTATIONS
for (group in 1:folds) {
  cat("Fold: ", group, "\n")
  dejar <- groups != group
  quitar <- !dejar

  # ODT
  cat("ODT\n")
  tree <- trainTree(exp_mut[dejar,], drug_response_matrix_both[dejar,], minbucket = 20)
  treatmentODTBoth[quitar] <- predictTree(tree,drug_response_matrix_both[quitar,], exp_mut[quitar,])

  # ODT sqrt
  cat("ODT sqrt\n")
  tree1 <- trainTree(exp_mut[dejar,], (drug_response_matrix_both[dejar,])^.5, minbucket = 20)
  treatmentODTSqrtBoth[quitar] <- predictTree(tree1,drug_response_matrix_both[quitar,], exp_mut[quitar,])

  # Using MnLasso
  cat("Multinomial Lasso\n")
  MnLasso <- trainMnLasso(exp_mut[dejar,], pseudop[dejar,])
  treatmentMnLassoBoth[quitar] <- predictMnLasso(MnLasso,drug_response_matrix_both[quitar,],exp_mut[quitar,])

  # Using Lasso
  cat("Lasso\n")
  Lasso <- trainLasso(exp_mut[dejar,], drug_response_matrix_both[dejar,])
  treatmentLassoBoth[quitar] <- predictMnLasso(Lasso,drug_response_matrix_both[quitar,],exp_mut[quitar,])
}

set.seed(2021)
folds <- 5
groups <- sample(1:folds, length(groups), replace = T)

# Use same groups as previously
treatmentBOSOBoth <- rep(NA, length(groups))

for (group in 1:folds) {
  dejar <- groups != group
  quitar <- !dejar

  BOSO <- trainBOSO(exp_mut[dejar,], drug_response_matrix_both[dejar,],maxVarsBlock = 10, standardize=F, IC = "eBIC")
  treatmentBOSOBoth[quitar] <- predictBOSO(BOSO,drug_response_matrix_both[quitar,], exp_mut[quitar,])
}

save(treatmentODTBoth, treatmentODTSqrtBoth, treatmentMnLassoBoth, treatmentLassoBoth, treatmentBOSOBoth, file = "C:/Users/ksada/OneDrive - Tecnun/Paper XAI Methods/Rdata/both.RData")

load("C:/Users/ksada/OneDrive - Tecnun/Paper XAI Methods/Rdata/both.RData")

treatmentOracle <- apply(drug_response_matrix_both, 1, which.min)

# Plot all results!!!

treatments2 <- cbind(
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentOracle)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentODTExp)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentODTExpSqrt)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentMnLassoExp)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentLassoExp)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentBOSOExp)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentODTMut2)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentODTMutSqrt2)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentMnLassoMut2)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentLassoMut2)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentBOSOMut2)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentODTBoth)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentODTSqrtBoth)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentMnLassoBoth)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentLassoBoth)],
  drug_response_matrix_both[cbind(1:nrow(drug_response_matrix_both), treatmentBOSOBoth)]
)

colnames(treatments2) <- c(
  "Oracle", "ODTExp", "ODTSqrtExp", "MnLassoExp", "LassoExp", "BOSOExp",
  "ODTMut", "ODTSqrtMut", "MnLassoMut", "LassoMut", "BOSOMut",
  "ODTBoth", "ODTSqrtBoth", "MnLassoBoth", "LassoBoth", "BOSOBoth"
)
boxplot(treatments2)

mycomparisons2 <- list(
  c("ODTExp", "ODTMut"), c("ODTMut", "ODTBoth"),
  c("ODTSqrtExp", "ODTSqrtMut"),c("ODTSqrtMut", "ODTSqrtBoth"),
  c("LassoExp", "LassoMut"),c("LassoMut", "LassoBoth"),
  c("BOSOExp", "BOSOMut"),c("BOSOMut", "BOSOBoth"),
  c("MnLassoExp", "MnLassoMut"),c("MnLassoMut", "MnLassoBoth")
)

mycomparisons3 <- list(
  c("ODTExp", "ODTBoth"),
  c("ODTSqrtExp", "ODTSqrtBoth"),
  c("LassoExp", "LassoBoth"),
  c("BOSOExp", "BOSOBoth"),
  c("MnLassoExp", "MnLassoBoth")
)


# Plotting Results.... 
plot_all<-data.frame(Method=rep(colnames(treatments2), each=nrow(treatments2)))
plot_all$IC50<-unlist(as.data.frame(treatments2))
plot_all$Type<-"GE"
plot_all$Type[grep("Mut", plot_all$Method)]<-"Mut" # add type "Mut" to mutation values
plot_all$Type[grep("Both", plot_all$Method)]<-"Both"


plot_all<-plot_all[order(plot_all$Method),]
plot_all$General_Method<-"unknown"
plot_all$General_Method[grep("Oracle", plot_all$Method)]<-"Oracle"
plot_all$General_Method[grep("BOSO", plot_all$Method)]<-"BOSO"
plot_all$General_Method[grep("ODT", plot_all$Method)]<-"ODT"
plot_all$General_Method[grep("ODTSqrt", plot_all$Method)]<-"ODTSqrt"
plot_all$General_Method[grep("Lasso", plot_all$Method)]<-"Lasso"
plot_all$General_Method[grep("MnLasso", plot_all$Method)]<-"MnLasso"

# plot_all<-rbind(plot_all, data.frame(Method=NA, IC50=5, Type="Mut", General_Method="MOM"))

plot_all$Method<-factor(plot_all$Method, levels=c("Oracle","BOSOExp","BOSOMut","BOSOBoth",
                                                        "ODTSqrtExp","ODTSqrtMut", "ODTSqrtBoth",
                                                        "LassoExp","LassoMut","LassoBoth",
                                                        "ODTExp","ODTMut","ODTBoth",
                                                        "MnLassoExp","MnLassoMut","MnLassoBoth"))
plot_all$General_Method<-factor(plot_all$General_Method, levels=c("Oracle", "BOSO", "ODTSqrt", "Lasso",
                                                                        "ODT","MnLasso"))
plot_all$Type<-factor(plot_all$Type, levels=c("GE", "Mut", "Both"))


mypal = pal_npg("nrc")(9)
mypal<-mypal[c(1,3:7)]
gg_ex_both<-ggplot(plot_all, aes(x=Method, y=IC50, fill=General_Method))+    # Different pattern for each group
  geom_boxplot_pattern(
                   pattern_color = "black",
                   pattern_fill = "black",
                   pattern_spacing = 0.02,
                   pattern_density      = 0.01,
                   aes(pattern = Type))+

  scale_fill_manual(values=mypal)+theme_bw()+ggtitle("5-fold CV in BeatAML")+ylab("IC50*")+
  labs(fill="General Method")+
  stat_compare_means(aes(Method),vjust=0, hjust=-1, comparisons=mycomparisons2,label.y = 5)+
  stat_compare_means(aes(Method),vjust=0, hjust=-1, comparisons=mycomparisons3,label.y = 5.4)+
  #scale_alpha_manual(values=c( 0.3,0.6,1))+
  theme(text = element_text(size = 20, family = "Roboto"),
  axis.text.x = element_text(size = 11),
  plot.background = element_rect(fill = "transparent", size = 0, color = "black"),
  panel.background = element_rect(fill = "transparent"))
gg_ex_both

ggsave("treatment_plot_all.png", gg_ex_both, width = 18, height = 6, dpi = 1000)

gg_ex_both<-ggplot(plot_all, aes(x=Method, y=IC50, fill=General_Method))+    # Different pattern for each group
  geom_boxplot_pattern(
                   pattern_color = "white",
                   pattern_fill = "black",
                   aes(pattern = Type))
gg_ex_both
ggsave("treatment_plot_all2.png", gg_ex_both, width = 18, height = 6, dpi = 1000)
