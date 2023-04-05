################################################################################
# MOM K-fold Cross Validation
# Author: Marian Gimeno email:mgimenoc@unav.es
# Supervisor: Angel Rubio
################################################################################



library(reticulate)
# use_python("C:/Users/mgimenoc/Anaconda3/envs/Tesis/python.exe")
use_python("C:/Users/ksada/Anaconda3/envs/AML_Patients_PHD/python.exe")

library(impute)
library(knitr)
library(reticulate)
library(readxl)
library(readr)
library(limma)
library(ggplot2)
library(ggpubr)
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

# Folder to save results
folder_dir <- "C:/Users/ksada/OneDrive - Tecnun/Paper XAI Methods/Rresults/"

# Loading Data -------------------------------------------- 
# Load clinical data 
clinical <- read_excel("./data/input/BeatAML2/beataml_wv1to4_clinical.xlsx",
                       sheet = "summary"
)
clinical$dbgap_subject_id <- as.character(clinical$dbgap_subject_id)
# keep only first specimen collected of each subject
clinical <- clinical[!clinical$timeOfSampleCollectionRelativeToInclusion>0,]
# keep subjects that have dna and rna
clinical <- clinical[complete.cases(clinical[, c("dbgap_dnaseq_sample", "dbgap_rnaseq_sample")]),]
# make sure no samples are duplicated: keep only one rna and dna sample per subject (some have more)
sum(duplicated(clinical$dbgap_subject_id)) # no samples are duplicated
# clinical <- clinical[!duplicated(clinical$dbgap_subject_id),]
clinical <- as.data.frame(clinical)
rownames(clinical) <- clinical$dbgap_subject_id

#Then, we will load all the data for the analysis.
#Data selected was the BeatAML2 cohort for which we had mutational data, drug sensitivity data and gene expression data.

# Load gene expression data
expression <- read_excel("./data/input/BeatAML2/beataml_waves1to4_norm_exp_dbgap.xlsx")
expression <- as.data.frame(expression)
# set symbols as row names 
rownames(expression) <- expression$display_label 
# set dbgap_subject_id as colnames
rna_ids <- clinical$dbgap_rnaseq_sample[clinical$dbgap_rnaseq_sample  %in% colnames(expression)]
expression <- expression[,rna_ids]
subject_ids <-na.omit(clinical$dbgap_subject_id[match(colnames(expression),clinical$dbgap_rnaseq_sample)])
colnames(expression) <- subject_ids

# Load mutations data
gene_variants <- read_excel("./data/input/BeatAML2/beataml_wes_wv1to4_mutations_dbgap.xlsx")
colnames(gene_variants)[1]="dbgap_dnaseq_sample"
# add dbgap_subject_id and keep only samples in clinical (have also dna)
gene_variants <-merge(gene_variants,clinical[,c("dbgap_subject_id","dbgap_dnaseq_sample")], by="dbgap_dnaseq_sample")

# Load drug sensitivity data
drug_response <- read_excel("./data/input/BeatAML2/beataml_probit_curve_fits_v4_dbgap.xlsx")
drug_response$dbgap_subject_id <- as.character(drug_response$dbgap_subject_id)

# Pre-processing drug and mutational information  --------------------------------------------

# Now it is time to pre-process the data in order to obtain two different matrices, a matrix containing all patients and their mutations available and a second matrix containing all patients and their sensitivity to the different drugs.
# For doing so, we will impute the missing values in the drug matrix.

# DRUG RESPONSE MATRIX 
drug_names <- unique(drug_response$inhibitor)
subjects_ids <- clinical$dbgap_subject_id

# Build Drug Matrix 
drug_matrix <- matrix( # empty matrix
  data = "NA", nrow = length(drug_names),
  ncol = length(subjects_ids)
) 

colnames(drug_matrix) <- subjects_ids # columns are subjects
rownames(drug_matrix) <- drug_names # rows are drugs

for (subject in subjects_ids) { # add the ic50
  inds_subject <- which(drug_response$dbgap_subject_id == subject) # get indices of subject
  drugs <- as.character(drug_response$inhibitor[inds_subject]) # get drug names of subject
  drug_matrix[drugs, subject] <- drug_response$ic50[inds_subject] # add ic50 of those drugs
}

# Identifying missing values
is.na(drug_matrix) <- drug_matrix == "NA"

# Change IC50 format to numeric
for (i in 2:(dim(drug_matrix)[2] - 1)) {
  drug_matrix[, i] <- as.numeric(drug_matrix[, i])
}

# Delete patients with few drug response and input missing drug response 
drug_matrix_2<-as.matrix(t(drug_matrix)) # transpose drug 
if (sum(is.na(drug_matrix_2)) > 0) { 
  drug_matrix_3 <- drug_matrix_2[
    # NAs percentage has to be less than 20% in rows and less than 30% in columns
    -which((rowSums(is.na(drug_matrix_2)) / dim(drug_matrix_2)[2]) > 0.8), # delete those tested in a few drugs, 
    -which((colSums(is.na(drug_matrix_2)) / dim(drug_matrix_2)[1]) > 0.7)] 
  out <- impute.knn(as.matrix(drug_matrix_3), # impute missing values
                    k = 10, rowmax = 0.8,
                    colmax = 0.75, maxp = 1500, rng.seed = 362436069
  )
  drug_matrix_4 <- out$data[, -which(colnames(out$data) %in% c("Elesclomol", "JNJ-7706621"))] # Los quitamos por ser toxicos, no los queremos, cuando hicimos la normalizacion se deberian de haber ido solitos
} else {
  drug_matrix_4 <- drug_matrix_2[, -which(colnames(drug_matrix_2) %in% c("Elesclomol", "JNJ-7706621"))]
} 

# MUTATIONS MATRIX 
# Build Mutations matrix
mutations <- matrix(0, # empty matrix
                    nrow = length(subjects_ids),
                    ncol = length(unique(gene_variants$symbol))
)

rownames(mutations) <- subjects_ids
colnames(mutations) <-unique(gene_variants$symbol)

for (i in as.character(unique(gene_variants$dbgap_subject_id))) {
  jj <- unique(gene_variants$symbol[which(gene_variants$dbgap_subject_id == i)])
  mutations[i, jj] <- 1 # add 1 if there is a mutation
}

# Add Translocations 
# look for translocations:  filtering the rows to only those that contain "AML with inv" or "AML with t" 
translocations_clinical<-clinical[,c("dbgap_subject_id","specificDxAtAcquisition")][grep("AML with inv", clinical$specificDxAtAcquisition),]
translocations_clinical<-rbind(translocations_clinical,clinical[,c("dbgap_subject_id","specificDxAtAcquisition")][grep("AML with t", clinical$specificDxAtAcquisition),])

trans<-as.character(unique(translocations_clinical$specificDxAtAcquisition))
trans_matrix<-matrix(0, nrow=nrow(mutations), ncol=length(trans)) # add columns for translocations
colnames(trans_matrix)<-trans
rownames(trans_matrix)<-rownames(mutations)

for(j in colnames(trans_matrix)){
  p<-as.character(translocations_clinical$dbgap_subject_id[which(translocations_clinical$specificDxAtAcquisition ==j)])
  p<-p[which(p %in% rownames(trans_matrix))]
  trans_matrix[p,j]<-1 # add 1 if there is a translocation
}

# Change column name and also add t: extract the second element of each translocation name (i.e., the specific type of translocation) and use it as the new column name in add. 
trans<-sapply(trans, function(X){unlist(strsplit(X, split = "; "))[2]})
colnames(trans_matrix)<-trans
mutations<-cbind(mutations,trans_matrix)
mutations<-as.matrix(mutations)

# Keep only samples in the 4 matrices (gene_variants,expression,drug_response,clinical)
drug_subjects <- rownames(drug_matrix_4)
mutation_subjects <- rownames(mutations)
expression_subjects <- colnames(expression)
clinical_subjects <- clinical$dbgap_subject_id

keep_subjects <- Reduce(intersect, list(drug_subjects,mutation_subjects,expression_subjects,clinical_subjects))
drug_matrix_4 <- drug_matrix_4[keep_subjects,]
mutations <- mutations[keep_subjects,]
expression <- expression[,keep_subjects]
clinical <- clinical[keep_subjects,]

table(clinical$cohort)

# IC50* Calculation 
drug<-log10(drug_matrix_4)-1
drug <- t(t(drug) - colMeans(drug, na.rm=T))  # unificar los farmacos por columnas, compensamos con la dosis

drug_response_matrix_MOM<-drug-max(drug,na.rm = T)  # MOM minimiza en vez de maximizar 

# Separate variables into train and test (Waves 1+2 or Both and Waves 3+4) REVISAR
# get ids of subjects
subject_id_w12 <- as.character(clinical$dbgap_subject_id[which(clinical$cohort=="Waves1+2" | clinical$cohort=="Both")])
subject_id_w34 <- as.character(clinical$dbgap_subject_id[which(clinical$cohort=="Waves3+4")])

drug_response_w12 <- drug_response_matrix_MOM[subject_id_w12,]
drug_response_w34 <- drug_response_matrix_MOM[subject_id_w34,]

# Filter: keep only mutations present in at least 1% of patients

keep_mut <- (colSums(mutations[subject_id_w34,]) > 0) & (colSums(mutations[subject_id_w12,]) > 0) & colSums(mutations)>0.01*nrow(mutations)

mutations<-mutations[, keep_mut] 

mutations_w12 <- mutations[subject_id_w12,]
mutations_w34 <- mutations[subject_id_w34,]

# Tidying up the variables (delete extra variables)
rm(list = c("drug","drug_matrix_2", "drug_matrix_3", "drug_matrix_4", "drug_response", "gene_variants"))
rm(list = c("p","j","out", "translocations_clinical", "trans_matrix", "drug_matrix","i","jj")) 

patient_names<-rownames(drug_response_matrix_MOM)

# drug_ori<-drug
# mut_ori<-mut
# 
# identical(colnames(drug), colnames(drug_ori))
# set.seed()

patient_names<-cbind(patient_names, kfold=sample(1:5,length(patient_names), replace=TRUE))
patient_names<-as.data.frame(patient_names)

# pat_kfold<-as.character(patient_names$patient_names[which(patient_names$kfold==11)])
# patient_final_Validation_drug<-drug[pat_kfold,]
# patient_final_Validation_mut<-mut[pat_kfold,]

drug_ori<-drug_response_matrix_MOM
mut_ori<-mutations

sensi<-vector(mode="list", length=10)
Results_list<-vector(mode="list", length=10)

Results_CV_pat<-NULL
kfold<-1
# CV Code-----------------------------------------------------------------------
for(kfold in 1:5){
  tic()
  
  cat(sprintf("##### \n \n FOLD %s/5 \n\n######\n\n", kfold))
  patient_names_kfold<-as.character(patient_names$patient_names[which(as.numeric(patient_names$kfold)==kfold)])
  
  kfold_validation_drug<-drug_ori[c(patient_names_kfold),]
  kfold_validation_mut<-mut_ori[patient_names_kfold,]
  
  drug<-drug_ori[-c(which(rownames(drug_ori) %in% patient_names_kfold)),]
  mut<-mut_ori[-c(which(rownames(mut_ori) %in% patient_names_kfold)),]
  mut<-mut[,colSums(mut)>0.01*nrow(mut)]
  
  Result<-getBMeffect_F(Mutations = mut, Drug=drug)
  output<-Get_Treatment_OptimizationMatrix(BM_Result = Result, Drug = drug, Mutations = mut)
  
  IC50<-output$`IC50*`
  PT<-output$PT
  P<-nrow(PT)
  Treat<-ncol(PT)
  S<-4
  sol_pat<-MOM_MILP(IC50, PT, S)
  
  MILP_classification<-Extract_results(sol_pat, PT)
  
  Prediction<-CV_Prediction2(MILP_classification = MILP_classification, 
                             kfold_validation_drug = kfold_validation_drug,
                             kfold_validation_mut = kfold_validation_mut )
  
  Prediction$Fold<-kfold
  Results_CV_pat<-rbind(Results_CV_pat, Prediction)
  toc()
}
treatmentMOM<-Results_CV_pat
treatmentMOM$IC50<-treatmentMOM$IC50-min(drug_ori)
MOM<-MILP_classification
save(MOM, file = paste(folder_dir,"Rdata/MOM_model.RData",sep=""))
save(treatmentMOM, file = paste(folder_dir,"Rdata/MOM_results.RData",sep=""))
# load(paste(folder_dir,"Rdata/MOM_model.RData",sep=""))
# load(paste(folder_dir,"Rdata/MOM_results.RData",sep=""))

### COMPARISON AGAINST ORACLE---------------------------------------------------
drug<-drug_ori-min(drug_ori)

ORACLE<- apply(drug, 1, which.min)
identical(names(ORACLE), rownames(drug))
ORACLE<-as.data.frame(ORACLE)
ORACLE$Drug<-colnames(drug)[ORACLE[,1]]
ORACLE$Patient<-rownames(ORACLE)

aux<-apply(ORACLE, 1, FUN=function(X){return(IC50<-drug[X[3], X[2]])})

ORACLE$IC50<-aux[rownames(ORACLE)]
MOM<-Results_CV_pat
MOM$IC50<-MOM$IC50-min(drug_ori)

boxplot(ORACLE$IC50, MOM$IC50, names = c("ORACLE", "MOM"), ylab="IC50*", ylim=c(0,5))

ORACLE_plot<-data.frame(Treatment=ORACLE$Drug, 
                        IC50=ORACLE$IC50, 
                        #Pat_name=ORACLE$Patient,
                        Method="ORACLE",
                        Fold=NA)

ORACLE_plot<-rbind(ORACLE_plot,cbind(MOM, Method="MOM"))


ggplot(ORACLE_plot, aes(x=Method, y=IC50, fill=Method))+geom_boxplot()+theme_bw()+ylim(0,5)+ylab("IC50*")



MOM$DeltaIC50<-abs(MOM$IC50-ORACLE$IC50[na.omit(match(ORACLE$Patient,MOM$Pat_name))])
MOM$Score<-median(MOM$DeltaIC50)




# save(MOM, ORACLE, file="./data/output/mom_COMPARISON.RData")


# Train the model with all the data---------------------------------------------
tic()
Result<-getBMeffect_F(Mutations = mutations_w12, Drug=drug_response_w12)
output<-Get_Treatment_OptimizationMatrix(BM_Result = Result, Drug = drug_response_w12, Mutations = mutations_w12)

IC50<-output$`IC50*`
PT<-output$PT
P<-nrow(PT)
Treat<-ncol(PT)
S<-4 # AÑADIIIIIIIIR
tic()
sol_pat<-MOM_MILP(IC50, PT)
toc()

MILP_classification<-Extract_results(sol_pat, PT)

Prediction34<-CV_Prediction2(MILP_classification = MILP_classification, 
                             kfold_validation_drug = drug_response_w34,
                             kfold_validation_mut = mutations_w34 )
toc()
treatmentMOM_w34<-Prediction34
treatmentMOM_w34$IC50<-treatmentMOM_w34$IC50-min(drug_response_w34)
MOM_All <- MILP_classification
save(MOM_All, file = paste(folder_dir,"Rdata/MOM_model_all.RData",sep=""))
save(treatmentMOM_w34, file = paste(folder_dir,"Rdata/MOM_results_test.RData",sep=""))
load(paste(folder_dir,"Rdata/MOM_model_all.RData",sep=""))
load(paste(folder_dir,"Rdata/MOM_results_test.RData",sep=""))


drug<-drug_ori-min(drug_response_w34)

ORACLE<- apply(drug, 1, which.min)
identical(names(ORACLE), rownames(drug))
ORACLE<-as.data.frame(ORACLE)
ORACLE$Drug<-colnames(drug)[ORACLE[,1]]
ORACLE$Patient<-rownames(ORACLE)

aux<-apply(ORACLE, 1, FUN=function(X){return(IC50<-drug[X[3], X[2]])})

ORACLE$IC50<-aux[rownames(ORACLE)]

MOM<-MOM_w34
MOM$IC50<-MOM$IC50-min(drug_ori)

boxplot(ORACLE$IC50, MOM$IC50, names = c("ORACLE", "MOM"), ylab="IC50*", ylim=c(0,5))

ORACLE_plot<-data.frame(Treatment=ORACLE$Drug, 
                        IC50=ORACLE$IC50, 
                        #Pat_name=ORACLE$Patient,
                        Method="ORACLE")

ORACLE_plot<-rbind(ORACLE_plot,cbind(MOM, Method="MOM"))


ggplot(ORACLE_plot, aes(x=Method, y=IC50, fill=Method))+geom_boxplot()+theme_bw()+ylim(0,5)+ylab("IC50*")



MOM$DeltaIC50<-abs(MOM$IC50-ORACLE$IC50[na.omit(match(ORACLE$Patient,MOM$Pat_name))])
MOM$Score<-median(MOM$DeltaIC50)