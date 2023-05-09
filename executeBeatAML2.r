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

windowsFonts("Roboto" = windowsFont("Roboto"))

# Folder to save results
folder_dir <- "./data/output"

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
drug_matrix_KRL <- drug_matrix_2 # KRL needs the NAs, we dont input data
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
drug_matrix_2 <- drug_matrix_2[keep_subjects,colnames(drug_matrix_4)] # same drugs and patients
mutations <- mutations[keep_subjects,]
expression <- expression[,keep_subjects]
clinical <- clinical[keep_subjects,]

table(clinical$cohort)

# Filter: keep only mutations present in at least 1% of patients 
mutations<-mutations[, colSums(mutations)>0.01*nrow(mutations)] 

# IC50* Calculation 
drug<-log10(drug_matrix_4)-1
drug <- t(t(drug) - colMeans(drug, na.rm=T))  # unify the drugs by columns, we compensate with the dosage

drug_response_matrix_MOM<-drug-max(drug,na.rm = T)  # MOM minimizes, the values must all be positive.
drug_response_matrix <- drug - min(drug) # 

drug_matrix_KRL <- apply(drug_matrix_2,2,as.numeric)
rownames(drug_matrix_KRL)<-rownames(drug_matrix_2[keep_subjects,])
drug_matrix_KRL<-log10(drug_matrix_KRL)-1
drug_matrix_KRL <- t(t(drug_matrix_KRL) - colMeans(drug_matrix_KRL, na.rm=T)) 
drug_matrix_KRL<-(drug_matrix_KRL-max(drug_matrix_KRL,na.rm = T))*(-1)

pseudop <- exp(-.5*drug_response_matrix) # Convert weights into pseudoprobabilities for Multinomial
pseudop <- pseudop /rowSums(pseudop) # Probabilities sum up one

# Separate variables into train and test (Waves 1+2 or Both and Waves 3+4) REVISAR
# get ids of subjects
subject_id_w12 <- as.character(clinical$dbgap_subject_id[which(clinical$cohort=="Waves1+2" | clinical$cohort=="Both")])
subject_id_w34 <- as.character(clinical$dbgap_subject_id[which(clinical$cohort=="Waves3+4")])

drug_response_w12 <- drug_response_matrix[subject_id_w12,]
drug_response_w34 <- drug_response_matrix[subject_id_w34,]

# drug_response_MOM_w12 <- drug_response_matrix_MOM[subject_id_w12,]
# drug_response_MOM_w34 <- drug_response_matrix_MOM[subject_id_w34,]

drug_matrix_KRL_w12 <- drug_matrix_KRL[subject_id_w12,]
drug_matrix_KRL_w34 <- drug_matrix_KRL[subject_id_w34,]

pseudop_w12 <- pseudop[subject_id_w12,]
pseudop_w34 <- pseudop[subject_id_w34,]

mutations_w12 <- mutations[subject_id_w12,]
mutations_w34 <- mutations[subject_id_w34,]

expression_w12 <- expression[,subject_id_w12]
expression_w34 <- expression[,subject_id_w34]

# FOR KRL
# write.csv(drug_matrix_KRL_w12, file="data/input/drug_matrix_KRL_w12_NAs.csv")
# write.csv(mutations_w12, file="data/input/mutations_w12L_NAs.csv")
# write.csv(drug_matrix_KRL_w34, file="data/input/drug_matrix_KRL_w34_NAs.csv")
# write.csv(mutations_w34, file="data/input/mutations_w34_NAs.csv")

# Tidying up the variables (delete extra variables)
rm(list = c("drug","drug_matrix_2", "drug_matrix_3", "drug_matrix_4", "drug_response", "gene_variants"))
rm(list = c("p","j","out", "translocations_clinical", "trans_matrix", "drug_matrix","i","jj")) 

# Save to use in GDSC testing
# save(drug_response_w12,pseudop_w12,expression_w12,mutations_w12, file = paste(folder_dir,"Rdata/training_data.RData",sep=""))
# load(file = paste(folder_dir,"Rdata/training_data.RData",sep=""))

#############################################################
# Mutation models
#############################################################

# Comparison of methods with Cross-Validation --------------------------------------------
"
We performed a 5-fold cross-validation using the BeatAML dataset.
We trained all models with genetic variants data from 319 patients, dividing the cohort between the training samples 4-folds and testing samples the selected 1-fold.
Each of the folds were tested, and the predicted IC50* for the 5-fold testing was compared for all the methods and compared against the Oracle -the drug with the optimum IC50.
We calculated the Oracle as the minimum IC50* value for each patient.
"

#Results from MOM and KRL were directly loaded into the workspace due to the computing limitations in R.
set.seed(2022)
folds <- 5
groups <- sample(1:folds, nrow(drug_response_w12), replace = T)
treatmentMnLassoMut <- treatmentLassoMut <- treatmentODTMut <- treatmentODTSqrtMut <- rep(NA, length(groups))

# Predict functions require the drug_response_matrix just to extract the names of the drugs
# Train and test
for (group in 1:folds) {
  cat("Fold: ", group, "\n")
  dejar <- groups != group
  quitar <- !dejar
  
  # ODTMut
  cat("ODTMut\n")
  ODTMut <- trainTreeMut(mutations_w12[dejar,], drug_response_w12[dejar,], minbucket = 10)
  treatmentODTMut[quitar] <- predictTreeMut(ODTMut,drug_response_w12[quitar,], mutations_w12[quitar,])
  
  # ODTMut sqrt  
  cat("ODTMut sqrt\n")
  ODTSqrtMut <- trainTreeMut(mutations_w12[dejar,], (drug_response_w12[dejar,])^.5, minbucket = 10)
  treatmentODTSqrtMut[quitar] <- predictTreeMut(ODTSqrtMut,drug_response_w12[quitar,], mutations_w12[quitar,])
  
  # Using MnLasso
  cat("Multinomial Lasso Mut\n")
  MnLassoMut <- trainMnLasso(mutations_w12[dejar,], pseudop_w12[dejar,])
  treatmentMnLassoMut[quitar] <- predictMnLasso(MnLassoMut,drug_response_w12[quitar,],mutations_w12[quitar,])
  
  # Using Lasso
  cat("Lasso Mut\n")
  LassoMut <- trainLasso(mutations_w12[dejar,], drug_response_w12[dejar,])
  treatmentLassoMut[quitar] <- predictLasso(LassoMut,drug_response_w12[quitar,],mutations_w12[quitar,])
  
}

# Use same groups as previously voy aqui
treatmentBOSOMut <- rep(NA, length(groups))
for (group in 1:folds) {
  cat("Fold: ", group, "\n")
  dejar <- groups != group
  quitar <- !dejar
  BOSOMut <- trainBOSO(mutations_w12[dejar,], drug_response_w12[dejar,],maxVarsBlock = 10, standardize=F, IC = "eBIC")
  treatmentBOSOMut[quitar] <- predictBOSO(BOSOMut, drug_response_w12[quitar,], mutations_w12[quitar,])
}

# save(ODTMut, ODTSqrtMut, MnLassoMut,LassoMut, BOSOMut, file = paste(folder_dir,"Rdata/mutations_models.RData",sep=""))
# save(treatmentODTMut, treatmentODTSqrtMut, treatmentMnLassoMut, treatmentLassoMut, treatmentBOSOMut, file = paste(folder_dir,"Rdata/mutations_results.RData",sep=""))
load(paste(folder_dir,"/mutations_models.RData",sep=""))
load(paste(folder_dir,"/mutations_results.RData",sep=""))

# load(paste(folder_dir,"/MOM_model.RData",sep="")) # model
load(paste(folder_dir,"/MOM_results.RData",sep="")) # predictions

# KRL
KRL_results_cv <- read_csv("data/output/KRL_results_5cv_BeatAML2.csv")
KRL_results_cv$Drug<-KRL_results_cv$Drug+1 
KRL_results_cv$IC50<-drug_response_w12[cbind(KRL_results_cv$Patients, KRL_results_cv$DrugName)]

treatmentOracle <- apply(drug_response_w12, 1, which.min)
### Plotting BeatAML Cross-validation Results and comparing means difference
treatment_plot<-data.frame(Method="ORACLE", IC50=drug_response_w12[cbind(1:nrow(drug_response_w12),treatmentOracle)])
treatment_plot<-rbind(treatment_plot,data.frame(Method="ODTMut", IC50=drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentODTMut)]))
treatment_plot<-rbind(treatment_plot,data.frame(Method="ODTSqrtMut", IC50=drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentODTSqrtMut)]))
treatment_plot<-rbind(treatment_plot,data.frame(Method="MnLassoMut", IC50=drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentMnLassoMut)]))
treatment_plot<-rbind(treatment_plot,data.frame(Method="LassoMut", IC50=drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentLassoMut)]))
treatment_plot<-rbind(treatment_plot,data.frame(Method="BOSOMut", IC50=drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentBOSOMut)]))
treatment_plot<-rbind(treatment_plot,data.frame(Method="MOM", IC50=na.omit(treatmentMOM$IC50)))
treatment_plot<-rbind(treatment_plot,data.frame(Method="KRL", IC50=KRL_results_cv$IC50)) # CAMBIAAAR

treatment_plot$Method<-factor(treatment_plot$Method, levels=c("ORACLE", "MOM","KRL", "BOSOMut", "ODTMut","ODTSqrtMut", "LassoMut", "MnLassoMut"))

colors <- c("#A83333","#367592", "#39A7AE", "#96D6B6", "#FDE5B0", "#F3908B", "#E36192", "#8E4884")
treatment_plot_gg <- ggplot(treatment_plot, aes(x=Method, y=IC50, fill=Method)) +
  geom_boxplot(lwd=0.8) +
  theme_bw()+
  scale_fill_manual(values = colors)+
  #scale_fill_npg()
  ylab("IC50*")+
  ggtitle("5-fold CV in BeatAML (Waves 1+2)")+
  theme(text = element_text(family = "Roboto"),
        axis.text.x = element_text(size = 24, angle=60,vjust=0.65),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 27),
        plot.title = element_text(size = 32, face = "bold"),
        plot.background = element_rect(fill = "transparent", size = 0, color = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        legend.position = "NULL")+
  ylim(0, 5)
treatment_plot_gg
ggsave(paste(folder_dir,"/images/boxplot_mutations.png",sep=""), treatment_plot_gg, width = 8, height = 7, dpi = 1000)

treatment_plot_violin <- ggplot(treatment_plot, aes(x=Method, y=IC50, fill=Method)) +
  geom_violin(lwd=0.8,position=position_dodge(0.8)) +
  theme_bw()+
  scale_fill_manual(values = colors)+
  #scale_fill_npg()
  ylab("IC50*")+
  ggtitle("5-fold CV in BeatAML (Waves 1+2)")+
  theme(text = element_text(family = "Roboto"),
        axis.text.x = element_text(size = 18, angle=60,vjust=0.65),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 21),
        plot.title = element_text(size = 26, face = "bold"),
        plot.background = element_rect(fill = "transparent", size = 0, color = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        )+
  ylim(0, 5)
treatment_plot_violin
ggsave(paste(folder_dir,"/images/violin_mutations.png",sep=""), treatment_plot_violin, width = 8, height = 7, dpi = 1000)

# Compute wilcoxon tests 
with(treatment_plot,                       # Order boxes by median
     reorder(Method,
             IC50,
             median))

with(treatment_plot,                       # Order boxes by median
     reorder(Method,
             IC50,
             var))

# library(Metrics)
# lapply(unique(treatment_plot$Method), function(method){
#   return(c(as.character(method),rmse(treatment_plot$IC50[treatment_plot$Method=="ORACLE"],treatment_plot$IC50[treatment_plot$Method==method])))
# })
# rmse(treatment_plot$IC50[treatment_plot$Method=="ORACLE"][-which(is.na(treatmentMOM$IC50))],na.omit(treatmentMOM$IC50))
# 

wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ODTMut"], 
            treatment_plot$IC50[treatment_plot$Method=="MnLassoMut"], alternative="less")

wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="MOM"], alternative="less")
wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="BOSOMut"], alternative="less")
wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="ODTSqrtMut"], alternative="less")
wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="LassoMut"], alternative="less")
wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="ODTMut"], alternative="less")
wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="MnLassoMut"], alternative="less")
wilcox.test(treatment_plot$IC50[treatment_plot$Method=="ORACLE"], 
            treatment_plot$IC50[treatment_plot$Method=="KRL"], alternative="less")

# Training models in complete BeatAML data and timing
"
We trained the models in BeatAML cohort and measured the training time for each model except KRL and MOM, that were meassured independently.
Technical implementation refers to the computational burden and software that the method requires.Despite it could be considered less important, some of the algorithms require hours of computing time for the BeatAML of the subset of AML samples in GDSC -that be considered to be small/medium size.
" 
# Train with ALL the data
# ODTMut
cat("ODTMut\n")
tic()
ODTMut_All <- trainTreeMut(mutations_w12, drug_response_w12, minbucket = 10)
treatmentODTMut_w34 <- predictTreeMut(ODTMut_All,drug_response_w34, mutations_w34)
toc()

# ODTMut sqrt
cat("ODTMut sqrt\n")
tic()
ODTSqrtMut_All <- trainTreeMut(mutations_w12, (drug_response_w12)^.5, minbucket = 10)
treatmentODTSqrtMut_w34 <- predictTreeMut(ODTSqrtMut_All,drug_response_w34, mutations_w34)
toc()

# Using MnLasso
cat("Multinomial Lasso Mut\n")
tic()
MnLassoMut_All <- trainMnLasso(mutations_w12, pseudop_w12)
treatmentMnLassoMut_w34 <- predictMnLasso(MnLassoMut_All,pseudop_w34,mutations_w34)
toc()

# Using Lasso
cat("Lasso Mut\n")
tic()
LassoMut_All <- trainLasso(mutations_w12, drug_response_w12)
treatmentLassoMut_w34 <- predictLasso(LassoMut_All,drug_response_w34,mutations_w34)
toc()

# Using Bosso
cat("Bosso Mut\n")
tic()
BOSOMut_All <- trainBOSO(mutations_w12, drug_response_w12,maxVarsBlock = 10, standardize=F, IC = "eBIC")
treatmentBOSOMut_w34 <- predictBOSO(BOSOMut_All,drug_response_w34,mutations_w34)
toc()

# save(treatmentODTMut_w34, treatmentODTSqrtMut_w34, treatmentMnLassoMut_w34, treatmentLassoMut_w34, treatmentBOSOMut_w34, file = paste(folder_dir,"Rdata/mutations_results_test.RData",sep=""))
# save(ODTMut_All, ODTSqrtMut_All, MnLassoMut_All, LassoMut_All,BOSOMut_All, file = paste(folder_dir,"Rdata/mutations_models_all.RData",sep=""))
load(paste(folder_dir,"/mutations_results_test.RData",sep=""))
load(paste(folder_dir,"/mutations_models_all.RData",sep=""))

load(paste(folder_dir,"/MOM_model_all.RData",sep="")) # model
load(paste(folder_dir,"/MOM_results_test.RData",sep="")) # predictions

# KRL
treatmentKRL_w34 <- read_csv(paste(folder_dir, "/KRL_results_w34_BeatAML2.csv",sep=""))
treatmentKRL_w34$`Best Drug`<-treatmentKRL_w34$`Best Drug`+1 
treatmentKRL_w34$IC50<-drug_response_w34[cbind(treatmentKRL_w34$Patients, treatmentKRL_w34$DrugName)]

treatmentOracle_w34 <- apply(drug_response_w34, 1, which.min)

### Plotting BeatAML Cross-validation Results and comparing means difference
treatment_plot_test<-data.frame(Method="ORACLE", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34),treatmentOracle_w34)])
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="ODTMut", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTMut_w34)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="ODTSqrtMut", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTSqrtMut_w34)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="MnLassoMut", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentMnLassoMut_w34)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="LassoMut", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentLassoMut_w34)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="BOSOMut", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentBOSOMut_w34)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="MOM", IC50=treatmentMOM_w34$IC50))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="KRL", IC50=treatmentKRL_w34$IC50)) # CAMBIAR

treatment_plot_test$Method<-factor(treatment_plot_test$Method, levels=c("ORACLE", "MOM","KRL", "BOSOMut", "ODTMut","ODTSqrtMut", "LassoMut", "MnLassoMut"))

treatment_plot_test_gg <- ggplot(treatment_plot_test, aes(x=Method, y=IC50, fill=Method)) +
  geom_boxplot(lwd=0.8,position=position_dodge(0.8)) +
  theme_bw()+
  scale_fill_manual(values = colors)+
  ylab("IC50*") +
  ggtitle("Test in BeatAML (Waves 3+4)")+
  theme(text = element_text(family = "Roboto"),
        axis.text.x = element_text(size = 24, angle=60,vjust=0.65),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 27),
        plot.title = element_text(size = 32, face = "bold"),
        plot.background = element_rect(fill = "transparent", size = 0, color = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        legend.position = "NULL")+
  ylim(0, 5)
treatment_plot_test_gg
ggsave(paste(folder_dir,"/images/boxplot_mutations_test.png",sep=""), treatment_plot_test_gg, width = 8, height = 7, dpi = 1000)

with(treatment_plot_test,                       # Order boxes by median
     reorder(Method,
             IC50,
             median))

with(treatment_plot_test,                       # Order boxes by median
     reorder(Method,
             IC50,
             var))

wilcox.test(treatment_plot_test$IC50[treatment_plot_test$Method=="MnLassoMut"], 
            treatment_plot_test$IC50[treatment_plot_test$Method=="BOSOMut"], alternative="less")

# Intragroup Validation  --------------------------------------------
"
We compared whether the IC50* of a drug in patients in whom it was recommended was lower than the IC50* in patients in whom it was not recommended.
Using this information we compared the sensitivity to a drug for a specific group against the sensitivity to that drug for the rest of the samples by using a 2-tailed Wilcoxon test.
This analysis was performed both for the BeatAML2 dataset (test dataset) and the GDSC AML cell lines cohort (predicted dataset).
"
### Obtaining MOM and KRL predicitions in BeatAML2

# treatmentKRL_All<- read_csv("data/output/KRL_results_all_samples_prediction.csv", show_col_types = FALSE)
# treatmentKRL_All<-treatmentKRL_All[which(treatmentKRL_All$Patients %in% rownames(mutations)),]
# treatmentKRL_All$DrugName[which(!(treatmentKRL_All$DrugName %in% colnames(drug_response_matrix)))]<-NA

### Intragroup Validation in BeatAML

# 1. BeatAML--------------------------------------------------------------------
number_methods <- 8
# number_methods <- nlevels(treatment_plot_test$Method)-1
Validation_Beat <- data.frame(
  patients = rep(rownames(drug_response_w34), number_methods),
  drug = c(
    colnames(drug_response_w34)[treatmentOracle_w34],
    colnames(drug_response_w34)[treatmentODTMut_w34],
    colnames(drug_response_w34)[treatmentODTSqrtMut_w34],
    colnames(drug_response_w34)[treatmentMnLassoMut_w34],
    colnames(drug_response_w34)[treatmentLassoMut_w34],
    colnames(drug_response_w34)[treatmentBOSOMut_w34],
    treatmentMOM_w34$Treatment,
    treatmentKRL_w34$DrugName # CAMBIAR
  ),
  IC50 = c(
    drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentOracle_w34)],
    drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTMut_w34)],
    drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTSqrtMut_w34)],
    drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentMnLassoMut_w34)],
    drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentLassoMut_w34)],
    drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentBOSOMut_w34)],
    treatmentMOM_w34$IC50,
    treatmentKRL_w34$IC50
  ),
  Method = rep(unique(treatment_plot_test$Method), each = nrow(drug_response_w34))
)

# CREATE PLOTS 
methods <- unique(treatment_plot_test$Method)[-1]
colors2 <-  c("#FDE5B0", "#F3908B", "#8E4884", "#E36192", "#96D6B6", "#367592","#39A7AE")

n <- 1 # choose method to plot
gg_drug_plots <- plotIntragroup(methods[n],colors2[n],Validation_Beat,drug_response_w34,"BeatAML2")
gg_drug_plots

# Plot all methods
for (i in 1:length(methods)) {
  gg_drug_plots <- plotIntragroup(methods[i],colors2[i],Validation_Beat,drug_response_w34,"BeatAML2")
}

############################################################
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
sdevs <- rowSds(as.matrix(expression_w12)) 
keep <- which(order(sdevs, decreasing = T) <= ngenes)
expression_w12 <- expression_w12[keep,]
expression_w12 <- t(expression_w12)

expression_w34 <- expression_w34[keep,] # keep same genes in both matrices
expression_w34 <- t(expression_w34)

### Train Models with expression data
library(robustbase)
set.seed(2022)
folds <- 5
groups <- sample(1:folds, nrow(drug_response_w12), replace = T)

treatmentMnLassoExp <- treatmentLassoExp <- treatmentODTExp <- treatmentODTSqrtExp <- rep(NA, length(groups))

for (group in 1:folds) {
  cat("Fold: ", group, "\n")
  dejar <- groups != group
  quitar <- !dejar
  
  # ODT
  cat("ODT\n")
  ODT <- trainTree(expression_w12[dejar,], drug_response_w12[dejar,], minbucket = 20)
  treatmentODTExp[quitar] <- predictTree(ODT,drug_response_w12[quitar,], expression_w12[quitar,])
  
  # ODT sqrt
  cat("ODT sqrt\n")
  ODTSqrt <- trainTree(expression_w12[dejar,], (drug_response_w12[dejar,])^.5, minbucket = 20)
  treatmentODTSqrtExp[quitar] <- predictTree(ODTSqrt,drug_response_w12[quitar,], expression_w12[quitar,])
  
  # Using MnLasso
  cat("Multinomial Lasso\n")
  MnLasso <- trainMnLasso(expression_w12[dejar,], pseudop_w12[dejar,])
  treatmentMnLassoExp[quitar] <- predictMnLasso(MnLasso,drug_response_w12[quitar,],expression_w12[quitar,])
  
  # Using Lasso
  cat("Lasso\n")
  Lasso <- trainLasso(expression_w12[dejar,], drug_response_w12[dejar,])
  treatmentLassoExp[quitar] <- predictMnLasso(Lasso,drug_response_w12[quitar,],expression_w12[quitar,])
  
}

# Use same groups as previously
treatmentBOSOExp <- rep(NA, length(groups))

for (group in 1:folds) {
  dejar <- groups != group
  quitar <- !dejar
  
  BOSO <- trainBOSO(expression_w12[dejar,], drug_response_w12[dejar,],maxVarsBlock = 10, standardize=F, IC = "eBIC")
  treatmentBOSOExp[quitar] <- predictBOSO(BOSO, drug_response_w12[quitar,], expression_w12[quitar,])
}
# save(ODT, ODTSqrt, MnLasso,Lasso, BOSO, file = paste(folder_dir,"Rdata/expression_models.RData",sep=""))
# save(treatmentODTExp, treatmentODTSqrtExp, treatmentMnLassoExp, treatmentLassoExp, treatmentBOSOExp, file = paste(folder_dir,"Rdata/expression_results.RData",sep=""))
load(paste(folder_dir,"/expression_models.RData",sep=""))
load(paste(folder_dir,"/expression_results.RData",sep=""))

treatmentOracle <- apply(drug_response_w12, 1, which.min)

treatments <- cbind(drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentOracle)],
                    drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentODTExp)],
                    drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentODTSqrtExp)],
                    drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentMnLassoExp)],
                    drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentLassoExp)],
                    drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentBOSOExp)],
                    drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentODTMut)],
                    drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentODTSqrtMut)],
                    drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentMnLassoMut)],
                    drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentLassoMut)],
                    drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentBOSOMut)]
)

colnames(treatments) <- c("Oracle","ODTExp","ODTSqrtExp","MnLassoExp","LassoExp","BOSOExp",
                          "ODTMut","ODTSqrtMut","MnLassoMut","LassoMut","BOSOMut")

# Plotting Results.... 
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

plot_Ex_Mut$Method<-factor(plot_Ex_Mut$Method, levels=c("Oracle",
                                                        "BOSOExp","BOSOMut",
                                                        "ODTSqrtExp","ODTSqrtMut",
                                                        "LassoExp","LassoMut",
                                                        "ODTExp","ODTMut",
                                                        "MnLassoExp","MnLassoMut"))
plot_Ex_Mut$General_Method<-factor(plot_Ex_Mut$General_Method, levels=c("Oracle", "BOSO", "ODTSqrt", "Lasso",
                                                                        "ODT","MnLasso"))

mycomparisons<-list(c("ODTExp", "ODTMut"), c("ODTSqrtExp", "ODTSqrtMut"), c("LassoExp", "LassoMut"),
                    c("BOSOExp", "BOSOMut"),
                    c("MnLassoExp", "MnLassoMut"))

colors3 <- c("#A83333", "#96D6B6","#FDE5B0", "#F3908B", "#E36192", "#8E4884") 
gg_ex_mut<-ggplot(plot_Ex_Mut, aes(x=Method, y=IC50, fill=General_Method))+
  geom_boxplot_pattern(
    pattern_color = "black",
    pattern_fill = "black",
    pattern_spacing = 0.02,
    pattern_density      = 0.01,
    aes(pattern = Type))+
  scale_fill_manual(values=colors3)+theme_bw()+
  ggtitle("5-fold CV in BeatAML (Waves 1+2)")+ylab("IC50*")+ 
  labs(fill="General Method")+
  stat_compare_means(aes(Method),vjust=0, hjust=-1, comparisons=mycomparisons,label.y = 5)+
  theme(text = element_text(size = 20, family = "Roboto"),
        axis.text.x = element_text(size = 16, angle=20,vjust=0.7),
        plot.background = element_rect(fill = "transparent", size = 0, color = "transparent"),
        panel.background = element_rect(fill = "transparent")) 
# ylim(0, 6)
gg_ex_mut

ggsave(paste(folder_dir,"/images/boxplot_both.png",sep=""), gg_ex_mut, width = 14, height = 6, dpi = 1000)

wilcox.test(plot_Ex_Mut$IC50[plot_Ex_Mut$Method=="MnLassoExp"], 
            plot_Ex_Mut$IC50[plot_Ex_Mut$Method=="ODTSqrtExp"], alternative="less")

lapply(unique(plot_Ex_Mut$Method), function(method){
  return(c(as.character(method),median(na.omit(plot_Ex_Mut$IC50[plot_Ex_Mut$Method==method]))))
})


# Train with all expression data and test in waves 3+4
# ODTMut
cat("ODTExp\n")
tic()
ODTExp_All <- trainTree(expression_w12, drug_response_w12, minbucket = 10)
treatmentODTExp_w34 <- predictTree(ODTExp_All,drug_response_w34, expression_w34)
toc()

# ODTExp sqrt
cat("ODTExp sqrt\n")
tic()
ODTSqrtExp_All <- trainTree(expression_w12, (drug_response_w12)^.5, minbucket = 10)
treatmentODTSqrtExp_w34 <- predictTree(ODTSqrtExp_All,drug_response_w34, expression_w34)
toc()

# Using MnLasso
cat("Multinomial Lasso Exp\n")
tic()
MnLassoExp_All <- trainMnLasso(expression_w12, pseudop_w12)
treatmentMnLassoExp_w34 <- predictMnLasso(MnLassoExp_All,pseudop_w34,expression_w34)
toc()

# Using Lasso
cat("Lasso Exp\n")
tic()
LassoExp_All <- trainLasso(expression_w12, drug_response_w12)
treatmentLassoExp_w34 <- predictLasso(LassoExp_All,drug_response_w34,expression_w34)
toc()

# Using Bosso
cat("Bosso Exp\n")
tic()
BOSOExp_All <- trainBOSO(expression_w12, drug_response_w12,maxVarsBlock = 10, standardize=F, IC = "eBIC")
treatmentBOSOExp_w34 <- predictBOSO(BOSOExp_All,drug_response_w34,expression_w34)
toc()

treatmentOracle_w34 <- apply(drug_response_w34, 1, which.min)

# save(ODTExp_All, ODTSqrtExp_All, MnLassoExp_All,LassoExp_All, BOSOExp_All, file = paste(folder_dir,"Rdata/expression_models_all.RData",sep=""))
# save(treatmentODTExp_w34, treatmentODTSqrtExp_w34, treatmentMnLassoExp_w34,treatmentLassoExp_w34, treatmentBOSOExp_w34, file = paste(folder_dir,"Rdata/expression_results_test.RData",sep=""))
load(paste(folder_dir,"/expression_models_all.RData",sep=""))
load(paste(folder_dir,"/expression_results_test.RData",sep=""))

### Plotting BeatAML Cross-validation Results and comparing means difference
treatment_plot_test_exp<-data.frame(Method="ORACLE", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34),treatmentOracle_w34)])
treatment_plot_test_exp<-rbind(treatment_plot_test_exp,data.frame(Method="ODTExp", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTExp_w34)]))
treatment_plot_test_exp<-rbind(treatment_plot_test_exp,data.frame(Method="ODTSqrtExp", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTSqrtExp_w34)]))
treatment_plot_test_exp<-rbind(treatment_plot_test_exp,data.frame(Method="MnLassoExp", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentMnLassoExp_w34)]))
treatment_plot_test_exp<-rbind(treatment_plot_test_exp,data.frame(Method="LassoExp", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentLassoExp_w34)]))
treatment_plot_test_exp<-rbind(treatment_plot_test_exp,data.frame(Method="BOSOExp", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentBOSOExp_w34)]))

treatment_plot_test_exp$Method<-factor(treatment_plot_test_exp$Method, levels=c("ORACLE", "BOSOExp", "ODTSqrtExp", "LassoExp", "ODTExp", "MnLassoExp"))

treatment_plot_test_exp_gg <- ggplot(treatment_plot_test_exp, aes(x=Method, y=IC50, fill=Method)) +
  geom_boxplot(lwd=0.8,position=position_dodge(0.8)) +
  theme_bw()+
  scale_fill_manual(values=colors3)+
  ylab("IC50*") +
  ggtitle("Test in BeatAML (Waves 3+4)")+
  theme(text = element_text(size = 20, family = "Roboto"),
        axis.text.x = element_text(size = 16, angle=20,vjust=0.7),
        plot.background = element_rect(fill = "transparent", size = 0, color = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  ylim(0, 4)
treatment_plot_test_exp_gg
ggsave(paste(folder_dir,"/images/boxplot_expression_test.png",sep=""), treatment_plot_test_exp_gg, width = 10, height = 6, dpi = 1000)


wilcox.test(treatment_plot_test_exp$IC50[treatment_plot_test_exp$Method=="MnLassoExp"], 
            treatment_plot_test_exp$IC50[treatment_plot_test_exp$Method=="ODTExp"], alternative="less")

lapply(unique(treatment_plot_test_exp$Method), function(method){
  return(c(as.character(method),median(na.omit(treatment_plot_test_exp$IC50[treatment_plot_test_exp$Method==method]))))
})

#############################################################
# Multiomics in BeatAML: Expression and Mutation models
#############################################################

# Join expression with mutations 
exp_mut_w12 <- merge(expression_w12, mutations_w12, by = "row.names", by.row = TRUE)
rownames(exp_mut_w12) <- exp_mut_w12$Row.names
exp_mut_w12$Row.names <- NULL
exp_mut_w12 <- data.matrix(exp_mut_w12, rownames.force = TRUE)

exp_mut_w34 <- merge(expression_w34, mutations_w34, by = "row.names", by.row = TRUE)
rownames(exp_mut_w34) <- exp_mut_w34$Row.names
exp_mut_w34$Row.names <- NULL
exp_mut_w34 <- data.matrix(exp_mut_w34, rownames.force = TRUE)

exp_mut_w12 <- as.matrix(exp_mut_w12)
exp_mut_w34 <- as.matrix(exp_mut_w34)
exp_mut_w12 <- exp_mut_w12[rownames(drug_response_w12),]
exp_mut_w34 <- exp_mut_w34[rownames(drug_response_w34),]

folds <- 5
groups <- sample(1:folds, nrow(drug_response_w12), replace = T)
treatmentMnLassoBoth <- treatmentLassoBoth <- treatmentODTBoth <- treatmentODTSqrtBoth <- rep(NA, length(groups))

# Train Models with EXPRESSION AND MUTATIONS
for (group in 1:folds) {
  cat("Fold: ", group, "\n")
  dejar <- groups != group
  quitar <- !dejar
  
  # ODT
  cat("ODT\n")
  ODTBoth <- trainTree(exp_mut_w12[dejar,], drug_response_w12[dejar,], minbucket = 20)
  treatmentODTBoth[quitar] <- predictTree(ODTBoth,drug_response_w12[quitar,], exp_mut_w12[quitar,])
  
  # ODT sqrt
  cat("ODT sqrt\n")
  ODTSqrtBoth <- trainTree(exp_mut_w12[dejar,], (drug_response_w12[dejar,])^.5, minbucket = 20)
  treatmentODTSqrtBoth[quitar] <- predictTree(ODTSqrtBoth,drug_response_w12[quitar,], exp_mut_w12[quitar,])
  
  # Using MnLasso
  cat("Multinomial Lasso\n")
  MnLassoBoth <- trainMnLasso(exp_mut_w12[dejar,], pseudop_w12[dejar,])
  treatmentMnLassoBoth[quitar] <- predictMnLasso(MnLassoBoth,drug_response_w12[quitar,],exp_mut_w12[quitar,])
  
  # Using Lasso
  cat("Lasso\n")
  LassoBoth <- trainLasso(exp_mut_w12[dejar,], drug_response_w12[dejar,])
  treatmentLassoBoth[quitar] <- predictMnLasso(LassoBoth,drug_response_w12[quitar,],exp_mut_w12[quitar,])
}

set.seed(2021)
folds <- 5
groups <- sample(1:folds, length(groups), replace = T)

# Use same groups as previously
treatmentBOSOBoth <- rep(NA, length(groups))

for (group in 1:folds) {
  dejar <- groups != group
  quitar <- !dejar
  
  BOSO <- trainBOSO(exp_mut_w12[dejar,], drug_response_w12[dejar,],maxVarsBlock = 10, standardize=F, IC = "eBIC")
  treatmentBOSOBoth[quitar] <- predictBOSO(BOSO,drug_response_w12[quitar,], exp_mut_w12[quitar,])
}
# save(ODTBoth, ODTSqrtBoth, MnLassoBoth,LassoBoth, BOSOBoth, file = paste(folder_dir,"Rdata/both_models.RData",sep=""))
# save(treatmentODTBoth, treatmentODTSqrtBoth, treatmentMnLassoBoth, treatmentLassoBoth, treatmentBOSOBoth, file = paste(folder_dir,"Rdata/both_results.RData",sep=""))
load(paste(folder_dir,"/both_models.RData",sep=""))
load(paste(folder_dir,"/both_results.RData",sep=""))

treatmentOracle <- apply(drug_response_w12, 1, which.min)

# Plot all results!!!

treatments2 <- cbind(
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentOracle)],
  treatmentMOM$IC50,
  KRL_results_cv$IC50,
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentODTExp)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentODTSqrtExp)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentMnLassoExp)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentLassoExp)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentBOSOExp)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentODTMut)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentODTSqrtMut)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentMnLassoMut)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentLassoMut)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentBOSOMut)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentODTBoth)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentODTSqrtBoth)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentMnLassoBoth)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentLassoBoth)],
  drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentBOSOBoth)]
)

colnames(treatments2) <- c(
  "Oracle", "MOM","KRL","ODTExp", "ODTSqrtExp", "MnLassoExp", "LassoExp", "BOSOExp",
  "ODTMut", "ODTSqrtMut", "MnLassoMut", "LassoMut", "BOSOMut",
  "ODTBoth", "ODTSqrtBoth", "MnLassoBoth", "LassoBoth", "BOSOBoth"
)

mycomparisons2 <- list(
  c("ODTExp", "ODTMut"), c("ODTMut", "ODTBoth"),
  c("ODTSqrtExp", "ODTSqrtMut"),c("ODTSqrtMut", "ODTSqrtBoth"),
  c("LassoExp", "LassoMut"),c("LassoMut", "LassoBoth"),
  c("BOSOExp", "BOSOMut"),
  c("BOSOMut", "BOSOBoth"),
  c("MnLassoExp", "MnLassoMut"),c("MnLassoMut", "MnLassoBoth"))

mycomparisons3 <- list(
  c("ODTExp", "ODTBoth"),
  c("ODTSqrtExp", "ODTSqrtBoth"),
  c("LassoExp", "LassoBoth"),
  c("BOSOExp", "BOSOBoth"),
  c("MnLassoExp", "MnLassoBoth"))

# Plotting Results.... 
plot_all<-data.frame(Method=rep(colnames(treatments2), each=nrow(treatments2)))
plot_all$IC50<-unlist(as.data.frame(treatments2))
plot_all$Type<-"GE"
plot_all$Type[grep("Mut", plot_all$Method)]<-"Mut" # add type "Mut" to mutation values
plot_all$Type[grep("Both", plot_all$Method)]<-"Both"
plot_all$Type[ which(plot_all$Method == "MOM")] <- "Mut"
plot_all$Type[ which(plot_all$Method == "KRL")] <- "Mut"

plot_all<-plot_all[order(plot_all$Method),]
plot_all$General_Method<-"unknown"
plot_all$General_Method[grep("Oracle", plot_all$Method)]<-"Oracle"
plot_all$General_Method[grep("MOM", plot_all$Method)]<-"MOM"
plot_all$General_Method[grep("KRL", plot_all$Method)]<-"KRL"
plot_all$General_Method[grep("BOSO", plot_all$Method)]<-"BOSO"
plot_all$General_Method[grep("ODT", plot_all$Method)]<-"ODT"
plot_all$General_Method[grep("ODTSqrt", plot_all$Method)]<-"ODTSqrt"
plot_all$General_Method[grep("Lasso", plot_all$Method)]<-"Lasso"
plot_all$General_Method[grep("MnLasso", plot_all$Method)]<-"MnLasso"

plot_all$Method<-factor(plot_all$Method, levels=c("Oracle","MOM","KRL","BOSOExp","BOSOMut","BOSOBoth",
                                                  "ODTExp","ODTMut","ODTBoth",
                                                  "ODTSqrtExp","ODTSqrtMut", "ODTSqrtBoth",
                                                  "LassoExp","LassoMut","LassoBoth",
                                                  "MnLassoExp","MnLassoMut","MnLassoBoth"))

plot_all$General_Method<-factor(plot_all$General_Method, levels=c("Oracle", "MOM","KRL","BOSO", "ODT","ODTSqrt", "Lasso", "MnLasso"))

plot_all$Type<-factor(plot_all$Type, levels=c("GE", "Mut", "Both"))

colors <- c("#A83333","#367592", "#39A7AE", "#96D6B6", "#FDE5B0", "#F3908B", "#E36192", "#8E4884")

gg_ex_both<-ggplot(plot_all, aes(x=Method, y=IC50, fill=General_Method))+    # Different pattern for each group
  geom_boxplot_pattern(
    pattern_color = "black",
    pattern_fill = "black",
    pattern_spacing = 0.02,
    pattern_density = 0.01,
    aes(pattern = Type))+
  scale_fill_manual(values=colors)+theme_bw()+ggtitle("5-fold CV in BeatAML (Waves 1+2)")+ylab("IC50*")+
  labs(fill="General Method")+
  stat_compare_means(aes(Method),vjust=0, hjust=-1, comparisons=mycomparisons2,label.y = 5)+
  stat_compare_means(aes(Method),vjust=0, hjust=-1, comparisons=mycomparisons3,label.y = 5.4)+
  theme(text = element_text(size = 20, family = "Roboto"),
        axis.text.x = element_text(size = 16, angle=20,vjust=0.7),
        plot.background = element_rect(fill = "transparent", size = 0, color = "transparent"),
        panel.background = element_rect(fill = "transparent"))
gg_ex_both
ggsave(paste(folder_dir,"/images/boxplot_all.png",sep=""), gg_ex_both, width = 17, height = 7, dpi = 1000)

# Train with BOTH expression AND mutations data and test in waves 3+4
# ODTMut
cat("ODTBoth\n")
tic()
ODTBoth_All <- trainTree(exp_mut_w12, drug_response_w12, minbucket = 10)
treatmentODTBoth_w34 <- predictTree(ODTBoth_All,drug_response_w34, exp_mut_w34)
toc()

# ODTBoth sqrt
cat("ODTBoth sqrt\n")
tic()
ODTSqrtBoth_All <- trainTree(exp_mut_w12, (drug_response_w12)^.5, minbucket = 10)
treatmentODTSqrtBoth_w34 <- predictTree(ODTSqrtBoth_All,drug_response_w34, exp_mut_w34)
toc()

# Using MnLasso
cat("Multinomial Lasso Both\n")
tic()
MnLassoBoth_All <- trainMnLasso(exp_mut_w12, pseudop_w12)
treatmentMnLassoBoth_w34 <- predictMnLasso(MnLassoBoth_All,pseudop_w34,exp_mut_w34)
toc()

# Using Lasso
cat("Lasso Both\n")
tic()
LassoBoth_All <- trainLasso(exp_mut_w12, drug_response_w12)
treatmentLassoBoth_w34 <- predictLasso(LassoBoth_All,drug_response_w34,exp_mut_w34)
toc()

# Using Bosso
cat("Bosso Both\n")
tic()
BOSOBoth_All <- trainBOSO(exp_mut_w12, drug_response_w12,maxVarsBlock = 10, standardize=F, IC = "eBIC")
treatmentBOSOBoth_w34 <- predictBOSO(BOSOBoth_All,drug_response_w34,exp_mut_w34)
toc()

treatmentOracle_w34 <- apply(drug_response_w34, 1, which.min)

# save(ODTBoth_All, ODTSqrtBoth_All, MnLassoBoth_All,LassoBoth_All, BOSOBoth_All, file = paste(folder_dir,"Rdata/both_models_all.RData",sep=""))
# save(treatmentODTBoth_w34, treatmentODTSqrtBoth_w34, treatmentMnLassoBoth_w34, treatmentLassoBoth_w34, treatmentBOSOBoth_w34, file = paste(folder_dir,"Rdata/both_results_test.RData",sep=""))
load(paste(folder_dir,"/both_models_all.RData",sep=""))
load(paste(folder_dir,"/both_results_test.RData",sep=""))

### Plotting BeatAML Cross-validation Results and comparing means difference
treatment_plot_test_both<-data.frame(Method="ORACLE", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34),treatmentOracle_w34)])
treatment_plot_test_both<-rbind(treatment_plot_test_both,data.frame(Method="ODTBoth", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTBoth_w34)]))
treatment_plot_test_both<-rbind(treatment_plot_test_both,data.frame(Method="ODTSqrtBoth", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTSqrtBoth_w34)]))
treatment_plot_test_both<-rbind(treatment_plot_test_both,data.frame(Method="MnLassoBoth", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentMnLassoBoth_w34)]))
treatment_plot_test_both<-rbind(treatment_plot_test_both,data.frame(Method="LassoBoth", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentLassoBoth_w34)]))
treatment_plot_test_both<-rbind(treatment_plot_test_both,data.frame(Method="BOSOBoth", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentBOSOBoth_w34)]))

treatment_plot_test_both$Method<-factor(treatment_plot_test_both$Method, levels=c("ORACLE", "BOSOBoth", "ODTBoth", "ODTSqrtBoth", "LassoBoth", "MnLassoBoth"))

treatment_plot_test_both_gg <- ggplot(treatment_plot_test_both, aes(x=Method, y=IC50, fill=Method)) +
  geom_boxplot(lwd=0.8,position=position_dodge(0.8)) +
  theme_bw()+
  scale_fill_manual(values=colors3)+ylab("IC50*") +
  ggtitle("Test in BeatAML (Waves 3+4)")+
  theme(text = element_text(size = 20, family = "Roboto"),
        axis.text.x = element_text(size = 16, angle=20,vjust=0.7),
        plot.background = element_rect(fill = "transparent", size = 0, color = "transparent"),
        panel.background = element_rect(fill = "transparent"))
treatment_plot_test_both_gg
ggsave(paste(folder_dir,"/images/boxplot_both_test.png",sep=""), treatment_plot_test_both_gg, width = 18, height = 7, dpi = 1000)

# Plot all TEST results!!! WAVES 3+4
treatmentOracle3 <- apply(drug_response_w34, 1, which.min)

treatments3 <- cbind(
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentOracle3)],
  treatmentMOM_w34$IC50,
  treatmentKRL_w34$IC50, 
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTExp_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTSqrtExp_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentMnLassoExp_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentLassoExp_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentBOSOExp_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTMut_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTSqrtMut_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentMnLassoMut_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentLassoMut_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentBOSOMut_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTBoth_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTSqrtBoth_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentMnLassoBoth_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentLassoBoth_w34)],
  drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentBOSOBoth_w34)]
)

colnames(treatments3) <- c(
  "Oracle", "MOM","KRL","ODTExp", "ODTSqrtExp", "MnLassoExp", "LassoExp", "BOSOExp",
  "ODTMut", "ODTSqrtMut", "MnLassoMut", "LassoMut", "BOSOMut",
  "ODTBoth", "ODTSqrtBoth", "MnLassoBoth", "LassoBoth", "BOSOBoth"
)

# Plotting Results.... 
plot_all<-data.frame(Method=rep(colnames(treatments3), each=nrow(treatments3)))
plot_all$IC50<-unlist(as.data.frame(treatments3))
plot_all$Type<-"GE"
plot_all$Type[grep("Mut", plot_all$Method)]<-"Mut" # add type "Mut" to mutation values
plot_all$Type[grep("Both", plot_all$Method)]<-"Both"
plot_all$Type[ which(plot_all$Method == "MOM")] <- "Mut"
plot_all$Type[ which(plot_all$Method == "KRL")] <- "Mut"

plot_all<-plot_all[order(plot_all$Method),]
plot_all$General_Method<-"unknown"
plot_all$General_Method[grep("Oracle", plot_all$Method)]<-"Oracle"
plot_all$General_Method[grep("BOSO", plot_all$Method)]<-"BOSO"
plot_all$General_Method[grep("MOM", plot_all$Method)]<-"MOM"
plot_all$General_Method[grep("KRL", plot_all$Method)]<-"KRL"
plot_all$General_Method[grep("ODT", plot_all$Method)]<-"ODT"
plot_all$General_Method[grep("ODTSqrt", plot_all$Method)]<-"ODTSqrt"
plot_all$General_Method[grep("Lasso", plot_all$Method)]<-"Lasso"
plot_all$General_Method[grep("MnLasso", plot_all$Method)]<-"MnLasso"

plot_all$Method<-factor(plot_all$Method, levels=c("Oracle","MOM","KRL","BOSOExp","BOSOMut","BOSOBoth",
                                                  "ODTExp","ODTMut","ODTBoth",
                                                  "ODTSqrtExp","ODTSqrtMut", "ODTSqrtBoth",
                                                  "LassoExp","LassoMut","LassoBoth",
                                                  "MnLassoExp","MnLassoMut","MnLassoBoth"))
plot_all$General_Method<-factor(plot_all$General_Method, levels=c("Oracle", "MOM","KRL","BOSO", "ODT","ODTSqrt", "Lasso", "MnLasso"))
plot_all$Type<-factor(plot_all$Type, levels=c("GE", "Mut", "Both"))

gg_ex_both<-ggplot(plot_all, aes(x=Method, y=IC50, fill=General_Method))+    # Different pattern for each group
  geom_boxplot_pattern(
    pattern_color = "black",
    pattern_fill = "black",
    pattern_spacing = 0.02,
    pattern_density      = 0.01,
    aes(pattern = Type))+
  
  scale_fill_manual(values=colors)+theme_bw()+ggtitle("Test in BeatAML (Waves 3+4)")+ylab("IC50*")+
  labs(fill="General Method")+
  stat_compare_means(aes(Method),vjust=0, hjust=-1, comparisons=mycomparisons2,label.y = 5)+
  stat_compare_means(aes(Method),vjust=0, hjust=-1, comparisons=mycomparisons3,label.y = 5.4)+
  #scale_alpha_manual(values=c( 0.3,0.6,1))+
  theme(text = element_text(size = 20, family = "Roboto"),
        axis.text.x = element_text(size = 16, angle=20,vjust=0.7),
        plot.background = element_rect(fill = "transparent", size = 0, color = "transparent" ),
        panel.background = element_rect(fill = "transparent"))
gg_ex_both
ggsave(paste(folder_dir,"/images/boxplot_all_test.png",sep=""), gg_ex_both, width = 17, height = 7, dpi = 1000)
