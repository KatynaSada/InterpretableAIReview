
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
library(impute)

source("./code/getBMeffect_F.R")
source("./code/Get_Treatment_OptimizationMatrix.R")
source("./code/MOM_MILP.R")
source("./code/Extact_results.R")
source("./code/CV_Prediction.R")

# Load mutational data----------------------------------------------------------

gene_variants<-read_excel("data/input/41586_2018_623_MOESM3_ESM.xlsx",
                          sheet="Table S7-Variants for Analysis")
clinical <- read_excel("data/input/41586_2018_623_MOESM3_ESM.xlsx",
                       sheet = "Tabe S5-Clinical Summary")

# Load drug sensitivity data----------------------------------------------------

Drug_response<-read_excel("data/input/41586_2018_623_MOESM3_ESM.xlsx",
                          sheet="Table S10-Drug Responses")




# Pre-processing Step-----------------------------------------------------------

# Build Drug Matrix

Drug<-matrix(data="NA",nrow=length(unique(Drug_response$inhibitor)), 
             ncol=length(unique(Drug_response$lab_id)))

colnames(Drug)<-as.character(unique(Drug_response$lab_id))
rownames(Drug)<-unique(Drug_response$inhibitor)


for(i in as.character(unique(Drug_response$lab_id))){
  ind<-which(Drug_response$lab_id==i)
  D<-as.character(Drug_response$inhibitor[ind])
  Drug[D, i]<-Drug_response$ic50[ind]
}

#Identifying missing values

is.na(Drug) <- Drug == "NA"

#change number formatting

for(i in 2:(dim(Drug)[2]-1)){
  Drug[,i]<-as.numeric(as.character(Drug[,i]))
}

# Building Mutation Matrix......................................................
Mutations<-matrix(0, nrow=length(unique(gene_variants$labId)), 
                  ncol=length(unique(gene_variants$symbol)))

rownames(Mutations)<-as.character(unique(gene_variants$labId))
colnames(Mutations)<-as.character(unique(gene_variants$symbol))


for(i in as.character(unique(gene_variants$labId))){
  jj<-unique(as.character(gene_variants$symbol[which(gene_variants$labId==i)]))
  Mutations[i, jj]<-1  
}



# Add Translocations
translocations_clinical<-clinical[,c("LabId","specificDxAtAcquisition")][grep("AML with inv", clinical$specificDxAtAcquisition),]
translocations_clinical<-rbind(translocations_clinical,clinical[,c("LabId","specificDxAtAcquisition")][grep("AML with t", clinical$specificDxAtAcquisition),])
trans<-as.character(unique(translocations_clinical$specificDxAtAcquisition))
add<-matrix(0, nrow=nrow(Mutations), ncol=length(trans))
colnames(add)<-trans
rownames(add)<-rownames(Mutations)
for(j in colnames(add)){
p<-translocations_clinical$LabId[which(translocations_clinical$specificDxAtAcquisition ==j)]
p<-p[which(p %in% rownames(add))]
add[p,j]<-1
}


trans<-sapply(trans, function(X){unlist(strsplit(X, split = "; "))[2]})
colnames(add)<-trans
Mutations<-cbind(Mutations,add)
print("Processing Mutations and Drugs Matrices")
Drug_2<-as.matrix(t(Drug))
Mutations<-as.matrix(Mutations)
patients<-rownames(Mutations)[which(rownames(Mutations) %in% rownames(Drug_2))]
mut<-Mutations[patients,colSums(Mutations[patients, ])>0]
mut_NA<-mut[, colSums(mut)>0.01*length(patients)]
sum(colSums(mut_NA)<1)




drug_NA<-Drug_2[patients,]
drug_NA<-apply(drug_NA,2,as.numeric)
rownames(drug_NA)<-rownames(Drug_2[patients,])

drug_NA<-log10(drug_NA)-1
drug_NA <- t(t(drug_NA) - colMeans(drug_NA, na.rm=T))
maximum_drug_correction<-max(drug_NA,na.rm = T)
drug_NA<-drug_NA-maximum_drug_correction
drug_NA<-drug_NA*(-1)

identical(rownames(mut_NA), rownames(Drug_2[patients,]))
identical(rownames(mut_NA), rownames(drug_NA))

#### LOAD original drug and mut2----------------------------------------------

gene_variants<-read_excel("data/input/41586_2018_623_MOESM3_ESM.xlsx",
                          sheet="Table S7-Variants for Analysis")
clinical <- read_excel("data/input/41586_2018_623_MOESM3_ESM.xlsx",
                       sheet = "Tabe S5-Clinical Summary")

# Load drug sensitivity data

Drug_response<-read_excel("data/input/41586_2018_623_MOESM3_ESM.xlsx",
                          sheet="Table S10-Drug Responses")




# Pre-processing Step

# Build Drug Matrix

Drug<-matrix(data="NA",nrow=length(unique(Drug_response$inhibitor)), 
             ncol=length(unique(Drug_response$lab_id)))

colnames(Drug)<-as.character(unique(Drug_response$lab_id))
rownames(Drug)<-unique(Drug_response$inhibitor)


for(i in as.character(unique(Drug_response$lab_id))){
  ind<-which(Drug_response$lab_id==i)
  D<-as.character(Drug_response$inhibitor[ind])
  Drug[D, i]<-Drug_response$ic50[ind]
}

#Identifying missing values

is.na(Drug) <- Drug == "NA"

#change number formatting

for(i in 2:(dim(Drug)[2]-1)){
  Drug[,i]<-as.numeric(as.character(Drug[,i]))
}

# Building Mutation Matrix......................................................
Mutations<-matrix(0, nrow=length(unique(gene_variants$labId)), 
                  ncol=length(unique(gene_variants$symbol)))

rownames(Mutations)<-as.character(unique(gene_variants$labId))
colnames(Mutations)<-as.character(unique(gene_variants$symbol))


for(i in as.character(unique(gene_variants$labId))){
  jj<-unique(as.character(gene_variants$symbol[which(gene_variants$labId==i)]))
  Mutations[i, jj]<-1  
}

# Add Translocations 

translocations_clinical<-clinical[,c("LabId","specificDxAtAcquisition")][grep("AML with inv", clinical$specificDxAtAcquisition),]

translocations_clinical<-rbind(translocations_clinical,clinical[,c("LabId","specificDxAtAcquisition")][grep("AML with t", clinical$specificDxAtAcquisition),])

trans<-as.character(unique(translocations_clinical$specificDxAtAcquisition))
add<-matrix(0, nrow=nrow(Mutations), ncol=length(trans))
colnames(add)<-trans
rownames(add)<-rownames(Mutations)

for(j in colnames(add)){
  p<-translocations_clinical$LabId[which(translocations_clinical$specificDxAtAcquisition ==j)]
  p<-p[which(p %in% rownames(add))]
  add[p,j]<-1
}

trans<-sapply(trans, function(X){unlist(strsplit(X, split = "; "))[2]})
colnames(add)<-trans
Mutations<-cbind(Mutations,add)

# Combining Matrices and Imputation
print("Processing Mutations and Drugs Matrices")
Drug_2<-as.matrix(t(Drug))
Mutations<-as.matrix(Mutations)

# impute missing values

if (sum(is.na(Drug_2))>0){
  Drug_3<-Drug_2[-which((rowSums(is.na(Drug_2))/dim(Drug_2)[2])>0.8),
                 -which((colSums(is.na(Drug_2))/dim(Drug_2)[1])>0.7)]
  Out<-impute.knn(as.matrix(Drug_3), k = 10, rowmax = 0.8, 
                  colmax = 0.75, maxp = 1500, rng.seed=362436069)
  Drug_4<-Out$data[, -which(colnames(Out$data) %in% c("Elesclomol","JNJ-7706621"))]
  
}else{
  Drug_4<-Drug_2[, -which(colnames(Drug_2) %in% c("Elesclomol","JNJ-7706621"))]
}


# Combining Drug and mutational Data
patients<-rownames(Mutations)[which(rownames(Mutations) %in% rownames(Drug_4))]

mut<-Mutations[patients,colSums(Mutations[patients, ])>0]
mut2<-mut[, colSums(mut)>0.01*length(patients)]

sum(colSums(mut2)<1)
rm(mut)

drug<-Drug_4[patients,]
drug<-log10(drug)-1
drug <- t(t(drug) - colMeans(drug, na.rm=T))
maximum_drug_correction<-max(drug,na.rm = T)
drug<-drug-maximum_drug_correction

#### JOIN with drug_NA and mut_NA-----------------------------------------------

drug_NA<-drug_NA[rownames(drug), colnames(drug)]
mut_NA<-mut_NA[rownames(drug), colnames(mut2)]


identical(rownames(mut_NA), rownames(drug_NA))

dim(drug_NA)
dim(mut_NA)

write.csv(drug_NA, file="../../data/input/drug_XAI_forKRL_NAs.csv")
write.csv(mut_NA, file="../../data/input/mut_XAI_forKRL_NAs.csv")
