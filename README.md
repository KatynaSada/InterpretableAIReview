
# **Precision Oncology: a review to assess interpretability in several explainable methods**
_Authors: Marian Gimeno, Katyna Sada and Angel Rubio_

_Date: 21-12-2022_



## Loading Libraries and Packages

The first part is to load all the libraries and external functions that are required for the analysis.

```{r}

library(readxl)
library(RColorBrewer)
library(matrixStats)
library(partykit)
library(glmnet)
library(BOSO)


source("Code_Analysis/XAIfunctions_NV.R")

library(impute)
library(knitr)
library(reticulate)
library(readxl)
library(readr)
library(limma)
library(ggplot2)
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

```




##  Load Data

Then, we will load all the data for the analysis.
Data selected was the BeatAML cohort for which we had mutational data, drug sensitivity data and gene expression data.

```{r}
# Load mutational data

gene_variants<-read_excel("./data/input/41586_2018_623_MOESM3_ESM.xlsx",
                          sheet="Table S7-Variants for Analysis")
clinical <- read_excel("./data/input/41586_2018_623_MOESM3_ESM.xlsx",
                       sheet = "Tabe S5-Clinical Summary")

# Load drug sensitivity data

Drug_response<-read_excel("./data/input/41586_2018_623_MOESM3_ESM.xlsx",
                          sheet="Table S10-Drug Responses")

# Load gene expression data

Expression <- read_excel("./data/input/41586_2018_623_MOESM3_ESM.xlsx",
                         sheet = "Table S8-Gene Counts RPKM")


```

## Pre-processing drug and mutational information 
Now it is time to pre-process the data in order to obtain two different matrices, a matrix containing all patients and their mutations available and a second matrix containing all patients and their sensitivity to the different drugs.
For doing so, we will impute the missing values in the drug matrix.

<details><summary>Click to expand</summary>


### Building the drug matrix

```{r}
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

```
### Building the mutations matrix

```{r}
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

```

### Unifying imputed Drug and Mutation matrices with common patients

```{r}
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


# Tidying up the variables

rm(list = c("Drug_2", "Drug_3", "Drug_4", "Drug_response", "Mutations", "gene_variants"))
rm(list = c("Out", "translocations_clinical", "clinical", "add", "Drug"))
```

### Pre-process gene expresison information

```{r}
# Keep genes with the largest variance
ngenes <- 1000
Expression <- as.data.frame(Expression)
rownames(Expression) <- Expression[,2]
Expression <- Expression[,-c(1,2)]
OldExpression <- as.matrix(Expression)
sdevs <- rowSds(OldExpression)
# sdevs <- rowMedians(OldExpression)

Keep <- which(order(sdevs, decreasing = T) <= ngenes)
Expression <- Expression[Keep,]

commonPatients <- intersect(rownames(drug), colnames(Expression))

```

### Create variables for model training

```{r}
C <- drug - min(drug)
X <- mut2

Y <- exp(-.5*C) # Convert weights into pseudoprobabilities for Multinomial
Y <- Y /rowSums(Y) # Probabilities sum up one
```

</details>





## Comparison of methods with Cross-Validation

We performed a 5-fold cross-validation using the BeatAML dataset.
We trained all models with genetic variants data from 319 patients, dividing the cohort between the training samples 4-folds and testing samples the selected 1-fold.
Each of the folds were tested, and the predicted IC50\* for the 5-fold testing was compared for all the methods and compared against the Oracle -the drug with the optimum IC50\*.

We calculated the Oracle as the minimum IC50\* value for each patient.

Results from MOM and KRL were directly loaded into the workspace due to the computing limitations in R.

<details><summary>Click to expand</summary>

```{r}
set.seed(2022)
Folds <- 5
Groups <- sample(1:Folds, nrow(C), replace = T)
TratamientoMnLassoMut <- TratamientoLassoMut <- TratamientoTreeMut <- TratamientoTreeMut1 <- rep(NA, length(Groups))

for (group in 1:Folds) {
  cat("Fold: ", group, "\n")
  Dejar <- Groups != group
  Quitar <- !Dejar
  
  # TreesMut
  cat("TreesMut\n")
  treeMut <- trainTreeMut(X[Dejar,], C[Dejar,], minbucket = 10)
  TratamientoTreeMut[Quitar] <- predictTreeMut(treeMut,C[Quitar,], X[Quitar,])

  # TreesMut sqrt
  cat("TreesMut sqrt\n")
  treeMut1 <- trainTreeMut(X[Dejar,], (C[Dejar,])^.5, minbucket = 10)
  TratamientoTreeMut1[Quitar] <- predictTreeMut(treeMut1,C[Quitar,], X[Quitar,])

  # Using MnLasso
  cat("Multinomial Lasso Mut\n")
  MnLassoMut <- trainMnLasso(X[Dejar,], Y[Dejar,])
  TratamientoMnLassoMut[Quitar] <- predictMnLasso(MnLassoMut,X[Quitar,])

  # Using Lasso
  cat("Lasso Mut\n")
  LassoMut <- trainLasso(X[Dejar,], C[Dejar,])
  TratamientoLassoMut[Quitar] <- predictLasso(LassoMut,X[Quitar,])

}

# Use same groups as previously
TratamientoBOSOMut <- rep(NA, length(Groups))

for (group in 1:Folds) {
  cat("Fold: ", group, "\n")
  Dejar <- Groups != group
  Quitar <- !Dejar
  BOSOMut <- trainBOSO(X[Dejar,], C[Dejar,],maxVarsBlock = 10, standardize=F, IC = "eBIC")
  TratamientoBOSOMut[Quitar] <- predictBOSO(BOSOMut, X[Quitar,])
}

load("./data/output/mom_COMPARISON.RData")

# KRL
KRL_results_cv <- read_csv("data/output/KRL_results_5cv_samples_prediction.csv")
KRL_all<-read_csv("data/output/KRL_results_all_samples_prediction.csv")

KRL_results_cv<-KRL_results_cv[which(KRL_results_cv$Patients %in% rownames(mut2)),]
KRL_all<-KRL_all[which(KRL_all$Patients %in% rownames(mut2)),]

KRL_results_cv<-KRL_results_cv[which(KRL_results_cv$DrugName %in% colnames(drug)),]
KRL_all<-KRL_all[which(KRL_all$DrugName %in% colnames(drug)),]


KRL_all$`Best Drug`<-KRL_all$`Best Drug`+1
KRL_results_cv$Drug<-KRL_results_cv$Drug+1


KRL_all$IC50<-C[cbind(KRL_all$Patients, KRL_all$DrugName)]
KRL_results_cv$IC50<-C[cbind(KRL_results_cv$Patients, KRL_results_cv$DrugName)]


TratamientoOracle <- apply(C, 1, which.min)
```

### Plotting BeatAML Cross-validation Results and comparing means difference

```{r}

Tratamiento_plot<-data.frame(Method="ORACLE", IC50=C[cbind(1:nrow(C),TratamientoOracle)])
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="TreeMut", IC50=C[cbind(1:nrow(C), TratamientoTreeMut)]))
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="TreeMutSqrt", IC50=C[cbind(1:nrow(C), TratamientoTreeMut1)]))
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="Multinomial", IC50=C[cbind(1:nrow(C), TratamientoMnLassoMut)]))
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="LassoMut", IC50=C[cbind(1:nrow(C), TratamientoLassoMut)]))
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="BOSOMut", IC50=C[cbind(1:nrow(C), TratamientoBOSOMut)]))
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="MOM", IC50=MOM$IC50))
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="KRL", IC50=KRL_results_cv$IC50))


Tratamiento_plot$Method<-factor(Tratamiento_plot$Method, levels=c("ORACLE", "MOM", "BOSOMut", "TreeMutSqrt", "LassoMut", "TreeMut", "Multinomial","KRL"))

ggplot(Tratamiento_plot, aes(x=Method, y=IC50, fill=Method))+geom_violin()+
  theme_bw()+scale_fill_npg()+ylab("IC50*")+ggtitle("BeatAML 5-fold Cross Validation")

ggplot(Tratamiento_plot, aes(x=Method, y=IC50, fill=Method))+geom_boxplot()+
  theme_bw()+scale_fill_npg()+ylab("IC50*")+theme(text = element_text(size = 15)) 


wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="MOM"], alternative="less")
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="BOSOMut"], alternative="less")
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="TreeMutSqrt"], alternative="less")
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="LassoMut"], alternative="less")
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="TreeMut"], alternative="less")
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="Multinomial"], alternative="less")
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="KRL"], alternative="less")

```


</details>


## Independent Cohort Validation

One of the main challenges of Machine Learning, including Precision Medicine, is generalization, i.e. the ability to adapt to new, previously unseen data.
All the methods were tested on the GDSC AML dataset to check their generalization ability.
The models were trained using the BeatAML dataset and were used to predict the optimal drug for AML cell lines from GDSC using its mutation files.
Each of the cell lines was recommended a drug, we compared the all-samples IC50 for all the models and against the Oracle (the drug with the minimum IC50 for each cell line).

<details><summary>Click to expand</summary>
### Training models in complete BeatAML data and timing

We trained the models in BeatAML cohort and measured the training time for each model except KRL and MOM, that were meassured independently.
Technical implementation refers to the computational burden and software that the method requires.Despite it could be considered less important, some of the algorithms require hours of computing time for the BeatAML of the subset of AML samples in GDSC -that be considered to be small/medium size.

```{r}
# Train with all the data

# TreesMut
cat("TreesMut\n")
tic()
treeMut_All <- trainTreeMut(X, C, minbucket = 10)
TratamientoTreeMut_ALL <- predictTreeMut(treeMut_All,C, X)
toc()

# TreesMut sqrt
cat("TreesMut sqrt\n")
tic()
treeMut1_All <- trainTreeMut(X, (C)^.5, minbucket = 10)
TratamientoTreeMut1_All <- predictTreeMut(treeMut1_All,C, X)
toc()

# Using MnLasso
cat("Multinomial Lasso Mut\n")
tic()
MnLassoMut_All <- trainMnLasso(X, Y)
TratamientoMnLassoMut_All <- predictMnLasso(MnLassoMut_All,X)
toc()

# Using Lasso
cat("Lasso Mut\n")
tic()
LassoMut_All <- trainLasso(X, C)
TratamientoLassoMut_All <- predictLasso(LassoMut_All,X)
toc()


# Use same groups as previously
TratamientoBOSOMut_All <- rep(NA, length(Groups))
tic()
BOSOMut_All <- trainBOSO(X, C,maxVarsBlock = 10, standardize=F, IC = "eBIC")
TratamientoBOSOMut_All <- predictBOSO(BOSOMut_All, X)
toc()

```

### Validation in GDSC

The models were trained using the BeatAML dataset and were used to predict the optimal drug for AML cell lines from GDSC using its mutation files.
Each of the cell lines was recommended a drug, we compared the all-samples IC50 for all the models and against the Oracle (the drug with the minimum IC50 for each cell line).
The function `predictMOM_ALL`.
Applies the decision tree generated by MOM in AML.

#### Load GDSC and CCLE data

```{r}
GDSC1_ESS<-read.csv("./data/input/LAML_IC_Wed Nov 24 10_31_42 2021.csv")
GDSC1_MUT<-read.csv("./data/input/LAML_Genetic_feature_variant_Wed Nov 24 GDSC1.csv")



load("./data/input/DEMETER2/Processed/new_CCLE_DEMETER2_20Q3_SYNLETHDB_data_v1.Rdata")
CCLE_mutation_data <- demeter2_mut

CCLE_mutation_data[c("MV411_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","NOMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","OCIAML2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"),"TP53"]   <- 1
CCLE_mutation_data["HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","KDM6A"]  <- 1
CCLE_mutation_data[c("MOLM13_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","MV411_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"),"FLT3"]   <- 1
CCLE_mutation_data["NOMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","KRAS"]   <- 1

rm(demeter2_mut, demeter2_scores, demeter2_fus)




```

#### Biomarker Pre-processing

We build the mutations matrix and complete that information with CCLE mutational data.

```{r}
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
                  ncol=ncol(mut2))

colnames(BM_matrix)<-colnames(mut2)
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


BM_matrix<-BM_matrix[,match( colnames(BM_matrix),colnames(mut2))]
# colnames(BM_matrix)<-sapply(colnames(BM_matrix), FUN=function(X){return(unlist(strsplit(X, split="_"))[1])})

Mutations_AML_GDSC<-BM_matrix


```

#### Drug information pre-processing

```{r}

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

drugsGDSC<-sapply(colnames(C), FUN=function(X){return(unlist(strsplit(X, split=" "))[1])})

GDSC_AML<-GDSC_AML[,na.omit(match(drugsGDSC, colnames(GDSC_AML)))]

colnamesGDSC<-colnames(C)[match(colnames(GDSC_AML),drugsGDSC)]


colnames(GDSC_AML)<-colnamesGDSC


GDSC_AML<-GDSC_AML[which(rownames(GDSC_AML) %in% rownames(Mutations_AML_GDSC)),]
GDSC_AML<-GDSC_AML[match(rownames(Mutations_AML_GDSC), rownames(GDSC_AML)),]


identical(rownames(GDSC_AML), rownames(Mutations_AML_GDSC))

```

#### Make Predictions with the models

We used the models trained in BeatAML to predict over GDSC.

```{r}
identical(colnames(Mutations_AML_GDSC), colnames(mut2)) #TRUE


Guideline_TreeMutSqrt<-predictTreeMut(treeMut1_All,C, Mutations_AML_GDSC)
Guideline_TreeMutSqrt<-colnames(C)[Guideline_TreeMutSqrt]


Guideline_BOSO<-predictBOSO(BOSOMut_All,Mutations_AML_GDSC)
Guideline_BOSO<-colnames(C)[Guideline_BOSO]

Guideline_TreeMut<-predictTreeMut(treeMut_All,C, Mutations_AML_GDSC)
Guideline_TreeMut<-colnames(C)[Guideline_TreeMut]

Guideline_Multinomial<-predictMnLasso(MnLassoMut_All, Mutations_AML_GDSC)
Guideline_Multinomial<-colnames(C)[Guideline_Multinomial]

Guideline_Lasso<-predictLasso(LassoMut_All, Mutations_AML_GDSC)
Guideline_Lasso<-colnames(C)[Guideline_Lasso]

Guideline_MOM<-predictMOM_ALL(Mutations_AML_GDSC)
Guideline_MOM<-colnames(C)[Guideline_MOM$Treatment]

Guideline_KRL<-read_csv("data/output/KRL_GDSC_pred.csv")
Guideline_KRL$DrugName[which(!(Guideline_KRL$DrugName %in% colnames(C)))]<-NA
Guideline_KRL<-Guideline_KRL$DrugName

```

#### Obtain Prediction IC50 values

We will now obtain the GDSC IC50 values corresponding to the predicted treatments

```{r}
# Validate the Guidelines..................................................

identical(rownames(GDSC_AML), rownames(Mutations_AML_GDSC)) # TRUE


GDSC_IC50<-data.frame(CL=rownames(Mutations_AML_GDSC),
                      TreeMut=NA,
                      TreeMutSqrt=NA,
                      MnLasso=NA,
                      BOSO=NA,
                      Lasso=NA,
                      MOM=NA,
                      KRL=NA,
                      ORACLE=NA)

for(i in 1:nrow(GDSC_IC50)){
  CL<-GDSC_IC50$CL[i]
  
  # TreeMut
  dtreeMut<-Guideline_TreeMut[i]
  if(dtreeMut %in% colnames(GDSC_AML)){
    GDSC_IC50$TreeMut[i]<-GDSC_AML[CL,dtreeMut]
  }
  
  
  # TreeMutSqrt
  dtreeMutSQrt<-Guideline_TreeMutSqrt[i]
  GDSC_IC50$TreeMutSqrt[i]<-GDSC_AML[CL,dtreeMutSQrt]
  
  # Multinomial
  dMultinomial<-Guideline_Multinomial[i]
  if(dMultinomial %in% colnames(GDSC_AML)){
    GDSC_IC50$MnLasso[i]<-GDSC_AML[CL,dMultinomial]
  }
  
  
  
  # BOSO
  dBOSO<-Guideline_BOSO[i]
  if(dBOSO %in% colnames(GDSC_AML)){
    GDSC_IC50$BOSO[i]<-GDSC_AML[CL,dBOSO]
  }
  
  
  
  # Lasso
  dLasso<-Guideline_Lasso[i]
  if(dLasso %in% colnames(GDSC_AML)){
    GDSC_IC50$Lasso[i]<-GDSC_AML[CL,dLasso]
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

```

#### Plotting GDSC Validation

```{r}
Validation_plot<-data.frame(Dataset=rep("GDSC", nrow(GDSC_IC50)*8), 
                                                   CellLine=rep(rownames(GDSC_IC50),8),
                                                   Method=rep(colnames(GDSC_IC50[,-1]), each=nrow(GDSC_IC50)), 
                                                   IC50=unlist(GDSC_IC50[,-1]))


rownames(Validation_plot)<-1:nrow(Validation_plot)

mycomparisons<-combn(unique(Validation_plot$Method),2)
mycomparisons<-as.list(data.frame(mycomparisons))
mycomparisons<-mycomparisons[-grep("ORACLE",mycomparisons)]


Plotting<-Validation_plot[Validation_plot$Dataset=="GDSC",]


library(ggsci)

Plotting$Method<-factor(Plotting$Method, levels=c("ORACLE", "MOM", "BOSO", "TreeMutSqrt", "Lasso", "TreeMut", "MnLasso", "KRL"))

g1_gdsc<-ggplot(Plotting, aes(x=Method, y=IC50, fill=Method))+geom_violin()+
  theme_bw()+scale_fill_npg()+ylab("GDSC IC50")+ggtitle("GDSC")

g2_gdsc<-ggplot(Plotting, aes(x=Method, y=IC50, fill=Method))+geom_boxplot()+
  theme_bw()+scale_fill_npg()+ylab("GDSC IC50")+ggtitle("GDSC")

g1_gdsc
g2_gdsc

```
</details>



## Intragroup Validation

We compared whether the IC50\* of a drug in patients in whom it was recommended was lower than the IC50\* in patients in whom it was not recommended.
Using this information we compared the sensitivity to a drug for a specific group against the sensitivity to that drug for the rest of the samples by using a 2-tailed Wilcoxon test.
This analysis was performed both for the BeatAML dataset (training dataset) and the GDSC AML cell lines cohort (predicted dataset).

<details><summary>Click to expand</summary>
### Obtaining MOM and KRL predicitions in all samples

```{r}
TratamientoMOM_All<-predictMOM_ALL(mut2)
identical(TratamientoMOM_All$CLs, rownames(C)) # TRUE

TratamientoKRL_All<- read_csv("data/output/KRL_results_all_samples_prediction.csv", show_col_types = FALSE)
TratamientoKRL_All<-TratamientoKRL_All[which(TratamientoKRL_All$Patients %in% rownames(mut2)),]
TratamientoKRL_All$DrugName[which(!(TratamientoKRL_All$DrugName %in% colnames(C)))]<-NA

```

</details>


### Intragroup Validation in BeatAML

<details><summary>Click to expand</summary>
```{r}
# 1. BeatAML--------------------------------------------------------------------

Validation_Beat<-data.frame(patients=rep(rownames(C),7),
                            drug=c(colnames(C)[TratamientoTreeMut_ALL],
                                   colnames(C)[TratamientoTreeMut1_All],
                                   colnames(C)[TratamientoMnLassoMut_All],
                                   colnames(C)[TratamientoLassoMut_All],
                                   colnames(C)[TratamientoBOSOMut_All],
                                   colnames(C)[TratamientoMOM_All$Treatment],
                                   TratamientoKRL_All$DrugName), 
                            IC50=c(C[cbind(1:nrow(C), TratamientoTreeMut_ALL)],
                                   C[cbind(1:nrow(C), TratamientoTreeMut1_All)],
                                   C[cbind(1:nrow(C), TratamientoMnLassoMut_All)],
                                   C[cbind(1:nrow(C), TratamientoLassoMut_All)],
                                   C[cbind(1:nrow(C), TratamientoBOSOMut_All)],
                                   C[cbind(1:nrow(C), TratamientoMOM_All$Treatment)],
                                   C[cbind(TratamientoKRL_All$Patients,TratamientoKRL_All$DrugName )]),
                            Method=rep(unique(Tratamiento_plot$Method[!(Tratamiento_plot$Method=="ORACLE")]), each=nrow(C)))


```

#### MOM

```{r}

Plot_aux<-Validation_Beat[Validation_Beat$Method=="MOM",]

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
  plotDrug$IC50[!ind]<-C[!ind, d]
  plotDrug$Recommended[!ind]<-"NO"
  
  MOM_plots<-rbind(MOM_plots,plotDrug)
  
}

gg1_validation_beat<-ggplot(MOM_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+
  scale_fill_npg()+stat_compare_means(vjust = 1)+ylab("IC50*")+ggtitle("MOM")

gg1_validation_beat
```

#### TreeMut

```{r}
# TREE MUT

Plot_aux<-Validation_Beat[Validation_Beat$Method=="TreeMut",]

Tree_plots<-NULL

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
  plotDrug$IC50[!ind]<-C[!ind, d]
  plotDrug$Recommended[!ind]<-"NO"
  
  Tree_plots<-rbind(Tree_plots,plotDrug)
  
}

gg3_validation_beat<-ggplot(Tree_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_npg()+
  stat_compare_means()+ylab("IC50*")+ggtitle("TreeMut")
gg3_validation_beat

```

#### TreeMutSqrt

```{r}
# TREE MUT SQRT

Plot_aux<-Validation_Beat[Validation_Beat$Method=="TreeMutSqrt",]

TreeSQR_plots<-NULL

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
  plotDrug$IC50[!ind]<-C[!ind, d]
  plotDrug$Recommended[!ind]<-"NO"
  
  TreeSQR_plots<-rbind(TreeSQR_plots,plotDrug)
  
}

gg4_validation_beat<-ggplot(TreeSQR_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_npg()+
  stat_compare_means()+ylab("IC50*")+ggtitle("TreeMutSqrt")

gg4_validation_beat
```

#### BOSO

```{r}
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
  plotDrug$IC50[!ind]<-C[!ind, d]
  plotDrug$Recommended[!ind]<-"NO"
  
  BOSO_plots<-rbind(BOSO_plots,plotDrug)
  
}

gg2_validation_beat<-ggplot(BOSO_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_npg()+
  stat_compare_means(vjust = 1)+ylab("IC50*")+ggtitle("BOSO")

gg2_validation_beat

```

#### Multinomial

```{r}
# MULTINOMIAL

Plot_aux<-Validation_Beat[Validation_Beat$Method=="Multinomial",]

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
  plotDrug$IC50[!ind]<-C[!ind, d]
  plotDrug$Recommended[!ind]<-"NO"
  
  MnLasso_plots<-rbind(MnLasso_plots,plotDrug)
  
}

gg5_validation_beat<-ggplot(MnLasso_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_npg()+
  stat_compare_means(vjust=1)+ylab("IC50*")+ggtitle("Multinomial")
gg5_validation_beat
```

#### KRL

```{r}
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
  plotDrug$IC50[!ind]<-C[!ind, d]
  plotDrug$Recommended[!ind]<-"NO"
  
 KRL_plots<-rbind(KRL_plots,plotDrug)
  
}

gg7_validation_beat<-ggplot(KRL_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_npg()+
  stat_compare_means(vjust=1)+ylab("IC50*")+ggtitle("KRL")

gg7_validation_beat
```
</details>


### Intragroup Validation in GDSC


<details><summary>Click to expand</summary>
```{r}
identical(colnames(Mutations_AML_GDSC), colnames(mut2)) #TRUE


Guideline_TreeMutSqrt<-predictTreeMut(treeMut1_All,C, Mutations_AML_GDSC)
Guideline_TreeMutSqrt<-colnames(C)[Guideline_TreeMutSqrt]


Guideline_BOSO<-predictBOSO(BOSOMut_All,Mutations_AML_GDSC)
Guideline_BOSO<-colnames(C)[Guideline_BOSO]

Guideline_TreeMut<-predictTreeMut(treeMut_All,C, Mutations_AML_GDSC)
Guideline_TreeMut<-colnames(C)[Guideline_TreeMut]

Guideline_Multinomial<-predictMnLasso(MnLassoMut_All, Mutations_AML_GDSC)
Guideline_Multinomial<-colnames(C)[Guideline_Multinomial]

Guideline_Lasso<-predictLasso(LassoMut_All, Mutations_AML_GDSC)
Guideline_Lasso<-colnames(C)[Guideline_Lasso]

Guideline_MOM<-predictMOM_ALL(Mutations_AML_GDSC)
Guideline_MOM<-colnames(C)[Guideline_MOM$Treatment]

Guideline_KRL<-read_csv("data/output/KRL_GDSC_pred.csv")
Guideline_KRL$DrugName[which(!(Guideline_KRL$DrugName %in% colnames(C)))]<-NA
Guideline_KRL<-Guideline_KRL$DrugName

Validation_GDSC<-data.frame(patients=rep(rownames(GDSC_AML),7),
                            drug=c(Guideline_TreeMut,
                                   Guideline_TreeMutSqrt,
                                   Guideline_Multinomial,
                                   Guideline_Lasso,
                                   Guideline_BOSO,
                                   Guideline_MOM, 
                                   Guideline_KRL), 
                            IC50=NA,
                            Method=rep(unique(Tratamiento_plot$Method[!(Tratamiento_plot$Method=="ORACLE")]), each=nrow(GDSC_AML)))




```

#### MOM

```{r}
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

```

#### TreeMut

```{r}
# TREE MUT

Plot_aux<-Validation_GDSC[Validation_GDSC$Method=="TreeMut",]

Tree_plots<-NULL

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
  
  
  Tree_plots<-rbind(Tree_plots,plotDrug)
  
}

gg3_validation_gdsc<-ggplot(Tree_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_manual(values=c("#999999", "#E69F00"))+
  stat_compare_means()+ylab("IC50*")+ggtitle("TreeMut")

gg3_validation_gdsc
```

#### TreeMutSqrt

```{r}
# TREE MUT SQRT

Plot_aux<-Validation_GDSC[Validation_GDSC$Method=="TreeMutSqrt",]

TreeSQR_plots<-NULL

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

  
  TreeSQR_plots<-rbind(TreeSQR_plots,plotDrug)
  
}

gg4_validation_gdsc<-ggplot(TreeSQR_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_manual(values=c("#999999", "#E69F00"))+
  stat_compare_means()+ylab("IC50*")+ggtitle("TreeMutSqrt")

gg4_validation_gdsc
```
</details>



## Multiomics in BeatAML

Some of the methods only accept as input binary variables.
Although, genomic variants can be transformed into binary variables, gene expression, methylation, or openness of the chromatin are intrinsically continuous variables.
We have included a table showing whether the algorithm accepts only binary inputs (genomic variants only) or whether it also accepts continuous data (gene expression, methylation, etc.).
For the methods that accept continuous variables, we assessed the performance of the predictions (5-fold cross-validation) in the BeatAML dataset using both data sources.
We state the statistical significance using a 2-tail Wilcoxon's test comparing the IC50\* using as input either genetic variants or gene expression.

<details><summary>Click to expand</summary>
### Load expression data

```{r}
# # Read expression data
Expression <- read_excel("./data/input/41586_2018_623_MOESM3_ESM.xlsx",
                         sheet = "Table S8-Gene Counts RPKM")

# Keep genes with the largest variance
ngenes <- 1000
Expression <- as.data.frame(Expression)
rownames(Expression) <- Expression[,2]
Expression <- Expression[,-c(1,2)]
OldExpression <- as.matrix(Expression)
sdevs <- rowSds(OldExpression)
# sdevs <- rowMedians(OldExpression)

Keep <- which(order(sdevs, decreasing = T) <= ngenes)
Expression <- Expression[Keep,]

commonPatients <- intersect(rownames(C), colnames(Expression))
Expression <- Expression[,commonPatients]
Expression <- t(Expression)

C <- C[commonPatients, ]
X<-mut2
X<-X[commonPatients, ]
Y <- exp(-.5*C) # Convert weights into pseudoprobabilities
Y <- Y /rowSums(Y) # Probabilities sum up one
```

### Train Models with Expression Data

```{r}
library(robustbase)
set.seed(2022)
Folds <- 5
Groups <- sample(1:Folds, nrow(C), replace = T)
TratamientoMnLasso <- TratamientoLasso <- TratamientoTree <- TratamientoTree1 <- rep(NA, length(Groups))
TratamientoMnLassoMut <- TratamientoLassoMut <- TratamientoTreeMut <- TratamientoTreeMut1 <- rep(NA, length(Groups))

for (group in 1:Folds) {
  cat("Fold: ", group, "\n")
  Dejar <- Groups != group
  Quitar <- !Dejar

  # Trees
  cat("Trees\n")
  tree <- trainTree(Expression[Dejar,], C[Dejar,], minbucket = 20)
  TratamientoTree[Quitar] <- predictTree(tree,C[Quitar,], Expression[Quitar,])

  # Trees sqrt
  cat("Trees sqrt\n")
  tree1 <- trainTree(Expression[Dejar,], (C[Dejar,])^.5, minbucket = 20)
  TratamientoTree1[Quitar] <- predictTree(tree1,C[Quitar,], Expression[Quitar,])

  # Using MnLasso
  cat("Multinomial Lasso\n")
  MnLasso <- trainMnLasso(Expression[Dejar,], Y[Dejar,])
  TratamientoMnLasso[Quitar] <- predictMnLasso(MnLasso,Expression[Quitar,])

  # Using Lasso
  cat("Lasso\n")
  Lasso <- trainLasso(Expression[Dejar,], C[Dejar,])
  TratamientoLasso[Quitar] <- predictMnLasso(Lasso,Expression[Quitar,])

  # TreesMut
  cat("TreesMut\n")
  treeMut <- trainTreeMut(X[Dejar,], C[Dejar,], minbucket = 20)
  TratamientoTreeMut[Quitar] <- predictTreeMut(treeMut,C[Quitar,], X[Quitar,])

  # TreesMut sqrt
  cat("TreesMut sqrt\n")
  treeMut1 <- trainTreeMut(X[Dejar,], (C[Dejar,])^.5, minbucket = 20)
  TratamientoTreeMut1[Quitar] <- predictTreeMut(treeMut1,C[Quitar,], X[Quitar,])

  # Using MnLasso
  cat("Multinomial Lasso Mut\n")
  MnLassoMut <- trainMnLasso(X[Dejar,], Y[Dejar,])
  TratamientoMnLassoMut[Quitar] <- predictMnLasso(MnLassoMut,X[Quitar,])

  # Using Lasso
  cat("Lasso Mut\n")
  LassoMut <- trainLasso(X[Dejar,], C[Dejar,])
  TratamientoLassoMut[Quitar] <- predictLasso(LassoMut,X[Quitar,])

}


# # BOSO should be included in the loop. However, the difference between train and test is not that large.
# # Let's see how it works out of the box

#BOSO <- trainBOSO(Expression, C)
#TratamientoBOSO <- predictBOSO(BOSO, Expression)

set.seed(2021)
Folds <- 5
Groups <- sample(1:Folds, 257, replace = T)

# Use same groups as previously
TratamientoBOSO <- TratamientoBOSOMutasso <- rep(NA, length(Groups))

for (group in 1:Folds) {
  BOSO <- trainBOSO(Expression[Dejar,], C[Dejar,],maxVarsBlock = 10, standardize=F, IC = "eBIC")
  TratamientoBOSO[Quitar] <- predictBOSO(BOSO, Expression[Quitar,])
  BOSOMut <- trainBOSO(X[Dejar,], C[Dejar,],maxVarsBlock = 10, standardize=F, IC = "eBIC")
  TratamientoBOSOMut[Quitar] <- predictBOSO(BOSOMut, X[Quitar,])
}

TratamientoOracle <- apply(C, 1, which.min)

Tratamientos <- cbind(C[cbind(1:nrow(C), TratamientoOracle)],
                      C[cbind(1:nrow(C), TratamientoTree)],
                      C[cbind(1:nrow(C), TratamientoTree1)],
                      C[cbind(1:nrow(C), TratamientoMnLasso)],
                      C[cbind(1:nrow(C), TratamientoLasso)],
                      C[cbind(1:nrow(C), TratamientoBOSO)],
                      C[cbind(1:nrow(C), TratamientoTreeMut)],
                      C[cbind(1:nrow(C), TratamientoTreeMut1)],
                      C[cbind(1:nrow(C), TratamientoMnLassoMut)],
                      C[cbind(1:nrow(C), TratamientoLassoMut)],
                      C[cbind(1:nrow(C), TratamientoBOSOMut)]
)
#
colnames(Tratamientos) <- c("Oracle","Tree","TreeSqrt","MnLasso","Lasso","BOSO",
                            "TreeMut","TreeSrqtMut","MnLassoMut","LassoMut","BOSOMut")


```

### Plotting Results

```{r}
boxplot(Tratamientos)
plot_Ex_Mut<-data.frame(Method=rep(colnames(Tratamientos), each=nrow(Tratamientos)))
plot_Ex_Mut$IC50<-unlist(as.data.frame(Tratamientos))
plot_Ex_Mut$Type<-"GE"
plot_Ex_Mut$Type[grep("Mut", plot_Ex_Mut$Method)]<-"Mut"


plot_Ex_Mut<-plot_Ex_Mut[order(plot_Ex_Mut$Method),]
plot_Ex_Mut$General_Method<-"unknown"
plot_Ex_Mut$General_Method[grep("Oracle", plot_Ex_Mut$Method)]<-"Oracle"
plot_Ex_Mut$General_Method[grep("BOSO", plot_Ex_Mut$Method)]<-"BOSO"
plot_Ex_Mut$General_Method[grep("Tree", plot_Ex_Mut$Method)]<-"Tree"
plot_Ex_Mut$General_Method[grep("TreeSqrt", plot_Ex_Mut$Method)]<-"TreeSqrt"
plot_Ex_Mut$General_Method[grep("TreeSrqtMut", plot_Ex_Mut$Method)]<-"TreeSqrt"
plot_Ex_Mut$General_Method[grep("Lasso", plot_Ex_Mut$Method)]<-"Lasso"
plot_Ex_Mut$General_Method[grep("MnLasso", plot_Ex_Mut$Method)]<-"MnLasso"

# plot_Ex_Mut<-rbind(plot_Ex_Mut, data.frame(Method=NA, IC50=5, Type="Mut", General_Method="MOM"))

plot_Ex_Mut$Method<-factor(plot_Ex_Mut$Method, levels=c("Oracle","BOSO","BOSOMut","TreeSqrt","TreeSrqtMut","Lasso","LassoMut",
                                                        "Tree","TreeMut",
                                                        "MnLasso","MnLassoMut"))
plot_Ex_Mut$General_Method<-factor(plot_Ex_Mut$General_Method, levels=c("Oracle", "BOSO", "TreeSqrt", "Lasso",
                                                                        "Tree","MnLasso"))


mycomparisons<-list(c("Tree", "TreeMut"), c("TreeSqrt", "TreeSrqtMut"), c("Lasso", "LassoMut"), c("BOSO", "BOSOMut"),
                    c("MnLasso", "MnLassoMut"))
mypal = pal_npg("nrc")(9)
mypal<-mypal[c(1,3:7)]
gg_ex_mut<-ggplot(plot_Ex_Mut, aes(x=Method, y=IC50, fill=General_Method, alpha=Type))+geom_boxplot()+
  scale_fill_manual(values=mypal)+theme_bw()+ggtitle("5-fold CV in BeatAML")+ylab("IC50*")+labs(fill="General Method")+
  stat_compare_means(aes(Method),vjust=0, hjust=-1, comparisons=mycomparisons,label.y = 5)+scale_alpha_manual(values=c( 0.5, 1))

gg_ex_mut

```
</details>


