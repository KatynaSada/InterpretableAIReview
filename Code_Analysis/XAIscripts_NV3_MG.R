####################################################################
#################################Script#############################
####################################################################
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
library(maftools)
library(qvalue)
library(tibble)
library(Rtsne)
library(IHW)
library(pheatmap)
library(DT)
library("dplyr")
library(tictoc)

source("./MOM/2022-05-24_MOM_source_Functions.R")


# Load mutational data----------------------------------------------------------

gene_variants<-read_excel("./data/input/41586_2018_623_MOESM3_ESM.xlsx",
                          sheet="Table S7-Variants for Analysis")
clinical <- read_excel("./data/input/41586_2018_623_MOESM3_ESM.xlsx",
                       sheet = "Tabe S5-Clinical Summary")

# Load drug sensitivity data----------------------------------------------------

Drug_response<-read_excel("./data/input/41586_2018_623_MOESM3_ESM.xlsx",
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

# Combining Matrices and Imputation---------------------------------------------
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


# Tidying up the variables

rm(list = c("Drug_2", "Drug_3", "Drug_4", "Drug_response", "Mutations", "gene_variants"))
rm(list = c("Out", "translocations_clinical", "clinical", "add", "Drug"))


# Read expression data
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

commonPatients <- intersect(rownames(drug), colnames(Expression))


C <- drug - min(drug)

# Prueba con tree mutations

treeMut <- trainTreeMut(mut2, C, minbucket = 10)
treeMutBonito <- niceTree(treeMut)
plot(treeMutBonito)
TratamientoTreeMut <- predictTreeMut(treeMut, C, mut2)
sort(table(TratamientoTreeMut))
TratamientoOracle <- apply(C, 1, which.min)

Tratamientos <- cbind(C[cbind(1:nrow(C), TratamientoOracle)],
                      C[cbind(1:nrow(C), TratamientoTreeMut)])
boxplot(Tratamientos)

# Prueba "en serio" con CV
X <- mut2

Y <- exp(-.5*C) # Convert weights into pseudoprobabilities
Y <- Y /rowSums(Y) # Probabilities sum up one


#####################################################
##### Comparison of methods with Cross-Validation
#####################################################
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


TratamientoOracle <- apply(C, 1, which.min)
Tratamientos <- cbind(C[cbind(1:nrow(C), TratamientoOracle)],
                      C[cbind(1:nrow(C), TratamientoTreeMut)],
                      C[cbind(1:nrow(C), TratamientoTreeMut1)],
                      C[cbind(1:nrow(C), TratamientoLassoMut)],
                      C[cbind(1:nrow(C), TratamientoMnLassoMut)],
                      C[cbind(1:nrow(C), TratamientoBOSOMut)], MOM$IC50,KRL)

colnames(Tratamientos) <- c("Oracle","TreeMut","TreeMut sqrt","MnLassoMut","LassoMut","BOSOMut",
                            "MOM", "KRL")

boxplot(Tratamientos)

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

# PLOTTING----------------------------------------------------------------------
load("./data/output/mom_COMPARISON.RData")
library(ggsci)

Tratamiento_plot<-data.frame(Method="ORACLE", IC50=C[cbind(1:nrow(C),TratamientoOracle)])
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="TreeMut", IC50=C[cbind(1:nrow(C), TratamientoTreeMut)]))
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="TreeMutSqrt", IC50=C[cbind(1:nrow(C), TratamientoTreeMut1)]))
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="Multinomial", IC50=C[cbind(1:nrow(C), TratamientoMnLassoMut)]))
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="LassoMut", IC50=C[cbind(1:nrow(C), TratamientoLassoMut)]))
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="BOSOMut", IC50=C[cbind(1:nrow(C), TratamientoBOSOMut)]))
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="MOM", IC50=MOM$IC50))
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="KRL", IC50=KRL_results_cv$IC50))
Tratamiento_plot<-rbind(Tratamiento_plot,data.frame(Method="UpperRef", IC50=C[cbind(1:nrow(C), which.min(colMeans(C)))]))

Tratamiento_plot$Method<-factor(Tratamiento_plot$Method, levels=c("ORACLE", "MOM", "BOSOMut", "TreeMutSqrt", "LassoMut", "TreeMut", "Multinomial","KRL", "UpperRef"))

ggplot(Tratamiento_plot, aes(x=Method, y=IC50, fill=Method))+geom_violin()+
  theme_bw()+scale_fill_npg()+ylab("IC50*")+ggtitle("BeatAML 5-fold Cross Validation")

ggplot(Tratamiento_plot, aes(x=Method, y=IC50, fill=Method))+geom_boxplot()+
  theme_bw()+scale_fill_npg()+ylab("IC50*")+theme(text = element_text(size = 15)) 


wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="UpperRef"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="MOM"], alternative="greater", paired = T)
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="UpperRef"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="BOSOMut"], alternative="greater", paired = T)
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="UpperRef"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="TreeMutSqrt"], alternative="greater", paired = T)
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="UpperRef"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="LassoMut"], alternative="greater", paired = T)
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="UpperRef"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="TreeMut"], alternative="greater", paired = T)
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="UpperRef"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="Multinomial"], alternative="greater", paired = T)
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="UpperRef"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="KRL"], alternative="greater", paired = T)



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

wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="Multinomial"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="MOM"], alternative="less")
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="Multinomial"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="BOSOMut"], alternative="less")
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="Multinomial"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="TreeMutSqrt"], alternative="less")
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="Multinomial"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="LassoMut"], alternative="less")
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="Multinomial"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="TreeMut"], alternative="less")
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="Multinomial"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="Multinomial"])
wilcox.test(Tratamiento_plot$IC50[Tratamiento_plot$Method=="Multinomial"], 
            Tratamiento_plot$IC50[Tratamiento_plot$Method=="KRL"], alternative="less")



# Method Precision lollipop plot------------------------------------------------
Tratamiento_plot$Score[Tratamiento_plot$Method=="ORACLE"]<-median(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"]-Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"])
Tratamiento_plot$Score[Tratamiento_plot$Method=="TreeMut"]<-median(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"]-Tratamiento_plot$IC50[Tratamiento_plot$Method=="TreeMut"])
Tratamiento_plot$Score[Tratamiento_plot$Method=="TreeMutSqrt"]<-median(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"]-Tratamiento_plot$IC50[Tratamiento_plot$Method=="TreeMutSqrt"])
Tratamiento_plot$Score[Tratamiento_plot$Method=="Multinomial"]<-median(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"]-Tratamiento_plot$IC50[Tratamiento_plot$Method=="Multinomial"])
Tratamiento_plot$Score[Tratamiento_plot$Method=="LassoMut"]<-median(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"]-Tratamiento_plot$IC50[Tratamiento_plot$Method=="LassoMut"])
Tratamiento_plot$Score[Tratamiento_plot$Method=="BOSOMut"]<-median(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"]-Tratamiento_plot$IC50[Tratamiento_plot$Method=="BOSOMut"])
Tratamiento_plot$Score[Tratamiento_plot$Method=="MOM"]<-median(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"]-Tratamiento_plot$IC50[Tratamiento_plot$Method=="MOM"])
Tratamiento_plot$Score[Tratamiento_plot$Method=="KRL"]<-median(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"]-Tratamiento_plot$IC50[Tratamiento_plot$Method=="KRL"])
Tratamiento_plot$Score[Tratamiento_plot$Method=="UpperRef"]<-median(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"]-Tratamiento_plot$IC50[Tratamiento_plot$Method=="UpperRef"])


Tratamiento_plot$Score<-abs(Tratamiento_plot$Score)


Lollipop_plot<-unique(data.frame(Tratamiento_plot$Method, Tratamiento_plot$Score))

Lollipop_plot<-Lollipop_plot[-which(Lollipop_plot$Tratamiento_plot.Method=="ORACLE"),]
mypal = pal_npg("nrc")(10)
mypal<-mypal[2:10]

gg_lolliplot<-ggplot(Lollipop_plot, aes(x=Tratamiento_plot.Score, y=Tratamiento_plot.Method))+
  geom_segment(aes(x = 0, xend = Tratamiento_plot.Score, y = Tratamiento_plot.Method, yend = Tratamiento_plot.Method),
               color = "gray", lwd = 1.5)+
  geom_point( shape=21, size=8, colour="black",stroke=2,aes(fill=Tratamiento_plot.Method)) +
  labs(fill="XAI Methods")+scale_fill_manual(values=mypal)+
  theme_bw(base_size = 18) +ylab("XAI Methods")+xlab(expression(paste("Oracle Loss=median(", Delta, "IC50*)")))

gg_lolliplot


#Method Loss Horizontal Boxplot
Tratamiento_plot$Loss[Tratamiento_plot$Method=="ORACLE"]<-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"])-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"])
Tratamiento_plot$Loss[Tratamiento_plot$Method=="TreeMut"]<-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"])-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="TreeMut"])
Tratamiento_plot$Loss[Tratamiento_plot$Method=="TreeMutSqrt"]<-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"])-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="TreeMutSqrt"])
Tratamiento_plot$Loss[Tratamiento_plot$Method=="Multinomial"]<-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"])-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="Multinomial"])
Tratamiento_plot$Loss[Tratamiento_plot$Method=="LassoMut"]<-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"])-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="LassoMut"])
Tratamiento_plot$Loss[Tratamiento_plot$Method=="BOSOMut"]<-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"])-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="BOSOMut"])
Tratamiento_plot$Loss[Tratamiento_plot$Method=="MOM"]<-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"])-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="MOM"])
Tratamiento_plot$Loss[Tratamiento_plot$Method=="KRL"]<-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"])-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="KRL"])
Tratamiento_plot$Loss[Tratamiento_plot$Method=="UpperRef"]<-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="ORACLE"])-abs(Tratamiento_plot$IC50[Tratamiento_plot$Method=="UpperRef"])

MethodLoss<-Tratamiento_plot[-which(Tratamiento_plot$Method=="ORACLE"),]
MethodLoss$Loss<-abs(MethodLoss$Loss)

mypal2 = pal_npg("nrc")(10)
mypal2<-mypal2[2:10]

ggplot(MethodLoss, aes(x=Loss, y=Method, fill=Method))+geom_boxplot()+theme_bw(base_size = 18)+
  scale_fill_manual(values=mypal2)
  


################################################################################
# VALIDATION IN GDSC, CERES AND DEMETER2
################################################################################

load("./data/input/Drug&Target.RData")
colnames(DT)<-c("Drug", "Target")
DT<-as.data.frame(DT)

# 1. DEMETER2--------------------------------------------------------------------

# Load DEMETER2 files

load("./data/input/DEMETER2/new_CCLE_sample_info_match.RData")
load("./data/input/DEMETER2/Processed/new_CCLE_DEMETER2_20Q3_SYNLETHDB_data_v1.Rdata")
CCLE_mutation_data <- demeter2_mut
CCLE_demeter<-demeter2_scores
CCLE_Fusions <-demeter2_fus


CCLE_mutation_data[c("MV411_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","NOMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","OCIAML2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"),"TP53"]   <- 1
CCLE_mutation_data["HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","KDM6A"]  <- 1
CCLE_mutation_data[c("MOLM13_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","MV411_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"),"FLT3"]   <- 1
CCLE_mutation_data["NOMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","KRAS"]   <- 1


rm(demeter2_mut, demeter2_scores, demeter2_fus)


# Pre-Process Drug Information to get biomarkers and essentiality scores in AML CLs


AML_CLs<-new_CCLE_sample_info_match$CCLE.name[new_CCLE_sample_info_match$Subtype=="acute_myeloid_leukaemia"]
CCLE_mutation_data_AML<-CCLE_mutation_data[(rownames(CCLE_mutation_data) %in% AML_CLs),]
CCLE_demeter_AML<-CCLE_demeter[(rownames(CCLE_demeter) %in% AML_CLs),]

CCLE_Fusions_AML<-CCLE_Fusions[(rownames(CCLE_Fusions) %in% AML_CLs),]

colnames(CCLE_Fusions_AML)<-sapply(colnames(CCLE_Fusions_AML), FUN=function(X){
  return(paste(unlist(strsplit(X, split="--")), collapse = "-"))
})

Combined_CCLE_Mut<-cbind(CCLE_mutation_data_AML,CCLE_Fusions_AML)

# Build Training Matrix

Mut_AML_DEMETER2<-Combined_CCLE_Mut[,match(colnames(mut2), colnames(Combined_CCLE_Mut))]
colnames(Mut_AML_DEMETER2)<-colnames(mut2)

CCLE_demeter_AML<-t(t(CCLE_demeter_AML)-colMeans(CCLE_demeter_AML))

Mut_AML_DEMETER2[is.na(Mut_AML_DEMETER2)]<-0

# write.csv(Mut_AML_DEMETER2, file="./KRL/mut_Matrix_DEMETER2.csv")
# Predict with the Models.........................................................
source("./MOM/predictMOM_ALL.R")

Guideline_TreeMutSqrt<-predictTreeMut(treeMut1_All,C, Mut_AML_DEMETER2)
Guideline_TreeMutSqrt<-colnames(C)[Guideline_TreeMutSqrt]


Guideline_BOSO<-predictBOSO(BOSOMut_All,Mut_AML_DEMETER2)
Guideline_BOSO<-colnames(C)[Guideline_BOSO]

Guideline_TreeMut<-predictTreeMut(treeMut_All,C, Mut_AML_DEMETER2)
Guideline_TreeMut<-colnames(C)[Guideline_TreeMut]

Guideline_Multinomial<-predictMnLasso(MnLassoMut_All, Mut_AML_DEMETER2)
Guideline_Multinomial<-colnames(C)[Guideline_Multinomial]

Guideline_Lasso<-predictLasso(LassoMut_All, Mut_AML_DEMETER2)
Guideline_Lasso<-colnames(C)[Guideline_Lasso]

Guideline_MOM<-predictMOM_ALL(Mut_AML_DEMETER2)
Guideline_MOM<-colnames(C)[Guideline_MOM$Treatment]

Guideline_KRL<-read_csv("data/output/KRL_DEMETER_pred.csv")
Guideline_KRL$DrugName[which(!(Guideline_KRL$DrugName %in% colnames(C)))]<-NA
Guideline_KRL<-Guideline_KRL$DrugName

# Validate the Guidelines..................................................

identical(rownames(CCLE_demeter_AML), rownames(Mut_AML_DEMETER2)) # TRUE


DEMETER2Score<-data.frame(CL=rownames(Mut_AML_DEMETER2),
                          TreeMut=NA,
                          TreeMutSqrt=NA,
                          MnLasso=NA,
                          BOSO=NA,
                          Lasso=NA,
                          MOM=NA, 
                          KRL=NA,
                          ORACLE=NA)

for(i in 1:nrow(DEMETER2Score)){
  CL<-DEMETER2Score$CL[i]
  
  # TreeMut
  dtreeMut<-Guideline_TreeMut[i]
  
  targets_dtreeMut<-unlist(strsplit(DT$Target[which(DT$Drug==dtreeMut)], split=", "))
  
  DEMETER2Score$TreeMut[i]<-mean(CCLE_demeter_AML[CL,na.omit(match(targets_dtreeMut,colnames(CCLE_demeter_AML)))])
  
  # TreeMutSqrt
  dtreeMutSQrt<-Guideline_TreeMutSqrt[i]
  
  targets_dtreeMutSQRT<-unlist(strsplit(DT$Target[which(DT$Drug==dtreeMutSQrt)], split=", "))
  
  DEMETER2Score$TreeMutSqrt[i]<-mean(CCLE_demeter_AML[CL,na.omit(match(targets_dtreeMutSQRT,colnames(CCLE_demeter_AML)))])
  
  
  # Multinomial
  dMultinomial<-Guideline_Multinomial[i]
  
  targets_dMultinomial<-unlist(strsplit(DT$Target[which(DT$Drug==dMultinomial)], split=", "))
  
  DEMETER2Score$MnLasso[i]<-mean(CCLE_demeter_AML[CL,na.omit(match(targets_dMultinomial,colnames(CCLE_demeter_AML)))])
  
  
  # BOSO
  dBOSO<-Guideline_BOSO[i]
  targets_dBOSO<-unlist(strsplit(DT$Target[which(DT$Drug==dBOSO)], split=", "))
  
  DEMETER2Score$BOSO[i]<-mean(CCLE_demeter_AML[CL,na.omit(match(targets_dBOSO,colnames(CCLE_demeter_AML)))])
  
  
  # Lasso
  dLasso<-Guideline_Lasso[i]
  targets_dLasso<-unlist(strsplit(DT$Target[which(DT$Drug==dLasso)], split=", "))
  
  DEMETER2Score$Lasso[i]<-mean(CCLE_demeter_AML[CL,na.omit(match(targets_dLasso,colnames(CCLE_demeter_AML)))])
  
  #MOM
  dMOM<-Guideline_MOM[i]
  targets_dMOM<-unlist(strsplit(DT$Target[which(DT$Drug==dMOM)], split=", "))
  
  DEMETER2Score$MOM[i]<-mean(CCLE_demeter_AML[CL,na.omit(match(targets_dMOM,colnames(CCLE_demeter_AML)))])
  
  #KRL
  dKRL<-Guideline_KRL[i]
  targets_dKRL<-unlist(strsplit(DT$Target[which(DT$Drug==dKRL)], split=", "))
  
  DEMETER2Score$KRL[i]<-mean(CCLE_demeter_AML[CL,na.omit(match(targets_dKRL,colnames(CCLE_demeter_AML)))])
  
  #ORACLE
  DEMETER2Score$ORACLE[i]<-min(CCLE_demeter_AML[CL,])
  
}

rownames(DEMETER2Score)<-DEMETER2Score$CL
boxplot(DEMETER2Score[,c(-1)],main="DEMETER 2", ylab="DEMETER 2 SCORE")



# 2. GDSC -----------------------------------------------------------------------------------------------

GDSC1_ESS<-read.csv("./data/input/LAML_IC_Wed Nov 24 10_31_42 2021.csv")
GDSC1_MUT<-read.csv("./data/input/LAML_Genetic_feature_variant_Wed Nov 24 GDSC1.csv")



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
                  ncol=ncol(Mut_AML_DEMETER2))

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



# Pre-process Drug Matrix

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

# write.csv(Mutations_AML_GDSC, file="./XAI methods/mut_Matrix_GDSC.csv")


# Predict with the Models.........................................................
source("./MOM/predictMOM_ALL.R")

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
                      ORACLE=NA, 
                      UpperRef=NA)

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
  
  #UpperRef
  GDSC_IC50$UpperRef[i]<-GDSC_AML[CL,which.min(colMeans(GDSC_AML))]
  
}

rownames(GDSC_IC50)<-GDSC_IC50$CL
boxplot(GDSC_IC50[,c(-1)],main="GDSC", ylab="GDSC IC50")




# 3. CERES----------------------------------------------------------------------

# LOADING CERES DATA

ceresAML<-read.csv("./data/input/wang_ceres_gene_effects.csv",row.names = 1)



# PRE-Processing Ceres Data

CCLE_rawmutation<-read.delim("./data/input/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", comment.char="#",row.names = 1)
# Update with CIMA mutations
# Matching with names
CCLE_rawmutation["TP53_MUT",c("MV411_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","NOMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","OCIAML2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE")]   <- 1
CCLE_rawmutation["KDM6A_MUT", "HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"]  <- 1
# CCLE_rawmutation["CDKN2A_MUT","HL60HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"] <- 1
CCLE_rawmutation["MLL_MUT","KASUMI1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"]    <- 1
CCLE_rawmutation["MLL_MUT",c("MOLM13_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","OCIAML2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE")]    <- 0
CCLE_rawmutation["FLT3_MUT",c("MOLM13_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","MV411_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE")]   <- 1
# CCLE_rawmutation["FLT3_MUT","OCIAML2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"]   <- 0
CCLE_rawmutation["KRAS_MUT","NOMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"]   <- 1

CCLE_rawmutation<-cbind(CCLE_rawmutation, 0)
colnames(CCLE_rawmutation)[ncol(CCLE_rawmutation)]<-"CBFB--MYH11"
CCLE_rawmutation["ME1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","CBFB--MYH11"]   <- 1

CCLE_mutation_data_ceres<-CCLE_rawmutation[grep("_MUT", rownames(CCLE_rawmutation)), na.omit(match(colnames(ceresAML), colnames(CCLE_rawmutation)))]


# CCLE_mutation_data_ceres<-CCLE_mutation_data_ceres[-c((nrow(CCLE_mutation_data_ceres)-7):nrow(CCLE_mutation_data_ceres)),]
rownames(CCLE_mutation_data_ceres)<-sapply(rownames(CCLE_mutation_data_ceres),FUN=function(X){return(unlist(strsplit(X, split="_"))[1])})
# rm(demeter2_mut, demeter2_scores, demeter2_fus)

CERES_Mutations<-t(CCLE_mutation_data_ceres)
CERES_Mutations<-cbind(CERES_Mutations,
                       Mut_AML_DEMETER2[match(rownames(CERES_Mutations), rownames(Mut_AML_DEMETER2)),
                                        which(!(colnames(Mut_AML_DEMETER2) %in% colnames(CERES_Mutations)))])


CERES_Mutations[is.na(CERES_Mutations)]<-0
CERES_Mutations<-CERES_Mutations[,match(colnames(mut2), colnames(CERES_Mutations))]





ceresAML<-t(ceresAML)
ceresAML<-t(t(ceresAML)-colMeans(ceresAML))

identical(rownames(CERES_Mutations), rownames(ceresAML))# TRUE

# write.csv(CERES_Mutations, file="./XAI methods/mut_Matrix_CERES.csv")

# MODEL PREDICTION..............................................................
identical(colnames(CERES_Mutations), colnames(mut2))
# [1] TRUE

source("./MOM/predictMOM_ALL.R")

Guideline_TreeMutSqrt<-predictTreeMut(treeMut1_All,C, CERES_Mutations)
Guideline_TreeMutSqrt<-colnames(C)[Guideline_TreeMutSqrt]


Guideline_BOSO<-predictBOSO(BOSOMut_All,CERES_Mutations)
Guideline_BOSO<-colnames(C)[Guideline_BOSO]

Guideline_TreeMut<-predictTreeMut(treeMut_All,C, CERES_Mutations)
Guideline_TreeMut<-colnames(C)[Guideline_TreeMut]

Guideline_Multinomial<-predictMnLasso(MnLassoMut_All, CERES_Mutations)
Guideline_Multinomial<-colnames(C)[Guideline_Multinomial]

Guideline_Lasso<-predictLasso(LassoMut_All, CERES_Mutations)
Guideline_Lasso<-colnames(C)[Guideline_Lasso]

Guideline_MOM<-predictMOM_ALL(CERES_Mutations)
Guideline_MOM<-colnames(C)[Guideline_MOM$Treatment]

Guideline_KRL<-read_csv("data/output/KRL_CERES_pred.csv")
Guideline_KRL$DrugName[which(!(Guideline_KRL$DrugName %in% colnames(C)))]<-NA
Guideline_KRL<-Guideline_KRL$DrugName

# CERES Comparison..............................................................



identical(rownames(ceresAML), rownames(CERES_Mutations)) # TRUE


CERESScore<-data.frame(CL=rownames(CERES_Mutations),
                       TreeMut=NA,
                       TreeMutSqrt=NA,
                       MnLasso=NA,
                       BOSO=NA,
                       Lasso=NA,
                       MOM=NA,
                       KRL=NA,
                       ORACLE=NA)

for(i in 1:nrow(CERESScore)){
  CL<-CERESScore$CL[i]
  
  # TreeMut
  dtreeMut<-Guideline_TreeMut[i]
  
  targets_dtreeMut<-unlist(strsplit(DT$Target[which(DT$Drug==dtreeMut)], split=", "))
  
  CERESScore$TreeMut[i]<-mean(ceresAML[CL,na.omit(match(targets_dtreeMut,colnames(ceresAML)))])
  
  # TreeMutSqrt
  dtreeMutSQrt<-Guideline_TreeMutSqrt[i]
  
  targets_dtreeMutSQRT<-unlist(strsplit(DT$Target[which(DT$Drug==dtreeMutSQrt)], split=", "))
  
  CERESScore$TreeMutSqrt[i]<-mean(ceresAML[CL,na.omit(match(targets_dtreeMutSQRT,colnames(ceresAML)))])
  
  
  # Multinomial
  dMultinomial<-Guideline_Multinomial[i]
  
  targets_dMultinomial<-unlist(strsplit(DT$Target[which(DT$Drug==dMultinomial)], split=", "))
  
  CERESScore$MnLasso[i]<-mean(ceresAML[CL,na.omit(match(targets_dMultinomial,colnames(ceresAML)))])
  
  
  # BOSO
  dBOSO<-Guideline_BOSO[i]
  targets_dBOSO<-unlist(strsplit(DT$Target[which(DT$Drug==dBOSO)], split=", "))
  
  CERESScore$BOSO[i]<-mean(ceresAML[CL,na.omit(match(targets_dBOSO,colnames(ceresAML)))])
  
  
  # Lasso
  dLasso<-Guideline_Lasso[i]
  targets_dLasso<-unlist(strsplit(DT$Target[which(DT$Drug==dLasso)], split=", "))
  
  CERESScore$Lasso[i]<-mean(ceresAML[CL,na.omit(match(targets_dLasso,colnames(ceresAML)))])
  
  #MOM
  dMOM<-Guideline_MOM[i]
  targets_dMOM<-unlist(strsplit(DT$Target[which(DT$Drug==dMOM)], split=", "))
  
  CERESScore$MOM[i]<-mean(ceresAML[CL,na.omit(match(targets_dMOM,colnames(ceresAML)))])
  
  #KRL
  dKRL<-Guideline_KRL[i]
  targets_dKRL<-unlist(strsplit(DT$Target[which(DT$Drug==dKRL)], split=", "))
  
  CERESScore$KRL[i]<-mean(ceresAML[CL,na.omit(match(targets_dKRL,colnames(ceresAML)))])
  
  
  #ORACLE
  CERESScore$ORACLE[i]<-as.numeric(min(ceresAML[CL,]))
  
}

rownames(CERESScore)<-CERESScore$CL
boxplot(CERESScore[,c(-1)],main="CERES", ylab="CERES SCORE")



# Final Plot Figure--------------------------------------------------------------

Validation_plot<-data.frame(Dataset=rep("DEMETER2", nrow(DEMETER2Score)*8), 
                            CellLine=rep(rownames(DEMETER2Score),8),
                            Method=rep(colnames(DEMETER2Score[,-1]), each=nrow(DEMETER2Score)), 
                            IC50=unlist(DEMETER2Score[,-1]))


Validation_plot<-rbind(Validation_plot, data.frame(Dataset=rep("CERES", nrow(CERESScore)*8), 
                                                   CellLine=rep(rownames(CERESScore),8),
                                                   Method=rep(colnames(CERESScore[,-1]), each=nrow(CERESScore)), 
                                                   IC50=unlist(CERESScore[,-1])))


Validation_plot<-rbind(Validation_plot, data.frame(Dataset=rep("GDSC", nrow(GDSC_IC50)*9), 
                                                   CellLine=rep(rownames(GDSC_IC50),9),
                                                   Method=rep(colnames(GDSC_IC50[,-1]), each=nrow(GDSC_IC50)), 
                                                   IC50=unlist(GDSC_IC50[,-1])))


rownames(Validation_plot)<-1:nrow(Validation_plot)

mycomparisons<-combn(unique(Validation_plot$Method),2)
mycomparisons<-as.list(data.frame(mycomparisons))
mycomparisons<-mycomparisons[-grep("ORACLE",mycomparisons)]


Plotting<-Validation_plot[Validation_plot$Dataset=="DEMETER2",]


library(ggsci)

Plotting$Method<-factor(Plotting$Method, levels=c("ORACLE", "MOM", "BOSO", "TreeMutSqrt", "Lasso", "TreeMut", "MnLasso", "KRL"))

g1_demeter<-ggplot(Plotting, aes(x=Method, y=IC50, fill=Method))+geom_violin()+
  theme_bw()+scale_fill_npg()+ylab("DEMETER 2 SCORE")+ggtitle("DEMETER 2")

g2_demeter<-ggplot(Plotting, aes(x=Method, y=IC50, fill=Method))+geom_boxplot()+
  theme_bw()+scale_fill_npg()+ylab("DEMETER 2 SCORE")+ggtitle("DEMETER 2")
#+stat_compare_means(method = "t.test",comparisons = mycomparisons[-grep("BOSO",mycomparisons)])



Plotting<-Validation_plot[Validation_plot$Dataset=="CERES",]


library(ggsci)

Plotting$Method<-factor(Plotting$Method, levels=c("ORACLE", "MOM", "BOSO", "TreeMutSqrt", "Lasso", "TreeMut", "MnLasso", "KRL"))

g1_ceres<-ggplot(Plotting, aes(x=Method, y=IC50, fill=Method))+geom_violin()+
  theme_bw()+scale_fill_npg()+ylab("CERES SCORE")+ggtitle("CERES")

g2_ceres<-ggplot(Plotting, aes(x=Method, y=IC50, fill=Method))+geom_boxplot()+
  theme_bw()+scale_fill_npg()+ylab("CERES SCORE")+ggtitle("CERES")









Plotting<-Validation_plot[Validation_plot$Dataset=="GDSC",]


library(ggsci)

Plotting$Method<-factor(Plotting$Method, levels=c("ORACLE", "MOM", "BOSO", "TreeMutSqrt", "Lasso", "TreeMut", "MnLasso", "KRL", "UpperRef"))

g1_gdsc<-ggplot(Plotting, aes(x=Method, y=IC50, fill=Method))+geom_violin()+
  theme_bw()+scale_fill_npg()+ylab("GDSC IC50")+ggtitle("GDSC")

g2_gdsc<-ggplot(Plotting, aes(x=Method, y=IC50, fill=Method))+geom_boxplot()+
  theme_bw()+scale_fill_npg()+ylab("GDSC IC50")+ggtitle("GDSC")

g1_gdsc
g2_gdsc+theme(legend.position="none", text = element_text(size = 15))



wilcox.test(Plotting$IC50[Plotting$Method=="UpperRef"], 
            Plotting$IC50[Plotting$Method=="MOM"], alternative="greater", paired = T)
wilcox.test(Plotting$IC50[Plotting$Method=="UpperRef"], 
            Plotting$IC50[Plotting$Method=="BOSO"], alternative="greater", paired = T)
wilcox.test(Plotting$IC50[Plotting$Method=="UpperRef"], 
            Plotting$IC50[Plotting$Method=="TreeMutSqrt"], alternative="greater", paired = T)
wilcox.test(Plotting$IC50[Plotting$Method=="UpperRef"], 
            Plotting$IC50[Plotting$Method=="Lasso"], alternative="greater", paired = T)
wilcox.test(Plotting$IC50[Plotting$Method=="UpperRef"], 
            Plotting$IC50[Plotting$Method=="TreeMut"], alternative="greater", paired = T)
wilcox.test(Plotting$IC50[Plotting$Method=="UpperRef"], 
            Plotting$IC50[Plotting$Method=="MnLasso"], alternative="greater", paired = T)
wilcox.test(Plotting$IC50[Plotting$Method=="UpperRef"], 
            Plotting$IC50[Plotting$Method=="KRL"], alternative="greater", paired = T)



wilcox.test(Plotting$IC50[Plotting$Method=="MOM"], 
            Plotting$IC50[Plotting$Method=="MOM"], alternative="less")
wilcox.test(Plotting$IC50[Plotting$Method=="MOM"], 
            Plotting$IC50[Plotting$Method=="BOSO"], alternative="less")
wilcox.test(Plotting$IC50[Plotting$Method=="MOM"], 
            Plotting$IC50[Plotting$Method=="TreeMutSqrt"], alternative="less")
wilcox.test(Plotting$IC50[Plotting$Method=="MOM"], 
            Plotting$IC50[Plotting$Method=="Lasso"], alternative="less")
wilcox.test(Plotting$IC50[Plotting$Method=="MOM"], 
            Plotting$IC50[Plotting$Method=="TreeMut"], alternative="less")
wilcox.test(Plotting$IC50[Plotting$Method=="MOM"], 
            Plotting$IC50[Plotting$Method=="MnLasso"], alternative="less")
wilcox.test(Plotting$IC50[Plotting$Method=="MOM"], 
            Plotting$IC50[Plotting$Method=="KRL"], alternative="less")

### 
# Join Figures Validation
####

library(ggpubr)


ggarrange(g1_demeter, g1_ceres, g1_gdsc, ncol=2, nrow=2, common.legend = T,  legend=c("right"))

ggarrange(g2_demeter, g2_ceres, g2_gdsc, ncol=2, nrow=2, common.legend = T, legend=c("right"))

ggarrange(g2_demeter,g2_ceres, g1_demeter,  g1_ceres, ncol=2,nrow=2,common.legend = T, legend=c("right") )

################################################################################
# Comparison Drug Recommended vs. Drug not recommended
################################################################################


TratamientoMOM_All<-predictMOM_ALL(mut2)
identical(TratamientoMOM_All$CLs, rownames(C)) # TRUE

TratamientoKRL_All<- read_csv("data/output/KRL_results_all_samples_prediction.csv")
TratamientoKRL_All<-TratamientoKRL_All[which(TratamientoKRL_All$Patients %in% rownames(mut2)),]
TratamientoKRL_All$DrugName[which(!(TratamientoKRL_All$DrugName %in% colnames(C)))]<-NA

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


### MOM

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


### BOSO

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

#### TREE MUT

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

### TREE MUT SQRT

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

#### MULTINOMIAL

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


#### LASSO

Plot_aux<-Validation_Beat[Validation_Beat$Method=="LassoMut",]

Lasso_plots<-NULL

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
  
  Lasso_plots<-rbind(Lasso_plots,plotDrug)
  
}

gg6_validation_beat<-ggplot(Lasso_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_npg()+
  stat_compare_means(vjust=1)+ylab("IC50*")+ggtitle("LASSO")


#### KRL

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




## COMBINED

ggarrange(gg1_validation_beat, gg2_validation_beat, gg3_validation_beat, 
          gg4_validation_beat, gg5_validation_beat, gg6_validation_beat, gg7_validation_beat,
          ncol=2, nrow=4)


gg1_validation_beat+theme(text=element_text(size=15), legend.position = "none")


gg2_validation_beat

gg3_validation_beat+theme(text=element_text(size=15), legend.position = "none")

gg4_validation_beat+theme(text=element_text(size=15), legend.position = "none")

gg5_validation_beat

gg6_validation_beat

gg7_validation_beat



ggarrange(gg1_validation_beat, gg3_validation_beat,ncols=2, nrows=0, common.legend = T)
# 2. GDSC-----------------------------------------------------------------------
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

# GDSC_AML<-Out_gdsc$data
# 
# 
# drugsGDSC<-sapply(colnames(C), FUN=function(X){return(unlist(strsplit(X, split=" "))[1])})
# 
# GDSC_AML<-GDSC_AML[,na.omit(match(drugsGDSC, colnames(GDSC_AML)))]
# 
# colnamesGDSC<-colnames(C)[match(colnames(GDSC_AML),drugsGDSC)]
# 
# 
# colnames(GDSC_AML)<-colnamesGDSC
# 
# 
# GDSC_AML<-GDSC_AML[which(rownames(GDSC_AML) %in% rownames(Mutations_AML_GDSC)),]
# GDSC_AML<-GDSC_AML[match(rownames(Mutations_AML_GDSC), rownames(GDSC_AML)),]
# 
# 
# identical(rownames(GDSC_AML), rownames(Mutations_AML_GDSC))
# 



# Validation_GDSC<-data.frame(patients=rep(rownames(GDSC_AML),6),
#                             drug=c(Guideline_TreeMut,
#                                    Guideline_TreeMutSqrt,
#                                    Guideline_Multinomial,
#                                    Guideline_Lasso,
#                                    Guideline_BOSO,
#                                    Guideline_MOM), 
#                           IC50=c(GDSC_AML[cbind(1:nrow(GDSC_AML), Guideline_TreeMut)],
#                                  GDSC_AML[cbind(1:nrow(GDSC_AML), Guideline_TreeMutSqrt)],
#                                  GDSC_AML[cbind(1:nrow(GDSC_AML), Guideline_Multinomial)],
#                                  GDSC_AML[cbind(1:nrow(GDSC_AML), Guideline_Lasso)],
#                                  GDSC_AML[cbind(1:nrow(GDSC_AML), Guideline_BOSO)],
#                                  GDSC_AML[cbind(1:nrow(GDSC_AML), Guideline_MOM)]),
#                             Method=rep(unique(Tratamiento_plot$Method[!(Tratamiento_plot$Method=="ORACLE")]), each=nrow(GDSC_AML)))


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


### MOM

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


### BOSO

Plot_aux<-Validation_GDSC[Validation_GDSC$Method=="BOSOMut",]

BOSO_plots<-NULL

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
  
  
  BOSO_plots<-rbind(BOSO_plots,plotDrug)
  
}

gg2_validation_gdsc<-ggplot(BOSO_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_manual(values=c("#999999", "#E69F00"))+
  stat_compare_means(vjust = 1)+ylab("IC50*")+ggtitle("BOSO")

#### TREE MUT

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

### TREE MUT SQRT

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

#### MULTINOMIAL

Plot_aux<-Validation_GDSC[Validation_GDSC$Method=="Multinomial",]

MnLasso_plots<-NULL

for(j in 1:length(unique(Plot_aux$drug))){
  d<-unique(Plot_aux$drug)[j]
  ind<-Plot_aux$drug==d
  plotDrug<-data.frame(Drug=j,
                       drugname=d,
                       Patients=Plot_aux$patients,
                       IC50=NA,
                       Recommended="NO")
  if(d %in% colnames(GDSC_AML)){
    plotDrug$Recommended[ind]<-"Yes"
    plotDrug$IC50[ind]<-GDSC_AML[ind,d]
    plotDrug$IC50[!ind]<-GDSC_AML[!ind, d]
    plotDrug$Recommended[!ind]<-"No"
  }
 
  
  MnLasso_plots<-rbind(MnLasso_plots,plotDrug)
  
}

gg5_validation_gdsc<-ggplot(MnLasso_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_manual(values=c("#999999","#E69F00"))+
  stat_compare_means(vjust=1)+ylab("IC50*")+ggtitle("Multinomial")


#### LASSO

Plot_aux<-Validation_GDSC[Validation_GDSC$Method=="LassoMut",]

Lasso_plots<-NULL

for(j in 1:length(unique(Plot_aux$drug))){
  d<-unique(Plot_aux$drug)[j]
  ind<-Plot_aux$drug==d
  plotDrug<-data.frame(Drug=j,
                       drugname=d,
                       Patients=Plot_aux$patients,
                       IC50=NA,
                       Recommended="NO")
  if(d %in% colnames(GDSC_AML)){
    plotDrug$Recommended[ind]<-"Yes"
    plotDrug$IC50[ind]<-GDSC_AML[ind,d]
    plotDrug$IC50[!ind]<-GDSC_AML[!ind, d]
    plotDrug$Recommended[!ind]<-"NO"
  }
  

  
  Lasso_plots<-rbind(Lasso_plots,plotDrug)
  
}

gg6_validation_gdsc<-ggplot(Lasso_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_manual(values=c("#999999", "#E69F00"))+
  stat_compare_means(vjust=1)+ylab("IC50*")+ggtitle("LASSO")


#### KRL

Plot_aux<-Validation_GDSC[Validation_GDSC$Method=="KRL",]

KRL_plots<-NULL

for(j in 1:length(unique(Plot_aux$drug))){
  d<-unique(Plot_aux$drug)[j]
  if(is.na(d)){next}
  ind<-na.omit(Plot_aux$drug==d)
  plotDrug<-data.frame(Drug=j,
                       drugname=d,
                       Patients=Plot_aux$patients,
                       IC50=NA,
                       Recommended="NO")
  if(d %in% colnames(GDSC_AML)){
    plotDrug$Recommended[ind]<-"Yes"
    plotDrug$IC50[ind]<-GDSC_AML[ind,d]
    plotDrug$IC50[!ind]<-GDSC_AML[!ind, d]
    plotDrug$Recommended[!ind]<-"NO"
  }
  
  
  
  KRL_plots<-rbind(KRL_plots,plotDrug)
  
}

gg7_validation_gdsc<-ggplot(KRL_plots, aes(x=Recommended, y=IC50, fill=Recommended))+
  geom_boxplot()+facet_wrap(~drugname)+theme_bw()+scale_fill_manual(values=c("#999999", "#E69F00"))+
  stat_compare_means(vjust=1)+ylab("IC50*")+ggtitle("KRL")


## COMBINED

ggarrange(gg1_validation_gdsc, gg2_validation_gdsc, gg3_validation_gdsc,
                gg4_validation_gdsc, gg5_validation_gdsc, gg6_validation_gdsc,gg7_validation_gdsc,
                ncol=2, nrow=4)


gg1_validation_gdsc+theme(text=element_text(size=15), legend.position = "none")


gg2_validation_gdsc

gg3_validation_gdsc+theme(text=element_text(size=15), legend.position = "none")

gg4_validation_gdsc+theme(text=element_text(size=15), legend.position = "none")

gg5_validation_gdsc

gg6_validation_gdsc

gg7_validation_gdsc


################################################################################
#Comparison Methods Mut vs. Expression in BeatAML and GDSC
################################################################################
# 1. Beat AML 5-CV CODE--------------------------------------------------------------------------
# 
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

BOSO <- trainBOSO(Expression, C)
TratamientoBOSO <- predictBOSO(BOSO, Expression)

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

summary(C[cbind(1:nrow(C), TratamientoOracle)])
summary(C[cbind(1:nrow(C), TratamientoTree)])
summary(C[cbind(1:nrow(C), TratamientoMnLasso)])
summary(C[cbind(1:nrow(C), TratamientoLasso)])
summary(C[cbind(1:nrow(C), TratamientoBOSO)])
summary(C[cbind(1:nrow(C), TratamientoTreeMut)])
summary(C[cbind(1:nrow(C), TratamientoMnLassoMut)])
summary(C[cbind(1:nrow(C), TratamientoLassoMut)])
summary(C[cbind(1:nrow(C), TratamientoBOSOMut)])

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
boxplot(Tratamientos)
# We can test each of the methods to check for statistiscal significance
wilcox.test(Tratamientos[,"BOSO"], Tratamientos[,"BOSOMut"], paired = T, alternative = "greater")
wilcox.test(Tratamientos[,"Tree"], Tratamientos[,"TreeMut"], paired = T, alternative = "greater")
wilcox.test(Tratamientos[,"MnLasso"], Tratamientos[,"MnLassoMut"], paired = T, alternative = "greater")
wilcox.test(Tratamientos[,"Lasso"], Tratamientos[,"LassoMut"], paired = T, alternative = "greater")




# Plotting----------------------------------------------------------------------

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


# 2. GDSC-----------------------------------------------------------------------

# Calculate Overall Treatment with mutated and Expression Models
# Trees
cat("Trees\n")
tic()
tree_Ex_ALL <- trainTree(Expression, C, minbucket = 20)
TratamientoTree_Ex_ALL <- predictTree(tree_Ex_ALL,C, Expression)
toc()

# Trees sqrt
cat("Trees sqrt\n")
tic()
tree1_Ex_ALL <- trainTree(Expression, (C)^.5, minbucket = 20)
TratamientoTree1_Ex_ALL<- predictTree(tree1_Ex_ALL,C, Expression)
toc()

# Using MnLasso
cat("Multinomial Lasso\n")
tic()
MnLasso_Ex_ALL <- trainMnLasso(Expression, Y)
TratamientoMnLasso_Ex_ALL<- predictMnLasso(MnLasso_Ex_ALL,Expression)
toc()

# Using Lasso
cat("Lasso\n")
tic()
Lasso_Ex_ALL <- trainLasso(Expression, C)
TratamientoLasso_Ex_ALL <- predictMnLasso(Lasso_Ex_ALL,Expression)
toc()

# TreesMut
cat("TreesMut\n")
treeMut_Ex_ALL <- trainTreeMut(X, C, minbucket = 20)
TratamientoTreeMut_Ex_ALL <- predictTreeMut(treeMut_Ex_ALL,C, X)


# TreesMut sqrt
cat("TreesMut sqrt\n")
treeMut1_Ex_ALL <- trainTreeMut(X, (C)^.5, minbucket = 20)
TratamientoTreeMut1_Ex_ALL <- predictTreeMut(treeMut1_Ex_ALL,C, X)

# Using MnLasso
cat("Multinomial Lasso Mut\n")
MnLassoMut_Ex_ALL <- trainMnLasso(X, Y)
TratamientoMnLassoMut_Ex_ALL <- predictMnLasso(MnLassoMut_Ex_ALL,X)

# Using Lasso
cat("Lasso Mut\n")
LassoMut_Ex_ALL <- trainLasso(X, C)
TratamientoLassoMut_Ex_ALL <- predictLasso(LassoMut_Ex_ALL,X)

cat("BOSO \n")
tic()
BOSO_Ex_ALL <- trainBOSO(Expression, C,maxVarsBlock = 10, standardize=F, IC = "eBIC")
TratamientoBOSO <- predictBOSO(BOSO_Ex_ALL, Expression)
toc()

BOSOMut_Ex_ALL <- trainBOSO(X, C,maxVarsBlock = 10, standardize=F, IC = "eBIC")
TratamientoBOSOMut_Ex_ALL <- predictBOSO(BOSOMut_Ex_ALL, X)


#### Predict in GDSC

load("data/input/CCLE_tpm_gn_allGenes.RData")
source("./MOM/predictMOM_ALL.R")
library(ensembldb)
library(EnsDb.Hsapiens.v79)

CCLE_tpm_gn<-as.data.frame(CCLE_tpm_gn)
rownames(CCLE_tpm_gn)<-CCLE_tpm_gn$Gene_ID
genes<-rownames(CCLE_tpm_gn)

symbol_genes<-ensembldb::select(EnsDb.Hsapiens.v79, keys= genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

dummy <- merge(symbol_genes, data.frame(GENEID = genes, EVENT1=genes))

identical(dummy$EVENT1[1:60245], rownames(CCLE_tpm_gn)[1:60245])# TRUE

CCLE_tpm_gn<-CCLE_tpm_gn[dummy$GENEID,]

CCLE_tpm_gn$Gene_ID<-dummy$SYMBOL
genes_to_select<-which(CCLE_tpm_gn$Gene_ID %in% colnames(Expression))

CCLE_tpm_GDSC<-CCLE_tpm_gn[genes_to_select,which(colnames(CCLE_tpm_gn)%in% rownames(Mutations_AML_GDSC))]
CCLE_tpm_GDSC<-as.data.frame(CCLE_tpm_GDSC)
CCLE_tpm_GDSC$SYMBOL<-CCLE_tpm_gn$Gene_ID[genes_to_select]


rowstoadd<-which(!(colnames(Expression) %in%  CCLE_tpm_gn$Gene_ID))
add_mat<-matrix(0, nrow=length(rowstoadd), ncol=ncol(CCLE_tpm_GDSC))
rownames(add_mat)<-colnames(Expression)[rowstoadd]
colnames(add_mat)<-colnames(CCLE_tpm_GDSC)
add_mat<-as.data.frame(add_mat)
add_mat$SYMBOL<-colnames(Expression)[rowstoadd]


CCLE_tpm_GDSC<-rbind(CCLE_tpm_GDSC, add_mat)

CCLE_tpm_GDSC<-CCLE_tpm_GDSC[match(colnames(Expression),CCLE_tpm_GDSC$SYMBOL),]

identical(CCLE_tpm_GDSC$SYMBOL, colnames(Expression)) # TRUE
rownames(CCLE_tpm_GDSC)<-CCLE_tpm_GDSC$SYMBOL

CCLE_tpm_GDSC<-CCLE_tpm_GDSC[,-which(colnames(CCLE_tpm_GDSC)=="SYMBOL")]


identical(rownames(CCLE_tpm_GDSC), colnames(Expression))# TRUE

CCLE_tpm_GDSC<-CCLE_tpm_GDSC[,na.omit(match(rownames(Mutations_AML_GDSC),colnames(CCLE_tpm_GDSC)))]
Mutations_AML_GDSC<-Mutations_AML_GDSC[match(colnames(CCLE_tpm_GDSC), rownames(Mutations_AML_GDSC)),]
GDSC_AML<-GDSC_AML[match(rownames(Mutations_AML_GDSC),rownames(GDSC_AML)),]

identical(colnames(Mutations_AML_GDSC), colnames(mut2)) #TRUE
identical(rownames(GDSC_AML), rownames(Mutations_AML_GDSC))#TRUE
identical(rownames(GDSC_AML), colnames(CCLE_tpm_GDSC)) # TRUE

CCLE_tpm_GDSC<-t(CCLE_tpm_GDSC)
CCLE_tpm_GDSC<-as.matrix(CCLE_tpm_GDSC)

# GET MUT and EX Guidelines in GDSC............................................

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

# Guideline_MOM<-predictMOM_ALL(Mutations_AML_GDSC)
# Guideline_MOM<-colnames(C)[Guideline_MOM$Treatment]



Guideline_Ex_TreeSqrt<-predictTree(tree1_Ex_ALL,C, CCLE_tpm_GDSC)
Guideline_Ex_TreeSqrt<-colnames(C)[Guideline_Ex_TreeSqrt]


Guideline_Ex_BOSO<-predictBOSO(BOSO_Ex_ALL,CCLE_tpm_GDSC)
Guideline_Ex_BOSO<-colnames(C)[Guideline_Ex_BOSO]

Guideline_Ex_Tree<-predictTree(tree_Ex_ALL,C, CCLE_tpm_GDSC)
Guideline_Ex_Tree<-colnames(C)[Guideline_Ex_Tree]

Guideline_Ex_Multinomial<-predictMnLasso(MnLasso_Ex_ALL, CCLE_tpm_GDSC)
Guideline_Ex_Multinomial<-colnames(C)[Guideline_Ex_Multinomial]

Guideline_Ex_Lasso<-predictLasso(Lasso_Ex_ALL, CCLE_tpm_GDSC)
Guideline_Ex_Lasso<-colnames(C)[Guideline_Ex_Lasso]

# Guideline_MOM<-predictMOM_ALL(CCLE_tpm_GDSC)
# Guideline_MOM<-colnames(C)[Guideline_MOM$Treatment]



# Validate the Guidelines..................................................

identical(rownames(GDSC_AML), rownames(Mutations_AML_GDSC)) # TRUE


GDSC_IC50<-data.frame(CL=rownames(Mutations_AML_GDSC),
                      TreeMut=NA,
                      TreeMutSqrt=NA,
                      MnLassoMut=NA,
                      BOSOMut=NA,
                      LassoMut=NA,
                      ORACLE=NA,
                      TreeEx=NA, 
                      TreeExSqrt=NA, 
                      BOSO=NA, 
                      MnLasso=NA, 
                      Lasso=NA)

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
  
  #TreeEX
  dTreeEx<-Guideline_Ex_Tree[i]
  if(dTreeEx %in% colnames(GDSC_AML)){
  GDSC_IC50$TreeEx[i]<-GDSC_AML[CL,dTreeEx]
  }
  
  #TreeSqrtEX
  dTreeExSqrt<-Guideline_Ex_TreeSqrt[i]
  if(dTreeExSqrt %in% colnames(GDSC_AML)){
  GDSC_IC50$TreeExSqrt[i]<-GDSC_AML[CL,dTreeExSqrt]
  }
  
  #BOSO
  dBOSOEx<-Guideline_Ex_BOSO[i]
  if(dBOSOEx %in% colnames(GDSC_AML)){
  GDSC_IC50$BOSO[i]<-GDSC_AML[CL,dBOSOEx]
  }
  
  #Lasso
  dLassoEx<-Guideline_Ex_Lasso[i]
  if(dLassoEx %in% colnames(GDSC_AML)){
  GDSC_IC50$Lasso[i]<-GDSC_AML[CL,dLassoEx]
  }
  
  #MnLasso
  dMnLassoEx<-Guideline_Ex_Multinomial[i]
  if(dMnLassoEx %in% colnames(GDSC_AML)){
  GDSC_IC50$MnLasso[i]<-GDSC_AML[CL,dMnLassoEx]
  }
  
  #ORACLE
  GDSC_IC50$ORACLE[i]<-min(GDSC_AML[CL,])
  
}

rownames(GDSC_IC50)<-GDSC_IC50$CL
boxplot(GDSC_IC50[,c(-1)],main="GDSC", ylab="GDSC IC50")


# Plotting----------------------------------------------------
# GDSC_IC50<-GDSC_IC50[,-which(colnames(GDSC_IC50)=="CL")]

plot_Ex_Mut<-data.frame(Method=rep(colnames(GDSC_IC50[,c(-1)]), each=nrow(GDSC_IC50[,c(-1)])))
plot_Ex_Mut$IC50<-unlist(as.data.frame(GDSC_IC50[,c(-1)]))
plot_Ex_Mut$Type<-"GE"
plot_Ex_Mut$Type[grep("Mut", plot_Ex_Mut$Method)]<-"Mut"


plot_Ex_Mut<-plot_Ex_Mut[order(plot_Ex_Mut$Method),]
plot_Ex_Mut$General_Method<-"unknown"
plot_Ex_Mut$General_Method[grep("ORACLE", plot_Ex_Mut$Method)]<-"ORACLE"
plot_Ex_Mut$General_Method[grep("BOSO", plot_Ex_Mut$Method)]<-"BOSO"
plot_Ex_Mut$General_Method[grep("Tree", plot_Ex_Mut$Method)]<-"Tree"
plot_Ex_Mut$General_Method[grep("Sqrt", plot_Ex_Mut$Method)]<-"TreeSqrt"
# plot_Ex_Mut$General_Method[grep("TreeSrqtMut", plot_Ex_Mut$Method)]<-"TreeSqrt"
plot_Ex_Mut$General_Method[grep("Lasso", plot_Ex_Mut$Method)]<-"Lasso"
plot_Ex_Mut$General_Method[grep("MnLasso", plot_Ex_Mut$Method)]<-"MnLasso"

# plot_Ex_Mut<-rbind(plot_Ex_Mut, data.frame(Method=NA, IC50=5, Type="Mut", General_Method="MOM"))

plot_Ex_Mut$Method<-factor(plot_Ex_Mut$Method, levels=c("ORACLE","BOSO","BOSOMut","TreeExSqrt","TreeMutSqrt","Lasso","LassoMut",
                                                        "TreeEx","TreeMut",
                                                        "MnLasso","MnLassoMut"))
plot_Ex_Mut$General_Method<-factor(plot_Ex_Mut$General_Method, levels=c("ORACLE", "BOSO", "TreeSqrt", "Lasso",
                                                                        "Tree","MnLasso"))


mycomparisons<-list(c("TreeEx", "TreeMut"), c("TreeExSqrt", "TreeMutSqrt"), c("Lasso", "LassoMut"), c("BOSO", "BOSOMut"),
                    c("MnLasso", "MnLassoMut"))
mypal = pal_npg("nrc")(9)
mypal<-mypal[c(1,3:7)]
gg_ex_mut<-ggplot(plot_Ex_Mut, aes(x=Method, y=IC50, fill=General_Method, alpha=Type))+geom_boxplot()+
  scale_fill_manual(values=mypal)+theme_bw()+ggtitle("GDSC")+ylab("IC50*")+labs(fill="General Method")+
  stat_compare_means(aes(Method),vjust=0, hjust=-1, comparisons=mycomparisons,label.y = 5)+scale_alpha_manual(values=c( 0.5, 1))

gg_ex_mut


################################################################################
# Variable Number
################################################################################

Variable_plot<-data.frame(Method=unique(Tratamiento_plot$Method[!(Tratamiento_plot$Method=="ORACLE")]),
                          Number_of_variables=NA)






### Lasso
Coeficientes <- coef(LassoMut_All, s = "lambda.min")
CoefM <- Matrix(0, nrow = nrow(Coeficientes[[1]]), ncol = length(Coeficientes))

for (n in 1:length(Coeficientes)) {
  CoefM[,n] <- Coeficientes[[n]]
  
}
colnames(CoefM) <- names(Coeficientes)
rownames(CoefM) <- rownames(Coeficientes[[1]])
image(CoefM)
nLasso<-sum(rowSums(abs(CoefM)) > 1e-4)-1
Variable_plot$Number_of_variables[Variable_plot$Method=="LassoMut"]<-nLasso

### MnLasso
Coeficientes <- coef(MnLassoMut_All, s = "lambda.min")
CoefM <- Matrix(0, nrow = nrow(Coeficientes[[1]]), ncol = length(Coeficientes))

for (n in 1:length(Coeficientes)) {
  CoefM[,n] <- Coeficientes[[n]]
  
}
colnames(CoefM) <- names(Coeficientes)
rownames(CoefM) <- rownames(Coeficientes[[1]])
image(CoefM)
nMultinomial<-sum(rowSums(abs(CoefM)) > 1e-4)
Variable_plot$Number_of_variables[Variable_plot$Method=="Multinomial"]<-nMultinomial


### BOSO
image(t(Matrix(BOSOMut_All)))
nBoso<-sum(rowSums(abs(CoefM)) > 1e-4)-1

Variable_plot$Number_of_variables[Variable_plot$Method=="BOSOMut"]<-nBoso



## MOM 
nMOM<-3
Variable_plot$Number_of_variables[Variable_plot$Method=="MOM"]<-nMOM

## Tree
nTree<-5
Variable_plot$Number_of_variables[Variable_plot$Method=="TreeMut"]<-nTree

## TreeSQRT
nTreeSqrt<-4
Variable_plot$Number_of_variables[Variable_plot$Method=="TreeMutSqrt"]<-nTreeSqrt

##KRL
nKRL<-69 #14 more relevant
Variable_plot$Number_of_variables[Variable_plot$Method=="KRL"]<-nKRL


##### Plotting-----------------------------------------------------------------
mypal = pal_npg("nrc")(9)
mypal<-mypal[2:8]

Variable_plot<-Variable_plot[order(Variable_plot$Number_of_variables, decreasing = T),]
Variable_plot$Method<-factor(Variable_plot$Method, levels=Variable_plot$Method)

pallete<-c("#DC0000FF", "#F39B7FFF",  "#91D1C2FF","#00A087FF", "#8491B4FF", "#3C5488FF", "#4DBBD5FF" ,"#4DBBD5FF")


gg_variable_number<-ggplot(Variable_plot, aes(x=Method, y=Number_of_variables, fill=Method))+
  geom_bar(stat = "identity")+scale_fill_manual(values=pallete)+theme_bw(base_size = 18)+
  ggtitle("Variable Number")+
  ylab("Number of required variables")+coord_flip()

gg_variable_number


ggarrange(gg_variable_number, gg_lolliplot, nrow=2, common.legend = T)

















