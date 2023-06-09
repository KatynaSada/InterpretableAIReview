# Loading Libraries and Packages --------------------------------------------
#The first part is to load all the libraries and external functions that are required for the analysis. # nolint # nolint

library(readxl)
library(RColorBrewer)
library(matrixStats)
library(partykit)
library(glmnet)
library(BOSO)

source("Code_Analysis/XAIfunctions_NV_new.R")
source("./MOM/MOM_source_Functions.R")

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

# Independent Cohort Validation --------------------------------------------
"
One of the main challenges of Machine Learning, including Precision Medicine, is generalization, i.e. the ability to adapt to new, previously unseen data.
All the methods were tested on the GDSC AML dataset to check their generalization ability.
The models were trained using the BeatAML dataset and were used to predict the optimal drug for AML cell lines from GDSC using its mutation files.
Each of the cell lines was recommended a drug, we compared the all-samples IC50 for all the models and against the Oracle (the drug with the minimum IC50 for each cell line).
"

# Load models and training data
load(paste(folder_dir,"/mutations_models_all.RData",sep=""))
load(paste(folder_dir,"/training_data.RData",sep=""))

load(paste(folder_dir,"/MOM_model_all.RData",sep="")) # model

# KRL
KRL_results_cv <- read_csv(paste(folder_dir,"/KRL_results_5cv_BeatAML2.csv",sep=""))
KRL_results_cv$Drug<-KRL_results_cv$Drug+1 
KRL_results_cv$IC50<-drug_response_w12[cbind(KRL_results_cv$Patients, KRL_results_cv$DrugName)]

### Validation in GDSC
"
The models were trained using the BeatAML dataset and were used to predict the optimal drug for AML cell lines from GDSC using its mutation files. 
Each of the cell lines was recommended a drug, we compared the all-samples IC50 for all the models and against the Oracle (the drug with the minimum IC50 for each cell line).
The function `predictMOM_ALL`.
Applies the decision tree generated by MOM in AML.

Expresion medida en GSDC puede que no se compare con la expresion de BeatAML, microArrays 
"
#### Load GDSC and CCLE data
GDSC1_ESS<-read.csv("./data/input/GDSC/LAML_IC_Wed Nov 24 10_31_42 2021.csv")
GDSC1_MUT<-read.csv("./data/input/GDSC/LAML_Genetic_feature_variant_Wed Nov 24 GDSC1.csv")

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
                  ncol=ncol(mutations_w12))

colnames(BM_matrix)<-colnames(mutations_w12)
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

BM_matrix<-BM_matrix[,match( colnames(BM_matrix),colnames(mutations_w12))]
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

drugsGDSC<-sapply(colnames(drug_response_w12), FUN=function(X){return(unlist(strsplit(X, split=" "))[1])})

GDSC_AML<-GDSC_AML[,na.omit(match(drugsGDSC, colnames(GDSC_AML)))]

colnamesGDSC<-colnames(drug_response_w12)[match(colnames(GDSC_AML),drugsGDSC)]

colnames(GDSC_AML)<-colnamesGDSC

GDSC_AML<-GDSC_AML[which(rownames(GDSC_AML) %in% rownames(mutations_AML_GDSC)),]
GDSC_AML<-GDSC_AML[match(rownames(mutations_AML_GDSC), rownames(GDSC_AML)),]

identical(rownames(GDSC_AML), rownames(mutations_AML_GDSC))

# FOR KRL
write.csv(mutations_AML_GDSC, file="data/input/GDSC/mutations_AML_GDSC_NAs.csv")

# Make Predictions with the models

# Predict with GDSC data
# We used the models trained in BeatAML to predict over GDSC. 
identical(rownames(GDSC_AML), rownames(mutations_AML_GDSC)) 
identical(colnames(mutations_AML_GDSC), colnames(mutations_w12)) #TRUE

Guideline_ODTSqrtMut<-predictTreeMut(ODTSqrtMut_All,drug_response_w12, mutations_AML_GDSC)
Guideline_ODTSqrtMut<-colnames(drug_response_w12)[Guideline_ODTSqrtMut]

Guideline_BOSO<-predictBOSO(BOSOMut_All,drug_response_w12,mutations_AML_GDSC)
Guideline_BOSO<-colnames(drug_response_w12)[Guideline_BOSO]

Guideline_ODTMut<-predictTreeMut(ODTMut_All,drug_response_w12, mutations_AML_GDSC)
Guideline_ODTMut<-colnames(drug_response_w12)[Guideline_ODTMut]

Guideline_Multinomial<-predictMnLasso(MnLassoMut_All,drug_response_w12, mutations_AML_GDSC)
Guideline_Multinomial<-colnames(drug_response_w12)[Guideline_Multinomial]

Guideline_Lasso<-predictLasso(LassoMut_All,drug_response_w12, mutations_AML_GDSC)
Guideline_Lasso<-colnames(drug_response_w12)[Guideline_Lasso]

Guideline_MOM<-CV_Prediction2(MILP_classification = MOM_All, 
                              kfold_validation_drug = NULL,
                              kfold_validation_mut = mutations_AML_GDSC )

Guideline_KRL<-read_csv("./data/output/KRL_results_GDSC.csv")
Guideline_KRL$DrugName[which(!(Guideline_KRL$DrugName %in% colnames(drug_response_w12)))]<-NA
Guideline_KRL<-Guideline_KRL$DrugName

#### Obtain Prediction IC50 values

#We will now obtain the GDSC IC50 values corresponding to the predicted treatments

# Validate the Guidelines..................................................

identical(rownames(GDSC_AML), rownames(mutations_AML_GDSC)) # TRUE

Guideline_ORACLE <- NULL

GDSC_IC50<-data.frame(CL=rownames(mutations_AML_GDSC),
                      ODTMut=NA,
                      ODTSqrtMut=NA,
                      LassoMut=NA,
                      MnLassoMut=NA,
                      BOSOMut=NA,
                      MOM=NA,
                      KRL=NA,
                      ORACLE=NA)
# ENTENDER BIEN
for(i in 1:nrow(GDSC_IC50)){
  CL<-GDSC_IC50$CL[i]
  
  # ODTMut
  dODTMut<-Guideline_ODTMut[i]
  if(dODTMut %in% colnames(GDSC_AML)){
    GDSC_IC50$ODTMut[i]<-GDSC_AML[CL,dODTMut]
  }
  
  # ODTSqrtMut
  dODTSqrtMut<-Guideline_ODTSqrtMut[i]
  if(dODTSqrtMut %in% colnames(GDSC_AML)){
    GDSC_IC50$ODTSqrtMut[i]<-GDSC_AML[CL,dODTSqrtMut]
  }
  
  # Lasso
  dLasso<-Guideline_Lasso[i]
  if(dLasso %in% colnames(GDSC_AML)){
    GDSC_IC50$LassoMut[i]<-GDSC_AML[CL,dLasso]
  } 
  
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
  
  #MOM
  dMOM<-Guideline_MOM[i]
  if(dMOM %in% colnames(GDSC_AML)){
    GDSC_IC50$MOM[i]<-GDSC_AML[CL,dMOM]
  }
  
  #KRL
  dKRL<-Guideline_KRL[i]
  if(dKRL %in% colnames(GDSC_AML)){
    GDSC_IC50$KRL[i]<-GDSC_AML[CL,dKRL]
  }
  
  #ORACLE
  Guideline_ORACLE[i] <- names(which(GDSC_AML[CL,] == min(GDSC_AML[CL,])))
  GDSC_IC50$ORACLE[i]<-min(GDSC_AML[CL,])
}

rownames(GDSC_IC50)<-GDSC_IC50$CL

# Plotting GDSC Validation 
number_methods <- ncol(GDSC_IC50)-1
Validation_plot<-data.frame(
  Dataset=rep("GDSC", nrow(GDSC_IC50)*(number_methods)), 
  CellLine=rep(rownames(GDSC_IC50),number_methods),
  Method=rep(colnames(GDSC_IC50[,-1]), each=nrow(GDSC_IC50)), 
  IC50=unlist(GDSC_IC50[,-1]))

rownames(Validation_plot)<-1:nrow(Validation_plot)

mycomparisons<-combn(unique(Validation_plot$Method),2)
mycomparisons<-as.list(data.frame(mycomparisons))
mycomparisons<-mycomparisons[-grep("ORACLE",mycomparisons)]

Plotting<-Validation_plot[Validation_plot$Dataset=="GDSC",]

library(ggsci)

Plotting$Method<-factor(Plotting$Method, levels=c("ORACLE", "MOM","KRL", "BOSOMut", "ODTMut", "ODTSqrtMut", "LassoMut", "MnLassoMut"))

treatment_plot_violin <- ggplot(Plotting, aes(x = Method, y = IC50, fill = Method)) +
  geom_violin(lwd = 0.8, position = position_dodge(0.8)) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ylab("GDSC IC50") +
  ggtitle("Test in GDSC (mutations)") +
  theme(text = element_text(family = "Roboto"),
        axis.text.x = element_text(size = 18, angle=60,vjust=0.65),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 21),
        plot.title = element_text(size = 26, face = "bold"),
        plot.background = element_rect(fill = "transparent", size = 0, color = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        )
treatment_plot_violin
ggsave(paste(folder_dir,"/images/violin_mutations_GDSC.png",sep=""), treatment_plot_violin, width = 8, height = 7, dpi = 1000)

colors <- c("#A83333","#367592", "#39A7AE", "#96D6B6", "#FDE5B0", "#F3908B", "#E36192", "#8E4884")

g2_gdsc <- ggplot(Plotting, aes(x = Method, y = IC50, fill = Method)) +
  geom_boxplot(lwd = 0.8, position = position_dodge(0.8)) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ylab("GDSC IC50") +
  ggtitle("Test in GDSC (mutations)") +
  theme(text = element_text(family = "Roboto"),
        axis.text.x = element_text(size = 24, angle=60,vjust=0.65),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 27),
        plot.title = element_text(size = 32, face = "bold"),
        plot.background = element_rect(fill = "transparent", size = 0, color = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        legend.position = "NULL")
g2_gdsc
ggsave(paste(folder_dir,"/images/boxplot_mutations_GDSC.png",sep=""), g2_gdsc, width = 8, height = 7, dpi = 1000)

lapply(unique(Plotting$Method), function(method){
  return(c(as.character(method),median(na.omit(Plotting$IC50[Plotting$Method==method]))))
})

wilcox.test(Plotting$IC50[Plotting$Method=="ODTSqrtMut"], 
            Plotting$IC50[Plotting$Method=="ODTMut"], alternative="less")


# Intragroup Validation in GDSC
Validation_GDSC<-data.frame(patients=rep(rownames(GDSC_AML),number_methods),
                            drug=c(
                              Guideline_ODTMut,
                              Guideline_ODTSqrtMut,
                              Guideline_Lasso,
                              Guideline_Multinomial,
                              Guideline_BOSO,
                              Guideline_MOM,
                              Guideline_KRL,
                              Guideline_ORACLE
                            ), 
                            IC50=c(
                              GDSC_IC50$ODTMut,
                              GDSC_IC50$ODTSqrtMut,
                              GDSC_IC50$LassoMut,
                              GDSC_IC50$MnLassoMut,
                              GDSC_IC50$BOSOMut,
                              GDSC_IC50$MOM,
                              GDSC_IC50$KRL,
                              GDSC_IC50$ORACLE
                            ),
                            Method=rep(unique(Validation_plot$Method), each=nrow(GDSC_AML)))

colors2 <-  c("#FDE5B0", "#F3908B", "#E36192", "#8E4884", "#96D6B6", "#367592","#39A7AE")

methods <- unique(Validation_GDSC$Method)[-8]

n <- 5
gg_drug_plots <- plotIntragroup(methods[n],colors2[n],Validation_GDSC,GDSC_AML,"GDSC")
gg_drug_plots

for (n in 1:length(methods)) {
  gg_drug_plots <- plotIntragroup(methods[n],colors2[n],Validation_GDSC,GDSC_AML,"GDSC")
}

chosen_method <- methods[n]
color <- colors2[n]
Validation_Beat <- Validation_GDSC
drug_response <- GDSC_AML
dataset<- "GDSC"
  