library(readr)
library(ggplot2)
library(ggbreak)
library(patchwork)
library(glmnet)
################################################################################
# Time
################################################################################
source("Code_Analysis/XAIfunctions_NV_new.R")
source("ODT/ODT functions.R")
folder_dir <- "C:/Users/ksada/OneDrive - Tecnun/Paper XAI Methods/Rresults/"

windowsFonts("Roboto" = windowsFont("Roboto"))


Table_Results <- read_delim("./data/output/Table_Results2.csv",
                            delim = ";", escape_double = FALSE, trim_ws = TRUE
)

colnames(Table_Results) <- c("Method", "Data_Type", "Computer_Time_sec")
Table_Results$Method_plot <- paste(Table_Results$Method, Table_Results$Data_Type, sep = "_")


Table_Results <- Table_Results[-which(Table_Results$Data_Type == "Both"), ]
Table_Results$Method<-factor(Table_Results$Method, levels=c("ORACLE", "MOM","KRL", "BOSO", "ODT", "ODTSqrt", "Lasso", "MnLasso"))

# KRL

# plot_computing_time<-ggplot(Table_Results, aes(x=Method, y=Computer_Time_sec, fill=Data_Type) )+
#   geom_bar(stat = "identity")+scale_fill_manual(values=pallete)+theme_bw(base_size = 18)+
#   ggtitle("Computing Time (s)")+facet_wrap(~Data_Type)+
#   ylab("Seconds(s)")
#
#
# ggplot(Table_Results, aes(x=Method_plot, y=Computer_Time_sec, fill=Method))+
#   geom_col()+scale_fill_manual(values=pallete)+theme_bw(base_size = 18)+
#   ggtitle("Computing Time (s)")+ylim(c(-1, 5100))+
#   scale_y_break(c(1, 10))+scale_y_break(c(11, 100))+scale_y_break(c(270, 3500))+
#   scale_y_break(c(4000, 4500))+scale_y_log10()+
#   ylab("Seconds(s)")+coord_flip()


pallete <- c("#367592", "#39A7AE", "#96D6B6", "#FDE5B0", "#F3908B", "#E36192", "#8E4884","#A83333")

plot_computing <- ggplot(Table_Results, aes(y = Method, x = Computer_Time_sec, fill = Method)) +
  geom_col() +
  scale_fill_manual(values = pallete) +
  # theme_bw(base_size = 25) +
  ggtitle("") +
  ylab("") +
  xlim(c(0, max(na.omit(Table_Results$Computer_Time_sec)) + 10)) +
  xlab("") +
  facet_wrap(~Data_Type,
             scales = "free_x",
             strip.position = "right",
             nrow = length(unique(Table_Results$Data_Type))
  ) +
  scale_x_cut(breaks = c(0.1, 18), which = c(1, 2, 3), scales = c(1, 2, 2)) +
  theme(
    axis.title = element_text(size = 35),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 30),
    text = element_text(size = 26, family = "Roboto"), # Increase font size and make drug names bold
    plot.background = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent")
  )
plot_computing

ggsave(paste(folder_dir, "Rimages/training_time.png", sep = ""), plot_computing, width = 25, height = 8, dpi = 1500)


################################################################################
# Variable Number
################################################################################
# Load models
load(paste(folder_dir, "Rdata/mutations_models_all.RData", sep = ""))
load(paste(folder_dir, "Rdata/MOM_model_all.RData", sep = "")) # model


Variable_plot <- data.frame(Method = unique(Table_Results$Method), Number_of_variables = NA)

### Lasso
Coeficientes <- coef(LassoMut_All, s = "lambda.min")
CoefM <- Matrix(0, nrow = nrow(Coeficientes[[1]]), ncol = length(Coeficientes))

for (n in 1:length(Coeficientes)) {
  CoefM[, n] <- Coeficientes[[n]]
}
colnames(CoefM) <- names(Coeficientes)
rownames(CoefM) <- rownames(Coeficientes[[1]])
image(CoefM)
nLasso <- sum(rowSums(abs(CoefM)) > 1e-4) - 1
Variable_plot$Number_of_variables[Variable_plot$Method == "Lasso"] <- nLasso
nLasso
### MnLasso
Coeficientes <- coef(MnLassoMut_All, s = "lambda.min")
CoefM <- Matrix(0, nrow = nrow(Coeficientes[[1]]), ncol = length(Coeficientes))

for (n in 1:length(Coeficientes)) {
  CoefM[, n] <- Coeficientes[[n]]
}
colnames(CoefM) <- names(Coeficientes)
rownames(CoefM) <- rownames(Coeficientes[[1]])
image(CoefM)
nMultinomial <- sum(rowSums(abs(CoefM)) > 1e-4)
Variable_plot$Number_of_variables[Variable_plot$Method == "MnLasso"] <- nMultinomial
nMultinomial

### BOSO
image(t(Matrix(BOSOMut_All)))
CoefM <- t(Matrix(BOSOMut_All))
nBoso <- sum(rowSums(abs(CoefM)) > 1e-4) - 1

Variable_plot$Number_of_variables[Variable_plot$Method == "BOSO"] <- nBoso

## MOM
nMOM <- length(MOM_All$TreatmentGuideline)
Variable_plot$Number_of_variables[Variable_plot$Method == "MOM"] <- nMOM

## Tree
nTree <- length(nodeids(ODTMut_All)) - length(nodeids(ODTMut_All, terminal = T))
plot(niceTree(ODTMut_All,folder))
nTree
Variable_plot$Number_of_variables[Variable_plot$Method == "ODT"] <- nTree

## TreeSQRT
nTreeSqrt <- length(nodeids(ODTSqrtMut_All)) - length(nodeids(ODTSqrtMut_All, terminal = T))
plot(niceTree(ODTSqrtMut_All))
nTreeSqrt
folder <- paste(folder_dir,"Rimages",sep="")
niceTree(ODTSqrtMut_All,folder)
Variable_plot$Number_of_variables[Variable_plot$Method == "ODTSqrt"] <- nTreeSqrt


plot(niceTree(ODTSqrt,folder))

## KRL
nKRL <- 70# 14 more relevant
Variable_plot$Number_of_variables[Variable_plot$Method == "KRL"] <- nKRL

Variable_plot <- Variable_plot[order(Variable_plot$Method), ]
##### Plotting-----------------------------------------------------------------

#Variable_plot <- Variable_plot[order(Variable_plot$Number_of_variables, decreasing = T), ]
Variable_plot$Method<-factor(Variable_plot$Method, levels=c("ORACLE", "MOM","KRL", "BOSO", "ODT", "ODTSqrt", "Lasso", "MnLasso"))

gg_variable_number <- ggplot(Variable_plot, aes(x = Method, y = Number_of_variables, fill = Method)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pallete) +
  ggtitle("") +
  ylab("") +
  xlab("") +
  coord_flip() +
  theme(
    axis.title = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    text = element_text(size = 26, family = "Roboto"), # Increase font size and make drug names bold
    #panel.background = element_rect(fill = "transparent"), # Set background to transparent
    plot.background = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent"),
    legend.position = "none"
  )

gg_variable_number

ggsave(paste(folder_dir, "Rimages/variable_number.png", sep = ""), gg_variable_number, width = 18, height = 4, dpi = 1500)


# MOM PLOT
GL <- MOM_All$TreatmentGuideline
for (s in 1:nrow(GL)) {
  selectedDrug <- unlist(strsplit(GL[s,2], split="_"))[1]
  mut_step <- unlist(strsplit(GL[s,2], split="_"))[2]
  status <- unlist(strsplit(GL[s,2], split="_"))[3]
  GL$Drug[s] <- selectedDrug
  GL$Gene[s] <- mut_step
  GL$Status[s] <- status
}
library(Rgraphviz)


swr = function(string, nwrap=12) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

# Create line breaks in DrugName
GL$Drug <- swr(GL$Drug)


node1<-GL$Gene[1]
node2<-GL$Gene[2]
node3<-GL$Drug[1]
node4<-GL$Gene[3]
node5<-GL$Drug[2]
node6<-GL$Drug[3]
node7<-GL$Drug[4]

nodeNames<-c(node1,node2,node3,node4, node5,node6, node7)

rEG <- new("graphNEL", nodes=nodeNames, edgemode="directed")

# Draw the "lines" or "branches" of the probability Tree
rEG <- addEdge(nodeNames[1], nodeNames[2], rEG, 1)
rEG <- addEdge(nodeNames[1], nodeNames[3], rEG, 1)
rEG <- addEdge(nodeNames[2], nodeNames[4], rEG, 1)
rEG <- addEdge(nodeNames[2], nodeNames[5], rEG, 1)
rEG <- addEdge(nodeNames[4], nodeNames[6], rEG, 1)
rEG <- addEdge(nodeNames[4], nodeNames[7], rEG, 1)

attributes<-list(node=list(label="foo", fillcolor="#F3908B", fontsize="16", font="Roboto",lwd=3),
                 edge=list(color="black",arrowsize=0.5,lwd=3),graph=list(rankdir="TB",font="Roboto"))

plot(rEG, attrs=attributes)

n <- network(rgraph(10, tprob = 0.2), directed = TRUE)
