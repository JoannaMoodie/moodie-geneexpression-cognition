# Script written by Joanna Moodie, July 2023. The gene expression components are calculated and validated, and correlations between gene expressions and g-morphometry profiles are calculated. 
library(ggplot2)
library(caret)
library(tidyverse)
library(psych)
library(tidyverse)
library(readr)
library(readxl)
library(ggplot2)
library(factoextra)
library(ggseg3d)
library(matrixStats)
library(reshape2)
library(GPArotation)
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")
library(cowplot)
library(ggpubr)
library(FSA)
library(lavaan)
library(lm.beta)


## load data
dataAllen8000 <- read.table("AllenHBA_DK_ExpressionMatrix.tsv", sep = "\t", header = T) # from French and Paus (2015)
tempdata <- dataAllen8000
tempdata[,2] = round(tempdata[,2], digits = 3) # subset to the between-donor consistent genes as defined in French and Paus (2015)
tempdata <- subset(tempdata, tempdata[,2] >= 0.446)
dataAllen8000 <- t(tempdata[,3:70])

# set up colnames and rownames
colnames(dataAllen8000) <- tempdata$X
dataAllen8000 <- dataAllen8000[order(rownames(dataAllen8000)),]
rownames(dataAllen8000) <- gsub("ctx-lh-", "", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("ctx-rh-", "", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("banks", "Bank s", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("caudalanteriorcin", "Caudal anterior cin", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("caudalmiddlefro", "Caudal middle fro", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("cuneus", "Cuneus", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("entorh", "Entorh", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("frontalpol", "Frontal pol", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("fusi", "Fusi", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("inferiorp", "Inferior p", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("inferiort", "Inferior t", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("insul", "Insul", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("isthmuscin", "Isthmus cin", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("lateralo", "Lateral o", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("lingual", "Lingual", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("medialorbitofro", "Medial orbito fro", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("middletem", "Middle tem", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("para", "Para", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("parso", "Pars o", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("parst", "Pars t", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("perical", "Perical", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("postcen", "Postcen", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("posteriorci", "Posterior ci", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("precen", "Precen", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("preCun", "Precun", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("rostralanteriorcing", "Rostral anterior cing", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("rostralmiddlefron", "Rostral middle fron", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("superiorfro", "Superior fro", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("superiorp", "Superior p", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("superiortemp", "Superior temp", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("suprama", "Suprama", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("temporalpo", "Temporal po", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("transversetem", "Transverse tem", rownames(dataAllen8000))
rownames(dataAllen8000) <- gsub("PreCun", "Precun", rownames(dataAllen8000))

# run PCA
pcaAllen <- prcomp(dataAllen8000, scale = T, center =T)
fviz_eig(pcaAllen, geom = c("line")) + geom_line(size = 1.5) + geom_point(size = 3) + theme_cowplot() + theme(text = element_text(size = 20, face = "bold"), axis.text.x = element_text(size = 20, face = "bold"), axis.text.y = element_text(size = 20, face = "bold")) + ylab("Variance explained (%)\n") + xlab("\nPrincipal component")
rawLoadings <- pcaAllen$rotation[,1:2] %*% diag(pcaAllen$sdev, 2, 2)
rotatedLoadings <- varimax(rawLoadings)$loadings
invLoadings     <- t(pracma::pinv(rotatedLoadings))
scores          <- scale(dataAllen8000) %*% invLoadings
scores <- as.data.frame(scores)
scores$hemisphere <- c(rep("lh", 34), rep("rh", 34))
scores$region <- rownames(scores)
scores$region <- gsub("ctx.lh.", "", scores$region)
scores$region <- gsub("ctx.rh.", "", scores$region)
colnames(scores)[1:2] <- c("unrotatedPC1", "unrotatedPC2")

write.table(rotatedLoadings, 'rotatedLoadings.csv', row.names = F)

cor(scores[,1:2])
cor(scores[1:34,1:2])
cor(scores[35:68,1:2])
scaledscores <- data.frame(region = scores$region, scaledscoresPC1 = c(scale(scores[1:34,1]), scale(scores[35:68,1])), scaledscoresPC2 = c(scale(scores[1:34,2]), scale(scores[35:68,2]))) # scale scores for left and right hemispheres

scaledscores$hemisphere <- NA
scaledscores$hemisphere[which(substr(rownames(scores), 1, 5) == "ctx.l")] <- "lh"
scaledscores$hemisphere[which(substr(rownames(scores), 1, 5) == "ctx.r")] <- "rh"

cor.test(scaledscores[1:34,2], scaledscores[35:68,2])
cor.test(scaledscores[1:34,3], scaledscores[35:68,3])
cor.test(scaledscores[,2], scaledscores[,3])

write.table(scaledscores, 'rotated_scaled_genexpression_PCscores.csv', sep = ",", row.names = F)
write.table(scores, 'rotated_unscaled_genexpression_PCscores.csv', sep = ",", row.names = F)

g_regional_cognition_metaout <- read.csv('data_metanalysis_output_estimates.csv')
g_regional_cognition_metaout$hemisphere <- NA
g_regional_cognition_metaout$hemisphere[which(substr(g_regional_cognition_metaout$region, 1, 3) == "lh_")] = "lh"
g_regional_cognition_metaout$hemisphere[which(substr(g_regional_cognition_metaout$region, 1, 3) == "rh_")] = "rh"
g_regional_cognition_metaout$region <- gsub("lh_vol_", "", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("rh_vol_", "", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("banks", "Bank s", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("caudalanteriorcin", "Caudal anterior cin", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("caudalmiddlefro", "Caudal middle fro", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("cuneus", "Cuneus", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("entorh", "Entorh", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("frontalpol", "Frontal pol", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("fusi", "Fusi", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("inferiorp", "Inferior p", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("inferiort", "Inferior t", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("insul", "Insul", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("isthmus", "Isthmus cingulate", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("lateralo", "Lateral o", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("lingual", "Lingual", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("medialorbitofro", "Medial orbito fro", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("middletem", "Middle tem", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("para", "Para", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("parso", "Pars o", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("parst", "Pars t", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("perical", "Perical", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("postcen", "Postcen", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("posteriorci", "Posterior ci", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("precen", "Precen", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("preCun", "Precun", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("rostralanteriorcing", "Rostral anterior cing", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("rostralmiddlefron", "Rostral middle fron", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("superiorfro", "Superior fro", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("superiorp", "Superior p", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("superiortemp", "Superior temp", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("suprama", "Suprama", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("temporalpo", "Temporal po", g_regional_cognition_metaout$region)
g_regional_cognition_metaout$region <- gsub("transversetem", "Transverse tem", g_regional_cognition_metaout$region)


testingcorrelations <- merge(scaledscores, g_regional_cognition_metaout, by = c("region", "hemisphere"))
testingcorrelations$absscaledscoresPC1 <- abs(testingcorrelations$scaledscoresPC1)
testingcorrelations$absscaledscoresPC2 <- abs(testingcorrelations$scaledscoresPC2)
testingcorrelations <- merge(testingcorrelations, scores, by = c("region", "hemisphere"))
write.table(testingcorrelations, 'testingcorrelations.csv', row.names =F)

# internal cross-validation -----------------------------------------------------
data <-  dataAllen8000
data<-as.data.frame(data)
data$index <- c(1:68)
a <- matrix(0, ncol = 250, nrow = 10) 
for (i in 1:250) {
  set.seed(1000)
  index <- createMultiFolds(data$index, k = 5, times = 50)[i]
  
  train_df <- data[unlist(index),]
  test_df <- data[-unlist(index),]
  dim(train_df)
  dim(test_df)
  
  train_pca1 <- stats::prcomp(train_df, center = T, scale = T)
  test_pca1 <- stats::prcomp(test_df, center = T, scale = T)
  
  a[,i] <- abs(diag(psych::factor.congruence(train_pca1$rotation,test_pca1$rotation)))[1:10]
  
}
boxplot(t(a), xlab = "PCs", ylab = "Absolute factor congruence coefficient")
plotk <- a
plotk <- as.data.frame(plotk)
plotk$PCNu <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
plotk <- melt(plotk, id = "PCNu") 
plotk$PCNu <- factor(plotk$PCNu, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
ggplot(data = plotk, aes(x = PCNu, y=value)) + xlab("Unrotated component") + ylab("\nAbsolute factor congruence coefficient") + geom_flat_violin(fill = "steel blue", position = position_nudge(x = .2, y = 0)) +geom_boxplot(width = .1, color = "BLACK", outlier.size = .5) +coord_flip() + theme_cowplot() + theme( text = element_text(size=16)) + scale_x_discrete(limits = rev(levels(plotk$PCNu))) + labs(subtitle ="68 regions, 5 folds, 50 times")
ggsave(filename = "kfold_Allenall.jpeg", bg = "white", width = 7, height = 7)

# external validation -----------------------------------------------------
# Brainspan external validation
loadBrainspan <- read_csv("/home/jmoodie/Documents/Expression/Protein/data/BSnew.csv")
STC <- loadBrainspan %>% dplyr:: select(starts_with("STC"))
dataBrainspanall <- as.data.frame(apply(STC, 1, median, na.rm = T))
dataBrainspanall$STC <- as.data.frame(apply(STC, 1, median, na.rm = T))
dataBrainspanall<- as.data.frame(dataBrainspanall[,1])
MFC <- loadBrainspan %>% dplyr:: select(starts_with("MFC"))
dataBrainspanall$MFC = as.data.frame(apply(MFC, 1, median, na.rm = T))
DFC <- loadBrainspan %>% dplyr:: select(starts_with("DFC"))
dataBrainspanall$DFC = as.data.frame(apply(DFC, 1, median, na.rm = T))
OFC <- loadBrainspan %>% dplyr:: select(starts_with("OFC"))
dataBrainspanall$OFC = as.data.frame(apply(OFC, 1, median, na.rm = T))
ITC <- loadBrainspan %>% dplyr:: select(starts_with("ITC"))
dataBrainspanall$ITC = as.data.frame(apply(ITC, 1, median, na.rm = T))
VFC <- loadBrainspan %>% dplyr:: select(starts_with("VFC"))
dataBrainspanall$VFC = as.data.frame(apply(VFC, 1, median, na.rm = T))
A1C <- loadBrainspan %>% dplyr:: select(starts_with("A1C"))
dataBrainspanall$A1C = as.data.frame(apply(A1C, 1, median, na.rm = T))
V1C <- loadBrainspan %>% dplyr:: select(starts_with("V1C"))
dataBrainspanall$V1C = as.data.frame(apply(V1C, 1, median, na.rm = T))
M1C <- loadBrainspan %>% dplyr:: select(starts_with("M1C"))
dataBrainspanall$M1C = as.data.frame(apply(M1C, 1, median, na.rm = T))
IPC <- loadBrainspan %>% dplyr:: select(starts_with("IPC"))
dataBrainspanall$IPC = as.data.frame(apply(IPC, 1, median, na.rm = T))
S1C <- loadBrainspan %>% dplyr:: select(starts_with("S1C"))
dataBrainspanall$S1C = as.data.frame(apply(S1C, 1, median, na.rm = T))
colnames(dataBrainspanall) <- cbind("STC", "MFC", "DFC", "OFC", "ITC", "VFC", "A1C", "V1C", "M1C", "IPC", "S1C")
allBS <- cbind(STC, MFC, DFC, OFC, ITC, VFC, A1C, V1C, M1C, IPC, S1C)
colnames(allBS)
test <- allBS %>% dplyr:: select(ends_with("C"))
dim(test)
a <- diag(cor(t(test), t(dataBrainspanall), method = "spearman"))
test <- allBS %>% dplyr:: select(contains("_1"))
dim(test)
b <- diag(cor(t(test), t(dataBrainspanall), method = "spearman"))
test <- allBS %>% dplyr:: select(contains("_2"))
dim(test)
c <- diag(cor(t(test), t(dataBrainspanall), method = "spearman"))
test <- allBS %>% dplyr:: select(contains("_3"))
dim(test)
d <- diag(cor(t(test), t(dataBrainspanall), method = "spearman"))
test <- allBS %>% dplyr:: select(contains("_4"))
dim(test)
e <- diag(cor(t(test), t(dataBrainspanall), method = "spearman"))
result <- cbind(a, b, c, d, e)
result <- rowMeans(result)
length(which(result > .446)) 
dataBrainspanall <- subset(dataBrainspanall, round(result, digits = 3) >= 0.446)
dim(dataBrainspanall)
dataBrainspanall <- t(dataBrainspanall)
colnames(dataBrainspanall) <- subset(loadBrainspan, round(result, digits = 3) >= 0.446)$Genes
dataBrainspanall <- dataBrainspanall[sort(row.names(dataBrainspanall)), ]

a <- dataBrainspanall[sort(row.names(dataBrainspanall)), ]
rownames(a)
names <- c("A1C", "DFC", "IPC", "ITC", "M1C", "MFC", "OFC", "S1C", "STC", "V1C", "VFC")
a<-melt(dataBrainspanall, measure.vars = rownames(dataBrainspanall)) 
colnames(a) <- c("region", "Gene", "value")
ggplot(data = a, aes(x = region, y = value, group = Gene)) +geom_line(aes(color=Gene), alpha = 0.03) + theme_classic() + theme(legend.position = "none", text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ xlab("\nRegion") + ylab("Expression Value") + scale_x_discrete(labels =names) +
  labs(subtitle="BrainSpan")+ expand_limits(y = c(3, 15)) +scale_y_continuous(breaks = c(3,6,9, 12 ,15))
ggsave("expressionplot_dataBrainspan.jpg")

pcaBrainspan <- stats::prcomp(dataBrainspanall, center=T, scale = T)
BSrawLoadings <- pcaBrainspan$rotation[,1:2] %*% diag(pcaBrainspan$sdev, 2, 2)
BSrotatedLoadings <- varimax(BSrawLoadings)$loadings
BSinvLoadings     <- t(pracma::pinv(BSrotatedLoadings))
BSscores          <- scale(dataBrainspanall) %*% BSinvLoadings

a <- pcaBrainspan$rotation %*% diag(pcaBrainspan$sdev,11,11)
b <- pcaAllen$rotation %*% diag(pcaAllen$sdev,68,68)
c <- merge(a, b, by = 0) 
dim(c) 
rownames(c) <-c[,1]
c[,1] <- NULL
d <- c[,1:11]
e <- c[,12:79]
f <- (psych::factor.congruence(d,e))
colnames(f) <- gsub("\\.y", "", colnames(f))
rownames(f) <- gsub("\\.x", "", rownames(f))
f<-f[1:5,1:5]
factorcong <- c(0.88, 0.24)
factorcong <- as.data.frame(t(factorcong))
colnames(factorcong) <- c("Component 1", "Component 2")
factorcong$validationset <- "BrainSpan"

# Kang et al. (2011) external validation
kangdata <-read_csv("/home/jmoodie/Documents/Expression/Protein/data/kangMatrix.csv")
which(a<-!(duplicated(kangdata) | duplicated(kangdata, fromLast = TRUE)) == FALSE)
OFC <- kangdata %>% dplyr:: select(contains("OFC.L"))
datakangleft <- as.data.frame(apply(OFC,1,median, nam = T))
datakangleft$OFC <-  as.data.frame(apply(OFC,1,median, nam = T))
datakangleft <- as.data.frame(datakangleft[,1])
DFC <- kangdata %>% dplyr:: select(contains("DFC.L"))
datakangleft$DFC <- as.data.frame(apply(DFC,1, median, nam = T))
VFC <- kangdata %>% dplyr:: select(contains("VFC.L"))
datakangleft$VFC <- as.data.frame(apply(VFC,1, median, nam = T))
MFC <- kangdata %>% dplyr:: select(contains("MFC.L"))
datakangleft$MFC <- as.data.frame(apply(MFC,1, median, nam = T))
M1C <- kangdata %>% dplyr:: select(contains("M1C.L"))
datakangleft$M1C <- as.data.frame(apply(M1C,1, median, nam = T))
S1C <- kangdata %>% dplyr:: select(contains("S1C.L"))
datakangleft$S1C <- as.data.frame(apply(S1C,1, median, nam = T))
IPC <- kangdata %>% dplyr:: select(contains("IPC.L"))
datakangleft$IPC <- as.data.frame(apply(IPC,1, median, nam = T))
A1C <- kangdata %>% dplyr:: select(contains("A1C.L"))
datakangleft$A1C <- as.data.frame(apply(A1C,1, median, nam = T))
STC <- kangdata %>% dplyr:: select(contains("STC.L"))
datakangleft$STC <- as.data.frame(apply(STC,1, median, nam = T))
ITC <- kangdata %>% dplyr:: select(contains("ITC.L"))
datakangleft$ITC <- as.data.frame(apply(ITC,1, median, nam = T))
V1C <- kangdata %>% dplyr:: select(contains("V1C.L"))
datakangleft$V1C <- as.data.frame(apply(V1C,1, median, nam = T))
colnames(datakangleft) <- c("OFC", "DFC", "VFC", "MFC", "M1C", "S1C", "IPC", "A1C", "STC", "ITC", "V1C")
allkang <- cbind(OFC, DFC, VFC, MFC, M1C, S1C, IPC, A1C, STC, ITC, V1C)
colnames(allkang)

# get consistent profiles
test <- allkang %>% dplyr:: select(contains("123"))
dim(test)
a <- diag(cor(t(test), t(datakangleft), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("125"))
dim(test)
colnames(test)
b <- diag(cor(t(test), t(datakangleft[,c(1, 4, 8)]), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("126"))
dim(test)
colnames(test)
c <- diag(cor(t(test), t(datakangleft[,-8]), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("133"))
dim(test)
colnames(test)
d <- diag(cor(t(test), t(datakangleft[,-c(3,5,6,8,9)]), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("135"))
dim(test)
e <- diag(cor(t(test), t(datakangleft), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("136"))
dim(test)
f <- diag(cor(t(test), t(datakangleft), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("144"))
dim(test)
g <- diag(cor(t(test), t(datakangleft), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("145"))
dim(test)
h <- diag(cor(t(test), t(datakangleft), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("182"))
dim(test)
colnames(test)
i <- diag(cor(t(test), t(datakangleft[,-c(5,6)]), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("183"))
dim(test)
j <- diag(cor(t(test), t(datakangleft), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("187"))
dim(test)
k <- diag(cor(t(test), t(datakangleft), method = "spearman"))
result <- cbind(a, b, c, d, e, f, g, h, i, j, k)
result <- rowMeans(result)
length(which(result > .446)) 
datakangleft <- subset(datakangleft, round(result, digits = 3) >= 0.446)
dim(datakangleft)
datakangleft <- t(datakangleft)
colnames(datakangleft) <- subset(kangdata, round(result, digits = 3) >= 0.446)$Genes
datakangleft <- datakangleft[sort(row.names(datakangleft)), ]
datakangleft <-as.data.frame(datakangleft) 
drops <- which(apply(datakangleft, 2, max) - apply(datakangleft, 2, min) > 4)
datakangleft$GABRQ <- NA
datakangleft$SLN <-NA
datakangleft <- datakangleft[,apply(datakangleft, 2, function(x) !any(is.na(x)))]
kangdata <-read_csv("/home/jmoodie/Documents/Expression/Protein/data/kangMatrix.csv")
OFC <- kangdata %>% dplyr:: select(contains("OFC.R"))
datakangright <- as.data.frame(apply(OFC,1,median, na.rm = T))
datakangright$OFC <-  as.data.frame(apply(OFC,1,median, na.rm = T))
datakangright <- as.data.frame(datakangright[,1])
DFC <- kangdata %>% dplyr:: select(contains("DFC.R"))
datakangright$DFC <- as.data.frame(apply(DFC,1, median, na.rm = T))
VFC <- kangdata %>% dplyr:: select(contains("VFC.R"))
datakangright$VFC <- as.data.frame(apply(VFC,1, median, na.rm = T))
MFC <- kangdata %>% dplyr:: select(contains("MFC.R"))
datakangright$MFC <- as.data.frame(apply(MFC,1, median, na.rm = T))
M1C <- kangdata %>% dplyr:: select(contains("M1C.R"))
datakangright$M1C <- as.data.frame(apply(M1C,1, median, na.rm = T))
S1C <- kangdata %>% dplyr:: select(contains("S1C.R"))
datakangright$S1C <- as.data.frame(apply(S1C,1, median, na.rm = T))
IPC <- kangdata %>% dplyr:: select(contains("IPC.R"))
datakangright$IPC <- as.data.frame(apply(IPC,1, median, na.rm = T))
A1C <- kangdata %>% dplyr:: select(contains("A1C.R"))
datakangright$A1C <- as.data.frame(apply(A1C,1, median, na.rm = T))
STC <- kangdata %>% dplyr:: select(contains("STC.R"))
datakangright$STC <- as.data.frame(apply(STC,1, median, na.rm = T))
ITC <- kangdata %>% dplyr:: select(contains("ITC.R"))
datakangright$ITC <- as.data.frame(apply(ITC,1, median, na.rm = T))
V1C <- kangdata %>% dplyr:: select(contains("V1C.R"))
datakangright$V1C <- as.data.frame(apply(V1C,1, median, na.rm = T))
colnames(datakangright) <- c("OFC", "DFC", "VFC", "MFC", "M1C", "S1C", "IPC", "A1C", "STC", "ITC", "V1C")
allkang <- cbind(OFC, DFC, VFC, MFC, M1C, S1C, IPC, A1C, STC, ITC, V1C)
# get consistent profiles
test <- allkang %>% dplyr:: select(contains("123"))
dim(test)
a <- diag(cor(t(test), t(datakangright), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("125"))
dim(test)
colnames(test)
b <- diag(cor(t(test), t(datakangright[,c(1, 4, 8)]), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("126"))
dim(test)
colnames(test)
c <- diag(cor(t(test), t(datakangright), method = "spearman"))
# test <- allkang %>% dplyr:: select(contains("133"))
# dim(test)
# colnames(test)
# d <- diag(cor(t(test), t(datakangright[,-c(3,5,6,8,9)]), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("135"))
dim(test)
e <- diag(cor(t(test), t(datakangright), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("136"))
dim(test)
f <- diag(cor(t(test), t(datakangright), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("144"))
dim(test)
g <- diag(cor(t(test), t(datakangright), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("145"))
dim(test)
colnames(test)
h <- diag(cor(t(test), t(datakangright[,-3]), method = "spearman"))
# test <- allkang %>% dplyr:: select(contains("182"))
# dim(test)
# colnames(test)
# i <- diag(cor(t(test), t(datakangright[,-c(5,6)]), method = "spearman"))
# test <- allkang %>% dplyr:: select(contains("183"))
# dim(test)
# j <- diag(cor(t(test), t(datakangright), method = "spearman"))
test <- allkang %>% dplyr:: select(contains("187"))
dim(test)
k <- diag(cor(t(test), t(datakangright), method = "spearman"))
result <- cbind(a, b, c, e, f, g, h, k)
result <- rowMeans(result)
length(which(result> .446))
datakangright <- subset(datakangright, round(result, digits = 3) >= 0.446)
dim(datakangright)
datakangright <- t(datakangright)
colnames(datakangright) <- subset(kangdata, round(result, digits = 3) >= 0.446)$Genes
datakangright <- datakangright[sort(row.names(datakangright)), ]
datakangright <- (as.data.frame(datakangright))
datakangright$MBP <- NA
datakangright <- datakangright[,apply(datakangright, 2, function(x) !any(is.na(x)))]
a <- dplyr::bind_rows(datakangleft, datakangright)
a <- na.omit(t(a))
a <- t(a)
dataKang <- a
dim(dataKang)

names <- c("A1C", "DFC", "IPC", "ITC", "M1C", "MFC", "OFC", "S1C", "STC", "V1C", "VFC")
datakang <- as.matrix(dataKang)
a<-melt(datakang[1:11,], measure.vars = rownames(dataKang)) 
colnames(a) <- c("region", "Gene", "value")
a<-melt(datakang[12:22,], measure.vars = rownames(dataKang)) 
colnames(a) <- c("region", "Gene", "value")
                                      
pcaKang <- stats::prcomp(dataKang, center=T, scale = T)
KangrawLoadings <- pcaKang$rotation[,1:2] %*% diag(pcaKang$sdev, 2, 2)
KangrotatedLoadings <- varimax(KangrawLoadings)$loadings
KanginvLoadings     <- t(pracma::pinv(KangrotatedLoadings))
Kangscores          <- scale(dataKang) %*% KanginvLoadings

cor(pcaKang$x[1:11,1:2], pcaKang$x[12:22,1:2])
cor(Kangscores[1:11,1:2], Kangscores[12:22,1:2])
t.test(Kangscores[1:11,1], Kangscores[12:22,1])
cor(pcaKang$x[1:22,1:2])
plot(pcaKang$x[1:11,1:2], pcaKang$x[12:22,1:2])

a <- pcaKang$rotation %*% diag(pcaKang$sdev,22,22)
b <- pcaAllen$rotation %*% diag(pcaAllen$sdev,68,68)
a<- KangrotatedLoadings[1:1046,]
b <- rotatedLoadings[1:8325,]
c <- merge(a, b, by = 0) 
dim(c) 
rownames(c) <-c[,1]
c[,1] <- NULL
d <- c[,1:2]
e <- c[,3:4]
f <- (psych::factor.congruence(e,d))
colnames(f) <- gsub("\\.x", "", colnames(f))
rownames(f) <- gsub("\\.y", "", rownames(f))
diag(psych::factor.congruence(d,e))
f <- f[1:5,1:5]
                                      
a <- pcaBrainspan$rotation %*% diag(pcaBrainspan$sdev,11,11)
b <- pcaKang$rotation %*% diag(pcaKang$sdev,22,22)
c <- merge(a, b, by = 0) 
dim(c) 
rownames(c) <-c[,1]
c[,1] <- NULL
d <- c[,1:11]
e <- c[,12:33]
f <- (psych::factor.congruence(e,d))
colnames(f) <- gsub("\\.x", "", colnames(f))
rownames(f) <- gsub("\\.y", "", rownames(f))
diag(psych::factor.congruence(d,e))
f <- f[1:5,1:5]

factorcong[2,] <- c(0.96, 0.63, "Kang et al.")
factorcong <- melt(factorcong, id.vars = "validationset")
factorcong$value <- as.numeric(factorcong$value)
factorcong$variable <-c("PC1","PC1","PC2","PC2")
factorcong$valcomp <- c("1", "1", "2", "3")

pca_kangright <- stats::prcomp(datakangright, scale = T, center = T)
pca_kangleft <- stats::prcomp(datakangleft, scale = T, center = T)
a <- pca_kangright$rotation %*% diag(pca_kangright$sdev,11,11)
b <- pca_kangleft$rotation %*% diag(pca_kangleft$sdev,11,11)
c <- merge(a, b, by = 0) 
dim(c) 
rownames(c) <-c[,1]
c[,1] <- NULL
d <- c[,1:11]
e <- c[,12:22]
f <- (psych::factor.congruence(e,d))
f <- f[1:5,1:5]

pca_kangright <- stats::prcomp(datakangright, scale = T, center = T)
pcaKang <- stats::prcomp(dataKang, scale = T, center = T)
a <- pca_kangright$rotation %*% diag(pca_kangright$sdev,11,11)
b <- pcaKang$rotation %*% diag(pcaKang$sdev,22,22)
c <- merge(a, b, by = 0) 
dim(c) 
rownames(c) <-c[,1]
c[,1] <- NULL
d <- c[,1:11]
e <- c[,12:33]
f <- (psych::factor.congruence(e,d))
f <- f[1:5,1:5]

pca_kangleft <- stats::prcomp(datakangleft, scale = T, center = T)
pcaKang <- stats::prcomp(dataKang, scale = T, center = T)
a <- pca_kangleft$rotation %*% diag(pca_kangleft$sdev,11,11)
b <- pcaKang$rotation %*% diag(pcaKang$sdev,22,22)
c <- merge(a, b, by = 0) 
dim(c) 
rownames(c) <-c[,1]
c[,1] <- NULL
d <- c[,1:11]
e <- c[,12:33]
f <- (psych::factor.congruence(e,d))
f <- f[1:5,1:5]


# abagen -------------------------------------------------
# are the PCs from French and Paus' matrix in line with the matrix produced by the abagen scripts?  (https://github.com/netneurolab/markello_transcriptome and abagen toolbox
abagenDK <- read.csv('abagen_DK.csv')
dim(abagenDK)
rownames(abagenDK) <- abagenDK[,1]
abagenDK <- abagenDK[,-1] 
abagenDK <- abagenDK[order(rownames(abagenDK)),]
abagenDK <- na.omit(abagenDK)
pca_allen <- stats::prcomp(allenabagen, scale = T, center = T)
abagenDK <- abagenDK[, colnames(abagenDK) %in% colnames(dataAllen8000)]
dim(abagenDK)
pca_abagenDK <- stats::prcomp(abagenDK, center =T, scale = T)
fviz_eig(pca_abagenDK, addlabels = T)
cor(pca_allen$x[c(1:66,68),1], pca_abagenDK$x[,1]) # no data was included for the right temporal pole in the matrix resulting from the abagen scripts
cor(pca_allen$x[c(1:66,68),2], pca_abagenDK$x[,2])
cor(pca_allen$x[c(1:66,68),2], pca_abagenDK$x[,3])
abagenrawLoadings <- pca_abagenDK$rotation[,1:2] %*% diag(pca_abagenDK$sdev, 2, 2)
abagenrotatedLoadings <- varimax(abagenrawLoadings)$loadings
abageninvLoadings     <- t(pracma::pinv(abagenrotatedLoadings))
abagenscores          <- scale(abagenDK) %*% abageninvLoadings
cor.test(scores[,1], abagenscores[,1])
cor.test(scores[,2], abagenscores[,2])
plot(scores[,1], abagenscores[,1])
plot(scores[,2], abagenscores[,2])
a <- pca_allen$rotation %*% diag(pca_allen$sdev,68,68)
b <- pca_abagenDK$rotation %*% diag(pca_abagenDK$sdev,67,67)
a <- rotatedLoadings[1:8235,1:2]
b <- abagenrotatedLoadings[1:6166,1:2]
c <- merge(a, b, by = 0) 
dim(c) 
rownames(c) <-c[,1]
c[,1] <- NULL
d <- c[,1:2]
e <- c[,3:4]
f <- (psych::factor.congruence(e,d))
colnames(f) <- gsub("\\.x", "", colnames(f))
rownames(f) <- gsub("\\.y", "", rownames(f))
diag(psych::factor.congruence(d,e))
f <- f[1:5,1:5]

                                      
# different parcellations vs Desikan-Killiany (with French and Paus' gene expression data processing pipeline)
yeo7 <- read.csv('abagen_yeo7.csv')
dim(yeo7)
yeo7 <- yeo7[, colnames(yeo7) %in% colnames(dataAllen8000)]
dim(yeo7)
pca_yeo7 <- stats::prcomp(yeo7, center =T, scale = T)
fviz_eig(pca_yeo7, addlabels = T)
yeo7rawLoadings <- pca_yeo7$rotation[,1:2] %*% diag(pca_yeo7$sdev, 2, 2)
yeo7rotatedLoadings <- varimax(yeo7rawLoadings)$loadings
yeo7invLoadings     <- t(pracma::pinv(yeo7rotatedLoadings))
yeo7scores          <- scale(yeo7) %*% yeo7invLoadings
a <- pca_allen$rotation %*% diag(pca_allen$sdev,68,68)
b <- pca_yeo7$rotation %*% diag(pca_yeo7$sdev,7,7)
#a <- rotatedLoadings[1:8235,]
#b <- yeo7rotatedLoadings[1:8108,]
c <- merge(a, b, by = 0) 
dim(c) 
rownames(c) <-c[,1]
c[,1] <- NULL
d <- c[,1:68]
e <- c[,69:75]
f <- (psych::factor.congruence(e,d))
colnames(f) <- gsub("\\.x", "", colnames(f))
rownames(f) <- gsub("\\.y", "", rownames(f))
diag(psych::factor.congruence(d,e))
f <- f[1:5,1:5]

yeo17 <- read.csv('abagen_yeo17.csv')
dim(yeo17)
yeo17 <- yeo17[, colnames(yeo17) %in% colnames(dataAllen8000)]
dim(yeo17)
pca_yeo17 <- stats::prcomp(yeo17, center =T, scale = T)
fviz_eig(pca_yeo17, addlabels = T)
yeo17rawLoadings <- pca_yeo17$rotation[,1:2] %*% diag(pca_yeo17$sdev, 2, 2)
yeo17rotatedLoadings <- varimax(yeo17rawLoadings)$loadings
yeo17invLoadings     <- t(pracma::pinv(yeo17rotatedLoadings))
yeo17scores          <- scale(yeo17) %*% yeo17invLoadings
a <- pca_allen$rotation %*% diag(pca_allen$sdev,68,68)
b <- pca_yeo17$rotation %*% diag(pca_yeo17$sdev,17,17)
#a <- rotatedLoadings[1:8235,]
#b <- yeo17rotatedLoadings[1:8108,]
c <- merge(a, b, by = 0) 
dim(c) 
rownames(c) <-c[,1]
c[,1] <- NULL
d <- c[,1:68]
e <- c[,69:85]
f <- (psych::factor.congruence(e,d))
colnames(f) <- gsub("\\.x", "", colnames(f))
rownames(f) <- gsub("\\.y", "", rownames(f))
diag(psych::factor.congruence(d,e))
f <- f[1:5,1:5]

a <- pca_yeo7$rotation %*% diag(pca_yeo7$sdev,7,7)
b <- pca_yeo17$rotation %*% diag(pca_yeo17$sdev,17,17)
#a <- rotatedLoadings[1:8235,]
#b <- yeo17rotatedLoadings[1:8108,]
c <- merge(a, b, by = 0) 
dim(c) 
rownames(c) <-c[,1]
c[,1] <- NULL
d <- c[,1::7]
e <- c[,8:24]
f <- (psych::factor.congruence(e,d))
colnames(f) <- gsub("\\.x", "", colnames(f))
rownames(f) <- gsub("\\.y", "", rownames(f))
diag(psych::factor.congruence(d,e))
f <- f[1:5,1:5]


destrieux <- read.csv('abagen_destrieux.csv')
dim(destrieux)
destrieux <- destrieux[, colnames(destrieux) %in% colnames(dataAllen8000)]
dim(destrieux)
pca_destrieux <- stats::prcomp(na.omit(destrieux), center =T, scale = T)
fviz_eig(pca_destrieux, addlabels = T)
destrieuxrawLoadings <- pca_destrieux$rotation[,1:2] %*% diag(pca_destrieux$sdev, 2, 2)
destrieuxrotatedLoadings <- varimax(destrieuxrawLoadings)$loadings
destrieuxinvLoadings     <- t(pracma::pinv(destrieuxrotatedLoadings))
destrieuxscores          <- scale(destrieux) %*% destrieuxinvLoadings
a <- pca_allen$rotation %*% diag(pca_allen$sdev,68,68)
b <- pca_destrieux$rotation %*% diag(pca_destrieux$sdev,134,134)
#a <- rotatedLoadings[1:8235,]
#b <- destrieuxrotatedLoadings[1:8108,]
c <- merge(a, b, by = 0) 
dim(c) 
rownames(c) <-c[,1]
c[,1] <- NULL
d <- c[,1:68]
e <- c[,69:202]
f <- (psych::factor.congruence(e,d))
colnames(f) <- gsub("\\.x", "", colnames(f))
rownames(f) <- gsub("\\.y", "", rownames(f))
diag(psych::factor.congruence(d,e))
f <- f[1:5,1:5]



                                      
#  interpretation: cell types -------------------------------------------
histology <- read_excel("histologygenes.xlsx") # (Zeisel et al. 2015)
rotatedLoadings <- read.csv('rotatedLoadings.csv')
colnames(rotatedLoadings) <- c("RC1", "RC2")
a <- print(rotatedLoadings[,1:2], cutoff = 0)
a <- as.data.frame(a)
a$Genes <- rownames(a)
histo <- merge(a, histology, by = "Genes")
histo <- as.data.frame(histo)
head(histo)
histo <- histo[,c(1, 2,3,4)]
rownames(histology) <- histology$Genes
rawhisto <- merge(t(dataAllen8000), histology, by = 0)
Microglia <- colMeans(rawhisto[which(rawhisto$Cell =="Microglia"),2:69])
Astrocyte <- colMeans(rawhisto[which(rawhisto$Cell =="Astrocyte"),2:69])
CA1Pyramidal <- colMeans(rawhisto[which(rawhisto$Cell =="CA1Pyramidal"),2:69])
Endothelial <- colMeans(rawhisto[which(rawhisto$Cell =="Endothelial"),2:69])
Ependymal <- colMeans(rawhisto[which(rawhisto$Cell =="Ependymal"),2:69])
Interneuron <- colMeans(rawhisto[which(rawhisto$Cell =="Interneuron"),2:69])
Mural <- colMeans(rawhisto[which(rawhisto$Cell =="Mural"),2:69])
Oligodendrocyte <- colMeans(rawhisto[which(rawhisto$Cell =="Oligodendrocyte"),2:69])
S1Pyramidal <- colMeans(rawhisto[which(rawhisto$Cell =="S1Pyramidal"),2:69])
Microglia <- c(scale(Microglia[1:34]), scale(Microglia[35:68]))
Astrocyte <- c(scale(Astrocyte[1:34]), scale(Astrocyte[35:68]))
CA1Pyramidal <- c(scale(CA1Pyramidal[1:34]), scale(CA1Pyramidal[35:68]))
Endothelial <- c(scale(Endothelial[1:34]), scale(Endothelial[35:68]))
Ependymal <- c(scale(Ependymal[1:34]), scale(Ependymal[35:68]))
Interneuron <- c(scale(Interneuron[1:34]), scale(Interneuron[35:68]))
Mural <- c(scale(Mural[1:34]), scale(Mural[35:68]))
Oligodendrocyte <- c(scale(Oligodendrocyte[1:34]), scale(Oligodendrocyte[35:68]))
S1Pyramidal <- c(scale(S1Pyramidal[1:34]), scale(S1Pyramidal[35:68]))
scalescore <- rbind(scale(scores[1:34,1:2]), scale(scores[35:68,1:2])) 
scalescore <- as.data.frame(scalescore)
colnames(scalescore) <- c("Component 1", "Component 2")

                                      
writeout <- cbind(Microglia, Astrocyte, Oligodendrocyte, CA1Pyramidal, Endothelial, Ependymal, Interneuron, Mural, S1Pyramidal)
write.table(writeout, 'meanCelltypes.csv', sep = ",")
testingcorrelations <- read.csv('writescores.csv')
                                      
Microgliacogvol <- lm(scale(Microglia) ~ scale(testingcorrelations$vol_beta) + scalescore[,1] +scalescore[,2])
Microgliacogsa <- lm( scale(Microglia) ~  scale(testingcorrelations$sa_beta)+ scalescore[,1] +scalescore[,2])
Microgliacogthk <- lm( scale(Microglia) ~ scale(testingcorrelations$thk_beta)+ scalescore[,1] +scalescore[,2])
Astrocytecogvol <- lm( scale(Astrocyte)~  scale(testingcorrelations$vol_beta)+ scalescore[,1] +scalescore[,2])
Astrocytecogsa <- lm(scale(Astrocyte)~   scale(testingcorrelations$sa_beta)+ scalescore[,1] +scalescore[,2])
Astrocytecogthk <- lm( scale(Astrocyte)~  scale(testingcorrelations$thk_beta)+ scalescore[,1] +scalescore[,2])
Oligodendrocytecogvol <- lm(scale(Oligodendrocyte) ~ scale(testingcorrelations$vol_beta)+ scalescore[,1] +scalescore[,2])
Oligodendrocytecogsa <- lm(scale(Oligodendrocyte)  ~ scale(testingcorrelations$sa_beta)+ scalescore[,1] +scalescore[,2])
Oligodendrocytecogthk <- lm(scale(Oligodendrocyte)~   scale(testingcorrelations$thk_beta)+ scalescore[,1] +scalescore[,2])
CA1Pyramidalcogvol <- lm( scale(CA1Pyramidal)~  scale(testingcorrelations$vol_beta)+ scalescore[,1] +scalescore[,2])
CA1Pyramidalcogsa <- lm(scale(CA1Pyramidal) ~  scale(testingcorrelations$sa_beta)+ scalescore[,1] +scalescore[,2])
CA1Pyramidalcogthk <- lm(scale(CA1Pyramidal)~  scale(testingcorrelations$thk_beta) + scalescore[,1] +scalescore[,2])
S1Pyramidalcogvol <- lm(scale(S1Pyramidal)  ~ scale(testingcorrelations$vol_beta)+ scalescore[,1] +scalescore[,2])
S1Pyramidalcogsa <- lm(scale(S1Pyramidal)~  scale(testingcorrelations$sa_beta) + scalescore[,1] +scalescore[,2])
S1Pyramidalcogthk <- lm( scale(S1Pyramidal)~ scale(testingcorrelations$thk_beta) + scalescore[,1] +scalescore[,2])
Interneuroncogvol <- lm( scale(Interneuron)~ scale(testingcorrelations$vol_beta)+ scalescore[,1] +scalescore[,2])
Interneuroncogsa <- lm(scale(Interneuron)~  scale(testingcorrelations$sa_beta) + scalescore[,1] +scalescore[,2])
Interneuroncogthk <- lm(scale(Interneuron)~  scale(testingcorrelations$thk_beta) + scalescore[,1] +scalescore[,2])
Ependymalcogvol <- lm(scale(Ependymal) ~  scale(testingcorrelations$vol_beta)+ scalescore[,1] +scalescore[,2])
Ependymalcogsa <- lm(scale(Ependymal) ~ scale(testingcorrelations$sa_beta) + scalescore[,1] +scalescore[,2])
Ependymalcogthk <- lm(scale(Ependymal) ~  scale(testingcorrelations$thk_beta)+ scalescore[,1] +scalescore[,2])
Endothelialcogvol <- lm(scale(Endothelial) ~ scale(testingcorrelations$vol_beta) + scalescore[,1] +scalescore[,2])
Endothelialcogsa <- lm(scale(Endothelial)~  scale(testingcorrelations$sa_beta) + scalescore[,1] +scalescore[,2])
Endothelialcogthk <- lm(scale(Endothelial) ~  scale(testingcorrelations$thk_beta)+ scalescore[,1] +scalescore[,2])
Muralcogvol <- lm( scale(Mural)  ~ scale(testingcorrelations$vol_beta)+ scalescore[,1] +scalescore[,2])
Muralcogsa <- lm(scale(Mural) ~ scale(testingcorrelations$sa_beta) + scalescore[,1] +scalescore[,2])
Muralcogthk <- lm(scale(Mural) ~  scale(testingcorrelations$thk_beta) + scalescore[,1] +scalescore[,2])

coef(summary(Microgliacogvol))[2,1:4] #
coef(summary(Microgliacogsa))[2,1:4] #
coef(summary(Microgliacogthk))[2,4]
coef(summary(Oligodendrocytecogvol))[2,4]
coef(summary(Oligodendrocytecogsa))[2,4]
coef(summary(Oligodendrocytecogthk))[2,4] # 
coef(summary(Muralcogvol))[2,4]
coef(summary(Muralcogsa))[2,4]
coef(summary(Muralcogthk))[2,4] 
coef(summary(CA1Pyramidalcogvol))[2,4] ##
coef(summary(CA1Pyramidalcogsa))[2,4]
coef(summary(CA1Pyramidalcogthk))[2,4] ##
coef(summary(S1Pyramidalcogvol))[2,4]
coef(summary(S1Pyramidalcogsa))[2,4]
coef(summary(S1Pyramidalcogthk))[2,4]
coef(summary(Endothelialcogvol))[2,4] ##
coef(summary(Endothelialcogsa))[2,4]
coef(summary(Endothelialcogthk))[2,4] 
coef(summary(Ependymalcogvol))[2,1:4] 
coef(summary(Ependymalcogsa))[2,4]
coef(summary(Ependymalcogthk))[2,1:4] 
coef(summary(Interneuroncogvol))[2,4]
coef(summary(Interneuroncogsa))[2,4]
coef(summary(Interneuroncogthk))[2,4]
coef(summary(Astrocytecogvol))[2,4]
coef(summary(Astrocytecogsa))[2,4]
coef(summary(Astrocytecogthk))[2,4]
all <- c(coef(summary(Microgliacogvol))[2,4],coef(summary(Microgliacogsa))[2,4] ,coef(summary(Microgliacogthk))[2,4],coef(summary(Oligodendrocytecogvol))[2,4],coef(summary(Oligodendrocytecogsa))[2,4],coef(summary(Oligodendrocytecogthk))[2,4] ,coef(summary(Muralcogvol))[2,4],coef(summary(Muralcogsa))[2,4],coef(summary(Muralcogthk))[2,4] ,coef(summary(CA1Pyramidalcogvol))[2,4],coef(summary(CA1Pyramidalcogsa))[2,4],coef(summary(CA1Pyramidalcogthk))[2,4],coef(summary(S1Pyramidalcogvol))[2,4], coef(summary(S1Pyramidalcogsa))[2,4],coef(summary(S1Pyramidalcogthk))[2,4],coef(summary(Endothelialcogvol))[2,4] ,coef(summary(Endothelialcogsa))[2,4] ,coef(summary(Endothelialcogthk))[2,4] ,coef(summary(Ependymalcogvol))[2,4] ,coef(summary(Ependymalcogsa))[2,4],coef(summary(Ependymalcogthk))[2,4],
coef(summary(Interneuroncogvol))[2,4],coef(summary(Interneuroncogsa))[2,4],coef(summary(Interneuroncogthk))[2,4],coef(summary(Astrocytecogvol))[2,4] ,coef(summary(Astrocytecogsa))[2,4],coef(summary(Astrocytecogthk))[2,4] )
all <- as.data.frame(all)
which(p.adjust(all$all, method = "BH") < .05)

all <- rbind(coef(summary(Microgliacogvol))[2,1:4],coef(summary(Microgliacogsa))[2,1:4] ,coef(summary(Microgliacogthk))[2,1:4],coef(summary(Oligodendrocytecogvol))[2,1:4],coef(summary(Oligodendrocytecogsa))[2,1:4],coef(summary(Oligodendrocytecogthk))[2,1:4] ,coef(summary(Muralcogvol))[2,1:4],coef(summary(Muralcogsa))[2,1:4],coef(summary(Muralcogthk))[2,1:4] ,coef(summary(CA1Pyramidalcogvol))[2,1:4],coef(summary(CA1Pyramidalcogsa))[2,1:4],coef(summary(CA1Pyramidalcogthk))[2,1:4],coef(summary(S1Pyramidalcogvol))[2,1:4], coef(summary(S1Pyramidalcogsa))[2,1:4],coef(summary(S1Pyramidalcogthk))[2,1:4],coef(summary(Endothelialcogvol))[2,1:4] ,coef(summary(Endothelialcogsa))[2,1:4] ,coef(summary(Endothelialcogthk))[2,1:4] ,coef(summary(Ependymalcogvol))[2,1:4] ,coef(summary(Ependymalcogsa))[2,1:4],coef(summary(Ependymalcogthk))[2,1:4],
coef(summary(Interneuroncogvol))[2,1:4],coef(summary(Interneuroncogsa))[2,1:4],coef(summary(Interneuroncogthk))[2,1:4],coef(summary(Astrocytecogvol))[2,1:4] ,coef(summary(Astrocytecogsa))[2,1:4],coef(summary(Astrocytecogthk))[2,1:4] )
all <- as.data.frame(all)
all$Cell <- c(rep("Microglia",3), rep("Oligodendrocyte",3), rep("Mural",3), rep("CA1Pyramidal",3), rep("S1Pyramidal",3), rep("Endothelial",3), rep("Ependymal",3), rep("Interneuron",3), rep("Astrocyte",3))
all$comparison <- rep(c("Volume", "Area", "Thickness"), 9)
all <- all[order(all$Cell),]
all <- as.data.frame(all)
colnames(all)[4] <- "p"
all$q <- p.adjust(all$p, method = "BH") 
write.table(all, 'celltypescog.csv', sep = ',', row.names = F, col.names = T)

histomelt <- melt(histo, id.vars = c("Genes", "Cell"))
head(histomelt)
colnames(histomelt) <- c("Genes", "Cell", "component", "loading")
str(histomelt)
histomelt$Cell <- as.factor(histomelt$Cell)
histomelt <- as.data.frame(histomelt)
incontext <- a
incontext <- merge(incontext, histology, by = "Genes", all.x = T)
incontext$Cell[is.na(incontext$Cell)] <- "Unclassified"
str(incontext)
incontext$Cell <- as.factor(incontext$Cell)
#plots
c1 <- histomelt[which(histomelt$component == "RC1"),]  
c1$Cell <- as.factor(c1$Cell)
c2 <- histomelt[which(histomelt$component == "RC2"),]  
c2$Cell <- as.factor(c2$Cell)

incontextmelt <- melt(incontext[,1:4], id.vars = c("Genes", "Cell"))
colnames(incontextmelt) <- c("Genes", "Cell","component", "loading")
c1 <- incontextmelt[which(incontextmelt$component == "RC1"),]  
c1$Cell <- as.factor(c1$Cell)
c2 <- incontextmelt[which(incontextmelt$component == "RC2"),]  
c2$Cell <- as.factor(c2$Cell)


#kruskal wallis 

kruskal.test(loading ~ Cell, data = c1)
dunntest1 <- dunnTest(c1$loading, c1$Cell)
kruskal.test(loading ~ Cell, data = c2)
dunntest2 <- dunnTest(c2$loading, c2$Cell)
dunntest1$res
dunntest2$res

                                    

# which individual gene expression profiles are associated with the g-morphometry profiles, after controlling for the first two principal components?
testingcorrelations <- read.csv('testingcorrelations.csv', sep = " ")
scaledlh <- scale(dataAllen8000[1:34,], scale =T)
scaledrh <- scale(dataAllen8000[35:68,], scale = T)
scaledgene <- rbind(scaledlh, scaledrh)
scaledgene <- as.data.frame(scaledgene)
scaledgene$region <- rownames(scaledgene)
scaledgene$region <- gsub("ctx.lh.", "", scaledgene$region )
scaledgene$region <- gsub("ctx.rh.", "", scaledgene$region )

testingcorrelations <- testingcorrelations[c(seq(1, 68, 2), seq(2, 68, 2)),]


cor(testingcorrelations$vol_beta, abs(testingcorrelations$scaledscoresPC1)) 
cor(testingcorrelations$vol_beta, (testingcorrelations$scaledscoresPC1)) 
cor(testingcorrelations$vol_beta, abs(testingcorrelations$scaledscoresPC2)) 
cor(testingcorrelations$vol_beta, (testingcorrelations$scaledscoresPC2)) 
cor(testingcorrelations$sa_beta, abs(testingcorrelations$scaledscoresPC1)) 
cor(testingcorrelations$sa_beta, (testingcorrelations$scaledscoresPC1)) 
cor(testingcorrelations$sa_beta, abs(testingcorrelations$scaledscoresPC2)) 
cor(testingcorrelations$sa_beta, (testingcorrelations$scaledscoresPC2)) 
cor(testingcorrelations$thk_beta, abs(testingcorrelations$scaledscoresPC1)) 
cor(testingcorrelations$thk_beta, (testingcorrelations$scaledscoresPC1)) 
cor(testingcorrelations$thk_beta, abs(testingcorrelations$scaledscoresPC2)) 
cor(testingcorrelations$thk_beta, (testingcorrelations$scaledscoresPC2)) 


volumemass <- matrix(0, 8235, 6)
for (i in 1:8235) {
  first <- lm(scale(scaledgene[,i]) ~ scale(testingcorrelations$vol_beta)+ testingcorrelations$scaledscoresPC1 + testingcorrelations$scaledscoresPC2)
  volumemass[i,1:4] <- summary(lm.beta(first))$coefficients[2,2:5]
	volumemass[i,5] <- cor(scale(testingcorrelations$vol_beta), scale(scaledgene[,i]))
	volumemass[i,6] <- summary(lm.beta(first))$r.squared
}


volumemass <- as.data.frame(volumemass)
volumemass$Genes <- colnames(dataAllen8000)
colnames(volumemass)[5] <- "cor"
volumemass <- volumemass[order(volumemass[,4]),]
notcorrectedsigvol <- volumemass[which(volumemass[,4] < .05),]
dim(notcorrectedsigvol)
q <- p.adjust(volumemass[,4], method = "BH")
q <- as.data.frame(q)
q$Genes <- volumemass$Genes
correctedsigvol <- merge(notcorrectedsigvol, q, by = "Genes")
notcorrectedsigvol <- correctedsigvol
correctedsigvol <- correctedsigvol[which(correctedsigvol$q < .05),]
dim(correctedsigvol)
                                      
areamass <- matrix(0, 8235, 6)
for (i in 1:8235) {
  first <- lm(scale(scaledgene[,i])~  scale(testingcorrelations$sa_beta)+ testingcorrelations$scaledscoresPC1 + testingcorrelations$scaledscoresPC2)
  areamass[i,1:4] <- summary(lm.beta(first))$coefficients[2,2:5]
	areamass[i,5] <- cor(scale(testingcorrelations$sa_beta), scale(scaledgene[,i]))
	areamass[i,6] <- summary(lm.beta(first))$r.squared
}

areamass <- as.data.frame(areamass)
areamass$Genes <- colnames(dataAllen8000)
colnames(areamass)[5] <- "cor"
areamass <- areamass[order(areamass[,4]),]
notcorrectedsigsa <- areamass[which(areamass[,4] < .05),]
dim(notcorrectedsigsa)
q <- p.adjust(areamass[,4], method = "BH")
q <- as.data.frame(q)
q$Genes <- areamass$Genes
correctedsigsa <- merge(notcorrectedsigsa, q, by = "Genes")
notcorrectedsigsa <- correctedsigsa
correctedsigsa <- correctedsigsa[which(correctedsigsa$q < .05),]
dim(correctedsigsa)

thickmass <- matrix(0, 8235, 6)
for (i in 1:8235) {
  first <- lm(scale(scaledgene[,i])~  scale(testingcorrelations$thk_beta)+ testingcorrelations$scaledscoresPC1 + testingcorrelations$scaledscoresPC2)
  thickmass[i,1:4] <- summary(lm.beta(first))$coefficients[2,2:5]
	thickmass[i,5] <- cor(scale(testingcorrelations$thk_beta), scale(scaledgene[,i]))
	thickmass[i,6] <- summary(lm.beta(first))$r.squared
}

thickmass <- as.data.frame(thickmass)
thickmass$Genes <- colnames(dataAllen8000)
colnames(thickmass)[5] <- "cor"
thickmass <- thickmass[order(thickmass[,4]),]
notcorrectedsigthk <- thickmass[which(thickmass[,4] < .05),]
dim(notcorrectedsigthk)
q <- p.adjust(thickmass[,4], method = "BH")
q <- as.data.frame(q)
q$Genes <- thickmass$Genes
correctedsigthk <- merge(notcorrectedsigthk, q, by = "Genes")
notcorrectedsigthk <- correctedsigthk
correctedsigthk <- correctedsigthk[which(correctedsigthk$q < .05),]
dim(correctedsigthk)


# which genes are have significant associations with g-volume, g-area AND g-thickness profiles?
inallthree <- rbind(volumemass,areamass, thickmass)
inallthree$q <- p.adjust(inallthree$V4, method = "BH")
backvolumemass <- inallthree[1:8235,]
backareamass <- inallthree[8236:16470,]
backthickmass <- inallthree[16471:24705,]
inallthree <- merge(backvolumemass, backareamass, by = "Genes")
inallthree <- merge(inallthree, backthickmass, by = "Genes")
dim(inallthree[which(inallthree$q < .05  & inallthree$q.y < .05),])
dim(inallthree[which(inallthree$q < .05  & inallthree$q.x < .05),])
dim(inallthree[which(inallthree$q.x < .05  & inallthree$q.y < .05),])
inallthree <- inallthree[which(inallthree$q < .05 & inallthree$q.x < .05 & inallthree$q.y < .05),]
write.table(inallthree, "inallthree.csv", sep = ",")




