# The Desikan-Killiany regional-g associations are found for each cohort: UKB, STRADL and LBC1936, and then meta-analysed. 
library(lavaan)
library(tidyverse)
library(ggplot2)
library(ggdist)
library(tidyquant) 
library(reshape2)
library(gapminder)
library(patchwork)
library(gridExtra)
library(cowplot)
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")
library(robumeta)
library(metafor)
library(dplyr)
library(ggforestplot)



outliers <- function(x, SD) {
a <- which(x> ((SD * sd(x, na.rm=T))+ (mean(x, na.rm=T))))
b <- which(x< ((-SD * sd(x, na.rm=T))+ (mean(x, na.rm=T))))
c <- c(a,b)
result <- c
return(result)
}

######## UKB


UKB <- read.csv("") # This file contains demographic information, the Desikan-Killiany 68 region parcellation values (columns) for participants who took part in MRI imaging in instance 2 (~40,000) and the cognitive test results from the 11 tests we included to estimate a latent g factor 

UKB$cog_trailB_log[which(UKB$cog_trailB_log == 0)] <- NA # prepare cognitive test data
UKB$cog_prosmem[which(UKB$cog_prosmem == 2)] <- 0
UKB$brokenletters_w2[which(UKB$brokenletters_w2 < 20)] <- NA
UKB$cog_numeric_memory[which(UKB$cog_numeric_memory == -1)] <- NA

neuroexclusions <- read.csv('UKB_neuroexclusions.csv', sep = " ")
colnames(neuroexclusions)[1] <- "ID"
UKB <- UKB[-which(UKB$ID %in% neuroexclusions[,1]),]

newUKB <- UKB 
for (i in 24:length(UKB)) { # index the regional brain data columns, and remove outliers for each column
	j = UKB[,colnames(UKB[i])]
	j[outliers(j, SD = 4)]=NA
	newUKB[,i] <- j
}
NAcountbeforeOutliers <- sapply(UKB, function(x) sum(is.na(x))) 
NAcountafterOutliers <- sapply(newUKB, function(x) sum(is.na(x)))
UKB <- newUKB


# create g scores 
theUKBCogdata = data.frame(ID = UKB$ID, cog_RT_log = UKB$cog_RT_log, cog_numeric_memory = UKB$cog_numeric_memory, cog_fluid_intelligence = UKB$cog_fluid_intelligence, cog_trailB_log = UKB$cog_trailB_log, cog_matrix_pattern = UKB$cog_matrix_pattern_correct,cog_tower = UKB$cog_tower,cog_digsym = UKB$cog_digsym, cog_pairedAss = UKB$cog_pairedAss, cog_prosmem = UKB$cog_prosmem, cog_pairsmatch_incorrect_log = UKB$cog_pairsmatch_incorrect_log, cog_picturevocab = UKB$picturevocab_w2, ageMRI = UKB$ageyears_MRI/100, sex = UKB$sex)
theUKBCogdata_models <- theUKBCogdata
theUKBCogdata_models[,c(8)] <- theUKBCogdata_models[,c(8)]/10 # allows variances of different cognitive tests to be on similar scales, so that lavaan runs
theUKBCogdata_models[,c(7)] <- theUKBCogdata_models[,c(7)]/10
par(mfrow = c(1, 1))
UKBcogmodel <- 'g =~cog_RT_log +cog_numeric_memory + cog_fluid_intelligence + cog_trailB_log + cog_matrix_pattern + cog_tower + cog_digsym + cog_pairsmatch_incorrect_log +cog_prosmem +cog_pairedAss + cog_picturevocab
cog_RT_log ~ ageMRI + sex
cog_numeric_memory ~ ageMRI + sex
cog_fluid_intelligence ~ ageMRI + sex
cog_trailB_log ~ ageMRI + sex 
cog_matrix_pattern ~ ageMRI + sex
cog_tower ~ ageMRI + sex
cog_digsym ~ ageMRI + sex
cog_pairsmatch_incorrect_log ~ ageMRI + sex
cog_prosmem ~ ageMRI + sex
cog_pairedAss ~ ageMRI + sex
cog_picturevocab ~ ageMRI + sex
'
UKBcogfit <- sem(UKBcogmodel, data=theUKBCogdata_models, missing="fiml.x")
summary(UKBcogfit, fit.measures=TRUE, standardized=T)

GpredictUKB <- lavPredict(UKBcogfit, theUKBCogdata_models)
GpredictUKB <- cbind(theUKBCogdata$ID, GpredictUKB)
colnames(GpredictUKB) <- c("ID", "g")
GpredictUKB <- as.data.frame(GpredictUKB)
GpredictUKB$g <- scale(GpredictUKB$g)
			       
# the desikan-killiany regions have slightly different names between freesurfer versions. Therefore, while this is a long write-out, it makes sure that the variable names will match up between cohorts (which is helpful later).
# Here, the variables are also scaled so that the variances are on similar scales, which is needed for lavaan to run successfully. 
theUKBdata<-data.frame(ID=UKB$ID,ageMRI=UKB$ageyears_MRI,dontuse=UKB$ageyears_MRI,lag=rep(0,dim(UKB)[1]),sex=UKB$sex,atrophy=UKB$ICV,cohort=rep("UKB",dim(UKB)[1]),freesurfer=rep("v6",dim(UKB)[1]),site=UKB$assCtr,X=UKB$headposX,Y=UKB$headposY,Z=UKB$headposZ/10,lh_vol_bankssts=UKB$lh_bankssts_volume/1000,lh_vol_caudalanteriorcingulate=UKB$lh_caudalanteriorcingulate_volume/1000,lh_vol_caudalmiddlefrontal=UKB$lh_caudalmiddlefrontal_volume/1000,lh_vol_cuneus=UKB$lh_cuneus_volume/1000,lh_vol_entorhinal=UKB$lh_entorhinal_volume/1000,lh_vol_fusiform=UKB$lh_fusiform_volume/1000,lh_vol_inferiorparietal=UKB$lh_inferiorparietal_volume/1000,lh_vol_inferiortemporal=UKB$lh_inferiortemporal_volume/1000,lh_vol_isthmus=UKB$lh_isthmuscingulate_volume/1000,lh_vol_lateraloccipital=UKB$lh_lateraloccipital_volume/1000,lh_vol_lateralorbitofrontal=UKB$lh_lateralorbitofrontal_volume/1000,lh_vol_lingual=UKB$lh_lingual_volume/1000,lh_vol_medialorbitofrontal=UKB$lh_medialorbitofrontal_volume/1000,lh_vol_middletemporal=UKB$lh_middletemporal_volume/1000,lh_vol_parahippocampal=UKB$lh_parahippocampal_volume/1000,lh_vol_paracentral=UKB$lh_paracentral_volume/1000,lh_vol_parsopercularis=UKB$lh_parsopercularis_volume/1000,lh_vol_parsorbitalis=UKB$lh_parsorbitalis_volume/1000,lh_vol_parstriangularis=UKB$lh_parstriangularis_volume/1000,lh_vol_pericalcarine=UKB$lh_pericalcarine_volume/1000,lh_vol_postcentral=UKB$lh_postcentral_volume/1000,lh_vol_posteriorcingulate=UKB$lh_posteriorcingulate_volume/1000,lh_vol_precentral=UKB$lh_precentral_volume/1000,lh_vol_precuneus=UKB$lh_precuneus_volume/1000,lh_vol_rostralanteriorcingulate=UKB$lh_rostralanteriorcingulate_volume/1000,lh_vol_rostralmiddlefrontal=UKB$lh_rostralmiddlefrontal_volume/1000,lh_vol_superiorfrontal=UKB$lh_superiorfrontal_volume/1000,lh_vol_superiorparietal=UKB$lh_superiorparietal_volume/1000,lh_vol_superiortemporal=UKB$lh_superiortemporal_volume/1000,lh_vol_supramarginal=UKB$lh_supramarginal_volume/1000,lh_vol_frontalpole=UKB$lh_frontalpole_volume/1000,lh_vol_transversetemporal=UKB$lh_transversetemporal_volume/1000,lh_vol_insula=UKB$lh_insula_volume/1000,lh_sa_bankssts=UKB$lh_bankssts_area/100,lh_sa_caudalanteriorcingulate=UKB$lh_caudalanteriorcingulate_area/100,lh_sa_caudalmiddlefrontal=UKB$lh_caudalmiddlefrontal_area/100,lh_sa_cuneus=UKB$lh_cuneus_area/100,lh_sa_entorhinal=UKB$lh_entorhinal_area/100,lh_sa_fusiform=UKB$lh_fusiform_area/100,lh_sa_inferiorparietal=UKB$lh_inferiorparietal_area/100,lh_sa_inferiortemporal=UKB$lh_inferiortemporal_area/100,lh_sa_isthmus=UKB$lh_isthmuscingulate_area/100,lh_sa_lateraloccipital=UKB$lh_lateraloccipital_area/100,lh_sa_lateralorbitofrontal=UKB$lh_lateralorbitofrontal_area/100,lh_sa_lingual=UKB$lh_lingual_area/100,lh_sa_medialorbitofrontal=UKB$lh_medialorbitofrontal_area/100,lh_sa_middletemporal=UKB$lh_middletemporal_area/100,lh_sa_parahippocampal=UKB$lh_parahippocampal_area/100,lh_sa_paracentral=UKB$lh_paracentral_area/100,lh_sa_parsopercularis=UKB$lh_parsopercularis_area/100,lh_sa_parsorbitalis=UKB$lh_parsorbitalis_area/100,lh_sa_parstriangularis=UKB$lh_parstriangularis_area/100,lh_sa_pericalcarine=UKB$lh_pericalcarine_area/100,lh_sa_postcentral=UKB$lh_postcentral_area/100,lh_sa_posteriorcingulate=UKB$lh_posteriorcingulate_area/100,lh_sa_precentral=UKB$lh_precentral_area/100,lh_sa_precuneus=UKB$lh_precuneus_area/100,lh_sa_rostralanteriorcingulate=UKB$lh_rostralanteriorcingulate_area/100,lh_sa_rostralmiddlefrontal=UKB$lh_rostralmiddlefrontal_area/100,lh_sa_superiorfrontal=UKB$lh_superiorfrontal_area/100,lh_sa_superiorparietal=UKB$lh_superiorparietal_area/100,lh_sa_superiortemporal=UKB$lh_superiortemporal_area/100,lh_sa_supramarginal=UKB$lh_supramarginal_area/100,lh_sa_frontalpole=UKB$lh_frontalpole_area/100,lh_sa_transversetemporal=UKB$lh_transversetemporal_area/100,lh_sa_insula=UKB$lh_insula_area/100,lh_thk_bankssts=UKB$lh_bankssts_thickness*10,lh_thk_caudalanteriorcingulate=UKB$lh_caudalanteriorcingulate_thickness*10,lh_thk_caudalmiddlefrontal=UKB$lh_caudalmiddlefrontal_thickness*10,lh_thk_cuneus=UKB$lh_cuneus_thickness*10,lh_thk_entorhinal=UKB$lh_entorhinal_thickness*10,lh_thk_fusiform=UKB$lh_fusiform_thickness*10,lh_thk_inferiorparietal=UKB$lh_inferiorparietal_thickness*10,lh_thk_inferiortemporal=UKB$lh_inferiortemporal_thickness*10,lh_thk_isthmus=UKB$lh_isthmuscingulate_thickness*10,lh_thk_lateraloccipital=UKB$lh_lateraloccipital_thickness*10,lh_thk_lateralorbitofrontal=UKB$lh_lateralorbitofrontal_thickness*10,lh_thk_lingual=UKB$lh_lingual_thickness*10,lh_thk_medialorbitofrontal=UKB$lh_medialorbitofrontal_thickness*10,lh_thk_middletemporal=UKB$lh_middletemporal_thickness*10,lh_thk_parahippocampal=UKB$lh_parahippocampal_thickness*10,lh_thk_paracentral=UKB$lh_paracentral_thickness*10,lh_thk_parsopercularis=UKB$lh_parsopercularis_thickness*10,lh_thk_parsorbitalis=UKB$lh_parsorbitalis_thickness*10,lh_thk_parstriangularis=UKB$lh_parstriangularis_thickness*10,lh_thk_pericalcarine=UKB$lh_pericalcarine_thickness*10,lh_thk_postcentral=UKB$lh_postcentral_thickness*10,lh_thk_posteriorcingulate=UKB$lh_posteriorcingulate_thickness*10,lh_thk_precentral=UKB$lh_precentral_thickness*10,lh_thk_precuneus=UKB$lh_precuneus_thickness*10,lh_thk_rostralanteriorcingulate=UKB$lh_rostralanteriorcingulate_thickness*10,lh_thk_rostralmiddlefrontal=UKB$lh_rostralmiddlefrontal_thickness*10,lh_thk_superiorfrontal=UKB$lh_superiorfrontal_thickness*10,lh_thk_superiorparietal=UKB$lh_superiorparietal_thickness*10,lh_thk_superiortemporal=UKB$lh_superiortemporal_thickness*10,lh_thk_supramarginal=UKB$lh_supramarginal_thickness*10,lh_thk_frontalpole=UKB$lh_frontalpole_thickness*10,lh_thk_transversetemporal=UKB$lh_transversetemporal_thickness*10,lh_thk_insula=UKB$lh_insula_thickness*10,rh_vol_bankssts=UKB$rh_bankssts_volume/1000,rh_vol_caudalanteriorcingulate=UKB$rh_caudalanteriorcingulate_volume/1000,rh_vol_caudalmiddlefrontal=UKB$rh_caudalmiddlefrontal_volume/1000,rh_vol_cuneus=UKB$rh_cuneus_volume/1000,rh_vol_entorhinal=UKB$rh_entorhinal_volume/1000,rh_vol_fusiform=UKB$rh_fusiform_volume/1000,rh_vol_inferiorparietal=UKB$rh_inferiorparietal_volume/1000,rh_vol_inferiortemporal=UKB$rh_inferiortemporal_volume/1000,rh_vol_isthmus=UKB$rh_isthmuscingulate_volume/1000,rh_vol_lateraloccipital=UKB$rh_lateraloccipital_volume/1000,rh_vol_lateralorbitofrontal=UKB$rh_lateralorbitofrontal_volume/1000,rh_vol_lingual=UKB$rh_lingual_volume/1000,rh_vol_medialorbitofrontal=UKB$rh_medialorbitofrontal_volume/1000,rh_vol_middletemporal=UKB$rh_middletemporal_volume/1000,rh_vol_parahippocampal=UKB$rh_parahippocampal_volume/1000,rh_vol_paracentral=UKB$rh_paracentral_volume/1000,rh_vol_parsopercularis=UKB$rh_parsopercularis_volume/1000,rh_vol_parsorbitalis=UKB$rh_parsorbitalis_volume/1000,rh_vol_parstriangularis=UKB$rh_parstriangularis_volume/1000,rh_vol_pericalcarine=UKB$rh_pericalcarine_volume/1000,rh_vol_postcentral=UKB$rh_postcentral_volume/1000,rh_vol_posteriorcingulate=UKB$rh_posteriorcingulate_volume/1000,rh_vol_precentral=UKB$rh_precentral_volume/1000,rh_vol_precuneus=UKB$rh_precuneus_volume/1000,rh_vol_rostralanteriorcingulate=UKB$rh_rostralanteriorcingulate_volume/1000,rh_vol_rostralmiddlefrontal=UKB$rh_rostralmiddlefrontal_volume/1000,rh_vol_superiorfrontal=UKB$rh_superiorfrontal_volume/1000,rh_vol_superiorparietal=UKB$rh_superiorparietal_volume/1000,rh_vol_superiortemporal=UKB$rh_superiortemporal_volume/1000,rh_vol_supramarginal=UKB$rh_supramarginal_volume/1000,rh_vol_frontalpole=UKB$rh_frontalpole_volume/1000,rh_vol_transversetemporal=UKB$rh_transversetemporal_volume/1000,rh_vol_insula=UKB$rh_insula_volume/1000,rh_sa_bankssts=UKB$rh_bankssts_area/100,rh_sa_caudalanteriorcingulate=UKB$rh_caudalanteriorcingulate_area/100,rh_sa_caudalmiddlefrontal=UKB$rh_caudalmiddlefrontal_area/100,rh_sa_cuneus=UKB$rh_cuneus_area/100,rh_sa_entorhinal=UKB$rh_entorhinal_area/100,rh_sa_fusiform=UKB$rh_fusiform_area/100,rh_sa_inferiorparietal=UKB$rh_inferiorparietal_area/100,rh_sa_inferiortemporal=UKB$rh_inferiortemporal_area/100,rh_sa_isthmus=UKB$rh_isthmuscingulate_area/100,rh_sa_lateraloccipital=UKB$rh_lateraloccipital_area/100,rh_sa_lateralorbitofrontal=UKB$rh_lateralorbitofrontal_area/100,rh_sa_lingual=UKB$rh_lingual_area/100,rh_sa_medialorbitofrontal=UKB$rh_medialorbitofrontal_area/100,rh_sa_middletemporal=UKB$rh_middletemporal_area/100,rh_sa_parahippocampal=UKB$rh_parahippocampal_area/100,rh_sa_paracentral=UKB$rh_paracentral_area/100,rh_sa_parsopercularis=UKB$rh_parsopercularis_area/100,rh_sa_parsorbitalis=UKB$rh_parsorbitalis_area/100,rh_sa_parstriangularis=UKB$rh_parstriangularis_area/100,rh_sa_pericalcarine=UKB$rh_pericalcarine_area/100,rh_sa_postcentral=UKB$rh_postcentral_area/100,rh_sa_posteriorcingulate=UKB$rh_posteriorcingulate_area/100,rh_sa_precentral=UKB$rh_precentral_area/100,rh_sa_precuneus=UKB$rh_precuneus_area/100,rh_sa_rostralanteriorcingulate=UKB$rh_rostralanteriorcingulate_area/100,rh_sa_rostralmiddlefrontal=UKB$rh_rostralmiddlefrontal_area/100,rh_sa_superiorfrontal=UKB$rh_superiorfrontal_area/100,rh_sa_superiorparietal=UKB$rh_superiorparietal_area/100,rh_sa_superiortemporal=UKB$rh_superiortemporal_area/100,rh_sa_supramarginal=UKB$rh_supramarginal_area/100,rh_sa_frontalpole=UKB$rh_frontalpole_area/100,rh_sa_transversetemporal=UKB$rh_transversetemporal_area/100,rh_sa_insula=UKB$rh_insula_area/100,rh_thk_bankssts=UKB$rh_bankssts_thickness*10,rh_thk_caudalanteriorcingulate=UKB$rh_caudalanteriorcingulate_thickness*10,rh_thk_caudalmiddlefrontal=UKB$rh_caudalmiddlefrontal_thickness*10,rh_thk_cuneus=UKB$rh_cuneus_thickness*10,rh_thk_entorhinal=UKB$rh_entorhinal_thickness*10,rh_thk_fusiform=UKB$rh_fusiform_thickness*10,rh_thk_inferiorparietal=UKB$rh_inferiorparietal_thickness*10,rh_thk_inferiortemporal=UKB$rh_inferiortemporal_thickness*10,rh_thk_isthmus=UKB$rh_isthmuscingulate_thickness*10,rh_thk_lateraloccipital=UKB$rh_lateraloccipital_thickness*10,rh_thk_lateralorbitofrontal=UKB$rh_lateralorbitofrontal_thickness*10,rh_thk_lingual=UKB$rh_lingual_thickness*10,rh_thk_medialorbitofrontal=UKB$rh_medialorbitofrontal_thickness*10,rh_thk_middletemporal=UKB$rh_middletemporal_thickness*10,rh_thk_parahippocampal=UKB$rh_parahippocampal_thickness*10,rh_thk_paracentral=UKB$rh_paracentral_thickness*10,rh_thk_parsopercularis=UKB$rh_parsopercularis_thickness*10,rh_thk_parsorbitalis=UKB$rh_parsorbitalis_thickness*10,rh_thk_parstriangularis=UKB$rh_parstriangularis_thickness*10,rh_thk_pericalcarine=UKB$rh_pericalcarine_thickness*10,rh_thk_postcentral=UKB$rh_postcentral_thickness*10,rh_thk_posteriorcingulate=UKB$rh_posteriorcingulate_thickness*10,rh_thk_precentral=UKB$rh_precentral_thickness*10,rh_thk_precuneus=UKB$rh_precuneus_thickness*10,rh_thk_rostralanteriorcingulate=UKB$rh_rostralanteriorcingulate_thickness*10,rh_thk_rostralmiddlefrontal=UKB$rh_rostralmiddlefrontal_thickness*10,rh_thk_superiorfrontal=UKB$rh_superiorfrontal_thickness*10,rh_thk_superiorparietal=UKB$rh_superiorparietal_thickness*10,rh_thk_superiortemporal=UKB$rh_superiortemporal_thickness*10,rh_thk_supramarginal=UKB$rh_supramarginal_thickness*10,rh_thk_frontalpole=UKB$rh_frontalpole_thickness*10,rh_thk_transversetemporal=UKB$rh_transversetemporal_thickness*10,rh_thk_insula=UKB$rh_insula_thickness*10,lh_vol_temporalpole=UKB$lh_temporalpole_volume/1000,rh_vol_temporalpole=UKB$rh_temporalpole_volume/1000,lh_sa_temporalpole=UKB$lh_temporalpole_area/100,rh_sa_temporalpole=UKB$rh_temporalpole_area/100,lh_thk_temporalpole=UKB$lh_temporalpole_thickness/10,rh_thk_temporalpole=UKB$rh_temporalpole_thickness/10)


select(theUKBdata, contains("_vol_")

theUKBdata <- merge(theUKBdata, theUKBCogdata[,1:12], by = "ID")
theUKBdata$g <- as.numeric(GpredictUKB$g)*-1 # so that a higher g means higher general cognitive functioning scores

UKB_g_regional_results <- matrix(NA, 68*3, 4)
colnames(UKB_g_regional_results) <- c("region", "beta", "se", "p")
index = 0
for (i in (1+12):((68*3)+12)) {
index = index + 1
loop_UKBmodeldata <- data.frame(theUKBdata[,1:12])
loop_UKBmodeldata$site1 <- 0
loop_UKBmodeldata$site1[which(theUKBdata$site==1)] = 1
loop_UKBmodeldata$site2 <- 0
loop_UKBmodeldata$site3[which(theUKBdata$site==2)] = 1
loop_UKBmodeldata$site3 <- 0
loop_UKBmodeldata$site3[which(theUKBdata$site==3)] = 1
loop_UKBmodeldata$site4 <- 0
loop_UKBmodeldata$site4[which(theUKBdata$site==4)] = 1
length(which(loop_UKBmodeldata$site1 == 1))
length(which(loop_UKBmodeldata$site2 == 1))
length(which(loop_UKBmodeldata$site3 == 1))
length(which(loop_UKBmodeldata$site4 == 1))
loop_UKBmodeldata$g <- theUKBdata$g
loop_UKBmodeldata <- cbind(loop_UKBmodeldata, theUKBdata[,i])
colnames(loop_UKBmodeldata)[dim(loop_UKBmodeldata)[2]] <- "loop_region"
regionname <- colnames(theUKBdata)[i]
loop_model <- 'loop_region ~ g + ageMRI + sex + site1 + site3 + site4 + X + Y + Z'
loop_fit <- sem(loop_model, data=loop_UKBmodeldata, missing = "fiml.x")
UKB_g_regional_results[index,1] <- regionname
UKB_g_regional_results[index,2] <- standardizedSolution(loop_fit)[1,c(4)] 
UKB_g_regional_results[index,3] <- standardizedSolution(loop_fit)[1,c(5)]
UKB_g_regional_results[index,4] <- standardizedSolution(loop_fit)[1,c(7)]
}

UKB_g_regional_results <- as.data.frame(UKB_g_regional_results)
UKB_g_regional_results[,2:4] <- sapply(UKB_g_regional_results[,2:4], as.numeric)
write.table(UKB_g_regional_results, 'UKB_g_regional_results.csv', row.names = F)
UKB_g_regional_results_vol <- UKB_g_regional_results[which(substr(UKB_g_regional_results$region, 4, 6) == "vol"),]
UKB_g_regional_results_sa <- UKB_g_regional_results[which(substr(UKB_g_regional_results$region, 4, 5) == "sa"),]
UKB_g_regional_results_thk <- UKB_g_regional_results[which(substr(UKB_g_regional_results$region, 4, 6) == "thk"),]


######## STRADL
STRADL <- read.csv("") # this file contains demographic information, the Desikan-Killiany regional volume, surface area, and thickness, and cognitive test scores for each STRADL participant.  
theSTRADLCogdata = data.frame(eid = STRADL$ID, age = STRADL$AgeFaceToFace, sex = STRADL$Sex, assCtr = STRADL$StudySite, mema = STRADL$mema, memdela = STRADL$memdela, digsym = STRADL$digsym, vftot = STRADL$vftot, mhv = STRADL$mhv, mrtotc = STRADL$mrtotc, logmem = STRADL$mema+STRADL$memdela)
cogmodel <- 'g =~ digsym + vftot + mhv + mrtotc + logmem 
digsym~age+sex
vftot ~age+sex
mhv~age+sex
mrtotc~age+sex
logmem~age+sex
'

cogfit <- sem(cogmodel, data = theSTRADLCogdata, missing = "fiml.x")
summary(cogfit, fit.measures=TRUE, standardized=T)

GpredictSTRADL <- lavPredict(cogfit, theSTRADLCogdata)
GpredictSTRADL <- cbind(STRADL$ID, GpredictSTRADL)
colnames(GpredictSTRADL) <- c("eid", "g")
GpredictSTRADL <- as.data.frame(GpredictSTRADL)
GpredictSTRADL$g <- as.numeric(GpredictSTRADL$g)
GpredictSTRADL$g <- scale(GpredictSTRADL$g)

#Organizing STRADL data for lavaan models  
theSTRADLdata=data.frame(ID=STRADL$ID,ageMRI=STRADL$AgeFaceToFace,dontuse=STRADL$AgeFaceToFace,lag=rep(0,dim(STRADL)[1]),sex=STRADL$Sex,atrophy=STRADL$atrophy,cohort=rep("STRADL",dim(STRADL)[1]),freesurfer=rep("unknown",dim(STRADL)[1]),site=STRADL$StudySite,X=STRADL$headposX,Y=STRADL$headposY,Z=STRADL$headposZ/10,lh_vol_bankssts=STRADL$lh_bankssts_volume/1000,lh_vol_caudalanteriorcingulate=STRADL$lh_caudalanteriorcingulate_volume/1000,lh_vol_caudalmiddlefrontal=STRADL$lh_caudalmiddlefrontal_volume/1000,lh_vol_cuneus=STRADL$lh_cuneus_volume/1000,lh_vol_entorhinal=STRADL$lh_entorhinal_volume/1000,lh_vol_fusiform=STRADL$lh_fusiform_volume/1000,lh_vol_inferiorparietal=STRADL$lh_inferiorparietal_volume/1000,lh_vol_inferiortemporal=STRADL$lh_inferiortemporal_volume/1000,lh_vol_isthmus=STRADL$lh_isthmuscingulate_volume/1000,lh_vol_lateraloccipital=STRADL$lh_lateraloccipital_volume/1000,lh_vol_lateralorbitofrontal=STRADL$lh_lateralorbitofrontal_volume/1000,lh_vol_lingual=STRADL$lh_lingual_volume/1000,lh_vol_medialorbitofrontal=STRADL$lh_medialorbitofrontal_volume/1000,lh_vol_middletemporal=STRADL$lh_middletemporal_volume/1000,lh_vol_parahippocampal=STRADL$lh_parahippocampal_volume/1000,lh_vol_paracentral=STRADL$lh_paracentral_volume/1000,lh_vol_parsopercularis=STRADL$lh_parsopercularis_volume/1000,lh_vol_parsorbitalis=STRADL$lh_parsorbitalis_volume/1000,lh_vol_parstriangularis=STRADL$lh_parstriangularis_volume/1000,lh_vol_pericalcarine=STRADL$lh_pericalcarine_volume/1000,lh_vol_postcentral=STRADL$lh_postcentral_volume/1000,lh_vol_posteriorcingulate=STRADL$lh_posteriorcingulate_volume/1000,lh_vol_precentral=STRADL$lh_precentral_volume/1000,lh_vol_precuneus=STRADL$lh_precuneus_volume/1000,lh_vol_rostralanteriorcingulate=STRADL$lh_rostralanteriorcingulate_volume/1000,lh_vol_rostralmiddlefrontal=STRADL$lh_rostralmiddlefrontal_volume/1000,lh_vol_superiorfrontal=STRADL$lh_superiorfrontal_volume/1000,lh_vol_superiorparietal=STRADL$lh_superiorparietal_volume/1000,lh_vol_superiortemporal=STRADL$lh_superiortemporal_volume/1000,lh_vol_supramarginal=STRADL$lh_supramarginal_volume/1000,lh_vol_frontalpole=STRADL$lh_frontalpole_volume/1000,lh_vol_transversetemporal=STRADL$lh_transversetemporal_volume/1000,lh_vol_insula=STRADL$lh_insula_volume/1000,lh_sa_bankssts=STRADL$lh_bankssts_area/100,lh_sa_caudalanteriorcingulate=STRADL$lh_caudalanteriorcingulate_area/100,lh_sa_caudalmiddlefrontal=STRADL$lh_caudalmiddlefrontal_area/100,lh_sa_cuneus=STRADL$lh_cuneus_area/100,lh_sa_entorhinal=STRADL$lh_entorhinal_area/100,lh_sa_fusiform=STRADL$lh_fusiform_area/100,lh_sa_inferiorparietal=STRADL$lh_inferiorparietal_area/100,lh_sa_inferiortemporal=STRADL$lh_inferiortemporal_area/100,lh_sa_isthmus=STRADL$lh_isthmuscingulate_area/100,lh_sa_lateraloccipital=STRADL$lh_lateraloccipital_area/100,lh_sa_lateralorbitofrontal=STRADL$lh_lateralorbitofrontal_area/100,lh_sa_lingual=STRADL$lh_lingual_area/100,lh_sa_medialorbitofrontal=STRADL$lh_medialorbitofrontal_area/100,lh_sa_middletemporal=STRADL$lh_middletemporal_area/100,lh_sa_parahippocampal=STRADL$lh_parahippocampal_area/100,lh_sa_paracentral=STRADL$lh_paracentral_area/100,lh_sa_parsopercularis=STRADL$lh_parsopercularis_area/100,lh_sa_parsorbitalis=STRADL$lh_parsorbitalis_area/100,lh_sa_parstriangularis=STRADL$lh_parstriangularis_area/100,lh_sa_pericalcarine=STRADL$lh_pericalcarine_area/100,lh_sa_postcentral=STRADL$lh_postcentral_area/100,lh_sa_posteriorcingulate=STRADL$lh_posteriorcingulate_area/100,lh_sa_precentral=STRADL$lh_precentral_area/100,lh_sa_precuneus=STRADL$lh_precuneus_area/100,lh_sa_rostralanteriorcingulate=STRADL$lh_rostralanteriorcingulate_area/100,lh_sa_rostralmiddlefrontal=STRADL$lh_rostralmiddlefrontal_area/100,lh_sa_superiorfrontal=STRADL$lh_superiorfrontal_area/100,lh_sa_superiorparietal=STRADL$lh_superiorparietal_area/100,lh_sa_superiortemporal=STRADL$lh_superiortemporal_area/100,lh_sa_supramarginal=STRADL$lh_supramarginal_area/100,lh_sa_frontalpole=STRADL$lh_frontalpole_area/100,lh_sa_transversetemporal=STRADL$lh_transversetemporal_area/100,lh_sa_insula=STRADL$lh_insula_area/100,lh_thk_bankssts=STRADL$lh_bankssts_thickness*10,lh_thk_caudalanteriorcingulate=STRADL$lh_caudalanteriorcingulate_thickness*10,lh_thk_caudalmiddlefrontal=STRADL$lh_caudalmiddlefrontal_thickness*10,lh_thk_cuneus=STRADL$lh_cuneus_thickness*10,lh_thk_entorhinal=STRADL$lh_entorhinal_thickness*10,lh_thk_fusiform=STRADL$lh_fusiform_thickness*10,lh_thk_inferiorparietal=STRADL$lh_inferiorparietal_thickness*10,lh_thk_inferiortemporal=STRADL$lh_inferiortemporal_thickness*10,lh_thk_isthmus=STRADL$lh_isthmuscingulate_thickness*10,lh_thk_lateraloccipital=STRADL$lh_lateraloccipital_thickness*10,lh_thk_lateralorbitofrontal=STRADL$lh_lateralorbitofrontal_thickness*10,lh_thk_lingual=STRADL$lh_lingual_thickness*10,lh_thk_medialorbitofrontal=STRADL$lh_medialorbitofrontal_thickness*10,lh_thk_middletemporal=STRADL$lh_middletemporal_thickness*10,lh_thk_parahippocampal=STRADL$lh_parahippocampal_thickness*10,lh_thk_paracentral=STRADL$lh_paracentral_thickness*10,lh_thk_parsopercularis=STRADL$lh_parsopercularis_thickness*10,lh_thk_parsorbitalis=STRADL$lh_parsorbitalis_thickness*10,lh_thk_parstriangularis=STRADL$lh_parstriangularis_thickness*10,lh_thk_pericalcarine=STRADL$lh_pericalcarine_thickness*10,lh_thk_postcentral=STRADL$lh_postcentral_thickness*10,lh_thk_posteriorcingulate=STRADL$lh_posteriorcingulate_thickness*10,lh_thk_precentral=STRADL$lh_precentral_thickness*10,lh_thk_precuneus=STRADL$lh_precuneus_thickness*10,lh_thk_rostralanteriorcingulate=STRADL$lh_rostralanteriorcingulate_thickness*10,lh_thk_rostralmiddlefrontal=STRADL$lh_rostralmiddlefrontal_thickness*10,lh_thk_superiorfrontal=STRADL$lh_superiorfrontal_thickness*10,lh_thk_superiorparietal=STRADL$lh_superiorparietal_thickness*10,lh_thk_superiortemporal=STRADL$lh_superiortemporal_thickness*10,lh_thk_supramarginal=STRAD$lh_supramarginal_thickness*10,lh_thk_frontalpole=STRADL$lh_frontalpole_thickness*10,lh_thk_transversetemporal=STRADL$lh_transversetemporal_thickness*10,lh_thk_insula=STRADL$lh_insula_thickness*10,rh_vol_bankssts=STRADL$rh_bankssts_volume/1000,rh_vol_caudalanteriorcingulate=STRADL$rh_caudalanteriorcingulate_volume/1000,rh_vol_caudalmiddlefrontal=STRADL$rh_caudalmiddlefrontal_volume/1000,rh_vol_cuneus=STRADL$rh_cuneus_volume/1000,rh_vol_entorhinal=STRADL$rh_entorhinal_volume/1000,rh_vol_fusiform=STRADL$rh_fusiform_volume/1000,rh_vol_inferiorparietal=STRADL$rh_inferiorparietal_volume/1000,rh_vol_inferiortemporal=STRADL$rh_inferiortemporal_volume/1000,rh_vol_isthmus=STRADL$rh_isthmuscingulate_volume/1000,rh_vol_lateraloccipital=STRADL$rh_lateraloccipital_volume/1000,rh_vol_lateralorbitofrontal=STRADL$rh_lateralorbitofrontal_volume/1000,rh_vol_lingual=STRADL$rh_lingual_volume/1000,rh_vol_medialorbitofrontal=STRADL$rh_medialorbitofrontal_volume/1000,rh_vol_middletemporal=STRADL$rh_middletemporal_volume/1000,rh_vol_parahippocampal=STRADL$rh_parahippocampal_volume/1000,rh_vol_paracentral=STRADL$rh_paracentral_volume/1000,rh_vol_parsopercularis=STRADL$rh_parsopercularis_volume/1000,rh_vol_parsorbitalis=STRADL$rh_parsorbitalis_volume/1000,rh_vol_parstriangularis=STRADL$rh_parstriangularis_volume/1000,rh_vol_pericalcarine=STRADL$rh_pericalcarine_volume/1000,rh_vol_postcentral=STRADL$rh_postcentral_volume/1000,rh_vol_posteriorcingulate=STRADL$rh_posteriorcingulate_volume/1000,rh_vol_precentral=STRADL$rh_precentral_volume/1000,rh_vol_precuneus=STRADL$rh_precuneus_volume/1000,rh_vol_rostralanteriorcingulate=STRADL$rh_rostralanteriorcingulate_volume/1000,rh_vol_rostralmiddlefrontal=STRADL$rh_rostralmiddlefrontal_volume/1000,rh_vol_superiorfrontal=STRADL$rh_superiorfrontal_volume/1000,rh_vol_superiorparietal=STRADL$rh_superiorparietal_volume/1000,rh_vol_superiortemporal=STRADL$rh_superiortemporal_volume/1000,rh_vol_supramarginal=STRADL$rh_supramarginal_volume/1000,rh_vol_frontalpole=STRADL$rh_frontalpole_volume/1000,rh_vol_transversetemporal=STRADL$rh_transversetemporal_volume/1000,rh_vol_insula=STRADL$rh_insula_volume/1000,rh_sa_bankssts=STRADL$rh_bankssts_area/100,rh_sa_caudalanteriorcingulate=STRADL$rh_caudalanteriorcingulate_area/100,rh_sa_caudalmiddlefrontal=STRADL$rh_caudalmiddlefrontal_area/100,rh_sa_cuneus=STRADL$rh_cuneus_area/100,rh_sa_entorhinal=STRADL$rh_entorhinal_area/100,rh_sa_fusiform=STRADL$rh_fusiform_area/100,rh_sa_inferiorparietal=STRADL$rh_inferiorparietal_area/100,rh_sa_inferiortemporal=STRADL$rh_inferiortemporal_area/100,rh_sa_isthmus=STRADL$rh_isthmuscingulate_area/100,rh_sa_lateraloccipital=STRADL$rh_lateraloccipital_area/100,rh_sa_lateralorbitofrontal=STRADL$rh_lateralorbitofrontal_area/100,rh_sa_lingual=STRADL$rh_lingual_area/100,rh_sa_medialorbitofrontal=STRADL$rh_medialorbitofrontal_area/100,rh_sa_middletemporal=STRADL$rh_middletemporal_area/100,rh_sa_parahippocampal=STRADL$rh_parahippocampal_area/100,rh_sa_paracentral=STRADL$rh_paracentral_area/100,rh_sa_parsopercularis=STRADL$rh_parsopercularis_area/100,rh_sa_parsorbitalis=STRADL$rh_parsorbitalis_area/100,rh_sa_parstriangularis=STRADL$rh_parstriangularis_area/100,rh_sa_pericalcarine=STRADL$rh_pericalcarine_area/100,rh_sa_postcentral=STRADL$rh_postcentral_area/100,rh_sa_posteriorcingulate=STRADL$rh_posteriorcingulate_area/100,rh_sa_precentral=STRADL$rh_precentral_area/100,rh_sa_precuneus=STRADL$rh_precuneus_area/100,rh_sa_rostralanteriorcingulate=STRADL$rh_rostralanteriorcingulate_area/100,rh_sa_rostralmiddlefrontal=STRADL$rh_rostralmiddlefrontal_area/100,rh_sa_superiorfrontal=STRADL$rh_superiorfrontal_area/100,rh_sa_superiorparietal=STRADL$rh_superiorparietal_area/100,rh_sa_superiortemporal=STRADL$rh_superiortemporal_area/100,rh_sa_supramarginal=STRADL$rh_supramarginal_area/100,rh_sa_frontalpole=STRADL$rh_frontalpole_area/100,rh_sa_transversetemporal=STRADL$rh_transversetemporal_area/100,rh_sa_insula=STRADL$rh_insula_area/100,rh_thk_bankssts=STRADL$rh_bankssts_thickness*10,rh_thk_caudalanteriorcingulate=STRADL$rh_caudalanteriorcingulate_thickness*10,rh_thk_caudalmiddlefrontal=STRADL$rh_caudalmiddlefrontal_thickness*10,rh_thk_cuneus=STRADL$rh_cuneus_thickness*10,rh_thk_entorhinal=STRADL$rh_entorhinal_thickness*10,rh_thk_fusiform=STRADL$rh_fusiform_thickness*10,rh_thk_inferiorparietal=STRADL$rh_inferiorparietal_thickness*10,rh_thk_inferiortemporal=STRADL$rh_inferiortemporal_thickness*10,rh_thk_isthmus=STRADL$rh_isthmuscingulate_thickness*10,rh_thk_lateraloccipital=STRADL$rh_lateraloccipital_thickness*10,rh_thk_lateralorbitofrontal=STRADL$rh_lateralorbitofrontal_thickness*10,rh_thk_lingual=STRADL$rh_lingual_thickness*10,rh_thk_medialorbitofrontal=STRADL$rh_medialorbitofrontal_thickness*10,rh_thk_middletemporal=STRADL$rh_middletemporal_thickness*10,rh_thk_parahippocampal=STRADL$rh_parahippocampal_thickness*10,rh_thk_paracentral=STRADL$rh_paracentral_thickness*10,rh_thk_parsopercularis=STRADL$rh_parsopercularis_thickness*10,rh_thk_parsorbitalis=STRADL$rh_parsorbitalis_thickness*10,rh_thk_parstriangularis=STRADL$rh_parstriangularis_thickness*10,rh_thk_pericalcarine=STRADL$rh_pericalcarine_thickness*10,rh_thk_postcentral=STRADL$rh_postcentral_thickness*10,rh_thk_posteriorcingulate=STRADL$rh_posteriorcingulate_thickness*10,rh_thk_precentral=STRADL$rh_precentral_thickness*10,rh_thk_precuneus=STRADL$rh_precuneus_thickness*10,rh_thk_rostralanteriorcingulate=STRADL$rh_rostralanteriorcingulate_thickness*10,rh_thk_rostralmiddlefrontal=STRADL$rh_rostralmiddlefrontal_thickness*10,rh_thk_superiorfrontal=STRADL$rh_superiorfrontal_thickness*10,rh_thk_superiorparietal=STRADL$rh_superiorparietal_thickness*10,rh_thk_superiortemporal=STRADL$rh_superiortemporal_thickness*10,rh_thk_supramarginal=STRADL$rh_supramarginal_thickness*10,rh_thk_frontalpole=STRADL$rh_frontalpole_thickness*10,rh_thk_transversetemporal=STRADL$rh_transversetemporal_thickness*10,rh_thk_insula=STRADL$rh_insula_thickness*10,lh_vol_temporalpole=STRADL$lh_temporalpole_volume/1000,rh_vol_temporalpole=STRADL$rh_temporalpole_volume/1000,lh_sa_temporalpole=STRADL$lh_temporalpole_area/100,rh_sa_temporalpole=STRADL$rh_temporalpole_area/100,lh_thk_temporalpole=STRADL$lh_temporalpole_thickness/10,rh_thk_temporalpole=STRADL$rh_temporalpole_thickness/10)

theSTRADLdata$Gpredict <- as.numeric(GpredictSTRADL$g)
colnames(theSTRADLdata)[217] <- "g"

       
STRADL_g_regional_results <- matrix(NA, 68*3, 4)
colnames(STRADL_g_regional_results) <- c("region", "beta", "se", "p")


index = 0
for (i in (1+12):((68*3)+12)) {
index = index + 1
loop_STRADLmodeldata <- data.frame(theSTRADLdata[,1:12])
loop_STRADLmodeldata$g <- theSTRADLdata$g
loop_STRADLmodeldata$ageMRIsquare <- loop_STRADLmodeldata$ageMRI^2 # including the square causes the model to crash because it's so highly correlated with ageMRI, so I've left this out for now
loop_STRADLmodeldata <- cbind(loop_STRADLmodeldata, theSTRADLdata[,i])
colnames(loop_STRADLmodeldata)[dim(loop_STRADLmodeldata)[2]] <- "loop_region"
regionname <- colnames(theSTRADLdata)[i]

loop_model <- 'loop_region ~ g + ageMRI + sex + site + X + Y + Z'
loop_fit <- sem(loop_model, data=loop_STRADLmodeldata, missing = "fiml.x")
STRADL_g_regional_results[index,1] <- regionname
STRADL_g_regional_results[index,2] <- standardizedSolution(loop_fit)[1,c(4)] 
STRADL_g_regional_results[index,3] <- standardizedSolution(loop_fit)[1,c(5)]
STRADL_g_regional_results[index,4] <- standardizedSolution(loop_fit)[1,c(7)]



}

STRADL_g_regional_results <- as.data.frame(STRADL_g_regional_results)
STRADL_g_regional_results[,2:4] <- sapply(STRADL_g_regional_results[,2:4], as.numeric)
write.table(STRADL_g_regional_results, 'STRADL_g_regional_results.csv', row.names = F)
STRADL_g_regional_results_vol <- STRADL_g_regional_results[which(substr(STRADL_g_regional_results$region, 4, 6) == "vol"),]
STRADL_g_regional_results_sa <- STRADL_g_regional_results[which(substr(STRADL_g_regional_results$region, 4, 5) == "sa"),]
STRADL_g_regional_results_thk <- STRADL_g_regional_results[which(substr(STRADL_g_regional_results$region, 4, 6) == "thk"),]

allbetas <- cbind(STRADL_g_regional_results_vol[,2], STRADL_g_regional_results_sa[,2], STRADL_g_regional_results_thk[,2], UKB_g_regional_results_vol[,2], UKB_g_regional_results_sa[,2], UKB_g_regional_results_thk[,2])
colnames(allbetas) <- c("vol_STRADL", "sa_STRADL", "thk_STRADL", "vol_UKB", "sa_UKB", "thk_UKB")
cor(allbetas)

######## LBC

LBC <- read.csv("") # this file contains demographic information, the Desikan-Killiany regional volume, surface area, and thickness, and cognitive test scores for each LBC1936 participant.
theLBCCogdata <- LBC[,1:25]
model2 <- 'g=~matrix_reasoning + block_design + spatial_span_total + NART + WTAR + verbal_fluency+verbal_paired_associates + logical_memory + digit_span_backward + symbol_search + digit_symbol + inspection_time + choice_reaction_time_reflected 

matrix_reasoning ~ AgeDays + sex
block_design~ AgeDays + sex
spatial_span_total~ AgeDays + sex
NART~ AgeDays + sex
WTAR ~ AgeDays + sex
verbal_fluency~ AgeDays + sex
verbal_paired_associates~ AgeDays + sex
logical_memory~ AgeDays + sex
digit_span_backward~ AgeDays + sex
symbol_search~ AgeDays + sex
digit_symbol~ AgeDays + sex
inspection_time~ AgeDays + sex
choice_reaction_time_reflected~ AgeDays + sex

# within-domain covariances
logical_memory~~verbal_paired_associates 
matrix_reasoning ~~ block_design
symbol_search ~~ digit_symbol
matrix_reasoning~~ spatial_span_total
block_design~~spatial_span_total
NART~~WTAR
NART~~verbal_fluency
WTAR~~verbal_fluency
verbal_paired_associates ~~ digit_span_backward
logical_memory ~~ digit_span_backward
symbol_search~~inspection_time
symbol_search~~choice_reaction_time_reflected
digit_symbol~~choice_reaction_time_reflected
inspection_time~~choice_reaction_time_reflected
digit_symbol~~inspection_time
'

fit1 <- sem(model2, data=theLBCCogdata, missing = "fiml.x")
summary(fit1, fit.measures=TRUE, standardized=T)
modificationIndices(fit1, sort = T)

GpredictLBC <- lavPredict(fit1, theLBCCogdata)
GpredictLBC <- as.data.frame(GpredictLBC)
GpredictLBC <- cbind(LBC$ID, GpredictLBC$g)
colnames(GpredictLBC) <- c("eid", "g")
GpredictLBC <- as.data.frame(GpredictLBC)
GpredictLBC$g <- as.numeric(GpredictLBC$g)
GpredictLBC$g <- scale(GpredictLBC$g)

#Organizing LBC data for lavaan models
theLBCdata = data.frame(ID = LBC$ID, ageMRI = LBC$AgeMRI/365.24, dontuse = LBC$AgeDays/365.24, lag = (LBC$AgeMRI-LBC$AgeDays), sex = LBC$sex, atrophy = LBC$atrophy, cohort = rep("LBC1936", dim(LBC)[1]), freesurfer = rep("NA", dim(LBC)[1]), site = rep("LBCsite", dim(LBC)[1]), X = LBC$headposX, Y = LBC$headposY, Z = LBC$headposZ/10, lh_vol_bankssts = LBC$lh_vol_bankssts/1000, lh_vol_caudalanteriorcingulate = LBC$lh_vol_caudalanteriorcingulate/1000, lh_vol_caudalmiddlefrontal = LBC$lh_vol_caudalmiddlefrontal/1000, lh_vol_cuneus = LBC$lh_vol_cuneus/1000, lh_vol_entorhinal = LBC$lh_vol_entorhinal/1000, lh_vol_fusiform = LBC$lh_vol_fusiform/1000, lh_vol_inferiorparietal = LBC$lh_vol_inferiorparietal/1000, lh_vol_inferiortemporal = LBC$lh_vol_inferiortemporal/1000, lh_vol_isthmus = LBC$lh_vol_isthmus/1000, lh_vol_lateraloccipital = LBC$lh_vol_lateraloccipital/1000, lh_vol_lateralorbitofrontal = LBC$lh_vol_lateralorbitofrontal/1000, lh_vol_lingual = LBC$lh_vol_lingual/1000, lh_vol_medialorbitofrontal = LBC$lh_vol_medialorbitofrontal/1000, lh_vol_middletemporal = LBC$lh_vol_middletemporal/1000, lh_vol_parahippocampal= LBC$lh_vol_parahippocampal/1000, lh_vol_paracentral= LBC$lh_vol_paracentral/1000, lh_vol_parsopercularis = LBC$lh_vol_parsopercularis/1000, lh_vol_parsorbitalis = LBC$lh_vol_parsorbitalis/1000, lh_vol_parstriangularis = LBC$lh_vol_parstriangularis/1000, lh_vol_pericalcarine = LBC$lh_vol_pericalcarine/1000, lh_vol_postcentral = LBC$lh_vol_postcentral/1000, lh_vol_posteriorcingulate= LBC$lh_vol_posteriorcingulate/1000, lh_vol_precentral = LBC$lh_vol_precentral/1000, lh_vol_precuneus= LBC$lh_vol_precuneus/1000, lh_vol_rostralanteriorcingulate = LBC$lh_vol_rostralanteriorcingulate/1000, lh_vol_rostralmiddlefrontal = LBC$lh_vol_rostralmiddlefrontal/1000,lh_vol_superiorfrontal = LBC$lh_vol_superiorfrontal/1000, lh_vol_superiorparietal = LBC$lh_vol_superiorparietal/1000, lh_vol_superiortemporal = LBC$lh_vol_superiortemporal/1000, lh_vol_supramarginal = LBC$lh_vol_supramarginal/1000, lh_vol_frontalpole = LBC$lh_vol_frontalpole/1000, lh_vol_transversetemporal = LBC$lh_vol_transversetemporal/1000, lh_vol_insula = LBC$lh_vol_insula/1000,lh_sa_bankssts = LBC$lh_sa_bankssts/100, lh_sa_caudalanteriorcingulate = LBC$lh_sa_caudalanteriorcingulate/100, lh_sa_caudalmiddlefrontal = LBC$lh_sa_caudalmiddlefrontal/100, lh_sa_cuneus = LBC$lh_sa_cuneus/100, lh_sa_entorhinal = LBC$lh_sa_entorhinal/100, lh_sa_fusiform= LBC$lh_sa_fusiform/100, lh_sa_inferiorparietal = LBC$lh_sa_inferiorparietal/100, lh_sa_inferiortemporal = LBC$lh_sa_inferiortemporal/100, lh_sa_isthmus = LBC$lh_sa_isthmus/100, lh_sa_lateraloccipital = LBC$lh_sa_lateraloccipital/100,  lh_sa_lateralorbitofrontal = LBC$lh_sa_lateralorbitofrontal/100, lh_sa_lingual =LBC$lh_sa_lingual/100, lh_sa_medialorbitofrontal = LBC$lh_sa_medialorbitofrontal/100, lh_sa_middletemporal = LBC$lh_sa_middletemporal/100, lh_sa_parahippocampal = LBC$lh_sa_parahippocampal/100, lh_sa_paracentral = LBC$lh_sa_paracentral/100, lh_sa_parsopercularis = LBC$lh_sa_parsopercularis/100,lh_sa_parsorbitalis = LBC$lh_sa_parsorbitalis/100, lh_sa_parstriangularis = LBC$lh_sa_parstriangularis/100, lh_sa_pericalcarine = LBC$lh_sa_pericalcarine/100, lh_sa_postcentral = LBC$lh_sa_postcentral/100, lh_sa_posteriorcingulate= LBC$lh_sa_posteriorcingulate/100, lh_sa_precentral = LBC$lh_sa_precentral/100,lh_sa_precuneus = LBC$lh_sa_precuneus/100, lh_sa_rostralanteriorcingulate = LBC$lh_sa_rostralanteriorcingulate/100, lh_sa_rostralmiddlefrontal = LBC$lh_sa_rostralmiddlefrontal/100, lh_sa_superiorfrontal = LBC$lh_sa_superiorfrontal/100, lh_sa_superiorparietal = LBC$lh_sa_superiorparietal/100, lh_sa_superiortemporal = LBC$lh_sa_superiortemporal/100, lh_sa_supramarginal = LBC$lh_sa_supramarginal/100, lh_sa_frontalpole = LBC$lh_sa_frontalpole/100, lh_sa_transversetemporal = LBC$lh_sa_transversetemporal/100, lh_sa_insula = LBC$lh_sa_insula/100, lh_thk_bankssts = LBC$lh_thk_bankssts*10, lh_thk_caudalanteriorcingulate = LBC$lh_thk_caudalanteriorcingulate*10, lh_thk_caudalmiddlefrontal = LBC$lh_thk_caudalmiddlefrontal*10, lh_thk_cuneus = LBC$lh_thk_cuneus*10, lh_thk_entorhinal = LBC$lh_thk_entorhinal*10, lh_thk_fusiform = LBC$lh_thk_fusiform*10, lh_thk_inferiorparietal = LBC$lh_thk_inferiorparietal*10, lh_thk_inferiortemporal = LBC$lh_thk_inferiortemporal*10, lh_thk_isthmus = LBC$lh_thk_isthmus*10, lh_thk_lateraloccipital = LBC$lh_thk_lateraloccipital*10, lh_thk_lateralorbitofrontal = LBC$lh_thk_lateralorbitofrontal*10,lh_thk_lingual = LBC$lh_thk_lingual*10, lh_thk_medialorbitofrontal = LBC$lh_thk_medialorbitofrontal*10, lh_thk_middletemporal = LBC$lh_thk_middletemporal*10, lh_thk_parahippocampal = LBC$lh_thk_parahippocampal*10, lh_thk_paracentral = LBC$lh_thk_paracentral*10, lh_thk_parsopercularis = LBC$lh_thk_parsopercularis*10, lh_thk_parsorbitalis = LBC$lh_thk_parsorbitalis*10, lh_thk_parstriangularis = LBC$lh_thk_parstriangularis*10, lh_thk_pericalcarine = LBC$lh_thk_pericalcarine*10,  lh_thk_postcentral = LBC$lh_thk_postcentral*10, lh_thk_posteriorcingulate = LBC$lh_thk_posteriorcingulate*10, lh_thk_precentral = LBC$lh_thk_precentral*10,lh_thk_precuneus = LBC$lh_thk_precuneus*10, lh_thk_rostralanteriorcingulate = LBC$lh_thk_rostralanteriorcingulate*10, lh_thk_rostralmiddlefrontal = LBC$lh_thk_rostralmiddlefrontal*10, lh_thk_superiorfrontal = LBC$lh_thk_superiorfrontal*10, lh_thk_superiorparietal = LBC$lh_thk_superiorparietal*10,  lh_thk_superiortemporal = LBC$lh_thk_superiortemporal*10, lh_thk_supramarginal = LBC$lh_thk_supramarginal*10, lh_thk_frontalpole = LBC$lh_thk_frontalpole*10, lh_thk_transversetemporal = LBC$lh_thk_transversetemporal*10, lh_thk_insula = LBC$lh_thk_insula*10, rh_vol_bankssts = LBC$rh_vol_bankssts/1000, rh_vol_caudalanteriorcingulate = LBC$rh_vol_caudalanteriorcingulate/1000, rh_vol_caudalmiddlefrontal = LBC$rh_vol_caudalmiddlefrontal/1000, rh_vol_cuneus = LBC$rh_vol_cuneus/1000, rh_vol_entorhinal = LBC$rh_vol_entorhinal/1000, rh_vol_fusiform = LBC$rh_vol_fusiform/1000, rh_vol_inferiorparietal = LBC$rh_vol_inferiorparietal/1000, rh_vol_inferiortemporal = LBC$rh_vol_inferiortemporal/1000, rh_vol_isthmus = LBC$rh_vol_isthmus/1000, rh_vol_lateraloccipital = LBC$rh_vol_lateraloccipital/1000, rh_vol_lateralorbitofrontal = LBC$rh_vol_lateralorbitofrontal/1000, rh_vol_lingual = LBC$rh_vol_lingual/1000, rh_vol_medialorbitofrontal = LBC$rh_vol_medialorbitofrontal/1000, rh_vol_middletemporal = LBC$rh_vol_middletemporal/1000, rh_vol_parahippocampal= LBC$rh_vol_parahippocampal/1000, rh_vol_paracentral= LBC$rh_vol_paracentral/1000, rh_vol_parsopercularis = LBC$rh_vol_parsopercularis/1000, rh_vol_parsorbitalis = LBC$rh_vol_parsorbitalis/1000, rh_vol_parstriangularis = LBC$rh_vol_parstriangularis/1000, rh_vol_pericalcarine = LBC$rh_vol_pericalcarine/1000, rh_vol_postcentral = LBC$rh_vol_postcentral/1000, rh_vol_posteriorcingulate= LBC$rh_vol_posteriorcingulate/1000, rh_vol_precentral = LBC$rh_vol_precentral/1000, rh_vol_precuneus= LBC$rh_vol_precuneus/1000, rh_vol_rostralanteriorcingulate = LBC$rh_vol_rostralanteriorcingulate/1000, rh_vol_rostralmiddlefrontal = LBC$rh_vol_rostralmiddlefrontal/1000,rh_vol_superiorfrontal = LBC$rh_vol_superiorfrontal/1000, rh_vol_superiorparietal = LBC$rh_vol_superiorparietal/1000, rh_vol_superiortemporal = LBC$rh_vol_superiortemporal/1000, rh_vol_supramarginal = LBC$rh_vol_supramarginal/1000, rh_vol_frontalpole = LBC$rh_vol_frontalpole/1000, rh_vol_transversetemporal = LBC$rh_vol_transversetemporal/1000, rh_vol_insula = LBC$rh_vol_insula/1000,rh_sa_bankssts = LBC$rh_sa_bankssts/100, rh_sa_caudalanteriorcingulate = LBC$rh_sa_caudalanteriorcingulate/100, rh_sa_caudalmiddlefrontal = LBC$rh_sa_caudalmiddlefrontal/100, rh_sa_cuneus = LBC$rh_sa_cuneus/100, rh_sa_entorhinal = LBC$rh_sa_entorhinal/100, rh_sa_fusiform= LBC$rh_sa_fusiform/100, rh_sa_inferiorparietal = LBC$rh_sa_inferiorparietal/100, rh_sa_inferiortemporal = LBC$rh_sa_inferiortemporal/100, rh_sa_isthmus = LBC$rh_sa_isthmus/100, rh_sa_lateraloccipital = LBC$rh_sa_lateraloccipital/100,  rh_sa_lateralorbitofrontal = LBC$rh_sa_lateralorbitofrontal/100, rh_sa_lingual =LBC$rh_sa_lingual/100, rh_sa_medialorbitofrontal = LBC$rh_sa_medialorbitofrontal/100, rh_sa_middletemporal = LBC$rh_sa_middletemporal/100, rh_sa_parahippocampal = LBC$rh_sa_parahippocampal/100, rh_sa_paracentral = LBC$rh_sa_paracentral/100,  rh_sa_parsopercularis = LBC$rh_sa_parsopercularis/100,rh_sa_parsorbitalis = LBC$rh_sa_parsorbitalis/100, rh_sa_parstriangularis = LBC$rh_sa_parstriangularis/100, rh_sa_pericalcarine = LBC$rh_sa_pericalcarine/100, rh_sa_postcentral = LBC$rh_sa_postcentral/100, rh_sa_posteriorcingulate= LBC$rh_sa_posteriorcingulate/100, rh_sa_precentral = LBC$rh_sa_precentral/100,rh_sa_precuneus = LBC$rh_sa_precuneus/100, rh_sa_rostralanteriorcingulate = LBC$rh_sa_rostralanteriorcingulate/100, rh_sa_rostralmiddlefrontal = LBC$rh_sa_rostralmiddlefrontal/100, rh_sa_superiorfrontal = LBC$rh_sa_superiorfrontal/100, rh_sa_superiorparietal = LBC$rh_sa_superiorparietal/100, rh_sa_superiortemporal = LBC$rh_sa_superiortemporal/100, rh_sa_supramarginal = LBC$rh_sa_supramarginal/100, rh_sa_frontalpole = LBC$rh_sa_frontalpole/100, rh_sa_transversetemporal = LBC$rh_sa_transversetemporal/100, rh_sa_insula = LBC$rh_sa_insula/100, rh_thk_bankssts = LBC$rh_thk_bankssts*10, rh_thk_caudalanteriorcingulate = LBC$rh_thk_caudalanteriorcingulate*10, rh_thk_caudalmiddlefrontal = LBC$rh_thk_caudalmiddlefrontal*10, rh_thk_cuneus = LBC$rh_thk_cuneus*10, rh_thk_entorhinal = LBC$rh_thk_entorhinal*10, rh_thk_fusiform = LBC$rh_thk_fusiform*10, rh_thk_inferiorparietal = LBC$rh_thk_inferiorparietal*10, rh_thk_inferiortemporal = LBC$rh_thk_inferiortemporal*10,  rh_thk_isthmus = LBC$rh_thk_isthmus*10, rh_thk_lateraloccipital = LBC$rh_thk_lateraloccipital*10, rh_thk_lateralorbitofrontal = LBC$rh_thk_lateralorbitofrontal*10,rh_thk_lingual = LBC$rh_thk_lingual*10, rh_thk_medialorbitofrontal = LBC$rh_thk_medialorbitofrontal*10, rh_thk_middletemporal = LBC$rh_thk_middletemporal*10, rh_thk_parahippocampal = LBC$rh_thk_parahippocampal*10, rh_thk_paracentral = LBC$rh_thk_paracentral*10, rh_thk_parsopercularis = LBC$rh_thk_parsopercularis*10, rh_thk_parsorbitalis = LBC$rh_thk_parsorbitalis*10, rh_thk_parstriangularis = LBC$rh_thk_parstriangularis*10, rh_thk_pericalcarine = LBC$rh_thk_pericalcarine*10,  rh_thk_postcentral = LBC$rh_thk_postcentral*10, rh_thk_posteriorcingulate = LBC$rh_thk_posteriorcingulate*10, rh_thk_precentral = LBC$rh_thk_precentral*10,rh_thk_precuneus = LBC$rh_thk_precuneus*10, rh_thk_rostralanteriorcingulate = LBC$rh_thk_rostralanteriorcingulate*10, rh_thk_rostralmiddlefrontal = LBC$rh_thk_rostralmiddlefrontal*10, rh_thk_superiorfrontal = LBC$rh_thk_superiorfrontal*10, rh_thk_superiorparietal = LBC$rh_thk_superiorparietal*10, rh_thk_superiortemporal = LBC$rh_thk_superiortemporal*10, rh_thk_supramarginal = LBC$rh_thk_supramarginal*10, rh_thk_frontalpole = LBC$rh_thk_frontalpole*10, rh_thk_transversetemporal = LBC$rh_thk_transversetemporal*10, rh_thk_insula = LBC$rh_thk_insula*10, lh_vol_temporalpole = LBC$lh_vol_temporalpole/1000, rh_vol_temporalpole = LBC$rh_vol_temporalpole/1000,lh_sa_temporalpole = LBC$lh_sa_temporalpole/100, rh_sa_temporalpole = LBC$rh_sa_temporalpole/100,lh_thk_temporalpole = LBC$lh_thk_temporalpole/10,rh_thk_temporalpole = LBC$rh_thk_temporal/10)

theLBCdata$Gpredict <- as.numeric(GpredictLBC$g)
colnames(theLBCdata)[217] <- "g"


LBC_g_regional_results <- matrix(NA, 68*3, 4)
colnames(LBC_g_regional_results) <- c("region", "beta", "se", "p")
index = 0
for (i in (1+12):((68*3)+12)) {
index = index + 1
loop_LBCmodeldata <- data.frame(theLBCdata[,1:12])
loop_LBCmodeldata$g <- theLBCdata$g
loop_LBCmodeldata$ageMRIsquare <- loop_LBCmodeldata$ageMRI^2 # including the square causes the model to crash because it's so highly correlated with ageMRI, so I've left this out for now
loop_LBCmodeldata <- cbind(loop_LBCmodeldata, theLBCdata[,i])
colnames(loop_LBCmodeldata)[dim(loop_LBCmodeldata)[2]] <- "loop_region"
regionname <- colnames(theLBCdata)[i]

loop_model <- 'loop_region ~ g + ageMRI + sex + lag + X + Y + Z'
loop_fit <- sem(loop_model, data=loop_LBCmodeldata, missing = "fiml.x")
LBC_g_regional_results[index,1] <- regionname
LBC_g_regional_results[index,2] <- standardizedSolution(loop_fit)[1,c(4)] 
LBC_g_regional_results[index,3] <- standardizedSolution(loop_fit)[1,c(5)]
LBC_g_regional_results[index,4] <- standardizedSolution(loop_fit)[1,c(7)]

}

LBC_g_regional_results <- as.data.frame(LBC_g_regional_results)
LBC_g_regional_results[,2:4] <- sapply(LBC_g_regional_results[,2:4], as.numeric)
write.table(LBC_g_regional_results, 'LBC_g_regional_results.csv', row.names = F)
LBC_g_regional_results_vol <- LBC_g_regional_results[which(substr(LBC_g_regional_results$region, 4, 6) == "vol"),]
LBC_g_regional_results_sa <- LBC_g_regional_results[which(substr(LBC_g_regional_results$region, 4, 5) == "sa"),]
LBC_g_regional_results_thk <- LBC_g_regional_results[which(substr(LBC_g_regional_results$region, 4, 6) == "thk"),]

# Meta-analyses
UKB_g_regional_results <- read.csv('UKB_g_regional_results.csv', sep = " ")
STRADL_g_regional_results <- read.csv('STRADL_g_regional_results.csv', sep = " ")
LBC_g_regional_results <- read.csv('LBC_g_regional_results.csv', sep = " ")

UKB_g_regional_results_vol <- UKB_g_regional_results[which(substr(UKB_g_regional_results$region, 4, 6) == "vol"),]
UKB_g_regional_results_sa <- UKB_g_regional_results[which(substr(UKB_g_regional_results$region, 4, 5) == "sa"),]
UKB_g_regional_results_thk <- UKB_g_regional_results[which(substr(UKB_g_regional_results$region, 4, 6) == "thk"),]
STRADL_g_regional_results_vol <- STRADL_g_regional_results[which(substr(STRADL_g_regional_results$region, 4, 6) == "vol"),]
STRADL_g_regional_results_sa <- STRADL_g_regional_results[which(substr(STRADL_g_regional_results$region, 4, 5) == "sa"),]
STRADL_g_regional_results_thk <- STRADL_g_regional_results[which(substr(STRADL_g_regional_results$region, 4, 6) == "thk"),]
LBC_g_regional_results_vol <- LBC_g_regional_results[which(substr(LBC_g_regional_results$region, 4, 6) == "vol"),]
LBC_g_regional_results_sa <- LBC_g_regional_results[which(substr(LBC_g_regional_results$region, 4, 5) == "sa"),]
LBC_g_regional_results_thk <- LBC_g_regional_results[which(substr(LBC_g_regional_results$region, 4, 6) == "thk"),]
allbetas <- cbind(STRADL_g_regional_results_vol[,2], STRADL_g_regional_results_sa[,2], STRADL_g_regional_results_thk[,2], UKB_g_regional_results_vol[,2], UKB_g_regional_results_sa[,2], UKB_g_regional_results_thk[,2], LBC_g_regional_results_vol[,2], LBC_g_regional_results_sa[,2], LBC_g_regional_results_thk[,2])
colnames(allbetas) <- c("vol_STRADL", "sa_STRADL", "thk_STRADL", "vol_UKB", "sa_UKB", "thk_UKB", "vol_LBC", "sa_LBC", "thk_LBC")

# test for between-cohort regional-g profile consistency
allbetas <- as.data.frame(allbetas)
cor.test(allbetas$vol_STRADL, allbetas$vol_LBC)
cor.test(allbetas$vol_STRADL, allbetas$vol_UKB)
cor.test(allbetas$vol_LBC, allbetas$vol_UKB)
cor.test(allbetas$sa_STRADL, allbetas$sa_LBC)
cor.test(allbetas$sa_STRADL, allbetas$sa_UKB)
cor.test(allbetas$sa_LBC, allbetas$sa_UKB)
cor.test(allbetas$thk_STRADL, allbetas$thk_LBC)
cor.test(allbetas$thk_STRADL, allbetas$thk_UKB)
cor.test(allbetas$thk_LBC, allbetas$thk_UKB)


UKBmetainfo <- data.frame(Cohort = rep("UKB", 68), n = rep(37840, 68), meanage =  rep(64.91, 68))
UKB_g_regional_results_vol <- cbind(UKB_g_regional_results_vol, UKBmetainfo)
UKB_g_regional_results_sa <- cbind(UKB_g_regional_results_sa, UKBmetainfo)
UKB_g_regional_results_thk <- cbind(UKB_g_regional_results_thk, UKBmetainfo)

STRADLmetainfo <- data.frame(Cohort = rep("STRADL", 68), n = rep(1043, 68), meanage = rep(59.29))
STRADL_g_regional_results_vol <- cbind(STRADL_g_regional_results_vol, STRADLmetainfo)
STRADL_g_regional_results_sa <- cbind(STRADL_g_regional_results_sa, STRADLmetainfo)
STRADL_g_regional_results_thk <- cbind(STRADL_g_regional_results_thk, STRADLmetainfo)

LBCmetainfo <- data.frame(Cohort = rep("LBC", 68), n = rep(636, 68), meanage = rep(72.67))
LBC_g_regional_results_vol <- cbind(LBC_g_regional_results_vol, LBCmetainfo)
LBC_g_regional_results_sa <- cbind(LBC_g_regional_results_sa, LBCmetainfo)
LBC_g_regional_results_thk <- cbind(LBC_g_regional_results_thk, LBCmetainfo)
			       
vol_outcomes <- matrix(0, 68, 11)
volage_outcomes <- matrix(0, 68, 11)
for (i in 1:68) {
loopdata <- rbind(UKB_g_regional_results_vol[i,], STRADL_g_regional_results_vol[i,], LBC_g_regional_results_vol[i,])
res <- rma(beta, se^2, data = loopdata, ni = n)
resmod <- rma(beta, se^2, mods  =  meanage, data = loopdata, ni = n)
vol_outcomes[i,1] <- as.numeric(res[17]) # q statistic 
vol_outcomes[i,2] <- as.numeric(res[18]) # p value for q
vol_outcomes[i,3] <- as.numeric(res[13]) # I^2
vol_outcomes[i,3] <- as.numeric(res[9]) # tau^2
vol_outcomes[i,5] <- as.numeric(res[10]) # setau^2
vol_outcomes[i,6] <- as.numeric(res[2]) #beta
vol_outcomes[i,7] <- as.numeric(res[3]) # se
vol_outcomes[i,8] <- as.numeric(res[4]) # z
vol_outcomes[i,9] <- as.numeric(res[5]) # p value
vol_outcomes[i,10] <- loopdata$region[1]
vol_outcomes[i,11] <- as.numeric(res[14])
colnames(vol_outcomes) <- c("vol_q", "vol_p_q", "vol_I^2", "vol_tau^2", "vol_setau^2", "vol_beta", "vol_se", "vol_z", "vol_p", "Region", "vol_H^2")

volage_outcomes[i,1] <- as.numeric(resmod[17]) # q statistic 
volage_outcomes[i,2] <- as.numeric(resmod[18]) # p value for q
volage_outcomes[i,3] <- as.numeric(resmod[13]) # I^2
volage_outcomes[i,3] <- as.numeric(resmod[9]) # tau^2
volage_outcomes[i,5] <- as.numeric(resmod[10]) # setau^2
volage_outcomes[i,6] <- as.numeric(as.data.frame(resmod[2])[2,1]) #beta
volage_outcomes[i,7] <- as.numeric(as.data.frame(resmod[3])[2,1])  # se
volage_outcomes[i,8] <- as.numeric(as.data.frame(resmod[4])[2,1])  # z
volage_outcomes[i,9] <- as.numeric(as.data.frame(resmod[5])[2,1])  # p value
volage_outcomes[i,10] <- loopdata$region[1]
	
vol_outcomes[i,11] <- as.numeric(res[14])
volage_outcomes[i,11] <- as.numeric(resmod[14])
colnames(volage_outcomes) <- c("vol_q", "vol_p_q", "vol_I^2", "vol_tau^2", "vol_setau^2", "vol_beta", "vol_se", "vol_z", "vol_p", "Region", "vol_H^2")
}

vol_outcomes <- as.data.frame(vol_outcomes)
vol_outcomes[,1:9] <-  sapply(vol_outcomes[,1:9], as.numeric)
vol_outcomes$bhQ <- p.adjust(vol_outcomes$vol_p, method = "BH")
vol_outcomes <- as.data.frame(vol_outcomes)
volage_outcomes <- as.data.frame(volage_outcomes)
volage_outcomes[,1:9] <-  sapply(volage_outcomes[,1:9], as.numeric)
volage_outcomes$bhQ <- p.adjust(volage_outcomes$vol_p, method = "BH")
volage_outcomes <- as.data.frame(volage_outcomes)
vol_outcomes$logQ <- log(vol_outcomes$bhQ)

sa_outcomes <- matrix(0, 68, 11)
saage_outcomes <- matrix(0, 68, 11)
for (i in 1:68) {
loopdata <- rbind(UKB_g_regional_results_sa[i,], STRADL_g_regional_results_sa[i,], LBC_g_regional_results_sa[i,])
res <- rma(beta, sa_va, data = loopdata, ni = n)
resmod <- rma(beta, sa_va, mods  =  meanage, data = loopdata, ni = n)
sa_outcomes[i,1] <- as.numeric(res[17]) 
sa_outcomes[i,2] <- as.numeric(res[18])
sa_outcomes[i,3] <- as.numeric(res[13]) 
sa_outcomes[i,3] <- as.numeric(res[9]) 
sa_outcomes[i,5] <- as.numeric(res[10])
sa_outcomes[i,6] <- as.numeric(res[2]) 
sa_outcomes[i,7] <- as.numeric(res[3])
sa_outcomes[i,8] <- as.numeric(res[4]) 
sa_outcomes[i,9] <- as.numeric(res[5]) 
sa_outcomes[i,10] <- loopdata$region[1]
sa_outcomes[i,11] <- as.numeric(res[14])
colnames(sa_outcomes) <- c("sa_q", "sa_p_q", "sa_I^2", "sa_tau^2", "sa_setau^2", "sa_beta", "sa_se", "sa_z", "sa_p", "Region", "sa_H^2")

saage_outcomes[i,1] <- as.numeric(resmod[17]) 
saage_outcomes[i,2] <- as.numeric(resmod[18]) 
saage_outcomes[i,3] <- as.numeric(resmod[13]) 
saage_outcomes[i,3] <- as.numeric(resmod[9]) 
saage_outcomes[i,5] <- as.numeric(resmod[10])
saage_outcomes[i,6] <- as.numeric(as.data.frame(resmod[2])[2,1])
saage_outcomes[i,7] <- as.numeric(as.data.frame(resmod[3])[2,1])  
saage_outcomes[i,8] <- as.numeric(as.data.frame(resmod[4])[2,1])  
saage_outcomes[i,9] <- as.numeric(as.data.frame(resmod[5])[2,1]) 
saage_outcomes[i,10] <- loopdata$region[1]
	
sa_outcomes[i,11] <- as.numeric(res[14])
saage_outcomes[i,11] <- as.numeric(resmod[14])
colnames(saage_outcomes) <- c("sa_q", "sa_p_q", "sa_I^2", "sa_tau^2", "sa_setau^2", "sa_beta", "sa_se", "sa_z", "sa_p", "Region", "sa_H^2")
}

sa_outcomes <- as.data.frame(sa_outcomes)
sa_outcomes[,1:9] <-  sapply(sa_outcomes[,1:9], as.numeric)
sa_outcomes$bhQ <- p.adjust(sa_outcomes$sa_p, method = "BH")
sa_outcomes <- as.data.frame(sa_outcomes)
saage_outcomes <- as.data.frame(saage_outcomes)
saage_outcomes[,1:9] <-  sapply(saage_outcomes[,1:9], as.numeric)
saage_outcomes$bhQ <- p.adjust(saage_outcomes$sa_p, method = "BH")
saage_outcomes <- as.data.frame(saage_outcomes)
sa_outcomes$logQ <- log(sa_outcomes$bhQ)

thk_outcomes <- matrix(0, 68, 11)
thkage_outcomes <- matrix(0, 68, 11)
for (i in 1:68) {
loopdata <- rbind(UKB_g_regional_results_thk[i,], STRADL_g_regional_results_thk[i,], LBC_g_regional_results_thk[i,])
res <- rma(beta, thk_va, data = loopdata, ni = n)
resmod <- rma(beta, thk_va, mods  =  meanage, data = loopdata, ni = n)
thk_outcomes[i,1] <- as.numeric(res[17]) 
thk_outcomes[i,2] <- as.numeric(res[18]) 
thk_outcomes[i,3] <- as.numeric(res[13]) 
thk_outcomes[i,3] <- as.numeric(res[9]) 
thk_outcomes[i,5] <- as.numeric(res[10]) 
thk_outcomes[i,6] <- as.numeric(res[2]) 
thk_outcomes[i,7] <- as.numeric(res[3]) 
thk_outcomes[i,8] <- as.numeric(res[4]) 
thk_outcomes[i,9] <- as.numeric(res[5]) 
thk_outcomes[i,10] <- loopdata$region[1]
thk_outcomes[i,11] <- as.numeric(res[14])
colnames(thk_outcomes) <- c("thk_q", "thk_p_q", "thk_I^2", "thk_tau^2", "thk_setau^2", "thk_beta", "thk_se", "thk_z", "thk_p", "Region", "thk_H^2")

thkage_outcomes[i,1] <- as.numeric(resmod[17])
thkage_outcomes[i,2] <- as.numeric(resmod[18]) 
thkage_outcomes[i,3] <- as.numeric(resmod[13])
thkage_outcomes[i,3] <- as.numeric(resmod[9]) 
thkage_outcomes[i,5] <- as.numeric(resmod[10])
thkage_outcomes[i,6] <- as.numeric(as.data.frame(resmod[2])[2,1]) 
thkage_outcomes[i,7] <- as.numeric(as.data.frame(resmod[3])[2,1]) 
thkage_outcomes[i,8] <- as.numeric(as.data.frame(resmod[4])[2,1]) 
thkage_outcomes[i,9] <- as.numeric(as.data.frame(resmod[5])[2,1]) 
thkage_outcomes[i,10] <- loopdata$region[1]
	
thk_outcomes[i,11] <- as.numeric(res[14])
thkage_outcomes[i,11] <- as.numeric(resmod[14])
colnames(thkage_outcomes) <- c("thk_q", "thk_p_q", "thk_I^2", "thk_tau^2", "thk_setau^2", "thk_beta", "thk_se", "thk_z", "thk_p", "Region", "thk_H^2")
}

thk_outcomes <- as.data.frame(thk_outcomes)
thk_outcomes[,1:9] <-  sapply(thk_outcomes[,1:9], as.numeric)
thk_outcomes$bhQ <- p.adjust(thk_outcomes$thk_p, method = "BH")
thk_outcomes <- as.data.frame(thk_outcomes)
thkage_outcomes <- as.data.frame(thkage_outcomes)
thkage_outcomes[,1:9] <-  sapply(thkage_outcomes[,1:9], as.numeric)
thkage_outcomes$bhQ <- p.adjust(thkage_outcomes$thk_p, method = "BH")
thkage_outcomes <- as.data.frame(thkage_outcomes)
thk_outcomes$logQ <- log(thk_outcomes$bhQ)

write.table(vol_outcomes, "voloutcomes.csv", sep = ",")
write.table(sa_outcomes, "saoutcomes.csv", sep = ",")
write.table(thk_outcomes, "thkoutcomes.csv", sep = ",")

write.table(volage_outcomes, "volage_outcomes.csv", sep = ",")
write.table(saage_outcomes, "saage_outcomes.csv", sep = ",")
write.table(thkage_outcomes, "thkage_outcomes.csv", sep = ",")

meta_out <- cbind(vol_outcomes$Region, vol_outcomes[,c(6,7,9)], sa_outcomes[,c(6,7,9)], thk_outcomes[,c(6,7,9)])
colnames(meta_out)[1] <- "region"
write.table(meta_out, "data_metanalysis_output_estimates.csv", sep = ",", row.names = F)


