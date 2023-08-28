# Script written by Joanna Moodie, May-June 2023. The Desikan-Killiany regional-g associations are found for each cohort: UKB, STRADL and LBC1936, and then meta-analysed. 
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

outliers <- function(x, SD) {
a <- which(x> ((SD * sd(x, na.rm=T))+ (mean(x, na.rm=T))))
b <- which(x< ((-SD * sd(x, na.rm=T))+ (mean(x, na.rm=T))))
c <- c(a,b)
result <- c
return(result)
}

######## UKB


UKB <- read.csv("") # This file contains demographic information, the Desikan-Killiany 68 region parcellation values (columns) for participants who took part in MRI imaging in instance 2 (~40,000) and the cognitive test results from the 11 tests we included to estimate a latent g factor 
# prepare cognitive test data
UKB$cog_trailB_log[which(UKB$cog_trailB_log == 0)] <- NA
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
			       
# the desikan-killiany regions have slightly different names between freesurfer versions. Therefore, while this is a long #write-out, it makes sure that the variable names will match up between cohorts (which is helpful later for entering the results into the meta-analysis etc.). 
theUKBdata <- data.frame(ID = UKB$ID,
 ageMRI = UKB$ageyears_MRI,
 dontuse = UKB$ageyears_MRI,
 lag = rep(0,
 dim(UKB)[1]),
 sex = UKB$sex,
 atrophy = UKB$ICV,
 cohort = rep("UKB",
 dim(UKB)[1]),
 freesurfer = rep("v6",
 dim(UKB)[1]),
 site = UKB$assCtr,
 X = UKB$headposX,
 Y = UKB$headposY,
 Z = UKB$headposZ/10,
 lh_vol_bankssts = UKB$lh_bankssts_volume/1000,
 lh_vol_caudalanteriorcingulate = UKB$lh_caudalanteriorcingulate_volume/1000,
 lh_vol_caudalmiddlefrontal = UKB$lh_caudalmiddlefrontal_volume/1000,
 lh_vol_cuneus = UKB$lh_cuneus_volume/1000,
 lh_vol_entorhinal = UKB$lh_entorhinal_volume/1000,
 lh_vol_fusiform = UKB$lh_fusiform_volume/1000,
 lh_vol_inferiorparietal = UKB$lh_inferiorparietal_volume/1000,
 lh_vol_inferiortemporal = UKB$lh_inferiortemporal_volume/1000,
 lh_vol_isthmus = UKB$lh_isthmuscingulate_volume/1000,
 lh_vol_lateraloccipital = UKB$lh_lateraloccipital_volume/1000,
 lh_vol_lateralorbitofrontal = UKB$lh_lateralorbitofrontal_volume/1000,
 lh_vol_lingual = UKB$lh_lingual_volume/1000,
 lh_vol_medialorbitofrontal = UKB$lh_medialorbitofrontal_volume/1000,
 lh_vol_middletemporal = UKB$lh_middletemporal_volume/1000,
 lh_vol_parahippocampal= UKB$lh_parahippocampal_volume/1000,
 lh_vol_paracentral= UKB$lh_paracentral_volume/1000, 
 lh_vol_parsopercularis = UKB$lh_parsopercularis_volume/1000,
 lh_vol_parsorbitalis = UKB$lh_parsorbitalis_volume/1000,
 lh_vol_parstriangularis = UKB$lh_parstriangularis_volume/1000,
 lh_vol_pericalcarine = UKB$lh_pericalcarine_volume/1000,
 lh_vol_postcentral = UKB$lh_postcentral_volume/1000,
 lh_vol_posteriorcingulate= UKB$lh_posteriorcingulate_volume/1000,
 lh_vol_precentral = UKB$lh_precentral_volume/1000,
 lh_vol_precuneus= UKB$lh_precuneus_volume/1000,
 lh_vol_rostralanteriorcingulate = UKB$lh_rostralanteriorcingulate_volume/1000,
 lh_vol_rostralmiddlefrontal = UKB$lh_rostralmiddlefrontal_volume/1000,
 lh_vol_superiorfrontal = UKB$lh_superiorfrontal_volume/1000,
 lh_vol_superiorparietal = UKB$lh_superiorparietal_volume/1000,
 lh_vol_superiortemporal = UKB$lh_superiortemporal_volume/1000,
 lh_vol_supramarginal = UKB$lh_supramarginal_volume/1000,
 lh_vol_frontalpole = UKB$lh_frontalpole_volume/1000,
 lh_vol_transversetemporal = UKB$lh_transversetemporal_volume/1000,
 lh_vol_insula = UKB$lh_insula_volume/1000,
 lh_sa_bankssts = UKB$lh_bankssts_area/100,
 lh_sa_caudalanteriorcingulate = UKB$lh_caudalanteriorcingulate_area/100,
 lh_sa_caudalmiddlefrontal = UKB$lh_caudalmiddlefrontal_area/100,
 lh_sa_cuneus = UKB$lh_cuneus_area/100,
 lh_sa_entorhinal = UKB$lh_entorhinal_area/100,
 lh_sa_fusiform= UKB$lh_fusiform_area/100,
 lh_sa_inferiorparietal = UKB$lh_inferiorparietal_area/100,
 lh_sa_inferiortemporal = UKB$lh_inferiortemporal_area/100,
 lh_sa_isthmus = UKB$lh_isthmuscingulate_area/100,
 lh_sa_lateraloccipital = UKB$lh_lateraloccipital_area/100,
 lh_sa_lateralorbitofrontal = UKB$lh_lateralorbitofrontal_area/100,
 lh_sa_lingual =UKB$lh_lingual_area/100,
 lh_sa_medialorbitofrontal = UKB$lh_medialorbitofrontal_area/100,
 lh_sa_middletemporal = UKB$lh_middletemporal_area/100,
 lh_sa_parahippocampal = UKB$lh_parahippocampal_area/100,
 lh_sa_paracentral = UKB$lh_paracentral_area/100,
 lh_sa_parsopercularis = UKB$lh_parsopercularis_area/100,
 lh_sa_parsorbitalis = UKB$lh_parsorbitalis_area/100,
 lh_sa_parstriangularis = UKB$lh_parstriangularis_area/100,
 lh_sa_pericalcarine = UKB$lh_pericalcarine_area/100,
 lh_sa_postcentral = UKB$lh_postcentral_area/100,
 lh_sa_posteriorcingulate= UKB$lh_posteriorcingulate_area/100,
 lh_sa_precentral = UKB$lh_precentral_area/100,
 lh_sa_precuneus = UKB$lh_precuneus_area/100,
 lh_sa_rostralanteriorcingulate = UKB$lh_rostralanteriorcingulate_area/100,
 lh_sa_rostralmiddlefrontal = UKB$lh_rostralmiddlefrontal_area/100,
 lh_sa_superiorfrontal = UKB$lh_superiorfrontal_area/100,
 lh_sa_superiorparietal = UKB$lh_superiorparietal_area/100,
 lh_sa_superiortemporal = UKB$lh_superiortemporal_area/100,
 lh_sa_supramarginal = UKB$lh_supramarginal_area/100,
 lh_sa_frontalpole = UKB$lh_frontalpole_area/100,
 lh_sa_transversetemporal = UKB$lh_transversetemporal_area/100,
 lh_sa_insula = UKB$lh_insula_area/100,
 lh_thk_bankssts = UKB$lh_bankssts_thickness*10,
 lh_thk_caudalanteriorcingulate = UKB$lh_caudalanteriorcingulate_thickness*10,
 lh_thk_caudalmiddlefrontal = UKB$lh_caudalmiddlefrontal_thickness*10,
 lh_thk_cuneus = UKB$lh_cuneus_thickness*10,
 lh_thk_entorhinal = UKB$lh_entorhinal_thickness*10,
 lh_thk_fusiform = UKB$lh_fusiform_thickness*10,
 lh_thk_inferiorparietal = UKB$lh_inferiorparietal_thickness*10,
 lh_thk_inferiortemporal = UKB$lh_inferiortemporal_thickness*10,
 lh_thk_isthmus = UKB$lh_isthmuscingulate_thickness*10,
 lh_thk_lateraloccipital = UKB$lh_lateraloccipital_thickness*10,
 lh_thk_lateralorbitofrontal = UKB$lh_lateralorbitofrontal_thickness*10,
 lh_thk_lingual = UKB$lh_lingual_thickness*10,
 lh_thk_medialorbitofrontal = UKB$lh_medialorbitofrontal_thickness*10,
 lh_thk_middletemporal = UKB$lh_middletemporal_thickness*10,
 lh_thk_parahippocampal = UKB$lh_parahippocampal_thickness*10,
 lh_thk_paracentral = UKB$lh_paracentral_thickness*10,
 lh_thk_parsopercularis = UKB$lh_parsopercularis_thickness*10,
 lh_thk_parsorbitalis = UKB$lh_parsorbitalis_thickness*10,
 lh_thk_parstriangularis = UKB$lh_parstriangularis_thickness*10,
 lh_thk_pericalcarine = UKB$lh_pericalcarine_thickness*10,
 lh_thk_postcentral = UKB$lh_postcentral_thickness*10,
 lh_thk_posteriorcingulate = UKB$lh_posteriorcingulate_thickness*10,
 lh_thk_precentral = UKB$lh_precentral_thickness*10,
 lh_thk_precuneus = UKB$lh_precuneus_thickness*10,
 lh_thk_rostralanteriorcingulate = UKB$lh_rostralanteriorcingulate_thickness*10,
 lh_thk_rostralmiddlefrontal = UKB$lh_rostralmiddlefrontal_thickness*10,
 lh_thk_superiorfrontal = UKB$lh_superiorfrontal_thickness*10,
 lh_thk_superiorparietal = UKB$lh_superiorparietal_thickness*10,
 lh_thk_superiortemporal = UKB$lh_superiortemporal_thickness*10,
 lh_thk_supramarginal = UKB$lh_supramarginal_thickness*10,
 lh_thk_frontalpole = UKB$lh_frontalpole_thickness*10,
 lh_thk_transversetemporal = UKB$lh_transversetemporal_thickness*10,
 lh_thk_insula = UKB$lh_insula_thickness*10,
 rh_vol_bankssts = UKB$rh_bankssts_volume/1000,
 rh_vol_caudalanteriorcingulate = UKB$rh_caudalanteriorcingulate_volume/1000,
 rh_vol_caudalmiddlefrontal = UKB$rh_caudalmiddlefrontal_volume/1000,
 rh_vol_cuneus = UKB$rh_cuneus_volume/1000,
 rh_vol_entorhinal = UKB$rh_entorhinal_volume/1000,
 rh_vol_fusiform = UKB$rh_fusiform_volume/1000,
 rh_vol_inferiorparietal = UKB$rh_inferiorparietal_volume/1000,
 rh_vol_inferiortemporal = UKB$rh_inferiortemporal_volume/1000,
 rh_vol_isthmus = UKB$rh_isthmuscingulate_volume/1000,
 rh_vol_lateraloccipital = UKB$rh_lateraloccipital_volume/1000,
 rh_vol_lateralorbitofrontal = UKB$rh_lateralorbitofrontal_volume/1000,
 rh_vol_lingual = UKB$rh_lingual_volume/1000,
 rh_vol_medialorbitofrontal = UKB$rh_medialorbitofrontal_volume/1000,
 rh_vol_middletemporal = UKB$rh_middletemporal_volume/1000,
 rh_vol_parahippocampal= UKB$rh_parahippocampal_volume/1000,
 rh_vol_paracentral= UKB$rh_paracentral_volume/1000,
 rh_vol_parsopercularis = UKB$rh_parsopercularis_volume/1000,
 rh_vol_parsorbitalis = UKB$rh_parsorbitalis_volume/1000,
 rh_vol_parstriangularis = UKB$rh_parstriangularis_volume/1000,
 rh_vol_pericalcarine = UKB$rh_pericalcarine_volume/1000,
 rh_vol_postcentral = UKB$rh_postcentral_volume/1000,
 rh_vol_posteriorcingulate= UKB$rh_posteriorcingulate_volume/1000,
 rh_vol_precentral = UKB$rh_precentral_volume/1000,
 rh_vol_precuneus= UKB$rh_precuneus_volume/1000,
 rh_vol_rostralanteriorcingulate = UKB$rh_rostralanteriorcingulate_volume/1000,
 rh_vol_rostralmiddlefrontal = UKB$rh_rostralmiddlefrontal_volume/1000,
 rh_vol_superiorfrontal = UKB$rh_superiorfrontal_volume/1000,
 rh_vol_superiorparietal = UKB$rh_superiorparietal_volume/1000,
 rh_vol_superiortemporal = UKB$rh_superiortemporal_volume/1000,
 rh_vol_supramarginal = UKB$rh_supramarginal_volume/1000,
 rh_vol_frontalpole = UKB$rh_frontalpole_volume/1000,
 rh_vol_transversetemporal = UKB$rh_transversetemporal_volume/1000,
 rh_vol_insula = UKB$rh_insula_volume/1000,
 rh_sa_bankssts = UKB$rh_bankssts_area/100,
 rh_sa_caudalanteriorcingulate = UKB$rh_caudalanteriorcingulate_area/100,
 rh_sa_caudalmiddlefrontal = UKB$rh_caudalmiddlefrontal_area/100,
 rh_sa_cuneus = UKB$rh_cuneus_area/100,
 rh_sa_entorhinal = UKB$rh_entorhinal_area/100, 
 rh_sa_fusiform= UKB$rh_fusiform_area/100,
 rh_sa_inferiorparietal = UKB$rh_inferiorparietal_area/100,
 rh_sa_inferiortemporal = UKB$rh_inferiortemporal_area/100,
 rh_sa_isthmus = UKB$rh_isthmuscingulate_area/100,
 rh_sa_lateraloccipital = UKB$rh_lateraloccipital_area/100,
 rh_sa_lateralorbitofrontal = UKB$rh_lateralorbitofrontal_area/100,
 rh_sa_lingual =UKB$rh_lingual_area/100,
 rh_sa_medialorbitofrontal = UKB$rh_medialorbitofrontal_area/100,
 rh_sa_middletemporal = UKB$rh_middletemporal_area/100,
 rh_sa_parahippocampal = UKB$rh_parahippocampal_area/100,
 rh_sa_paracentral = UKB$rh_paracentral_area/100,
 rh_sa_parsopercularis = UKB$rh_parsopercularis_area/100,
 rh_sa_parsorbitalis = UKB$rh_parsorbitalis_area/100,
 rh_sa_parstriangularis = UKB$rh_parstriangularis_area/100,
 rh_sa_pericalcarine = UKB$rh_pericalcarine_area/100,
 rh_sa_postcentral = UKB$rh_postcentral_area/100,
 rh_sa_posteriorcingulate= UKB$rh_posteriorcingulate_area/100,
 rh_sa_precentral = UKB$rh_precentral_area/100,
 rh_sa_precuneus = UKB$rh_precuneus_area/100,
 rh_sa_rostralanteriorcingulate = UKB$rh_rostralanteriorcingulate_area/100,
 rh_sa_rostralmiddlefrontal = UKB$rh_rostralmiddlefrontal_area/100,
 rh_sa_superiorfrontal = UKB$rh_superiorfrontal_area/100,
 rh_sa_superiorparietal = UKB$rh_superiorparietal_area/100,
 rh_sa_superiortemporal = UKB$rh_superiortemporal_area/100,
 rh_sa_supramarginal = UKB$rh_supramarginal_area/100,
 rh_sa_frontalpole = UKB$rh_frontalpole_area/100,
 rh_sa_transversetemporal = UKB$rh_transversetemporal_area/100,
 rh_sa_insula = UKB$rh_insula_area/100,
 rh_thk_bankssts = UKB$rh_bankssts_thickness*10,
 rh_thk_caudalanteriorcingulate = UKB$rh_caudalanteriorcingulate_thickness*10,
 rh_thk_caudalmiddlefrontal = UKB$rh_caudalmiddlefrontal_thickness*10,
 rh_thk_cuneus = UKB$rh_cuneus_thickness*10,
 rh_thk_entorhinal = UKB$rh_entorhinal_thickness*10,
 rh_thk_fusiform = UKB$rh_fusiform_thickness*10,
 rh_thk_inferiorparietal = UKB$rh_inferiorparietal_thickness*10,
 rh_thk_inferiortemporal = UKB$rh_inferiortemporal_thickness*10,
 rh_thk_isthmus = UKB$rh_isthmuscingulate_thickness*10,
 rh_thk_lateraloccipital = UKB$rh_lateraloccipital_thickness*10,
 rh_thk_lateralorbitofrontal = UKB$rh_lateralorbitofrontal_thickness*10,
 rh_thk_lingual = UKB$rh_lingual_thickness*10, 
 rh_thk_medialorbitofrontal = UKB$rh_medialorbitofrontal_thickness*10,
 rh_thk_middletemporal = UKB$rh_middletemporal_thickness*10,
 rh_thk_parahippocampal = UKB$rh_parahippocampal_thickness*10,
 rh_thk_paracentral = UKB$rh_paracentral_thickness*10,
 rh_thk_parsopercularis = UKB$rh_parsopercularis_thickness*10,
 rh_thk_parsorbitalis = UKB$rh_parsorbitalis_thickness*10,
 rh_thk_parstriangularis = UKB$rh_parstriangularis_thickness*10,
 rh_thk_pericalcarine = UKB$rh_pericalcarine_thickness*10,
 rh_thk_postcentral = UKB$rh_postcentral_thickness*10,
 rh_thk_posteriorcingulate = UKB$rh_posteriorcingulate_thickness*10,
 rh_thk_precentral = UKB$rh_precentral_thickness*10,
 rh_thk_precuneus = UKB$rh_precuneus_thickness*10,
 rh_thk_rostralanteriorcingulate = UKB$rh_rostralanteriorcingulate_thickness*10,
 rh_thk_rostralmiddlefrontal = UKB$rh_rostralmiddlefrontal_thickness*10,
 rh_thk_superiorfrontal = UKB$rh_superiorfrontal_thickness*10,
 rh_thk_superiorparietal = UKB$rh_superiorparietal_thickness*10,
 rh_thk_superiortemporal = UKB$rh_superiortemporal_thickness*10,
 rh_thk_supramarginal = UKB$rh_supramarginal_thickness*10,
 rh_thk_frontalpole = UKB$rh_frontalpole_thickness*10,
 rh_thk_transversetemporal = UKB$rh_transversetemporal_thickness*10,
 rh_thk_insula = UKB$rh_insula_thickness*10,
 lh_vol_temporalpole = UKB$lh_temporalpole_volume/1000,
 rh_vol_temporalpole = UKB$rh_temporalpole_volume/1000,
 lh_sa_temporalpole = UKB$lh_temporalpole_area/100,
 rh_sa_temporalpole = UKB$rh_temporalpole_area/100,
 lh_thk_temporalpole = UKB$lh_temporalpole_thickness/10,
 rh_thk_temporalpole = UKB$rh_temporalpole_thickness/10)

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
STRADL <- read.csv("") # this file contains demographic information, the Desikan-Killiany regional volume, surface area, and thickness, and cognitive test scores for each participant.  
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

# make sure the regional 
theSTRADLdata = data.frame(ID = STRADL$ID, ageMRI = STRADL$AgeFaceToFace, dontuse = STRADL$AgeFaceToFace, lag = rep(0, dim(STRADL)[1]), sex = STRADL$Sex, atrophy = STRADL$atrophy, cohort = rep("STRADL", dim(STRADL)[1]), freesurfer = rep("unknown", dim(STRADL)[1]), site = STRADL$StudySite, X = STRADL$headposX, Y = STRADL$headposY, Z = STRADL$headposZ/10, lh_vol_bankssts = STRADL$lh_bankssts_volume/1000,lh_vol_caudalanteriorcingulate = STRADL$lh_caudalanteriorcingulate_volume/1000, lh_vol_caudalmiddlefrontal = STRADL$lh_caudalmiddlefrontal_volume/1000, lh_vol_cuneus = STRADL$lh_cuneus_volume/1000, lh_vol_entorhinal = STRADL$lh_entorhinal_volume/1000, lh_vol_fusiform = STRADL$lh_fusiform_volume/1000, lh_vol_inferiorparietal = STRADL$lh_inferiorparietal_volume/1000, lh_vol_inferiortemporal = STRADL$lh_inferiortemporal_volume/1000, lh_vol_isthmus = STRADL$lh_isthmuscingulate_volume/1000, lh_vol_lateraloccipital = STRADL$lh_lateraloccipital_volume/1000, lh_vol_lateralorbitofrontal = STRADL$lh_lateralorbitofrontal_volume/1000, lh_vol_lingual = STRADL$lh_lingual_volume/1000, lh_vol_medialorbitofrontal = STRADL$lh_medialorbitofrontal_volume/1000, lh_vol_middletemporal = STRADL$lh_middletemporal_volume/1000, lh_vol_parahippocampal= STRADL$lh_parahippocampal_volume/1000, lh_vol_paracentral= STRADL$lh_paracentral_volume/1000, 
lh_vol_parsopercularis = STRADL$lh_parsopercularis_volume/1000, lh_vol_parsorbitalis = STRADL$lh_parsorbitalis_volume/1000, 
                      lh_vol_parstriangularis = STRADL$lh_parstriangularis_volume/1000, lh_vol_pericalcarine = STRADL$lh_pericalcarine_volume/1000, 
                      lh_vol_postcentral = STRADL$lh_postcentral_volume/1000, lh_vol_posteriorcingulate= STRADL$lh_posteriorcingulate_volume/1000, 
                      lh_vol_precentral = STRADL$lh_precentral_volume/1000, lh_vol_precuneus= STRADL$lh_precuneus_volume/1000, 
                      lh_vol_rostralanteriorcingulate = STRADL$lh_rostralanteriorcingulate_volume/1000, lh_vol_rostralmiddlefrontal = STRADL$lh_rostralmiddlefrontal_volume/1000,lh_vol_superiorfrontal = STRADL$lh_superiorfrontal_volume/1000, lh_vol_superiorparietal = STRADL$lh_superiorparietal_volume/1000, lh_vol_superiortemporal = STRADL$lh_superiortemporal_volume/1000, lh_vol_supramarginal = STRADL$lh_supramarginal_volume/1000, lh_vol_frontalpole = STRADL$lh_frontalpole_volume/1000, 
                      lh_vol_transversetemporal = STRADL$lh_transversetemporal_volume/1000, lh_vol_insula = STRADL$lh_insula_volume/1000,lh_sa_bankssts = STRADL$lh_bankssts_area/100, lh_sa_caudalanteriorcingulate = STRADL$lh_caudalanteriorcingulate_area/100, 
                      lh_sa_caudalmiddlefrontal = STRADL$lh_caudalmiddlefrontal_area/100, lh_sa_cuneus = STRADL$lh_cuneus_area/100, 
                      lh_sa_entorhinal = STRADL$lh_entorhinal_area/100, 
                      lh_sa_fusiform= STRADL$lh_fusiform_area/100, lh_sa_inferiorparietal = STRADL$lh_inferiorparietal_area/100, 
                      lh_sa_inferiortemporal = STRADL$lh_inferiortemporal_area/100,
                      lh_sa_isthmus = STRADL$lh_isthmuscingulate_area/100, lh_sa_lateraloccipital = STRADL$lh_lateraloccipital_area/100, 
                      lh_sa_lateralorbitofrontal = STRADL$lh_lateralorbitofrontal_area/100, lh_sa_lingual =STRADL$lh_lingual_area/100, 
                      lh_sa_medialorbitofrontal = STRADL$lh_medialorbitofrontal_area/100, lh_sa_middletemporal = STRADL$lh_middletemporal_area/100, 
                      lh_sa_parahippocampal = STRADL$lh_parahippocampal_area/100, lh_sa_paracentral = STRADL$lh_paracentral_area/100, 
                      lh_sa_parsopercularis = STRADL$lh_parsopercularis_area/100,lh_sa_parsorbitalis = STRADL$lh_parsorbitalis_area/100, 
                      lh_sa_parstriangularis = STRADL$lh_parstriangularis_area/100, lh_sa_pericalcarine = STRADL$lh_pericalcarine_area/100, 
                      lh_sa_postcentral = STRADL$lh_postcentral_area/100, lh_sa_posteriorcingulate= STRADL$lh_posteriorcingulate_area/100, 
                      lh_sa_precentral = STRADL$lh_precentral_area/100,lh_sa_precuneus = STRADL$lh_precuneus_area/100, 
                      lh_sa_rostralanteriorcingulate = STRADL$lh_rostralanteriorcingulate_area/100,
                      lh_sa_rostralmiddlefrontal = STRADL$lh_rostralmiddlefrontal_area/100, lh_sa_superiorfrontal = STRADL$lh_superiorfrontal_area/100, 
                      lh_sa_superiorparietal = STRADL$lh_superiorparietal_area/100, lh_sa_superiortemporal = STRADL$lh_superiortemporal_area/100, 
                      lh_sa_supramarginal = STRADL$lh_supramarginal_area/100, lh_sa_frontalpole = STRADL$lh_frontalpole_area/100, 
                      lh_sa_transversetemporal = STRADL$lh_transversetemporal_area/100, lh_sa_insula = STRADL$lh_insula_area/100, 
                      lh_thk_bankssts = STRADL$lh_bankssts_thickness*10,
                      lh_thk_caudalanteriorcingulate = STRADL$lh_caudalanteriorcingulate_thickness*10, lh_thk_caudalmiddlefrontal = STRADL$lh_caudalmiddlefrontal_thickness*10,
                      lh_thk_cuneus = STRADL$lh_cuneus_thickness*10, lh_thk_entorhinal = STRADL$lh_entorhinal_thickness*10, lh_thk_fusiform = STRADL$lh_fusiform_thickness*10, 
                      lh_thk_inferiorparietal = STRADL$lh_inferiorparietal_thickness*10, lh_thk_inferiortemporal = STRADL$lh_inferiortemporal_thickness*10, 
                      lh_thk_isthmus = STRADL$lh_isthmuscingulate_thickness*10, lh_thk_lateraloccipital = STRADL$lh_lateraloccipital_thickness*10, 
                      lh_thk_lateralorbitofrontal = STRADL$lh_lateralorbitofrontal_thickness*10,lh_thk_lingual = STRADL$lh_lingual_thickness*10, 
                      lh_thk_medialorbitofrontal = STRADL$lh_medialorbitofrontal_thickness*10, lh_thk_middletemporal = STRADL$lh_middletemporal_thickness*10,
                      lh_thk_parahippocampal = STRADL$lh_parahippocampal_thickness*10, lh_thk_paracentral = STRADL$lh_paracentral_thickness*10, 
                      lh_thk_parsopercularis = STRADL$lh_parsopercularis_thickness*10, lh_thk_parsorbitalis = STRADL$lh_parsorbitalis_thickness*10, 
                      lh_thk_parstriangularis = STRADL$lh_parstriangularis_thickness*10, lh_thk_pericalcarine = STRADL$lh_pericalcarine_thickness*10, 
                      lh_thk_postcentral = STRADL$lh_postcentral_thickness*10, lh_thk_posteriorcingulate = STRADL$lh_posteriorcingulate_thickness*10, 
                      lh_thk_precentral = STRADL$lh_precentral_thickness*10,lh_thk_precuneus = STRADL$lh_precuneus_thickness*10, 
                      lh_thk_rostralanteriorcingulate = STRADL$lh_rostralanteriorcingulate_thickness*10, lh_thk_rostralmiddlefrontal = STRADL$lh_rostralmiddlefrontal_thickness*10,
                      lh_thk_superiorfrontal = STRADL$lh_superiorfrontal_thickness*10, lh_thk_superiorparietal = STRADL$lh_superiorparietal_thickness*10, 
                      lh_thk_superiortemporal = STRADL$lh_superiortemporal_thickness*10, lh_thk_supramarginal = STRADL$lh_supramarginal_thickness*10, 
                      lh_thk_frontalpole = STRADL$lh_frontalpole_thickness*10,
                      lh_thk_transversetemporal = STRADL$lh_transversetemporal_thickness*10, lh_thk_insula = STRADL$lh_insula_thickness*10,rh_vol_bankssts = STRADL$rh_bankssts_volume/1000, rh_vol_caudalanteriorcingulate = STRADL$rh_caudalanteriorcingulate_volume/1000, rh_vol_caudalmiddlefrontal = STRADL$rh_caudalmiddlefrontal_volume/1000, 
                      rh_vol_cuneus = STRADL$rh_cuneus_volume/1000, rh_vol_entorhinal = STRADL$rh_entorhinal_volume/1000, rh_vol_fusiform = STRADL$rh_fusiform_volume/1000, rh_vol_inferiorparietal = STRADL$rh_inferiorparietal_volume/1000, rh_vol_inferiortemporal = STRADL$rh_inferiortemporal_volume/1000, rh_vol_isthmus = STRADL$rh_isthmuscingulate_volume/1000, rh_vol_lateraloccipital = STRADL$rh_lateraloccipital_volume/1000, rh_vol_lateralorbitofrontal = STRADL$rh_lateralorbitofrontal_volume/1000, rh_vol_lingual = STRADL$rh_lingual_volume/1000, rh_vol_medialorbitofrontal = STRADL$rh_medialorbitofrontal_volume/1000, rh_vol_middletemporal = STRADL$rh_middletemporal_volume/1000, rh_vol_parahippocampal= STRADL$rh_parahippocampal_volume/1000, rh_vol_paracentral= STRADL$rh_paracentral_volume/1000, 
                      rh_vol_parsopercularis = STRADL$rh_parsopercularis_volume/1000, rh_vol_parsorbitalis = STRADL$rh_parsorbitalis_volume/1000, 
                      rh_vol_parstriangularis = STRADL$rh_parstriangularis_volume/1000, rh_vol_pericalcarine = STRADL$rh_pericalcarine_volume/1000, 
                      rh_vol_postcentral = STRADL$rh_postcentral_volume/1000, rh_vol_posteriorcingulate= STRADL$rh_posteriorcingulate_volume/1000, 
                      rh_vol_precentral = STRADL$rh_precentral_volume/1000, rh_vol_precuneus= STRADL$rh_precuneus_volume/1000, 
                      rh_vol_rostralanteriorcingulate = STRADL$rh_rostralanteriorcingulate_volume/1000, rh_vol_rostralmiddlefrontal = STRADL$rh_rostralmiddlefrontal_volume/1000,rh_vol_superiorfrontal = STRADL$rh_superiorfrontal_volume/1000, rh_vol_superiorparietal = STRADL$rh_superiorparietal_volume/1000, rh_vol_superiortemporal = STRADL$rh_superiortemporal_volume/1000, rh_vol_supramarginal = STRADL$rh_supramarginal_volume/1000, rh_vol_frontalpole = STRADL$rh_frontalpole_volume/1000, 
                      rh_vol_transversetemporal = STRADL$rh_transversetemporal_volume/1000, rh_vol_insula = STRADL$rh_insula_volume/1000,rh_sa_bankssts = STRADL$rh_bankssts_area/100, rh_sa_caudalanteriorcingulate = STRADL$rh_caudalanteriorcingulate_area/100, 
                      rh_sa_caudalmiddlefrontal = STRADL$rh_caudalmiddlefrontal_area/100, rh_sa_cuneus = STRADL$rh_cuneus_area/100, 
                      rh_sa_entorhinal = STRADL$rh_entorhinal_area/100, 
                      rh_sa_fusiform= STRADL$rh_fusiform_area/100, rh_sa_inferiorparietal = STRADL$rh_inferiorparietal_area/100, 
                      rh_sa_inferiortemporal = STRADL$rh_inferiortemporal_area/100,
                      rh_sa_isthmus = STRADL$rh_isthmuscingulate_area/100, rh_sa_lateraloccipital = STRADL$rh_lateraloccipital_area/100, 
                      rh_sa_lateralorbitofrontal = STRADL$rh_lateralorbitofrontal_area/100, rh_sa_lingual =STRADL$rh_lingual_area/100, 
                      rh_sa_medialorbitofrontal = STRADL$rh_medialorbitofrontal_area/100, rh_sa_middletemporal = STRADL$rh_middletemporal_area/100, 
                      rh_sa_parahippocampal = STRADL$rh_parahippocampal_area/100, rh_sa_paracentral = STRADL$rh_paracentral_area/100, 
                      rh_sa_parsopercularis = STRADL$rh_parsopercularis_area/100,rh_sa_parsorbitalis = STRADL$rh_parsorbitalis_area/100, 
                      rh_sa_parstriangularis = STRADL$rh_parstriangularis_area/100, rh_sa_pericalcarine = STRADL$rh_pericalcarine_area/100, 
                      rh_sa_postcentral = STRADL$rh_postcentral_area/100, rh_sa_posteriorcingulate= STRADL$rh_posteriorcingulate_area/100, 
                      rh_sa_precentral = STRADL$rh_precentral_area/100,rh_sa_precuneus = STRADL$rh_precuneus_area/100, 
                      rh_sa_rostralanteriorcingulate = STRADL$rh_rostralanteriorcingulate_area/100,
                      rh_sa_rostralmiddlefrontal = STRADL$rh_rostralmiddlefrontal_area/100, rh_sa_superiorfrontal = STRADL$rh_superiorfrontal_area/100, 
                      rh_sa_superiorparietal = STRADL$rh_superiorparietal_area/100, rh_sa_superiortemporal = STRADL$rh_superiortemporal_area/100, 
                      rh_sa_supramarginal = STRADL$rh_supramarginal_area/100, rh_sa_frontalpole = STRADL$rh_frontalpole_area/100, 
                      rh_sa_transversetemporal = STRADL$rh_transversetemporal_area/100, rh_sa_insula = STRADL$rh_insula_area/100, 
                      rh_thk_bankssts = STRADL$rh_bankssts_thickness*10,    rh_thk_caudalanteriorcingulate = STRADL$rh_caudalanteriorcingulate_thickness*10, rh_thk_caudalmiddlefrontal = STRADL$rh_caudalmiddlefrontal_thickness*10,
                      rh_thk_cuneus = STRADL$rh_cuneus_thickness*10, rh_thk_entorhinal = STRADL$rh_entorhinal_thickness*10, rh_thk_fusiform = STRADL$rh_fusiform_thickness*10, 
                      rh_thk_inferiorparietal = STRADL$rh_inferiorparietal_thickness*10, rh_thk_inferiortemporal = STRADL$rh_inferiortemporal_thickness*10, 
                      rh_thk_isthmus = STRADL$rh_isthmuscingulate_thickness*10, rh_thk_lateraloccipital = STRADL$rh_lateraloccipital_thickness*10, 
                      rh_thk_lateralorbitofrontal = STRADL$rh_lateralorbitofrontal_thickness*10,rh_thk_lingual = STRADL$rh_lingual_thickness*10, 
                      rh_thk_medialorbitofrontal = STRADL$rh_medialorbitofrontal_thickness*10, rh_thk_middletemporal = STRADL$rh_middletemporal_thickness*10,
                      rh_thk_parahippocampal = STRADL$rh_parahippocampal_thickness*10, rh_thk_paracentral = STRADL$rh_paracentral_thickness*10, 
                      rh_thk_parsopercularis = STRADL$rh_parsopercularis_thickness*10, rh_thk_parsorbitalis = STRADL$rh_parsorbitalis_thickness*10, 
                      rh_thk_parstriangularis = STRADL$rh_parstriangularis_thickness*10, rh_thk_pericalcarine = STRADL$rh_pericalcarine_thickness*10, 
                      rh_thk_postcentral = STRADL$rh_postcentral_thickness*10, rh_thk_posteriorcingulate = STRADL$rh_posteriorcingulate_thickness*10, 
                      rh_thk_precentral = STRADL$rh_precentral_thickness*10,rh_thk_precuneus = STRADL$rh_precuneus_thickness*10, 
                      rh_thk_rostralanteriorcingulate = STRADL$rh_rostralanteriorcingulate_thickness*10, rh_thk_rostralmiddlefrontal = STRADL$rh_rostralmiddlefrontal_thickness*10,
                      rh_thk_superiorfrontal = STRADL$rh_superiorfrontal_thickness*10, rh_thk_superiorparietal = STRADL$rh_superiorparietal_thickness*10, 
                      rh_thk_superiortemporal = STRADL$rh_superiortemporal_thickness*10, rh_thk_supramarginal = STRADL$rh_supramarginal_thickness*10, 
                      rh_thk_frontalpole = STRADL$rh_frontalpole_thickness*10,
                      rh_thk_transversetemporal = STRADL$rh_transversetemporal_thickness*10, rh_thk_insula = STRADL$rh_insula_thickness*10, 
	              lh_vol_temporalpole = STRADL$lh_temporalpole_volume/1000, rh_vol_temporalpole = STRADL$rh_temporalpole_volume/1000,
		      lh_sa_temporalpole = STRADL$lh_temporalpole_area/100, rh_sa_temporalpole = STRADL$rh_temporalpole_area/100, lh_thk_temporalpole = STRADL$lh_temporalpole_thickness/10,
		      rh_thk_temporalpole = STRADL$rh_temporalpole_thickness/10)

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
write.table(STRADL_g_regional_results, '/home/jmoodie/Documents/Expression/Scripts_toshare/STRADL_g_regional_results.csv', row.names = F)
STRADL_g_regional_results_vol <- STRADL_g_regional_results[which(substr(STRADL_g_regional_results$region, 4, 6) == "vol"),]
STRADL_g_regional_results_sa <- STRADL_g_regional_results[which(substr(STRADL_g_regional_results$region, 4, 5) == "sa"),]
STRADL_g_regional_results_thk <- STRADL_g_regional_results[which(substr(STRADL_g_regional_results$region, 4, 6) == "thk"),]

allbetas <- cbind(STRADL_g_regional_results_vol[,2], STRADL_g_regional_results_sa[,2], STRADL_g_regional_results_thk[,2], UKB_g_regional_results_vol[,2], UKB_g_regional_results_sa[,2], UKB_g_regional_results_thk[,2])
colnames(allbetas) <- c("vol_STRADL", "sa_STRADL", "thk_STRADL", "vol_UKB", "sa_UKB", "thk_UKB")
cor(allbetas)

######## LBC



fsl <- read.csv('/home/jmoodie/Documents/fsl_XYZ_LBC1936.csv')

LBC <- read.csv("/home/jmoodie/Documents/Expression/LBC/data_LBC_expression.csv")

x <- read.csv("/home/jmoodie/Documents/Expression/All_together_now/xyz/data_LBC_xyz.csv", sep = ",", header = T)
LBC <- cbind(LBC, x)

theLBCCogdata <- LBC[,1:25]
corrplot::corrplot(cor(theLBCCogdata[,8:25], use = "complete.obs"))

theLBCCogdata$AgeDays <- theLBCCogdata$AgeDays
model2 <- 'g=~matrix_reasoning + block_design + spatial_span_total + NART + WTAR + verbal_fluency+verbal_paired_associates + logical_memory + digit_span_backward + symbol_search + digit_symbol + inspection_time + choice_reaction_time_reflected 

#matrix_reasoning ~ AgeDays + sex
#block_design~ AgeDays + sex
#spatial_span_total~ AgeDays + sex
#NART~ AgeDays + sex
#WTAR ~ AgeDays + sex
#verbal_fluency~ AgeDays + sex
#verbal_paired_associates~ AgeDays + sex
#logical_memory~ AgeDays + sex
#digit_span_backward~ AgeDays + sex
#symbol_search~ AgeDays + sex
#digit_symbol~ AgeDays + sex
#inspection_time~ AgeDays + sex
#choice_reaction_time_reflected~ AgeDays + sex

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

ggplot(data = GpredictLBC, aes(x = g)) + geom_density(color = "gold", fill = "gold", alpha = 0.5) + xlab("Scaled cognitive ability score") + ylab("Density\n") +theme_cowplot() + theme( text = element_text(size=29), axis.text=element_text(size=29), legend.title=element_blank()) + theme(legend.position = "none") +scale_y_continuous(expand = c(0, 0)) +scale_x_continuous(breaks = c(-4,-2, 0,2,4), limits = c(-4,4))


theLBCdata = data.frame(ID = LBC$ID, ageMRI = LBC$AgeMRI/365.24, dontuse = LBC$AgeDays/365.24, lag = (LBC$AgeMRI-LBC$AgeDays), sex = LBC$sex, atrophy = LBC$atrophy, cohort = rep("LBC1936", dim(LBC)[1]), freesurfer = rep("unknown", dim(LBC)[1]), site = rep("LBCsite", dim(LBC)[1]), X = LBC$headposX, Y = LBC$headposY, Z = LBC$headposZ/10, lh_vol_bankssts = LBC$lh_vol_bankssts/1000, lh_vol_caudalanteriorcingulate = LBC$lh_vol_caudalanteriorcingulate/1000, lh_vol_caudalmiddlefrontal = LBC$lh_vol_caudalmiddlefrontal/1000, 
 lh_vol_cuneus = LBC$lh_vol_cuneus/1000, lh_vol_entorhinal = LBC$lh_vol_entorhinal/1000, lh_vol_fusiform = LBC$lh_vol_fusiform/1000, lh_vol_inferiorparietal = LBC$lh_vol_inferiorparietal/1000, lh_vol_inferiortemporal = LBC$lh_vol_inferiortemporal/1000, lh_vol_isthmus = LBC$lh_vol_isthmus/1000, lh_vol_lateraloccipital = LBC$lh_vol_lateraloccipital/1000, lh_vol_lateralorbitofrontal = LBC$lh_vol_lateralorbitofrontal/1000, lh_vol_lingual = LBC$lh_vol_lingual/1000, lh_vol_medialorbitofrontal = LBC$lh_vol_medialorbitofrontal/1000, lh_vol_middletemporal = LBC$lh_vol_middletemporal/1000, lh_vol_parahippocampal= LBC$lh_vol_parahippocampal/1000, lh_vol_paracentral= LBC$lh_vol_paracentral/1000, 
 lh_vol_parsopercularis = LBC$lh_vol_parsopercularis/1000, lh_vol_parsorbitalis = LBC$lh_vol_parsorbitalis/1000, 
 lh_vol_parstriangularis = LBC$lh_vol_parstriangularis/1000, lh_vol_pericalcarine = LBC$lh_vol_pericalcarine/1000, 
 lh_vol_postcentral = LBC$lh_vol_postcentral/1000, lh_vol_posteriorcingulate= LBC$lh_vol_posteriorcingulate/1000, 
 lh_vol_precentral = LBC$lh_vol_precentral/1000, lh_vol_precuneus= LBC$lh_vol_precuneus/1000, 
 lh_vol_rostralanteriorcingulate = LBC$lh_vol_rostralanteriorcingulate/1000, lh_vol_rostralmiddlefrontal = LBC$lh_vol_rostralmiddlefrontal/1000,lh_vol_superiorfrontal = LBC$lh_vol_superiorfrontal/1000, lh_vol_superiorparietal = LBC$lh_vol_superiorparietal/1000, lh_vol_superiortemporal = LBC$lh_vol_superiortemporal/1000, lh_vol_supramarginal = LBC$lh_vol_supramarginal/1000, lh_vol_frontalpole = LBC$lh_vol_frontalpole/1000, 
 lh_vol_transversetemporal = LBC$lh_vol_transversetemporal/1000, lh_vol_insula = LBC$lh_vol_insula/1000,lh_sa_bankssts = LBC$lh_sa_bankssts/100, lh_sa_caudalanteriorcingulate = LBC$lh_sa_caudalanteriorcingulate/100, 
lh_sa_caudalmiddlefrontal = LBC$lh_sa_caudalmiddlefrontal/100, lh_sa_cuneus = LBC$lh_sa_cuneus/100, 
lh_sa_entorhinal = LBC$lh_sa_entorhinal/100, 
lh_sa_fusiform= LBC$lh_sa_fusiform/100, lh_sa_inferiorparietal = LBC$lh_sa_inferiorparietal/100, 
lh_sa_inferiortemporal = LBC$lh_sa_inferiortemporal/100,
lh_sa_isthmus = LBC$lh_sa_isthmus/100, lh_sa_lateraloccipital = LBC$lh_sa_lateraloccipital/100, 
lh_sa_lateralorbitofrontal = LBC$lh_sa_lateralorbitofrontal/100, lh_sa_lingual =LBC$lh_sa_lingual/100, 
lh_sa_medialorbitofrontal = LBC$lh_sa_medialorbitofrontal/100, lh_sa_middletemporal = LBC$lh_sa_middletemporal/100, 
lh_sa_parahippocampal = LBC$lh_sa_parahippocampal/100, lh_sa_paracentral = LBC$lh_sa_paracentral/100, 
lh_sa_parsopercularis = LBC$lh_sa_parsopercularis/100,lh_sa_parsorbitalis = LBC$lh_sa_parsorbitalis/100, 
lh_sa_parstriangularis = LBC$lh_sa_parstriangularis/100, lh_sa_pericalcarine = LBC$lh_sa_pericalcarine/100, 
lh_sa_postcentral = LBC$lh_sa_postcentral/100, lh_sa_posteriorcingulate= LBC$lh_sa_posteriorcingulate/100, 
lh_sa_precentral = LBC$lh_sa_precentral/100,lh_sa_precuneus = LBC$lh_sa_precuneus/100, 
lh_sa_rostralanteriorcingulate = LBC$lh_sa_rostralanteriorcingulate/100,
lh_sa_rostralmiddlefrontal = LBC$lh_sa_rostralmiddlefrontal/100, lh_sa_superiorfrontal = LBC$lh_sa_superiorfrontal/100, 
lh_sa_superiorparietal = LBC$lh_sa_superiorparietal/100, lh_sa_superiortemporal = LBC$lh_sa_superiortemporal/100, 
lh_sa_supramarginal = LBC$lh_sa_supramarginal/100, lh_sa_frontalpole = LBC$lh_sa_frontalpole/100, 
lh_sa_transversetemporal = LBC$lh_sa_transversetemporal/100, lh_sa_insula = LBC$lh_sa_insula/100, 
 lh_thk_bankssts = LBC$lh_thk_bankssts*10,
 lh_thk_caudalanteriorcingulate = LBC$lh_thk_caudalanteriorcingulate*10, lh_thk_caudalmiddlefrontal = LBC$lh_thk_caudalmiddlefrontal*10,
 lh_thk_cuneus = LBC$lh_thk_cuneus*10, lh_thk_entorhinal = LBC$lh_thk_entorhinal*10, lh_thk_fusiform = LBC$lh_thk_fusiform*10, 
 lh_thk_inferiorparietal = LBC$lh_thk_inferiorparietal*10, lh_thk_inferiortemporal = LBC$lh_thk_inferiortemporal*10, 
 lh_thk_isthmus = LBC$lh_thk_isthmus*10, lh_thk_lateraloccipital = LBC$lh_thk_lateraloccipital*10, 
 lh_thk_lateralorbitofrontal = LBC$lh_thk_lateralorbitofrontal*10,lh_thk_lingual = LBC$lh_thk_lingual*10, 
 lh_thk_medialorbitofrontal = LBC$lh_thk_medialorbitofrontal*10, lh_thk_middletemporal = LBC$lh_thk_middletemporal*10,
 lh_thk_parahippocampal = LBC$lh_thk_parahippocampal*10, lh_thk_paracentral = LBC$lh_thk_paracentral*10, 
 lh_thk_parsopercularis = LBC$lh_thk_parsopercularis*10, lh_thk_parsorbitalis = LBC$lh_thk_parsorbitalis*10, 
 lh_thk_parstriangularis = LBC$lh_thk_parstriangularis*10, lh_thk_pericalcarine = LBC$lh_thk_pericalcarine*10, 
 lh_thk_postcentral = LBC$lh_thk_postcentral*10, lh_thk_posteriorcingulate = LBC$lh_thk_posteriorcingulate*10, 
 lh_thk_precentral = LBC$lh_thk_precentral*10,lh_thk_precuneus = LBC$lh_thk_precuneus*10, 
 lh_thk_rostralanteriorcingulate = LBC$lh_thk_rostralanteriorcingulate*10, lh_thk_rostralmiddlefrontal = LBC$lh_thk_rostralmiddlefrontal*10,
 lh_thk_superiorfrontal = LBC$lh_thk_superiorfrontal*10, lh_thk_superiorparietal = LBC$lh_thk_superiorparietal*10, 
 lh_thk_superiortemporal = LBC$lh_thk_superiortemporal*10, lh_thk_supramarginal = LBC$lh_thk_supramarginal*10, 
 lh_thk_frontalpole = LBC$lh_thk_frontalpole*10, 
 lh_thk_transversetemporal = LBC$lh_thk_transversetemporal*10, lh_thk_insula = LBC$lh_thk_insula*10, rh_vol_bankssts = LBC$rh_vol_bankssts/1000, rh_vol_caudalanteriorcingulate = LBC$rh_vol_caudalanteriorcingulate/1000, rh_vol_caudalmiddlefrontal = LBC$rh_vol_caudalmiddlefrontal/1000, 
 rh_vol_cuneus = LBC$rh_vol_cuneus/1000, rh_vol_entorhinal = LBC$rh_vol_entorhinal/1000, rh_vol_fusiform = LBC$rh_vol_fusiform/1000, rh_vol_inferiorparietal = LBC$rh_vol_inferiorparietal/1000, rh_vol_inferiortemporal = LBC$rh_vol_inferiortemporal/1000, rh_vol_isthmus = LBC$rh_vol_isthmus/1000, rh_vol_lateraloccipital = LBC$rh_vol_lateraloccipital/1000, rh_vol_lateralorbitofrontal = LBC$rh_vol_lateralorbitofrontal/1000, rh_vol_lingual = LBC$rh_vol_lingual/1000, rh_vol_medialorbitofrontal = LBC$rh_vol_medialorbitofrontal/1000, rh_vol_middletemporal = LBC$rh_vol_middletemporal/1000, rh_vol_parahippocampal= LBC$rh_vol_parahippocampal/1000, rh_vol_paracentral= LBC$rh_vol_paracentral/1000, 
 rh_vol_parsopercularis = LBC$rh_vol_parsopercularis/1000, rh_vol_parsorbitalis = LBC$rh_vol_parsorbitalis/1000, 
 rh_vol_parstriangularis = LBC$rh_vol_parstriangularis/1000, rh_vol_pericalcarine = LBC$rh_vol_pericalcarine/1000, 
 rh_vol_postcentral = LBC$rh_vol_postcentral/1000, rh_vol_posteriorcingulate= LBC$rh_vol_posteriorcingulate/1000, 
 rh_vol_precentral = LBC$rh_vol_precentral/1000, rh_vol_precuneus= LBC$rh_vol_precuneus/1000, 
 rh_vol_rostralanteriorcingulate = LBC$rh_vol_rostralanteriorcingulate/1000, rh_vol_rostralmiddlefrontal = LBC$rh_vol_rostralmiddlefrontal/1000,rh_vol_superiorfrontal = LBC$rh_vol_superiorfrontal/1000, rh_vol_superiorparietal = LBC$rh_vol_superiorparietal/1000, rh_vol_superiortemporal = LBC$rh_vol_superiortemporal/1000, rh_vol_supramarginal = LBC$rh_vol_supramarginal/1000, rh_vol_frontalpole = LBC$rh_vol_frontalpole/1000, 
 rh_vol_transversetemporal = LBC$rh_vol_transversetemporal/1000, rh_vol_insula = LBC$rh_vol_insula/1000,rh_sa_bankssts = LBC$rh_sa_bankssts/100, rh_sa_caudalanteriorcingulate = LBC$rh_sa_caudalanteriorcingulate/100, 
rh_sa_caudalmiddlefrontal = LBC$rh_sa_caudalmiddlefrontal/100, rh_sa_cuneus = LBC$rh_sa_cuneus/100, 
rh_sa_entorhinal = LBC$rh_sa_entorhinal/100,
rh_sa_fusiform= LBC$rh_sa_fusiform/100, rh_sa_inferiorparietal = LBC$rh_sa_inferiorparietal/100, 
rh_sa_inferiortemporal = LBC$rh_sa_inferiortemporal/100,
rh_sa_isthmus = LBC$rh_sa_isthmus/100, rh_sa_lateraloccipital = LBC$rh_sa_lateraloccipital/100, 
rh_sa_lateralorbitofrontal = LBC$rh_sa_lateralorbitofrontal/100, rh_sa_lingual =LBC$rh_sa_lingual/100, 
rh_sa_medialorbitofrontal = LBC$rh_sa_medialorbitofrontal/100, rh_sa_middletemporal = LBC$rh_sa_middletemporal/100, 
rh_sa_parahippocampal = LBC$rh_sa_parahippocampal/100, rh_sa_paracentral = LBC$rh_sa_paracentral/100, 
rh_sa_parsopercularis = LBC$rh_sa_parsopercularis/100,rh_sa_parsorbitalis = LBC$rh_sa_parsorbitalis/100, 
rh_sa_parstriangularis = LBC$rh_sa_parstriangularis/100, rh_sa_pericalcarine = LBC$rh_sa_pericalcarine/100, 
rh_sa_postcentral = LBC$rh_sa_postcentral/100, rh_sa_posteriorcingulate= LBC$rh_sa_posteriorcingulate/100, 
rh_sa_precentral = LBC$rh_sa_precentral/100,rh_sa_precuneus = LBC$rh_sa_precuneus/100, 
rh_sa_rostralanteriorcingulate = LBC$rh_sa_rostralanteriorcingulate/100,
rh_sa_rostralmiddlefrontal = LBC$rh_sa_rostralmiddlefrontal/100, rh_sa_superiorfrontal = LBC$rh_sa_superiorfrontal/100, 
rh_sa_superiorparietal = LBC$rh_sa_superiorparietal/100, rh_sa_superiortemporal = LBC$rh_sa_superiortemporal/100, 
rh_sa_supramarginal = LBC$rh_sa_supramarginal/100, rh_sa_frontalpole = LBC$rh_sa_frontalpole/100, 
rh_sa_transversetemporal = LBC$rh_sa_transversetemporal/100, rh_sa_insula = LBC$rh_sa_insula/100, 
 rh_thk_bankssts = LBC$rh_thk_bankssts*10,
 rh_thk_caudalanteriorcingulate = LBC$rh_thk_caudalanteriorcingulate*10, rh_thk_caudalmiddlefrontal = LBC$rh_thk_caudalmiddlefrontal*10,
 rh_thk_cuneus = LBC$rh_thk_cuneus*10, rh_thk_entorhinal = LBC$rh_thk_entorhinal*10, rh_thk_fusiform = LBC$rh_thk_fusiform*10, 
 rh_thk_inferiorparietal = LBC$rh_thk_inferiorparietal*10, rh_thk_inferiortemporal = LBC$rh_thk_inferiortemporal*10, 
 rh_thk_isthmus = LBC$rh_thk_isthmus*10, rh_thk_lateraloccipital = LBC$rh_thk_lateraloccipital*10, 
 rh_thk_lateralorbitofrontal = LBC$rh_thk_lateralorbitofrontal*10,rh_thk_lingual = LBC$rh_thk_lingual*10, 
 rh_thk_medialorbitofrontal = LBC$rh_thk_medialorbitofrontal*10, rh_thk_middletemporal = LBC$rh_thk_middletemporal*10,
 rh_thk_parahippocampal = LBC$rh_thk_parahippocampal*10, rh_thk_paracentral = LBC$rh_thk_paracentral*10, 
 rh_thk_parsopercularis = LBC$rh_thk_parsopercularis*10, rh_thk_parsorbitalis = LBC$rh_thk_parsorbitalis*10, 
 rh_thk_parstriangularis = LBC$rh_thk_parstriangularis*10, rh_thk_pericalcarine = LBC$rh_thk_pericalcarine*10, 
 rh_thk_postcentral = LBC$rh_thk_postcentral*10, rh_thk_posteriorcingulate = LBC$rh_thk_posteriorcingulate*10, 
 rh_thk_precentral = LBC$rh_thk_precentral*10,rh_thk_precuneus = LBC$rh_thk_precuneus*10, 
 rh_thk_rostralanteriorcingulate = LBC$rh_thk_rostralanteriorcingulate*10, rh_thk_rostralmiddlefrontal = LBC$rh_thk_rostralmiddlefrontal*10,
 rh_thk_superiorfrontal = LBC$rh_thk_superiorfrontal*10, rh_thk_superiorparietal = LBC$rh_thk_superiorparietal*10, 
 rh_thk_superiortemporal = LBC$rh_thk_superiortemporal*10, rh_thk_supramarginal = LBC$rh_thk_supramarginal*10, 
 rh_thk_frontalpole = LBC$rh_thk_frontalpole*10, 
 rh_thk_transversetemporal = LBC$rh_thk_transversetemporal*10, rh_thk_insula = LBC$rh_thk_insula*10, lh_vol_temporalpole = LBC$lh_vol_temporalpole/1000, rh_vol_temporalpole = LBC$rh_vol_temporalpole/1000,lh_sa_temporalpole = LBC$lh_sa_temporalpole/100, rh_sa_temporalpole = LBC$rh_sa_temporalpole/100,lh_thk_temporalpole = LBC$lh_thk_temporalpole/10,rh_thk_temporalpole = LBC$rh_thk_temporal/10)

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
write.table(LBC_g_regional_results, '/home/jmoodie/Documents/Expression/Scripts_toshare/LBC_g_regional_results.csv', row.names = F)
LBC_g_regional_results_vol <- LBC_g_regional_results[which(substr(LBC_g_regional_results$region, 4, 6) == "vol"),]
LBC_g_regional_results_sa <- LBC_g_regional_results[which(substr(LBC_g_regional_results$region, 4, 5) == "sa"),]
LBC_g_regional_results_thk <- LBC_g_regional_results[which(substr(LBC_g_regional_results$region, 4, 6) == "thk"),]



allbetas <- cbind(LBC_g_regional_results_vol[,2], LBC_g_regional_results_sa[,2], LBC_g_regional_results_thk[,2])
colnames(allbetas) <- c("vol", "sa", "thk")
cor(allbetas)
allbetas <- cbind(STRADL_g_regional_results_vol[,2], STRADL_g_regional_results_sa[,2], STRADL_g_regional_results_thk[,2], UKB_g_regional_results_vol[,2], UKB_g_regional_results_sa[,2], UKB_g_regional_results_thk[,2], LBC_g_regional_results_vol[,2], LBC_g_regional_results_sa[,2], LBC_g_regional_results_thk[,2])
allbetas <- as.data.frame(allbetas)
colnames(allbetas) <- c("vol_STRADL", "sa_STRADL", "thk_STRADL", "vol_UKB", "sa_UKB", "thk_UKB", "vol_LBC", "sa_LBC", "thk_LBC")
cor(allbetas)



# plot the cognitive tests


LBCcogplot <- data.frame(`Age` = theLBCCogdata$AgeMRI/365.25, `Sex` = theLBCCogdata$sex, `WTAR` = theLBCCogdata$WTAR, `NART` = theLBCCogdata$NART, `Verbal fluency` = theLBCCogdata$verbal_fluency, `Logical memory` = theLBCCogdata$logical_memory, `Verbal paired associates` = theLBCCogdata$verbal_paired_associates, `Digit symbol` = theLBCCogdata$digit_symbol, `Symbol search`= theLBCCogdata$symbol_search, `Choice RT` = theLBCCogdata$choice_reaction_time_reflected,  `Inspection time` = theLBCCogdata$inspection_time, `Block design` = theLBCCogdata$block_design,  `Matrix reasoning` = theLBCCogdata$matrix_reasoning, `Spatial span` = theLBCCogdata$spatial_span_total,   `Digit span` = theLBCCogdata$digit_span_backward)

STRADLcogplot <- data.frame(`Age` = theSTRADLCogdata$age, `Sex` = theSTRADLCogdata$sex, `Logical memory` = theSTRADLCogdata$logmem,  `Digit symbol` = theSTRADLCogdata$digsym, `Matrix reasoning` = theSTRADLCogdata$mrtotc, `Verbal fluency` = theSTRADLCogdata$vftot, `Mill Hill` = theSTRADLCogdata$mhv) 


UKBcogplot <- data.frame(`Age` = theUKBCogdata$ageMRI*100, `Sex` = theUKBCogdata$sex, `Paired associates` = theUKBCogdata$cog_pairedAss, `Fluid intelligence` = theUKBCogdata$cog_fluid_intelligence, `Picture vocabulary` = theUKBCogdata$cog_picturevocab, `Prospective memory` = theUKBCogdata$cog_prosmem, `Numeric memory` = theUKBCogdata$cog_numeric_memory,  `Tower` = theUKBCogdata$cog_tower, `Matrix pattern` = theUKBCogdata$cog_matrix_pattern, `Trail B log` = theUKBCogdata$cog_trailB_log, `Digit symbol` = theUKBCogdata$cog_digsym, `RT log` = theUKBCogdata$cog_RT_log,     `Pairs matching log` = theUKBCogdata$cog_pairsmatch_incorrect_log)

pdf('/home/jmoodie/Documents/Expression/Scripts_toshare/chartcorrelation_LBC1936.pdf')
jo_chart.Correlation(LBCcogplot, label.pos = 1.5, cex.size = 1, lwd.size = 0.1)
dev.off()

jpeg('/home/jmoodie/Documents/Expression/Scripts_toshare/chartcorrelationcorr_LBC.jpeg', 800, 800)
corrplot::corrplot(cor(LBCcogplot, use = "complete.obs"), tl.col = "black", method = "number", cl.ratio = 0.1)
dev.off()

pdf('/home/jmoodie/Documents/Expression/Scripts_toshare/chartcorrelation_STRADL.pdf')
jo_chart.Correlation(STRADLcogplot, label.pos = 1.5, cex.size = 1, lwd.size = 0.1)
dev.off()

jpeg('/home/jmoodie/Documents/Expression/Scripts_toshare/chartcorrelationcorr_STRADL.jpeg', 700,700)
corrplot::corrplot(cor(STRADLcogplot, use = "complete.obs"), tl.col = "black", method = "number", cl.ratio = 0.1)
dev.off()

pdf('/home/jmoodie/Documents/Expression/Scripts_toshare/chartcorrelation_UKB.pdf')
jo_chart.Correlation(UKBcogplot, label.pos = 1.5, cex.size = 1, lwd.size = 0.1)
dev.off()

jpeg('/home/jmoodie/Documents/Expression/Scripts_toshare/chartcorrelationcorr_UKB.jpeg', 800, 800)
corrplot::corrplot(cor(UKBcogplot, use = "complete.obs"), tl.col = "black", method = "number", cl.ratio = 0.1)
dev.off()



# Plotting brain by cohort

plotdataUKB = data.frame(lh_vol_bankssts = UKB$lh_bankssts_volume,lh_vol_caudalanteriorcingulate = UKB$lh_caudalanteriorcingulate_volume, lh_vol_caudalmiddlefrontal = UKB$lh_caudalmiddlefrontal_volume, lh_vol_cuneus = UKB$lh_cuneus_volume, lh_vol_entorhinal = UKB$lh_entorhinal_volume, lh_vol_fusiform = UKB$lh_fusiform_volume, lh_vol_inferiorparietal = UKB$lh_inferiorparietal_volume, lh_vol_inferiortemporal = UKB$lh_inferiortemporal_volume, lh_vol_isthmus = UKB$lh_isthmuscingulate_volume, lh_vol_lateraloccipital = UKB$lh_lateraloccipital_volume, lh_vol_lateralorbitofrontal = UKB$lh_lateralorbitofrontal_volume, lh_vol_lingual = UKB$lh_lingual_volume, lh_vol_medialorbitofrontal = UKB$lh_medialorbitofrontal_volume, lh_vol_middletemporal = UKB$lh_middletemporal_volume, lh_vol_parahippocampal= UKB$lh_parahippocampal_volume, lh_vol_paracentral= UKB$lh_paracentral_volume, 
lh_vol_parsopercularis = UKB$lh_parsopercularis_volume, lh_vol_parsorbitalis = UKB$lh_parsorbitalis_volume, 
                      lh_vol_parstriangularis = UKB$lh_parstriangularis_volume, lh_vol_pericalcarine = UKB$lh_pericalcarine_volume, 
                      lh_vol_postcentral = UKB$lh_postcentral_volume, lh_vol_posteriorcingulate= UKB$lh_posteriorcingulate_volume, 
                      lh_vol_precentral = UKB$lh_precentral_volume, lh_vol_precuneus= UKB$lh_precuneus_volume, 
                      lh_vol_rostralanteriorcingulate = UKB$lh_rostralanteriorcingulate_volume, lh_vol_rostralmiddlefrontal = UKB$lh_rostralmiddlefrontal_volume,lh_vol_superiorfrontal = UKB$lh_superiorfrontal_volume, lh_vol_superiorparietal = UKB$lh_superiorparietal_volume, lh_vol_superiortemporal = UKB$lh_superiortemporal_volume, lh_vol_supramarginal = UKB$lh_supramarginal_volume, lh_vol_frontalpole = UKB$lh_frontalpole_volume, 
                      lh_vol_transversetemporal = UKB$lh_transversetemporal_volume, lh_vol_insula = UKB$lh_insula_volume,lh_sa_bankssts = UKB$lh_bankssts_area, lh_sa_caudalanteriorcingulate = UKB$lh_caudalanteriorcingulate_area, 
                      lh_sa_caudalmiddlefrontal = UKB$lh_caudalmiddlefrontal_area, lh_sa_cuneus = UKB$lh_cuneus_area, 
                      lh_sa_entorhinal = UKB$lh_entorhinal_area, 
                      lh_sa_fusiform= UKB$lh_fusiform_area, lh_sa_inferiorparietal = UKB$lh_inferiorparietal_area, 
                      lh_sa_inferiortemporal = UKB$lh_inferiortemporal_area,
                      lh_sa_isthmus = UKB$lh_isthmuscingulate_area, lh_sa_lateraloccipital = UKB$lh_lateraloccipital_area, 
                      lh_sa_lateralorbitofrontal = UKB$lh_lateralorbitofrontal_area, lh_sa_lingual =UKB$lh_lingual_area, 
                      lh_sa_medialorbitofrontal = UKB$lh_medialorbitofrontal_area, lh_sa_middletemporal = UKB$lh_middletemporal_area, 
                      lh_sa_parahippocampal = UKB$lh_parahippocampal_area, lh_sa_paracentral = UKB$lh_paracentral_area, 
                      lh_sa_parsopercularis = UKB$lh_parsopercularis_area,lh_sa_parsorbitalis = UKB$lh_parsorbitalis_area, 
                      lh_sa_parstriangularis = UKB$lh_parstriangularis_area, lh_sa_pericalcarine = UKB$lh_pericalcarine_area, 
                      lh_sa_postcentral = UKB$lh_postcentral_area, lh_sa_posteriorcingulate= UKB$lh_posteriorcingulate_area, 
                      lh_sa_precentral = UKB$lh_precentral_area,lh_sa_precuneus = UKB$lh_precuneus_area, 
                      lh_sa_rostralanteriorcingulate = UKB$lh_rostralanteriorcingulate_area,
                      lh_sa_rostralmiddlefrontal = UKB$lh_rostralmiddlefrontal_area, lh_sa_superiorfrontal = UKB$lh_superiorfrontal_area, 
                      lh_sa_superiorparietal = UKB$lh_superiorparietal_area, lh_sa_superiortemporal = UKB$lh_superiortemporal_area, 
                      lh_sa_supramarginal = UKB$lh_supramarginal_area, lh_sa_frontalpole = UKB$lh_frontalpole_area, 
                      lh_sa_transversetemporal = UKB$lh_transversetemporal_area, lh_sa_insula = UKB$lh_insula_area, 
                      lh_thk_bankssts = UKB$lh_bankssts_thickness,
                      lh_thk_caudalanteriorcingulate = UKB$lh_caudalanteriorcingulate_thickness, lh_thk_caudalmiddlefrontal = UKB$lh_caudalmiddlefrontal_thickness,
                      lh_thk_cuneus = UKB$lh_cuneus_thickness, lh_thk_entorhinal = UKB$lh_entorhinal_thickness, lh_thk_fusiform = UKB$lh_fusiform_thickness, 
                      lh_thk_inferiorparietal = UKB$lh_inferiorparietal_thickness, lh_thk_inferiortemporal = UKB$lh_inferiortemporal_thickness, 
                      lh_thk_isthmus = UKB$lh_isthmuscingulate_thickness, lh_thk_lateraloccipital = UKB$lh_lateraloccipital_thickness, 
                      lh_thk_lateralorbitofrontal = UKB$lh_lateralorbitofrontal_thickness,lh_thk_lingual = UKB$lh_lingual_thickness, 
                      lh_thk_medialorbitofrontal = UKB$lh_medialorbitofrontal_thickness, lh_thk_middletemporal = UKB$lh_middletemporal_thickness,
                      lh_thk_parahippocampal = UKB$lh_parahippocampal_thickness, lh_thk_paracentral = UKB$lh_paracentral_thickness, 
                      lh_thk_parsopercularis = UKB$lh_parsopercularis_thickness, lh_thk_parsorbitalis = UKB$lh_parsorbitalis_thickness, 
                      lh_thk_parstriangularis = UKB$lh_parstriangularis_thickness, lh_thk_pericalcarine = UKB$lh_pericalcarine_thickness, 
                      lh_thk_postcentral = UKB$lh_postcentral_thickness, lh_thk_posteriorcingulate = UKB$lh_posteriorcingulate_thickness, 
                      lh_thk_precentral = UKB$lh_precentral_thickness,lh_thk_precuneus = UKB$lh_precuneus_thickness, 
                      lh_thk_rostralanteriorcingulate = UKB$lh_rostralanteriorcingulate_thickness, lh_thk_rostralmiddlefrontal = UKB$lh_rostralmiddlefrontal_thickness,
                      lh_thk_superiorfrontal = UKB$lh_superiorfrontal_thickness, lh_thk_superiorparietal = UKB$lh_superiorparietal_thickness, 
                      lh_thk_superiortemporal = UKB$lh_superiortemporal_thickness, lh_thk_supramarginal = UKB$lh_supramarginal_thickness, 
                      lh_thk_frontalpole = UKB$lh_frontalpole_thickness,
                      lh_thk_transversetemporal = UKB$lh_transversetemporal_thickness, lh_thk_insula = UKB$lh_insula_thickness,rh_vol_bankssts = UKB$rh_bankssts_volume, rh_vol_caudalanteriorcingulate = UKB$rh_caudalanteriorcingulate_volume, rh_vol_caudalmiddlefrontal = UKB$rh_caudalmiddlefrontal_volume, 
                      rh_vol_cuneus = UKB$rh_cuneus_volume, rh_vol_entorhinal = UKB$rh_entorhinal_volume, rh_vol_fusiform = UKB$rh_fusiform_volume, rh_vol_inferiorparietal = UKB$rh_inferiorparietal_volume, rh_vol_inferiortemporal = UKB$rh_inferiortemporal_volume, rh_vol_isthmus = UKB$rh_isthmuscingulate_volume, rh_vol_lateraloccipital = UKB$rh_lateraloccipital_volume, rh_vol_lateralorbitofrontal = UKB$rh_lateralorbitofrontal_volume, rh_vol_lingual = UKB$rh_lingual_volume, rh_vol_medialorbitofrontal = UKB$rh_medialorbitofrontal_volume, rh_vol_middletemporal = UKB$rh_middletemporal_volume, rh_vol_parahippocampal= UKB$rh_parahippocampal_volume, rh_vol_paracentral= UKB$rh_paracentral_volume, 
                      rh_vol_parsopercularis = UKB$rh_parsopercularis_volume, rh_vol_parsorbitalis = UKB$rh_parsorbitalis_volume, 
                      rh_vol_parstriangularis = UKB$rh_parstriangularis_volume, rh_vol_pericalcarine = UKB$rh_pericalcarine_volume, 
                      rh_vol_postcentral = UKB$rh_postcentral_volume, rh_vol_posteriorcingulate= UKB$rh_posteriorcingulate_volume, 
                      rh_vol_precentral = UKB$rh_precentral_volume, rh_vol_precuneus= UKB$rh_precuneus_volume, 
                      rh_vol_rostralanteriorcingulate = UKB$rh_rostralanteriorcingulate_volume, rh_vol_rostralmiddlefrontal = UKB$rh_rostralmiddlefrontal_volume,rh_vol_superiorfrontal = UKB$rh_superiorfrontal_volume, rh_vol_superiorparietal = UKB$rh_superiorparietal_volume, rh_vol_superiortemporal = UKB$rh_superiortemporal_volume, rh_vol_supramarginal = UKB$rh_supramarginal_volume, rh_vol_frontalpole = UKB$rh_frontalpole_volume, 
                      rh_vol_transversetemporal = UKB$rh_transversetemporal_volume, rh_vol_insula = UKB$rh_insula_volume,rh_sa_bankssts = UKB$rh_bankssts_area, rh_sa_caudalanteriorcingulate = UKB$rh_caudalanteriorcingulate_area, 
                      rh_sa_caudalmiddlefrontal = UKB$rh_caudalmiddlefrontal_area, rh_sa_cuneus = UKB$rh_cuneus_area, 
                      rh_sa_entorhinal = UKB$rh_entorhinal_area, 
                      rh_sa_fusiform= UKB$rh_fusiform_area, rh_sa_inferiorparietal = UKB$rh_inferiorparietal_area, 
                      rh_sa_inferiortemporal = UKB$rh_inferiortemporal_area,
                      rh_sa_isthmus = UKB$rh_isthmuscingulate_area, rh_sa_lateraloccipital = UKB$rh_lateraloccipital_area, 
                      rh_sa_lateralorbitofrontal = UKB$rh_lateralorbitofrontal_area, rh_sa_lingual =UKB$rh_lingual_area, 
                      rh_sa_medialorbitofrontal = UKB$rh_medialorbitofrontal_area, rh_sa_middletemporal = UKB$rh_middletemporal_area, 
                      rh_sa_parahippocampal = UKB$rh_parahippocampal_area, rh_sa_paracentral = UKB$rh_paracentral_area, 
                      rh_sa_parsopercularis = UKB$rh_parsopercularis_area,rh_sa_parsorbitalis = UKB$rh_parsorbitalis_area, 
                      rh_sa_parstriangularis = UKB$rh_parstriangularis_area, rh_sa_pericalcarine = UKB$rh_pericalcarine_area, 
                      rh_sa_postcentral = UKB$rh_postcentral_area, rh_sa_posteriorcingulate= UKB$rh_posteriorcingulate_area, 
                      rh_sa_precentral = UKB$rh_precentral_area,rh_sa_precuneus = UKB$rh_precuneus_area, 
                      rh_sa_rostralanteriorcingulate = UKB$rh_rostralanteriorcingulate_area,
                      rh_sa_rostralmiddlefrontal = UKB$rh_rostralmiddlefrontal_area, rh_sa_superiorfrontal = UKB$rh_superiorfrontal_area, 
                      rh_sa_superiorparietal = UKB$rh_superiorparietal_area, rh_sa_superiortemporal = UKB$rh_superiortemporal_area, 
                      rh_sa_supramarginal = UKB$rh_supramarginal_area, rh_sa_frontalpole = UKB$rh_frontalpole_area, 
                      rh_sa_transversetemporal = UKB$rh_transversetemporal_area, rh_sa_insula = UKB$rh_insula_area, 
                      rh_thk_bankssts = UKB$rh_bankssts_thickness,
                      rh_thk_caudalanteriorcingulate = UKB$rh_caudalanteriorcingulate_thickness, rh_thk_caudalmiddlefrontal = UKB$rh_caudalmiddlefrontal_thickness,
                      rh_thk_cuneus = UKB$rh_cuneus_thickness, rh_thk_entorhinal = UKB$rh_entorhinal_thickness, rh_thk_fusiform = UKB$rh_fusiform_thickness, 
                      rh_thk_inferiorparietal = UKB$rh_inferiorparietal_thickness, rh_thk_inferiortemporal = UKB$rh_inferiortemporal_thickness, 
                      rh_thk_isthmus = UKB$rh_isthmuscingulate_thickness, rh_thk_lateraloccipital = UKB$rh_lateraloccipital_thickness, 
                      rh_thk_lateralorbitofrontal = UKB$rh_lateralorbitofrontal_thickness,rh_thk_lingual = UKB$rh_lingual_thickness, 
                      rh_thk_medialorbitofrontal = UKB$rh_medialorbitofrontal_thickness, rh_thk_middletemporal = UKB$rh_middletemporal_thickness,
                      rh_thk_parahippocampal = UKB$rh_parahippocampal_thickness, rh_thk_paracentral = UKB$rh_paracentral_thickness, 
                      rh_thk_parsopercularis = UKB$rh_parsopercularis_thickness, rh_thk_parsorbitalis = UKB$rh_parsorbitalis_thickness, 
                      rh_thk_parstriangularis = UKB$rh_parstriangularis_thickness, rh_thk_pericalcarine = UKB$rh_pericalcarine_thickness, 
                      rh_thk_postcentral = UKB$rh_postcentral_thickness, rh_thk_posteriorcingulate = UKB$rh_posteriorcingulate_thickness, 
                      rh_thk_precentral = UKB$rh_precentral_thickness,rh_thk_precuneus = UKB$rh_precuneus_thickness, 
                      rh_thk_rostralanteriorcingulate = UKB$rh_rostralanteriorcingulate_thickness, rh_thk_rostralmiddlefrontal = UKB$rh_rostralmiddlefrontal_thickness,
                      rh_thk_superiorfrontal = UKB$rh_superiorfrontal_thickness, rh_thk_superiorparietal = UKB$rh_superiorparietal_thickness, 
                      rh_thk_superiortemporal = UKB$rh_superiortemporal_thickness, rh_thk_supramarginal = UKB$rh_supramarginal_thickness, 
                      rh_thk_frontalpole = UKB$rh_frontalpole_thickness,
                      rh_thk_transversetemporal = UKB$rh_transversetemporal_thickness, rh_thk_insula = UKB$rh_insula_thickness, lh_vol_temporalpole = UKB$lh_temporalpole_volume, rh_vol_temporalpole = UKB$rh_temporalpole_volume,lh_sa_temporalpole = UKB$lh_temporalpole_area, rh_sa_temporalpole = UKB$rh_temporalpole_are,lh_thk_temporalpole = UKB$lh_temporalpole_thickness, rh_thk_temporalpole = UKB$rh_temporalpole_thickness, sex = UKB$sex, age = UKB$ageyears_MRI, ID = UKB$ID)


plotdataSTRADL = data.frame(lh_vol_bankssts = STRADL$lh_bankssts_volume,lh_vol_caudalanteriorcingulate = STRADL$lh_caudalanteriorcingulate_volume, lh_vol_caudalmiddlefrontal = STRADL$lh_caudalmiddlefrontal_volume, lh_vol_cuneus = STRADL$lh_cuneus_volume, lh_vol_entorhinal = STRADL$lh_entorhinal_volume, lh_vol_fusiform = STRADL$lh_fusiform_volume, lh_vol_inferiorparietal = STRADL$lh_inferiorparietal_volume, lh_vol_inferiortemporal = STRADL$lh_inferiortemporal_volume, lh_vol_isthmus = STRADL$lh_isthmuscingulate_volume, lh_vol_lateraloccipital = STRADL$lh_lateraloccipital_volume, lh_vol_lateralorbitofrontal = STRADL$lh_lateralorbitofrontal_volume, lh_vol_lingual = STRADL$lh_lingual_volume, lh_vol_medialorbitofrontal = STRADL$lh_medialorbitofrontal_volume, lh_vol_middletemporal = STRADL$lh_middletemporal_volume, lh_vol_parahippocampal= STRADL$lh_parahippocampal_volume, lh_vol_paracentral= STRADL$lh_paracentral_volume, 
lh_vol_parsopercularis = STRADL$lh_parsopercularis_volume, lh_vol_parsorbitalis = STRADL$lh_parsorbitalis_volume, 
                      lh_vol_parstriangularis = STRADL$lh_parstriangularis_volume, lh_vol_pericalcarine = STRADL$lh_pericalcarine_volume, 
                      lh_vol_postcentral = STRADL$lh_postcentral_volume, lh_vol_posteriorcingulate= STRADL$lh_posteriorcingulate_volume, 
                      lh_vol_precentral = STRADL$lh_precentral_volume, lh_vol_precuneus= STRADL$lh_precuneus_volume, 
                      lh_vol_rostralanteriorcingulate = STRADL$lh_rostralanteriorcingulate_volume, lh_vol_rostralmiddlefrontal = STRADL$lh_rostralmiddlefrontal_volume,lh_vol_superiorfrontal = STRADL$lh_superiorfrontal_volume, lh_vol_superiorparietal = STRADL$lh_superiorparietal_volume, lh_vol_superiortemporal = STRADL$lh_superiortemporal_volume, lh_vol_supramarginal = STRADL$lh_supramarginal_volume, lh_vol_frontalpole = STRADL$lh_frontalpole_volume, 
                      lh_vol_transversetemporal = STRADL$lh_transversetemporal_volume, lh_vol_insula = STRADL$lh_insula_volume,lh_sa_bankssts = STRADL$lh_bankssts_area, lh_sa_caudalanteriorcingulate = STRADL$lh_caudalanteriorcingulate_area, 
                      lh_sa_caudalmiddlefrontal = STRADL$lh_caudalmiddlefrontal_area, lh_sa_cuneus = STRADL$lh_cuneus_area, 
                      lh_sa_entorhinal = STRADL$lh_entorhinal_area, 
                      lh_sa_fusiform= STRADL$lh_fusiform_area, lh_sa_inferiorparietal = STRADL$lh_inferiorparietal_area, 
                      lh_sa_inferiortemporal = STRADL$lh_inferiortemporal_area,
                      lh_sa_isthmus = STRADL$lh_isthmuscingulate_area, lh_sa_lateraloccipital = STRADL$lh_lateraloccipital_area, 
                      lh_sa_lateralorbitofrontal = STRADL$lh_lateralorbitofrontal_area, lh_sa_lingual =STRADL$lh_lingual_area, 
                      lh_sa_medialorbitofrontal = STRADL$lh_medialorbitofrontal_area, lh_sa_middletemporal = STRADL$lh_middletemporal_area, 
                      lh_sa_parahippocampal = STRADL$lh_parahippocampal_area, lh_sa_paracentral = STRADL$lh_paracentral_area, 
                      lh_sa_parsopercularis = STRADL$lh_parsopercularis_area,lh_sa_parsorbitalis = STRADL$lh_parsorbitalis_area, 
                      lh_sa_parstriangularis = STRADL$lh_parstriangularis_area, lh_sa_pericalcarine = STRADL$lh_pericalcarine_area, 
                      lh_sa_postcentral = STRADL$lh_postcentral_area, lh_sa_posteriorcingulate= STRADL$lh_posteriorcingulate_area, 
                      lh_sa_precentral = STRADL$lh_precentral_area,lh_sa_precuneus = STRADL$lh_precuneus_area, 
                      lh_sa_rostralanteriorcingulate = STRADL$lh_rostralanteriorcingulate_area,
                      lh_sa_rostralmiddlefrontal = STRADL$lh_rostralmiddlefrontal_area, lh_sa_superiorfrontal = STRADL$lh_superiorfrontal_area, 
                      lh_sa_superiorparietal = STRADL$lh_superiorparietal_area, lh_sa_superiortemporal = STRADL$lh_superiortemporal_area, 
                      lh_sa_supramarginal = STRADL$lh_supramarginal_area, lh_sa_frontalpole = STRADL$lh_frontalpole_area, 
                      lh_sa_transversetemporal = STRADL$lh_transversetemporal_area, lh_sa_insula = STRADL$lh_insula_area, 
                      lh_thk_bankssts = STRADL$lh_bankssts_thickness,
                      lh_thk_caudalanteriorcingulate = STRADL$lh_caudalanteriorcingulate_thickness, lh_thk_caudalmiddlefrontal = STRADL$lh_caudalmiddlefrontal_thickness,
                      lh_thk_cuneus = STRADL$lh_cuneus_thickness, lh_thk_entorhinal = STRADL$lh_entorhinal_thickness, lh_thk_fusiform = STRADL$lh_fusiform_thickness, 
                      lh_thk_inferiorparietal = STRADL$lh_inferiorparietal_thickness, lh_thk_inferiortemporal = STRADL$lh_inferiortemporal_thickness, 
                      lh_thk_isthmus = STRADL$lh_isthmuscingulate_thickness, lh_thk_lateraloccipital = STRADL$lh_lateraloccipital_thickness, 
                      lh_thk_lateralorbitofrontal = STRADL$lh_lateralorbitofrontal_thickness,lh_thk_lingual = STRADL$lh_lingual_thickness, 
                      lh_thk_medialorbitofrontal = STRADL$lh_medialorbitofrontal_thickness, lh_thk_middletemporal = STRADL$lh_middletemporal_thickness,
                      lh_thk_parahippocampal = STRADL$lh_parahippocampal_thickness, lh_thk_paracentral = STRADL$lh_paracentral_thickness, 
                      lh_thk_parsopercularis = STRADL$lh_parsopercularis_thickness, lh_thk_parsorbitalis = STRADL$lh_parsorbitalis_thickness, 
                      lh_thk_parstriangularis = STRADL$lh_parstriangularis_thickness, lh_thk_pericalcarine = STRADL$lh_pericalcarine_thickness, 
                      lh_thk_postcentral = STRADL$lh_postcentral_thickness, lh_thk_posteriorcingulate = STRADL$lh_posteriorcingulate_thickness, 
                      lh_thk_precentral = STRADL$lh_precentral_thickness,lh_thk_precuneus = STRADL$lh_precuneus_thickness, 
                      lh_thk_rostralanteriorcingulate = STRADL$lh_rostralanteriorcingulate_thickness, lh_thk_rostralmiddlefrontal = STRADL$lh_rostralmiddlefrontal_thickness,
                      lh_thk_superiorfrontal = STRADL$lh_superiorfrontal_thickness, lh_thk_superiorparietal = STRADL$lh_superiorparietal_thickness, 
                      lh_thk_superiortemporal = STRADL$lh_superiortemporal_thickness, lh_thk_supramarginal = STRADL$lh_supramarginal_thickness, 
                      lh_thk_frontalpole = STRADL$lh_frontalpole_thickness,
                      lh_thk_transversetemporal = STRADL$lh_transversetemporal_thickness, lh_thk_insula = STRADL$lh_insula_thickness,rh_vol_bankssts = STRADL$rh_bankssts_volume, rh_vol_caudalanteriorcingulate = STRADL$rh_caudalanteriorcingulate_volume, rh_vol_caudalmiddlefrontal = STRADL$rh_caudalmiddlefrontal_volume, 
                      rh_vol_cuneus = STRADL$rh_cuneus_volume, rh_vol_entorhinal = STRADL$rh_entorhinal_volume, rh_vol_fusiform = STRADL$rh_fusiform_volume, rh_vol_inferiorparietal = STRADL$rh_inferiorparietal_volume, rh_vol_inferiortemporal = STRADL$rh_inferiortemporal_volume, rh_vol_isthmus = STRADL$rh_isthmuscingulate_volume, rh_vol_lateraloccipital = STRADL$rh_lateraloccipital_volume, rh_vol_lateralorbitofrontal = STRADL$rh_lateralorbitofrontal_volume, rh_vol_lingual = STRADL$rh_lingual_volume, rh_vol_medialorbitofrontal = STRADL$rh_medialorbitofrontal_volume, rh_vol_middletemporal = STRADL$rh_middletemporal_volume, rh_vol_parahippocampal= STRADL$rh_parahippocampal_volume, rh_vol_paracentral= STRADL$rh_paracentral_volume, 
                      rh_vol_parsopercularis = STRADL$rh_parsopercularis_volume, rh_vol_parsorbitalis = STRADL$rh_parsorbitalis_volume, 
                      rh_vol_parstriangularis = STRADL$rh_parstriangularis_volume, rh_vol_pericalcarine = STRADL$rh_pericalcarine_volume, 
                      rh_vol_postcentral = STRADL$rh_postcentral_volume, rh_vol_posteriorcingulate= STRADL$rh_posteriorcingulate_volume, 
                      rh_vol_precentral = STRADL$rh_precentral_volume, rh_vol_precuneus= STRADL$rh_precuneus_volume, 
                      rh_vol_rostralanteriorcingulate = STRADL$rh_rostralanteriorcingulate_volume, rh_vol_rostralmiddlefrontal = STRADL$rh_rostralmiddlefrontal_volume,rh_vol_superiorfrontal = STRADL$rh_superiorfrontal_volume, rh_vol_superiorparietal = STRADL$rh_superiorparietal_volume, rh_vol_superiortemporal = STRADL$rh_superiortemporal_volume, rh_vol_supramarginal = STRADL$rh_supramarginal_volume, rh_vol_frontalpole = STRADL$rh_frontalpole_volume, 
                      rh_vol_transversetemporal = STRADL$rh_transversetemporal_volume, rh_vol_insula = STRADL$rh_insula_volume,rh_sa_bankssts = STRADL$rh_bankssts_area, rh_sa_caudalanteriorcingulate = STRADL$rh_caudalanteriorcingulate_area, 
                      rh_sa_caudalmiddlefrontal = STRADL$rh_caudalmiddlefrontal_area, rh_sa_cuneus = STRADL$rh_cuneus_area, 
                      rh_sa_entorhinal = STRADL$rh_entorhinal_area, 
                      rh_sa_fusiform= STRADL$rh_fusiform_area, rh_sa_inferiorparietal = STRADL$rh_inferiorparietal_area, 
                      rh_sa_inferiortemporal = STRADL$rh_inferiortemporal_area,
                      rh_sa_isthmus = STRADL$rh_isthmuscingulate_area, rh_sa_lateraloccipital = STRADL$rh_lateraloccipital_area, 
                      rh_sa_lateralorbitofrontal = STRADL$rh_lateralorbitofrontal_area, rh_sa_lingual =STRADL$rh_lingual_area, 
                      rh_sa_medialorbitofrontal = STRADL$rh_medialorbitofrontal_area, rh_sa_middletemporal = STRADL$rh_middletemporal_area, 
                      rh_sa_parahippocampal = STRADL$rh_parahippocampal_area, rh_sa_paracentral = STRADL$rh_paracentral_area, 
                      rh_sa_parsopercularis = STRADL$rh_parsopercularis_area,rh_sa_parsorbitalis = STRADL$rh_parsorbitalis_area, 
                      rh_sa_parstriangularis = STRADL$rh_parstriangularis_area, rh_sa_pericalcarine = STRADL$rh_pericalcarine_area, 
                      rh_sa_postcentral = STRADL$rh_postcentral_area, rh_sa_posteriorcingulate= STRADL$rh_posteriorcingulate_area, 
                      rh_sa_precentral = STRADL$rh_precentral_area,rh_sa_precuneus = STRADL$rh_precuneus_area, 
                      rh_sa_rostralanteriorcingulate = STRADL$rh_rostralanteriorcingulate_area,
                      rh_sa_rostralmiddlefrontal = STRADL$rh_rostralmiddlefrontal_area, rh_sa_superiorfrontal = STRADL$rh_superiorfrontal_area, 
                      rh_sa_superiorparietal = STRADL$rh_superiorparietal_area, rh_sa_superiortemporal = STRADL$rh_superiortemporal_area, 
                      rh_sa_supramarginal = STRADL$rh_supramarginal_area, rh_sa_frontalpole = STRADL$rh_frontalpole_area, 
                      rh_sa_transversetemporal = STRADL$rh_transversetemporal_area, rh_sa_insula = STRADL$rh_insula_area, 
                      rh_thk_bankssts = STRADL$rh_bankssts_thickness,
                      rh_thk_caudalanteriorcingulate = STRADL$rh_caudalanteriorcingulate_thickness, rh_thk_caudalmiddlefrontal = STRADL$rh_caudalmiddlefrontal_thickness,
                      rh_thk_cuneus = STRADL$rh_cuneus_thickness, rh_thk_entorhinal = STRADL$rh_entorhinal_thickness, rh_thk_fusiform = STRADL$rh_fusiform_thickness, 
                      rh_thk_inferiorparietal = STRADL$rh_inferiorparietal_thickness, rh_thk_inferiortemporal = STRADL$rh_inferiortemporal_thickness, 
                      rh_thk_isthmus = STRADL$rh_isthmuscingulate_thickness, rh_thk_lateraloccipital = STRADL$rh_lateraloccipital_thickness, 
                      rh_thk_lateralorbitofrontal = STRADL$rh_lateralorbitofrontal_thickness,rh_thk_lingual = STRADL$rh_lingual_thickness, 
                      rh_thk_medialorbitofrontal = STRADL$rh_medialorbitofrontal_thickness, rh_thk_middletemporal = STRADL$rh_middletemporal_thickness,
                      rh_thk_parahippocampal = STRADL$rh_parahippocampal_thickness, rh_thk_paracentral = STRADL$rh_paracentral_thickness, 
                      rh_thk_parsopercularis = STRADL$rh_parsopercularis_thickness, rh_thk_parsorbitalis = STRADL$rh_parsorbitalis_thickness, 
                      rh_thk_parstriangularis = STRADL$rh_parstriangularis_thickness, rh_thk_pericalcarine = STRADL$rh_pericalcarine_thickness, 
                      rh_thk_postcentral = STRADL$rh_postcentral_thickness, rh_thk_posteriorcingulate = STRADL$rh_posteriorcingulate_thickness, 
                      rh_thk_precentral = STRADL$rh_precentral_thickness,rh_thk_precuneus = STRADL$rh_precuneus_thickness, 
                      rh_thk_rostralanteriorcingulate = STRADL$rh_rostralanteriorcingulate_thickness, rh_thk_rostralmiddlefrontal = STRADL$rh_rostralmiddlefrontal_thickness,
                      rh_thk_superiorfrontal = STRADL$rh_superiorfrontal_thickness, rh_thk_superiorparietal = STRADL$rh_superiorparietal_thickness, 
                      rh_thk_superiortemporal = STRADL$rh_superiortemporal_thickness, rh_thk_supramarginal = STRADL$rh_supramarginal_thickness, 
                      rh_thk_frontalpole = STRADL$rh_frontalpole_thickness,
                      rh_thk_transversetemporal = STRADL$rh_transversetemporal_thickness, rh_thk_insula = STRADL$rh_insula_thickness, lh_vol_temporalpole = STRADL$lh_temporalpole_volume, rh_vol_temporalpole = STRADL$rh_temporalpole_volume,lh_sa_temporalpole = STRADL$lh_temporalpole_area, rh_sa_temporalpole = STRADL$rh_temporalpole_are,lh_thk_temporalpole = STRADL$lh_temporalpole_thickness, rh_thk_temporalpole = STRADL$rh_temporalpole_thickness, sex = STRADL$Sex, age = STRADL$AgeFaceToFace, ID = STRADL$ID)


########################### load LBC ##########

plotdataLBC = data.frame(lh_vol_bankssts = LBC$lh_vol_bankssts, lh_vol_caudalanteriorcingulate = LBC$lh_vol_caudalanteriorcingulate, lh_vol_caudalmiddlefrontal = LBC$lh_vol_caudalmiddlefrontal, 
                      lh_vol_cuneus = LBC$lh_vol_cuneus, lh_vol_entorhinal = LBC$lh_vol_entorhinal, lh_vol_fusiform = LBC$lh_vol_fusiform, lh_vol_inferiorparietal = LBC$lh_vol_inferiorparietal, lh_vol_inferiortemporal = LBC$lh_vol_inferiortemporal, lh_vol_isthmus = LBC$lh_vol_isthmus, lh_vol_lateraloccipital = LBC$lh_vol_lateraloccipital, lh_vol_lateralorbitofrontal = LBC$lh_vol_lateralorbitofrontal, lh_vol_lingual = LBC$lh_vol_lingual, lh_vol_medialorbitofrontal = LBC$lh_vol_medialorbitofrontal, lh_vol_middletemporal = LBC$lh_vol_middletemporal, lh_vol_parahippocampal= LBC$lh_vol_parahippocampal, lh_vol_paracentral= LBC$lh_vol_paracentral, 
                      lh_vol_parsopercularis = LBC$lh_vol_parsopercularis, lh_vol_parsorbitalis = LBC$lh_vol_parsorbitalis, 
                      lh_vol_parstriangularis = LBC$lh_vol_parstriangularis, lh_vol_pericalcarine = LBC$lh_vol_pericalcarine, 
                      lh_vol_postcentral = LBC$lh_vol_postcentral, lh_vol_posteriorcingulate= LBC$lh_vol_posteriorcingulate, 
                      lh_vol_precentral = LBC$lh_vol_precentral, lh_vol_precuneus= LBC$lh_vol_precuneus, 
                      lh_vol_rostralanteriorcingulate = LBC$lh_vol_rostralanteriorcingulate, lh_vol_rostralmiddlefrontal = LBC$lh_vol_rostralmiddlefrontal,lh_vol_superiorfrontal = LBC$lh_vol_superiorfrontal, lh_vol_superiorparietal = LBC$lh_vol_superiorparietal, lh_vol_superiortemporal = LBC$lh_vol_superiortemporal, lh_vol_supramarginal = LBC$lh_vol_supramarginal, lh_vol_frontalpole = LBC$lh_vol_frontalpole, 
                      lh_vol_transversetemporal = LBC$lh_vol_transversetemporal, lh_vol_insula = LBC$lh_vol_insula,lh_sa_bankssts = LBC$lh_sa_bankssts, lh_sa_caudalanteriorcingulate = LBC$lh_sa_caudalanteriorcingulate, 
                      lh_sa_caudalmiddlefrontal = LBC$lh_sa_caudalmiddlefrontal, lh_sa_cuneus = LBC$lh_sa_cuneus, 
                      lh_sa_entorhinal = LBC$lh_sa_entorhinal, 
                      lh_sa_fusiform= LBC$lh_sa_fusiform, lh_sa_inferiorparietal = LBC$lh_sa_inferiorparietal, 
                      lh_sa_inferiortemporal = LBC$lh_sa_inferiortemporal,
                      lh_sa_isthmus = LBC$lh_sa_isthmus, lh_sa_lateraloccipital = LBC$lh_sa_lateraloccipital, 
                      lh_sa_lateralorbitofrontal = LBC$lh_sa_lateralorbitofrontal, lh_sa_lingual =LBC$lh_sa_lingual, 
                      lh_sa_medialorbitofrontal = LBC$lh_sa_medialorbitofrontal, lh_sa_middletemporal = LBC$lh_sa_middletemporal, 
                      lh_sa_parahippocampal = LBC$lh_sa_parahippocampal, lh_sa_paracentral = LBC$lh_sa_paracentral, 
                      lh_sa_parsopercularis = LBC$lh_sa_parsopercularis,lh_sa_parsorbitalis = LBC$lh_sa_parsorbitalis, 
                      lh_sa_parstriangularis = LBC$lh_sa_parstriangularis, lh_sa_pericalcarine = LBC$lh_sa_pericalcarine, 
                      lh_sa_postcentral = LBC$lh_sa_postcentral, lh_sa_posteriorcingulate= LBC$lh_sa_posteriorcingulate, 
                      lh_sa_precentral = LBC$lh_sa_precentral,lh_sa_precuneus = LBC$lh_sa_precuneus, 
                      lh_sa_rostralanteriorcingulate = LBC$lh_sa_rostralanteriorcingulate,
                      lh_sa_rostralmiddlefrontal = LBC$lh_sa_rostralmiddlefrontal, lh_sa_superiorfrontal = LBC$lh_sa_superiorfrontal, 
                      lh_sa_superiorparietal = LBC$lh_sa_superiorparietal, lh_sa_superiortemporal = LBC$lh_sa_superiortemporal, 
                      lh_sa_supramarginal = LBC$lh_sa_supramarginal, lh_sa_frontalpole = LBC$lh_sa_frontalpole, 
                      lh_sa_transversetemporal = LBC$lh_sa_transversetemporal, lh_sa_insula = LBC$lh_sa_insula, 
                      lh_thk_bankssts = LBC$lh_thk_bankssts,
                      lh_thk_caudalanteriorcingulate = LBC$lh_thk_caudalanteriorcingulate, lh_thk_caudalmiddlefrontal = LBC$lh_thk_caudalmiddlefrontal,
                      lh_thk_cuneus = LBC$lh_thk_cuneus, lh_thk_entorhinal = LBC$lh_thk_entorhinal, lh_thk_fusiform = LBC$lh_thk_fusiform, 
                      lh_thk_inferiorparietal = LBC$lh_thk_inferiorparietal, lh_thk_inferiortemporal = LBC$lh_thk_inferiortemporal, 
                      lh_thk_isthmus = LBC$lh_thk_isthmus, lh_thk_lateraloccipital = LBC$lh_thk_lateraloccipital, 
                      lh_thk_lateralorbitofrontal = LBC$lh_thk_lateralorbitofrontal,lh_thk_lingual = LBC$lh_thk_lingual, 
                      lh_thk_medialorbitofrontal = LBC$lh_thk_medialorbitofrontal, lh_thk_middletemporal = LBC$lh_thk_middletemporal,
                      lh_thk_parahippocampal = LBC$lh_thk_parahippocampal, lh_thk_paracentral = LBC$lh_thk_paracentral, 
                      lh_thk_parsopercularis = LBC$lh_thk_parsopercularis, lh_thk_parsorbitalis = LBC$lh_thk_parsorbitalis, 
                      lh_thk_parstriangularis = LBC$lh_thk_parstriangularis, lh_thk_pericalcarine = LBC$lh_thk_pericalcarine, 
                      lh_thk_postcentral = LBC$lh_thk_postcentral, lh_thk_posteriorcingulate = LBC$lh_thk_posteriorcingulate, 
                      lh_thk_precentral = LBC$lh_thk_precentral,lh_thk_precuneus = LBC$lh_thk_precuneus, 
                      lh_thk_rostralanteriorcingulate = LBC$lh_thk_rostralanteriorcingulate, lh_thk_rostralmiddlefrontal = LBC$lh_thk_rostralmiddlefrontal,
                      lh_thk_superiorfrontal = LBC$lh_thk_superiorfrontal, lh_thk_superiorparietal = LBC$lh_thk_superiorparietal, 
                      lh_thk_superiortemporal = LBC$lh_thk_superiortemporal, lh_thk_supramarginal = LBC$lh_thk_supramarginal, 
                      lh_thk_frontalpole = LBC$lh_thk_frontalpole,
                      lh_thk_transversetemporal = LBC$lh_thk_transversetemporal, lh_thk_insula = LBC$lh_thk_insula, rh_vol_bankssts = LBC$rh_vol_bankssts, rh_vol_caudalanteriorcingulate = LBC$rh_vol_caudalanteriorcingulate, rh_vol_caudalmiddlefrontal = LBC$rh_vol_caudalmiddlefrontal, 
                      rh_vol_cuneus = LBC$rh_vol_cuneus, rh_vol_entorhinal = LBC$rh_vol_entorhinal, rh_vol_fusiform = LBC$rh_vol_fusiform, rh_vol_inferiorparietal = LBC$rh_vol_inferiorparietal, rh_vol_inferiortemporal = LBC$rh_vol_inferiortemporal, rh_vol_isthmus = LBC$rh_vol_isthmus, rh_vol_lateraloccipital = LBC$rh_vol_lateraloccipital, rh_vol_lateralorbitofrontal = LBC$rh_vol_lateralorbitofrontal, rh_vol_lingual = LBC$rh_vol_lingual, rh_vol_medialorbitofrontal = LBC$rh_vol_medialorbitofrontal, rh_vol_middletemporal = LBC$rh_vol_middletemporal, rh_vol_parahippocampal= LBC$rh_vol_parahippocampal, rh_vol_paracentral= LBC$rh_vol_paracentral, 
                      rh_vol_parsopercularis = LBC$rh_vol_parsopercularis, rh_vol_parsorbitalis = LBC$rh_vol_parsorbitalis, 
                      rh_vol_parstriangularis = LBC$rh_vol_parstriangularis, rh_vol_pericalcarine = LBC$rh_vol_pericalcarine, 
                      rh_vol_postcentral = LBC$rh_vol_postcentral, rh_vol_posteriorcingulate= LBC$rh_vol_posteriorcingulate, 
                      rh_vol_precentral = LBC$rh_vol_precentral, rh_vol_precuneus= LBC$rh_vol_precuneus, 
                      rh_vol_rostralanteriorcingulate = LBC$rh_vol_rostralanteriorcingulate, rh_vol_rostralmiddlefrontal = LBC$rh_vol_rostralmiddlefrontal,rh_vol_superiorfrontal = LBC$rh_vol_superiorfrontal, rh_vol_superiorparietal = LBC$rh_vol_superiorparietal, rh_vol_superiortemporal = LBC$rh_vol_superiortemporal, rh_vol_supramarginal = LBC$rh_vol_supramarginal, rh_vol_frontalpole = LBC$rh_vol_frontalpole, 
                      rh_vol_transversetemporal = LBC$rh_vol_transversetemporal, rh_vol_insula = LBC$rh_vol_insula,rh_sa_bankssts = LBC$rh_sa_bankssts, rh_sa_caudalanteriorcingulate = LBC$rh_sa_caudalanteriorcingulate, 
                      rh_sa_caudalmiddlefrontal = LBC$rh_sa_caudalmiddlefrontal, rh_sa_cuneus = LBC$rh_sa_cuneus, 
                      rh_sa_entorhinal = LBC$rh_sa_entorhinal, 
                      rh_sa_fusiform= LBC$rh_sa_fusiform, rh_sa_inferiorparietal = LBC$rh_sa_inferiorparietal, 
                      rh_sa_inferiortemporal = LBC$rh_sa_inferiortemporal,
                      rh_sa_isthmus = LBC$rh_sa_isthmus, rh_sa_lateraloccipital = LBC$rh_sa_lateraloccipital, 
                      rh_sa_lateralorbitofrontal = LBC$rh_sa_lateralorbitofrontal, rh_sa_lingual =LBC$rh_sa_lingual, 
                      rh_sa_medialorbitofrontal = LBC$rh_sa_medialorbitofrontal, rh_sa_middletemporal = LBC$rh_sa_middletemporal, 
                      rh_sa_parahippocampal = LBC$rh_sa_parahippocampal, rh_sa_paracentral = LBC$rh_sa_paracentral, 
                      rh_sa_parsopercularis = LBC$rh_sa_parsopercularis,rh_sa_parsorbitalis = LBC$rh_sa_parsorbitalis, 
                      rh_sa_parstriangularis = LBC$rh_sa_parstriangularis, rh_sa_pericalcarine = LBC$rh_sa_pericalcarine, 
                      rh_sa_postcentral = LBC$rh_sa_postcentral, rh_sa_posteriorcingulate= LBC$rh_sa_posteriorcingulate, 
                      rh_sa_precentral = LBC$rh_sa_precentral,rh_sa_precuneus = LBC$rh_sa_precuneus, 
                      rh_sa_rostralanteriorcingulate = LBC$rh_sa_rostralanteriorcingulate,
                      rh_sa_rostralmiddlefrontal = LBC$rh_sa_rostralmiddlefrontal, rh_sa_superiorfrontal = LBC$rh_sa_superiorfrontal, 
                      rh_sa_superiorparietal = LBC$rh_sa_superiorparietal, rh_sa_superiortemporal = LBC$rh_sa_superiortemporal, 
                      rh_sa_supramarginal = LBC$rh_sa_supramarginal, rh_sa_frontalpole = LBC$rh_sa_frontalpole, 
                      rh_sa_transversetemporal = LBC$rh_sa_transversetemporal, rh_sa_insula = LBC$rh_sa_insula, 
                      rh_thk_bankssts = LBC$rh_thk_bankssts,
                      rh_thk_caudalanteriorcingulate = LBC$rh_thk_caudalanteriorcingulate, rh_thk_caudalmiddlefrontal = LBC$rh_thk_caudalmiddlefrontal,
                      rh_thk_cuneus = LBC$rh_thk_cuneus, rh_thk_entorhinal = LBC$rh_thk_entorhinal, rh_thk_fusiform = LBC$rh_thk_fusiform, 
                      rh_thk_inferiorparietal = LBC$rh_thk_inferiorparietal, rh_thk_inferiortemporal = LBC$rh_thk_inferiortemporal, 
                      rh_thk_isthmus = LBC$rh_thk_isthmus, rh_thk_lateraloccipital = LBC$rh_thk_lateraloccipital, 
                      rh_thk_lateralorbitofrontal = LBC$rh_thk_lateralorbitofrontal,rh_thk_lingual = LBC$rh_thk_lingual, 
                      rh_thk_medialorbitofrontal = LBC$rh_thk_medialorbitofrontal, rh_thk_middletemporal = LBC$rh_thk_middletemporal,
                      rh_thk_parahippocampal = LBC$rh_thk_parahippocampal, rh_thk_paracentral = LBC$rh_thk_paracentral, 
                      rh_thk_parsopercularis = LBC$rh_thk_parsopercularis, rh_thk_parsorbitalis = LBC$rh_thk_parsorbitalis, 
                      rh_thk_parstriangularis = LBC$rh_thk_parstriangularis, rh_thk_pericalcarine = LBC$rh_thk_pericalcarine, 
                      rh_thk_postcentral = LBC$rh_thk_postcentral, rh_thk_posteriorcingulate = LBC$rh_thk_posteriorcingulate, 
                      rh_thk_precentral = LBC$rh_thk_precentral,rh_thk_precuneus = LBC$rh_thk_precuneus, 
                      rh_thk_rostralanteriorcingulate = LBC$rh_thk_rostralanteriorcingulate, rh_thk_rostralmiddlefrontal = LBC$rh_thk_rostralmiddlefrontal,
                      rh_thk_superiorfrontal = LBC$rh_thk_superiorfrontal, rh_thk_superiorparietal = LBC$rh_thk_superiorparietal, 
                      rh_thk_superiortemporal = LBC$rh_thk_superiortemporal, rh_thk_supramarginal = LBC$rh_thk_supramarginal, 
                      rh_thk_frontalpole = LBC$rh_thk_frontalpole, 
                      rh_thk_transversetemporal = LBC$rh_thk_transversetemporal, rh_thk_insula = LBC$rh_thk_insula, lh_vol_temporalpole = LBC$lh_vol_temporalpole, rh_vol_temporalpole = LBC$rh_vol_temporalpole,lh_sa_temporalpole = LBC$lh_sa_temporalpole, rh_sa_temporalpole = LBC$rh_sa_temporalpole,lh_thk_temporalpole = LBC$lh_thk_temporalpole, rh_thk_temporalpole = LBC$rh_thk_temporalpole,sex = LBC$sex, age = LBC$AgeMRI/365.25, ID = LBC$ID)


# merge 3 cohorts for plots

dim(plotdataUKB)
dim(plotdataSTRADL)
dim(plotdataLBC) # all have the same columns in the same order!


plot3data <- rbind(plotdataUKB, plotdataSTRADL, plotdataLBC)
dim(plot3data) # just to check
plot3data$cohort <- c(rep("UKB", dim(plotdataUKB)[1]), rep("STRADL", dim(plotdataSTRADL)[1]), rep("LBC1936", dim(plotdataLBC)[1]))
plot3data$cohort <- as.factor(plot3data$cohort)

regionlabelsrev <- rev(c("Bank ssts", "Caudal anterior cingulate", "Caudal middle frontal", "Cuneus", "Entorhinal",  "Fusiform", "Inferior parietal", "Inferior temporal", "Isthmus cingulate", "Lateral occipital", "Lateral orbitofrontal", "Lingual", "Middle temporal", "Medial orbitofrontal", "Paracentral", "Parahippocampal", "Pars opercularis", "Pars orbitalis", "Pars triangularis", "Pericalcarine", "Postcentral", "Posterior cingulate", "Precentral", "Precuneus", "Rostral anterior cingulate", "Rostral middle frontal", "Superior frontal", "Superior parietal", "Superior temporal", "Supramarginal","Frontal pole",  "Transverse temporal", "Insula", "Temporal pole"))


plotvol3_lh <- plot3data[ ,c(1:33, 199, 206, 207, 208)] # check this 

plotvol3_lh <-   melt(plot3data[ ,c(1:33, 199, 206, 207, 208)], id.vars = c("age","ID", "cohort")) 
plotvol3_lh$variable <- as.character(plotvol3_lh$variable)

plotvol3_lh$variable <- gsub("lh_vol_bankssts", "Bank ssts", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_caudalanteriorcingulate", "Caudal anterior cingulate", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_caudalmiddlefront", "Caudal middle front", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_cuneus", "Cuneus", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_entorhin", "Entorhin", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_fusi", "Fusi", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_inferiorparietal", "Inferior parietal", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_inferiortemp", "Inferior temp", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_insul", "Insul", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_isthmus", "Isthmus cingulate", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_lateraloccipital", "Lateral occipital", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_lateralorbito", "Lateral orbito", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_lingual", "Lingual", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_medialorbi", "Medial orbi", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_middletemporal", "Middle temporal", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_paracent", "Paracent", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_parahipp", "Parahipp", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_parsopercularis", "Pars opercularis", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_parsorb", "Pars orb", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_parstriang", "Pars triang", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_perical", "Perical", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_postcent", "Postcent", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_posteriorcingulate", "Posterior cingulate", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_precent", "Precent", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_precuneu", "Precuneu", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_rostralanteriorcingulate", "Rostral anterior cingulate", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_rostralmiddlefro", "Rostral middle fro", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_superiorfro", "Superior fro", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_superiorpari", "Superior pari", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_superiortemp", "Superior temp", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_supra", "Supra", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_caudalmiddlefront", "Caudal middle front", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_temporalpol", "Temporal pol", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_frontalpol", "Frontal pol", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_temporalpol", "Temporal pol", plotvol3_lh$variable)
plotvol3_lh$variable <- gsub("lh_vol_transversetemp", "Transverse temp", plotvol3_lh$variable)
plotthk3_lh$variable <- gsub("lh_vol_temporalpol", "Temporal pol, plotthk3_lh$variable)




plotvol3_lh$variable <- as.factor(plotvol3_lh$variable)
levels(plotvol3_lh$variable) <-levels(plotvol3_lh$variable[order(levels(plotvol3_lh$variable))])





ggplot(data = plotvol3_lh, aes(x = variable, y = value, color = cohort)) + geom_flat_violin(position = position_nudge(x = 0, y = 0), width = 1.8, size = 0.3, alpha = 0.2) + xlab("Region") + ylab(expression("Volume, mm" ^3)) +coord_flip() +theme_cowplot() + theme( text = element_text(size=16), axis.text.y = element_text(size = 11), legend.title=element_blank()) +labs(subtitle="Left Hemisphere") + scale_y_continuous(breaks = c(0, 20000, 40000), limits = c(0,40000)) + scale_x_discrete(limits = rev(levels(plotvol3_lh$variable))) + scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE))
ggsave("/home/jmoodie/Documents/Expression/Scripts_toshare/plot_lhvol.jpeg", bg = "white", width = 8, height = 11)



###############


plotvol3_rh <- plot3data[ ,c(100:132, 200, 206, 207, 208)] # can check this

plotvol3_rh <-   melt(plot3data[ ,c(100:132, 200, 206, 207, 208)], id.vars = c("age", "ID", "cohort")) 
plotvol3_rh$variable <- as.character(plotvol3_rh$variable)

plotvol3_rh$variable <- gsub("rh_vol_bankssts", "Bank ssts", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_caudalanteriorcingulate", "Caudal anterior cingulate", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_caudalmiddlefront", "Caudal middle front", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_cuneus", "Cuneus", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_entorhin", "Entorhin", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_fusi", "Fusi", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_inferiorparietal", "Inferior parietal", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_inferiortemp", "Inferior temp", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_insul", "Insul", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_isthmus", "Isthmus cingulate", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_lateraloccipital", "Lateral occipital", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_lateralorbito", "Lateral orbito", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_lingual", "Lingual", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_medialorbi", "Medial orbi", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_middletemporal", "Middle temporal", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_paracent", "Paracent", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_parahipp", "Parahipp", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_parsopercularis", "Pars opercularis", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_parsorb", "Pars orb", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_parstriang", "Pars triang", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_perical", "Perical", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_postcent", "Postcent", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_posteriorcingulate", "Posterior cingulate", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_precent", "Precent", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_precuneu", "Precuneu", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_rostralanteriorcingulate", "Rostral anterior cingulate", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_rostralmiddlefro", "Rostral middle fro", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_superiorfro", "Superior fro", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_superiorpari", "Superior pari", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_superiortemp", "Superior temp", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_supra", "Supra", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_caudalmiddlefront", "Caudal middle front", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_temporalpol", "Temporal pol", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_frontalpol", "Frontal pol", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_temporalpol", "Temporal pol", plotvol3_rh$variable)
plotvol3_rh$variable <- gsub("rh_vol_transversetemp", "Transverse temp", plotvol3_rh$variable)



plotvol3_rh$variable <- as.factor(plotvol3_rh$variable)
levels(plotvol3_rh$variable) <- levels(plotvol3_rh$variable)[order(levels(plotvol3_rh$variable))]


ggplot(data = plotvol3_rh, aes(x = variable, y = value, color = cohort)) + geom_flat_violin(position = position_nudge(x = 0, y = 0), width = 1.8, size = 0.3, alpha = 0.2) + xlab("Region") +ylab(expression("Volume, mm" ^3)) +coord_flip() +theme_cowplot() + theme( text = element_text(size=16), axis.text.y = element_text(size = 11), legend.title=element_blank()) +labs(subtitle="Right Hemisphere") + scale_y_continuous(breaks = c(0, 20000, 40000), limits = c(0,40000)) + scale_x_discrete(limits = rev(levels(plotvol3_rh$variable))) + scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) 
ggsave("/home/jmoodie/Documents/Expression/Scripts_toshare/plot_rhvol.jpeg", bg = "white", width = 8, height = 11)

###############

plotthk3_lh <- plot3data[ ,c(67:99, 203, 206, 207, 208)] # can check this

plotthk3_lh <-   melt(plot3data[ ,c(67:99, 203, 206, 207, 208)], id.vars = c("age", "ID", "cohort")) 
plotthk3_lh$variable <- as.character(plotthk3_lh$variable)

plotthk3_lh$variable <- gsub("lh_thk_bankssts", "Bank ssts", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_caudalanteriorcingulate", "Caudal anterior cingulate", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_caudalmiddlefront", "Caudal middle front", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_cuneus", "Cuneus", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_entolhin", "Entolhin", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_fusi", "Fusi", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_inferiorparietal", "Inferior parietal", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_inferiortemp", "Inferior temp", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_insul", "Insul", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_isthmus", "Isthmus cingulate", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_lateraloccipital", "Lateral occipital", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_lateralorbito", "Lateral orbito", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_lingual", "Lingual", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_medialorbi", "Medial orbi", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_middletemporal", "Middle temporal", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_paracent", "Paracent", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_parahipp", "Parahipp", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_parsopercularis", "Pars opercularis", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_parsorb", "Pars orb", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_parstriang", "Pars triang", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_perical", "Perical", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_postcent", "Postcent", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_posteriorcingulate", "Posterior cingulate", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_precent", "Precent", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_precuneu", "Precuneu", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_rostralanteriorcingulate", "Rostral anterior cingulate", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_rostralmiddlefro", "Rostral middle fro", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_superiorfro", "Superior fro", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_superiorpari", "Superior pari", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_superiortemp", "Superior temp", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_supra", "Supra", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_caudalmiddlefront", "Caudal middle front", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_temporalpol", "Temporal pol", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_frontalpol", "Frontal pol", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_temporalpol", "Temporal pol", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_transversetemp", "Transverse temp", plotthk3_lh$variable)
plotthk3_lh$variable <- gsub("lh_thk_temporalpol", "Temporal pol, plotthk3_lh$variable)


plotthk3_lh$variable <- as.factor(plotthk3_lh$variable)
levels(plotthk3_lh$variable) <- levels(plotthk3_lh$variable)[order(levels(plotthk3_lh$variable))]



ggplot(data = plotthk3_lh, aes(x = variable, y = value, color = cohort)) + geom_flat_violin(position = position_nudge(x = 0, y = 0), width = 1.8, size = 0.3, alpha = 0.2) + xlab("Region") +ylab(expression("Thickness, mm")) +coord_flip() +theme_cowplot() + theme( text = element_text(size=16), axis.text.y = element_text(size = 11), legend.title=element_blank()) +labs(subtitle="Left Hemisphere") + scale_y_continuous(breaks = c(1, 2, 3, 4, 5), limits = c(0,5)) + scale_x_discrete(limits = rev(levels(plotthk3_lh$variable))) + scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) 
ggsave("/home/jmoodie/Documents/Expression/Scripts_toshare/plot_lhthk.jpeg", bg = "white", width = 8, height = 11)




#########################################

plotthk3_rh <- plot3data[ ,c(166:198, 204, 206, 207, 208)] # can check this

plotthk3_rh <-   melt(plot3data[ ,c(166:198, 204, 206, 207, 208)], id.vars = c("age", "ID", "cohort")) 
plotthk3_rh$variable <- as.character(plotthk3_rh$variable)

plotthk3_rh$variable <- gsub("rh_thk_bankssts", "Bank ssts", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_caudalanteriorcingulate", "Caudal anterior cingulate", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_caudalmiddlefront", "Caudal middle front", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_cuneus", "Cuneus", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_entorhin", "Entorhin", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_fusi", "Fusi", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_inferiorparietal", "Inferior parietal", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_inferiortemp", "Inferior temp", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_insul", "Insul", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_isthmus", "Isthmus cingulate", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_lateraloccipital", "Lateral occipital", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_lateralorbito", "Lateral orbito", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_lingual", "Lingual", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_medialorbi", "Medial orbi", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_middletemporal", "Middle temporal", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_paracent", "Paracent", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_parahipp", "Parahipp", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_parsopercularis", "Pars opercularis", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_parsorb", "Pars orb", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_parstriang", "Pars triang", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_perical", "Perical", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_postcent", "Postcent", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_posteriorcingulate", "Posterior cingulate", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_precent", "Precent", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_precuneu", "Precuneu", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_rostralanteriorcingulate", "Rostral anterior cingulate", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_rostralmiddlefro", "Rostral middle fro", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_superiorfro", "Superior fro", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_superiorpari", "Superior pari", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_superiortemp", "Superior temp", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_supra", "Supra", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_caudalmiddlefront", "Caudal middle front", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_temporalpol", "Temporal pol", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_frontalpol", "Frontal pol", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_temporalpol", "Temporal pol", plotthk3_rh$variable)
plotthk3_rh$variable <- gsub("rh_thk_transversetemp", "Transverse temp", plotthk3_rh$variable)



plotthk3_rh$variable <- as.factor(plotthk3_rh$variable)
levels(plotthk3_rh$variable) <- levels(plotthk3_rh$variable)[order(levels(plotthk3_rh$variable))]



ggplot(data = plotthk3_rh, aes(x = variable, y = value, color = cohort)) + geom_flat_violin(position = position_nudge(x = 0, y = 0), width = 1.8, size = 0.3, alpha = 0.2) + xlab("Region") +ylab(expression("Thickness, mm")) +coord_flip() +theme_cowplot() + theme( text = element_text(size=16), axis.text.y = element_text(size = 11), legend.title=element_blank()) +labs(subtitle="Right Hemisphere") + scale_y_continuous(breaks = c(1, 2, 3, 4, 5), limits = c(0,5)) + scale_x_discrete(limits = rev(levels(plotthk3_rh$variable))) + scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) 
ggsave("/home/jmoodie/Documents/Expression/Scripts_toshare/plot_rhthk.jpeg", bg = "white", width = 8, height = 11)

###############



plotsa3_lh <- plot3data[ ,c(34:66, 201, 206, 207, 208)] # can check this

plotsa3_lh <-   melt(plot3data[ ,c(34:66, 201, 206, 207, 208)], id.vars = c("age", "ID", "cohort")) 
plotsa3_lh$variable <- as.character(plotsa3_lh$variable)

plotsa3_lh$variable <- gsub("lh_sa_bankssts", "Bank ssts", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_caudalanteriorcingulate", "Caudal anterior cingulate", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_caudalmiddlefront", "Caudal middle front", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_cuneus", "Cuneus", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_entolhin", "Entolhin", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_fusi", "Fusi", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_inferiorparietal", "Inferior parietal", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_inferiortemp", "Inferior temp", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_insul", "Insul", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_isthmus", "Isthmus cingulate", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_lateraloccipital", "Lateral occipital", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_lateralorbito", "Lateral orbito", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_lingual", "Lingual", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_medialorbi", "Medial orbi", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_middletemporal", "Middle temporal", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_paracent", "Paracent", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_parahipp", "Parahipp", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_parsopercularis", "Pars opercularis", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_parsorb", "Pars orb", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_parstriang", "Pars triang", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_perical", "Perical", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_postcent", "Postcent", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_posteriorcingulate", "Posterior cingulate", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_precent", "Precent", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_precuneu", "Precuneu", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_rostralanteriorcingulate", "Rostral anterior cingulate", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_rostralmiddlefro", "Rostral middle fro", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_superiorfro", "Superior fro", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_superiorpari", "Superior pari", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_superiortemp", "Superior temp", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_supra", "Supra", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_caudalmiddlefront", "Caudal middle front", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_temporalpol", "Temporal pol", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_frontalpol", "Frontal pol", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_temporalpol", "Temporal pol", plotsa3_lh$variable)
plotsa3_lh$variable <- gsub("lh_sa_transversetemp", "Transverse temp", plotsa3_lh$variable)



plotsa3_lh$variable <- as.factor(plotsa3_lh$variable)
levels(plotsa3_lh$variable) <- levels(plotsa3_lh$variable)[order(levels(plotsa3_lh$variable))]



ggplot(data = plotsa3_lh, aes(x = variable, y = value, color = cohort)) + geom_flat_violin(position = position_nudge(x = 0, y = 0), width = 1.8, size = 0.3, alpha = 0.2) + xlab("Region") + ylab(expression("Surface area, mm" ^2))  +coord_flip() +theme_cowplot() + theme( text = element_text(size=16), axis.text.y = element_text(size = 11), legend.title=element_blank()) +labs(subtitle="Left Hemisphere") + scale_y_continuous(breaks = c(0, 4000, 8000, 12000), limits = c(0,12000)) + scale_x_discrete(limits = rev(levels(plotsa3_lh$variable))) + scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) 
ggsave("/home/jmoodie/Documents/Expression/Scripts_toshare/plot_lhsa.jpeg", bg = "white", width = 8, height = 11)

###############


plotsa3_rh <- plot3data[ ,c(133:165, 202, 206, 207, 208)] # can check this

plotsa3_rh <-   melt(plot3data[ ,c(133:165, 202, 206, 207, 208)], id.vars = c("age", "ID", "cohort")) 
plotsa3_rh$variable <- as.character(plotsa3_rh$variable)

plotsa3_rh$variable <- gsub("rh_sa_bankssts", "Bank ssts", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_caudalanteriorcingulate", "Caudal anterior cingulate", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_caudalmiddlefront", "Caudal middle front", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_cuneus", "Cuneus", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_entorhin", "Entorhin", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_fusi", "Fusi", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_inferiorparietal", "Inferior parietal", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_inferiortemp", "Inferior temp", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_insul", "Insul", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_isthmus", "Isthmus cingulate", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_lateraloccipital", "Lateral occipital", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_lateralorbito", "Lateral orbito", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_lingual", "Lingual", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_medialorbi", "Medial orbi", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_middletemporal", "Middle temporal", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_paracent", "Paracent", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_parahipp", "Parahipp", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_parsopercularis", "Pars opercularis", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_parsorb", "Pars orb", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_parstriang", "Pars triang", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_perical", "Perical", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_postcent", "Postcent", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_posteriorcingulate", "Posterior cingulate", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_precent", "Precent", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_precuneu", "Precuneu", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_rostralanteriorcingulate", "Rostral anterior cingulate", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_rostralmiddlefro", "Rostral middle fro", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_superiorfro", "Superior fro", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_superiorpari", "Superior pari", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_superiortemp", "Superior temp", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_supra", "Supra", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_caudalmiddlefront", "Caudal middle front", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_temporalpol", "Temporal pol", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_frontalpol", "Frontal pol", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_temporalpol", "Temporal pol", plotsa3_rh$variable)
plotsa3_rh$variable <- gsub("rh_sa_transversetemp", "Transverse temp", plotsa3_rh$variable)



plotsa3_rh$variable <- as.factor(plotsa3_rh$variable)
levels(plotsa3_rh$variable) <- levels(plotsa3_rh$variable)[order(levels(plotsa3_rh$variable))]




ggplot(data = plotsa3_rh, aes(x = variable, y = value, color = cohort)) + geom_flat_violin(position = position_nudge(x = 0, y = 0), width = 1.8, size = 0.3, alpha = 0.2) + xlab("Region") + ylab(expression("Surface area, mm" ^2))  +coord_flip() +theme_cowplot() + theme( text = element_text(size=16), axis.text.y = element_text(size = 11), legend.title=element_blank()) +labs(subtitle="Right Hemisphere") + scale_y_continuous(breaks = c(0, 4000, 8000, 12000), limits = c(0,12000)) + scale_x_discrete(limits = rev(levels(plotsa3_rh$variable))) + scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) 
ggsave("/home/jmoodie/Documents/Expression/Scripts_toshare/plot_rhsa.jpeg", bg = "white", width = 8, height = 11)


UKBtotvol <- rowSums(plotdataUKB[,c(1:33, 100:132, 199:200)])
UKBtotvol <- as.data.frame(UKBtotvol)
UKBtotvol$age <- plotdataUKB$age
UKBtotvol$ID <- plotdataUKB$ID
colnames(UKBtotvol) <- c("Total Volume", "Age", "ID")



STRADLtotvol <- rowSums(plotdataSTRADL[,c(1:33, 100:132, 199:200)])
STRADLtotvol <- as.data.frame(STRADLtotvol)
STRADLtotvol$age <- plotdataSTRADL$age
STRADLtotvol$ID <- plotdataSTRADL$ID
colnames(STRADLtotvol) <- c("Total Volume", "Age", "ID")


LBCtotvol <- rowSums(plotdataLBC[,c(1:33, 100:132, 199:200)])
LBCtotvol <- as.data.frame(LBCtotvol)
LBCtotvol$age <- plotdataLBC$age
LBCtotvol$ID <- plotdataLBC$ID
colnames(LBCtotvol) <- c("Total Volume", "Age", "ID")

cohort <- c(rep("UKB", dim(UKBtotvol)[1]), rep("STRADL", dim(STRADLtotvol)[1]), rep("LBC", dim(LBCtotvol)[1]))

totvol <- rbind(UKBtotvol, STRADLtotvol, LBCtotvol)
totvol <- as.data.frame(totvol)
totvol$cohort <- cohort
totvol$cohort <- as.factor(totvol$cohort)
totvol$ID <- as.factor(totvol$ID)


ggplot() +  geom_point(data = totvol, aes(y = `Total Volume`, x=Age, color = cohort), size = .2, alpha = .06) + geom_density_2d(data = totvol[which(totvol$cohort == "LBC"),], aes(y = `Total Volume`, x=Age, color = cohort))+ geom_density_2d(data = totvol[which(totvol$cohort == "STRADL"),], aes(y = `Total Volume`, x=Age, color = cohort))+ geom_density_2d(data = totvol[which(totvol$cohort == "UKB"),], aes(y = `Total Volume`, x=Age, color = cohort))  + ylab(expression("Total volume, mm" ^3)) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size=16),  plot.margin=unit(c(1, 1, 0, 1), "cm"))  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + scale_y_continuous(breaks = c(300000, 400000, 500000, 600000, 700000))
ggsave(filename = "/home/jmoodie/Documents/Expression/Scripts_toshare/plot_agetotvol.jpeg", bg = "white",  width = 8.8, height = 7)





UKBtotsa <- rowSums(plotdataUKB[,c(34:66, 133:165, 201:202)])
UKBtotsa <- as.data.frame(UKBtotsa)
UKBtotsa$age <- plotdataUKB$age
UKBtotsa$ID <- plotdataUKB$ID
colnames(UKBtotsa) <- c("Total area", "Age", "ID")



STRADLtotsa <- rowSums(plotdataSTRADL[,c(34:66, 133:165, 201:202)])
STRADLtotsa <- as.data.frame(STRADLtotsa)
STRADLtotsa$age <- plotdataSTRADL$age
STRADLtotsa$ID <- plotdataSTRADL$ID
colnames(STRADLtotsa) <- c("Total area", "Age", "ID")


LBCtotsa <- rowSums(plotdataLBC[,c(34:66, 133:165, 201:202)])
LBCtotsa <- as.data.frame(LBCtotsa)
LBCtotsa$age <- plotdataLBC$age
LBCtotsa$ID <- plotdataLBC$ID
colnames(LBCtotsa) <- c("Total area", "Age", "ID")

cohort <- c(rep("UKB", dim(UKBtotsa)[1]), rep("STRADL", dim(STRADLtotsa)[1]), rep("LBC", dim(LBCtotsa)[1]))

totsa <- rbind(UKBtotsa, STRADLtotsa, LBCtotsa)
totsa <- as.data.frame(totsa)
totsa$cohort <- cohort
totsa$cohort <- as.factor(totsa$cohort)
totsa$ID <- as.factor(totsa$ID)



ggplot() + geom_point(data = totsa, aes(y = `Total area`, x=Age, color = cohort), size = .2, alpha = .06) + geom_density_2d(data = totsa[which(totsa$cohort == "LBC"),], aes(y = `Total area`, x=Age, color = cohort))+ geom_density_2d(data = totsa[which(totsa$cohort == "STRADL"),], aes(y = `Total area`, x=Age, color = cohort))+ geom_density_2d(data = totsa[which(totsa$cohort == "UKB"),], aes(y = `Total area`, x=Age, color = cohort)) +  ylab(expression("Total surface area, mm" ^2)) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size=16),  plot.margin=unit(c(1, 1, 0, 1), "cm"))  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + scale_y_continuous(breaks = c(120000, 140000, 160000, 180000, 200000, 220000))
ggsave(filename = "/home/jmoodie/Documents/Expression/Scripts_toshare/plot_agetotsa.jpeg", bg = "white",  width = 8.8, height = 7)



UKBmeanthk <- rowMeans(plotdataUKB[,c(67:99, 166:198, 203:204)])
UKBmeanthk <- as.data.frame(UKBmeanthk)
UKBmeanthk$age <- plotdataUKB$age
UKBmeanthk$ID <- plotdataUKB$ID
colnames(UKBmeanthk) <- c("mean thickness", "Age", "ID")



STRADLmeanthk <- rowMeans(plotdataSTRADL[,c(67:99, 166:198, 203:204)])
STRADLmeanthk <- as.data.frame(STRADLmeanthk)
STRADLmeanthk$age <- plotdataSTRADL$age
STRADLmeanthk$ID <- plotdataSTRADL$ID
colnames(STRADLmeanthk) <- c("mean thickness", "Age", "ID")


LBCmeanthk <- rowMeans(plotdataLBC[,c(67:99, 166:198, 203:204)])
LBCmeanthk <- as.data.frame(LBCmeanthk)
LBCmeanthk$age <- plotdataLBC$age
LBCmeanthk$ID <- plotdataLBC$ID
colnames(LBCmeanthk) <- c("mean thickness", "Age", "ID")

cohort <- c(rep("UKB", dim(UKBmeanthk)[1]), rep("STRADL", dim(STRADLmeanthk)[1]), rep("LBC", dim(LBCmeanthk)[1]))

meanthk <- rbind(UKBmeanthk, STRADLmeanthk, LBCmeanthk)
meanthk <- as.data.frame(meanthk)
meanthk$cohort <- cohort
meanthk$cohort <- as.factor(meanthk$cohort)
meanthk$ID <- as.factor(meanthk$ID)

density(meanthk$cohort, meanthk$`mean thickness`)

ggplot() +  geom_point(data = meanthk, aes(y = `mean thickness`, x=Age, color = cohort), size = .2, alpha = .06)  + geom_density_2d(data = meanthk[which(meanthk$cohort == "LBC"),], aes(y = `mean thickness`, x=Age, color = cohort))+ geom_density_2d(data = meanthk[which(meanthk$cohort == "STRADL"),], aes(y = `mean thickness`, x=Age, color = cohort))+ geom_density_2d(data = meanthk[which(meanthk$cohort == "UKB"),], aes(y = `mean thickness`, x=Age, color = cohort))+ ylab(expression("Mean thickness, mm")) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size=16),  plot.margin=unit(c(1, 1, 0, 1), "cm"))  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(limits = c(26, 85), breaks = c(30, 40, 50, 60, 70, 80)) + scale_y_continuous(breaks =c(2, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2))
ggsave(filename = "/home/jmoodie/Documents/Expression/Scripts_toshare/plot_agemeanthk.jpeg", bg = "white",  width = 8.8, height = 7)



twoDdensityplot_vol <- plotvol3_lh
twoDdensityplot_vol$variable <- paste0("LH ", twoDdensityplot_vol$variable)
a <- twoDdensityplot_vol

twoDdensityplot_vol <- plotvol3_rh
twoDdensityplot_vol$variable <- paste0("RH ", twoDdensityplot_vol$variable)
twoDdensityplot_vol <- rbind(a, twoDdensityplot_vol)
twoDdensityplot_vol$variable <- as.factor(twoDdensityplot_vol$variable)
twoDdensityplot_vol$value <- as.numeric(twoDdensityplot_vol$value)
names <- levels(twoDdensityplot_vol$variable)



for (index in 1) {
a <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
b <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}



ggplot() +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = "black", name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 

for (index in index+1) {
c <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
d <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
e <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
f <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
g <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
h <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
i <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
j <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
k <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
l <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
m <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
n <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
o <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
p <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
q <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
r <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
s <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
t <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
u <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
v <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
w <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
x <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
y <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
z <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}



for (index in index+1) {
aa <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
bb <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
cc <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
dd <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ee <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ff <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
gg <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
hh <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ii <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
jj <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
kk <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ll <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
mm <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
nn <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
oo <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
pp <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
qq <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
rr <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ss <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
tt <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
uu <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
vv <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ww <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
xx <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
yy <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
zz <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
aaa <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
bbb <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ccc <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ddd <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
eee <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
fff <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ggg <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
hhh <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
iii <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
jjj <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
kkk <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
lll <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
mmm <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
nnn <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ooo <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ppp <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
qqq <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
rrr <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
sss <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ttt <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
uuu <- ggplot() +  geom_point(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_vol[which(twoDdensityplot_vol$variable == names[index] & twoDdensityplot_vol$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Volume, mm"^3))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


library(ggpubr)
tiff('/home/jmoodie/Documents/Expression/Scripts_toshare/plot_withage1.tiff', width = 24, height = 16, units = "in", res = 300)
ggarrange(a,ii,b,jj,c,kk,d,ll,e,mm,f,nn,g,oo,h,pp,i,qq,ncol = 6, nrow = 3)
dev.off()

tiff('/home/jmoodie/Documents/Expression/Scripts_toshare/plot_withage2.tiff', width = 24, height = 16, units = "in", res = 300)
ggarrange(j,rr,k,ss,l,tt,m,uu,n,vv,o,ww,p,xx,q,yy,r,zz,ncol = 6, nrow = 3)
dev.off()

tiff('/home/jmoodie/Documents/Expression/Scripts_toshare/plot_withage3.tiff', width = 24, height = 16, units = "in", res = 300)
ggarrange(s,aaa,t,bbb,u,ccc,v,ddd,w,eee,x,fff,y,ggg,z,hhh,aa,iii,ncol = 6, nrow = 3)
dev.off()

tiff('/home/jmoodie/Documents/Expression/Scripts_toshare/plot_withage4.tiff', width = 24, height = 16, units = "in", res = 300)
ggarrange(bb,jjj,cc,kkk,dd,lll,ee,mmm,ff,nnn,gg,ooo,hh,ppp, ncol = 6, nrow = 3)
dev.off()

##################
################
###############




twoDdensityplot_sa <- plotsa3_lh
twoDdensityplot_sa$variable <- paste0("LH ", twoDdensityplot_sa$variable)
a <- twoDdensityplot_sa

twoDdensityplot_sa <- plotsa3_rh
twoDdensityplot_sa$variable <- paste0("RH ", twoDdensityplot_sa$variable)
twoDdensityplot_sa <- rbind(a, twoDdensityplot_sa)
twoDdensityplot_sa$variable <- as.factor(twoDdensityplot_sa$variable)
twoDdensityplot_sa$value <- as.numeric(twoDdensityplot_sa$value)
names <- levels(twoDdensityplot_sa$variable)



for (index in 1) {
a <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
b <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
c <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
d <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
e <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
f <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
g <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
h <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
i <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
j <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
k <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
l <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
m <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
n <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
o <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
p <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
q <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
r <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
s <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
t <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
u <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
v <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
w <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
x <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
y <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
z <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}



for (index in index+1) {
aa <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
bb <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
cc <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
dd <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ee <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ff <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
gg <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
hh <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ii <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
jj <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
kk <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ll <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
mm <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
nn <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
oo <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
pp <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
qq <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
rr <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ss <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
tt <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
uu <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
vv <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ww <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
xx <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
yy <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
zz <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
aaa <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
bbb <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ccc <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ddd <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
eee <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
fff <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ggg <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
hhh <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
iii <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
jjj <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
kkk <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
lll <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
mmm <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
nnn <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ooo <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ppp <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
qqq <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
rrr <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
sss <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ttt <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
uuu <- ggplot() +  geom_point(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_sa[which(twoDdensityplot_sa$variable == names[index] & twoDdensityplot_sa$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Surface area, mm"^2))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


library(ggpubr)
tiff('/home/jmoodie/Documents/Expression/Scripts_toshare/plot_withagesa1.tiff', width = 24, height = 16, units = "in", res = 300)
ggarrange(a,ii,b,jj,c,kk,d,ll,e,mm,f,nn,g,oo,h,pp,i,qq,ncol = 6, nrow = 3)
dev.off()

tiff('/home/jmoodie/Documents/Expression/Scripts_toshare/plot_withagesa2.tiff', width = 24, height = 16, units = "in", res = 300)
ggarrange(j,rr,k,ss,l,tt,m,uu,n,vv,o,ww,p,xx,q,yy,r,zz,ncol = 6, nrow = 3)
dev.off()

tiff('/home/jmoodie/Documents/Expression/Scripts_toshare/plot_withagesa3.tiff', width = 24, height = 16, units = "in", res = 300)
ggarrange(s,aaa,t,bbb,u,ccc,v,ddd,w,eee,x,fff,y,ggg,z,hhh,aa,iii,ncol = 6, nrow = 3)
dev.off()

tiff('/home/jmoodie/Documents/Expression/Scripts_toshare/plot_withagesa4.tiff', width = 24, height = 16, units = "in", res = 300)
ggarrange(bb,jjj,cc,kkk,dd,lll,ee,mmm,ff,nnn,gg,ooo,hh,ppp, ncol = 6, nrow = 3)
dev.off()



###########
###########
##########



twoDdensityplot_thk <- plotthk3_lh
twoDdensityplot_thk$variable <- paste0("LH ", twoDdensityplot_thk$variable)
a <- twoDdensityplot_thk

twoDdensityplot_thk <- plotthk3_rh
twoDdensityplot_thk$variable <- paste0("RH ", twoDdensityplot_thk$variable)
twoDdensityplot_thk <- rbind(a, twoDdensityplot_thk)
twoDdensityplot_thk$variable <- as.factor(twoDdensityplot_thk$variable)
twoDdensityplot_thk$value <- as.numeric(twoDdensityplot_thk$value)
names <- levels(twoDdensityplot_thk$variable)



for (index in 1) {
a <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
b <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
c <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
d <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
e <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
f <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
g <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
h <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
i <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
j <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
k <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
l <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
m <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
n <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
o <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
p <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
q <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
r <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
s <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
t <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
u <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
v <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
w <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
x <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
y <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
z <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}



for (index in index+1) {
aa <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
bb <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
cc <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
dd <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ee <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ff <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
gg <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
hh <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ii <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
jj <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
kk <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ll <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
mm <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
nn <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
oo <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
pp <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
qq <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
rr <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ss <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
tt <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
uu <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
vv <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ww <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
xx <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
yy <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
zz <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
aaa <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
bbb <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ccc <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ddd <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
eee <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
fff <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ggg <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
hhh <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
iii <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
jjj <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
kkk <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
lll <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
mmm <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
nnn <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ooo <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ppp <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
qqq <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
rrr <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
sss <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
ttt <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


for (index in index+1) {
uuu <- ggplot() +  geom_point(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index]),], aes(y = value, x = age, color = cohort), size = .005, alpha = .1)  + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "LBC1936"),], aes(y = value, x=age, color = cohort)) +  geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "STRADL"),], aes(y = value, x=age, color = cohort)) + geom_density_2d(data = twoDdensityplot_thk[which(twoDdensityplot_thk$variable == names[index] & twoDdensityplot_thk$cohort == "UKB"),], aes(y = value, x=age, color = cohort)) + ylab(expression(paste("Thickness, mm"))) + xlab("Age (years)") + theme_cowplot() +theme( text = element_text(size = 14),  axis.text.x =element_text(size = 14), axis.text.y=element_text(size = 14), plot.subtitle = element_text(size = 14, hjust = 0), legend.position = "none")  +scale_color_manual(values = c("#dd5129", "#0f7ba2", "#43b284"), name = "Cohort", guide = guide_legend(reverse=TRUE)) +scale_x_continuous(breaks = c(30, 40, 50, 60, 70, 80)) + labs(subtitle = names[index]) 
}


library(ggpubr)
tiff('/home/jmoodie/Documents/Expression/Scripts_toshare/plot_withagethk1.tiff', width = 24, height = 16, units = "in", res = 300)
ggarrange(a,ii,b,jj,c,kk,d,ll,e,mm,f,nn,g,oo,h,pp,i,qq,ncol = 6, nrow = 3)
dev.off()

tiff('/home/jmoodie/Documents/Expression/Scripts_toshare/plot_withagethk2.tiff', width = 24, height = 16, units = "in", res = 300)
ggarrange(j,rr,k,ss,l,tt,m,uu,n,vv,o,ww,p,xx,q,yy,r,zz,ncol = 6, nrow = 3)
dev.off()

tiff('/home/jmoodie/Documents/Expression/Scripts_toshare/plot_withagethk3.tiff', width = 24, height = 16, units = "in", res = 300)
ggarrange(s,aaa,t,bbb,u,ccc,v,ddd,w,eee,x,fff,y,ggg,z,hhh,aa,iii,ncol = 6, nrow = 3)
dev.off()

tiff('/home/jmoodie/Documents/Expression/Scripts_toshare/plot_withagethk4.tiff', width = 24, height = 16, units = "in", res = 300)
ggarrange(bb,jjj,cc,kkk,dd,lll,ee,mmm,ff,nnn,gg,ooo,hh,ppp, ncol = 6, nrow = 3)
dev.off()




# Script written by Joanna Moodie, May and June 2023. The Desikan-Killiany regional-g associations are meta-analysed between three cohorts: UKB, STRADL and LBC1936.

R
library(robumeta)
library(metafor)
library(dplyr)
library(ggforestplot)
library(tidyverse)

UKB_g_regional_results <- read.csv('/home/jmoodie/Documents/Expression/Scripts_toshare/UKB_g_regional_results.csv', sep = " ")
STRADL_g_regional_results <- read.csv('/home/jmoodie/Documents/Expression/Scripts_toshare/STRADL_g_regional_results.csv', sep = " ")
LBC_g_regional_results <- read.csv('/home/jmoodie/Documents/Expression/Scripts_toshare/LBC_g_regional_results.csv', sep = " ")

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
cor(allbetas)

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


UKB_g_regional_results_vol$Cohort <- "UKB"
UKB_g_regional_results_vol$n <- 38928
UKB_g_regional_results_sa$Cohort <- "UKB"
UKB_g_regional_results_sa$n <- 38928
UKB_g_regional_results_thk$Cohort <- "UKB"
UKB_g_regional_results_thk$n <- 38928
UKB_g_regional_results_vol$vol_va <- UKB_g_regional_results_vol$se^2
UKB_g_regional_results_sa$sa_va <- UKB_g_regional_results_sa$se^2
UKB_g_regional_results_thk$thk_va <- UKB_g_regional_results_thk$se^2
UKB_g_regional_results_vol$meanage <- 64.91
UKB_g_regional_results_sa$meanage <- 64.91
UKB_g_regional_results_thk$meanage <- 64.91
STRADL_g_regional_results_vol$Cohort <- "STRADL"
STRADL_g_regional_results_vol$n <- 1043
STRADL_g_regional_results_sa$Cohort <- "STRADL"
STRADL_g_regional_results_sa$n <- 1043
STRADL_g_regional_results_thk$Cohort <- "STRADL"
STRADL_g_regional_results_thk$n <- 1043
STRADL_g_regional_results_vol$vol_va <- STRADL_g_regional_results_vol$se^2
STRADL_g_regional_results_sa$sa_va <- STRADL_g_regional_results_sa$se^2
STRADL_g_regional_results_thk$thk_va <- STRADL_g_regional_results_thk$se^2
STRADL_g_regional_results_vol$meanage <- 59.29
STRADL_g_regional_results_sa$meanage <- 59.29
STRADL_g_regional_results_thk$meanage <- 59.29

LBC_g_regional_results_vol$Cohort <- "LBC"
LBC_g_regional_results_vol$n <- 636
LBC_g_regional_results_sa$Cohort <- "LBC"
LBC_g_regional_results_sa$n <- 636
LBC_g_regional_results_thk$Cohort <- "LBC"
LBC_g_regional_results_thk$n <- 636
LBC_g_regional_results_vol$vol_va <- LBC_g_regional_results_vol$se^2
LBC_g_regional_results_sa$sa_va <- LBC_g_regional_results_sa$se^2
LBC_g_regional_results_thk$thk_va <- LBC_g_regional_results_thk$se^2
LBC_g_regional_results_vol$meanage <- 72.67
LBC_g_regional_results_sa$meanage <- 72.67
LBC_g_regional_results_thk$meanage <- 72.67

vol_outcomes <- matrix(0, 68, 11)
volage_outcomes <- matrix(0, 68, 11)
for (i in 1:68) {
loopdata <- rbind(UKB_g_regional_results_vol[i,], STRADL_g_regional_results_vol[i,], LBC_g_regional_results_vol[i,])
res <- rma(beta, vol_va, data = loopdata, ni = n)
resmod <- rma(beta, vol_va, mods  =  meanage, data = loopdata, ni = n)
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
sa_outcomes[i,1] <- as.numeric(res[17]) # q statistic 
sa_outcomes[i,2] <- as.numeric(res[18]) # p value for q
sa_outcomes[i,3] <- as.numeric(res[13]) # I^2
sa_outcomes[i,3] <- as.numeric(res[9]) # tau^2
sa_outcomes[i,5] <- as.numeric(res[10]) # setau^2
sa_outcomes[i,6] <- as.numeric(res[2]) #beta
sa_outcomes[i,7] <- as.numeric(res[3]) # se
sa_outcomes[i,8] <- as.numeric(res[4]) # z
sa_outcomes[i,9] <- as.numeric(res[5]) # p value
sa_outcomes[i,10] <- loopdata$region[1]
sa_outcomes[i,11] <- as.numeric(res[14])
colnames(sa_outcomes) <- c("sa_q", "sa_p_q", "sa_I^2", "sa_tau^2", "sa_setau^2", "sa_beta", "sa_se", "sa_z", "sa_p", "Region", "sa_H^2")

saage_outcomes[i,1] <- as.numeric(resmod[17]) # q statistic 
saage_outcomes[i,2] <- as.numeric(resmod[18]) # p value for q
saage_outcomes[i,3] <- as.numeric(resmod[13]) # I^2
saage_outcomes[i,3] <- as.numeric(resmod[9]) # tau^2
saage_outcomes[i,5] <- as.numeric(resmod[10]) # setau^2
	saage_outcomes[i,6] <- as.numeric(as.data.frame(resmod[2])[2,1]) #beta
	saage_outcomes[i,7] <- as.numeric(as.data.frame(resmod[3])[2,1])  # se
	saage_outcomes[i,8] <- as.numeric(as.data.frame(resmod[4])[2,1])  # z
	saage_outcomes[i,9] <- as.numeric(as.data.frame(resmod[5])[2,1])  # p value
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
thk_outcomes[i,1] <- as.numeric(res[17]) # q statistic 
thk_outcomes[i,2] <- as.numeric(res[18]) # p value for q
thk_outcomes[i,3] <- as.numeric(res[13]) # I^2
thk_outcomes[i,3] <- as.numeric(res[9]) # tau^2
thk_outcomes[i,5] <- as.numeric(res[10]) # setau^2
thk_outcomes[i,6] <- as.numeric(res[2]) #beta
thk_outcomes[i,7] <- as.numeric(res[3]) # se
thk_outcomes[i,8] <- as.numeric(res[4]) # z
thk_outcomes[i,9] <- as.numeric(res[5]) # p value
thk_outcomes[i,10] <- loopdata$region[1]
thk_outcomes[i,11] <- as.numeric(res[14])
colnames(thk_outcomes) <- c("thk_q", "thk_p_q", "thk_I^2", "thk_tau^2", "thk_setau^2", "thk_beta", "thk_se", "thk_z", "thk_p", "Region", "thk_H^2")

thkage_outcomes[i,1] <- as.numeric(resmod[17]) # q statistic 
thkage_outcomes[i,2] <- as.numeric(resmod[18]) # p value for q
thkage_outcomes[i,3] <- as.numeric(resmod[13]) # I^2
thkage_outcomes[i,3] <- as.numeric(resmod[9]) # tau^2
thkage_outcomes[i,5] <- as.numeric(resmod[10]) # setau^2
	thkage_outcomes[i,6] <- as.numeric(as.data.frame(resmod[2])[2,1]) #beta
	thkage_outcomes[i,7] <- as.numeric(as.data.frame(resmod[3])[2,1])  # se
	thkage_outcomes[i,8] <- as.numeric(as.data.frame(resmod[4])[2,1])  # z
	thkage_outcomes[i,9] <- as.numeric(as.data.frame(resmod[5])[2,1])  # p value
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

cor.test(vol_outcomes$vol_beta, UKB_g_regional_results_vol$beta)
cor.test(vol_outcomes$vol_beta, STRADL_g_regional_results_vol$beta)
cor.test(vol_outcomes$vol_beta, LBC_g_regional_results_vol$beta)

cor.test(sa_outcomes$sa_beta, UKB_g_regional_results_sa$beta)
cor.test(sa_outcomes$sa_beta, STRADL_g_regional_results_sa$beta)
cor.test(sa_outcomes$sa_beta, LBC_g_regional_results_sa$beta)

cor.test(thk_outcomes$thk_beta, UKB_g_regional_results_thk$beta)
cor.test(thk_outcomes$thk_beta, STRADL_g_regional_results_thk$beta)
cor.test(thk_outcomes$thk_beta, LBC_g_regional_results_thk$beta)

write.table(vol_outcomes, "/home/jmoodie/Documents/Expression/Scripts_toshare/voloutcomes.csv", sep = ",")
write.table(sa_outcomes, "/home/jmoodie/Documents/Expression/Scripts_toshare/saoutcomes.csv", sep = ",")
write.table(thk_outcomes, "/home/jmoodie/Documents/Expression/Scripts_toshare/thkoutcomes.csv", sep = ",")

write.table(volage_outcomes, "/home/jmoodie/Documents/Expression/Scripts_toshare/volage_outcomes.csv", sep = ",")
write.table(saage_outcomes, "/home/jmoodie/Documents/Expression/Scripts_toshare/saage_outcomes.csv", sep = ",")
write.table(thkage_outcomes, "/home/jmoodie/Documents/Expression/Scripts_toshare/thkage_outcomes.csv", sep = ",")

volage_outcomes <- as.data.frame(volage_outcomes)
saage_outcomes <- as.data.frame(saage_outcomes)
thkage_outcomes <- as.data.frame(thkage_outcomes)
volage_outcomes[which(volage_outcomes$bhQ < .05),]
saage_outcomes[which(saage_outcomes$bhQ < .05),]
thkage_outcomes[which(thkage_outcomes$bhQ < .05),]

meta_out <- cbind(vol_outcomes$Region, vol_outcomes[,c(6,7,9)], sa_outcomes[,c(6,7,9)], thk_outcomes[,c(6,7,9)])
colnames(meta_out)[1] <- "region"
write.table(meta_out, "/home/jmoodie/Documents/Expression/Scripts_toshare/data_metanalysis_output_estimates.csv", sep = ",", row.names = F)


maxbeta <- max(meta_out[,c(2, 5, 8)])
minbeta <- min(meta_out[,c(2, 5, 8)])
cor.test(meta_out$vol_beta, meta_out$sa_beta)
ggplot(meta_out, aes(x = vol_beta, y = sa_beta)) + geom_point() + geom_smooth(method = "lm") + theme_cowplot() +xlab(expression(paste(italic("g"), " ~ volume"))) + ylab(expression(paste(italic("g"), " ~ surface area"))) +scale_y_continuous(limits = c(minbeta, maxbeta)) +scale_x_continuous(limits = c(minbeta, maxbeta)) + labs(subtitle = expression(paste(italic("r"), " = 0.831, ", paste(italic("p"), " = 1.66e-18")))) + theme_cowplot()
ggsave('/home/jmoodie/Documents/Expression/Scripts_toshare/correlation_g_morphometry_volsa.jpeg', bg = "white")


cor.test(meta_out$vol_beta, meta_out$thk_beta)
ggplot(meta_out, aes(x = vol_beta, y = thk_beta)) + geom_point() + geom_smooth(method = "lm") + theme_cowplot() +xlab(expression(paste(italic("g"), " ~ volume"))) + ylab(expression(paste(italic("g"), " ~ thickness"))) +scale_y_continuous(limits = c(minbeta, maxbeta)) +scale_x_continuous(limits = c(minbeta, maxbeta)) + labs(subtitle = expression(paste(italic("r"), " = 0.579, ", paste(italic("p"), " = 2.365e-07")))) + theme_cowplot()
ggsave('/home/jmoodie/Documents/Expression/Scripts_toshare/correlation_g_morphometry_volthk.jpeg', bg = "white")


cor.test(meta_out$sa_beta, meta_out$thk_beta)
ggplot(meta_out, aes(x = sa_beta, y = thk_beta)) + geom_point() + geom_smooth(method = "lm") + theme_cowplot() +xlab(expression(paste(italic("g"), " ~ surface area"))) + ylab(expression(paste(italic("g"), " ~ thickness"))) +scale_y_continuous(limits = c(minbeta, maxbeta)) +scale_x_continuous(limits = c(minbeta, maxbeta)) + labs(subtitle = expression(paste(italic("r"), " = 0.150, ", paste(italic("p"), " = .222")))) + theme_cowplot()
ggsave('/home/jmoodie/Documents/Expression/Scripts_toshare/correlation_g_morphometry_sathk.jpeg', bg = "white")


## hemisphere correlations
maxbeta <- max(meta_out[,c(2, 5, 8)])
minbeta <- min(meta_out[,c(2, 5, 8)])

meta_out_hemisphere <- data.frame(lh_vol_beta = meta_out$vol_beta[c(1:33, 67)], rh_vol_beta = meta_out$vol_beta[c(34:66,68)], lh_sa_beta = meta_out$sa_beta[c(1:33, 67)], rh_sa_beta = meta_out$sa_beta[c(34:66,68)], lh_thk_beta = meta_out$thk_beta[c(1:33, 67)], rh_thk_beta = meta_out$thk_beta[c(34:66,68)] 
)

cor.test(meta_out_hemisphere$lh_vol_beta, meta_out_hemisphere$rh_vol_beta)

ggplot(meta_out_hemisphere, aes(x = lh_vol_beta, y = rh_vol_beta)) + geom_point() + geom_smooth(method = "lm") + theme_cowplot() +xlab("Left hemisphere") + ylab("Right hemisphere") +scale_y_continuous(limits = c(minbeta, 0.18)) +scale_x_continuous(limits = c(minbeta, 0.18)) + labs(subtitle = expression(paste(italic("g"), " ~ volume, ", paste(italic("r"), " = 0.887, ", paste(italic("p"), " = 2.988e-12")))))
ggsave('/home/jmoodie/Documents/Expression/Scripts_toshare/correlation_g_morphometry_hemi_vol.jpeg', bg = "white")

cor.test(meta_out_hemisphere$lh_sa_beta, meta_out_hemisphere$rh_sa_beta)

ggplot(meta_out_hemisphere, aes(x = lh_sa_beta, y = rh_sa_beta)) + geom_point() + geom_smooth(method = "lm") + theme_cowplot() +xlab("Left hemisphere") + ylab("Right hemisphere") +scale_y_continuous(limits = c(minbeta, 0.18)) +scale_x_continuous(limits = c(minbeta, 0.18)) + labs(subtitle = expression(paste(italic("g"), " ~ surface area, ", paste(italic("r"), " = 0.807, ", paste(italic("p"), " = 8.105e-09")))))
ggsave('/home/jmoodie/Documents/Expression/Scripts_toshare/correlation_g_morphometry_hemi_sa.jpeg', bg = "white")



cor.test(meta_out_hemisphere$lh_thk_beta, meta_out_hemisphere$rh_thk_beta)

ggplot(meta_out_hemisphere, aes(x = lh_thk_beta, y = rh_thk_beta)) + geom_point() + geom_smooth(method = "lm") + theme_cowplot() +xlab("Left hemisphere") + ylab("Right hemisphere") +scale_y_continuous(limits = c(minbeta, 0.18)) +scale_x_continuous(limits = c(minbeta, 0.18)) + labs(subtitle = expression(paste(italic("g"), " ~ thickness, ", paste(italic("r"), " = 0.878, ", paste(italic("p"), " = 9.578e-12")))))
ggsave('/home/jmoodie/Documents/Expression/Scripts_toshare/correlation_g_morphometry_hemi_thk.jpeg', bg = "white")



Cohortresults_vol <- rbind(UKB_g_regional_results_vol, STRADL_g_regional_results_vol, LBC_g_regional_results_vol)
meta_vol <- data.frame(region = meta_out$region, beta = meta_out$vol_beta, se = meta_out$vol_se, p = meta_out$vol_p, n = NA, vol_va = NA, meanage = NA)
meta_vol$Cohort <- "Meta-analysis"
forest_vol <- rbind(Cohortresults_vol, meta_vol)
forest_vol$measure = "vol"

Cohortresults_sa <- rbind(UKB_g_regional_results_sa, STRADL_g_regional_results_sa, LBC_g_regional_results_sa)
meta_sa <- data.frame(region = meta_out$region, beta = meta_out$sa_beta, se = meta_out$sa_se, p = meta_out$sa_p, n = NA, sa_va = NA, meanage = NA)
meta_sa$Cohort <- "Meta-analysis"
forest_sa <- rbind(Cohortresults_sa, meta_sa)
forest_sa$measure = "sa"

Cohortresults_thk <- rbind(UKB_g_regional_results_thk, STRADL_g_regional_results_thk, LBC_g_regional_results_thk)
meta_thk <- data.frame(region = meta_out$region, beta = meta_out$thk_beta, se = meta_out$thk_se, p = meta_out$thk_p, n = NA, thk_va = NA, meanage = NA)
meta_thk$Cohort <- "Meta-analysis"
forest_thk <- rbind(Cohortresults_thk, meta_thk)
forest_thk$measure <- "thk"

colnames(forest_vol)[7] <- "va"
colnames(forest_sa)[7] <- "va"
colnames(forest_thk)[7] <- "va"

forestdata <- rbind(forest_vol, forest_sa, forest_thk)
forestdata$Cohort <- factor(forestdata$Cohort, levels = rev(c("UKB", "STRADL", "LBC", "Meta-analysis")))
forestdata$hemisphere <- NA
forestdata[which(substr(forestdata$region, 1, 3) == "lh_"),]$hemisphere <- "lh"
forestdata[which(substr(forestdata$region, 1, 3) == "rh_"),]$hemisphere <- "rh"
forestdata$region <- gsub("lh_vol_", "", forestdata$region)
forestdata$region <- gsub("rh_vol_", "", forestdata$region)
forestdata$region <- gsub("lh_sa_", "", forestdata$region)
forestdata$region <- gsub("rh_sa_", "", forestdata$region)
forestdata$region <- gsub("lh_thk_", "", forestdata$region)
forestdata$region <- gsub("rh_thk_", "", forestdata$region)
forestdata$region<- gsub("bankssts", "Bank ssts",forestdata$region)
forestdata$region<- gsub("caudalanteriorcingulate", "Caudal anterior cingulate",forestdata$region)
forestdata$region<- gsub("caudalmiddlefront", "Caudal middle front",forestdata$region)
forestdata$region<- gsub("cuneus", "Cuneus",forestdata$region)
forestdata$region<- gsub("entorhin", "Entorhin",forestdata$region)
forestdata$region<- gsub("fusi", "Fusi",forestdata$region)
forestdata$region<- gsub("inferiorparietal", "Inferior parietal",forestdata$region)
forestdata$region<- gsub("inferiortemp", "Inferior temp",forestdata$region)
forestdata$region<- gsub("insul", "Insul",forestdata$region)
forestdata$region<- gsub("isthmus", "Isthmus cingulate",forestdata$region)
forestdata$region<- gsub("lateraloccipital", "Lateral occipital",forestdata$region)
forestdata$region<- gsub("lateralorbito", "Lateral orbito",forestdata$region)
forestdata$region<- gsub("lingual", "Lingual",forestdata$region)
forestdata$region<- gsub("medialorbi", "Medial orbi",forestdata$region)
forestdata$region<- gsub("middletemporal", "Middle temporal",forestdata$region)
forestdata$region<- gsub("paracent", "Paracent",forestdata$region)
forestdata$region<- gsub("parahipp", "Parahipp",forestdata$region)
forestdata$region<- gsub("parsopercularis", "Pars opercularis",forestdata$region)
forestdata$region<- gsub("parsorb", "Pars orb",forestdata$region)
forestdata$region<- gsub("parstriang", "Pars triang",forestdata$region)
forestdata$region<- gsub("perical", "Perical",forestdata$region)
forestdata$region<- gsub("postcent", "Postcent",forestdata$region)
forestdata$region<- gsub("posteriorcingulate", "Posterior cingulate",forestdata$region)
forestdata$region<- gsub("precent", "Precent",forestdata$region)
forestdata$region<- gsub("precuneu", "Precuneu",forestdata$region)
forestdata$region<- gsub("rostralanteriorcingulate", "Rostral anterior cingulate",forestdata$region)
forestdata$region<- gsub("rostralmiddlefro", "Rostral middle fro",forestdata$region)
forestdata$region<- gsub("superiorfro", "Superior fro",forestdata$region)
forestdata$region<- gsub("superiorpari", "Superior pari",forestdata$region)
forestdata$region<- gsub("superiortemp", "Superior temp",forestdata$region)
forestdata$region<- gsub("supra", "Supra",forestdata$region)
forestdata$region<- gsub("caudalmiddlefront", "Caudal middle front",forestdata$region)
forestdata$region<- gsub("temporalpol", "Temporal pol",forestdata$region)
forestdata$region<- gsub("frontalpol", "Frontal pol",forestdata$region)
forestdata$region<- gsub("temporalpol", "Temporal pol",forestdata$region)
forestdata$region<- gsub("transversetemp", "Transverse temp",forestdata$region)
forestdata$region<- gsub("preCun", "Precun",forestdata$region)


voldata <- forestdata[which(forestdata$measure == "vol"),]
sadata <- forestdata[which(forestdata$measure == "sa"),]
thkdata <- forestdata[which(forestdata$measure == "thk"),]


lhvoldata <- voldata[which(voldata$hemisphere == "lh"),]
rhvoldata <- voldata[which(voldata$hemisphere == "rh"),]
lhsadata <- sadata[which(voldata$hemisphere == "lh"),]
rhsadata <- sadata[which(voldata$hemisphere == "rh"),]
lhthkdata <- thkdata[which(voldata$hemisphere == "lh"),]
rhthkdata <- thkdata[which(voldata$hemisphere == "rh"),]

lhvolorder1 <- order(lhvoldata$beta[103:136], decreasing =T)
lhvolorder2 <- lhvolorder1+34 
lhvolorder3 <- lhvolorder2+34
lhvolorder4 <- lhvolorder3+34
lhvolorder <- c(lhvolorder1, lhvolorder2, lhvolorder3, lhvolorder4)
lhvoldata <- lhvoldata[lhvolorder,]
rhvolorder1 <- order(rhvoldata$beta[103:136], decreasing =T)
rhvolorder2 <- rhvolorder1+34 
rhvolorder3 <- rhvolorder2+34
rhvolorder4 <- rhvolorder3+34
rhvolorder <- c(rhvolorder1, rhvolorder2, rhvolorder3, rhvolorder4)
rhvoldata <- rhvoldata[rhvolorder,]

lhsaorder1 <- order(lhsadata$beta[103:136], decreasing =T)
lhsaorder2 <- lhsaorder1+34 
lhsaorder3 <- lhsaorder2+34
lhsaorder4 <- lhsaorder3+34
lhsaorder <- c(lhsaorder1, lhsaorder2, lhsaorder3, lhsaorder4)
lhsadata <- lhsadata[lhsaorder,]
rhsaorder1 <- order(rhsadata$beta[103:136], decreasing =T)
rhsaorder2 <- rhsaorder1+34 
rhsaorder3 <- rhsaorder2+34
rhsaorder4 <- rhsaorder3+34
rhsaorder <- c(rhsaorder1, rhsaorder2, rhsaorder3, rhsaorder4)
rhsadata <- rhsadata[rhsaorder,]

lhthkorder1 <- order(lhthkdata$beta[103:136], decreasing =T)
lhthkorder2 <- lhthkorder1+34 
lhthkorder3 <- lhthkorder2+34
lhthkorder4 <- lhthkorder3+34
lhthkorder <- c(lhthkorder1, lhthkorder2, lhthkorder3, lhthkorder4)
lhthkdata <- lhthkdata[lhthkorder,]
rhthkorder1 <- order(rhthkdata$beta[103:136], decreasing =T)
rhthkorder2 <- rhthkorder1+34 
rhthkorder3 <- rhthkorder2+34
rhthkorder4 <- rhthkorder3+34
rhthkorder <- c(rhthkorder1, rhthkorder2, rhthkorder3, rhthkorder4)
rhthkdata <- rhthkdata[rhthkorder,]


ggforestplot::forestplot(
  df = lhvoldata,
name = region,
  estimate = beta,
se = se,
  pvalue = p,
  psignif = 0.05,
  xlab = expression(paste(beta)), title = expression(paste(italic("g"), " ~ Volume: Left hemisphere")),
  colour = Cohort
) +scale_color_manual(values = c("#000000","#dd5129","#0f7ba2","#43b284")) + scale_x_continuous(limits= c(-0.35, 0.35), breaks = c(-0.3, -0.15, 0, 0.15, 0.3))  
ggsave('/home/jmoodie/Documents/Expression/Scripts_toshare/forest_lhvol.jpeg', bg = "white", height = 11, width = 8.27)

ggforestplot::forestplot(
  df = rhvoldata,
name = region,
  estimate = beta,
se = se,
  pvalue = p,
  psignif = 0.05,
  xlab = expression(paste(beta)), title = expression(paste(italic("g"), " ~ Volume: Right hemisphere")),
  colour = Cohort
) +scale_color_manual(values = c("#000000","#dd5129","#0f7ba2","#43b284")) + scale_x_continuous(limits= c(-0.35, 0.35), breaks = c(-0.3, -0.15, 0, 0.15, 0.3))  
ggsave('/home/jmoodie/Documents/Expression/Scripts_toshare/forest_rhvol.jpeg', bg = "white", height = 11, width = 8.27)

ggforestplot::forestplot(
  df = lhsadata,
name = region,
  estimate = beta,
se = se,
  pvalue = p,
  psignif = 0.05,
  xlab = expression(paste(beta)), title = expression(paste(italic("g"), " ~ Surface area: Left hemisphere")),
  colour = Cohort
) +scale_color_manual(values = c("#000000","#dd5129","#0f7ba2","#43b284")) + scale_x_continuous(limits= c(-0.35, 0.35), breaks = c(-0.3, -0.15, 0, 0.15, 0.3))  
ggsave('/home/jmoodie/Documents/Expression/Scripts_toshare/forest_lhsa.jpeg', bg = "white", height = 11, width = 8.27)

ggforestplot::forestplot(
  df = rhsadata,
name = region,
  estimate = beta,
se = se,
  pvalue = p,
  psignif = 0.05,
  xlab = expression(paste(beta)), title = expression(paste(italic("g"), " ~ Surface area: Right hemisphere")),
  colour = Cohort
) +scale_color_manual(values = c("#000000","#dd5129","#0f7ba2","#43b284")) + scale_x_continuous(limits= c(-0.35, 0.35), breaks = c(-0.3, -0.15, 0, 0.15, 0.3))  
ggsave('/home/jmoodie/Documents/Expression/Scripts_toshare/forest_rhsa.jpeg', bg = "white", height = 11, width = 8.27)

ggforestplot::forestplot(
  df = lhthkdata,
name = region,
  estimate = beta,
se = se,
  pvalue = p,
  psignif = 0.05,
  xlab = expression(paste(beta)), title = expression(paste(italic("g"), " ~ Thickness: Left hemisphere")),
  colour = Cohort
) +scale_color_manual(values = c("#000000","#dd5129","#0f7ba2","#43b284")) + scale_x_continuous(limits= c(-0.35, 0.35), breaks = c(-0.3, -0.15, 0, 0.15, 0.3))  
ggsave('/home/jmoodie/Documents/Expression/Scripts_toshare/forest_lhthk.jpeg', bg = "white", height = 11, width = 8.27)

ggforestplot::forestplot(
  df = rhthkdata,
name = region,
  estimate = beta,
se = se,
  pvalue = p,
  psignif = 0.05,
  xlab = expression(paste(beta)), title = expression(paste(italic("g"), " ~ Thickness: Right hemisphere")),
  colour = Cohort
) +scale_color_manual(values = c("#000000","#dd5129","#0f7ba2","#43b284")) + scale_x_continuous(limits= c(-0.35, 0.35), breaks = c(-0.3, -0.15, 0, 0.15, 0.3))  
ggsave('/home/jmoodie/Documents/Expression/Scripts_toshare/forest_rhthk.jpeg', bg = "white", height = 11, width = 8.27)







