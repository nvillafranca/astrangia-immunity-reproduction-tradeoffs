#Setting working directory
setwd("~/Desktop/RI_21") 

#Importing data and making row 1 the row names
ThesisData <- read.csv("ThesisData.csv",row.names = 1)

#averaging sperm into new column
ThesisData$SpermAvg <- rowMeans(subset(ThesisData, select = c(SpermCount1, SpermCount2, SpermCount3)), na.rm = TRUE)

#averaging polyp counts
ThesisData$PolypAvg <- rowMeans(subset(ThesisData, select = c(PC1, PC2, PC3, na.rm = TRUE)))

#deleting sperm 1, 2, 3 count and symbiont 1 2 3  columns 
ThesisData <- ThesisData[-c(7:12)]

#deleting polyp count 1, 2, 3 
ThesisData <- ThesisData[-c(8:10)]

#standardizing sperm/size 
ThesisData$SpermSizeStand <- ThesisData$SpermAvg/ThesisData$Size.cm.

#standardizing sperm/polyp 
ThesisData$SpermPolyp <- ThesisData$SpermAvg/ThesisData$PolypAvg

#normalizing symbiont data 
ThesisData$nsymbiont_per_cm2 <- sqrt(ThesisData$symbiont_per_cm2)


#Getting rid of outlier from POX
RI_21_DATANoPOXOutliers <- subset(ThesisData, rownames(ThesisData) != c("RI55"))

#Getting rid of outlier from melanin
RI_21_DATANoMELOutliers <- subset(ThesisData, rownames(ThesisData) != c("RI51"))

#getting rid of outliers from lipid assay
ThesisDataLipidOutliers <- ThesisData[!(row.names(ThesisData) %in% c("RI58", "RI53")),]

#Getting rid of outliers from lipid assay + POX 
LipidPOXOutliers <- ThesisData[!(row.names(ThesisData) %in% c("RI55", "RI58", "RI53")),]

#subset with no melanin or lipid outliers
MelaninLipNoOutliers <- ThesisData[!(row.names(ThesisData) %in% c("RI51", "RI58", "RI53")),]


#doing t test on sperm polyp relationship 
Sperm <- ThesisData$SpermAvg
Polyp <- ThesisData$PolypAvg
t.test(Sperm, Polyp, paired = TRUE)
####SIGNIFICANT, p = 0.01022

library(dplyr)

#doing t test on symbiont white and brown (one way anova because type is a categorical variable)
one.way <- aov(nsymbiont_per_cm2 ~ Type, data = ThesisData)
summary(one.way) 
####SIGNIFICANT p = 2.43e-09 ***

#doing a linear regression on sperm count and symbiont 
linreg <- lm(formula = nsymbiont_per_cm2 ~ SpermPolyp, data = ThesisData)
print(summary(linreg))
####not significant 

###################################################################
#Collinearity tests to determine whether there is a difference in 
##lipid and carb content between browns and whites
lmcl <- (lm(mg.carb.mg.tissue + mg.lipid.mg.tissue, data = ThesisDataLipidOutliers))
summary(lmcl)

lmslc <- (lm(TypeN ~ mg.lipid.mg.tissue + mg.carb.mg.tissue, data = ThesisDataLipidOutliers))
summary(lmslc)

library("olsrr")
#testing for collinearity 
ols_vif_tol(lmslc)
#no collinearity

ols_eigen_cindex(lmslc)
#no collinearity


###General linear models to determine the relationshhip between symbiotic state and lipid/carb content
#GLM testing the relationship between symbiont count and carb
#testing for outleirs 
library(EnvStats)
(rosnerTest(ThesisData$mg.carb.mg.tissue))
#Plotting histograms and testing for normality 
#If data are not normal, transform and test for normality again
#Normal data = p = > 0.05
hist(ThesisData$mg.carb.mg.tissue, main = "carb", breaks = 10)
shapiro.test(ThesisData$mg.carb.mg.tissue)
#normal, no outliers 
glmcs <- (glm(ThesisData$mg.carb.mg.tissue ~ ThesisData$nsymbiont_per_cm2, family=gaussian))
summary(glmcs)


#GLM testing relationship between symbiont count and lipid 
#testing for outliers 
library(EnvStats)
rosnerTest(ThesisData$mg.lipid.mg.tissue, k = 3)
#2 outliers, RI 58 and 53, removed above 

#Plotting histograms and testing for normality 
#If data are not normal, transform and test for normality again
#Normal data = p = > 0.05
hist(ThesisDataLipidOutliers$mg.lipid.mg.tissue, main = "lipid", breaks = 10)
shapiro.test(ThesisDataLipidOutliers$mg.lipid.mg.tissue)
#normal

glmls <- glm(mg.lipid.mg.tissue ~ nsymbiont_per_cm2, data=ThesisDataLipidOutliers, family = gaussian)
summary(glmls)

###GLMs between assays and lipid and carb
catlip <- (glm(Cat.min.mg ~ mg.lipid.mg.tissue, data = ThesisDataLipidOutliers, family=gaussian))
summary(catlip)

catcar <- (glm(Cat.min.mg ~ mg.carb.mg.tissue, data = ThesisData, family=gaussian))
summary(catcar)

poxlip <- (glm(POX.min.mg ~ mg.lipid.mg.tissue, data = LipidPOXOutliers, family=gaussian))
summary(poxlip)

poxcarb <- (glm(POX.min.mg ~ mg.carb.mg.tissue, data = RI_21_DATANoPOXOutliers, family=gaussian))
summary(poxcarb)

ablip <- (glm(Doubling.Time ~ mg.lipid.mg.tissue, data = ThesisDataLipidOutliers, family = gaussian))
summary(ablip)

abcarb <- (glm(Doubling.Time ~ mg.carb.mg.tissue, data = ThesisData, family = gaussian))
summary(abcarb)

ppolip <- (glm(TransformedPPO ~ mg.lipid.mg.tissue, data = ThesisDataLipidOutliers, family=gaussian))
summary(ppolip)

ppocarb <- (glm(TransformedPPO ~ mg.carb.mg.tissue, data = ThesisData, family=gaussian))
summary(ppocarb)

ablip <- glm(Doubling.Time ~ mg.lipid.mg.tissue, data = ThesisDataLipidOutliers, family=gaussian)
summary(ablip)

abcarb <- glm(Doubling.Time ~ mg.carb.mg.tissue, data = ThesisData, family=gaussian)
summary(abcarb)

mellip <- (glm(mg.melanin.mg.tissue ~ mg.lipid.mg.tissue, data = MelaninLipNoOutliers, family=gaussian))
summary(mellip)
#significant

melcarb <- (lm(mg.melanin.mg.tissue ~ mg.carb.mg.tissue, data = RI_21_DATANoMELOutliers, family=gaussian))
summary(melcarb)
#significant 

#graphing melanin/lipid melanin/carb
library(ggplot2)
library(ggpubr)

melaninlipid<- ggplot(MelaninLipNoOutliers, aes(mg.melanin.mg.tissue, mg.lipid.mg.tissue)) + geom_smooth(method = "lm", color = "black", lty = 2) + geom_point(pch = 21, color = "black", size = 3, aes(fill = mg.melanin.mg.tissue)) +
  scale_fill_viridis_c(direction= -1, option = "A") + theme_bw() + theme(panel.grid = element_blank()) + xlab("Melanin Concentration\n(mg melanin/mg tissue)") + ylab("Lipid Concentration\n(mg lipid/mg tissue)")
melaninlipid    
ggsave("melaninlipid.jpg", melaninlipid, height = 6, width = 6)

melanincarb<- ggplot(RI_21_DATANoMELOutliers, aes(mg.melanin.mg.tissue, mg.carb.mg.tissue)) + geom_smooth(method = "lm", color = "black", lty = 2) + geom_point(pch = 21, color = "black", size = 3, aes(fill = mg.melanin.mg.tissue)) +
  scale_fill_viridis_c(direction= -1, option = "A") + theme_bw() + theme(panel.grid = element_blank()) + xlab("Melanin Concentration\n(mg melanin/mg tissue)") + ylab("Carbohydrate Concentration\n(mg carb/mg tissue)")
melanincarb    
ggsave("melanincarb.jpg", melanincarb, height = 6, width = 6)

#graphing comparison of melanin vs lipid/carb in browns and whites
#to get R and y values for each line, use facet_grid(~Type). Pearson correlation used to get the R value 
library(ggplot2)
lipmelbw <- ggplot(MelaninLipNoOutliers, aes(mg.melanin.mg.tissue, mg.lipid.mg.tissue)) + geom_smooth(aes(color = Type), method = "lm", se = FALSE) + geom_point(aes(fill = Type), pch = 21, color = "black", size = 3) + 
  scale_fill_manual(values = c('Brown' = "burlywood4", 'White' = "white")) + scale_color_manual(values = c('Brown' = "burlywood4", 'White' = "grey")) + theme_bw() + theme(panel.grid = element_blank()) +
  xlab("Melanin Concentration\n(mg melanin/mg tissue)") + ylab("Lipid Concentration\n(mg lipid/mg tissue)")
lipmelbw
ggsave("lipmelbw.pdf", lipmelbw, height = 6, width = 6)

#determining if symbiont type affects melanin/lipid relationship
mellipsymb <- glm(mg.melanin.mg.tissue ~ mg.lipid.mg.tissue + TypeN + mg.lipid.mg.tissue*TypeN, data = MelaninLipNoOutliers, family=gaussian)
summary(mellipsymb)
#not significant

carbmelbw <- ggplot(RI_21_DATANoMELOutliers, aes(mg.melanin.mg.tissue, mg.carb.mg.tissue)) + geom_smooth(aes(color = Type), method = "lm", se = FALSE) + geom_point(aes(fill = Type), pch = 21, color = "black", size = 3) + 
  scale_fill_manual(values = c('Brown' = "burlywood4", 'White' = "white")) + scale_color_manual(values = c('Brown' = "burlywood4", 'White' = "grey")) + theme_bw() + theme(panel.grid = element_blank()) + 
  xlab("Melanin Concentration\n(mg melanin/mg tissue)") + ylab("Carbohydrate Concentration\n(mg carb/mg tissue)")
carbmelbw
ggsave("carbmelbw.pdf", carbmelbw, height = 6, width = 6)

#determining if symbiont type affects melanin/carb relationship
melcarbsymb <- glm(mg.melanin.mg.tissue ~ mg.carb.mg.tissue + TypeN + mg.carb.mg.tissue*TypeN, data = RI_21_DATANoMELOutliers, family=gaussian)
summary(melcarbsymb)
#not significant

Griddy <- ggarrange(melaninlipid, melanincarb, lipmelbw, carbmelbw, nrow = 2, ncol = 2, align = "hv", legend = "right")
Griddy
ggsave("Griddy.pdf", Griddy, height = 10, width = 12)

#correlation test
cor(ThesisData$mg.lipid.mg.tissue, ThesisData$SymbioAvg)
#negative correlation

cor(ThesisData$mg.carb.mg.tissue, ThesisData$SymbioAvg)
#very weak correlation

########CATALASE#########
#CATALASE testing for outliers
library(EnvStats)
rosnerTest(ThesisData$Cat.min.mg, k = 3)
#no outliers

#CATALASE Plotting histograms and testing for normality 
#If data are not normal, transform and test for normality again
#Normal data = p = > 0.05
hist(ThesisData$Cat.min.mg, main = "catalase", breaks = 10)
shapiro.test(ThesisData$Cat.min.mg)
boxplot(ThesisData$Cat.min.mg)

#Running glm for CATALASE = sperm/polyp + symbiont, summary shows the results of the GLM 
CATALASE.GLM <- glm(ThesisData$Cat.min.mg ~ ThesisData$SpermPolyp + ThesisData$nsymbiont_per_cm2, family=gaussian())
summary(CATALASE.GLM)

##########POX########### was transformed and had one outlier
#testing for outliers 
rosnerTest(ThesisData$POX.min.mg)

#POX Plotting histograms and testing for normality 
#If data are not normal, transform and test for normality again
#Normal data = p = > 0.05
hist(RI_21_DATANoPOXOutliers$POX.min.mg, main = "pox", breaks = 10)
shapiro.test(RI_21_DATANoPOXOutliers$POX.min.mg)
boxplot(ThesisData$POX.min.mg)
#normal

#testing for outliers
rosnerTest(ThesisData$POX.min.mg, k = 1)
#RI55 was the outlier, removed outlier at the beginning

#Running glm for POX = sperm/polyp + symbiont, summary shows the results of the GLM 
POX.GLM <- glm(RI_21_DATANoPOXOutliers$POX.min.mg ~ RI_21_DATANoPOXOutliers$SpermPolyp + RI_21_DATANoPOXOutliers$nsymbiont_per_cm2, family=gaussian())
summary(POX.GLM)

##############PPO############ was transformed
#testing for outliers
rosnerTest(ThesisData$PPO.min.mg, k = 1)
#no outliers found

#PPO Plotting histograms and testing for normality of transformed PPO 
#If data are not normal, transform and test for normality again
#Normal data = p = > 0.05
hist(ThesisData$PPO.min.mg, main = "ppo", breaks = 10)
shapiro.test(ThesisData$PPO.min.mg)
boxplot(ThesisData$PPO.min.mg)
#not normal

#Transforming PPO data 
ThesisData$TransformedPPO <- log(ThesisData$PPO.min.mg)

#Confirming normality
shapiro.test(ThesisData$TransformedPPO)
#normal

#Running glm for ppo = sperm/polyp + symbiont, summary shows the results of the GLM 
PPO.GLM <- glm(ThesisData$TransformedPPO ~ ThesisData$SpermPolyp + ThesisData$nsymbiont_per_cm2, family=gaussian())
summary(PPO.GLM)


############ANTIBACTERIAL###########
#testing for outliers
rosnerTest(ThesisData$Doubling.Time, k = 1)

#ANTIBACTERIAL Plotting histograms and testing for normality of transformed ab 
#If data are not normal, transform and test for normality again
#Normal data = p = > 0.05
shapiro.test(ThesisData$Doubling.Time)

#Running glm for ab = sperm/polyp + symbiont, summary shows the results of the GLM 
ab.GLM <- glm(ThesisData$Doubling.Time ~ ThesisData$SpermPolyp + ThesisData$nsymbiont_per_cm2, family=gaussian())
summary(ab.GLM)


############MELANIN############ one outlier 
#testing for outliers
rosnerTest(ThesisData$mg.melanin.mg.tissue, k = 3)
#1 outlier, RI51

#Melanin Plotting histograms and testing for normality of transformed melanin 
#If data are not normal, transform and test for normality again
#Normal data = p = > 0.05
shapiro.test(RI_21_DATANoMELOutliers$mg.melanin.mg.tissue)
#normal, more than 0.05

#Running glm for melanin = sperm/polyp + symbiont, summary shows the results of the GLM 
melanin.GLM <- glm(RI_21_DATANoMELOutliers$mg.melanin.mg.tissue ~ RI_21_DATANoMELOutliers$SpermPolyp + RI_21_DATANoMELOutliers$nsymbiont_per_cm2, family=gaussian())
summary(melanin.GLM)


