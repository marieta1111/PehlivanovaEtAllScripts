##############################################################################################################
#  Analysis code used in "Diminished Cortical Thickness is Associated with Impulsive Choice in Adolescence"  #
##############################################################################################################

require(mgcv)
library(visreg)
library(ggplot2)
library(reshape2)
library(corrplot)
library(voxel)
library(ppcor)

# add R functions to run models in individual networks and extract relevant statistics
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/regDataGam.R") 
source("/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/Rfunctions/regDataGamMultOut.R")

# load data
load("/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/DataK_n427_112917.Rda")

### demographic + cognitive analyses for manuscript
# descriptive statistics for k
mean(data_k$k)
sd(data_k$k)

# associations with age and sex in GAM models
gmAgeSex1 <- gam(logk ~ s(ageAtScan) + sex.o, data=data_k)
summary(gmAgeSex1)
gmAgeSex2 <- gam(logk ~ s(ageAtScan) + sex.o + s(ageAtScan, by=sex.o), data=data_k)
summary(gmAgeSex2)

# association with cognitive variables 
pcor.test(data_k$logk,data_k$Overall_Accuracy,data_k[,c("ageAtScan","sex.o","ageSq")],method="pearson")
pcor.test(data_k$logk,data_k$F1_Exec_Comp_Res_Accuracy,data_k[,c("ageAtScan","sex.o","ageSq")],method="pearson")
pcor.test(data_k$logk,data_k$F2_Social_Cog_Accuracy,data_k[,c("ageAtScan","sex.o","ageSq")],method="pearson")
pcor.test(data_k$logk,data_k$F3_Memory_Accuracy,data_k[,c("ageAtScan","sex.o","ageSq")],method="pearson")

### main analyses of NMF networks
# create vector of model statements
nmf20Names <-names(data_k)[1467:1486] # NMF names
ms<-paste0(nmf20Names[-17],"~logk+sex.o+s(ageAtScan)")

# running basic models for association between NMF networks and logk
modOut<-regDataGam(ms,data_k,nmf20Names[-17],1) # w/o the noise component
sigComp<-rownames(modOut[which(modOut$fdr.p<.05),])

# sample code for correlations between logk and CT networks 
pcor.test(data_k$logk,data_k[,nmf20Names[1]],data_k[,c("ageAtScan","sex.o","ageSq")],method="pearson")

### visualizations of CT-logk effects (scatterplots in Figure 4)
gm.comp14 <- gam(Nmf20C14 ~ logk + s(ageAtScan) + sex.o, method="REML", data=data_k)      
gm.comp15 <- gam(Nmf20C15 ~ logk + s(ageAtScan) + sex.o, method="REML", data=data_k)

# network 14 with labels
visreg(gm.comp14, 'logk', overlay=T,ylab="CT Scores in Network 14",xlab="Discount Rate (log k)",cex.axis=1.5,cex.lab=2,cex.main=2,points=list(cex=1.01))
# network 15 with labels
visreg(gm.comp15, 'logk', overlay=T,ylab="CT Scores in Network 15",xlab="Discount Rate (log k)",cex.axis=1.5,cex.lab=1.5,cex.main=2,points=list(cex=1.01))

### age effects 
# interaction models in individual NMF networks
mods<-paste0(sigComp,"~ ageScanY*logk + logk + sex.o + ageScanY + ageSq")
pvals = NA
for (i in 1:length(mods)){
        foo<-lm(as.formula(mods[i]),data=data_k)
        tab<-summary(foo)
        pvals[i]<-tab$coefficients[6,4]
}
median(pvals)
min(pvals)
max(pvals)

### interaction plots (Figure 5)
# dataset with top and bottom logk quartiles
dataQ14<-data_k[data_k$logkQ.o %in% c('Q1','Q4'),]
# network 14
gm14 <- gam(Nmf20C14 ~ s(ageScanY) + logkQ.o, data=dataQ14, method="REML")
plot14<-plotGAM(gm14, "ageScanY", "logkQ.o", orderedAsFactor = T, rawOrFitted = "raw", plotCI = T)

plot14 + theme_bw() + ylab("CT Scores in Network 14") + xlab("Age") +
       ggtitle("") + scale_size(range=c(8,20))+
       theme(plot.title = element_text(hjust = 0.5,size=22,face="bold"), panel.background = element_blank(), panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),axis.title=element_text(size=24),legend.text = element_text(size = 16, hjust = 3, vjust = 3),
       legend.title=element_text(size=16),axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)) +
       scale_colour_discrete(name="DD")+scale_x_continuous(breaks=seq(10,24,2))

# network 15 
gm15 <- gam(Nmf20C15 ~ s(ageScanY) + logkQ.o, data=dataQ14, method="REML")
plot15<-plotGAM(gm15, "ageScanY", "logkQ.o", orderedAsFactor = T, rawOrFitted = "raw", plotCI = T)

plot15 + theme_bw() + ylab("CT Scores in Network 15") + xlab("Age") +
       ggtitle("") + scale_size(range=c(8,20))+
       theme(plot.title = element_text(hjust = 0.5,size=22,face="bold"), panel.background = element_blank(), panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),axis.title=element_text(size=24),legend.text = element_text(size = 16, hjust = 3, vjust = 3),
       legend.title=element_text(size=16),axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)) +
       scale_colour_discrete(name="DD")+scale_x_continuous(breaks=seq(10,24,2))


### sensitivity analyses
# adding Medu as a covariate 
data_medu <- data_k[!is.na(data_k$meduCnbGo1),]
pcor.test(data_medu$logk,data_medu$meduCnbGo1,data_medu[,c("ageAtScan","sex.o","ageSq")],method="pearson")

ms6 <-paste0(nmf20Names[-17],"~logk+meduCnbGo1+sex.o+s(ageAtScan)")
ms6_model<-regDataGamMultOut(ms6,data_k,nmf20Names[-17],1,2)

# adding TBV rating as covariate
ms7 <-paste0(nmf20Names[-17],"~logk+tbv+sex.o+s(ageAtScan)")
ms7_model<-regDataGamMultOut(ms7,data_k,nmf20Names[-17],1,2)

# adding T1 rating as covariate
pcor.test(data_k$logk,data_k$newRating,data_k[,c("ageAtScan","sex.o","ageSq")],method="spearman")
ms8 <-paste0(nmf20Names[-17],"~logk+newRating+sex.o+s(ageAtScan)")
ms8_model<-regDataGamMultOut(ms8,data_k,nmf20Names[-17],1,2)

# adding general cognitive factor as covariate 
ms9 <-paste0(nmf20Names[-17],"~logk+Overall_Accuracy+sex.o+s(ageAtScan)")
ms9_model<-regDataGamMultOut(ms9,data_k,nmf20Names[-17],1,2)

# sensitivity analysis on subjects not using medication
data_k_noMed <-data_k[which(data_k$medsComb==0),] # subsetting data with subjects NOT using medications
modOutNoMeds<-regDataGam(ms,data_k_noMed,nmf20Names[-17],1) # w/o the noise component
medSigComp<-rownames(modOutNoMeds[which(modOutNoMeds$fdr.p<.05),])


### analyses with regional JLF cortical thickness
# read in relevant JLF CT names for analysis
load("/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/jlfCtnames.Rda")
n<-length(jlfCtNames1)

# model statements
ms_k <-paste0(jlfCtNames1,"~logk+sex.o+s(ageAtScan)")
# running models
jlfCtMods <- regDataGam(ms_k,data_k,jlfCtNames1,1)

### multivariate prediction with NMF components 
# baseline covariates
covs1<-"sex.o+ageAtScan"
covs2<-"sex.o+ageAtScan+F1_Exec_Comp_Res_Accuracy+F3_Memory_Accuracy"
nmfcovs<-"Nmf20C1+Nmf20C2+Nmf20C3+Nmf20C4+Nmf20C5+Nmf20C6+Nmf20C7+Nmf20C8+Nmf20C9+Nmf20C10+Nmf20C11+Nmf20C12+Nmf20C13+Nmf20C14+Nmf20C15+Nmf20C16+Nmf20C18+Nmf20C19+Nmf20C20+"
covsMedu<-"sex.o+ageAtScan+meduCnbGo1"

# linear models
lmBase1<-lm(as.formula(paste0("logk~",covs1)),data=data_k)	          # model with demographic covariates
lmBase2<-lm(as.formula(paste0("logk~",covs2)),data=data_k)	          # model with demographics + cognitive variables
lmBase3<-lm(as.formula(paste0("logk~",covs1)),data=data_medu)             # model with demographic covariates in dataset with maternal education
lmBase4<-lm(as.formula(paste0("logk~",covsMedu)),data=data_medu)          # model with demographics + maternal in dataset with maternal education
lmCt1<-lm(as.formula(paste0("logk~",nmfcovs,covs1)),data=data_k)          # model with demographics + NMF networks
lmCt2<-lm(as.formula(paste0("logk~",nmfcovs,covs2)),data=data_k)          # model with demographics + cognitive variables + NMF networks
lmCtMedu<-lm(as.formula(paste0("logk~",nmfcovs,covsMedu)),data=data_medu) # model with demographics + maternal education + NMF networks

# model with just demographics
cor.test(predict(lmBase1), data_k$logk)
# model with just demographics + cognition
cor.test(predict(lmBase2),data_k$logk) 
# CT model against baseline of just age and sex
anova(lmCt1,lmBase1)
cor.test(predict(lmCt1),data_k$logk)
# CT model against baseline of just age and sex + cognition
anova(lmCt2,lmBase2)
cor.test(predict(lmCt2),data_k$logk)
# model with just demographics + medu
cor.test(predict(lmBase4),data_k$logk[-which(is.na(data_k$meduCnbGo1))])
# CT model against baseline of just age and sex + medu
anova(lmCtMedu,lmBase4)
cor.test(predict(lmCtMedu),data_medu$logk)

### plotting multivariate prediction (Figure 7)
# plotting values predicted from model with NMF networks + demographics against actual log k values
gamCovs1<-"sex.o+s(ageAtScan)"
nmfcovs<-"Nmf20C1+Nmf20C2+Nmf20C3+Nmf20C4+Nmf20C5+Nmf20C6+Nmf20C7+Nmf20C8+Nmf20C9+Nmf20C10+Nmf20C11+Nmf20C12+Nmf20C13+Nmf20C14+Nmf20C15+Nmf20C16+Nmf20C18+Nmf20C19+Nmf20C20+"
nmfModel<-gam(as.formula(paste0("logk~",nmfcovs,gamCovs1)),data=data_k,method="REML")
data_k$nmfPred<-predict(nmfModel)

plotNmfDemoModel <- lm(logk~nmfPred,data=data_k)
cor.test(data_k$logk,data_k$nmfPred,method="pearson")

# plotting
par("las"=1)
visreg(plotNmfDemoModel, 'nmfPred',xlab="Predicted Discount Rate (Log k)",ylab="Actual Discount Rate (Log k)",cex.axis=1.5,cex.lab=1.5,cex.main=1.8,points=list(cex=0.9),line.par=list(col="blue"))
legend(x='bottomright', legend='r = 0.33, p < 0.0001',bty = "n",cex=1.235)

