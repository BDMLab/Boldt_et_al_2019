# This analysis code was written by Annika Boldt, 2015-2019
# Boldt, A., Blundell, C., & De Martino, B. (2019). Confidence
# modulates exploration and exploitation in value-based learning. Neuroscience
# of Consciousness, 2019(1), 1–12. https://doi.org/10.1093/nc/niz004
# All analyses are listed in the same order as in the supplement.


rm(list=ls(all=TRUE)) #clears the workspace
load('Data/Boldt_et_al_2019_data')


# EXPERIMENT 1

# Model Comparisons: Predicting Decision Confidence

datanow=subset(data1,filt==1 & part==4 & metatype==2 & !is.nan(prob_c) & !is.nan(prob_uc))
datanow$subnew=as.factor(datanow$sub)

M9b <- lmer(scale(cj2) ~ scale(withinblocktrial) + scale(cj1_c) + scale(cj1_uc) + 
              scale(prob_c) * scale(prob_uc) + scale(cornew) + scale(log(rt)) + 
              ((scale(withinblocktrial) + scale(cj1_c) + scale(cj1_uc) + 
                  scale(prob_c) * scale(prob_uc) + scale(cornew) + scale(log(rt)))|subnew), data = datanow, REML=FALSE)
unlist(summary(M9b))$AICtab.BIC

M9 <- lmer(scale(cj2) ~ scale(withinblocktrial) + scale(cj1_c) * scale(cj1_uc) + 
             scale(prob_c) * scale(prob_uc) + scale(cornew) + scale(log(rt)) + 
             ((scale(withinblocktrial) + scale(cj1_c) * scale(cj1_uc) + 
                 scale(prob_c) * scale(prob_uc) + scale(cornew) + scale(log(rt)))|subnew), data = datanow, REML=FALSE)
unlist(summary(M9))$AICtab.BIC

M8b <- lmer(scale(cj2) ~ scale(cj1_c) + scale(cj1_uc) + 
              scale(prob_c) * scale(prob_uc) + scale(cornew) + scale(log(rt)) + 
              ((scale(cj1_c) + scale(cj1_uc) + 
                  scale(prob_c) * scale(prob_uc) + scale(cornew) + scale(log(rt)))|subnew), data = datanow, REML=FALSE)
unlist(summary(M8b))$AICtab.BIC

M7b <- lmer(scale(cj2) ~ scale(cj1_c) + scale(cj1_uc) + 
              scale(prob_c) * scale(prob_uc) + scale(cornew) + 
              ((scale(cj1_c) + scale(cj1_uc) + 
                  scale(prob_c) * scale(prob_uc) + scale(cornew))|subnew), data = datanow, REML=FALSE)
unlist(summary(M7b))$AICtab.BIC

M6 <- lmer(scale(cj2) ~ scale(cj1_c) * scale(cj1_uc) + 
             scale(prob_c) * scale(prob_uc) + 
             ((scale(cj1_c) * scale(cj1_uc) + 
                 scale(prob_c) * scale(prob_uc))|subnew), data = datanow, REML=FALSE)
unlist(summary(M6))$AICtab.BIC

M5 <- lmer(scale(cj2) ~ scale(cj1_c) + scale(cj1_uc) + 
             scale(prob_c) * scale(prob_uc) + 
             ((scale(cj1_c) + scale(cj1_uc) + 
                 scale(prob_c) * scale(prob_uc))|subnew), data = datanow, REML=FALSE)
unlist(summary(M5))$AICtab.BIC

M4 <- lmer(scale(cj2) ~ scale(cj1_c) + 
             scale(prob_c) * scale(prob_uc) + 
             ((scale(cj1_c) + 
                 scale(prob_c) * scale(prob_uc))|subnew), data = datanow, REML=FALSE)
unlist(summary(M4))$AICtab.BIC

M3 <- lmer(scale(cj2) ~ scale(prob_c) * scale(prob_uc) + 
             ((scale(prob_c) * scale(prob_uc))|subnew), data = datanow, REML=FALSE)
unlist(summary(M3))$AICtab.BIC

M2 <- lmer(scale(cj2) ~ scale(prob_c) + scale(prob_uc) + 
             ((scale(prob_c) + scale(prob_uc))|subnew), data = datanow, REML=FALSE)
unlist(summary(M2))$AICtab.BIC

M1 <- lmer(scale(cj2) ~ scale(prob_c) + 
             ((scale(prob_c))|subnew), data = datanow, REML=FALSE)
unlist(summary(M1))$AICtab.BIC



plotBICs = c(unlist(summary(M1))$AICtab.BIC,
             unlist(summary(M2))$AICtab.BIC,
             unlist(summary(M3))$AICtab.BIC,
             unlist(summary(M4))$AICtab.BIC,
             unlist(summary(M5))$AICtab.BIC,
             unlist(summary(M6))$AICtab.BIC,
             unlist(summary(M7b))$AICtab.BIC,
             unlist(summary(M8b))$AICtab.BIC,
             unlist(summary(M9b))$AICtab.BIC,
             unlist(summary(M9))$AICtab.BIC)
which(min(plotBICs)==plotBICs)
plotBICs


layout(1)
par(cex.main = 4.2, mar = c(10, 6, 1.5, 0), mgp = c(3, 1, 0), cex.lab = 1.6, font.lab = 1.6, cex.axis = 1.4, bty = "n", lwd=6, pch=19, las=1)
mp <- barplot(plotBICs-5000, col=c(rep("white",7),"grey","white","white"), offset=5000, beside = TRUE, names.arg=c("Model 1","Model 2","Model 3","Model 4","Model 5","Model 6","Model 7","Model 8","Model 9","Model 10"), ylim=c(5000,6500), main="", ylab="", axes=FALSE, las=2, cex.names=2.4)
axis(2, c(5000,5500,6000,6500), c(5000,5500,6000,"BIC"), cex.axis=2.4, lwd=6)






# EXPERIMENT 2

# Model Comparisons: Predicting Exploration

datanow=data2       
datanow=subset(datanow, filt==1 & !is.nan(prevcj11))
datanow$highv=(datanow$prevprob1<datanow$prevprob2)+1    
datanow$cj1_h=NA
datanow$cj1_l=NA
datanow$cj1_h[datanow$highv==1]=datanow$prevcj11[datanow$highv==1]
datanow$cj1_h[datanow$highv==2]=datanow$prevcj12[datanow$highv==2]
datanow$cj1_l[datanow$highv==1]=datanow$prevcj12[datanow$highv==1]
datanow$cj1_l[datanow$highv==2]=datanow$prevcj11[datanow$highv==2]    
datanow$pickedlowv=abs((datanow$prevprob1<datanow$prevprob2)+1-datanow$resp)

datanow$prob_h=NA
datanow$prob_l=NA
datanow$prob_h[datanow$highv==1]=datanow$prevprob1[datanow$highv==1]
datanow$prob_h[datanow$highv==2]=datanow$prevprob2[datanow$highv==2]
datanow$prob_l[datanow$highv==1]=datanow$prevprob2[datanow$highv==1]
datanow$prob_l[datanow$highv==2]=datanow$prevprob1[datanow$highv==2]  
datanow$dv=datanow$prob_h-datanow$prob_l

for(isub in 1:npp2) {
  cutoffs=median(datanow$cj1_h[datanow$sub==filenames2[isub]])
  datanow$bins_h[datanow$sub==filenames2[isub]]=2
  datanow$bins_h[datanow$cj1_h<=cutoffs & datanow$sub==filenames2[isub]]=1
  cutoffs=median(datanow$cj1_l[datanow$sub==filenames2[isub]])
  datanow$bins_l[datanow$sub==filenames2[isub]]=2
  datanow$bins_l[datanow$cj1_l<=cutoffs & datanow$sub==filenames2[isub]]=1    
}

datanow$subnew=as.factor(datanow$sub)

M7 <- glmer(pickedlowv ~ scale(cj1_h) +
              ((scale(cj1_h))|subnew), 
            family = binomial(link="logit"), data = datanow)
unlist(summary(M7))$AICtab.BIC

M12 <- glmer(pickedlowv ~ scale(cj1_h) + scale(dv) +
               ((scale(cj1_h) + scale(dv))|subnew), 
             family = binomial(link="logit"), data = datanow)
unlist(summary(M12))$AICtab.BIC

M13 <- glmer(pickedlowv ~ scale(cj1_h) * scale(dv) +
               ((scale(cj1_h) * scale(dv))|subnew), 
             family = binomial(link="logit"), data = datanow)
unlist(summary(M13))$AICtab.BIC

M14 <- glmer(pickedlowv ~ scale(cj1_h) * scale(dv) + scale(cj1_l) +
               ((scale(cj1_h) * scale(dv) + scale(cj1_l))|subnew), 
             family = binomial(link="logit"), data = datanow)
unlist(summary(M14))$AICtab.BIC

M5 <- glmer(pickedlowv ~ scale(cj1_l) * scale(cj1_h) * scale(dv) +
              ((scale(cj1_l) * scale(cj1_h) * scale(dv))|subnew), 
            family = binomial(link="logit"), data = datanow)
unlist(summary(M5))$AICtab.BIC
coefs <- data.frame(coef(summary(M5)))
M5coefs <- data.frame(coef(summary(M5)))




plotmyBICs = c(unlist(summary(M7))$AICtab.BIC,
               unlist(summary(M12))$AICtab.BIC,
               unlist(summary(M13))$AICtab.BIC,
               unlist(summary(M14))$AICtab.BIC,
               unlist(summary(M5))$AICtab.BIC)


layout(matrix(c(1,2), 2, 1, byrow = TRUE), widths=c(1), heights=c(1,1))
par(cex.main = 4.2, mar = c(10, 10, 3, 0), mgp = c(3, 1, 0), cex.lab = 1.6, font.lab = 1.6, cex.axis = 1.4, bty = "n", lwd=6, pch=19, las=1)
mp <- barplot(plotmyBICs-10500, col=c(rep("white",2),"grey","white","black"), offset=10500, beside = TRUE, names.arg=c("Model 1","Model 2","Model 3","Model 4","Model 5"), ylim=c(10500,11500), main="", ylab="", axes=FALSE, las=2, cex.names=2.4)
axis(2, c(10500,10700,10900,11100,11300,11500), c(10500,10700,10900,11100,11300,"BIC"), cex.axis=2.4, lwd=6)
mtext("A)", side=3, las=1, line=0.5, at=-1.05, cex=3.6)

par(mar = c(1, 12, 4, 4))
plot(1,coefs$Estimate[1],type="n",xlim=c(0.5,7.5),ylim=c(-4.6,0.5),xlab="",ylab="",
     main="", axes=F)
axis(2, c(-0.5,0,0.5), c(-0.5,"0.0",0.5), cex.axis=2.0, lwd=4)
mtext(text="Exploration", side=2, las=3, line=7.5, cex=2.2, at=0)  
mtext(text="less", side=2, las=3, line=4.0, cex=1.8, at=-0.45)  
mtext(text="more", side=2, las=3, line=4.0, cex=1.8, at=0.40)  
arrows(-0.85, -0.6, -0.85, 0.6, xpd = TRUE, code=3)
lines(seq(0.5,7.5,0.5),rep(0,length(seq(0.5,7.5,0.5))),col="grey",lty="dashed")
neworder=c(4,3,2,7,6,5,8)
sign=c("*","*","","*","","*","")
for (ivar in 1:47) {
  points(ivar, coefs$Estimate[neworder[ivar]], cex=3, col="black")
  arrows(ivar, coefs$Estimate[neworder[ivar]]-coefs$Std..Error[neworder[ivar]], ivar,
         coefs$Estimate[neworder[ivar]]+coefs$Std..Error[neworder[ivar]],
         length = 0, angle = 30, code = 2, lwd=7, col="black")
}
ttext=c("Difference in Value\n(DV)", "Higher-value\nBelief Confidence", "Lower-value\nBelief Confidence", "Interaction\nDV and Higher-value\nBelief Confidence","Interaction\nDV and Lower-value\nBelief Confidence","Interaction\nBeliefs Confidences", "Interaction\nDV and Belief Confidences")
for (i in 1:7) {
  mtext(text=ttext[i], side=1, las=2, line=-16, cex=1.6, at=i)
  mtext(text=sign[i], side=1, las=2, line=-18, cex=4.0, at=i+0.13)
}
mtext("B)", side=3, las=1, line=1.4, at=-1.55, cex=3.6)





# Model Comparisons: Predicting Belief Confidence

datanow1=subset(data2,part==3 & metatype==1)
datanow1$ttemp=1
datanow2=subset(data2,part==3 & metatype==1)
datanow2$ttemp=2
datanow=rbind(datanow1,datanow2)
datanow$subnew=as.factor(datanow$sub)
datanow$cj1=c(datanow1$cj11[datanow1$ttemp==1],datanow2$cj12[datanow2$ttemp==2])
datanow$sig=c(datanow1$banditsoutcomevar1[datanow1$ttemp==1],datanow2$banditsoutcomevar2[datanow2$ttemp==2])
datanow$mu=c(datanow1$banditsoutcomemean1[datanow1$ttemp==1],datanow2$banditsoutcomemean2[datanow2$ttemp==2])
datanow$arm=c(datanow1$ttemp[datanow1$ttemp==1],datanow2$ttemp[datanow2$ttemp==2])

theBICs = rep(NA,8)

require(lme4)
M1 <- lmer(scale(cj1) ~ scale(sig) + 
             ((scale(sig))|subnew), data = datanow, REML=FALSE)
theBICs[1]=unlist(summary(M1))$AICtab.BIC
M1coefs <- data.frame(coef(summary(M1)))
require(lmerTest)
M1.semTest <- lmer(scale(cj1) ~ scale(sig) + 
                     ((scale(sig))|subnew), data = datanow, REML=FALSE)
M1coefs$df.Satt <- coef(summary(M1.semTest))[, 3]
M1coefs$p.Satt <- coef(summary(M1.semTest))[, 5]
M1coefs

require(lme4)
M2 <- lmer(scale(cj1) ~ scale(sig) + scale(mu) + 
             ((scale(sig) + scale(mu))|subnew), data = datanow, REML=FALSE)
theBICs[2]=unlist(summary(M2))$AICtab.BIC
M2coefs <- data.frame(coef(summary(M2)))
require(lmerTest)
M2.semTest <- lmer(scale(cj1) ~ scale(sig) + scale(mu) + 
                     ((scale(sig) + scale(mu))|subnew), data = datanow, REML=FALSE)
M2coefs$df.Satt <- coef(summary(M2.semTest))[, 3]
M2coefs$p.Satt <- coef(summary(M2.semTest))[, 5]
M2coefs

require(lme4)
M3 <- lmer(scale(cj1) ~ scale(sig) * scale(mu) + 
             ((scale(sig) * scale(mu))|subnew), data = datanow, REML=FALSE)
theBICs[3]=unlist(summary(M3))$AICtab.BIC
M3coefs <- data.frame(coef(summary(M3)))
require(lmerTest)
M3.semTest <- lmer(scale(cj1) ~ scale(sig) * scale(mu) + 
                     ((scale(sig) * scale(mu))|subnew), data = datanow, REML=FALSE)
M3coefs$df.Satt <- coef(summary(M3.semTest))[, 3]
M3coefs$p.Satt <- coef(summary(M3.semTest))[, 5]
M3coefs

require(lme4)
M4 <- lmer(scale(cj1) ~ scale(log(withinblocktrial)) + scale(sig) * scale(mu) + 
             ((scale(log(withinblocktrial)) + scale(sig) * scale(mu))|subnew), data = datanow, REML=FALSE)
theBICs[4]=unlist(summary(M4))$AICtab.BIC
M4coefs <- data.frame(coef(summary(M4)))
require(lmerTest)
M4.semTest <- lmer(scale(cj1) ~ scale(log(withinblocktrial)) + scale(sig) * scale(mu) + 
                     ((scale(log(withinblocktrial)) + scale(sig) * scale(mu))|subnew), data = datanow, REML=FALSE)
M4coefs$df.Satt <- coef(summary(M4.semTest))[, 3]
M4coefs$p.Satt <- coef(summary(M4.semTest))[, 5]
M4coefs

require(lme4)
M5 <- lmer(scale(cj1) ~ scale(log(withinblocktrial)) * scale(sig) * scale(mu) + 
             ((scale(log(withinblocktrial)) * scale(sig) * scale(mu))|subnew), data = datanow, REML=FALSE)
theBICs[5]=unlist(summary(M5))$AICtab.BIC
M5coefs <- data.frame(coef(summary(M5)))
require(lmerTest)
M5.semTest <- lmer(scale(cj1) ~ scale(log(withinblocktrial)) * scale(sig) * scale(mu) + 
                     ((scale(log(withinblocktrial)) * scale(sig) * scale(mu))|subnew), data = datanow, REML=FALSE)
M5coefs$df.Satt <- coef(summary(M5.semTest))[, 3]
M5coefs$p.Satt <- coef(summary(M5.semTest))[, 5]
M5coefs

require(lme4)
M6 <- lmer(scale(cj1) ~ scale(log(withinblocktrial)) * scale(sig) * scale(mu) + scale(arm) + 
             ((scale(log(withinblocktrial)) * scale(sig) * scale(mu) + scale(arm))|subnew), data = datanow, REML=FALSE)
theBICs[6]=unlist(summary(M6))$AICtab.BIC
M6coefs <- data.frame(coef(summary(M6)))
require(lmerTest)
M6.semTest <- lmer(scale(cj1) ~ scale(log(withinblocktrial)) * scale(sig) * scale(mu) + scale(arm) + 
                     ((scale(log(withinblocktrial)) * scale(sig) * scale(mu) + scale(arm))|subnew), data = datanow, REML=FALSE)
M6coefs$df.Satt <- coef(summary(M6.semTest))[, 3]
M6coefs$p.Satt <- coef(summary(M6.semTest))[, 5]
M6coefs

require(lme4)
M7 <- lmer(scale(cj1) ~ scale(log(withinblocktrial)) * scale(sig) * scale(mu) * scale(arm) + 
             ((scale(log(withinblocktrial)) * scale(sig) * scale(mu) * scale(arm))|subnew), data = datanow, REML=FALSE)
theBICs[7]=unlist(summary(M7))$AICtab.BIC
M7coefs <- data.frame(coef(summary(M7)))
require(lmerTest)
M7coefs



layout(1)
par(cex.main = 4.2, mar = c(10, 7, 1.5, 0), mgp = c(3, 1, 0), cex.lab = 1.6, font.lab = 1.6, cex.axis = 1.4, bty = "n", lwd=6, pch=19, las=1)
mp <- barplot(theBICs[1:7]-15000, col=c(rep("white",5),"grey","white"), offset=15000, beside = TRUE, names.arg=c("Model 1","Model 2","Model 3","Model 4","Model 5","Model 6","Model 7"), ylim=c(15000,20000), main="", ylab="", axes=FALSE, las=2, cex.names=2.4)
axis(2, c(15000,16000,17000,18000,19000,20000), c(15000,16000,17000,18000,19000,"BIC"), cex.axis=2.4, lwd=6)

