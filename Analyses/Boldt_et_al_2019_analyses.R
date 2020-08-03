# This analysis code was written by Annika Boldt, 2015-2019
# Boldt, A., Blundell, C., & De Martino, B. (2019). Confidence
# modulates exploration and exploitation in value-based learning. Neuroscience
# of Consciousness, 2019(1), 1–12. https://doi.org/10.1093/nc/niz004
# All analyses are listed in the same order as in the paper.


rm(list=ls(all=TRUE)) #clears the workspace
load('Data/Boldt_et_al_2019_data')

library(plyr)
library(viridis)
library(Hmisc)
library(lme4)
library(fields)
library(ez)
library(reshape2)

head(data1)
head(data2)


# EXPERIMENT 1
# Participants formed value beliefs over time

errsplit = array(NA,c(npp1,2))
cj1split = array(NA,c(npp1,2))

for (isub in 1:npp1) {
  tmean=data1$banditsmean1[data1$sub==filenames1[isub] & data1$metatype==1 & data1$part==4]
  tmean2=data1$banditsmean2[data1$sub==filenames1[isub] & data1$metatype==1 & data1$part==4]    
  tobs=data1$whichobs[data1$sub==filenames1[isub] & data1$metatype==1 & data1$part==4]
  tmean[tobs==2]=tmean2[tobs==2]
  terr=abs(tmean*100 -
             data1$prob[data1$sub==filenames1[isub] & data1$metatype==1 & data1$part==4])
  tcj1=data1$cj1[data1$sub==filenames1[isub] & data1$metatype==1 & data1$part==4]    
  ttrial=data1$withinblocktrial[data1$sub==filenames1[isub] & data1$metatype==1 & data1$part==4]
  cutoffs=median(ttrial)
  errsplit[isub,1] = mean(terr[ttrial<=cutoffs], na.rm=T)
  errsplit[isub,2] = mean(terr[ttrial>cutoffs], na.rm=T)  
  cj1split[isub,1] = mean(tcj1[ttrial<=cutoffs], na.rm=T)
  cj1split[isub,2] = mean(tcj1[ttrial>cutoffs], na.rm=T) 
}

show(colMeans(errsplit))
show(t.test(errsplit[,1],errsplit[,2],paired=T))

show(colMeans(cj1split))
show(t.test(cj1split[,1],cj1split[,2],paired=T))




corrupdate = array(NA,c(npp1,2))

for(isub in 1:npp1) {
  datanow=data1[data1$sub==filenames1[isub] & data1$part==4 & data1$metatype==1,]
  twin=datanow$outcome[datanow$whichobs==1]
  tprob=datanow$prevprob1[datanow$whichobs==1]
  tprevprob=c(NA,tprob[1:(length(tprob)-1)])
  tcj1=datanow$prevcj11[datanow$whichobs==1]
  tprevcj1=c(NA,tcj1[1:(length(tcj1)-1)])
  twin=twin[2:length(twin)]
  tprob=tprob[2:length(tprob)]
  tprevprob=tprevprob[2:length(tprevprob)]
  tcj1=tcj1[2:length(tcj1)]
  tprevcj1=tprevcj1[2:length(tprevcj1)]
  tcj1diff=tcj1-tprevcj1
  tvalupdate=abs(twin-tprevprob)
  corrupdate[isub,1]=rcorr(tvalupdate,tcj1diff)$r[1,2]
  corrupdate[isub,2]=rcorr(tvalupdate,tcj1diff)$P[1,2]
}

length(which(corrupdate[,1]<0))
min(corrupdate[corrupdate[,1]<0,1])
max(corrupdate[corrupdate[,1]<0,1])
show(length(which(corrupdate[,2]<0.05)))
show(max(corrupdate[which(corrupdate[,2]<0.05),2]))






# FIGURE 2

aggtraces=ddply(data1[data1$part==4,], .(blocktype,withinblocktrial), summarise, meanzcj11=mean(zcj11, na.rm=T), meanzcj12=mean(zcj12, na.rm=T), meanzprob1=mean(zprob1, na.rm=T), meanzprob2=mean(zprob2, na.rm=T))

layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=rep(1,2), heights=rep(1,1))
par(cex.main = 2.4, mar = c(6, 6, 1, 1), mgp = c(3, 1, 0), cex.lab = 1.6,
    font.lab = 1.6, cex.axis = 1.4, bty = "n", lwd=4, pch=19, las=1)

iblock=9
tcj1 = aggtraces$meanzcj11[aggtraces$blocktype==iblock]
tcj1 = tcj1[!is.na(tcj1)]
if (length(tcj1)>0) {
  tprob = aggtraces$meanzprob1[aggtraces$blocktype==iblock]
  tprob = tprob[!is.na(tprob)]
  N1=length(tcj1) 
  tchoice = aggtraces$withinblocktrial[aggtraces$blocktype==iblock & is.na(aggtraces$meanzcj11) & is.na(aggtraces$meanzcj12)]
  ttemp = which(!is.na(aggtraces$meanzcj11[aggtraces$blocktype==iblock]))
  plot(tcj1, tprob, type="n", xlim=c(-1.7,1.7), ylim=c(-2.3,2.4),xlab="",ylab="", main="", axes=F)
  axis(2, c(-2,0,2), c(-2,0,2), cex.axis=2.0, lwd=4)
  mtext(text="Value Belief", side=2, las=3, line=3.3, cex=2.0, at=0)   
  axis(1, c(-1.5,0,1.5), c(-1.5,"0.0",1.5), cex.axis=2.0, lwd=4, mgp=c(3,1.5,0))   
  mtext(text="Belief Confidence", side=1, las=1, line=4.3, cex=2.0, at=0)
  lines(tcj1,tprob)
  
  tactprob=data1$zbanditsoutcomemean1new[data1$blocktype==iblock & data1$whichobs==1]
  for (itrial in 1:N1) {
    arrows(tcj1[itrial], tprob[itrial], tcj1[itrial], tactprob[itrial], length = 0.03, angle = 90, code = 2, lwd=1.5, col="darkgrey")
  }   
  if (unique(data1$banditsdiff[data1$blocktype==iblock]) < 0) {
    points(tcj1[c(1:4,7:length(tcj1))], tprob[c(1:4,7:length(tprob))], col=viridis(N1*2)[N1:1][c(1:4,7:length(tprob))], pch=19, cex=1.4)
    points(tcj1[5:6], tprob[5:6], col=viridis(N1*2)[N1:1][5:6], pch=8, cex=1.4)
  } else {
    points(tcj1, tprob, col=viridis(N1*2)[(N1*2):(N1*2-N1)], pch=19, cex=1.4)
  } 
  
  tcj1 = aggtraces$meanzcj12[aggtraces$blocktype==iblock]
  tcj1 = tcj1[!is.na(tcj1)]
  tprob = aggtraces$meanzprob2[aggtraces$blocktype==iblock]
  tprob = tprob[!is.na(tprob)]
  N2=length(tcj1)
  ttemp = which(!is.na(aggtraces$meanzcj12[aggtraces$blocktype==iblock]))             
  lines(tcj1,tprob)         
  tactprob=data1$zbanditsoutcomemean2new[data1$blocktype==iblock & data1$whichobs==2]
  for (itrial in 1:N2) {
    arrows(tcj1[itrial], tprob[itrial], tcj1[itrial], tactprob[itrial], length = 0.03, angle = 90, code = 2, lwd=1.5, col="darkgrey")
  }
  if (unique(data1$banditsdiff[data1$blocktype==iblock]) > 0) {                
    points(tcj1, tprob, col=viridis(N2*2)[N2:1], pch=19, cex=1.4)
  } else {
    points(tcj1, tprob, col=viridis(N2*2)[(N2*2):(N2*2-N2)], pch=19, cex=1.4)
  } 
  xcoord = 0
  for (i in (N2*2):(N2*2-N2+1)) {
    for (j in 1:10) {
      xcoord = xcoord + 0.01			
      lines(rep(xcoord,2)-0.5,c(1.8,2.0),col=viridis(N2*2)[i],lwd=1.3)
    }
  }	
  xcoord = 0
  for (i in N2:1) {
    for (j in 1:10) {
      xcoord = xcoord + 0.01			
      lines(rep(xcoord,2)-0.5,c(1.3,1.5),col=viridis(N2*2)[i],lwd=1.3)
    }
  }	
  text(-0.3,2.4,"Arm of the Bandit", cex=1.8, font=2)
  text(-0.8,1.9,"better", cex=1.8)
  text(-0.8,1.4,"worse", cex=1.8)
  text(0.1-0.5,1.1,"early", cex=1.0)
  text(0.93-0.5,1.1,"late", cex=1.0)
  
  mtext("A", side=3, las=1, line=-1.2, cex=3.0, at=-2.5, font=2)	
}	

zcj1quant = array(NA,c(npp1,5))
trialquant = array(NA,c(npp1,5))

for (isub in 1:npp1) {
  
  tzcj1=data1$zcj1[data1$sub==filenames1[isub] & data1$metatype==1 & data1$part==4]
  ttrial=data1$withinblocktrial[data1$sub==filenames1[isub] & data1$metatype==1 & data1$part==4]
  cutoffs=quantile(ttrial,probs = seq(0, 1, 0.2))
  
  for (iquant in 2:6) {
    zcj1quant[isub,iquant-1] = mean(tzcj1[ttrial>cutoffs[iquant-1] & ttrial<=cutoffs[iquant]])
    trialquant[isub,iquant-1] = mean(ttrial[ttrial>cutoffs[iquant-1] & ttrial<=cutoffs[iquant]])
  }
  
}
plot(1:5,zcj1quant[1,],type="n",xlim=c(0.5,5.5),ylim=c(-1.5,1.5),xlab="",ylab="", main="", axes=F)
axis(1, 1:5, 1:5, cex.axis=2.0, lwd=4, mgp=c(3,1.5,0))            
mtext(text="Trial Quantile", side=1, las=1, line=4.0, cex=2.0, at=3)
axis(2, c(-1.5,0,1.5), c(-1.5,"0.0",1.5), cex.axis=2.0, lwd=4)
mtext(text="Belief Confidence", side=2, las=3, line=3.8, cex=2.0, at=0)   
points(1:5,colMeans(zcj1quant), pch=1, col="black", cex=2.4, lwd=3)
for (i in 1:5) {
  arrows(i,mean(zcj1quant[,i])-sd(zcj1quant[,i]), 
         i, mean(zcj1quant[,i])+sd(zcj1quant[,i]), 
         length = 0.03, angle = 90, code = 3, lwd=3)
}
mtext("B", side=3, las=1, line=-1.2, cex=3.0, at=-0.9, font=2)	







# Linking belief confidence to decision confidence

datanow=subset(data1,filt==1 & part==4 & metatype==2 & !is.nan(prob_c) & !is.nan(prob_uc))
datanow$subnew=as.factor(datanow$sub)
M8b <- lmer(scale(cj2) ~ scale(cj1_c) + scale(cj1_uc) + 
              scale(prob_c) * scale(prob_uc) + scale(cornew) + scale(log(rt)) + 
              ((scale(cj1_c) + scale(cj1_uc) + 
                  scale(prob_c) * scale(prob_uc) + scale(cornew) + scale(log(rt)))|subnew), data = datanow, REML=FALSE)

M8bcoefs <- data.frame(coef(summary(M8b)))
library(lmerTest)
M8b.semTest <- lmer(scale(cj2) ~ scale(cj1_c) + scale(cj1_uc) + 
                      scale(prob_c) * scale(prob_uc) + scale(cornew) + scale(log(rt)) + 
                      ((scale(cj1_c) + scale(cj1_uc) + 
                          scale(prob_c) * scale(prob_uc) + scale(cornew) + scale(log(rt)))|subnew), data = datanow, REML=FALSE)
M8bcoefs$df.Satt <- coef(summary(M8b.semTest))[, 3]
M8bcoefs$p.Satt <- coef(summary(M8b.semTest))[, 5]
M8bcoefs
library(lme4)


unlist(summary(M8b))$AICtab.BIC

M8bcoefs[4,c(1,5)]
M8bcoefs[5,c(1,5)]
M8bcoefs[8,c(1,5)]
M8bcoefs[2,c(1,5)]
M8bcoefs[3,c(1,5)]
M8bcoefs[6,c(1,5)]
M8bcoefs[7,c(1,5)]


# FIGURE 3

# First we prepare the two colour grids at the bottom

# Fitting the model
zprob_c = scale(datanow$prob_c)
zprob_uc = scale(datanow$prob_uc)
zcj2 = scale(datanow$cj2)
smoothed_mod = mgcv::gam(zcj2~s(zprob_c, zprob_uc),data=datanow)
# Generaring the predictor grid
z_seq = seq(-2, 2, 0.1)
zprob_c = rep(z_seq, each=41)
zprob_uc = rep(z_seq, times=41)
sim_data = data.frame(zprob_c, zprob_uc)
# Using the smooth spline to simuluate outcome values for every point on the grid
smoothed_pred = predict(smoothed_mod, newdata=sim_data)
x1new = matrix(smoothed_pred,nrow=41,byrow=T)

# Fitting the model
zcj1_c = scale(datanow$cj1_c)
zcj1_uc = scale(datanow$cj1_uc)
zcj2 = scale(datanow$cj2)
smoothed_mod = mgcv::gam(zcj2~s(zcj1_c, zcj1_uc),data=datanow)
# Generaring the predictor grid
z_seq = seq(-2, 2, 0.1)
zcj1_c = rep(z_seq, each=41)
zcj1_uc = rep(z_seq, times=41)
sim_data = data.frame(zcj1_c, zcj1_uc)
# Using the smooth spline to simuluate outcome values for every point on the grid
smoothed_pred = predict(smoothed_mod, newdata=sim_data)
x2new = matrix(smoothed_pred,nrow=41,byrow=T)

layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE), widths=c(1,1), heights=c(1,1,1))
par(cex.main = 4.2, mar = c(1, 1, 4, 1), mgp = c(3, 1.3, 0), cex.lab = 1.6, font.lab = 1.6, 
    cex.axis = 1.4, bty = "n", lwd=6, pch=19, las=1)     

x=seq(-0.5,0.5,length=1000)
y=dnorm(x,mean=-0.1, sd=.135)
plot(x, y, type="n", xlab="", ylab="",xlim=c(-1.0,1.0),ylim=c(-1,6),main="",axes=FALSE) 
lines(x, y, col="darkgrey")
text(-0.2,3.34,"Value Belief A", cex=2.8, col="darkgrey")
y=dnorm(x,mean=0.07, sd=.10)
lines(x, y, col="azure4")
text(0.078,4.4,"Value Belief B", cex=2.8, col="azure4")
arrows(-0.5,0,0.53,0,col ="black",code=2,angle=30,length=0.2)
mtext(text="Value", side=1, las=1, line=-3.5, cex=1.8, at=0.43)
lines(c(-0.1-0.15*2,-0.1+0.15*2),c(6.2,6.2),col="darkgrey")
mtext("Belief Confidence A",side=3, las=1, line=0, at=-0.1, cex=1.8, col="darkgrey")
lines(c(0.07-0.135*2,0.07+0.135*2),c(5.4,5.4),col="azure4")
text(0.07,5.7,"Belief Confidence B",cex=2.8,col="azure4")
arrows(-0.1,-1,0.07,-1,col ="black",code=3,angle=90,length=0.1)
text(-0.015,-0.6,"DV", cex=2.8, col="black")
mtext("A", side=3, las=1, line=-0.5, at=-0.55, cex=3.6, font=2)


par(mar = c(4, 16, 3, 1))
plot(1,M8bcoefs$Estimate[1],type="n",xlim=c(0.5,7.5),ylim=c(-2.4,0.5),xlab="",ylab="",
     main="", axes=F)
axis(2, c(-0.5,0,0.5), c(-0.5,"0.0",0.5), cex.axis=2.0, lwd=4)
mtext(text="Decision\nConfidence", side=2, las=3, line=9.0, cex=2.2, at=0)  
mtext(text="low", side=2, las=3, line=4.3, cex=1.8, at=-0.45)  
mtext(text="high", side=2, las=3, line=4.3, cex=1.8, at=0.41)  
arrows(-0.55, -0.6, -0.55, 0.6, xpd = TRUE, code=3)
legend(90, 1.2, title="Participant", legend = filenames1, bty="n", col=rainbow(npp1), pch=19, cex=1.4)
lines(seq(0.5,7.5,0.5),rep(0,15),col="grey",lty="dashed")
colours=c(rep("darkgrey",3), rep("darkslategrey",3), rep("black",3))
neworder=c(4,5,8,2,3,6,7)
for (ivar in 1:7) {
  points(ivar,M8bcoefs$Estimate[neworder[ivar]], cex=3, col=colours[ivar])
  arrows(ivar,M8bcoefs$Estimate[neworder[ivar]]-M8bcoefs$Std..Error[neworder[ivar]], ivar,
         M8bcoefs$Estimate[neworder[ivar]]+M8bcoefs$Std..Error[neworder[ivar]],
         length = 0, angle = 30, code = 2, lwd=7, col=colours[ivar])
}
ttext=c("Value\n(chosen)","Value\n(unchosen)", "Interaction Value", 
        "Belief Confidence\n(chosen)", "Belief Confidence\n(unchosen)",
        "Accuracy", "Response Time")
sign=c("*","*","*","*","*","*","*")        
for (i in 1:7) {
  mtext(text=ttext[i], side=1, las=2, line=-17, cex=2.0, at=i)
  mtext(text=sign[i], side=1, las=2, line=-20, cex=4.0, at=i+0.16)
}
mtext("B", side=3, las=1, line=3.0, at=-1.1, cex=3.6, font=2)


par(mar = c(8, 22, 5, 5))
image.plot(x1new, col=viridis(300), axes=F, legend.mar=6, legend.width=1.5)
contour(x1new, add=T, vfont = c("sans serif", "plain"), lwd=0.5, lty="dashed")
axis(2, c(0,0.25,0.5,0.75,1), c(-2,-1,0,1,2), cex.axis=2.0, lwd=4)
mtext(text="Value Belief\n(unchosen; quintile)", side=2, las=3, line=3.8, cex=1.8, at=0.5)  
axis(1, c(0,0.25,0.5,0.75,1), c(-2,-1,0,1,2), cex.axis=2.0, lwd=4)
mtext(text="Value Belief\n(chosen; quintile)", side=1, las=1, line=6.3, cex=1.8, at=0.5)  
mtext("C", side=3, las=1, line=0, at=-0.31, cex=3.6, font=2)

image.plot(x2new, col=viridis(300), axes=F, legend.mar=6, legend.width=1.5)
contour(x2new, add=T, vfont = c("sans serif", "plain"), lwd=0.5, lty="dashed")
axis(2, c(0,0.25,0.5,0.75,1), c(-2,-1,0,1,2), cex.axis=2.0, lwd=4)
mtext(text="Belief Confidence\n(unchosen; quintile)", side=2, las=3, line=3.8, cex=1.8, at=0.5)  
axis(1, c(0,0.25,0.5,0.75,1), c(-2,-1,0,1,2), cex.axis=2.0, lwd=4)
mtext(text="Belief Confidence\n(chosen; quintile)", side=1, las=1, line=6.3, cex=1.8, at=0.5) 
mtext("D", side=3, las=1, line=0, at=-0.31, cex=3.6, font=2)





# EXPERIMENT 2

# Percent trials in which lower-value option was chosen
# Exp1:
datanow=data1[data1$part==4,]
datanow$pickedlowv=abs((datanow$prevprob1<datanow$prevprob2)+1-datanow$resp)
overview=ddply(datanow, .(sub), summarise, pickedlowv=mean(pickedlowv, na.rm=T))
mean(overview$pickedlowv)
# Exp2:
datanow=data2[data2$part==3,]
datanow$pickedlowv=abs((datanow$prevprob1<datanow$prevprob2)+1-datanow$resp)
overview=ddply(datanow, .(sub), summarise, pickedlowv=mean(pickedlowv, na.rm=T))
mean(overview$pickedlowv)



# Participants have meaningful insight into their gambling behavior

# Average choice RT
datanow=data2[data2$part==3,]
overview=ddply(datanow, .(sub), summarise, meanrt=mean(rt[rt>0 & datanow$filt==1], na.rm=T),
               meanerrnew=mean(errnew[filt==1], na.rm=T))
mean(overview$meanrt)
1-mean(overview$meanerrnew)





errnewquant03 = array(NA,c(npp2,5))
winquant03 = array(NA,c(npp2,5))

for (isub in 1:npp2) {
  
  terr=data2$err[data2$sub==filenames2[isub] & data2$metatype==2 & data2$filt==1]
  terrnew=data2$errnew[data2$sub==filenames2[isub] & data2$metatype==2 & data2$filt==1]
  tzcj2=data2$zcj2[data2$sub==filenames2[isub] & data2$metatype==2 & data2$filt==1]
  twin=data2$win[data2$sub==filenames2[isub] & data2$metatype==2 & data2$filt==1]    
  cutoffs=quantile(tzcj2,probs = seq(0, 1, 0.2))
  
  for (iquant in 2:6) {
    errnewquant03[isub,iquant-1] = mean(terrnew[tzcj2>cutoffs[iquant-1] & tzcj2<=cutoffs[iquant]], na.rm=T)        
    winquant03[isub,iquant-1] = mean(twin[tzcj2>cutoffs[iquant-1] & tzcj2<=cutoffs[iquant]], na.rm=T)
  }
  
}                                              




errquantlong = melt(cbind(filenames2,data.frame(errnewquant03)), id.vars=c("filenames2"))
errquantlong$filenames2 = as.factor(errquantlong$filenames2)
anovares=aov(value~(variable)+Error(filenames2/(variable)),errquantlong)
summary(anovares)
ezANOVA(data=errquantlong,within=.(variable),wid=.(filenames2),dv=.(value),type=2)

winquantlong = melt(cbind(filenames2,data.frame(winquant03)), id.vars=c("filenames2"))
winquantlong$filenames2 = as.factor(winquantlong$filenames2)
anovares=aov(value~(variable)+Error(filenames2/(variable)),winquantlong)
summary(anovares)
ezANOVA(data=winquantlong,within=.(variable),wid=.(filenames2),dv=.(value),type=2)



# please note: the analyses reported in the paper have been calculated using SPSS to allow easy calculation of the GG correction and linear trends
# to reproduce these analyses please either correct by hand using the GG epsilon
# or export data (wide format: errnewquant03) into SPSS and run this script:

# DATASET ACTIVATE DataSet0.
# GLM ERq1 ERq2 ERq3 ERq4 ERq5
# /WSFACTOR=quantile 5 Polynomial 
# /METHOD=SSTYPE(3)
# /PRINT=ETASQ 
# /CRITERIA=ALPHA(.05)
# /WSDESIGN=quantile.
# 
# GLM WINq1 WINq2 WINq3 WINq4 WINq5
# /WSFACTOR=quantile 5 Polynomial 
# /METHOD=SSTYPE(3)
# /PRINT=ETASQ 
# /CRITERIA=ALPHA(.05)
# /WSDESIGN=quantile.




# FIGURE 4

layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=rep(1,2), heights=rep(1,1))
par(cex.main = 2.4, mar = c(7, 6, 1, 1), mgp = c(3, 1, 0), cex.lab = 1.6,
    font.lab = 1.6, cex.axis = 1.4, bty = "n", lwd=4, pch=19, las=1)

plot(1:5,errnewquant03[1,],type="n",xlim=c(0.5,5.5),ylim=c(-0.03,0.7),xlab="",ylab="",
     main="", axes=F)
axis(1, 1:5, 1:5, cex.axis=2.0, lwd=4, mgp=c(3,1.4,0))
mtext(text="Decision Confidence\nQuantile", side=1, las=1, line=5.5, cex=2.0, at=3)
axis(2, c(0.0,0.2,0.4,0.6), c("0.0",0.2,0.4,0.6), cex.axis=2.0, lwd=4)
mtext(text="Error Rate", side=2, las=3, line=4.2, cex=2.0, at=0.3)
points(1:5,colMeans(errnewquant03), pch=1, col="black", cex=2.4, lwd=3)
for (i in 1:5) {
  arrows(i,mean(errnewquant03[,i])-std.error(errnewquant03[,i]), 
         i, mean(errnewquant03[,i])+std.error(errnewquant03[,i]), 
         length = 0.03, angle = 90, code = 3, lwd=3)
}
mtext("A", side=3, las=1, line=-1.5, cex=3.0, at=-0.65, font=2)

plot(1:5,winquant03[1,],type="n",xlim=c(0.5,5.5),ylim=c(40,73),xlab="",ylab="",
     main="", axes=F)
mtext(text="Rewards", side=2, las=3, line=3.8, cex=2.0, at=55)              
axis(1, 1:5, 1:5, cex.axis=2.0, lwd=4, mgp=c(3,1.4,0))
mtext(text="Decision Confidence\nQuantile", side=1, las=1, line=5.5, cex=2.0, at=3)
axis(2, c(40,50,60,70), c(40,50,60,70), cex.axis=2.0, lwd=4)
points(1:5,colMeans(winquant03), pch=1, col="black", cex=2.4, lwd=3)
for (i in 1:5) {
  arrows(i,mean(winquant03[,i])-std.error(winquant03[,i]), 
         i, mean(winquant03[,i])+std.error(winquant03[,i]), 
         length = 0.03, angle = 90, code = 3, lwd=3)
}       
mtext("B", side=3, las=1, line=-1.5, cex=3.0, at=-0.65, font=2) 






# Confidence-guided exploration

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



cj1bins_hl=ddply(datanow, .(sub,bins_h,bins_l), summarise, meanpickedlowv=mean(pickedlowv))
agg_cj1bins_hl=ddply(cj1bins_hl, .(bins_h,bins_l), summarise, gmeanpickedlowv=mean(meanpickedlowv), 
                     sepickedlowv=std.error(meanpickedlowv)) 




datanow$subnew=as.factor(datanow$sub)

M13 <- glmer(pickedlowv ~ scale(cj1_h) * scale(dv) +
               ((scale(cj1_h) * scale(dv))|subnew), 
             family = binomial(link="logit"), data = datanow)
M13coefs <- data.frame(coef(summary(M13)))
unlist(summary(M13))$AICtab.BIC

M5 <- glmer(pickedlowv ~ scale(cj1_l) * scale(cj1_h) * scale(dv) +
              ((scale(cj1_l) * scale(cj1_h) * scale(dv))|subnew), 
            family = binomial(link="logit"), data = datanow)
unlist(summary(M5))$AICtab.BIC


M13coefs[2,c(1,4)]
M13coefs[3,c(1,4)]
M13coefs[4,c(1,4)]





# FIGURE 5

# Prepare image plot
# Fitting the model
zcj1_h = scale(datanow$cj1_h)
zdv = scale(datanow$dv)
pickedlowv = scale(datanow$pickedlowv)
smoothed_mod = mgcv::gam(pickedlowv~s(zcj1_h, zdv),data=datanow)
# Generaring the predictor grid
z_seq = seq(-2, 2, 0.1)
zcj1_h = rep(z_seq, each=41)
zdv = rep(z_seq, times=41)
sim_data = data.frame(zcj1_h, zdv)
# Using the smooth spline to simuluate outcome values for every point on the grid
smoothed_pred = predict(smoothed_mod, newdata=sim_data)
x1new = matrix(smoothed_pred,nrow=41,byrow=T)


layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), widths=c(1,1), heights=c(1,1))
par(cex.main = 4.2, mar = c(10, 10, 4, 1), mgp = c(3, 1.3, 0), cex.lab = 1.6, font.lab = 1.6, 
    cex.axis = 1.4, bty = "n", lwd=6, pch=19, las=1)     
plot(1:2,agg_cj1bins_hl$gmeanpickedlowv[agg_cj1bins_hl$bins_l==1],type="n",xlim=c(0.8,2.2),ylim=c(0.08,0.42),axes=F,
     xlab="",ylab="")
lines(1:2,agg_cj1bins_hl$gmeanpickedlowv[agg_cj1bins_hl$bins_l==1],col="black", type="o", pch=19, cex=1.6, lty="solid")
for (i in 1:2) {
  arrows(i, agg_cj1bins_hl$gmeanpickedlowv[agg_cj1bins_hl$bins_l==1][i]-
           agg_cj1bins_hl$sepickedlowv[agg_cj1bins_hl$bins_l==1][i], 
         i, agg_cj1bins_hl$gmeanpickedlowv[agg_cj1bins_hl$bins_l==1][i]+
           agg_cj1bins_hl$sepickedlowv[agg_cj1bins_hl$bins_l==1][i], 
         col ="black", code = 3, angle = 90, length = 0.1)
}
lines(1:2,agg_cj1bins_hl$gmeanpickedlowv[agg_cj1bins_hl$bins_l==2],col="black", type="o", pch=19, cex=1.6, lty="dotted")
for (i in 1:2) {
  arrows(i, agg_cj1bins_hl$gmeanpickedlowv[agg_cj1bins_hl$bins_l==2][i]-
           agg_cj1bins_hl$sepickedlowv[agg_cj1bins_hl$bins_l==2][i], 
         i, agg_cj1bins_hl$gmeanpickedlowv[agg_cj1bins_hl$bins_l==2][i]+
           agg_cj1bins_hl$sepickedlowv[agg_cj1bins_hl$bins_l==2][i], 
         col ="black", code = 3, angle = 90, length = 0.1)
}
axis(1,1:2,c("low","high"), lwd=5, cex.axis=2.6, mgp = c(3, 1.8, 0))
mtext("Higher-value\nBelief Confidence", side=1, line=6.5, at=1.5, cex=2)
axis(2,c(0.1,0.2,0.3,0.4),c(10,20,30,40), lwd=5, cex.axis=2.6)   
mtext("Percentage of Exploration Trials", side=2, line=5, at=0.25, cex=2, las=3)
legend(1.4,0.31, title="Lower-value\nBelief Confidence", legend=c("low","high"), 
       bty="n", lty=c("solid","dotted"), col=c("black","black"), cex=2.2)
mtext("A", side=3, las=1, line=0.8, at=0.45, cex=3.6, font=2)




par(mar = c(11, 7, 9, 5))
image.plot(x1new, col=viridis(300), axes=F, legend.mar=6, legend.width=1.5)
contour(x1new, add=T, vfont = c("sans serif", "plain"), lwd=0.5, lty="dashed")
axis(2, c(0,0.25,0.5,0.75,1), c(-2,-1,0,1,2), cex.axis=2.0, lwd=4)
mtext(text="Difference in Value (Quintile)", side=2, las=3, line=3.8, cex=1.8, at=0.5)  
axis(1, c(0,0.25,0.5,0.75,1), c(-2,-1,0,1,2), cex.axis=2.0, lwd=4)
mtext(text="Higher-value\nBelief Confidence (Quintile)", side=1, las=1, line=5.3, cex=1.8, at=0.5)  
mtext(text="low     % exploration     high", side=4, las=3, line=1.3, cex=1.4, at=0.5)  
mtext("B", side=3, las=1, line=5.8, at=-0.25, cex=3.6, font=2)



par(mar = c(1, 26, 6, 16))
plot(1,M13coefs$Estimate[1],type="n",xlim=c(0.5,3.5),ylim=c(-2.6,0.5),xlab="",ylab="",
     main="", axes=F)
axis(2, c(-0.5,0,0.5), c(-0.5,"0.0",0.5), cex.axis=2.0, lwd=4)
mtext(text="Exploration", side=2, las=3, line=8.5, cex=2.2, at=0)  
mtext(text="less", side=2, las=3, line=4.8, cex=1.8, at=-0.4)  
mtext(text="more", side=2, las=3, line=4.8, cex=1.8, at=0.35)  
arrows(-0.12, -0.45, -0.12, 0.45, xpd = TRUE, code=3)
legend(90, 1.2, title="Participant", legend = filenames2, bty="n", col=rainbow(npp2), pch=19, cex=1.4)
lines(seq(0.5,3.5,0.5),rep(0,length(seq(0.5,3.5,0.5))),col="grey",lty="dashed")
neworder=c(3,2,4)
sign=c("*","*","*")
for (ivar in 1:3) {
  points(ivar, M13coefs$Estimate[neworder[ivar]], cex=3, col="black")
  arrows(ivar, M13coefs$Estimate[neworder[ivar]]-M13coefs$Std..Error[neworder[ivar]], ivar,
         M13coefs$Estimate[neworder[ivar]]+M13coefs$Std..Error[neworder[ivar]],
         length = 0, angle = 30, code = 2, lwd=7, col="black")
}
ttext=c("Difference in Value\n(DV)", "Higher-value\nBelief Confidence", "Interaction Higher-\nvalue Belief\nConfidence and DV")
for (i in 1:3) {
  mtext(text=ttext[i], side=1, las=2, line=-17, cex=2.0, at=i)
  mtext(text=sign[i], side=1, las=2, line=-21, cex=4.0, at=i+0.11)
}
mtext("C", side=3, las=1, line=0.3, at=-0.8, cex=3.6, font=2)





# Inter-individual differences in participants’ capability to track uncertainty

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


M6 <- lmer(scale(cj1) ~ scale(log(withinblocktrial)) * scale(sig) * scale(mu) + scale(arm) + 
             ((scale(log(withinblocktrial)) * scale(sig) * scale(mu) + scale(arm))|subnew), data = datanow, REML=FALSE)
M6coefs <- data.frame(coef(summary(M6)))
require(lmerTest)
M6.semTest <- lmer(scale(cj1) ~ scale(log(withinblocktrial)) * scale(sig) * scale(mu) + scale(arm) + 
                     ((scale(log(withinblocktrial)) * scale(sig) * scale(mu) + scale(arm))|subnew), data = datanow, REML=FALSE)
M6coefs$df.Satt <- coef(summary(M6.semTest))[, 3]
M6coefs$p.Satt <- coef(summary(M6.semTest))[, 5]

unlist(summary(M6))$AICtab.BIC


M6coefs[3,c(1,7)]
M6coefs[4,c(1,7)]
M6coefs[2,c(1,7)]
M6coefs[5,c(1,7)]
M6coefs[8,c(1,7)]
M6coefs[9,c(1,7)]
M6coefs[6,c(1,7)]
M6coefs[7,c(1,7)]


# these da and meta-da values were calculated using the scripts by Maniscalco & Lau (2012), see http://www.columbia.edu/~bsm2105/type2sdt/
da = c(0.1887584,1.994162,2.1731,1.680478,3.135439,2.651225,0.9291517,1.784013,1.914013,1.317915,
       2.87685,2.237916,1.630057,3.458307,0.8236186,1.905082,1.3894,0.7365487,1.783772,1.283469,
       1.673106,1.851337,1.939298,2.005032,0.7473835,1.16974,1.507387,2.497787,2.402886,1.292816)
metada=c(0.2555760,1.8502687,3.5466544,3.6581343,2.9490425,3.4551174,4.1011440,3.3465060,0.8587106,3.5783166,
         3.1823545,1.4198827,2.9011860,2.2373448,2.7532595,2.2676003,2.7032654,1.0070821,1.0994642,3.3547256,
         2.4147144,2.6739457,3.1057162,3.2472819,0.7916017,2.6723194,1.9246497,3.0526514,3.5502326,1.1945134)
metadadf=data.frame(sub=1:npp2,da=matrix(da,npp2,1),metada=matrix(metada,npp2,1))
names(metadadf)=list("sub","da","metada")
metadadf$Mratio=metadadf$metada/metadadf$da
metadadf$logMratio=log(metadadf$Mratio,10)
agg_metadadf=ddply(metadadf, .(), summarise, meanda=mean(da), meanmetada=mean(metada),
                   meanMratio=mean(Mratio), meanlogMratio=mean(logMratio,na.rm=T))

rcorr(coef(M6)$subnew[,3],metadadf$logMratio)
rcorr(coef(M6)$subnew[,4],metadadf$logMratio)




# FIGURE 6

sign=c("","*","*","","","","*","*") 
neworder=c(3,4,2,5,8,6,7,9)

layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(2,1), heights=rep(1,1))
par(cex.main = 2.4, mar = c(5.5, 10, 5, 1), mgp = c(3, 1, 0), cex.lab = 1.6,
    font.lab = 1.6, cex.axis = 1.4, bty = "n", lwd=4, pch=19, las=1)

plot(1,M6coefs$Estimate[neworder][1],type="n",xlim=c(0.5,8.5),ylim=c(-2.2,0.8),xlab="",ylab="",
     main="", axes=F)           #xlim=c(0,100)
axis(2, c(-0.5,0,0.5), c(-0.5,"0.0",0.5), cex.axis=1.8, lwd=4)

mtext(text="Belief Confidence", side=2, las=3, line=7.5, cex=1.8, at=0)  
mtext(text="low", side=2, las=3, line=3.8, cex=1.8, at=-0.4)  
mtext(text="high", side=2, las=3, line=3.8, cex=1.8, at=0.4)  
arrows(-1.5, -0.6, -1.5, 0.6, xpd = TRUE, code=3)

#abline(h=0,col="grey",lty="dashed")
xes = seq(0,9,0.1)
lines(xes, rep(0,length(xes)),col="grey",lty="dashed")
colours=c("darkgrey","black","black","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")
for (ivar in 1:8) {
  points(ivar, M6coefs$Estimate[neworder][ivar], cex=2, col=colours[neworder[ivar]-1])
  arrows(ivar, M6coefs$Estimate[neworder][ivar]-M6coefs$Std..Error[neworder][ivar], ivar, M6coefs$Estimate[neworder][ivar]+M6coefs$Std..Error[neworder][ivar],
         length = 0, angle = 30, code = 2, lwd=5, col=colours[neworder[ivar]-1])
}
ttext=c("log(Trial)", expression(paste("Outcome ",sigma,sep="")), expression(paste("Outcome ",mu,sep="")), "Arm", expression(paste("Interaction\nOutcome ",sigma,", log(Trial)",sep="")), expression(paste("Interaction\nOutcome ",mu,", log(Trial)",sep="")), expression(paste("Interaction\nOutcome ",mu,", ",sigma,sep="")), expression(paste("Interaction\nOutcome ",mu,", ",sigma,", log(Trial)",sep="")))
for (i in 1:8) {
  mtext(text=ttext[neworder[i]-1], side=1, las=2, line=-8.5, cex=1.5, at=i-0.02)
  mtext(text=sign[neworder[i]-1], side=1, las=2, line=-10.5, cex=4.0, at=i+0.18)
}
mtext("A", side=3, las=1, line=2.2, at=-2.0, cex=3.6, font=2)

par(mar = c(10, 6, 4, 1))    
plot(coef(M6)$subnew[,3],metadadf$logMratio,type="n",xlim=c(-1.1,1.0),ylim=c(-0.4,0.9),xlab="",ylab="",main="", axes=F)
axis(2, c(-0.4,0,0.4,0.8), c(-0.4,"0.0",0.4,0.8), cex.axis=1.8, lwd=4)
mtext(text="Metacognitive Efficiency", side=2, las=3, line=4.3, cex=1.8, at=0.2)
axis(1, c(-1.0,0,1.0), c(-1.0,0,1.0), cex.axis=1.8, lwd=4, mgp=c(3,1.5,0))                
mtext(text=expression(paste("Beta Outcome ",sigma,sep="")), side=1, las=1, line=3.5, cex=1.8, at=0)
points(coef(M6)$subnew[,3],metadadf$logMratio, pch=19, cex=1.8, col="black")
mtext("B", side=3, las=1, line=1.2, at=-1.9, cex=3.6, font=2)
xes = seq(-1.05,1.05,0.1)
tempres=unlist(lm(metadadf$logMratio ~ coef(M6)$subnew[,3]))
lines(xes, as.numeric(tempres[1]) + xes * as.numeric(tempres[2]))

