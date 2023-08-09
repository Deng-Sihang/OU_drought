rm(list=ls())
library(lme4)
library(car)
library(plyr)
library(vegan)
library(dplyr)
library(ieggr)
library(reshape2)
library(picante)
library(NST)
save.wd <- iwd(choose.dir())

### Calculation of MST values --------------------------------------------------------
com.file = "Bacteria-zotu-resampled.csv"
treat.file = "Treatment info-for MST.csv"

comm = read.csv(com.file, row.names = 1, header = T, sep = ",")
comm = comm[,1:48]
treat = read.csv(treat.file, row.names = 1, header = T, sep = ",")
comm = as.data.frame(t(comm))

dim(comm)
comm = comm[, colSums(comm) > 0]
dim(comm)

dim(comm)
comm = comm[rowSums(comm) > 0,]
dim(comm)

for(i in 1:6) ###6 represents the 6-year prolonged drought.
{
  comm1 <- comm[(8*i-7):(8*i),]
  
  samp = match.name(rn.list = list(treat = treat, comm1 = comm1))
  comm1 = samp$comm1
  treat1 = samp$treat
  
  set.seed(1028)
  ### taxonomic (Sorensen) diversity
  tMST <- tNST(comm = comm1, group = treat1, dist.method = "bray",
               null.model = 'PF', abundance.weighted = FALSE,
               rand = 1000, nworker = 80)
  if(i ==1)
  {
    MST_group <- tMST$index.pair.grp
    MST_result <- tMST$index.pair.grp
  }
  else{
    MST_group <- tMST$index.pair.grp
    MST_result <-rbind(MST_result, MST_group)
  }
  NST_result
}

### Temporal succession of MST ----------------------------------------------------------
### linear mixed-effects models
dat <- read.table("MSTvalue.csv",sep=",",header = TRUE)
treat.file = "Treatment info-Microbial biomass.csv"
treat = read.csv(treat.file, row.names = 1, header = T, sep = ",")
dat$treatment <- treat$Treatment[match(dat$name1,rownames(treat))]
dat$block1 <- treat$block[match(dat$name1,rownames(treat))]
dat$block2 <- treat$block[match(dat$name2,rownames(treat))]
dat$block <- paste(dat$block1, dat$block2)
dat$year <- treat$year[match(dat$name1,rownames(treat))]
dat$year <- dat$year-2009
dat[,5:7]<-scale(dat[,5:7])

tpNST.lmm <- function(tpNST,treat)
{ 
  drought=treat$Treatment
  drought.lev=unique(drought)
  out=list()
  for(j in 1:length(drought.lev))
  {
    tpNST1 = tpNST[(36*j-35):(36*j),]
    lmij=lmer(value~year+((1+year)|block),data=tpNST1)
    lmijsm=summary(lmij)
    AIC1=AIC(lmij)
    r2ij=rsquared.glmm(lmij)
    lmijCS=Anova(lmij,type = "II")
    out[[j]]=c(slop.fix=lmijsm$coefficients[2,1],slop.sd=coef(summary(lmij))[ , "Std. Error"][2],
               R2M=r2ij$Marginal,R2C=r2ij$Conditional,AIC1=AIC1,AIC2=r2ij$AIC,
               P.typeII=lmijCS[[3]],Chisq=lmijCS[[1]])
  }
  outs=Reduce(cbind,out)
  colnames(outs)=drought.lev
  outs
}

for (i in 5:7)
{
  tpNST <-data.frame(value=dat[,i],dat)
  tpNST.result =tpNST.lmm(tpNST= tpNST,treat = treat)
  if(i == 5)
  {
    output=tpNST.result
  }else{
    output=cbind(output,tpNST.result)
  }
  output
}

### Relationships between microbial diversity indices and relative abundances of C/N cycling genes --------------------------------------------
### function to get the R, slope, and p value from the linear mixed model
corenvs<-function(divtest,genes){
  sapply(colnames(genes),function(x){
    message(x)
    if(sum(abs(divtest-genes[,x]),na.rm = T)<0.01 | length(unique(genes[,x])) < 7 | length(unique(divtest)) < 7) { 
      result=list(r=NA,pvalue=NA)
    }else{
      div<-data.frame(divtest=divtest,genes=genes[,x],dat)
      div<-div[(!is.na(div$divtest)) & (!is.na(div$genes)),]
      if (length(unique(div$year))<2) {
        fm1<-lmer(scale(genes)~scale(divtest)+(1|block),data=div)
      } else {
        fm1<-lmer(scale(genes)~scale(divtest)+(1|year)+(1|block),data=div)
      }
      
      presult<-car::Anova(fm1,type=2)
      coefs<-coef(summary(fm1))[ , "Estimate"]
      pvalue=presult[,3]
      r2<-rsquared.glmm(fm1)
      r=r2[1,5]
      coefs <- coefs[[2]]
      result=list(r=r,pvalue=pvalue,slope=coefs)
    }
    result
  },simplify = F)
}

### Importing data
dat<-read.csv("Alpha-gene.csv",row.names = 1)
### Drought and control groups should be calculated separately
dat <- dat[dat$treatment == "Drought",]

divs<-dat[,1:9]  
genes<-dat[,10:144]   

divs1<-sapply(colnames(divs),function(y){
  message(y)
  divtest=divs[,y]
  corenvs(divtest=divtest,genes = genes)
},simplify = F)

### slope of the linear mixed model
test<-sapply(divs1, function(x){
  sapply(x, function(y){
    y$slope
  })
})

### p-value of the linear mixed model
test1<-sapply(divs1, function(x){
  sapply(x, function(y){
    y$pvalue
  })
})
