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
### Taking bacterial communities as an example for calculation
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

### PLS --------------------------------------------
############
#install.packages("OmicsPLS")
library(OmicsPLS)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ropls")
library(ropls)
#################
#################
R2ff<-function(xm,ym,o2p)
{
  c(R2x=mean(sapply(1:ncol(xm),function(i){1-sum((o2p$X_hat[,i]-xm[,i])^2)/sum((xm[,i]-mean(xm[,i]))^2)})),
    R2y=mean(sapply(1:ncol(ym),function(i){1-sum((o2p$Y_hat[,i]-ym[,i])^2)/sum((ym[,i]-mean(ym[,i]))^2)})))
}
R2sep<-function(xm,pls){(((getVipVn(pls))^2)/ncol(xm))*getSummaryDF(pls)[,'R2Y(cum)']}

plstest<-function(xm,ym,rand=100)
{
  nx=ncol(xm)
  combs=list()
  for(i in 1:nx)
  {
    message("i=",i," ",date())
    cbni=combn(nx,i)
    combs=c(combs,lapply(1:ncol(cbni),function(i){cbni[,i]}))
  }
  message("Total of ",length(combs)," combinations. ",date())
  #trac=seq(from=1,to=length(combs),by=50)
  id=t(rep(0,nx))
  colnames(id)=colnames(xm)
  
  dfna=t(rep(NA,8))
  colnames(dfna)=c('R2X(cum)','R2Y(cum)','Q2(cum)','RMSEE','pre','ort','pR2Y','pQ2')
  
  res=lapply(1:length(combs),
             function(i)
             {
               message("----- PLS i=",i," ",date())
               xmi=xm[,combs[[i]],drop=FALSE]
               pls=try(opls(x=xmi,y=ym,predI=NA,orthoI=0,permI=rand,fig.pdfC='none',info.txtC='none'))
               if(class(pls)=="try-error")
               {
                 pls=try(opls(x=xmi,y=ym,predI=1,orthoI=0,permI=rand,fig.pdfC='none',info.txtC='none'))
               }
               out=dfna
               if(class(pls)!="try-error"){dfi=getSummaryDF(pls);out[1,match(colnames(dfi),colnames(dfna))]=as.vector(as.matrix(dfi[1,]))}
               idni=id
               idni[combs[[i]]]=1
               cbind(idni,out)
             })
  res
}

plsfw<-function(xm,ym,r2buf=0.98,Q2ck=FALSE,SEEck=FALSE,rand=100)
{
  nx=ncol(xm)
  fs<-list(integer(0))
  dfna=t(rep(NA,8))
  colnames(dfna)=c('R2X(cum)','R2Y(cum)','Q2(cum)','RMSEE','pre','ort','pR2Y','pQ2')
  id=t(rep(0,nx))
  colnames(id)=colnames(xm)
  R2Yn=list()
  fijn=list()
  outrc=list()
  if(Q2ck){Q2n=list()}
  if(SEEck){SEEn=list()}
  k=1
  kn=1
  for(fn in 1:nx)
  {
    for(i in 1:length(fs))
    {
      fsi=fs[[i]]
      fai=which(!((1:nx) %in% fsi))
      for(j in 1:length(fai))
      {
        fij=c(fsi,fai[j])
        
        message("----- PLS fn=",fn," in ",nx,", i=",i," in ",length(fs),", j=",j," in ",length(fai),". ",date())
        xmi=xm[,fij,drop=FALSE]
        pls=try(opls(x=xmi,y=ym,predI=NA,orthoI=0,permI=rand,fig.pdfC='none',info.txtC='none'))
        if(class(pls)=="try-error")
        {
          pls=try(opls(x=xmi,y=ym,predI=1,orthoI=0,permI=rand,fig.pdfC='none',info.txtC='none'))
        }
        out=dfna
        if(class(pls)!="try-error"){dfi=getSummaryDF(pls);out[1,match(colnames(dfi),colnames(dfna))]=as.vector(as.matrix(dfi[1,]))}
        idni=id
        idni[fij]=1
        outrc[[k]]=cbind(idni,out)
        k=k+1
        R2Yn[[kn]]=out[,"R2Y(cum)"][[1]]
        if(Q2ck){Q2n[[kn]]=out[,"Q2(cum)"][[1]]}
        if(SEEck){SEEn[[kn]]=out[,"RMSEE"][[1]]}
        fijn[[kn]]=fij
        kn=kn+1
      }
    }
    maxR2n=max(unlist(R2Yn),na.rm = TRUE)
    kns=which(unlist(R2Yn)>=(maxR2n*r2buf))
    if(Q2ck)
    {
      if(length(kns)>1)
      {
        Q2nk=Q2n[kns]
        maxQ2nk=max(unlist(Q2nk),na.rm = TRUE)
        if(maxQ2nk>=0){maxQ2nkb=maxQ2nk*r2buf}else{maxQ2nkb=maxQ2nk*(1+(1-r2buf))}
        knsk1=which(unlist(Q2nk)>=maxQ2nkb)
        kns=kns[knsk1]
      }
    }
    if(SEEck)
    {
      if(length(kns)>1)
      {
        SEEnk=SEEn[kns]
        minSEE=min(unlist(SEEnk),na.rm = TRUE)
        knsk2=which(unlist(SEEnk)<=(minSEE*(1+(1-r2buf))))
        kns=kns[knsk2]
      }
    }
    EPS <- (.Machine$double.eps)
    if((max(kns)<=length(fs)) | sum(is.na(kns))>0 | maxR2n >= (1-EPS)){break}else{
      fs=fijn[kns[which(kns>length(fs))]]
      R2Yn=R2Yn[kns]
      fijn=fijn[kns]
      kn=length(R2Yn)+1
    }
  }
  outrcm=Reduce(rbind,outrc)
}

rp.pls<-function(xm,ym,rand=100)
{
  R2sepf<-function(xm,ym,predI)
  {
    pls=try(opls(x=xm,y=ym,predI=predI,orthoI=0,permI=100,info.txtC="none",fig.pdfC="none"))
    if(class(pls)=="try-error"){pls=try(opls(x=xm,y=ym,predI=1,orthoI=0,permI=100,info.txtC="none",fig.pdfC="none"))}
    if(class(pls)=="try-error"){out1=rep(NA,1+ncol(xm))}else{
      out1=c(R2Y=getSummaryDF(pls)[,"R2Y(cum)"],(((getVipVn(pls))^2)/ncol(xm))*getSummaryDF(pls)[,'R2Y(cum)'])
    }
    out1
  }
  pls=try(opls(x=xm,y=ym,predI=NA,orthoI=0,permI=100,info.txtC="none",fig.pdfC="none"))
  predI=NA
  if(class(pls)=="try-error"){predI=1;pls=try(opls(x=xm,y=ym,predI=1,orthoI=0,permI=100,info.txtC="none",fig.pdfC="none"))}
  if(class(pls)=="try-error"){out=rep(NA,2*(1+ncol(xm)))}else{
    R2obs=R2sepf(xm,ym,predI)
    perm=permute::shuffleSet(nrow(xm),nset = rand)
    tracs=seq(from=1,to=nrow(perm),by=20)
    R2rm=sapply(1:nrow(perm),
                function(i)
                {
                  if(i %in% tracs){message("i=",i,". ",date())}
                  ymri=ym[perm[i,],,drop=FALSE]
                  rownames(ymri)=rownames(ym)
                  R2sepf(xm,ymri,predI)
                })
    EPS <- (.Machine$double.eps)
    dR2=((R2rm-R2obs)>=(-EPS))
    out=c(R2obs,rowSums(dR2,na.rm = TRUE)/rand)
  }
  names(out)=c(paste0("R2.",c('Y',colnames(xm))),paste0("P.",c('R2Y',colnames(xm))))
  out
}

r2adj<-function(r2,n,p)
{
  idx=which((n-p-1)<0)
  out=1-((1-r2)*((n-1)/(n-p-1)))
  out[idx]=NA
  out
}

### Read the datasets
xms1=read.table("PLS-average.csv",sep=",",header = TRUE,row.names=1)
head(xms1)
xms1[,2:25]<-scale(xms1[,2:25])

### Take some environmental factors as an example
#####
# 1 # for Soil total organic C
be=xms1[,'TOC',drop=FALSE]
ym=as.matrix(be)
Yname=paste(colnames(ym))

colnames(xms1)
xms2=as.matrix(xms1[,c("Drought","Temperature","moisture","pH","TN","NH4.N","NO3.N","FlTotl","FlC3","PR"), drop=FALSE])
head(xms2)
plstm=plsfw(xms2,ym,r2buf=0.98,Q2ck=TRUE,SEEck=FALSE,rand=100) 
ieggr::save.file(plstm,prefix="PLS",filename = paste0(Yname,".PLSFWtest"))

# Optimuize result: R2Y larger than 98% of maximum R2Y, and P values for R2Y and Q2 should be <0.05; the least factor number; if multiple hits, use minimum RMSEE; if still multiple hits, choose the one with better biological sense.
sel.id=c(1,2,3,5,9)
xmi=xms2[,sel.id,drop=FALSE]
head(xmi)

pls=try(opls(x=xmi,y=ym,predI=NA,orthoI=0,permI=1000,fig.pdfC="none"))
if(class(pls)=="try-error"){pls=try(opls(x=xmi,y=ym,predI=1,orthoI=0,permI=1000,fig.pdfC="none"))}

vip=getVipVn(pls)
names(vip)=paste0("VIP.",names(vip))
rpi=R2sep(xmi,pls)
names(rpi)=paste0("R2.",names(rpi))

R2sig<-Psig<-rep(NA,ncol(xmi))
names(R2sig)<-paste0("R2.Single.",colnames(xmi))
names(Psig)<-paste0("P.Single.",colnames(xmi))
for(j in 1:ncol(xmi))
{
  message("-----Single j=",j,". ",date())
  plsj=try(opls(x=xmi[,j,drop=FALSE],y=ym,predI=NA,orthoI=0,permI=1000,fig.pdfC="none"))
  if(class(plsj)=="try-error"){plsj=try(opls(x=xmi[,j,drop=FALSE],y=ym,predI=1,orthoI=0,permI=1000,fig.pdfC="none"))}
  sdfj=getSummaryDF(plsj)
  R2sig[j]=sdfj[,'R2Y(cum)'][[1]]
  if('pR2Y' %in% colnames(sdfj)){Psig[j]=sdfj[,'pR2Y'][[1]]}
  if(is.na(Psig[j])){rpj=rp.pls(xm=xmi[,j,drop=FALSE],ym=ym,rand = 100);Psig[j]=rpj['P.R2Y']}
}

#####
# 2 # for Soil bacterial PC1
be=xms1[,'B.PC1',drop=FALSE]
ym=as.matrix(be)
Yname=paste(colnames(ym))

colnames(xms1)
xms2=as.matrix(xms1[,c("Drought","Temperature","moisture","pH",
                       "TOC","TN","NH4.N","NO3.N","FlTotl",
                       "FlC3","PR","X16S","ITS","X18S","F.PC1",
                       "P.PC1","Carbon.degradation","Carbon.Fixation",
                       "Nitrification","Assimilatory.N..reduction","Dissimilatory.N.reduction","Stress"), drop=FALSE])
head(xms2)
plstm=plsfw(xms2,ym,r2buf=0.98,Q2ck=TRUE,SEEck=FALSE,rand=100) 
ieggr::save.file(plstm,prefix="PLS",filename = paste0(Yname,".PLSFWtest"))

# Optimuize result: R2Y larger than 98% of maximum R2Y, and P values for R2Y and Q2 should be <0.05; the least factor number; if multiple hits, use minimum RMSEE; if still multiple hits, choose the one with better biological sense.
sel.id=c(2,4,16,19,21,22)
xmi=xms2[,sel.id,drop=FALSE]
head(xmi)

pls=try(opls(x=xmi,y=ym,predI=NA,orthoI=0,permI=1000,fig.pdfC="none"))
if(class(pls)=="try-error"){pls=try(opls(x=xmi,y=ym,predI=1,orthoI=0,permI=1000,fig.pdfC="none"))}

vip=getVipVn(pls)
names(vip)=paste0("VIP.",names(vip))
rpi=R2sep(xmi,pls)
names(rpi)=paste0("R2.",names(rpi))

R2sig<-Psig<-rep(NA,ncol(xmi))
names(R2sig)<-paste0("R2.Single.",colnames(xmi))
names(Psig)<-paste0("P.Single.",colnames(xmi))
for(j in 1:ncol(xmi))
{
  message("-----Single j=",j,". ",date())
  plsj=try(opls(x=xmi[,j,drop=FALSE],y=ym,predI=NA,orthoI=0,permI=1000,fig.pdfC="none"))
  if(class(plsj)=="try-error"){plsj=try(opls(x=xmi[,j,drop=FALSE],y=ym,predI=1,orthoI=0,permI=1000,fig.pdfC="none"))}
  sdfj=getSummaryDF(plsj)
  R2sig[j]=sdfj[,'R2Y(cum)'][[1]]
  if('pR2Y' %in% colnames(sdfj)){Psig[j]=sdfj[,'pR2Y'][[1]]}
  if(is.na(Psig[j])){rpj=rp.pls(xm=xmi[,j,drop=FALSE],ym=ym,rand = 100);Psig[j]=rpj['P.R2Y']}
}

#####
# 3 # for ER
be=xms1[,'ER',drop=FALSE]
ym=as.matrix(be)
Yname=paste(colnames(ym))

colnames(xms1)
xms2=as.matrix(xms1[,c("Drought","Temperature","moisture","pH",
                       "TOC","TN","NH4.N","NO3.N","FlTotl",
                       "FlC3","PR","X16S","ITS","X18S","B.PC1",
                       "F.PC1","P.PC1","Carbon.degradation",
                       "Carbon.Fixation","Nitrification","Assimilatory.N..reduction",
                       "Dissimilatory.N.reduction",'Stress',"NEE"), drop=FALSE])
head(xms2)
plstm=plsfw(xms2,ym,r2buf=0.98,Q2ck=TRUE,SEEck=TRUE,rand=100) 
ieggr::save.file(plstm,prefix="PLS",filename = paste0(Yname,".PLSFWtest"))

# Optimuize result: R2Y larger than 98% of maximum R2Y, and P values for R2Y and Q2 should be <0.05; the least factor number; if multiple hits, use minimum RMSEE; if still multiple hits, choose the one with better biological sense.
sel.id=c(1,2,6,8,12,17,19,20,21,22,23)
sel.id=c(3,8,11,13,18,20,21,22,24)
xmi=xms2[,sel.id,drop=FALSE]
head(xmi)

pls=try(opls(x=xmi,y=ym,predI=NA,orthoI=0,permI=1000,fig.pdfC="none"))
if(class(pls)=="try-error"){pls=try(opls(x=xmi,y=ym,predI=1,orthoI=0,permI=1000,fig.pdfC="none"))}

vip=getVipVn(pls)
names(vip)=paste0("VIP.",names(vip))
rpi=R2sep(xmi,pls)
names(rpi)=paste0("R2.",names(rpi))

R2sig<-Psig<-rep(NA,ncol(xmi))
names(R2sig)<-paste0("R2.Single.",colnames(xmi))
names(Psig)<-paste0("P.Single.",colnames(xmi))
for(j in 1:ncol(xmi))
{
  message("-----Single j=",j,". ",date())
  plsj=try(opls(x=xmi[,j,drop=FALSE],y=ym,predI=NA,orthoI=0,permI=1000,fig.pdfC="none"))
  if(class(plsj)=="try-error"){plsj=try(opls(x=xmi[,j,drop=FALSE],y=ym,predI=1,orthoI=0,permI=1000,fig.pdfC="none"))}
  sdfj=getSummaryDF(plsj)
  R2sig[j]=sdfj[,'R2Y(cum)'][[1]]
  if('pR2Y' %in% colnames(sdfj)){Psig[j]=sdfj[,'pR2Y'][[1]]}
  if(is.na(Psig[j])){rpj=rp.pls(xm=xmi[,j,drop=FALSE],ym=ym,rand = 100);Psig[j]=rpj['P.R2Y']}
}

#############Repeat the above analyses for the remaining parameters