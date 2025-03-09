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

### Temporal succession of MST ----------------------------------------------------------
### linear mixed-effects models + permutation tests
treatused<-read.csv("Treatment info-Microbial biomass.csv",header = T,row.names = 1)
divindex<-read.csv("Bacteria-Fungi-Protists-MSTvalue.csv",header = T)
treatused$year<-treatused$year-2009 

### Match the MST data with the corresponding sample information
divindex$treatment <- treatused$Treatment[match(divindex$name1,rownames(treatused))]
divindex$block1 <- treatused$block[match(divindex$name1,rownames(treatused))]
divindex$block2 <- treatused$block[match(divindex$name2,rownames(treatused))]
divindex$block <- paste(divindex$block1,divindex$block2)
divindex[,5:7]<-scale(divindex[,5:7])

### Columns 5 to 7 correspond to the MST values for bacteria, fungi, and protists, respectively
for(i in 5:7)
{
  div<-data.frame(divtest=divindex[,i],divindex)
  
  ### MST.lmm: function used to calculate the MST trend over time
  MST.lmm <- function(MSTi,treat)
  {
    drought=MSTi$treatment
    drought.lev=unique(drought)
    out=list()
    for(j in 1:length(drought.lev))
    {
      idj=which(drought==drought.lev[j])
      MSTi1 = MSTi[idj,]
      lmij=lmer(divtest~year+((1+year)|block),data=MSTi1)
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
  
  MSTi=MST.lmm(MSTi=div,treat = treatused)
  slope.obs = as.vector(MSTi[1,])
  r2.obs=as.vector(MSTi[3:4,])
  aic.obs=as.vector(MSTi[5:6,])
  ds.obs=(-MSTi[1,1])-(-MSTi[1,2])
  
  ### Permutation test
  year.lev=unique(treatused$year)
  rand =1000
  year.perm=vegan:::getPermuteMatrix(rand,length(year.lev))
  trace.seq=seq(from=1,to=rand,by = 100)
  ind.rand=lapply(1:nrow(year.perm),
                  function(k)
                  {
                    if(k %in% trace.seq) message("-------Now randomizing k=",k,". ",date())
                    out=list()
                    idi=year.perm[k,]
                    div1=div
                    div1[,"year"]=year.lev[idi[match(div$year,year.lev)]]
                    MSTr=MST.lmm(MSTi=div1,treat = treatused)
                    out$slop.fix=as.vector(MSTr[1,])
                    out$r2=as.vector(MSTr[3:4,])
                    out$aic=as.vector(MSTr[5:6,])
                    out$ds=(-MSTr[1,1])-(-MSTr[1,2])
                    out
                  })
  slop.ran =sapply(1:length(ind.rand),function(k){ind.rand[[k]]$slop.fix})
  r2.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$r2})
  aic.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$aic})
  ds.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$ds})
  p.ran = rbind(slop.ran,r2.ran,aic.ran,ds.ran) 
  
  if(slope.obs[1]<0){for (k in 1:719){if(p.ran[1,k]>0){p.ran[3,k]<-NA;p.ran[4,k]<-NA;p.ran[7,k]<-NA;p.ran[8,k]<-NA}}
  }else{for (k in 1:719){if(p.ran[1,k]<0){p.ran[3,k]<-NA;p.ran[4,k]<-NA;p.ran[7,k]<-NA;p.ran[8,k]<-NA}}}
  if(slope.obs[2]<0){for (l in 1:719){if(p.ran[2,l]>0){p.ran[5,l]<-NA;p.ran[6,l]<-NA;p.ran[9,l]<-NA;p.ran[10,l]<-NA}}
  }else{for (l in 1:719){if(p.ran[2,l]<0){p.ran[5,l]<-NA;p.ran[6,l]<-NA;p.ran[9,l]<-NA;p.ran[10,l]<-NA}}}
  if((slope.obs[1]<0)&(slope.obs[2]<0)){for (m in 1:719){if((p.ran[1,m]>0)|(p.ran[2,m]>0)){p.ran[11,m]<-NA}}}
  if((slope.obs[1]<0)&(slope.obs[2]>0)){for (m in 1:719){if((p.ran[1,m]>0)|(p.ran[2,m]<0)){p.ran[11,m]<-NA}}}
  if((slope.obs[1]>0)&(slope.obs[2]<0)){for (m in 1:719){if((p.ran[1,m]<0)|(p.ran[2,m]>0)){p.ran[11,m]<-NA}}}
  if((slope.obs[1]>0)&(slope.obs[2]>0)){for (m in 1:719){if((p.ran[1,m]<0)|(p.ran[2,m]<0)){p.ran[11,m]<-NA}}}
  
  EPS <- sqrt(.Machine$double.eps)
  slop.ran.1 = p.ran[1,]
  slop.ran.2 = p.ran[2,]
  r2.ran = p.ran[3:6,]
  aic.ran = p.ran[7:10,]
  ds.ran = p.ran[11,]
  
  if(slope.obs[1]<0){
    p.slope.1=(sum(slop.ran.1<=(slope.obs[1]+EPS))+1)/(length(slop.ran.1)+1)
  }else{
    p.slope.1=(sum(slop.ran.1>=(slope.obs[1]-EPS))+1)/(length(slop.ran.1)+1)
  }
  if(slope.obs[2]<0){
    p.slope.2=(sum(slop.ran.2<=(slope.obs[2]+EPS))+1)/(length(slop.ran.2)+1)
  }else{
    p.slope.2=(sum(slop.ran.2>=(slope.obs[2]-EPS))+1)/(length(slop.ran.2)+1)
  }
  
  p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS),na.rm=TRUE)+1)/(ncol(r2.ran)+1)
  p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS),na.rm=TRUE)+1)/(ncol(aic.ran)+1)
  if(ds.obs>0)
  {
    p.perm=(sum(ds.ran>=(ds.obs-EPS),na.rm=TRUE)+1)/(length(ds.ran)+1)
  }else{
    p.perm=(sum(ds.ran<=(ds.obs+EPS),na.rm=TRUE)+1)/(length(ds.ran)+1)
  }
  
  ### Results output
  p.values=rbind(c(p.slope.1,p.slope.2),matrix(p.r2,2,2),matrix(p.aic,2,2),c(p.perm,NA))
  rownames(p.values)=c("P.Slope","P.R2M","P.R2C","P.AIC1","P.AIC2","P.dS.perm")
  
  output=rbind(MSTi,p.values)
  colnames(output)=c((paste0(colnames(divindex)[i],".Drought")),(paste0(colnames(divindex)[i],".Control")))
  if (i == 5){
    result = NULL
    result <- cbind(result, output)
  }else{
    result <- cbind(result, output)
  }
  
  result <- as.data.frame(result)
  result
}
write.csv(result,"AllassemblyMST-LMM-Permutation.csv")

### Correlations between deterministic processes and environmental variables --------------------------------------------
Env<-read.csv("Env-factor.csv",header = T,row.names = 1)
MST<-read.csv("Bacteria-Fungi-Protists-MSTvalue.csv",header = T)
treatused<-read.csv("Treatment info-Microbial biomass.csv",header = T,row.names = 1)
### Maintain key environmental factors
Env<-Env[,-(13:15)]
Env<-Env[,-(16:18)]

### Datasets pre-treatment
Env<-log10(Env)
treatused$year <- treatused$year-2009
### Proportion of deterministic processes = 1 - proportion of stochastic processes
MST$Bacteria <- 1-MST$Bacteria
MST$Fungi <- 1-MST$Fungi
MST$Protists <- 1-MST$Protists

### Calculate the correlations
for (i in 1:ncol(Env))
{
### Calculation of Euclidean distance between samples of environmental factor
env.dist <- vegdist(Env[,i],"euclidean")
env.dist <- as.matrix(env.dist)
row.names(env.dist) <- row.names(Env)
colnames(env.dist) <- row.names(Env)
env.distij3=dist.3col(env.dist)
env.distij3$treatment1 <- treatused$Treatment[match(env.distij3$name1,rownames(treatused))]
env.distij3$treatment2 <- treatused$Treatment[match(env.distij3$name2,rownames(treatused))]
env.distij3$year1 <- treatused$year[match(env.distij3$name1,rownames(treatused))]
env.distij3$year2 <- treatused$year[match(env.distij3$name2,rownames(treatused))]

### only drought with same year
env.distij3.drought <- env.distij3[(env.distij3$treatment1=="Drought")&(env.distij3$treatment2=="Drought")&(env.distij3$year1==env.distij3$year2),]
MST.drought <- MST[MST$group=="Drought",]
MST.drought$env <- env.distij3.drought$dis

### linear mixed-effect regression
MST.drought$block1 <- treatused$block[match(MST.drought$name1,rownames(treatused))]
MST.drought$block2 <- treatused$block[match(MST.drought$name2,rownames(treatused))]
MST.drought$block <- paste(MST.drought$block1,MST.drought$block2)

for (j in 5:7)
{
  lmij<-lmer(MST.drought[,j]~env+(1|block)+(1|year),data=MST.drought)
  lmijsm=summary(lmij)
  AIC1=AIC(lmij)
  r2ij=rsquared.glmm(lmij)
  lmijCS=Anova(lmij,type = "II")
  out=c(slop.fix=lmijsm$coefficients[2,1],slop.sd=coef(summary(lmij))[ , "Std. Error"][2],
        R2M=r2ij$Marginal,R2C=r2ij$Conditional,AIC1=AIC1,AIC2=r2ij$AIC,
        P.typeII=lmijCS[[3]],Chisq=lmijCS[[1]])
  outs=as.data.frame(out)
  colnames(outs)<-paste0(colnames(Env)[i],".",colnames(MST.drought)[j])
  if (j == 5){
    result = NULL
    result <- outs
  }else{
    result <- cbind(result, outs)
  }
}

if (i == 1){
  output = NULL
  output <- result
}else{
  output <- cbind(output, result)
}

output
}

write.csv(output,"AllassemblyMST-LMM-CorrelationEnv-Control.csv")

### Geochip analysis --------------------------------------------
treatused<-read.csv("Treatment info-Microbial biomass.csv",header = T,row.names = 1)
divindex<-read.csv("GeoChip-rawdata-2009-2014-Gene.csv",header = T,row.names = 1)

### Dataset pre-treatment
divindex<-divindex[,1:48]
counts<-rowSums(divindex>0)
divindex<-divindex[counts>=24,]
divindex<-as.data.frame(t(divindex))
treatused$year<-treatused$year-2009 

### match names
divindex<-divindex[match(row.names(treatused),row.names(divindex)),] 
### check if names are matched
sum(row.names(divindex)==row.names(treatused)) 
### rescale the diversities
divindex<-scale(divindex) 

### Calculate the effect sizes of treatment on relative abundance of specific probes
divs1<-sapply(1:ncol(divindex),function(j){
  message("Now j=",j," in ",ncol(divindex),". ",date())
  if (length(unique(divindex[,j]))<3){
    result<-rep(NA,8)
  } else {
    div<-data.frame(divtest=divindex[,j],treatused)
    
    fm1<-lmer(divtest~Drought+(1|block)+(1|year),data=div)
    
    presult<-car::Anova(fm1,type=2)
    coefs<-coef(summary(fm1))[ , "Estimate"]  ##four coefs
    names(coefs)<-paste0(names(coefs),".mean")
    
    SEvalues<-coef(summary(fm1))[ , "Std. Error"] ##standard errors
    names(SEvalues)<-paste0(names(SEvalues),".se")
    
    tvalues<-coef(summary(fm1))[ , "t value"] ##t values
    names(tvalues)<-paste0(names(tvalues),".t")
    
    chisqP<-c(presult[,1],presult[,3])
    names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
    
    result<-c(coefs,tvalues,SEvalues,chisqP)}
  result
})

colnames(divs1)<-colnames(divindex)
divs1 <- as.data.frame(t(divs1))

### fdr adjustment
colSums(divs1<0.05)
divs1$p.adjust <- p.adjust(divs1$Drought.P,method="fdr")
colSums(divs1<0.05)
write.csv(divs1,"GeoChip-rawdata-2009-2014-Gene-LMM.csv")

### PLS --------------------------------------------
############
### install.packages("OmicsPLS")
library(OmicsPLS)
### if (!requireNamespace("BiocManager", quietly = TRUE))
### install.packages("BiocManager")
### BiocManager::install("ropls")
library(ropls)
#################
#################
### Loading functions
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
xms1[,2:24]<-scale(xms1[,2:24])

### Take some environmental factors as an example
#####
### 1 ### for Soil total organic C
be=xms1[,'TOC',drop=FALSE]
ym=as.matrix(be)
Yname=paste(colnames(ym))

colnames(xms1)
xms2=as.matrix(xms1[,c("Drought","Temperature","moisture","pH","TN","NH4.N","NO3.N","FlTotl","FlC3"), drop=FALSE])
head(xms2)
plstm=plsfw(xms2,ym,r2buf=0.98,Q2ck=TRUE,SEEck=FALSE,rand=100) 
ieggr::save.file(plstm,prefix="PLS",filename = paste0(Yname,".PLSFWtest"))

### Optimuize result: R2Y larger than 98% of maximum R2Y, and P values for R2Y and Q2 should be <0.05; the least factor number; if multiple hits, use minimum RMSEE; if still multiple hits, choose the one with better biological sense.
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
### 2 ### for Soil bacterial PC1
be=xms1[,'B.PC1',drop=FALSE]
ym=as.matrix(be)
Yname=paste(colnames(ym))

colnames(xms1)
xms2=as.matrix(xms1[,c("Drought","Temperature","moisture","pH",
                       "TOC","TN","NH4.N","NO3.N","FlTotl",
                       "FlC3","X16S","ITS","X18S","F.PC1",
                       "P.PC1","Carbon.degradation","Carbon.Fixation",
                       "Nitrification","Assimilatory.N..reduction","Dissimilatory.N.reduction","Stress"), drop=FALSE])
head(xms2)
plstm=plsfw(xms2,ym,r2buf=0.98,Q2ck=TRUE,SEEck=TRUE,rand=100) 
ieggr::save.file(plstm,prefix="PLS",filename = paste0(Yname,".PLSFWtest"))

### Optimuize result: R2Y larger than 98% of maximum R2Y, and P values for R2Y and Q2 should be <0.05; the least factor number; if multiple hits, use minimum RMSEE; if still multiple hits, choose the one with better biological sense.
sel.id=c(1,4,8,13,15,17,18,21)
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
### 3 ### for ER
be=xms1[,'ER',drop=FALSE]
ym=as.matrix(be)
Yname=paste(colnames(ym))

colnames(xms1)
xms2=as.matrix(xms1[,c("Drought","Temperature","moisture","pH",
                       "TOC","TN","NH4.N","NO3.N","FlTotl",
                       "FlC3","X16S","ITS","X18S","B.PC1",
                       "F.PC1","P.PC1","Carbon.degradation",
                       "Carbon.Fixation","Nitrification","Assimilatory.N..reduction",
                       "Dissimilatory.N.reduction",'Stress',"NEE"), drop=FALSE])
head(xms2)
plstm=plsfw(xms2,ym,r2buf=0.98,Q2ck=TRUE,SEEck=TRUE,rand=100) 
ieggr::save.file(plstm,prefix="PLS",filename = paste0(Yname,".PLSFWtest"))

### Optimuize result: R2Y larger than 98% of maximum R2Y, and P values for R2Y and Q2 should be <0.05; the least factor number; if multiple hits, use minimum RMSEE; if still multiple hits, choose the one with better biological sense.
sel.id=c(1,2,3,8,10,11,12,13,17,18,19,20,23)
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