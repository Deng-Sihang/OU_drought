rm(list=ls())
library(lme4)
library(car)
library(plyr)
library(vegan)
library(dplyr)
library(ieggr)
library(reshape2)
save.wd <- iwd(choose.dir())

### Calculation of beta diversity indices -------------------------------------------------------------------
### Taking bacterial communities as an example for calculation
### For the phyla/classes/guilds obtained after annotation, the calculation process is the same
com.file = "Bacteria-zotu-resampled.csv"
treat.file = "Treatment info-Microbial biomass.csv"
tree.file = "Bacteria-tree.nwk"

comm = read.csv(com.file, row.names = 1, header = T, sep = ",")
### Columns 49-55 are annotated information for species
comm = comm[,1:48]
treat = read.csv(treat.file, row.names = 1, header = T, sep = ",")
comm = as.data.frame(t(comm))

dim(comm)
comm = comm[, colSums(comm) > 0]
dim(comm)

dim(comm)
comm = comm[rowSums(comm) > 0,]
dim(comm)

samp = match.name(rn.list = list(treat = treat, comm = comm))
comm = samp$comm
treat = samp$treat

tree = lazyopen(tree.file)
spc = match.name(cn.list = list(comm = comm), tree.list = list(tree = tree))
comm = spc$comm
tree = spc$tree

### Sorensen distances
dist.soren = vegdist(comm, binary = TRUE)
### unweighted Unifrac distances
unifc = unifrac.g(otu.tab = comm, tree = tree,alpha = c(0, 0.5, 1))
dist.uw = unifc$d_UW %>%
  as.dist()

### Non-parametric test based on Sorensen distances ------------------------------------------------------
### whole community ###
adonis.soren = adonis2(dist.soren ~ Treatment+ block*Year, data = treat, permutations = 999)
anosim.soren = anosim(dist.soren,treat$Treatment,permutations = 999,
                     strata = as.factor(paste(treat$Year,treat$block,sep = ".")))
summary(anosim.soren)
mrpp.soren = mrpp(dist.soren,treat$Treatment,permutations = 999,
                 strata = as.factor(paste(treat$Year,treat$block,sep = ".")))

### For a specific year, you can use the following code###
### Taking the year 2009 as an example:
dist.soren = vegdist(comm[1:8,], binary = TRUE)
treat = treat[1:8,]
adonis.soren = adonis2(dist.soren ~ Treatment+ block, data = treat, permutations = 999)
anosim.soren = anosim(dist.soren,treat$Treatment,permutations = 999,
                      strata = as.factor(paste(treat$block,sep = ".")))
summary(anosim.soren)
mrpp.soren = mrpp(dist.soren,treat$Treatment,permutations = 999,
                  strata = as.factor(paste(treat$block,sep = ".")))

### NMDS analysis based on Sorensen distances--------------------------------------------------------------------
nmds_dis <- metaMDS(dist.soren, k = 2)
nmds_dis$stress
nmds_dis_site <- data.frame(nmds_dis$points)
nmds_dis_site$Treatment <- treat$Treatment[match(rownames(nmds_dis_site),rownames(treat))]
nmds_dis_site$Year <- treat$year[match(rownames(nmds_dis_site),rownames(treat))]
write.csv(nmds_dis_site,"NMDS-score-Bacteria.csv")

### Temporal succession of microbial communities' dispersion --------------------------------------------------------------
### Calculation of dispersion based on Sorensen distances
for (i in 1:6) ### 6 represents the 6-year prolonged drought.
{
    time = i+2008
    dist.soren1 = vegdist(comm[(8*i-7):(8*i),], binary = TRUE)
    list.beta <- betadisper(dist.soren1, group = subset(treat, year == time)$Treatment,type = c("centroid"))
    disper1 = as.data.frame(list.beta$distances)
    if (i == 1)
    {
      disper.soren = NULL
      disper.soren = rbind(disper.soren,disper1)
    }else {
      disper.soren = rbind(disper.soren,disper1)
    }
}

### Calculation of dispersion based on unweighted Unifrac distances
for (i in 1:6)
{
    time = i+2008
    comm1 = comm[(8*i-7):(8*i),]
    dim(comm1)
    comm1 = comm1[, colSums(comm1) > 0]
    dim(comm1)
    
    spc1 = match.name(cn.list = list(comm1 = comm1), tree.list = list(tree = tree))
    comm1 = spc1$comm1
    tree1 = spc1$tree
    
    unifc = unifrac.g(otu.tab = comm1, tree = tree1,alpha = c(0, 0.5, 1))
    dist.uw1 = unifc$d_UW %>%
      as.dist()
    list.beta <- betadisper(dist.uw1, group = subset(treat, year == time)$Treatment,type = c("centroid"))
    disper1 = as.data.frame(list.beta$distances)
    if (i == 1)
    {
      disper.uw = NULL
      disper.uw = rbind(disper.uw,disper1)
    }else {
      disper.uw = rbind(disper.uw,disper1)
    }
}

### Check if the sample names are matched  
### based on Sorensen distances
samp = match.name(rn.list = list(treat = treat, disper.soren = disper.soren))
disper.soren = samp$disper.soren
treat = samp$treat
dispersion.soren <-data.frame(disper=disper.soren[,1],treat)

### based on unweighted Unifrac distances
samp = match.name(rn.list = list(treat = treat, disper.uw = disper.uw))
disper.uw = samp$disper.uw
treat = samp$treat
dispersion.uw <-data.frame(disper=disper.uw[,1],treat)

### Temporal succession (permutation test)
disper.lmm <- function(dispersion,treat)
{
    drought=treat$Treatment
    drought.lev=unique(drought)
    out=list()
    for(j in 1:length(drought.lev))
    {
      idj=which(drought==drought.lev[j])
      dispersion1 = dispersion[idj,]
      lmij=lmer(disper~year+((1+year)|plot),data=dispersion1)
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

### Taking dispersion.soren as an example
disperi =disper.lmm(dispersion = dispersion.soren,treat = treat)
slope.obs = as.vector(disperi[1,])
r2.obs=as.vector(disperi[3:4,])
aic.obs=as.vector(disperi[5:6,])
ds.obs=(-disperi[1,1])-(-disperi[1,2])

### randomize time points and other the same as observed. 
year.lev=unique(treat$year)
rand =1000
year.perm=vegan:::getPermuteMatrix(rand,length(year.lev))
trace.seq=seq(from=1,to=rand,by = 100)
ind.rand=lapply(1:nrow(year.perm),
                  function(k)
                  {
                    if(k %in% trace.seq) message("-------Now randomizing k=",k,". ",date())
                    out=list()
                    idi=year.perm[k,]
                    dispersionr=dispersion.soren
                    dispersionr[,"year"]=year.lev[idi[match(dispersion.soren$year,year.lev)]]
                    disperc=disper.lmm(dispersion = dispersionr,treat = treat)
                    out$slop.fix=as.vector(disperc[1,])
                    out$r2=as.vector(disperc[3:4,])
                    out$aic=as.vector(disperc[5:6,])
                    out$ds=(-disperc[1,1])-(-disperc[1,2])
                    out
                  })
slop.ran =sapply(1:length(ind.rand),function(k){ind.rand[[k]]$slop.fix})
r2.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$r2})
aic.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$aic})
ds.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$ds})
p.ran = rbind(slop.ran,r2.ran,aic.ran,ds.ran) 
EPS <- sqrt(.Machine$double.eps)
p.slope=(rowSums(slop.ran>=(matrix(slope.obs,nr=nrow(slop.ran),nc=ncol(slop.ran))+EPS))+1)/(ncol(slop.ran)+1)
p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS))+1)/(ncol(r2.ran)+1)
p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS))+1)/(ncol(aic.ran)+1)
if(ds.obs>0)
  {
    p.perm=(sum(ds.ran>=(ds.obs-EPS))+1)/(length(ds.ran)+1)
  }else{
    p.perm=(sum(ds.ran<=(ds.obs+EPS))+1)/(length(ds.ran)+1)
  }
p.values=rbind(matrix(p.slope,1,2),matrix(p.r2,2,2),matrix(p.aic,2,2),c(p.perm,NA))
rownames(p.values)=c("P.Slope","P.R2M","P.R2C","P.AIC1","P.AIC2","P.dS.perm")
output=rbind(disperi,p.values)

### Time-decay relationships (TDRs) -------------------------------------------------
### Taking the calculation based on Sorensen distances as an example
betai <- as.matrix(dist.soren)
prefixi = "Bacteria.Sorensen"
tdc.lmm<-function(betai,treat,save.output=FALSE,prefixi=NULL)
{
  drought=treat$Treatment
  drought.lev=unique(drought)
  out=list()
  for(j in 1:length(drought.lev))
  {
    idj=which(drought==drought.lev[j])
    sampj=rownames(treat)[idj]
    betaij=betai[idj,idj]
    betaij3=dist.3col(betaij)
    dtj=abs(treat$year[match(betaij3[,1],rownames(treat))]-treat$year[match(betaij3[,2],rownames(treat))])
    plotj1=treat$plot[match(betaij3[,1],rownames(treat))]
    plotj2=treat$plot[match(betaij3[,2],rownames(treat))]
    idj.use=which(plotj1==plotj2)
    cij.use=1-betaij3[idj.use,3]
    dtj.use=dtj[idj.use]
    plotj.use=plotj1[idj.use]
    logCij=log(cij.use)
    logCij[logCij==-Inf]<-NA
    logdtj=log(dtj.use)
    logdtj[logdtj==-Inf]<-NA
    lmij=lmer(logCij~logdtj+((1+logdtj)|plotj.use))
    lmijsm=summary(lmij)
    AIC1=AIC(lmij)
    r2ij=rsquared.glmm(lmij)
    lmijCS=Anova(lmij,type = "II")
    if(save.output)
    {
      dataij=data.frame(treat=drought.lev[j],betaij3[idj.use,,drop=FALSE],
                        time.diff=dtj.use,similarity=cij.use,
                        ln.timediff=logdtj,ln.S=logCij,block=blockj.use)
      save.file(dataij,prefix = prefixi,filename = paste0("TimeDecay.Data.",drought.lev[j]))
    }
    out[[j]]=c(slop.fix=lmijsm$coefficients[2,1],slop.sd=coef(summary(lmij))[ , "Std. Error"][2],
               R2M=r2ij$Marginal,R2C=r2ij$Conditional,AIC1=AIC1,AIC2=r2ij$AIC,
               P.typeII=lmijCS[[3]],Chisq=lmijCS[[1]])
  }
  outs=Reduce(cbind,out)
  colnames(outs)=drought.lev
  outs
}

tdci=tdc.lmm(betai=betai,treat = treat,save.output = FALSE,prefixi = prefixi)
slope.obs = as.vector(tdci[1,])
r2.obs=as.vector(tdci[3:4,])
aic.obs=as.vector(tdci[5:6,])
ds.obs=(-tdci[1,1])-(-tdci[1,2])

### randomize time points and other the same as observed.
year.lev=unique(treat$year)
rand = 1000
year.perm=vegan:::getPermuteMatrix(rand,length(year.lev))
trace.seq=seq(from=1,to=rand,by = 100)
ind.rand=lapply(1:nrow(year.perm),
                function(k)
                {
                  if(k %in% trace.seq) message("-------Now randomizing k=",k,". ",date())
                  out=list()
                  idi=year.perm[k,]
                  perm.treat=treat
                  perm.treat[,"year"]=year.lev[idi[match(treat$year,year.lev)]]
                  tdcr=tdc.lmm(betai = betai, treat = perm.treat,save.output = FALSE,prefixi = NULL)
                  out$slop.fix=as.vector(tdcr[1,])
                  out$r2=as.vector(tdcr[3:4,])
                  out$aic=as.vector(tdcr[5:6,])
                  out$ds=(-tdcr[1,1])-(-tdcr[1,2])
                  out
                })
slop.ran =sapply(1:length(ind.rand),function(k){ind.rand[[k]]$slop.fix})
r2.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$r2})
aic.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$aic})
ds.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$ds})
p.ran = rbind(slop.ran,r2.ran,aic.ran,ds.ran) 
EPS <- sqrt(.Machine$double.eps)
p.slope=(rowSums(slop.ran<=(matrix(slope.obs,nr=nrow(slop.ran),nc=ncol(slop.ran))+EPS))+1)/(ncol(slop.ran)+1)
p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS))+1)/(ncol(r2.ran)+1)
p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS))+1)/(ncol(aic.ran)+1)
if(ds.obs>0)
{
  p.perm=(sum(ds.ran>=(ds.obs-EPS))+1)/(length(ds.ran)+1)
}else{
  p.perm=(sum(ds.ran<=(ds.obs+EPS))+1)/(length(ds.ran)+1)
}
p.values=rbind(matrix(p.slope,1,2),matrix(p.r2,2,2),matrix(p.aic,2,2),c(p.perm,NA))
rownames(p.values)=c("P.Slope","P.R2M","P.R2C","P.AIC1","P.AIC2","P.dS.perm")
output=rbind(disperi,p.values)

### Temporal change in community differences (dissimilarity) ------------------------------------------
### Taking the calculation based on Sorensen distances as an example
betai <- as.matrix(dist.soren)
dv.lmm<-function(betai,treat,rand=1000,prefixi,save.data=FALSE)
  {
  betai3=dist.3col(betai)
  id1=match(betai3[,1],rownames(treat))
  id2=match(betai3[,2],rownames(treat))
  w1=treat$Treatment[id1]
  w2=treat$Treatment[id2]
  y1=treat$year[id1]
  y2=treat$year[id2]
  p1=treat$block[id1]
  p2=treat$block[id2]
  bys=paste0(treat$block,".",treat$year)
  by1=bys[id1]
  by2=bys[id2]
  id.use=which((by1==by2)&(w1!=w2))
  d.use=betai3[id.use,3]
  t.use=y1[id.use]
  pp.use=sapply(id.use,function(x){paste(sort(c(p1[x],p2[x])),collapse = "")})
  data.out=data.frame(betai3[id.use,1:2,drop=FALSE],t=t.use,D.A=d.use,
                      rand.eff=pp.use,stringsAsFactors = FALSE)
  if(save.data)  save.file(data.out,prefix=prefixi,filename = "Dissimilarity")
  lmi=list()
  lmi[[1]]=lmer(D.A~t+((1+t)|rand.eff),data = data.out)
  out1=sapply(1:length(lmi),
              function(k)
              {
                lmism=summary(lmi[[k]])
                AIC1=AIC(lmi[[k]])
                r2i=rsquared.glmm(lmi[[k]])
                lmiCS=Anova(lmi[[k]],type = "II")
                c(slop.fix=lmism$coefficients[2,1],slop.sd=coef(lmism)[ , "Std. Error"][2],
                  R2M=r2i$Marginal,R2C=r2i$Conditional,AIC1=AIC1,AIC2=r2i$AIC,
                  P.typeII=lmiCS[[3]],Chisq=lmiCS[[1]])
              })
  out1
}
  
dvi = dv.lmm(betai=betai,treat=treat,rand=1000,prefixi = "Bacteria.Sorensen",save.data=FALSE) 
slope.obs = as.vector(dvi[1,])
r2.obs=as.vector(dvi[3:4,])
aic.obs=as.vector(dvi[5:6,])

rand = 1000
year.lev=unique(treat$year)
year.perm=vegan:::getPermuteMatrix(rand,length(year.lev))
trace.seq=seq(from=1,to=rand,by = 100)
ind.rand=lapply(1:nrow(year.perm),
                  function(k)
                  {
                    if(k %in% trace.seq) message("-------Now randomizing k=",k,". ",date())
                    out=list()
                    idi=year.perm[k,]
                    perm.treat=treat
                    perm.treat[,"year"]=year.lev[idi[match(treat$year,year.lev)]]
                    dvr=dv.lmm(betai=betai,treat=perm.treat,rand=1000,
                               prefixi = NULL,save.data=FALSE)
                    out$slop.fix=as.vector(dvr[1,])
                    out$r2=as.vector(as.matrix(dvr[3:4,]))
                    out$aic=as.vector(as.matrix(dvr[5:6,]))
                    out
                  })
slop.ran =t(as.matrix(sapply(1:length(ind.rand),function(k){ind.rand[[k]]$slop.fix})))
r2.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$r2})
aic.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$aic})
p.ran = rbind(slop.ran,r2.ran,aic.ran)
EPS <- sqrt(.Machine$double.eps)
p.slope=(rowSums(slop.ran>=(matrix(slope.obs,nr=nrow(slop.ran),nc=ncol(slop.ran))-EPS))+1)/(ncol(slop.ran)+1)
p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS))+1)/(ncol(r2.ran)+1)
p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS))+1)/(ncol(aic.ran)+1)
p.values=rbind(matrix(p.slope,1,1),matrix(p.r2,2,1),matrix(p.aic,2,1))
rownames(p.values)=c("P.slope","P.R2M","P.R2C","P.AIC1","P.AIC2")
output=rbind(dvi,p.values)