rm(list=ls())
library(lme4)
library(car)
library(plyr)
library(vegan)
library(dplyr)
library(ieggr)
library(reshape2)
library(boot)
save.wd <- iwd(choose.dir())

### Calculation of beta diversity indices -------------------------------------------------------------------
### Taking bacterial communities as an example for calculation
### For the phyla/classes/guilds obtained after annotation, the calculation process is the same
com.file = "Bacteria-zotu-resampled.csv"
treat.file = "Treatment info-Microbial biomass.csv"
tree.file ="Bacteria-tree.nwk"
comm = read.csv(com.file, row.names = 1, header = T, sep = ",")

### For different phylogenetic lineages and functional guilds
### After column 49, Listed as annotated information for species
### Bacteria
#comm <- comm[comm$Phylum=="D_1__Verrucomicrobia",1:48]
#comm <- comm[comm$Class=="D_2__Deltaproteobacteria",1:48]
#comm <- comm[comm$Order=="D_3__Betaproteobacteriales",1:48]
#comm <- comm[(comm$Class == "D_2__Gammaproteobacteria")&(comm$Order!= "D_3__Betaproteobacteriales"),1:48]
### Fungi
#comm <- comm[comm$Phylum=="p__Mortierellomycota",1:48]
### Protists
#comm <- comm[comm$Class=="Ochrophyta",1:48]
#comm <- comm[comm$Function=="Phototroph",1:48]

### After column 49, Listed as annotated information for species
comm = comm[,1:48]
treat = read.csv(treat.file, row.names = 1, header = T, sep = ",")
comm = as.data.frame(t(comm))

dim(comm)
comm = comm[, colSums(comm) > 0]
dim(comm)

row.names(comm) <- row.names(treat)
dim(comm)
comm = comm[rowSums(comm) > 0,]
dim(comm)

### Align OTU tables, treatment information, and phylogenetic tree files
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

### For a specific year, you can use the following code ###
### Taking the year 2009 as an example:
dist.soren = vegdist(comm[1:8,], binary = TRUE) ### We have a total of 8 samples per year
treat = treat[1:8,] ### We have a total of 8 samples per year
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
    dist.soren1 = vegdist(comm[(8*i-7):(8*i),], binary = TRUE) ### We have a total of 8 samples per year
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
write.csv(disper.soren,"Dispersion-Sorensen-Bacteria.csv")

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
write.csv(disper.uw,"Dispersion-Unifrac-Bacteria.csv")

### Temporal succession
### linear mixed-effects models + permutation tests
treatused<-read.csv("F:\\1-drought\\All analyses\\UNOISE3-0914-New\\GitHub\\Dataset\\Treatment info-Microbial biomass.csv",header = T,row.names = 1)
divindex<-read.csv("F:\\1-drought\\All analyses\\UNOISE3-0914-New\\GitHub\\Dataset\\Dispersion-Sorensen.csv",header = T,row.names = 1)
treatused$year<-treatused$year-2009  

### match names
divindex<-divindex[match(row.names(treatused),row.names(divindex)),] 
### check if names are matched
sum(row.names(divindex)==row.names(treatused)) 

output1 <- c()
output.slope <- c()
for(i in 1:ncol(divindex))
{
  div<-data.frame(disper=divindex[,i],treatused)
  
  ### disper.lmm: function used to calculate the dispersion trend over time
  disper.lmm <- function(dispersion,treat)
  {
    drought=treat$Treatment
    drought.lev=unique(drought)
    out=list()
    for(j in 1:length(drought.lev))
    {
      idj=which(drought==drought.lev[j])
      dispersion1 = dispersion[idj,]
      lmij=lmer(disper~year+((1+year)|block/plot),data=dispersion1)
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
  
  disperi =disper.lmm(dispersion = div,treat = treatused)
  slope.obs = as.vector(disperi[1,])
  r2.obs=as.vector(disperi[3:4,])
  aic.obs=as.vector(disperi[5:6,])
  ds.obs=(disperi[1,1])-(disperi[1,2])
  
  ### Permutation tests
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
                    dispersionr=div
                    dispersionr[,"year"]=year.lev[idi[match(div$year,year.lev)]]
                    disperc=disper.lmm(dispersion = dispersionr,treat = treatused)
                    out$slop.fix=as.vector(disperc[1,])
                    out$r2=as.vector(disperc[3:4,])
                    out$aic=as.vector(disperc[5:6,])
                    out$ds=(disperc[1,1])-(disperc[1,2])
                    out
                  })
  slop.ran =sapply(1:length(ind.rand),function(k){ind.rand[[k]]$slop.fix})
  r2.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$r2})
  aic.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$aic})
  ds.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$ds})
  p.ran = rbind(slop.ran,r2.ran,aic.ran,ds.ran) 

  EPS <- sqrt(.Machine$double.eps)
  slop.ran.1 = p.ran[1,]
  slop.ran.2 = p.ran[2,]
  r2.ran = p.ran[3:6,]
  aic.ran = p.ran[7:10,]
  ds.ran = p.ran[11,]
  
  p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS),na.rm=TRUE)+1)/(ncol(r2.ran)+1)
  p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS),na.rm=TRUE)+1)/(ncol(aic.ran)+1)
  p.values=rbind(matrix(p.r2,2,2),matrix(p.aic,2,2))
  rownames(p.values)=c("P.R2M","P.R2C","P.AIC1","P.AIC2")
  output=rbind(disperi,p.values)
  
  ### Permutation tests (Bootstrapping)
  m1 <- div[div$Treatment=="Drought",]
  m1 <- arrange(m1,desc(block))
  m2 <- div[div$Treatment=="Control",]
  m2 <- arrange(m2,desc(block))
  
  boot.slope <- function(dat,indices) {
    slope <- as.numeric(summary(lmer(disper~year+((1+year)|block/plot),data=dat[indices,]))$coefficients[2,1])
    slope
  }
  set.seed(1028)
  boot.slope.drought <- boot(data=m1,statistic=boot.slope,R=100)
  boot.slope.control <- boot(data=m2,statistic=boot.slope,R=100)
  slope <- as.data.frame(cbind(boot.slope.drought$t,boot.slope.control$t)) 
  colnames(slope)[1] <- paste("drought.slope",colnames(divindex)[i],sep="_")
  colnames(slope)[2] <- paste("control.slope",colnames(divindex)[i],sep="_")
  output.slope <- rbind(output.slope,as.data.frame(t(slope)))
  rowMeans(output.slope)
  result <- wilcox.test(as.numeric(slope[,1]),as.numeric(slope[,2]))
  slope.signif <- result$p.value
  
  ### Results output
  output=as.data.frame(t(output))
  output$p.dis.slope = slope.signif
  rownames(output)=c((paste0(colnames(divindex)[i],".drought")),(paste0(colnames(divindex)[i],".control")))
  output1 <- rbind(output1,output)
}

output1 <- as.data.frame(t(output1))
write.csv(output1,"AllbetaTaxa-LMM-Permutation-Dispersion-Sorensen.csv")

### Time-decay relationships (TDRs) -------------------------------------------------
### Taking the calculation based on Sorensen distances as an example
betai <- as.matrix(dist.soren)

### tdc.lmm: function used to calculate the TDRs
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
    blockj=treat$block[match(betaij3[,1],rownames(treat))]
    idj.use=which(plotj1==plotj2)
    cij.use=1-betaij3[idj.use,3]
    dtj.use=dtj[idj.use]
    plotj.use=plotj1[idj.use]
    blockj.use=blockj[idj.use]
    logCij=log(cij.use)
    logCij[logCij==-Inf]<-NA
    logdtj=log(dtj.use)
    logdtj[logdtj==-Inf]<-NA
    lmij=lmer(logCij~logdtj+((1+logdtj)|blockj.use/plotj.use))
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

### Permutation test
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
slop.ran.1 = p.ran[1,]
slop.ran.2 = p.ran[2,]
r2.ran = p.ran[3:6,]
aic.ran = p.ran[7:10,]
ds.ran = p.ran[11,]

p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS),na.rm=TRUE)+1)/(ncol(r2.ran)+1)
p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS),na.rm=TRUE)+1)/(ncol(aic.ran)+1)
p.values=rbind(matrix(p.r2,2,2),matrix(p.aic,2,2))
rownames(p.values)=c("P.R2M","P.R2C","P.AIC1","P.AIC2")
output=rbind(tdci,p.values)

### Permutation tests (Bootstrapping)
betaij3=dist.3col(betai)
betaij3$dtj=abs(treat$year[match(betaij3[,1],rownames(treat))]-treat$year[match(betaij3[,2],rownames(treat))])
betaij3$plotj1=treat$plot[match(betaij3[,1],rownames(treat))]
betaij3$plotj2=treat$plot[match(betaij3[,2],rownames(treat))]
betaij3$blockj=treat$block[match(betaij3[,1],rownames(treat))]
idj.use=which(betaij3$plotj1==betaij3$plotj2)
betaij3 = betaij3[idj.use,]
betaij3$Treatment=treat$Treatment[match(betaij3[,1],rownames(treat))]
betaij3$logCij = log(1-betaij3[,3])
betaij3[betaij3$logCij==-Inf]$logCij<-NA
betaij3$logdtj = log(betaij3[,4])

m1 <- betaij3[betaij3$Treatment=="Drought",]
m1 <- arrange(m1,desc(blockj))
m2 <- betaij3[betaij3$Treatment=="Control",]
m2 <- arrange(m2,desc(blockj))

boot.slope <- function(dat,indices) {
  slope <- as.numeric(summary(lmer(logCij~logdtj+((1+logdtj)|blockj/plotj1),data=dat[indices,]))$coefficients[2,1])
  slope
}
set.seed(1028)
boot.slope.drought <- boot(data=m1,statistic=boot.slope,R=100)
boot.slope.control <- boot(data=m2,statistic=boot.slope,R=100)
slope <- as.data.frame(cbind(boot.slope.drought$t,boot.slope.control$t)) 
output.slope <- c()
output.slope <- rbind(output.slope,as.data.frame(t(slope)))
rowMeans(output.slope)
result <- wilcox.test(as.numeric(slope[,1]),as.numeric(slope[,2]))
slope.signif <- result$p.value

### Results output
output=as.data.frame(t(tdci))
output$p.dis.slope = slope.signif
prefixi = "Bacteria"
rownames(output)=c((paste0(prefixi,".Sorensen.drought")),(paste0(prefixi,".Sorensen.control")))
#output1 <- c()
output1 <- rbind(output1,output)

output1 <- as.data.frame(t(output1))
write.csv(output1,"AllbetaTaxa-LMM-Permutation-TDR-Sorensen.csv")

### Temporal change in community differences (dissimilarity) ------------------------------------------
### Taking the calculation based on Sorensen distances as an example
betai <- as.matrix(dist.soren)

### dv.lmm: function used to calculate the dissimilarity trend over time
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

### Permutation test
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

### Results output
EPS <- sqrt(.Machine$double.eps)
p.slope=(rowSums(slop.ran>=(matrix(slope.obs,nr=nrow(slop.ran),nc=ncol(slop.ran))-EPS))+1)/(ncol(slop.ran)+1)
p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS))+1)/(ncol(r2.ran)+1)
p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS))+1)/(ncol(aic.ran)+1)
p.values=rbind(matrix(p.slope,1,1),matrix(p.r2,2,1),matrix(p.aic,2,1))
rownames(p.values)=c("P.slope","P.R2M","P.R2C","P.AIC1","P.AIC2")
output=rbind(dvi,p.values)
write.csv(output,"AllbetaTaxa-LMM-Permutation-Dissimilarity-Sorensen.csv")
