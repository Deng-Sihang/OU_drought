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

### Calculation of microbial community dispersion --------------------------------------------------------------
### Sorensen distances
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

### Time-decay relationships (TDRs) -------------------------------------------------
### Taking the calculation based on Sorensen distances as an example
betai <- as.matrix(dist.soren)

output1 <- c()
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
tdci=as.data.frame(t(tdci))

### Permutation tests
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
p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS),na.rm=TRUE)+1)/(ncol(r2.ran)+1)
p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS),na.rm=TRUE)+1)/(ncol(aic.ran)+1)
p.values=rbind(matrix(p.r2,2,2),matrix(p.aic,2,2))
rownames(p.values)=c("P.R2M","P.R2C","P.AIC1","P.AIC2")
tdci$p.slope.perm = c(p.values[4,1],p.values[4,2])

### Interaction term
betai <- as.matrix(dist.soren)
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

lmij=lmer(logCij~logdtj*Treatment+((1+logdtj)|blockj/plotj1),data=betaij3)
lmijsm=summary(lmij)
lmijCS=Anova(lmij,type = "II")
beta.year = lmijsm$coefficients[2,1]
beta.drought = lmijsm$coefficients[3,1]
betai = lmijsm$coefficients[4,1]
tdci$beta = betai
tdci$p.beta = lmijCS[3,3]

### Permutation tests
year.lev=unique(betaij3$dtj)
rand =1000
year.perm=vegan:::getPermuteMatrix(rand,length(year.lev))
trace.seq=seq(from=1,to=rand,by = 100)
betar <- c()
for (j in 1:nrow(year.perm)){
  idi=year.perm[j,]
  betaij31=betaij3
  betaij31[,"dtj"]=year.lev[idi[match(betaij3$dtj,year.lev)]]
  betaij31$logdtj <- log(betaij31$dtj)

  lmij=lmer(logCij~logdtj*Treatment+((1+logdtj)|blockj/plotj1),data=betaij31)
  lmijsm=summary(lmij)
  result <- data.frame(beta.year.r = lmijsm$coefficients[2,1],
                       beta.drought.r = lmijsm$coefficients[3,1],
                       beta = lmijsm$coefficients[4,1])
  betar = rbind(betar,result)
}

p.ds.perm <- sum((if(tdci["Drought",]$slop.fix < 0) out$Drought < 0 else out$Drought > 0) & 
                   (if(tdci["Control",]$slop.fix < 0) out$Control < 0 else out$Control > 0) & 
                   (if(ds.obs > 0) out$ds.r > ds.obs else out$ds.r < ds.obs), na.rm = TRUE)/(nrow(out)+1)

p.beta.perm <- sum((if(beta.year < 0) betar$beta.year.r < 0 else betar$beta.year.r > 0) & 
                     (if(beta.drought > 0) betar$beta.drought.r > 0 else betar$beta.drought.r < 0) & 
                     (if(betai < 0) betar$beta < betai else betar$beta > betai), na.rm = TRUE)/(nrow(out)+1)
tdci$p.ds.perm = p.ds.perm
tdci$p.beta.perm = p.beta.perm

### Results output
prefixi = "Bacteria.Sorensen"
output=data.frame(type=prefixi,treatment = row.names(tdci),tdci)
row.names(output) <- NULL
output1 <- rbind(output1,output)
write.csv(output1,"AllbetaTaxa-LMM-TDR-Sorensen.csv")
