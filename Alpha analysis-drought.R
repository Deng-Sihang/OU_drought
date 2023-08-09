rm(list=ls())
library(lme4)
library(car)
library(plyr)
library(vegan)
library(dplyr)
library(picante)
library(ieggr)
save.wd <- iwd(choose.dir())

### Calculation of alpha diversity indices ---------------------------------------
### For bacterial communities as an example
### For phyla/classes/guilds, is the same
com.file = "Bacteria-zotu-resampled.csv"
tree.file = "Bacteria-tree.nwk"

comm = read.csv(com.file, row.names = 1, header = T, sep = ",")
comm = comm[,1:48]
dim(comm)
comm <- comm[, colSums(comm) > 0]
dim(comm)
comm = comm[rowSums(comm) > 0,]
dim(comm)

comm <- as.data.frame(t(comm))
alpha.tax <- alpha.g(comm, td.method = c("richness","shannon","invsimpson"))
tree <- lazyopen(tree.file)
is.rooted(tree)
tree = root(tree,1,r=T) 
is.rooted(tree)
spc <- match.name(cn.list = list(comm = comm), tree.list = list(tree = tree))
comm <- spc$comm
tree <- spc$tree
Faith.PD <- picante::pd(comm, tree, include.root = FALSE)
result <- cbind(alpha.tax,Faith.PD[,1])
result <- as.data.frame(result)
colnames(result) <- c("richness","shannon","invsimpson","PD")
write.csv(result,"Bacteria-alpha.csv")

### Temporal succession of soil microbial alpha-diversity ------------------------------------------------
### linear mixed-effects models + permutation tests
### It can also be used to calculate the temporal successions of biomass, relative abundance, and network topological properties
treatused<-read.csv("Treatment info-Microbial biomass.csv",header = T,row.names = 1)
divindex<-read.csv("Bacteria-alpha.csv",header = T,row.names = 1)
treatused$year<-treatused$year-2009  
divindex<-divindex[match(row.names(treatused),row.names(divindex)),] #match names
sum(row.names(divindex)==row.names(treatused)) #check if names are mached

divindex<-scale(divindex) ##recale the diversities

for(i in 1:ncol(divindex))
{
  div<-data.frame(divtest=divindex[,i],treatused)
  alpha.lmm <- function(alphai,treat)
  {
    drought=treat$Treatment
    drought.lev=unique(drought)
    out=list()
    for(j in 1:length(drought.lev))
    {
      idj=which(drought==drought.lev[j])
      alphai1 = alphai[idj,]
      lmij=lmer(divtest~year+((1+year)|plot),data=alphai1)
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
  
  alphai=alpha.lmm(alphai=div,treat = treatused)
  slope.obs = as.vector(alphai[1,])
  r2.obs=as.vector(alphai[3:4,])
  aic.obs=as.vector(alphai[5:6,])
  ds.obs=(-alphai[1,1])-(-alphai[1,2])
  
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
                    alphar=alpha.lmm(alphai=div1,treat = treatused)
                    out$slop.fix=as.vector(alphar[1,])
                    out$r2=as.vector(alphar[3:4,])
                    out$aic=as.vector(alphar[5:6,])
                    out$ds=(-alphar[1,1])-(-alphar[1,2])
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
  
  output=rbind(alphai,p.values)
  if (i == 1)
  {
    result = NULL
    result <- cbind(result, output)
  }else{
    result <- cbind(result, output)
  }
  result <- as.data.frame(result)
  result
}
write.csv(result,"Bacteria-alpha-LMM-Permutation.csv")


### Yearly impact of drought on microbial alpha diversity -------------------
### Effect size
treatused<-read.csv("Treatment info-Microbial biomass.csv",header = T,row.names = 1)
divindex<-read.csv("Fungi-alpha.csv",header = T,row.names = 1)
treatused$year<-treatused$year-2009  
divindex<-divindex[match(row.names(treatused),row.names(divindex)),] #match names
sum(row.names(divindex)==row.names(treatused)) #check if names are mached

### Using the following function, you can directly output the results of the effect size for each year
for(i in 1:6) ###6 represents the 6-year prolonged drought.
{
  dat = divindex[(8*i-7):(8*i),]
  dat <-dat[match(row.names(treatused),row.names(dat)),] #match names
  sum(row.names(dat)==row.names(treatused))
  dat <-scale(dat) 
  
  divs1<-sapply(1:ncol(dat),function(j){
    message("Now j=",j," in ",ncol(dat),". ",date())
    if (length(unique(dat[,j]))<3){
      result<-rep(NA,8)
    } else {
      div<-data.frame(divtest=dat[,j],treatused)
      
      fm1<-lmer(divtest~Drought+(1|block),data=div)
      
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
  colnames(divs1)<-colnames(dat)
  tag<-as.character(rownames(divindex)[8*i])
  write.csv(divs1,paste0(tag,".csv"))
}

### Identification of environmental drivers influencing the microbial biomass---------------------------------------------------
### Random forest
library(randomForest)
library(rfPermute)
library(rfUtilities)
dat <- read.table("Random forest-biomass.csv",sep=",",check.names = FALSE,row.names =1,header = TRUE)
### Taking total microbial biomass as an example
dat <- dat[,1:7] 
dat <- as.data.frame(scale(dat))

set.seed(1028)
biomass.forest <- randomForest(total~.,data = dat,importance = TRUE,proximity=TRUE)
### Calculate the overall significance and explanatory power (r-values) of the model
biomass.pvalue <- rf.significance(biomass.forest,dat[,-1],nperm = 1000,ntree=500)
biomass.pvalue
### Calculate the IncMSE(%) and significance of each variable
biomass.rfP <- rfPermute(total~.,data = dat,importance = TRUE,ntree = 500,nrep = 1000,num.cores =1,proximity=TRUE)
index.importance.rfp <- data.frame(importance(biomass.rfP,scale=TRUE),check.names = FALSE)