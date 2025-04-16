rm(list=ls())
library(lme4)
library(car)
library(plyr)
library(vegan)
library(dplyr)
library(picante)
library(ieggr)
library(boot)
save.wd <- iwd(choose.dir())

### Calculation of alpha diversity indices ---------------------------------------
### Taking bacterial communities as an example for calculation
com.file = "Bacteria-zotu-resampled.csv"
tree.file = "Bacteria-tree.nwk"
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
#comm <- comm[comm$Function=="Consumer",1:48]

### After column 49, Listed as annotated information for species
comm = comm[,1:48]
dim(comm)
comm <- comm[, colSums(comm) > 0]
dim(comm)
comm = comm[rowSums(comm) > 0,]
dim(comm)

### Calculate taxonomic diversity indices 
comm <- as.data.frame(t(comm))
alpha.tax <- alpha.g(comm, td.method = c("richness","shannon","invsimpson"))

### Import the phylogenetic tree file and align it with the OTU table
tree <- lazyopen(tree.file)
is.rooted(tree)
tree = root(tree,1,r=T) 
is.rooted(tree)
spc <- match.name(cn.list = list(comm = comm), tree.list = list(tree = tree))
comm <- spc$comm
tree <- spc$tree

### Calculate phylogenetic diversity index 
Faith.PD <- picante::pd(comm, tree, include.root = FALSE)
result <- cbind(alpha.tax,Faith.PD[,1])
result <- as.data.frame(result)
colnames(result) <- c("richness","shannon","invsimpson","PD")
write.csv(result,"Bacteria-alpha.csv")

### Temporal succession of soil microbial alpha-diversity ------------------------------------------------
### linear mixed-effects models + permutation tests (Bootstrapping)
### It can also be used to calculate the temporal successions of microbial biomass, relative abundance, and network topological properties.
treatused<-read.csv("Treatment info-Microbial biomass.csv",header = T,row.names = 1)
divindex<-read.csv("Bacteria-Fungi-Protists-alphadiversity.csv",header = T,row.names = 1)
treatused$year<-treatused$year-2009  

### match names
divindex<-divindex[match(row.names(treatused),row.names(divindex)),] 
### check if names are matched
sum(row.names(divindex)==row.names(treatused)) 
### rescale the diversities
divindex<-scale(divindex) 

output1 <- c()
output.slope <- c()
for(i in 1:ncol(divindex))
{
  div<-data.frame(divtest=divindex[,i],treatused)
  
  ### alpha.lmm: function used to calculate the diversity trend over time
  alpha.lmm <- function(alphai,treat)
  {
    drought=treat$Treatment
    drought.lev=unique(drought)
    out=list()
    for(j in 1:length(drought.lev))
    {
      idj=which(drought==drought.lev[j])
      alphai1 = alphai[idj,]
      lmij=lmer(divtest~year+((1+year)|block/plot),data=alphai1)
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
  
  ### Permutation tests (Bootstrapping)
  m1 <- div[div$Treatment=="Drought",]
  m1 <- arrange(m1,desc(block))
  m2 <- div[div$Treatment=="Control",]
  m2 <- arrange(m2,desc(block))

  boot.slope <- function(dat,indices) {
    slope <- as.numeric(summary(lmer(divtest~year+((1+year)|block/plot),data=dat[indices,]))$coefficients[2,1])
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
  output=as.data.frame(t(alphai))
  output$p.dis.slope = slope.signif
  rownames(output)=c((paste0(colnames(divindex)[i],".drought")),(paste0(colnames(divindex)[i],".control")))
  output1 <- rbind(output1,output)
}

output1 <- as.data.frame(t(output1))
write.csv(output1,"Allnetwork-LMM-Permutation.csv")

### Yearly impact of drought on fungal alpha diversity -------------------
### Effect size
treatused<-read.csv("Treatment info-Microbial biomass.csv",header = T,row.names = 1)
divindex<-read.csv("Bacteria-Fungi-Protists-alphadiversity.csv",header = T,row.names = 1)
treatused$year<-treatused$year-2009  

### match names
divindex<-divindex[match(row.names(treatused),row.names(divindex)),] 
### check if names are matched
sum(row.names(divindex)==row.names(treatused)) 

### Using the following function, you can directly output the results of the effect size for each year
for(i in 1:6) ### 6 represents the 6-year prolonged drought.
{
  dat = divindex[(8*i-7):(8*i),] ### We have a total of 8 samples per year
  dat <-dat[match(row.names(treatused),row.names(dat)),] ### match names
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
      coefs<-coef(summary(fm1))[ , "Estimate"]  ### four coefs
      names(coefs)<-paste0(names(coefs),".mean")
      
      SEvalues<-coef(summary(fm1))[ , "Std. Error"] ### standard errors
      names(SEvalues)<-paste0(names(SEvalues),".se")
      
      tvalues<-coef(summary(fm1))[ , "t value"] ### t values
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

### Taking total microbial biomass as an example (1-7 total microbial biomass, 8-14 total bacterial biomass)
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
