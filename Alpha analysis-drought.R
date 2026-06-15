rm(list=ls())
library(lme4)
library(car)
library(plyr)
library(vegan)
library(dplyr)
library(picante)
library(ieggr)
library(boot)
library(nlme)
library(tibble)
library(effsize)
library(reshape2)
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
#comm <- comm[comm$Strategy=="Oligotroph",1:48]
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
### It can also be used to calculate the temporal successions of microbial biomass, relative abundance, community dispersion, and network topological properties.
treatused<-read.csv("Treatment info-Microbial biomass.csv",header = T,row.names = 1)
divindex<-read.csv("Bacteria-Fungi-Protists-alphadiversity.csv",header = T,row.names = 1)
treatused$year<-treatused$year-2009  

### match names
divindex<-divindex[match(row.names(treatused),row.names(divindex)),] 
### check if names are matched
sum(row.names(divindex)==row.names(treatused)) 
### rescale the diversities
divindex <- divindex %>% mutate(across(everything(), as.numeric))
divindex<-scale(divindex) 

output1 <- c()
for(i in 1:ncol(divindex))
{
  div<-data.frame(divtest=divindex[,i],treatused)
  div <- na.omit(div)
  
  result1 <- c()
  for (treat in unique(div$Treatment)){
    result <- try(lme(divtest~year, data = div[div$Treatment==treat,],
                      random = ~ 1|block/plot,
                      correlation = corAR1(form = ~ year|block/plot)), silent = TRUE)
    
    if (inherits(result, "try-error")) {
      lmij=lme(divtest~year, data = div[div$Treatment==treat,],
                        random = ~ 1|block,
                        correlation = corAR1(form = ~ year|block))
      } else {
        lmij=lme(divtest~year, data = div[div$Treatment==treat,],
                 random = ~ 1|block/plot,
                 correlation = corAR1(form = ~ year|block/plot))
      }
    
    lmijsm=summary(lmij)
    AIC1=AIC(lmij)
    result <- try(rsquared.glmm(lmij), silent = TRUE)
    
    if (inherits(result, "try-error")) {
      r2ij=performance::r2(lmij)
      r2ij=data.frame(Marginal=r2ij$R2_marginal,Conditional=r2ij$R2_conditional,AIC=AIC1)
    } else {
      r2ij=rsquared.glmm(lmij)
    }
    lmijCS=Anova(lmij,type = "II")
    result=data.frame(Class=colnames(divindex)[i],Treatment = treat,slop.fix=lmijsm$coefficients$fixed[2],slop.sd=coef(summary(lmij))[2,2],
               R2M=r2ij$Marginal,R2C=r2ij$Conditional,AIC1=AIC1,AIC2=r2ij$AIC,
               P.typeII=lmijCS[[3]],Chisq=lmijCS[[1]])
    result1 <- rbind(result1,result)
    row.names(result1) <- NULL
  }
  
  result1=as.data.frame(t(result1))
  slope.obs = as.vector(as.numeric(result1[3,]))
  r2.obs=as.vector(as.numeric(result1[6,]))
  aic.obs=as.vector(as.numeric(result1[8,]))
  ds.obs=(-as.numeric(result1[3,1]))-(-as.numeric(result1[3,2]))
  result1=as.data.frame(t(result1))
  row.names(result1) <- result1$Treatment
  
  ### Interaction term
  result <- try(lme(divtest~year*Drought, data = div,
                    random = ~ 1|block/plot,
                    correlation = corAR1(form = ~ year|block/plot)), silent = TRUE)
  
  if (inherits(result, "try-error")) {
    lmij=lme(divtest~year*Drought, data = div,
             random = ~ 1|block,
             correlation = corAR1(form = ~ year|block))
  } else {
    lmij=lme(divtest~year*Drought, data = div,
             random = ~ 1|block/plot,
             correlation = corAR1(form = ~ year|block/plot))
  }
  
  lmijsm=summary(lmij)
  beta.year = lmijsm$coefficients$fixed[2]
  beta.drought = lmijsm$coefficients$fixed[3]
  betai = lmijsm$coefficients$fixed[4]
  result1$beta = betai
  
  ### Permutation tests
  year.lev=unique(div$year)
  rand =1000
  year.perm=vegan:::getPermuteMatrix(rand,length(year.lev))
  trace.seq=seq(from=1,to=rand,by = 100)
  out <- c()
  betar <- c()
  for (j in 1:nrow(year.perm)){
   idi=year.perm[j,]
   div1=div
   div1[,"year"]=year.lev[idi[match(div$year,year.lev)]]
     
    resultp1 <- c()
    for (treat in unique(div1$Treatment)){
      result <- try(lme(divtest~year, data = div1[div1$Treatment==treat,],
                       random = ~ 1|block/plot,
                        correlation = corAR1(form = ~ year|block/plot)), silent = TRUE)
   
     if (inherits(result, "try-error")) {
       result <- try(lme(divtest~year, data = div1[div1$Treatment==treat,],
                         random = ~ 1|block,
                          correlation = corAR1(form = ~ year|block)), silent = TRUE)
       
        if (inherits(result, "try-error")) {
          next
        } else {
         lmij=lme(divtest~year, data = div1[div1$Treatment==treat,],
                 random = ~ 1|block,
                 correlation = corAR1(form = ~ year|block))
        }
      } else {
        lmij=lme(divtest~year, data = div1[div1$Treatment==treat,],
                 random = ~ 1|block/plot,
                 correlation = corAR1(form = ~ year|block/plot))
      }
      
      lmijsm=summary(lmij)
     AIC1=AIC(lmij)
      
      result <- try(rsquared.glmm(lmij), silent = TRUE)
      if (inherits(result, "try-error")) {
        next
      } else {
        r2ij=rsquared.glmm(lmij)
      }
      
      lmijCS=Anova(lmij,type = "II")
      result=data.frame(Class=colnames(divindex)[i],Treatment = treat,slop.fix=lmijsm$coefficients$fixed[2],slop.sd=coef(summary(lmij))[2,2],
                        R2M=r2ij$Marginal,R2C=r2ij$Conditional,AIC1=AIC1,AIC2=r2ij$AIC,
                        P.typeII=lmijCS[[3]],Chisq=lmijCS[[1]])
      resultp1 <- rbind(resultp1,result)
      row.names(resultp1) <- NULL
    }
    
    if (nrow(resultp1)!=2){
      next
    }else{
      resultp1=as.data.frame(t(resultp1[,3]))
     colnames(resultp1) <- c("Drought","Control")
      resultp1$ds.r <- (-resultp1$Drought)-(-resultp1$Control)
      row.names(resultp1) <- NULL
      out = rbind(out,resultp1)
    }
     
     ## Interaction term
    result <- try(lme(divtest~year*Drought, data = div1,
                      random = ~ 1|block/plot,
                      correlation = corAR1(form = ~ year|block/plot)), silent = TRUE)
    
    if (inherits(result, "try-error")) {
      result <- try(lme(divtest~year*Drought, data = div1,
                        random = ~ 1|block,
                        correlation = corAR1(form = ~ year|block)), silent = TRUE)
      if (inherits(result, "try-error")) {
        next
      }else{
        lmij=lme(divtest~year*Drought, data = div1,
                 random = ~ 1|block,
                 correlation = corAR1(form = ~ year|block))
      }
    } else {
      lmij=lme(divtest~year*Drought, data = div1,
               random = ~ 1|block/plot,
              correlation = corAR1(form = ~ year|block/plot))
    }
    
    lmijsm=summary(lmij)
    result <- data.frame(beta.year.r = lmijsm$coefficients$fixed[2],
                         beta.drought.r = lmijsm$coefficients$fixed[3],
                         beta = lmijsm$coefficients$fixed[4])
    row.names(result) <- NULL
    betar = rbind(betar,result)
  }
    
    if (result1["Drought",]$slop.fix<0){
     p.slope.perm.drought <- nrow(subset(out,(Drought < as.numeric(result1["Drought",]$slop.fix))))/(nrow(out)+1)
    }else {p.slope.perm.drought <- nrow(subset(out,(Drought > as.numeric(result1["Drought",]$slop.fix))))/(nrow(out)+1)}
    
    if (result1["Control",]$slop.fix<0){
     p.slope.perm.control <- nrow(subset(out,(Control < as.numeric(result1["Control",]$slop.fix))))/(nrow(out)+1)
    }else {p.slope.perm.control <- nrow(subset(out,(Control > as.numeric(result1["Control",]$slop.fix))))/(nrow(out)+1)}
    result1$p.slope.perm = c(p.slope.perm.drought,p.slope.perm.control)
    
    p.beta.perm <- sum((if(beta.year < 0) betar$beta.year.r < 0 else betar$beta.year.r > 0) & 
                         (if(beta.drought > 0) betar$beta.drought.r > 0 else betar$beta.drought.r < 0) & 
                         (if(betai < 0) betar$beta < betai else betar$beta > betai), na.rm = TRUE)/(nrow(out)+1)
    result1$p.beta.perm = p.beta.perm
    
    ## Results output
    row.names(result1) <- NULL
    output1 <- rbind(output1,result1)
}
