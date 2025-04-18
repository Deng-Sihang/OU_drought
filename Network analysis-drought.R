rm(list=ls())
library(lme4)
library(car)
library(plyr)
library(vegan)
library(dplyr)
library(ieggr)
library(reshape2)
library(picante)
library(igraph)
library(sna)
save.wd <- iwd(choose.dir())

### Extracting subset networks from global MENs --------
### The matrix "assmc" represents the adjacency matrix, which is obtained as follows: 
### First, a correlation matrix is generated based on Pearson correlations of centered log-transformed OTU abundances. 
### The Random Matrix Theory (RMT) method is then applied to determine the optimal threshold value, and a cutoff value of 0.8 is selected. 
### The resulting adjacency matrix includes only the correlations with absolute coefficient values equal to or greater than 0.8.
### All the calculations mentioned above were conducted on the MENAP website (http://ieg4.rccc.ou.edu/MENA/)
### Note: Drought and control groups should be separately considered for network construction
assmc = read.csv("Bacteria-drought-corrletion0.8.csv",header = T,row.names = 1)

### The function "net.pro" can be used to calculate the main topological properties of a network
### including: Total nodes,Total links,Average degree (avgK),Average weighted degree (avgKw), Average clustering coefficient (avgCC)
### including: Average path distance (GD), Density (D),Transitivity (Trans),Krackhardt Connectedness (Con),Power-law distribution
net.pro = function (assmc)
{
  require(igraph)
  amg = abs(as.matrix(assmc))
  amg[is.na(amg)] = 0
  diag(amg) = 0
  ag = igraph::graph_from_adjacency_matrix(adjmatrix = amg, 
                                           mode = "undirected", weighted = TRUE, diag = FALSE)
  nod.name = V(ag)$name
  centr.deg = igraph::centr_degree(ag, mode = "all", loops = FALSE, 
                                   normalized = TRUE)
  nod.deg = centr.deg$res
  net.deg = mean(nod.deg)
  nod.wdeg = igraph::strength(ag, mode = "all", loops = FALSE)
  net.wdeg = mean(nod.wdeg)
  require(sna)
  nod.transit = igraph::transitivity(ag, type = "local")
  net.noden = nrow(assmc)
  net.edgen = sum(abs(amg[upper.tri(amg)]) > 0)
  net.meandis = igraph::mean_distance(ag, directed = FALSE)
  net.density = igraph::edge_density(ag, loops = FALSE)
  net.cc = mean(nod.transit, na.rm = TRUE)
  net.transit = igraph::transitivity(ag, type = "global")
  net.connect = sna::connectedness(amg)
  fitpowerlaw = function(graph) {
    d = igraph::degree(graph, mode = "all")
    dd = igraph::degree_distribution(graph, mode = "all", 
                                     loops = FALSE, cumulative = FALSE)
    degree = 1:max(d)
    probability = dd[-1]
    nonzero.position = which(probability != 0)
    probability = probability[nonzero.position]
    degree = degree[nonzero.position]
    reg = lm(log(probability) ~ log(degree))
    cozf = coef(reg)
    alpha = -cozf[[2]]
    R.square = summary(reg)$r.squared
    lmijCS=Anova(reg,type = "II")
    p = lmijCS$`Pr(>F)`[[1]]
    list(alpha = alpha, R.square = R.square,p = p)
  }
  net.plr2 = fitpowerlaw(ag)$R.square
  net.alpha = fitpowerlaw(ag)$alpha
  net.p = fitpowerlaw(ag)$p
  net.att = data.frame(NetworkIndex = c("Total nodes", "Total links", "Average degree (avgK)", "Average weighted degree (avgKw)", 
                                        "Average clustering coefficient (avgCC)", "Average path distance (GD)", 
                                        "Density (D)","Transitivity (Trans)", "Krackhardt Connectedness (Con)",
                                        "alpha","r","P"), 
                       Value = c(net.noden, net.edgen,net.deg, net.wdeg,  
                                 net.cc, net.meandis, net.density,net.transit, net.connect,
                                 net.alpha,net.plr2,net.p), 
                       stringsAsFactors = FALSE)
  net.att
}

### Taking the global MEN as an example
a = net.pro(assmc)

### The function "net.rand" can be used to generate random networks and compare them with empirical networks
net.rand = function (assmc)
{
  require(igraph)
  amg = abs(as.matrix(assmc))
  amg[is.na(amg)] = 0
  diag(amg) = 0
  ag = igraph::graph_from_adjacency_matrix(adjmatrix = amg, 
                                           mode = "undirected", weighted = TRUE, diag = FALSE)
  neta.obs = net.pro(assmc)
  modul.obs = module(assmc = assmc, absolute = TRUE, methods = c("greedy"))
  mod.obs = modul.obs$module
  mod.node = modul.obs$node
  
  netp.obs = rbind(neta.obs, data.frame(NetworkIndex = c(paste0("Module.number.", mod.obs$Method), paste0("Modularity.", mod.obs$Method)), 
                                        Value = c(mod.obs$Module.number, mod.obs$Modularity), 
                                        stringsAsFactors = FALSE))
  rand.method = c("swap")
  randone <- function(ag, rand.method, swap.iter = 100, 
                      endpoint.p = 0.5) {
    require(igraph)
    if (rand.method == "swap") {
      randm = igraph::rewire(ag, with = igraph::keeping_degseq(loops = FALSE, 
                                                               niter = igraph::vcount(ag) * 10))
      randm = igraph::set_edge_attr(graph = randm, 
                                    name = "weight", value = E(ag)$weight)
    }
    else if (rand.method == "endpoints") {
      randm = igraph::rewire(ag, with = igraph::each_edge(prob = endpoint.p, 
                                                          loops = FALSE, multiple = FALSE))
    }
    else {
      stop("rand.method is not valid.")
    }
    igraph::as_adj(randm, sparse = FALSE, names = TRUE, 
                   attr = "weight")
  }
  netp.rand <- sapply(1:100, function(i) {
    message("Now i=", i, " in ", 100, ". ", date())
    ramgi = randone(ag = ag, rand.method = rand.method, 
                    swap.iter = swap.iter, endpoint.p = endpoint.p)
    netai = net.pro(ramgi)
    moduli = module(assmc = ramgi, absolute = TRUE,methods = c("greedy"))
    modi = moduli$module
    as.numeric(c(netai[,2], modi$Module.number, 
                 modi$Modularity))
  })
  z.test = function(a, mu, two.tail = FALSE) {
    zeta = (mean(a, na.rm = TRUE) - mu)/sd(a, na.rm = TRUE)
    if (two.tail) {
      p = 2 * pnorm(-abs(zeta))
    }
    else {
      p = pnorm(-abs(zeta))
    }
    list(z = zeta, p = p)
  }
  random.mean = rowMeans(netp.rand, na.rm = TRUE)
  random.sd = apply(netp.rand, 1, sd, na.rm = TRUE)
  z.value = (random.mean - as.numeric(netp.obs$Value))/random.sd
  p.ztest = pnorm(-abs(z.value))
  EPS <- sqrt(.Machine$double.eps)
  p.lower = (rowSums(netp.rand >= (matrix(as.numeric(netp.obs$Value), 
                                          nrow = nrow(netp.rand), ncol = ncol(netp.rand)) - EPS), 
                     na.rm = TRUE) + 1)/(1 + rowSums(!is.na(netp.rand)))
  p.higher = (rowSums(netp.rand <= (matrix(as.numeric(netp.obs$Value), 
                                           nrow = nrow(netp.rand), ncol = ncol(netp.rand)) + EPS), 
                      na.rm = TRUE) + 1)/(1 + rowSums(!is.na(netp.rand)))
  p.count = apply(cbind(p.lower, p.higher), 1, min)
  out = data.frame(netp.obs, random.mean = random.mean, random.sd = random.sd, 
                   z.value = z.value, p.ztest = p.ztest, p.count = p.count,stringsAsFactors = FALSE)
  colnames(out) = c("NetworkIndex", "Empirical.Results", "Random.Mean", 
                    "Random.Stdev", "Z.value", "P.Ztest", "P.count")
  out
}

### The function “module” can be used to calculate network modularity-related parameters and identify key nodes
module = function (assmc, methods = c("greedy", "walk", "eigen", "betweenness","infomap", "MLopt"), absolute = TRUE, zi.threshold = 2.5, 
          Pi.threshold = 0.62) 
{
  require(igraph)
  if (absolute) {
    amg = abs(as.matrix(assmc))
  }
  else {
    amg = as.matrix(assmc)
  }
  amg[is.na(amg)] = 0
  diag(amg) = 0
  ag = igraph::graph_from_adjacency_matrix(adjmatrix = amg, 
                                           mode = "undirected", weighted = TRUE, diag = FALSE)
  mods = list()
  if ("greedy" %in% methods) {
    mods$Greedy.optimization = try(igraph::cluster_fast_greedy(graph = ag, 
                                                               merges = TRUE, modularity = TRUE, membership = TRUE, 
                                                               weights = E(ag)$weight))
  }
  if ("eigen" %in% methods) {
    mods$Leading.eigenvector = try(igraph::cluster_leading_eigen(ag))
  }
  if ("walk" %in% methods) {
    mods$Short.random.walk = try(igraph::cluster_walktrap(ag, 
                                                          weights = E(ag)$weight, steps = 4, merges = TRUE, 
                                                          modularity = TRUE, membership = TRUE))
  }
  if ("betweenness" %in% methods) {
    mods$Edge.betweenness = try(igraph::cluster_edge_betweenness(ag, 
                                                                 weights = E(ag)$weight, directed = FALSE, edge.betweenness = TRUE, 
                                                                 merges = TRUE, bridges = TRUE, modularity = TRUE, 
                                                                 membership = TRUE))
  }
  if ("infomap" %in% methods) {
    mods$infomap = try(igraph::cluster_infomap(ag, nb.trials = 10, 
                                               modularity = TRUE))
  }
  if ("MLopt" %in% methods) {
    mods$Multi.level.optimization = try(igraph::cluster_louvain(ag))
  }
  if (length(mods) == 0) 
    stop("All method names are not valid.")
  module.att = data.frame(Method = names(mods), Module.number = sapply(1:length(mods), 
                                                                       function(i) {
                                                                         if (inherits(mods[[i]], "try-error")) {
                                                                           out = NA
                                                                         }
                                                                         else {
                                                                           out = length(mods[[i]])
                                                                         }
                                                                         out
                                                                       }), Modularity = sapply(1:length(mods), function(i) {
                                                                         if (inherits(mods[[i]], "try-error")) {
                                                                           out = NA
                                                                         }
                                                                         else {
                                                                           out = modularity(mods[[i]])
                                                                         }
                                                                         out
                                                                       }), stringsAsFactors = FALSE)
  rownames(module.att) = c()
  zpvalue <- function(mod, amg, method.name) {
    if (inherits(mod, "try-error")) {
      out = data.frame(Method = method.name, ID = rownames(amg), 
                       ModuleID = NA, zi = NA, Pi = NA, Classification = NA, 
                       stringsAsFactors = FALSE)
    }
    else {
      nod.modid = igraph::membership(mod)
      modlen = length(mod)
      zi <- kis2 <- rep(0, length(nod.modid))
      ki = rowSums(amg > 0)
      for (i in 1:modlen) {
        idi = which(nod.modid == i)
        nodi = names(nod.modid)[idi]
        idii = match(nodi, rownames(amg))
        kib = rowSums(amg[idii, idii, drop = FALSE] > 
                        0)
        zib = (kib - mean(kib))/sd(kib)
        zib[which(is.na(zib))] = 0
        zi[idi] = zib
        kisi = rowSums(amg[, idii, drop = FALSE] > 0)
        kis2 = kis2 + (kisi^2)
      }
      Pi = 1 - (kis2/(ki^2))
      nod.class = rep("Peripheral", length(nod.modid))
      nod.class[which(zi <= zi.threshold & Pi > Pi.threshold)] = "Connectors"
      nod.class[which(zi > zi.threshold & Pi <= Pi.threshold)] = "Module.hubs"
      nod.class[which(zi > zi.threshold & Pi > Pi.threshold)] = "Network.hubs"
      out = data.frame(Method = method.name, ID = names(nod.modid), 
                       ModuleID = as.vector(nod.modid), zi = zi, Pi = Pi, 
                       Classification = nod.class, stringsAsFactors = FALSE)
    }
    rownames(out) = c()
    out
  }
  node.mod.att = Reduce(rbind, lapply(1:length(mods), function(i) {
    zpvalue(mod = mods[[i]], amg = amg, method.name = names(mods)[i])
  }))
  list(module = module.att, node = node.mod.att)
}

### Taking the global MEN as an example to calculate network properties and identify key nodes
result = net.rand(assmc) 
modul.obs = module(assmc = assmc, absolute = TRUE, methods = c("greedy"))
mod.node = modul.obs$node

### Extracting: Compute the topological properties of each subset network
otu.file="Bacteria-drought-oturesampled.txt"
comm = lazyopen(otu.file)

for (i in 1:ncol(comm)){
  ### Extract the corresponding OTU table associated with the subset network
  comms = comm[comm[,i]!=0,]
  match.rowname=match.name(both.list = list(assmc=assmc),rn.list = list(comms=comms))
  assmc1 =  match.rowname$assmc
  assmc2 <- assmc1[rowSums(assmc1,na.rm = TRUE)!=1,colSums(assmc1,na.rm = TRUE)!=1]
  
  ### Compute the properties of the subset network
  result = net.pro(assmc2)
  mod = module(assmc = assmc2, absolute = TRUE, methods = c("greedy"))
  modu = mod$module
  
  ### Results output
  result.all = rbind(result, data.frame(NetworkIndex = c(paste0("Module.number.", modu$Method), paste0("Modularity.", modu$Method)), 
                                        Value = c(modu$Module.number, modu$Modularity), 
                                        stringsAsFactors = FALSE))
  if (i ==1){
    output = result.all
  }else{
    output = cbind(output,result.all[,2])
  }
}

### Calculate the changes in properties of key nodes -------------------
treatused<-read.csv("Treatment info-Microbial biomass.csv",header = T,row.names = 1)
divindex<-read.csv("keystone nodes-properties.csv",header = T,row.names = 1)
### It can also be applied for calculating changes in environmental factors

treatused$year<-treatused$year-2009 
### match names
divindex<-divindex[match(row.names(treatused),row.names(divindex)),] 
### check if names are matched
sum(row.names(divindex)==row.names(treatused)) 
### rescale the diversities
divindex<-scale(divindex) 

### Calculate the effect sizes of treatment on the properties of key nodes
divs1<-sapply(1:ncol(divindex),function(j){
  message("Now j=",j," in ",ncol(divindex),". ",date())
  if (length(unique(divindex[,j]))<3){
    result<-rep(NA,8)
  } else {
    div<-data.frame(divtest=divindex[,j],treatused)
    
    fm1<-lmer(divtest~Drought+(1|block)+(1|year),data=div)
    
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
colnames(divs1)<-colnames(divindex)

### Calculate the temporal succession of network dissimilarity -------------------
treatused<-read.csv("Treatment info-Microbial biomass.csv",header = T,row.names = 1)
divindex<-read.csv("Bacteria-Fungi-Protists-networkdissimlarity.csv",header = T,row.names = 1)
treatused$year<-treatused$year-2009  

### match names
divindex<-divindex[match(row.names(treatused),row.names(divindex)),] 
### check if names are matched
sum(row.names(divindex)==row.names(treatused)) 
### rescale the diversities
divindex<-scale(divindex) 

for(i in 1:ncol(divindex))
{
  div<-data.frame(divtest=divindex[,i],treatused)
  
  ### nd.lmm: function used to calculate the network dissimilarity trend over time
  nd.lmm <- function(ndi,treat)
  {
     lmij=lmer(divtest~year+((1+year)|block),data=ndi)
     lmijsm=summary(lmij)
     AIC1=AIC(lmij)
     r2ij=rsquared.glmm(lmij)
     lmijCS=Anova(lmij,type = "II")
     out=c(slop.fix=lmijsm$coefficients[2,1],slop.sd=coef(summary(lmij))[ , "Std. Error"][2],
                 R2M=r2ij$Marginal,R2C=r2ij$Conditional,AIC1=AIC1,AIC2=r2ij$AIC,
                 P.typeII=lmijCS[[3]],Chisq=lmijCS[[1]])
     outs=as.data.frame(out)
     colnames(outs)=colnames(divindex)[i]
     outs
  }
  
  ndi=nd.lmm(ndi=div,treat = treatused)
  slope.obs = as.vector(ndi[1,])
  r2.obs=as.vector(ndi[3:4,])
  aic.obs=as.vector(ndi[5:6,])
  
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
                    ndr=nd.lmm(ndi=div1,treat = treatused)
                    out$slop.fix=as.vector(ndr[1,])
                    out$r2=as.vector(ndr[3:4,])
                    out$aic=as.vector(ndr[5:6,])
                    out
                  })
  slop.ran =sapply(1:length(ind.rand),function(k){ind.rand[[k]]$slop.fix})
  r2.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$r2})
  aic.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$aic})
  p.ran = rbind(slop.ran,r2.ran,aic.ran) 
  
  if(slope.obs<0){for (k in 1:719){if(p.ran[1,k]>0){p.ran[2,k]<-NA;p.ran[3,k]<-NA;p.ran[4,k]<-NA;p.ran[5,k]<-NA}}
  }else{for (k in 1:719){if(p.ran[1,k]<0){p.ran[2,k]<-NA;p.ran[3,k]<-NA;p.ran[4,k]<-NA;p.ran[5,k]<-NA}}}
 
  EPS <- sqrt(.Machine$double.eps)
  slop.ran = p.ran[1,]
  r2.ran = p.ran[2:3,]
  aic.ran = p.ran[4:5,]
  
  if(slope.obs<0){
    p.slope=(sum(slop.ran<=(slope.obs[1]+EPS))+1)/(length(slop.ran)+1)
  }else{
    p.slope=(sum(slop.ran>=(slope.obs[1]-EPS))+1)/(length(slop.ran)+1)
  }
  
  p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS),na.rm=TRUE)+1)/(ncol(r2.ran)+1)
  p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS),na.rm=TRUE)+1)/(ncol(aic.ran)+1)
  p.values=rbind(c(p.slope),matrix(p.r2,2,1),matrix(p.aic,2,1))
  rownames(p.values)=c("P.Slope","P.R2M","P.R2C","P.AIC1","P.AIC2")
  colnames(p.values)=colnames(divindex)[i]
  
  ### Results output
  output=rbind(ndi,p.values)
  colnames(output)=paste0(colnames(divindex)[i],".Dissimilarity")
  if (i == 1){
    result = NULL
    result <- output
  }else{
    result <- cbind(result, output)
  }
  
  result <- as.data.frame(result)
  result
}

write.csv(result,"Allnetwork-LMM-Permutation-sharednodes.csv")
