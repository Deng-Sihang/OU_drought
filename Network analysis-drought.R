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
