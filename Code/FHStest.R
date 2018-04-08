library(igraph)
library(ColorPalette)
library(RColorBrewer)
##
source("Code/MoranI.R")
source("Code/Phi.R")
##
RdPalette = colorRampPalette(brewer.pal(9,"Reds"))(20)
BlPalette = colorRampPalette(brewer.pal(9,"Blues"))(20)
GnPalette = colorRampPalette(brewer.pal(9,"Greens"))(20)
## read data
#network_c2 = read.table("data/phs000153_c2.txt", sep = "\t", header = TRUE)
#psy.data = read.table("data/phs000007.v29.pht000100.v6.p10.c2.psych1_3s.HMB-IRB-NPU-MDS.txt", sep = "\t", header = TRUE)
#network2 = network_c2[(network_c2$SPELLBEGIN == 1),]
#pheno_c2_ex1_5 = read.table("data/phs000007.v29.pht000034.v7.p10.c2.ex1_5s.HMB-IRB-NPU-MDS.txt",sep = "\t", header = TRUE)
#pheno.data = read.table("data/phs000007.v29.pht000034.v7.p10.c2.ex1_5s.HMB-IRB-NPU-MDS.txt",sep = "\t", header = TRUE)
pheno.info = pheno.data
network_c2_1.5 = network_c2[(network_c2$SPELLBEGIN <= 12*(95-71+1)) & (network_c2$SPELLEND > 12*(91-71)),  ]
Adj = matrix(0, nrow = nrow(pheno.info), ncol = nrow(pheno.info))
ids = pheno.info$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network_c2_1.5[which(network_c2_1.5$shareid == ids[i]), 34])
  friends2 = which(ids %in% network_c2_1.5[which(network_c2_1.5$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 
  Adj[friends2, i ] = 1
}
sum(Adj)/2

## Continuous
# blood pressure : systolic - 1st md reading  
focus.data = cbind(pheno.info$shareid, pheno.info$E485, pheno.info$E486)
focus.data = as.data.frame(na.omit(focus.data))
names(focus.data) = c("shareid", "systolic", "diastolic")
focus.data = as.data.frame(na.omit(focus.data))
network_c2_1.5 = network_c2[(network_c2$SPELLBEGIN <= 12*(95-71+1)) & (network_c2$SPELLEND > 12*(91-71)),  ]
Adj = matrix(0, nrow = nrow(focus.data), ncol = nrow(focus.data))
ids = focus.data$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network_c2_1.5[which(network_c2_1.5$shareid == ids[i]), 34])
  friends2 = which(ids %in% network_c2_1.5[which(network_c2_1.5$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 
  Adj[friends2, i ] = 1
}
bp.moran = make.permute.moran(Adj, focus.data$systolic, 500)

## Figure
Y = focus.data$systolic
A = Adj
G = graph.adjacency(A, "undirected")
V(G)$Y = Y
subgraph = induced.subgraph(G, components(G)$membership == 5)
members = which(components(G)$membership == 5)
for(i in 1:20){
  lower = quantile(V(subgraph)$Y, 0.05*(i-1))
  upper = quantile(V(subgraph)$Y, 0.05*i)
  V(subgraph)[V(subgraph)$Y >= lower & V(subgraph)$Y < upper]$color = BlPalette[i]
}
V(subgraph)[which.max(V(subgraph)$Y)]$color = BlPalette[20]


igraph.options(vertex.size = 4, edge.arrow.size = 0.1,
               vertex.label = NULL)
pdf("Figure/bloodpressure.pdf", height = 14, width = 10)
par(mar=c(5,0,5,0), cex.main = 3)
set.seed(123)
plot(subgraph, layout = layout.fruchterman.reingold,
     vertex.label= "", main = "Systolic blood pressure")
mtext(paste("Moran's I: ", formatC(bp.moran[1], digits = 2, format = "f")), 
      side = 1, cex = 2, line = -1)
mtext(paste("p-value (permutation): ", formatC(bp.moran[3], digits = 2, format = "f")), 
      side = 1, cex = 2, line = 1)
dev.off()


# blood pressure : diastolic - 1st md reading  
dia.moran = make.permute.moran(Adj, focus.data$diastolic, 500)
## Figure
Y = focus.data$diastolic
A = Adj
G = graph.adjacency(A, "undirected")
V(G)$Y = Y
subgraph = induced.subgraph(G, members)
for(i in 1:20){
  lower = quantile(V(subgraph)$Y, 0.05*(i-1))
  upper = quantile(V(subgraph)$Y, 0.05*i)
  V(subgraph)[V(subgraph)$Y >= lower & V(subgraph)$Y < upper]$color = GnPalette[i]
}
V(subgraph)[which.max(V(subgraph)$Y)]$color = GnPalette[20]

igraph.options(vertex.size = 4, edge.arrow.size = 0.1,
               vertex.label = NULL)
pdf("Figure/diastolic.pdf", height = 14, width = 10)
par(mar=c(5,0,5,0), cex.main = 3)
set.seed(123)
plot(subgraph, layout = layout.fruchterman.reingold,
     vertex.label= "", main = "Diastolic blood pressure")
mtext(paste("Moran's I: ", formatC(dia.moran[1], digits = 2, format = "f")), 
      side = 1, cex = 2, line = -1)
mtext(paste("p-value (permutation): ", formatC(dia.moran[3], digits = 2, format = "f")), 
      side = 1, cex = 2, line = 1)
dev.off()

## weight
focus.data = cbind(pheno.info$shareid, pheno.info$E024)
focus.data = as.data.frame(na.omit(focus.data))
names(focus.data) = c("shareid", "pheno")
focus.data = as.data.frame(na.omit(focus.data))
network_c2_1.5 = network_c2[(network_c2$SPELLBEGIN <= 12*(95-71+1)) & (network_c2$SPELLEND > 12*(91-71)),  ]
Adj = matrix(0, nrow = nrow(focus.data), ncol = nrow(focus.data))
ids = focus.data$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network_c2_1.5[which(network_c2_1.5$shareid == ids[i]), 34])
  friends2 = which(ids %in% network_c2_1.5[which(network_c2_1.5$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 
  Adj[friends2, i ] = 1
}
weight.moran = make.permute.moran(Adj, focus.data$pheno, 500)

Y = focus.data$pheno
A = Adj
G = graph.adjacency(A, "undirected")
V(G)$Y = Y
subgraph = induced.subgraph(G, components(G)$membership == 5)
for(i in 1:20){
  lower = quantile(V(subgraph)$Y, 0.05*(i-1))
  upper = quantile(V(subgraph)$Y, 0.05*i)
  V(subgraph)[V(subgraph)$Y >= lower & V(subgraph)$Y < upper]$color = RdPalette[i]
}
V(subgraph)[which.max(V(subgraph)$Y)]$color = RdPalette[20]


igraph.options(vertex.size = 4, edge.arrow.size = 0.1,
               vertex.label = NULL)
pdf("Figure/weight.pdf", height = 14, width = 10)
par(mar=c(5,0,5,0), cex.main = 3)
set.seed(123)
plot(subgraph, layout = layout.fruchterman.reingold,
     vertex.label= "", main = "Weight")
mtext(paste("Moran's I: ", formatC(weight.moran[1], digits = 2, format = "f")), 
      side = 1, cex = 2, line = -1)
mtext(paste("p-value (permutation): ", formatC(weight.moran[3], digits = 2, format = "f")), 
      side = 1, cex = 2, line = 1)
dev.off()

##################################
## employment status
# 1 : full time, 2: part time, 0 : no
focus.data = cbind(pheno.info$shareid, pheno.info$E067)
focus.data = as.data.frame(na.omit(focus.data))
names(focus.data) = c("shareid", "pheno")
focus.data = as.data.frame(na.omit(focus.data))
network_c2_1.5 = network_c2[(network_c2$SPELLBEGIN <= 12*(95-71+1)) & (network_c2$SPELLEND > 12*(91-71)),  ]
Adj = matrix(0, nrow = nrow(focus.data), ncol = nrow(focus.data))
ids = focus.data$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network_c2_1.5[which(network_c2_1.5$shareid == ids[i]), 34])
  friends2 = which(ids %in% network_c2_1.5[which(network_c2_1.5$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 
  Adj[friends2, i ] = 1
}
employ = make.permute.Phi(Adj, as.integer(focus.data$pheno + 1), 500) # 3.400979 0.000000 0.002000

Y = focus.data$pheno
A = Adj
G = graph.adjacency(A, "undirected")
V(G)$Y = Y
subgraph = induced.subgraph(G, components(G)$membership == 5)
for(i in 1:length(V(subgraph))){
  V(subgraph)$color[i] =  ifelse(V(subgraph)$Y[i] == 0, "firebrick1", "dodgerblue1")
  V(subgraph)$color[i] = ifelse(V(subgraph)$Y[i] == 1, "seagreen", V(subgraph)$color[i])
}

igraph.options(vertex.size = 4, edge.arrow.size = 0.1,
               vertex.label = NULL)
pdf("Figure/employed_cate.pdf", height = 14, width = 10)
par(mar=c(5,0,5,0), cex.main = 3)
set.seed(123)
plot(subgraph, layout = layout.fruchterman.reingold,
     vertex.label= "", main = "Employment status")
mtext(expression(paste(Phi,": 3.40")), 
      side = 1, cex = 2, line = -1)
mtext(paste("p-value (permutation): ", formatC(employ[3], digits = 4, format = "f")), 
      side = 1, cex = 2, line = 1)
legend("topright", c("Fulltime", "Parttime", "Not employed"), 
       cex = 2, pch = 20, col = c("dodgerblue1", "seagreen", "firebrick1"),
       bty = 'n', pt.cex = 2.5)
dev.off()

## COFFEE/CAFFEINATED - PREDOMINANT METHOD
focus.data = cbind(pheno.info$shareid, pheno.info$E303)
focus.data = as.data.frame(na.omit(focus.data))
names(focus.data) = c("shareid", "pheno")
focus.data = as.data.frame(na.omit(focus.data))
focus.data = focus.data[focus.data$pheno != 3 & focus.data$pheno != 8,]
network_c2_1.5 = network_c2[(network_c2$SPELLBEGIN <= 12*(95-71+1)) & (network_c2$SPELLEND > 12*(91-71)),  ]
Adj = matrix(0, nrow = nrow(focus.data), ncol = nrow(focus.data))
ids = focus.data$shareid
for(i in 1:nrow(Adj)){
  friends1 = which(ids %in% network_c2_1.5[which(network_c2_1.5$shareid == ids[i]), 34])
  friends2 = which(ids %in% network_c2_1.5[which(network_c2_1.5$sharealterid == ids[i]), 33])
  Adj[i, friends1] = 1
  Adj[i, friends2] = 1
  Adj[friends1, i] = 1 
  Adj[friends2, i ] = 1
}
coffee.method = make.permute.Phi(Adj, focus.data$pheno, 500) 

Y = focus.data$pheno
A = Adj
G = graph.adjacency(A, "undirected")
V(G)$Y = Y
subgraph = induced.subgraph(G, components(G)$membership == 5)
for(i in 1:length(V(subgraph))){
  V(subgraph)$color[i] =  ifelse(V(subgraph)$Y[i] == 0, "mediumspringgreen", "dodgerblue1")
  V(subgraph)$color[i] = ifelse(V(subgraph)$Y[i] == 1, "royalblue", V(subgraph)$color[i])
  V(subgraph)$color[i] =  ifelse(V(subgraph)$Y[i] == 2, "mediumvioletred", V(subgraph)$color[i])
  V(subgraph)$color[i] =  ifelse(V(subgraph)$Y[i] == 4, "yellow", V(subgraph)$color[i])
}

igraph.options(vertex.size = 4, edge.arrow.size = 0.1,
               vertex.label = NULL)
pdf("Figure/coffee_cate.pdf", height = 14, width = 10)
par(mar=c(5,0,5,0), cex.main = 3)
set.seed(321)
plot(subgraph, layout = layout.fruchterman.reingold,
     vertex.label= "", main = "Ways to make coffee")
mtext(expression(paste(Phi,": 3.21")), 
      side = 1, cex = 2, line = -1)
mtext(paste("p-value (permutation): ", formatC(coffee.method[3], digits = 4, format = "f")), 
      side = 1, cex = 2, line = 1)
legend("topright", c("Non-drinker", "Filter", "Perc", "Instant"), 
       cex = 2, pch = 20, 
       col = c("mediumspringgreen", "royalblue", "mediumvioletred", "yellow"),
       bty = 'n', pt.cex = 2.5)
dev.off()