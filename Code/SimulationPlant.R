library(igraph)
library(MASS)
library(geosphere)
library(spdep)
##
source("Code/MoranI.R")
source("Code/Phi.R")
## read data
load("Data/analysis_dat.RData")
latitude = analysis_dat$Fac.Latitude
longitude = analysis_dat$Fac.Longitude
dist.matrix = distm(cbind(longitude, latitude), cbind(longitude, latitude), fun = distHaversine)
weights = max(dist.matrix) / dist.matrix
diag(weights) = 0
summary(as.numeric(weights))
weights = ifelse(weights > 10, 10, weights)
## neighborhood choice (needed for join count analysis)
#
nb_24nn = knn2nb(knearneigh(cbind(analysis_dat$Fac.Longitude, analysis_dat$Fac.Latitude), k = 24))
nblist24 = nb2listw(nb_24nn)
#
nb_30nn = knn2nb(knearneigh(cbind(analysis_dat$Fac.Longitude, analysis_dat$Fac.Latitude), k = 30))
nblist30 = nb2listw(nb_30nn)
#
nb_48nn = knn2nb(knearneigh(cbind(analysis_dat$Fac.Longitude, analysis_dat$Fac.Latitude), k = 48))
nblist48 = nb2listw(nb_48nn)
#
nb_95nn = knn2nb(knearneigh(cbind(analysis_dat$Fac.Longitude, analysis_dat$Fac.Latitude), k = 95))
nblist95 = nb2listw(nb_95nn)
#
nb_142nn = knn2nb(knearneigh(cbind(analysis_dat$Fac.Longitude, analysis_dat$Fac.Latitude), k = 142))
nblist142 = nb2listw(nb_142nn)
##
nominal_neighbor = function(dist.mat, time_point, mprob, multip, effect.quantile){
  ### input
  # dist.mat : distance matrix
  # time_point : a vector of time point you want to observe the outcomes.
  # mprob : a maximum susceptibability probability.
  # multip : multinomial probability.
  # effect.quantile : quantile of distances under the influence 
  ### output
  #  an array of observations at time_point.
  
  popn = nrow(dist.mat)
  # popn : the number of subjects (n)
  
  max_t = time_point[length(time_point)] + 1  # maximum number of time for observation
  max_t = as.integer(max_t)
  
  outcome = matrix(0, ncol = popn, nrow = max_t) # initialize the outcome data frame
  
  # Generate observations at initial stage
  k = length(multip)
  for(i in 1:popn){
    dumb = rmultinom(1, size = 1, multip)
    outcome[1,i] = which(dumb == 1) 
  }
  
  for (t in 2:max_t){
    for (i in 1:popn){  
      p = runif(1, 0, mprob) # Generate a new susceptibility probability
      dummy = rbinom(1,1,p)
      
      if(dummy == 0){
        outcome[t,i] = outcome[(t-1),i]
      }else{
        ind = dist.mat[,i] <= quantile(dist.mat[,i], effect.quantile)
        new.prop = table(outcome[t-1, ind])
        if(length(new.prop) == 0){
          outcome[t,i] = outcome[(t-1), i]
        }else{
          outcome[t,i] = as.integer(sample(names(new.prop), 1, prob = new.prop))
        }
      }
      
    }
  }
  
  return(outcome[time_point + 1,])
}

######### Specified neighborhood ##############
Phisim1_Phi = Phisim2_Phi = list()
PhiSim1_neighbors = PhiSim2_neighbors = list()
for(mm in 1:500){
  set.seed(mm) 
  
  multip =  c(0.1, 0.2, 0.3, 0.25, 0.15) 
  popn = nrow(dist.matrix)
  outcome0 = c()
  for(i in 1:popn){
    dumb = rmultinom(1, size = 1, multip)
    outcome0[i] = which(dumb == 1) 
  }
  
  outcomes = nominal_neighbor(dist.matrix, c(1:20), 0.5, c(0.1, 0.2, 0.3, 0.25, 0.15), 
                              effect.quantile = 0.1)
  # outcomes = nominal_neighbor(dist.matrix, c(1:20), 0.5, c(0.1, 0.2, 0.3, 0.25, 0.15), 
  # effect.quantile = 0.2) # under 20% influence range
  
  ## neighbor = 30
  # time = 0
  mc.test = joincount.mc(as.factor(outcome0), nblist30, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors30.time0 = list(mc.test, BB.stat, BB.pval)
  
  # time = 3
  mc.test = joincount.mc(as.factor(outcomes[3,]), nblist30, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors30.time3 = list(mc.test, BB.stat, BB.pval)
  
  
  # time = 5
  mc.test = joincount.mc(as.factor(outcomes[5,]), nblist30, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors30.time5 = list(mc.test, BB.stat, BB.pval)
  
  # time = 10
  mc.test = joincount.mc(as.factor(outcomes[10,]), nblist30, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors30.time10 = list(mc.test, BB.stat, BB.pval)
  
  # time = 20
  mc.test = joincount.mc(as.factor(outcomes[20,]), nblist30, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors30.time20 = list(mc.test, BB.stat, BB.pval)
  
  neighbors30 = list(neighbors30.time0, neighbors30.time3, neighbors30.time5, neighbors30.time10, neighbors30.time20)
  
  
  ## neighbor = 48
  # time = 0
  mc.test = joincount.mc(as.factor(outcome0), nblist48, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors48.time0 = list(mc.test, BB.stat, BB.pval)
  
  # time = 3
  mc.test = joincount.mc(as.factor(outcomes[3,]), nblist48, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors48.time3 = list(mc.test, BB.stat, BB.pval)
  
  
  # time = 5
  mc.test = joincount.mc(as.factor(outcomes[5,]), nblist48, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors48.time5 = list(mc.test, BB.stat, BB.pval)
  
  # time = 10
  mc.test = joincount.mc(as.factor(outcomes[10,]), nblist48, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors48.time10 = list(mc.test, BB.stat, BB.pval)
  
  # time = 20
  mc.test = joincount.mc(as.factor(outcomes[20,]), nblist48, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors48.time20 = list(mc.test, BB.stat, BB.pval)
  
  neighbors48 = list(neighbors48.time0, neighbors48.time3, neighbors48.time5, neighbors48.time10, neighbors48.time20)
  
  
  ## neighbor = 95
  # time = 0
  mc.test = joincount.mc(as.factor(outcome0), nblist95, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors95.time0 = list(mc.test, BB.stat, BB.pval)
  
  # time = 3
  mc.test = joincount.mc(as.factor(outcomes[3,]), nblist95, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors95.time3 = list(mc.test, BB.stat, BB.pval)
  
  
  # time = 5
  mc.test = joincount.mc(as.factor(outcomes[5,]), nblist95, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors95.time5 = list(mc.test, BB.stat, BB.pval)
  
  # time = 10
  mc.test = joincount.mc(as.factor(outcomes[10,]), nblist95, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors95.time10 = list(mc.test, BB.stat, BB.pval)
  
  # time = 20
  mc.test = joincount.mc(as.factor(outcomes[20,]), nblist95, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors95.time20 = list(mc.test, BB.stat, BB.pval)
  
  
  neighbors95 = list(neighbors95.time0, neighbors95.time3, neighbors95.time5, neighbors95.time10, neighbors95.time20) 
  
  
  ## neighbor = 142
  # time = 0
  mc.test = joincount.mc(as.factor(outcome0), nblist142, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors142.time0 = list(mc.test, BB.stat, BB.pval)
  
  # time = 3
  mc.test = joincount.mc(as.factor(outcomes[3,]), nblist142, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors142.time3 = list(mc.test, BB.stat, BB.pval)
  
  
  # time = 5
  mc.test = joincount.mc(as.factor(outcomes[5,]), nblist142, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors142.time5 = list(mc.test, BB.stat, BB.pval)
  
  # time = 10
  mc.test = joincount.mc(as.factor(outcomes[10,]), nblist142, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors142.time10 = list(mc.test, BB.stat, BB.pval)
  
  # time = 20
  mc.test = joincount.mc(as.factor(outcomes[20,]), nblist142, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors142.time20 = list(mc.test, BB.stat, BB.pval)
  
  neighbors142 = list(neighbors142.time0, neighbors142.time3, neighbors142.time5, neighbors142.time10, neighbors142.time20)   
  
  Phisim1_Phi[[mm]] = list(Phi0, Phi5, Phi10, Phi50, Phi100)
  
  PhiSim1_neighbors[[mm]] = list(neighbors30, neighbors48, neighbors95, neighbors142)
  ## PhiSim2_neighbors[[mm]] = list(neighbors30, neighbors48, neighbors95, neighbors142)
}

save(PhiSim1_neighbors, file = "Data/PhiSim1_neighbors.RData")
save(PhiSim1_Phi, file = "Data/PhiSim1_Phi.RData")



########## Autocorrelated error model ##############
autoerror0_neighbors = autoerror1_neighbors = autoerror2_neighbors = autoerror3_neighbors = list()
cor0_Phi = cor1_Phi = cor2_Phi = cor3_Phi = list()
cor0_truePhi = cor1_truePhi = cor2_truePhi = cor3_truePhi = list()


for(mm in 1:500){
  set.seed(mm) 
  # weights matrix
  p = 1/max((dist.matrix)^(1/5))
  # p = 5/max((dist.matrix)^(1/2)) : autoerror2_neighbors
  # p = 10/max((dist.matrix)^(1/2)) : autoerror3_neighbors
  w = exp(-p*(dist.matrix)^(1/5))
  Ww = chol(w)
  
  # errors
  z = t(Ww) %*% rnorm(nrow(dist.matrix),0, 1) 
  z = as.numeric(t(z))
  quantile.point = quantile(z,c(0.1, 0.2+0.1, 0.3+0.1+0.2, 0.25+0.1+0.2+0.3, 0.15+0.1+0.2+0.3+0.25)) 
  # (0.1, 0.2, 0.3, 0.25, 0.15)
  cate.outcome = ifelse(z < quantile.point[1], 1, 5)
  cate.outcome = ifelse(z >= quantile.point[1] & z < quantile.point[2], 2, cate.outcome)
  cate.outcome = ifelse(z >= quantile.point[2] & z < quantile.point[3], 3, cate.outcome)
  cate.outcome = ifelse(z >= quantile.point[3] & z < quantile.point[4], 4, cate.outcome)
  ## null case : ind = sample(1:length(cate.outcome), replace = FALSE)
  ## null case : cate.outcome = cate.outcome[ind]
  
  ## neighbor = 24
  mc.test = joincount.mc(as.factor(cate.outcome), nblist24, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors24 = list(mc.test, BB.stat, BB.pval)
  
  
  ## neighbor = 48
  mc.test = joincount.mc(as.factor(cate.outcome), nblist48, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors48 = list(mc.test, BB.stat, BB.pval)
  
  ## neighbor = 95
  mc.test = joincount.mc(as.factor(cate.outcome), nblist95, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors95 = list(mc.test, BB.stat, BB.pval)
  
  ## neighbor = 142
  mc.test = joincount.mc(as.factor(cate.outcome), nblist142, nsim = 500)
  BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
              mc.test[[3]]$statistic, mc.test[[4]]$statistic, mc.test[[5]]$statistic)
  BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
              mc.test[[3]]$p.value,  mc.test[[4]]$p.value,  mc.test[[5]]$p.value)
  neighbors142 = list(mc.test, BB.stat, BB.pval)
  
  
  co1_Phi[[mm]] = make.permute.Phi(weights, cate.outcome, 500)
  cor1_truePhi[[mm]] = Phi.result = make.permute.Phi(w, cate.outcome, 500)
  
  # autoerror1_neighbors[[mm]] = list(neighbors24, neighbors48, neighbors95, neighbors142) : nullcase
  autoerror1_neighbors[[mm]] = list(neighbors24, neighbors48, neighbors95, neighbors142)
  #autoerror2_neighbors[[mm]] = list(neighbors24, neighbors48, neighbors95, neighbors142)
  #autoerror3_neighbors[[mm]] = list(neighbors24, neighbors48, neighbors95, neighbors142)
}

save(cor1_Phi, file = "Data/cor1_Phi.RData")
save(cor1_truePhi, file = "Data/cor1_truePhi.RData")
save(autoerror1_neighbors, file = "Data/autoerror1_neighbors.RData")