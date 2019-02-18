# Author: Emma E. Glennon
# Last edited: 15 February 2019
# Purpose: estimate number of undetected Ebola virus spillover events and small
#          outbreaks based on three estimated offspring distributions
# Contact: eeg31@cam.ac.uk/eeglennon@gmail.com
#
#*******************************************************************************

#SETTINGS***********************************************************************
#options: 1) setting from.scratch to true will do the full analysis from
#         scratch, starting from randomly sampling values of R0 and k for the
#         secondary infection distributions
#         2) setting generate.outbreaks to true (with from.scratch=FALSE) 
#         will use pre-generated values of R0 and k but will simulate
#         outbreaks to generate outbreak size distributions (this is
#         the slowest step of analysis; set nsim carefully)
#         3) setting fit.models to TRUE (with above two settings FALSE) will fit
#         the observation models with pre-loaded simulated outbreak size
#         distributions and estimates of R0, k
#         4) sensitivity.analysis, when set to TRUE, will cause the following 
#         sensitivity analyses to be repeated: validation of expectation step 
#         from maximum likelihood estimation (via perturbations of final estimates);
#         simulating outbreaks and observation process until 13 outbreaks have been 
#         observed (as in the data); extended model fitting with decay constants
#         and wider R0 and k ranges (as set below; also requires fit.models=TRUE)
# note: with all settings set to FALSE, existing results will be plotted, but no
#       analysis will be completed

generate.outbreaks <- FALSE
fit.models <- FALSE 
from.scratch <- FALSE 
sensitivity.analysis <- FALSE #whether to repeat sensitivity analysis

#ADDITIONAL CONTROL PARAMETERS AND CONSTANTS
nsim <- 10000         #number of simulations
nsets <- 500          #number of parameter sets

cutoff <- 57          #outbreak size cutoff
maxgen <- 50          #maximum number of generations for which to simulate an outbreak
maxsize <- cutoff + 1 #maximum number of individuals for which to simulate an outbreak

#settings for sensitivity analyses:
nToSim <- 10000     #number of times to simulate outbreaks and fit observation process
SAsets <- 100       #number of offpsring dists to use when repeating fitting across cutoff values
cuts <- seq(5,55,5) #cutoff values for use in cutoff sensitivity analysis
decays <- seq(0,.3,0.1)    #generational decay constants for sensitivity analysis
R0s.SA <- seq(0.5,3,0.25)  #expanded R0 values for sensitivity analysis
ks.SA <- seq(0.1,2.1,0.25) #expanded R0 values for sensitivity analysis


#*******************************************************************************
# EDITING BELOW THIS LINE WILL CHANGE THE ANALYSIS
#*******************************************************************************

#setup and data loading*********************************************************
#packages required for plotting (minimum functionality)
require(ggplot2)
require(reshape2)
require(cowplot)
require(scales)

if(from.scratch) {
  generate.outbreaks <- TRUE
  fit.models <- TRUE
}
if(generate.outbreaks){
  fit.models <- TRUE
  require(igraph)
}
if(fit.models){
  require(foreach)  #allows for easy parallelization if needed
  require(numDeriv)
}

#setup by reading and transforming data and defining plotting function
source('setup.R')

#get sizes of observed epidemic sizes smaller than cutoff
obs.dist <- known$cases
obs.freq <- numeric(cutoff+1)
for (i in 1:cutoff){
  obs.freq[i] <- length(which(obs.dist==i))
}
#number of reported outbreaks larger than cutoff:
obs.inf <- length(which(obs.dist>cutoff))
#sizes of reported outbreaks smaller than cutoff:
obs.dist <- obs.dist[which(obs.dist < cutoff)]

xmax <- 14
spillseen <- obs.freq[1]
allseen <- sum(obs.freq)+obs.inf
obs.freq.all <- obs.freq
obs.freq <- obs.freq[1:cutoff]

#****************************************************************************************
distnames <- c('FO','SL','Guin')

if(from.scratch) {
  #Full 2013-16 outbreak (WHO Ebola Response Team, PLOS Medicine)
  R0s <- runif(nsets, 1, 1.5)
  sizes <- runif(nsets, .03, 0.52)
  save(R0s, file='results/R0s-FO.Rda')
  save(sizes, file='results/sizes-FO.Rda')
  
  #Western Area, Sierra Leone (Lau et al., PNAS)
  #estimated to match credible interval
  R0s <- rlnorm(nsets,log(2.39),0.09) 
  #variance estimated from figure; not given in paper
  sizes <- rnorm(nsets, 0.37, 0.025)
  save(R0s, file='results/R0s-SL.Rda')
  save(sizes, file='results/sizes-SL.Rda')
  
  #Conakry, Guinea (Faye et al., Lancet ID; Althaus Lancet ID)
  #R0s and sizes estimated to match confidence intervals
  R0s <- rlnorm(nsets,log(0.95),0.09) #95% CI 0·57–1·34
  sizes <- rnorm(nsets, 0.18, 0.025) #95% CI 0·10–0·26
  save(R0s, file='results/R0s-Guin.Rda')
  save(sizes, file='results/sizes-Guin.Rda')
} 

if(generate.outbreaks) {
  for(d in distnames) {
    load(eval(paste('results/R0s-',d,'.Rda',sep="")))
    load(eval(paste('results/sizes-',d,'.Rda',sep="")))
    distname <- d
    source('sim-by-name.R')
  }
}

#fit observation models with logistic linking function
if(fit.models) {
  source('define-model-fit.R')
  for(d in distnames) {
    load(eval(paste('results/sim-freqs-',d,'.Rda',sep="")))
    attempt <- tryFit(beta.init=5,alpha.init=1,f=all.freqs,obs.fun="logistic")
    estimates <- attempt$estimates
    fits <- attempt$fits
    save(estimates, file=paste('results/estimates-',d,'.Rda',sep=""))
    save(fits, file=paste('results/fits-',d,'.Rda',sep=""))
  }
}
#plot estimates (with logistic linking function)
source('plot01.R')

#fit observation models with geometric linking function
if(fit.models) {
  source('define-model-fit.R')
  for(d in distnames) {
    load(eval(paste('results/sim-freqs-',d,'.Rda',sep="")))
    attempt <- tryFit(p.init=.01,f=all.freqs,obs.fun="geometric")
    estimates <- attempt$estimates
    fits <- attempt$fits
    save(estimates, file=paste('results/estimates-geom-',d,'.Rda',sep=""))
    save(fits, file=paste('results/fits-geom-',d,'.Rda',sep=""))
  }
}

#compare AICc of different observation functions
source('plot02.R')

#STATISTICS****************************************************************************
source('print-statistics.R')

#Cutoff sensitivity analysis************************************************************

if(fit.models) {
  source('cutoff-sensitivity.R')
} else {
  load(file='./results/SA.estlist.Rda')
}
source('plot03.R')

#Reff decay and R0/k expansion sensitivity analysis**************************************
R0s <- R0s.SA
sizes <- ks.SA
obs.freq <- obs.freq[1:cutoff]

if(sensitivity.analysis){
  if(generate.outbreaks) source('simulate-outbreaks-SA.R')
  
  if(fit.models){
    est.matrix <- matrix(numeric(length(sizes)*length(R0s)),
                         nrow=length(sizes))
    nsets <- length(sizes)
    source('define-model-fit.R')
    est.list <- list()
    for(d in 1:length(decays)){
      decay <- decays[d]
      load(file=paste('results/sim-freqs-SA-decay',decay,'.Rda',sep=""))
      sub.list <- list()
      for(i in 1:length(R0s)){
        est <- tryFit(beta.init=5,alpha.init=1,f=SA.freqs[[i]])
        overall <- allseen/rowSums(est$estimates)
        est.matrix[,i] <- overall
        sub.list[[i]] <- est
        
        print(paste(round(((d-1)*length(R0s)+i)/(length(R0s)*length(decays))*100,2),
                    '% of sensitivity analysis complete',sep=""))
      }
      
      est.list[[d]] <- sub.list
    }
    save(est.list,file='results/fit-SA-all.Rda')    
  }
}

#Goodness of fit test and expectation step validation************************************
if(sensitivity.analysis) {
  source('validate-fitting.R')
  source('simulate-13.R')
}

source('plot04.R') #plot results of all sensitivity analyses
