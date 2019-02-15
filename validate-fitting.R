ntotest <- 50
perturbations <- 1000
strengths <- 1:500

source('define-model-fit.R')

lik <- function(true.f, t.inf, estimate, prob){
  if (any(true.f==0)){
    zeros <- which(true.f==0)
    true.f[zeros] <- min(true.f[-zeros])
  }
  
  return(-sum(dbinom(obs.freq, size=estimate, prob=prob, log=TRUE)) -
           sum(dmultinom(c(estimate,obs.inf), size=sum(estimate)+obs.inf, prob=c(true.f,t.inf), log=TRUE)))
}

deltas <- matrix(numeric(perturbations*length(strengths)),nrow=perturbations)

i <- 1
for(d in distnames){
  
  load(file=paste('results/fits-',d,'.Rda',sep=""))
  load(eval(paste('results/sim-freqs-',d,'.Rda',sep="")))
  
  for(i in 1:ntotest){
    print(i)
    
    est <- fits[[i]]$est
    obs <- fits[[i]]$prob
    
    f <- all.freqs[[i]]
    l <- sum(f)
    true.freq <- f/l
    true.inf <- max(1/nsim,sum(true.freq[(cutoff+1):maxsize]))
    true.freq <- true.freq[1:cutoff]
    
    likelihood <- lik(true.freq, true.inf, est, obs)
    for(k in 1:length(strengths)){
    strength <- strengths[k]
    perturbed.lik <- numeric(perturbations)
      for(j in 1:perturbations){
        new.est <- est + 
                   rbinom(cutoff,size=strength,prob=true.freq) - 
                   rbinom(cutoff,size=strength,prob=true.freq) 
        new.est[which(new.est<obs.freq)] <- obs.freq[which(new.est<obs.freq)]
        perturbed.lik[j] <- lik(true.freq, true.inf, new.est, obs)
        if(perturbed.lik[j] - likelihood < 0) {
          print(new.est)
          print(strength)
        }
      }
    deltas[,k] <- perturbed.lik - likelihood
    }
    ans <- any(deltas < 0)
    print(paste('weighted, both directions: any better? ', ans))
    if(ans) print(deltas[which(deltas<0)])
   
    for(k in 1:length(strengths)){
      strength <- strengths[k]
      perturbed.lik <- numeric(perturbations)
      for(j in 1:perturbations){
        new.est <- est + 
          rbinom(cutoff,size=strength,prob=true.freq) 
        new.est[which(new.est<obs.freq)] <- obs.freq[which(new.est<obs.freq)]
        perturbed.lik[j] <- lik(true.freq, true.inf, new.est, obs)
      }
      deltas[,k] <- perturbed.lik - likelihood
    }
    ans <- any(deltas < 0)
    print(paste('weighted, only larger: any better? ', ans))
    if(ans) print(deltas[which(deltas<0)])
    
    for(k in 1:length(strengths)){
      strength <- strengths[k]
      perturbed.lik <- numeric(perturbations)
      for(j in 1:perturbations){
        new.est <- est - 
          rbinom(cutoff,size=strength,prob=true.freq) 
        new.est[which(new.est<obs.freq)] <- obs.freq[which(new.est<obs.freq)]
        perturbed.lik[j] <- lik(true.freq, true.inf, new.est, obs)
      }
      deltas[,k] <- perturbed.lik - likelihood
    }
    ans <- any(deltas < 0)
    print(paste('weighted, only smaller: any better? ', ans))
    if(ans) print(deltas[which(deltas<0)])
    
    for(k in 1:length(strengths)){
      strength <- strengths[k]
      perturbed.lik <- numeric(perturbations)
      for(j in 1:perturbations){
        new.est <- est - 
          rbinom(cutoff,size=strength,prob=1/cutoff) +
          rbinom(cutoff,size=strength,prob=1/cutoff) 
        new.est[which(new.est<obs.freq)] <- obs.freq[which(new.est<obs.freq)]
        perturbed.lik[j] <- lik(true.freq, true.inf, new.est, obs)
      }
      deltas[,k] <- perturbed.lik - likelihood
    }
    ans <- any(deltas < 0)
    print(paste('unweighted, both directions: any better? ', ans))
    if(ans) print(deltas[which(deltas<0)])
  }
}