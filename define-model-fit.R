#likelihood of observing reported outbreaks give a 'true' distribution of epidemic sizes
getObs <- function(true.f,     #frequencies of outbreak sizes (from simulation)
                   obs.freq,   #frequencies of outbreak sizes (from reported outbreaks)
                   beta=NULL,       #estimate of observation function parameter 1
                   alpha=NULL,      #estimate of observation function parameter 2
                   p=NULL,
                   cutoff=cutoff, 
                   t.inf=NULL, #proportion of simulated outbreaks larger than cutoff
                   obs.inf=obs.inf, #number of outbreaks larger than cutoff (100% detection assumed)
                   obs.fun="logistic",
                   npert=1000,
                   strengthmax=100,
                   strengths=10
){
  init.est <- NA
  
  #assume there are no true zero frequencies, just not enough simulations, 
  #and set all zeros to the minimum nonzero frequency 
  #(rarely necessary given high number of simulations; edge case handling)
  if (any(true.f==0)){
    zeros <- which(true.f==0)
    if(min(true.f) > 0) {
      true.f[zeros] <- min(true.f[-zeros])
    } else{
      true.f[zeros] <- 1/nsim
    }
  }
  if (t.inf==0) t.inf <- min(true.f)

  lik <- function(estimate, prob){
    out <- -sum(dbinom(obs.freq, size=estimate, prob=prob, log=TRUE)) -
      sum(dmultinom(c(estimate,obs.inf), size=sum(estimate)+obs.inf, prob=c(true.f,t.inf), log=TRUE))

    return(out)
  }

  EM <- function(par){
    tol <- 1e-10
    its <- 1000

    #to store 'true' outbreaks smaller than cutoff
    index <- which(obs.freq>0)
    true.size <- numeric(cutoff)
    
    loglik <- rep(Inf,its)
    loglik[1] <- Inf
    loglik[2] <- 1e50
    
    k <- 3
    while(k <= its && (loglik[k-2]-loglik[k-1])>=tol){
      #probability of observation as function of epidemic size
      obs <- NA
      if(obs.fun=="logistic") {
        beta <- par[1]
        alpha <- exp(par[2])/(1+exp(par[2]))*100
        obs <- (1+exp(beta-1:cutoff))^(-alpha)
      }
      if(obs.fun=="geometric"){
        p <- par[1]
        obs <- 1-(1-p)^(1:cutoff)
      }
      
      #initial estimate of T to be adjusted to maximize likelihood
      N.min <- sum(obs.freq) + obs.inf
      N.est <- round(obs.inf/t.inf)
      init.est <- qbinom(0.5,size=N.est,prob=true.f)
      init.est[which(init.est<obs.freq)] <- obs.freq[which(init.est<obs.freq)]
      init.est[which(obs==1)] <- obs.freq[which(obs==1)]
      
      N.est2 <- round(obs.inf/t.inf)
      init.est2 <- qbinom(0.5,size=N.est2,prob=true.f)
      init.est2[which(init.est2<obs.freq)] <- obs.freq[which(init.est2<obs.freq)]
      init.est2[which(obs==1)] <- obs.freq[which(obs==1)]
      if(lik(init.est2,obs) < lik(init.est,obs)) init.est <- init.est2
      rm(init.est2)
      
      N.est <- sum(init.est) + obs.inf
      
      #EXPECTATION STEP
      #adjust estimates until best estimate (T) is reached
      stop <- FALSE
      while(!stop){
        increment <- rep(-Inf, cutoff)
        decrement <- rep(-Inf, cutoff)
        
        for(i in 1:cutoff){
          vals <- lik(init.est, obs)

          new.est <- init.est
          new.est[i] <- init.est[i] + 1
          increment[i] <- vals - lik(new.est, obs)

          if(init.est[i] > obs.freq[i] & N.est > N.min){
            new.est[i] <- init.est[i] - 1
            decrement[i] <- vals - lik(new.est, obs)
          }
        }
        
        if(any(increment > 0)){
          if(max(increment) > max(decrement)){
            init.est[which(increment==max(increment))[1]] <- 
              init.est[which(increment==max(increment))[1]] + 1
            N.est <- N.est + 1
          } else {
            init.est[which(decrement==max(decrement))[1]] <- 
              init.est[which(decrement==max(decrement))[1]] - 1
            N.est <- N.est - 1
          }
        } else if (any(decrement > 0)){
          init.est[which(decrement==max(decrement))[1]] <- 
            init.est[which(decrement==max(decrement))[1]] - 1
          N.est <- N.est - 1
        } else {
          #if at a local likelihood maximum, perturb to see if any changes greater than
          #a single change in the current estimate lead to higher likelihoods
          best.est <- init.est
          for(strength in round(seq(1,strengthmax,length.out=strengths))){
            for(pert in 1:npert) {
            perturbed <- init.est + 
              rbinom(cutoff,size=strength,prob=true.f) - 
              rbinom(cutoff,size=strength,prob=true.f)
            perturbed[which(perturbed<obs.freq)] <- obs.freq[which(perturbed<obs.freq)]
            if(lik(perturbed, obs) - lik(best.est, obs) < 0) best.est <- perturbed
            
            perturbed <- init.est + 
              rbinom(cutoff,size=strength,prob=true.f)
            perturbed[which(perturbed<obs.freq)] <- obs.freq[which(perturbed<obs.freq)]
            if(lik(perturbed, obs) - lik(best.est, obs) < 0) best.est <- perturbed
            }
          }
          
          if(all(best.est == init.est)) {
            #if perturbations do not change the likelihood, progress to maximization step
            stop <- TRUE
          } else {
            #if perturbations lead to a higher likelihood, search again for incremental changes
            init.est <- best.est
            N.est <- sum(init.est)+obs.inf
            print('perturbed')
          }
        }
      }
      
      true.size <- init.est
      N <- sum(true.size)+obs.inf
      
      #MAXIMIZATION STEP: given T, find par to maximize L
      logL <- function(par.new){
        obs <- NA
        if(obs.fun=="logistic") {
          beta <- par.new[1]
          alpha <- exp(par.new[2])/(1+exp(par.new[2]))*100
          obs <- (1+exp(beta-1:cutoff))^(-alpha)
        }
        if(obs.fun=="geometric"){
          p <- par.new[1]
          obs <- 1-(1-p)^(1:cutoff)
        }
        
        return(lik(true.size,obs))
      }
      
      maxed <- NA
      if(obs.fun=="logistic") maxed <- optim(fn=logL, par=par)
      if(obs.fun=="geometric")maxed <- optim(fn=logL, par=par, method="Brent", lower=0, upper=.9)
      par <- maxed$par
      if(loglik[k] != 1e50) loglik[k] <- maxed$value
      k <- k + 1
    }
    
    par <- maxed$par
    obs <- NA
    beta <- NA
    alpha <- NA
    p <- NA
    if(obs.fun=="logistic") {
      beta <- par[1]
      alpha <- exp(par[2])/(1+exp(par[2]))*100
      obs <- (1+exp(beta-1:cutoff))^(-alpha)
    }
    if(obs.fun=="geometric"){
      p <- par[1]
      obs <- 1-(1-p)^(1:cutoff)
    }
    
    N.min <- sum(obs.freq) + obs.inf
    N.est <- round(obs.inf/t.inf)
    init.est <- qbinom(0.5,size=N.est,prob=true.f)
    init.est[which(init.est<obs.freq)] <- obs.freq[which(init.est<obs.freq)]
    init.est[which(obs==1)] <- obs.freq[which(obs==1)]
    
    N.est <- sum(init.est) + obs.inf
    stop <- FALSE
    while(!stop){
      increment <- rep(-Inf, cutoff)
      decrement <- rep(-Inf, cutoff)
      
      for(i in 1:cutoff){
        vals <- lik(init.est, obs)
        new.est <- init.est
        new.est[i] <- init.est[i] + 1
        increment[i] <- vals - lik(new.est, obs)
        
        if(init.est[i] > obs.freq[i] & N.est > N.min){
          new.est[i] <- init.est[i] - 1
          decrement[i] <- vals - lik(new.est, obs)
        }
      }
      
      if(any(increment > 0)){
        if(max(increment) > max(decrement)){
          init.est[which(increment==max(increment))[1]] <- 
            init.est[which(increment==max(increment))[1]] + 1
          N.est <- N.est + 1
        } else {
          init.est[which(decrement==max(decrement))[1]] <- 
            init.est[which(decrement==max(decrement))[1]] - 1
          N.est <- N.est - 1
        }
      } else if (any(decrement > 0)){
        init.est[which(decrement==max(decrement))[1]] <- 
          init.est[which(decrement==max(decrement))[1]] - 1
        N.est <- N.est - 1
      } else {
        #if at a local likelihood maximum, perturb to see if any changes greater than
        #a single change in the current estimate lead to higher likelihoods
        best.est <- init.est
        for(strength in round(seq(1,strengthmax,length.out=strengths))){
          for(pert in 1:npert) {
            perturbed <- init.est + 
              rbinom(cutoff,size=strength,prob=true.f) - 
              rbinom(cutoff,size=strength,prob=true.f)
            perturbed[which(perturbed<obs.freq)] <- obs.freq[which(perturbed<obs.freq)]
            if(lik(perturbed, obs) - lik(best.est, obs) < 0) best.est <- perturbed
            
            perturbed <- init.est + 
              rbinom(cutoff,size=strength,prob=true.f)
            perturbed[which(perturbed<obs.freq)] <- obs.freq[which(perturbed<obs.freq)]
            if(lik(perturbed, obs) - lik(best.est, obs) < 0) best.est <- perturbed
          }
        }
        
        if(all(best.est == init.est)) {
          #if perturbations do not change the likelihood, progress to maximization step
          stop <- TRUE
        } else {
          #if perturbations lead to a higher likelihood, search again for incremental changes
          init.est <- best.est
          N.est <- sum(init.est)+obs.inf
          print('perturbed')
        }
      }
    }
    
    true.size <- init.est
  
    return(list(beta=beta, alpha=alpha, p=p, value=min(loglik), est=true.size, prob=obs))
  }

  fit <- NA
  if(obs.fun=="logistic") {
    fit <- EM(par=c(beta,log(alpha/(1-alpha))))
  }
  if(obs.fun=="geometric"){
    fit <- EM(par=p)
  }
  
  est <- fit$est
  N.est <- sum(est) + obs.inf
 
  return(list(beta=fit$beta,alpha=fit$alpha,p=fit$p,
              N=N.est,val=fit$value,est=est,prob=fit$prob))
}

#fit observation function and get corresponding estimates of true outbreak size dist
tryFit <- function(beta.init=NULL,   #initial estimate of observation function parameter 1
                   alpha.init=NULL,  #initial estimate of observation function parameter 2
                   p.init=NULL,
                   values=rep(Inf,nsets),   #-log(L) values of all attempts
                   fits=list(),             #additional fit information
                   estimates=matrix(numeric(nsets*cutoff),nrow=nsets), #estimated true outbreaks
                   f, #outbreak size probabilties from simulation (may be list or vector)
                   cutoff=57,
                   obs.fun="logistic"
)#simulated outbreak distributions
{
  sets <- nsets
  if(is.vector(f))  sets <- 1
  if(is.list(f))  sets <- length(f)

  if(sets!=nsets){
    values=rep(Inf,sets)
    estimates=matrix(numeric(sets*cutoff),nrow=sets)
  }
  
  fitsAll <- foreach(i=1:sets,
                    .export=c('getObs','obs.freq','obs.inf')) %do% {
    fit <- list()
    
    l <- sum(f[[i]])
    true.freq <- f[[i]]/l
    if(sets==1) {
      l <- sum(f)
      true.freq <- f/l
    }
    true.inf <- max(1/nsim,sum(true.freq[(cutoff+1):maxsize]))
    true.freq <- true.freq[1:cutoff]
    
    par <- NA
    if(obs.fun=='logistic'){
      par <- getObs(true.freq, obs.freq, beta=beta.init, alpha=alpha.init/100,
                    obs.inf=obs.inf, t.inf=true.inf, cutoff=cutoff, obs.fun=obs.fun)
    }
    if(obs.fun=='geometric'){
      par <- getObs(true.freq, obs.freq, p=p.init,
                    obs.inf=obs.inf, t.inf=true.inf, cutoff=cutoff, obs.fun=obs.fun)
    }
    print(paste('fit ',i,sep=""))
    
    fit <- list(freq=true.freq, length=l,
                par=par, prob=par$prob,
                values=par$val,
                estimates=par$est)
    fit
  }
  
  for(i in 1:sets){
    fits[[i]] <- fitsAll[[i]]$par
    values[i] <- fitsAll[[i]]$values
    estimates[i,] <- fitsAll[[i]]$estimates
  }
  
  return(list(fits=fits, values=values, estimates=estimates))
}
