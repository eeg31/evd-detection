load('results/sim-freqs-FO.Rda')
freqs <- all.freqs

fit.list <- list()
est.list <- list()
j <- 1
for (q in cuts){
    cutoff <- q
    
    obs.dist <- known$cases
    obs.freq <- numeric(cutoff)
    for (i in 1:cutoff){
        obs.freq[i] <- length(which(obs.dist==i))
    }
    obs.inf <- length(which(obs.dist>cutoff))
    obs.dist <- obs.dist[which(obs.dist < cutoff)]
    source('define-model-fit.R')
    
    attempt <- tryFit(beta.init=5,alpha.init=1,f=freqs[1:SAsets],cutoff=q)
    estimates <- attempt$estimates
    fits <- attempt$fits
    
    median <- apply(estimates, 2, median)
    print(paste('cutoff:',q))
    print(paste('all outbreaks est',sum(median)+length(which(obs.dist>cutoff)),
                "(",100*nrow(known)/(sum(median)+length(which(obs.dist>cutoff))),"%)"))
    print(paste('dead end',median[1],"(",100*obs.freq[1]/median[1],"%)"))
    est.list[[j]] <- estimates
    fit.list[[j]] <- fits
    
    j <- j+1
}
save(est.list,file='results/SA.estlist.Rda')
save(fit.list,file='results/SA.fitlist.Rda')
