load('results/sim-freqs-FO.Rda')
load('./results/fits-FO.Rda')

simPlot <- data.frame(size=1:cutoff,count=obs.freq,obs=rep(TRUE,cutoff))
observed <- data.frame(size=1:cutoff,count=obs.freq,obs=rep(TRUE,cutoff))

for(i in 1:nToSim){
    parSet <- sample(1:nsets,1)
    
    sim.all <- all.freqs[[parSet]]
    sim <- sim.all[1:cutoff]
    sim <- rep(1:length(sim),sim)
    
    sample <- table(sample(sim,22))
    x <- 1:cutoff
    y <- rep(0, cutoff)
    y[as.numeric(names(sample))] <- sample
    seen <- numeric(cutoff)
    obs <- fits[[parSet]]$prob
    for(j in 1:cutoff){
        seen[j] <- rbinom(1,size=y[j],prob=obs[j])
    }
    
    while(sum(seen) < 13){
        sample <- sample(sim,1)
        y[sample] <- y[sample] + 1
        seen[sample] <- seen[sample] + rbinom(1,size=1,prob=obs[sample])
    }
    
    simPlot <- rbind(simPlot,data.frame(size=x, count=y, obs=rep(FALSE,length(y))))    
    observed <- rbind(observed,data.frame(size=x, count=seen, obs=rep(FALSE,length(y))))
    print(i)
}
Tmedian <- numeric(cutoff) 
TIQRmin <- numeric(cutoff)
T95min <- numeric(cutoff)
TIQRmax <- numeric(cutoff)
T95max <- numeric(cutoff)
Smedian <- numeric(cutoff) 
SIQRmin <- numeric(cutoff)
S95min <- numeric(cutoff)
SIQRmax <- numeric(cutoff)
S95max <- numeric(cutoff)

simPlot2 <- simPlot[which(!simPlot$obs),]
observed <- observed[which(!observed$obs),]

for(x in 1:cutoff){
    T95min[x] <- quantile(simPlot2[which(simPlot2$size==x),]$count,0.025)
    T95max[x] <- quantile(simPlot2[which(simPlot2$size==x),]$count,0.975)
    TIQRmin[x] <-  quantile(simPlot2[which(simPlot2$size==x),]$count,0.25)
    TIQRmax[x] <-  quantile(simPlot2[which(simPlot2$size==x),]$count,0.75)
    Tmedian[x] <- quantile(simPlot2[which(simPlot2$size==x),]$count,0.5)
    
    S95min[x] <- quantile(observed[which(observed$size==x),]$count,0.025)
    S95max[x] <- quantile(observed[which(observed$size==x),]$count,0.975)
    SIQRmin[x] <-  quantile(observed[which(observed$size==x),]$count,0.25)
    SIQRmax[x] <-  quantile(observed[which(observed$size==x),]$count,0.75)
    Smedian[x] <- quantile(observed[which(observed$size==x),]$count,0.5)
}

df.Sim <- data.frame(x=c(1:cutoff,cutoff:1),freq=obs.freq, 
                     IQR=c(TIQRmax[1:cutoff],rev(TIQRmin[1:cutoff])),
                     CI95=c(T95max[1:cutoff],rev(T95min[1:cutoff])),
                     median=c(Tmedian[1:cutoff],rev(Tmedian[1:cutoff])),
                     SIQR=c(SIQRmax[1:cutoff],rev(SIQRmin[1:cutoff])),
                     SCI95=c(S95max[1:cutoff],rev(S95min[1:cutoff])),
                     Smedian=c(Smedian[1:cutoff],rev(Smedian[1:cutoff])))
save(df.Sim,file='./results/fit-13sim.Rda')
