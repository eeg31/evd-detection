probs <- matrix(numeric(cutoff*nsets),nrow=nsets)
for (i in 1:nsets){
  prob <- fits[[i]]$prob
  probs[i,] <- prob
}

median <- numeric(cutoff)
min <- numeric(cutoff)
max <- numeric(cutoff)
minQ <- numeric(cutoff)
maxQ <- numeric(cutoff)
medianprob <- numeric(cutoff)
maxprob <- numeric(cutoff)
minprob <- numeric(cutoff)

for(x in 1:cutoff){
  median[x] <- quantile(estimates[,x],0.5)
  min[x] <- min(estimates[,x])
  max[x] <- max(estimates[,x])
  maxQ[x] <- quantile(estimates[,x],0.75)
  minQ[x] <- quantile(estimates[,x],0.25)
  
  
  medianprob[x] <- quantile(probs[,x],0.5) 
  minprob[x] <- min(probs[,x])
  maxprob[x] <- max(probs[,x])
}

probscale <- 260
df.new <- data.frame(x=1:xmax,freq=obs.freq[1:xmax], 
                    median=median[1:xmax],
                    min=min[1:xmax],
                    max=max[1:xmax],
                    minQ=minQ[1:xmax],
                    maxQ=maxQ[1:xmax],
                    probmedian=medianprob[1:xmax]*probscale,
                    probmax=maxprob[1:xmax]*probscale,
                    probmin=minprob[1:xmax]*probscale)
df.new$dataset <- rep(distname,nrow(df.new))

df.new.full <- data.frame(x=1:cutoff,freq=obs.freq[1:cutoff], 
                     median=median[1:cutoff],
                     min=min[1:cutoff],
                     max=max[1:cutoff],
                     minQ=minQ[1:cutoff],
                     maxQ=maxQ[1:cutoff],
                     probmedian=medianprob[1:cutoff]*probscale,
                     probmax=maxprob[1:cutoff]*probscale,
                     probmin=minprob[1:cutoff]*probscale)
df.new.full$dataset <- rep(distname,nrow(df.new.full))

if(is.null(df.all)){
  df.all <- df.new
  df.all.full <- df.new.full
} else{
  df.all <- rbind(df.all,df.new)
  df.all.full <- rbind(df.all.full,df.new.full)
}

