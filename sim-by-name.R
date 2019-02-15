all.freqs <- list()

for (i in 1:nsets){
  #number of secondary infections caused by index case in each simulation
  children <- rnbinom(nsim, mu=R0s[i], size=sizes[i])
    
  dist.epi.sizes <- rep(1,nsim)
  for (j in 1:nsim){
    #if any secondary infections, calculate their offspring and store in tree (igraph)
     if (children[j] > 0){
          edges <- matrix(numeric(2*children[j]),ncol=2)
          edges[,1] <- rep(1,children[j])
          edges[,2] <- 2:(children[j]+1)
          tree <- graph_from_edgelist(edges)
            
          branches <- rnbinom(children[j], mu=R0s[i], size=sizes[i])
          thisGen <- 2:(children[j]+1)
          ID <- children[j] + 2
          gen <- 2
          #continue adding leaves and gens until epidemic dies out or reaches time limit
          while(gen <= maxgen & sum(branches) > 0 & ID<maxsize){
              newGen <- NULL
              for(l in 1:length(thisGen)){
                 if(branches[l] > 0) {
                      kids <- ID:(ID + branches[l] - 1)
                      ID <- ID + branches[l]
                      newGen <- c(newGen, kids)
                      
                      tree <- add_vertices(tree, nv=length(kids), names=kids)
                      for (k in 1:length(kids)){
                      tree <- add_edges(tree, c(thisGen[l], kids[k]))
                      }
                  } 
              }
              
              if (length(newGen)>0) {
                  branches <- rnbinom(length(newGen), mu=R0s[i], size=sizes[i])
              } else {
                 branches <- 0
              }
              thisGen <- newGen
              gen <- gen + 1
            }
          #store final epidemic size
          dist.epi.sizes[j] <- ID - 1
          #assume outbreak size is infinite if it exceeds max cases/generations
          if (ID >= maxsize | gen >= maxgen) dist.epi.sizes[j] <- Inf
      }
  }
    
  freqs <- numeric(maxsize)
  for(j in 1:maxsize){
    freqs[j] <- length(which(dist.epi.sizes==j))
  }
  freqs[maxsize] <- freqs[maxsize] + length(which(is.infinite(dist.epi.sizes)))
  
  #store distribution of epidemic sizes
  all.freqs[[i]] <- freqs
  print(paste(i*nsim,' of ',nsets*nsim,' (',round(i/nsets*100,2),'%) simulations complete',sep=""))
}
save(all.freqs, file=paste('results/sim-freqs-',distname,'.Rda',sep=""))
