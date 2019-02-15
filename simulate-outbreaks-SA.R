for(decay in decays){

  SA.freqs <- list()
  for (i in 1:length(R0s)){
    freqs.list <- list()
    
    for (k in 1:length(sizes)){
      size <- sizes[k]
      
      ##number of secondary infections caused by index case in each simulation
      children <- rnbinom(nsim, size=size, mu=R0)

      dist.epi.sizes <- rep(1,nsim)
      for (j in 1:nsim){
        R0 <- R0s[i]
        gen <- 1
        
        #if any secondary infections, calculate their offspring and store in tree (igraph)
        if (children[j] > 0){
          edges <- matrix(numeric(2*children[j]),ncol=2)
          edges[,1] <- rep(1,children[j])
          edges[,2] <- 2:(children[j]+1)
          tree <- graph_from_edgelist(edges)
          
          R0 <- R0*(1-decay)
          branches <- rnbinom(children[j], size=size, mu=R0)
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
                for (m in 1:length(kids)){
                  tree <- add_edges(tree, c(thisGen[l], kids[m]))
                }
              } 
            }
            
            if (length(newGen)>0) {
              branches <- rnbinom(length(newGen), size=size, mu=R0)
            } else {
              branches <- 0
            }
            thisGen <- newGen
            gen <- gen + 1
          }
          
          rm(tree)
          #store final epidemic size
          dist.epi.sizes[j] <- ID - 1
          #assume outbreak size is infinite if it exceeds max cases/generations           
          if (ID >= maxsize | gen >= maxgen) dist.epi.sizes[j] <- Inf
        }
        
      }
      #store distribution of epidemic sizes
      freqs <- numeric(maxsize)
      for(j in 1:maxsize){
        freqs[j] <- length(which(dist.epi.sizes==j))
      }
      if(any((is.infinite(dist.epi.sizes)))) freqs[maxsize] <- freqs[maxsize] + length(which(is.infinite(dist.epi.sizes)))
      rm(dist.epi.sizes)    
      
      freqs.list[[k]] <- freqs
      rm(freqs)
      
      #to track progress
      print(paste(round(((i-1)*length(sizes)+k)/(length(R0s)*length(sizes))*100,2),'%',sep=""))
      
      }
      SA.freqs[[i]] <- freqs.list
      
  }
  save(SA.freqs,file=paste('results/sim-freqs-SA-decay',decay,'.Rda',sep=""))
}
