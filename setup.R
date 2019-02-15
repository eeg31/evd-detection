#read in outbreak sizes for reported EVD outbreaks
known <- read.csv('data/detections.csv')
#exclude outbreaks which have been combined, which occur in labs, or are of RESTV
known <- known[which(known$include=='Y'),]

AICc <- function(nLogLik, n, k){
  lik <- exp(-nLogLik)
  AIC <- 2*k - 2*n*log(lik)
  out <- AIC + (2*k^2+2*k)/(n-k-1)
  return(out)
}
