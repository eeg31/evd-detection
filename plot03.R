SA <- data.frame(cutoff=c(cuts,rev(cuts)))
j <- 1

obs.dist <- known$cases
obs.freq <- numeric(max(known$cases))
for (i in 1:max(known$cases)){
  obs.freq[i] <- length(which(obs.dist==i))
}
largeseen <- numeric(length(cuts))
for(i in 1:length(cuts)){
  largeseen[i] <- length(which(obs.dist>cuts[i]))
}


median.outbreaks=numeric(length(cuts))
median.dead=numeric(length(cuts))
prop.all.outbreaks=numeric(length(cuts))
prop.DE.outbreaks=numeric(length(cuts))
median.detect.rate=numeric(length(cuts))
prop.DE.outbreaks.min95=numeric(length(cuts))
prop.DE.outbreaks.max95=numeric(length(cuts))
prop.all.outbreaks.min95=numeric(length(cuts))
prop.all.outbreaks.max95=numeric(length(cuts))
prop.DE.outbreaks.minIQR=numeric(length(cuts))
prop.DE.outbreaks.maxIQR=numeric(length(cuts))
prop.all.outbreaks.minIQR=numeric(length(cuts))
prop.all.outbreaks.maxIQR=numeric(length(cuts))

for (q in cuts){
  estimates <- est.list[[j]]
  median <- apply(estimates, 2, median)
  max95 <- apply(estimates, 2, quantile, probs=0.975)
  min95 <- apply(estimates, 2, quantile, probs=0.025)
  maxIQR <- apply(estimates, 2, quantile, probs=0.75)
  minIQR <- apply(estimates, 2, quantile, probs=0.25)
  
  median.outbreaks[j] <- sum(median) + largeseen[j]
  median.detect.rate[j] <- (sum(obs.freq[1:q])+ largeseen[j])/(sum(median)+ largeseen[j])
  median.dead[j] <- sum(obs.freq[1])/median[1]
  prop.DE.outbreaks[j] <- obs.freq[1]/median[1]
  prop.all.outbreaks[j] <- nrow(known)/(sum(median[1:q])+largeseen[j])
  prop.DE.outbreaks.min95[j] <- obs.freq[1]/max95[1]
  prop.all.outbreaks.min95[j] <- nrow(known)/(sum(max95[1:q])+largeseen[j])
  prop.DE.outbreaks.minIQR[j] <- obs.freq[1]/maxIQR[1]
  prop.all.outbreaks.minIQR[j] <- nrow(known)/(sum(maxIQR[1:q])+largeseen[j])  
  prop.DE.outbreaks.max95[j] <- obs.freq[1]/min95[1]
  prop.all.outbreaks.max95[j] <- nrow(known)/(sum(min95[1:q])+largeseen[j]) 
  prop.DE.outbreaks.maxIQR[j] <- obs.freq[1]/minIQR[1]
  prop.all.outbreaks.maxIQR[j] <- nrow(known)/(sum(minIQR[1:q])+largeseen[j])
  j <- j+1
}

SA$median.detect.rate <- c(median.detect.rate,rev(median.detect.rate))
SA$median.dead <- c(median.dead,rev(median.dead))
SA$median.outbreaks <- c(median.outbreaks,rev(median.outbreaks))
SA$DE.IQR <- c(prop.DE.outbreaks.maxIQR,rev(prop.DE.outbreaks.minIQR))
SA$all.IQR <- c(prop.all.outbreaks.maxIQR,rev(prop.all.outbreaks.minIQR))
SA$DE.CI95 <- c(prop.DE.outbreaks.max95,rev(prop.DE.outbreaks.min95))
SA$all.CI95 <- c(prop.all.outbreaks.max95,rev(prop.all.outbreaks.min95))

pdf(file='./figures/cutoff-sensitivity.pdf',width=5,height=4)
print({
  ggplot(SA,aes(x=cutoff)) +
  geom_polygon(aes(y=all.CI95), stat="identity", fill=rgb(.8,.8,.8)) +
  geom_polygon(aes(y=all.IQR), stat="identity", fill=rgb(.5,.5,.5)) +
  scale_x_continuous(name="cutoff for stochastic outbreaks",limits=c(0,55),expand=c(0,0)) +
  geom_line(aes(y=median.detect.rate,x=cutoff)) +
  scale_y_continuous(name="proportion of all outbreaks reported",limits=c(0,.4),expand=c(0,0)) +
  theme_classic()+
  geom_point(data=data.frame(x=known$cases,y=rep(0.05,length(known$cases))),aes(x=x,y=y),
             pch=4,stroke=1.5)
})
dev.off()
