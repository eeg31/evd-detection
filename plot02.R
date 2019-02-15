AIC <- data.frame(dataset=NULL, type=NULL, AICc=NULL)
alphas <- numeric(nsets)
betas <- numeric(nsets)
ps <- numeric(nsets)
for(d in distnames){
  load(file=paste('results/fits-geom-',d,'.Rda',sep=""))
  new.data <- data.frame(dataset=rep(d,nsets), type=rep('geometric',nsets), AICc=numeric(nsets))
  for(i in 1:nsets){
    new.data$AICc[i] <- AICc(fits[[i]]$val, n=allseen, k=1)
    if(d=='FO') ps[i] <- fits[[i]]$p
  }
  
  AIC <- rbind(AIC,new.data)
  load(file=paste('results/fits-',d,'.Rda',sep=""))
  new.data <- data.frame(dataset=rep(d,nsets), type=rep('logistic',nsets), AICc=numeric(nsets))
  for(i in 1:nsets){
    new.data$AICc[i] <- AICc(fits[[i]]$val, n=allseen, k=2)
    if(d=='FO') {
      alphas[i] <- fits[[i]]$alpha
      betas[i] <- fits[[i]]$beta
    }
    
  }  
  AIC <- rbind(AIC,new.data)
}

geometricObs <- function(p,i){
  out <- 1-(1-p)^i
  return(out)
}
convertToLogistic <- function(alpha,beta,i){
  p <- 1-(1-(1+exp(beta-i))^(-alpha))^(1/i)
  return(p)
}

obsXmax <- 20
obsCompare <- data.frame(x=NULL, type=NULL, func=NULL, val=NULL)
j <- 1
for (i in 1:nsets){
  toAdd <- data.frame(x=rep(1:obsXmax, 4),
                      type=c(rep('individual', 2*obsXmax),
                             rep('cumulative', 2*obsXmax)),
                      func=c(rep('geometric',obsXmax),
                             rep('logistic',obsXmax), 
                             rep('geometric',obsXmax),
                             rep('logistic',obsXmax)),
                      val=c(rep(ps[i],obsXmax),
                            convertToLogistic(alphas[i],betas[i],1:obsXmax),
                            geometricObs(ps[i],1:obsXmax),
                            probs[i,1:obsXmax]))
  obsCompare <- rbind(obsCompare,toAdd)
                      
}

figInd <- ggplot(subset(obsCompare,func=='logistic'&type=='individual')) +
  geom_boxplot(aes(x=x,y=val,group=x,color=func,fill=func),alpha=0.5)+
  geom_boxplot(data=subset(obsCompare,func=='geometric'&type=='individual'),
               aes(x=x,y=val,group=x,color=func,fill=func),alpha=0.5)+
  scale_fill_manual(values=cols[1:2],name="observation\nfunction")+
  scale_color_manual(values=cols[1:2],name="observation\nfunction")+
  scale_y_continuous(name="individual probability of detection")+
  scale_x_continuous(name="outbreak/cluster size")

figA <- ggplot(subset(obsCompare,func=='logistic'&type=='cumulative')) +
  geom_boxplot(aes(x=x,y=val,group=x,color=func,fill=func),alpha=0.5)+
  geom_boxplot(data=subset(obsCompare,func=='geometric'&type=='cumulative'),
               aes(x=x,y=val,group=x,color=func,fill=func),alpha=0.5)+
  scale_fill_manual(values=cols[1:2],name="observation\nfunction")+
  scale_color_manual(values=cols[1:2],name="observation\nfunction")+
  scale_y_continuous(name="probability of cluster detection")+
  scale_x_continuous(name="outbreak/cluster size")+
  guides(color=FALSE,alpha = guide_legend(override.aes = list(fill = 'black')))


figB <- ggplot() +
        geom_boxplot(aes(x=dataset,fill=type,y=AICc), 
                     data=subset(AIC), position=position_dodge())+
        scale_x_discrete(breaks=c("FO","SL","Guin"),
                         labels=c("West Africa", "Sierra Leone", "Guinea"))+
        guides(fill=FALSE) +
        scale_fill_manual(values=cols[1:2])
pAIC <- plot_grid(figA,figB,labels=c('A','B'),ncol=2,rel_widths = c(0.6,0.4))

pdf('figures/individual-prob.pdf',width=6,height=4)
print(figInd)
dev.off()

pdf('figures/AIC-compare.pdf',width=12,height=4)
print(pAIC)
dev.off()
