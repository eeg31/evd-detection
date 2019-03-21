load(file='results/fit-SA-all.Rda') 
figs <- list()
d.all <- NULL
for(d in 1:length(decays)){
  
  decay <- decays[d]
  ests <- est.list[[d]]
  est.matrix <- matrix(numeric(length(sizes)*length(R0s)),
                       nrow=length(sizes))
  for(i in 1:length(R0s)){
    est <- ests[[i]]
    overall <- allseen/(obs.inf+rowSums(est$estimates))
    est.matrix[,i] <- overall
  }    
  colnames(est.matrix) <- R0s
  rownames(est.matrix) <- sizes
  est.data <- melt(est.matrix, varnames=c('sizes','R0'), value.name='estimate',as.is=T)
  if(is.null(d.all)) {d.all <- est.data
  } else d.all <- rbind(d.all,est.data)
  
  figs[[d]] <- ggplot(est.data)+
    geom_tile(aes(x=sizes, y=R0, fill=estimate))+
    scale_fill_gradient2(name='prop. observed',
                         high='black', mid=cols[1], low='white',
                         midpoint = 0.5, limits=c(0,1))+ 
    geom_rect(data=data.frame(),
              aes(xmin=0.5, xmax=3.5, ymin=2.5, ymax=5.5), colour='black', alpha=0, lwd=2
    ) + theme_classic() + theme(legend.position="none") +
    scale_x_discrete(name='dispersion (k)')+ggtitle(bquote(delta==.(decay)))
}
leg.plot <- ggplot(est.data)+
  geom_tile(aes(x=sizes, y=R0, fill=estimate))+
  scale_fill_gradient2(name='prop. observed',
                       high='black', mid=cols[1], low='white',
                       midpoint = 0.5, limits=c(0,1))
legend <- get_legend(leg.plot)
rm(leg.plot)

f1 <- plot_grid(figs[[1]],figs[[3]],nrow=2, labels=c('A','C'))
f2 <- plot_grid(figs[[2]],figs[[4]],nrow=2, labels=c('B','D'))

pdf('figures/SA-R0-decay.pdf',width=7.5,height=6.5)
print(plot_grid(f1,f2,legend, nrow=1, rel_widths = c(.4,.4,.2)))
dev.off()


load(file='./results/fit-13sim.Rda')

xmax <- 13

df.Sim.extra <- df.Sim[which(df.Sim$x>xmax),]
df.Sim.extra <- rbind(colSums(df.Sim.extra[1:(cutoff-xmax),]), colSums(df.Sim.extra[(cutoff-xmax+1):nrow(df.Sim.extra),]))
df.Sim.extra <- as.data.frame(rbind(df.Sim.extra,df.Sim.extra[2,],df.Sim.extra[1,]))
df.Sim.extra$x <- c(xmax+0.5,xmax+0.5,xmax+1.5,xmax+1.5)
df.Sim.extra$logCI <- log(df.Sim.extra$CI95)
df.Sim.extra$logCI[which(df.Sim.extra$logCI==-Inf)] <- -1
df.Sim.extra$logIQR <- log(df.Sim.extra$IQR)
df.Sim.extra$logIQR[which(df.Sim.extra$logIQR==-Inf)] <- -1


df.Sim <- df.Sim[which(df.Sim$x<=xmax),]

df.Sim$logCI <- log(df.Sim$CI95)
df.Sim$logCI[which(df.Sim$logCI==-Inf)] <- -1
df.Sim$logIQR <- log(df.Sim$IQR)
df.Sim$logIQR[which(df.Sim$logIQR==-Inf)] <- -1
df.Sim$logSCI <- log(df.Sim$SCI95)
df.Sim$logSCI[which(df.Sim$logSCI==-Inf)] <- -1
df.Sim$logSIQR <- log(df.Sim$SIQR)
df.Sim$logSIQR[which(df.Sim$logSIQR==-Inf)] <- -1
df.Sim$logmedian <- log(df.Sim$median)
df.Sim$logmedian[which(df.Sim$logmedian==-Inf)] <- -1
df.Sim$logSmedian <- log(df.Sim$Smedian)
df.Sim$logSmedian[which(df.Sim$logSmedian==-Inf)] <- -1


logX <- log(obs.freq[1:xmax])
logX[which(logX==-Inf)] <- -1

pdf(file='./figures/simulate13.pdf',width=5,height=4)
  print({
    ggplot(df.Sim,aes(x=x))+
    geom_polygon(aes(y=log(CI95+1)), stat="identity", fill=cols[1],alpha=0.3) +
    geom_polygon(aes(y=log(IQR+1)), stat="identity", fill=cols[1], alpha=0.3) +
    scale_x_continuous(name="outbreak size", limits = c(1,xmax+1.5),breaks=1:(xmax+1),labels = c(1:xmax,paste(xmax+1,cutoff,sep="-"))) +
    scale_y_continuous(name="number of outbreaks", limits = c(0,log(400)), 
                       breaks=log(c(seq(1,11,by=2),seq(21,272,by=50))), 
                       labels=c(seq(1,11,by=2),seq(21,272,by=50))-1) +
    geom_line(aes(y=log(median+1)),col=cols[1]) +
    geom_polygon(aes(y=log(SCI95+1)), stat="identity", fill=rgb(0,0,0,0.4)) +
    geom_polygon(aes(y=log(SIQR+1)), stat="identity", fill=rgb(0,0,0,0.4)) +
    geom_line(aes(y=logSmedian)) +
    geom_point(data=data.frame(x=1:xmax,y=log(obs.freq[1:xmax]+1)), 
               aes(y=y, x=x), stat="identity", shape=4,size=1.5,stroke=1.5) +
    theme_classic()+
    geom_vline(xintercept=xmax,linetype='1F')+
    geom_polygon(data=df.Sim.extra, aes(x=x,y=log(CI95+1)), fill=rgb(0,0,0,0.4))+
    geom_polygon(data=df.Sim.extra, aes(x=x,y=log(IQR+1)), fill=rgb(0,0,0,0.4))+
    geom_point(data=data.frame(x=xmax+1,y=log(allseen+1)), aes(x=x,y=y), shape=4, size=1.5, stroke=1.5)+
    scale_fill_manual(values=c(rgb(0,0,0,0.4),rgb(1,0,0,0.4)),
                      labels=c("True", "Observed"))+
    theme(text = element_text(size=12), legend.position="right")
  })
dev.off()
