df.all <- NULL
df.all.full <- NULL
df.poly <- NULL
for(d in distnames) {
  distname <- d
  load(file=paste('results/estimates-',d,'.Rda',sep=""))
  load(file=paste('results/fits-',d,'.Rda',sep=""))
  source('compile-results.R')
  
  if(is.null(df.poly)){
    df.poly <- data.frame(x=c(1:xmax,xmax:1),
                          dataset=rep(d,xmax*2),
                          estimates=c(df.new$min,rev(df.new$max)),
                          probs=c(df.new$probmin,rev(df.new$probmax)))
  } else {
    df.poly <- rbind(df.poly,
                     data.frame(x=c(1:xmax,xmax:1),
                                dataset=rep(d,xmax*2),
                                estimates=c(df.new$min,rev(df.new$max)),
                                probs=c(df.new$probmin,rev(df.new$probmax))))
  }
}

cols <- brewer_pal(type='qual', palette='Dark2')(3)

largeseen <- allseen - sum(obs.freq)

figC <- ggplot(data=subset(df.all,dataset=="FO")) + 
  scale_x_continuous(name="outbreak size", expand=c(0.01,0.01), limits = c(1,xmax+1.5)) +
  geom_polygon(data=subset(df.poly,dataset=="FO"), aes(group=dataset, x=x, y=estimates), 
               alpha=0.5, fill=cols[1])+
  geom_polygon(data=subset(df.poly,dataset=="FO"), aes(group=dataset, x=x, y=probs), 
               alpha=0.3)+
  geom_line(aes(x=x,y=median,group=dataset),
            color=cols[1], lwd=1.1) +
  geom_line(aes(x=x,y=probmedian,group=dataset),
            lwd=1.1, linetype='dashed') +
  geom_point(data=data.frame(x=1:xmax,y=obs.freq[1:xmax]), aes(y=y, x=x), stat="identity", shape=4,size=2) +
  geom_bar(data=data.frame(x=xmax+1,y=largeseen), aes(x=x,y=y), width=1,stat='identity')+
  scale_y_continuous(name="", expand=c(0.01,0.01),
                     sec.axis=sec_axis(~./probscale, name="")) +
  theme_classic()+  theme_classic() +theme(axis.text = element_text(size=16),
                                           axis.title =  element_text(size=16))

figD <- ggplot(data=subset(df.all,dataset=="SL")) + 
  scale_x_continuous(name="outbreak size", limits = c(1,xmax+1.5), expand=c(0.01,0.01)) +
  geom_polygon(data=subset(df.poly,dataset=="SL"), aes(group=dataset, x=x, y=estimates), 
               alpha=0.6, fill=cols[2])+
  geom_polygon(data=subset(df.poly,dataset=="SL"), aes(group=dataset, x=x, y=probs), 
               alpha=0.3)+
  geom_line(aes(x=x,y=median,group=dataset),
            color=cols[2],lwd=1.1) +
  geom_line(aes(x=x,y=probmedian,group=dataset),
            lwd=1.1, linetype='dashed') +
  geom_point(data=data.frame(x=1:xmax,y=obs.freq[1:xmax]), aes(y=y, x=x), stat="identity", shape=4,size=2) +
  geom_bar(data=data.frame(x=xmax+1,y=largeseen), aes(x=x,y=y), width=1,stat='identity')+
  scale_y_continuous(name="number of outbreaks", expand=c(0.01,0.01),
                     sec.axis=sec_axis(~./probscale, name="probability of observation")) +
  theme_classic()+  theme_classic() +theme(axis.text = element_text(size=16), 
                                           axis.title =  element_text(size=16))

figE <- ggplot(data=subset(df.all,dataset=="Guin")) + 
  scale_x_continuous(name="outbreak size", limits = c(1,xmax+1.5),expand=c(0.01,0.01)) +
  geom_polygon(data=subset(df.poly,dataset=="Guin"), aes(group=dataset, x=x, y=estimates), 
               alpha=0.6, fill=cols[3])+
  geom_polygon(data=subset(df.poly,dataset=="Guin"), aes(group=dataset, x=x, y=probs), 
               alpha=0.3)+
  geom_line(aes(x=x,y=median,group=dataset),
            color=cols[3],lwd=1.1) +
  geom_line(aes(x=x,y=probmedian,group=dataset),
            lwd=1.1,linetype='dashed') +
  geom_point(data=data.frame(x=1:xmax,y=obs.freq[1:xmax]), aes(y=y, x=x), stat="identity", shape=4,size=2) +
  geom_bar(data=data.frame(x=xmax+1,y=largeseen), aes(x=x,y=y), width=1,stat='identity')+
  scale_y_continuous(name="", expand=c(0.01,0.01),
                     sec.axis=sec_axis(~./probscale, name="")) +
  theme_classic() +theme(axis.text = element_text(size=16),
                         axis.title =  element_text(size=16))

sums <- data.frame(distname = c('FO','SL','Guin'),
                   location=c('West Africa','Sierra Leone','Guinea'))
sums$type <- 'all sizes'
sums$med <- NA
sums$min <- NA
sums$max <- NA
sums$minQ <- NA
sums$maxQ <- NA
for(i in 1:nrow(sums)){
  d <- sums$distname[i]
  sums$med[i] <- sum(subset(df.all, dataset==d)$median) + obs.inf
  sums$min[i] <- sum(subset(df.all, dataset==d)$min) + obs.inf
  sums$max[i] <- sum(subset(df.all, dataset==d)$max) + obs.inf
  sums$minQ[i] <- sum(subset(df.all, dataset==d)$minQ) + obs.inf
  sums$maxQ[i] <- sum(subset(df.all, dataset==d)$maxQ) + obs.inf
}
spills <- sums
spills$type <- 'single case'
for(i in 1:nrow(sums)){
  d <- spills$distname[i]
  spills$med[i] <- subset(df.all, dataset==d)$median[1]
  spills$min[i] <- subset(df.all, dataset==d)$min[1]
  spills$max[i] <- subset(df.all, dataset==d)$max[1]
  spills$minQ[i] <- subset(df.all, dataset==d)$minQ[1]
  spills$maxQ[i] <- subset(df.all, dataset==d)$maxQ[1]
}
sums <- rbind(sums,spills)
rm(spills)

figB <- ggplot(sums)+
  geom_boxplot(aes(x=location, group=type, middle=med, lower=minQ, upper=maxQ, ymin=min,ymax=max,
                   color=location, fill=location, alpha=type),
               stat='identity',lwd=1.1)+
  geom_hline(yintercept=allseen, lty='dashed')+
  geom_hline(yintercept=spillseen, lty='dotted')+
  scale_y_continuous(name='number of outbreaks',limits=c(0,300),expand=c(0,0))+
  scale_x_discrete(name="dataset")+
  scale_alpha_manual(values=c(0.5,0))+ 
  scale_fill_manual(values=cols[3:1])+
  scale_color_manual(values=cols[3:1])+
  theme(legend.position = c(.2,.8), legend.text = element_text(size=16),
        legend.title = element_text(size=16), axis.text = element_text(size=16),
        axis.title =  element_text(size=16))+
  guides(alpha = guide_legend(override.aes = list(fill = 'black')), fill=FALSE, color=FALSE)



offsprings <- NULL
for(i in 1:length(distnames)){
  d <- distnames[i]
  load(eval(paste('results/R0s-',d,'.Rda',sep="")))
  load(eval(paste('results/sizes-',d,'.Rda',sep="")))
  
  min <- numeric(xmax)
  max <- numeric(xmax)
  med <- numeric(xmax)
  
  for (j in 1:xmax){
    min[j] <- min(dnbinom(j, size=sizes, mu=R0s))
    max[j] <- max(dnbinom(j, size=sizes, mu=R0s))
    med[j] <- median(dnbinom(j, size=sizes, mu=R0s))
  }
  
  new.offs <- data.frame(dataset=rep(d,xmax*2),
                         x=c(1:xmax,xmax:1),
                         med=c(med,rev(med)),
                         bounds=c(min,rev(max)))
  if(is.null(offsprings)) {
    offsprings <- new.offs
  } else {
    offsprings <- rbind(offsprings, new.offs)
  }
}

figA <- ggplot(offsprings) +
  geom_polygon(aes(x=x,y=bounds,fill=dataset, color=dataset, group=dataset), alpha=0.5)+
  scale_fill_manual(values=c('FO'=cols[1],
                            'SL'=cols[2],
                            'Guin'=cols[3]),
                    labels=c('West Africa', 'Sierra Leone','Guinea'))+
  scale_color_manual(values=c('FO'=cols[1],
                              'SL'=cols[2],
                              'Guin'=cols[3]),
                     labels=c('West Africa', 'Sierra Leone','Guinea'))+
  scale_y_continuous('probability',expand=c(0,0))+
  scale_x_continuous('number of secondary infections',expand=c(0.01,0.01)) +
  theme(legend.position = c(.6,.6), legend.text = element_text(size=16),
        legend.title = element_text(size=16), axis.text = element_text(size=16),
        axis.title =  element_text(size=16))
library(cowplot)

p1 <- plot_grid(figA,figB,ncol=1,labels = c('A','B'), label_size=18)
p2 <- plot_grid(figC,figD,figE,ncol=1,labels=c('C','D','E'), label_size=18, label_x=0.15)
pF <- plot_grid(p1,p2,ncol=2)

pdf(file='figures/main-results.pdf',width=12,height=8)
print(pF)
dev.off()


