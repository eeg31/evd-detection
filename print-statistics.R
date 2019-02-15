file <- './results/stats.txt'

for(d in distnames){
  set <- df.all.full[which(df.all.full$dataset==d),]
  sink(file=file,append=TRUE)
  print(paste('Dataset: ',d,"* * * * * * * * * * * * * * * * * * * * * *",sep=""))
  print(paste('number of observed outbreaks: ',allseen))
  print(paste('number of observed dead end spillovers: ',spillseen))
  print(paste('median number of expected outbreaks: ',sum(set$median)+obs.inf))
  print(paste('range of expected outbreaks: ',sum(set$min)+obs.inf,' to ',sum(set$max)+obs.inf,sep=""))
  print(paste('median number of expected dead-end spillovers: ',set$median[which(set$x==1)]))
  print(paste('range of expected dead-end spillovers: ',set$min[which(set$x==1)],' to ',set$max[which(set$x==1)],sep=""))
  sink(file=NULL)
}
