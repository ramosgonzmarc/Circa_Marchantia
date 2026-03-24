#Function for drawing the expression profile of a gene through LD and SD
plot.ld.sd <- function(gene.id, gene.name, genes.expression)
{
  current.gene.expression.ld <- 0
  current.gene.expression.sd <- 0
  
  ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.ld <- genes.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"), paste(ld.zt,3,sep="_"))]
  
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.sd <- genes.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"), paste(sd.zt,3,sep="_"))]
  
  max.expr <- max(c(current.gene.expression.ld, current.gene.expression.sd))
  min.expr <- min(c(current.gene.expression.ld, current.gene.expression.sd))
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="FPKM",
       ylim=c((min.expr-(max.expr- min.expr)/5), max.expr),xlim=c(0,72),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  axis(side=2,lwd=3)
  zt.points<-paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4)),sep="_")
  zt.points.rep<-c(paste(zt.points,1,sep="_"),paste(zt.points,2,sep="_"), paste(zt.points,3,sep="_"))
  axis(side = 1,lwd = 3, at=seq(from=0,to=68,by=4),labels = zt.points.rep,las=2)
  
  lines(x = seq(from=0,by=4,to=68),current.gene.expression.ld,type="o",lwd=3,col="blue")
  lines(x = seq(from=0,by=4,to=68),current.gene.expression.sd,type="o",lwd=3,col="red")
  
  
  for(i in 0:2)
  {
    current.line <- (max.expr- min.expr)/20
    
    polygon(x = c(24*i, 24*i+16, 24*i+16, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(3/2*current.line), min.expr-(3/2*current.line)),
            lwd=2,border="blue")
    polygon(x = c(24*i+16,24*(i+1),24*(i+1),24*i+16),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(3/2*current.line), min.expr-(3/2*current.line)),
            lwd=2,border="blue",col="blue")
    current.line <- 2*current.line
    
    polygon(x = c(24*i, 24*i+8, 24*i+8, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red")
    polygon(x = c(24*i+8,24*(i+1),24*(i+1),24*i+8),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red",col="red")
  }
  
  return(0)  
}

# Function for plotting the expression profile of a gene through SD and LL conditions
plot.sd.ll <- function(gene.id, gene.name, genes.expression)
{
  current.gene.expression.sd <- 0
  
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  sd.zt.i <- paste(rep(sd.zt,3), rep(1:3, each=6), sep="_")
  ll.sd.zt <- paste("ll_sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  ll.sd.zt.i <- paste(rep(ll.sd.zt,2), rep(1:2, each=6), sep="_")
  
  current.gene.expression.sd <- genes.expression[gene.id,c(sd.zt.i, ll.sd.zt.i)]
  
  max.expr <- max(current.gene.expression.sd)
  min.expr <- min(current.gene.expression.sd)
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="",
       ylim=c((min.expr-(max.expr- min.expr)/5), max.expr),xlim=c(0,120),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  abline(v=72, lwd=2, col= "green")
  axis(side=2,lwd=3)
  zt.points<-paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4)),sep="_")
  zt.points.rep<-c(paste(zt.points,1,sep="_"),paste(zt.points,2,sep="_"),paste(zt.points,3,sep="_"), paste(zt.points,4,sep="_"), paste(zt.points,5,sep="_"))
  axis(side = 1,lwd = 3, at=seq(from=0,to=116,by=4),labels = zt.points.rep,las=2)
  
  lines(x = seq(from=0,by=4,to=116),current.gene.expression.sd,type="o",lwd=3,col="red")
  
  
  for(i in 0:2)
  {
    current.line <- (max.expr- min.expr)/10
    
    polygon(x = c(24*i, 24*i+8, 24*i+8, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red")
    polygon(x = c(24*i+8,24*(i+1),24*(i+1),24*i+8),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red",col="red")
  }
  
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(72, 96, 96, 72),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red")
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(96, 120, 120, 96),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red")
  
  return(0)  
}

# Function for plotting the expression profile of a gene through SD and DD conditions
plot.sd.dd <- function(gene.id, gene.name, genes.expression)
{
  current.gene.expression.sd <- 0
  
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  sd.zt.i <- paste(rep(sd.zt,3), rep(1:3, each=6), sep="_")
  dd.sd.zt <- paste("dd_sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  dd.sd.zt.i <- paste(rep(dd.sd.zt,2), rep(1:2, each=6), sep="_")
  
  current.gene.expression.sd <- genes.expression[gene.id,c(sd.zt.i, dd.sd.zt.i)]
  
  max.expr <- max(current.gene.expression.sd)
  min.expr <- min(current.gene.expression.sd)
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="",
       ylim=c((min.expr-(max.expr- min.expr)/5), max.expr),xlim=c(0,120),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  abline(v=72, lwd=2, col= "green")
  axis(side=2,lwd=3)
  zt.points<-paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4)),sep="_")
  zt.points.rep<-c(paste(zt.points,1,sep="_"),paste(zt.points,2,sep="_"),paste(zt.points,3,sep="_"), paste(zt.points,4,sep="_"), paste(zt.points,5,sep="_"))
  axis(side = 1,lwd = 3, at=seq(from=0,to=116,by=4),labels = zt.points.rep,las=2)
  
  lines(x = seq(from=0,by=4,to=116),current.gene.expression.sd,type="o",lwd=3,col="red")
  
  
  for(i in 0:2)
  {
    current.line <- (max.expr- min.expr)/10
    
    polygon(x = c(24*i, 24*i+8, 24*i+8, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red")
    polygon(x = c(24*i+8,24*(i+1),24*(i+1),24*i+8),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red",col="red")
  }
  
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(72, 96, 96, 72),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red", col="red")
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(96, 120, 120, 96),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red", col="red")
  
  return(0)  
}

# Function for plotting the expression profile of a gene through LD and LL conditions
plot.ld.ll <- function(gene.id, gene.name, genes.expression)
{
  current.gene.expression.ld <- 0
  
  ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  ld.zt.i <- paste(rep(ld.zt,3), rep(1:3, each=6), sep="_")
  ll.ld.zt <- paste("ll_ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  ll.ld.zt.i <- paste(rep(ll.ld.zt,2), rep(1:2, each=6), sep="_")
  
  current.gene.expression.ld <- genes.expression[gene.id,c(ld.zt.i, ll.ld.zt.i)]
  
  max.expr <- max(current.gene.expression.ld)
  min.expr <- min(current.gene.expression.ld)
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="",
       ylim=c((min.expr-(max.expr- min.expr)/5), max.expr),xlim=c(0,120),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  abline(v=72, lwd=2, col= "green")
  axis(side=2,lwd=3)
  zt.points<-paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4)),sep="_")
  zt.points.rep<-c(paste(zt.points,1,sep="_"),paste(zt.points,2,sep="_"),paste(zt.points,3,sep="_"), paste(zt.points,4,sep="_"), paste(zt.points,5,sep="_"))
  axis(side = 1,lwd = 3, at=seq(from=0,to=116,by=4),labels = zt.points.rep,las=2)
  
  lines(x = seq(from=0,by=4,to=116),current.gene.expression.ld,type="o",lwd=3,col="blue")
  
  
  for(i in 0:2)
  {
    current.line <- (max.expr- min.expr)/10
    
    polygon(x = c(24*i, 24*i+16, 24*i+16, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="blue")
    polygon(x = c(24*i+16,24*(i+1),24*(i+1),24*i+16),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="blue",col="blue")
  }
  
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(72, 96, 96, 72),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="blue")
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(96, 120, 120, 96),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="blue")
  
  return(0)  
}

# Function for plotting the expression profile of a gene through LD and DD conditions
plot.ld.dd <- function(gene.id, gene.name, genes.expression)
{
  current.gene.expression.ld <- 0
  
  ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  ld.zt.i <- paste(rep(ld.zt,3), rep(1:3, each=6), sep="_")
  dd.ld.zt <- paste("dd_ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  dd.ld.zt.i <- paste(rep(dd.ld.zt,2), rep(1:2, each=6), sep="_")
  
  current.gene.expression.ld <- genes.expression[gene.id,c(ld.zt.i, dd.ld.zt.i)]
  
  max.expr <- max(current.gene.expression.ld)
  min.expr <- min(current.gene.expression.ld)
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="",
       ylim=c((min.expr-(max.expr- min.expr)/5), max.expr),xlim=c(0,120),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  abline(v=72, lwd=2, col= "green")
  axis(side=2,lwd=3)
  zt.points<-paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4)),sep="_")
  zt.points.rep<-c(paste(zt.points,1,sep="_"),paste(zt.points,2,sep="_"),paste(zt.points,3,sep="_"), paste(zt.points,4,sep="_"), paste(zt.points,5,sep="_"))
  axis(side = 1,lwd = 3, at=seq(from=0,to=116,by=4),labels = zt.points.rep,las=2)
  
  lines(x = seq(from=0,by=4,to=116),current.gene.expression.ld,type="o",lwd=3,col="blue")
  
  
  for(i in 0:2)
  {
    current.line <- (max.expr- min.expr)/10
    
    polygon(x = c(24*i, 24*i+16, 24*i+16, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="blue")
    polygon(x = c(24*i+16,24*(i+1),24*(i+1),24*i+16),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="blue",col="blue")
  }
  
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(72, 96, 96, 72),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="blue", col="blue")
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(96, 120, 120, 96),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="blue", col="blue")
  
  return(0)  
}

#Function to determine the expression peak and troughs of a gene
classify.genes <- function(rain.results, gene.set)
{
  circadian.rain.results <- rain.results[gene.set,]
  
  gene.peaks <- circadian.rain.results$phase - 4
  gene.troughs <- (gene.peaks + circadian.rain.results$peak.shape) %% 24
  
  return(data.frame(gene=gene.set,peak=gene.peaks,trough=gene.troughs,stringsAsFactors = F))
}

#Function to make clusters according to the ZT at which the genes show their peak expression
#and another dividing those previous clusters by the troughs of the genes
clusters <- function(rain.results, gene.set, photoperiod)
{
  classified.genes <- classify.genes(rain.results,gene.set)
  
  peaks <- seq(from=0,to=20,by=4)
  
  for(i in 1:length(peaks))
  {
    file.name <- paste(c(photoperiod,"peak_clusters/cluster_ZT",peaks[i],".txt"),collapse = "")
    write.table(x = subset(classified.genes, peak == peaks[i])$gene,file = file.name,quote = F,row.names = F,col.names = F)
  }
  
  peaks <- troughs <- seq(from=0,to=20,by=4)
  
  indices <- 1:length(peaks)
  
  for(i in indices)
  {
    current.peak <- peaks[i]
    for (j in indices[-i])
    {
      current.trough <- troughs[j]
      file.name <- paste(c(photoperiod,"peak_trough_clusters/cluster_peak_ZT",current.peak,"_trough_ZT",
                           current.trough,".txt"),collapse="")
      write.table(x = subset(classified.genes, peak == current.peak & trough == current.trough)$gene,
                  file = file.name,quote = F,row.names = F,col.names = F)
    }
  }
}

scale.and.interpolate.mean.expression <- function(mean.expression)
{
  time.points <- seq(from=0,to=24,4)
  sorted.genes <- rownames(mean.expression)
  
  scaled.gene.expression <- matrix(nrow=nrow(mean.expression),ncol=49)
  rownames(scaled.gene.expression) <- sorted.genes
  colnames(scaled.gene.expression) <- paste("ZT",seq(from=0,to=24,by=0.5),sep="_")
  
  i <- 1
  for(i in 1:length(sorted.genes))
  {
    current.gene <- sorted.genes[i]
    
    expression.current.gene.i <- mean.expression[current.gene,]
    #expression.current.gene.i <- short.day.mean.expression[current.gene,]
    
    expression.current.gene.i <- c(expression.current.gene.i, expression.current.gene.i[1])
    #scaled.gene.expression[i,] <- approx(x=time.points,expression.current.gene.i,xout=0:24)[[2]]
    
    scaled.gene.expression[i,] <- approx(x=time.points,expression.current.gene.i,xout=seq(from=0,to=24,by=0.5))[[2]]
    scaled.gene.expression[i,] <- scaled.gene.expression[i,]/max(scaled.gene.expression[i,])
  }
  
  return(scaled.gene.expression)
}

plot.ld.sd.wave <- function(gene.id, gene.name, genes.expression, period=1,
                            wave_parameters)
{
  current.gene.expression.ld <- 0
  current.gene.expression.sd <- 0
  
  ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.ld <- genes.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"), paste(ld.zt,3,sep="_"))]
  
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.sd <- genes.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"), paste(sd.zt,3,sep="_"))]
  
  max.expr <- max(c(current.gene.expression.ld, current.gene.expression.sd))
  min.expr <- min(c(current.gene.expression.ld, current.gene.expression.sd))
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="FPKM",
       ylim=c((min.expr-(max.expr- min.expr)/5), max.expr),xlim=c(0,72),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  axis(side=2,lwd=3)
  zt.points<-paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4)),sep="_")
  zt.points.rep<-c(paste(zt.points,1,sep="_"),paste(zt.points,2,sep="_"), paste(zt.points,3,sep="_"))
  axis(side = 1,lwd = 3, at=seq(from=0,to=68,by=4),labels = zt.points.rep,las=2)
  
  points(x = seq(from=0,by=4,to=68),current.gene.expression.ld,lwd=3,col="blue", pch=16)
  points(x = seq(from=0,by=4,to=68),current.gene.expression.sd,lwd=3,col="red", pch=16)
  
  wave.ld <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["ld.mesor.estimate"]+wave_parameters["ld.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["ld.peak.time.hours"]*pi/12)))
  lines(x = seq(from=0,by=0.1,to=68),wave.ld,lwd=3,col="blue")
  
  wave.sd <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["sd.mesor.estimate"]+wave_parameters["sd.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["sd.peak.time.hours"]*pi/12)))
  lines(x = seq(from=0,by=0.1,to=68),wave.sd,lwd=3,col="red")
  
  
  for(i in 0:2)
  {
    current.line <- (max.expr- min.expr)/20
    
    polygon(x = c(24*i, 24*i+16, 24*i+16, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(3/2*current.line), min.expr-(3/2*current.line)),
            lwd=2,border="blue")
    polygon(x = c(24*i+16,24*(i+1),24*(i+1),24*i+16),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(3/2*current.line), min.expr-(3/2*current.line)),
            lwd=2,border="blue",col="blue")
    current.line <- 2*current.line
    
    polygon(x = c(24*i, 24*i+8, 24*i+8, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red")
    polygon(x = c(24*i+8,24*(i+1),24*(i+1),24*i+8),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red",col="red")
  }
  
  return(0)  
}

plot.sd.ll.wave <- function(gene.id, gene.name, genes.expression, period=1, wave_parameters)
{
  current.gene.expression.sd <- 0
  
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  sd.zt.i <- paste(rep(sd.zt,3), rep(1:3, each=6), sep="_")
  ll.sd.zt <- paste("ll_sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  ll.sd.zt.i <- paste(rep(ll.sd.zt,2), rep(1:2, each=6), sep="_")
  
  current.gene.expression.sd <- genes.expression[gene.id,c(sd.zt.i)]
  current.gene.expression.ll <- genes.expression[gene.id,c(ll.sd.zt.i)]
  
  max.expr <- max(c(current.gene.expression.sd, current.gene.expression.ll))
  min.expr <- min(c(current.gene.expression.sd, current.gene.expression.ll))
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="",
       ylim=c((min.expr-(max.expr- min.expr)/5), max.expr),xlim=c(0,120),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  #abline(v=72, lwd=2, col= "green")
  axis(side=2,lwd=3)
  zt.points<-paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4)),sep="_")
  zt.points.rep<-c(paste(zt.points,1,sep="_"),paste(zt.points,2,sep="_"),paste(zt.points,3,sep="_"), paste(zt.points,4,sep="_"), paste(zt.points,5,sep="_"))
  axis(side = 1,lwd = 3, at=seq(from=0,to=116,by=4),labels = zt.points.rep,las=2)
  
  #lines(x = seq(from=0,by=4,to=116),current.gene.expression.sd,type="o",lwd=3,col="red")
  points(x = seq(from=0,by=4,to=68),current.gene.expression.sd,lwd=3,col="firebrick2", pch=16)
  points(x = seq(from=72,by=4,to=116),current.gene.expression.ll,lwd=3,col="salmon", pch=16)
  
  wave.sd <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["sd.mesor.estimate"]+wave_parameters["sd.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["sd.peak.time.hours"]*pi/12)))
  lines(x = seq(from=0,by=0.1,to=68),wave.sd,lwd=3,col="firebrick2")
  
  wave.ll <- sapply(seq(from=72,by=0.1,to=116), function(x) 
    wave_parameters["ll.mesor.estimate"]+wave_parameters["ll.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["ll.peak.time.hours"]*pi/12)))
  lines(x = seq(from=72,by=0.1,to=116),wave.ll,lwd=3,col="salmon")
  
  
  for(i in 0:2)
  {
    current.line <- (max.expr- min.expr)/10
    
    polygon(x = c(24*i, 24*i+8, 24*i+8, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="firebrick2")
    polygon(x = c(24*i+8,24*(i+1),24*(i+1),24*i+8),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="firebrick2",col="firebrick2")
  }
  
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(72, 96, 96, 72),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="firebrick2", col="salmon")
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(96, 120, 120, 96),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="firebrick2", col="salmon")
  
  return(0)  
}

plot.ld.ll.wave <- function(gene.id, gene.name, genes.expression, period=1, wave_parameters)
{
  current.gene.expression.ld <- 0
  
  ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  ld.zt.i <- paste(rep(ld.zt,3), rep(1:3, each=6), sep="_")
  ll.ld.zt <- paste("ll_ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  ll.ld.zt.i <- paste(rep(ll.ld.zt,2), rep(1:2, each=6), sep="_")
  
  current.gene.expression.ld <- genes.expression[gene.id,c(ld.zt.i)]
  current.gene.expression.ll <- genes.expression[gene.id,c(ll.ld.zt.i)]
  
  max.expr <- max(c(current.gene.expression.ld, current.gene.expression.ll))
  min.expr <- min(c(current.gene.expression.ld, current.gene.expression.ll))
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="",
       ylim=c((min.expr-(max.expr- min.expr)/5), max.expr),xlim=c(0,120),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  #abline(v=72, lwd=2, col= "green")
  axis(side=2,lwd=3)
  zt.points<-paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4)),sep="_")
  zt.points.rep<-c(paste(zt.points,1,sep="_"),paste(zt.points,2,sep="_"),paste(zt.points,3,sep="_"), paste(zt.points,4,sep="_"), paste(zt.points,5,sep="_"))
  axis(side = 1,lwd = 3, at=seq(from=0,to=116,by=4),labels = zt.points.rep,las=2)
  
  #lines(x = seq(from=0,by=4,to=116),current.gene.expression.ld,type="o",lwd=3,col="blue")
  points(x = seq(from=0,by=4,to=68),current.gene.expression.ld,lwd=3,col="royalblue3", pch=16)
  points(x = seq(from=72,by=4,to=116),current.gene.expression.ll,lwd=3,col="skyblue", pch=16)
  
  wave.ld <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["ld.mesor.estimate"]+wave_parameters["ld.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["ld.peak.time.hours"]*pi/12)))
  lines(x = seq(from=0,by=0.1,to=68),wave.ld,lwd=3,col="royalblue3")
  
  wave.ll <- sapply(seq(from=72,by=0.1,to=116), function(x) 
    wave_parameters["ll.mesor.estimate"]+wave_parameters["ll.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["ll.peak.time.hours"]*pi/12)))
  lines(x = seq(from=72,by=0.1,to=116),wave.ll,lwd=3,col="skyblue")
  
  
  for(i in 0:2)
  {
    current.line <- (max.expr- min.expr)/10
    
    polygon(x = c(24*i, 24*i+8, 24*i+8, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="royalblue3")
    polygon(x = c(24*i+8,24*(i+1),24*(i+1),24*i+8),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="royalblue3",col="royalblue3")
  }
  
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(72, 96, 96, 72),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="royalblue3", col="lightcyan")
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(96, 120, 120, 96),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="royalblue3", col="lightcyan")
  
  return(0)  
}

plot.ld.dd.wave <- function(gene.id, gene.name, genes.expression, period=1, wave_parameters)
{
  current.gene.expression.ld <- 0
  
  ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  ld.zt.i <- paste(rep(ld.zt,3), rep(1:3, each=6), sep="_")
  dd.ld.zt <- paste("dd_ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  dd.ld.zt.i <- paste(rep(dd.ld.zt,2), rep(1:2, each=6), sep="_")
  
  current.gene.expression.ld <- genes.expression[gene.id,c(ld.zt.i)]
  current.gene.expression.dd <- genes.expression[gene.id,c(dd.ld.zt.i)]
  
  max.expr <- max(c(current.gene.expression.ld, current.gene.expression.dd))
  min.expr <- min(c(current.gene.expression.ld, current.gene.expression.dd))
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="",
       ylim=c((min.expr-(max.expr- min.expr)/5), max.expr),xlim=c(0,120),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  #abline(v=72, lwd=2, col= "green")
  axis(side=2,lwd=3)
  zt.points<-paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4)),sep="_")
  zt.points.rep<-c(paste(zt.points,1,sep="_"),paste(zt.points,2,sep="_"),paste(zt.points,3,sep="_"), paste(zt.points,4,sep="_"), paste(zt.points,5,sep="_"))
  axis(side = 1,lwd = 3, at=seq(from=0,to=116,by=4),labels = zt.points.rep,las=2)
  
  #lines(x = seq(from=0,by=4,to=116),current.gene.expression.ld,type="o",lwd=3,col="blue")
  points(x = seq(from=0,by=4,to=68),current.gene.expression.ld,lwd=3,col="royalblue3", pch=16)
  points(x = seq(from=72,by=4,to=116),current.gene.expression.dd,lwd=3,col="slateblue4", pch=16)
  
  wave.ld <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["ld.mesor.estimate"]+wave_parameters["ld.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["ld.peak.time.hours"]*pi/12)))
  lines(x = seq(from=0,by=0.1,to=68),wave.ld,lwd=3,col="royalblue3")
  
  wave.dd <- sapply(seq(from=72,by=0.1,to=116), function(x) 
    wave_parameters["dd.mesor.estimate"]+wave_parameters["dd.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["dd.peak.time.hours"]*pi/12)))
  lines(x = seq(from=72,by=0.1,to=116),wave.dd,lwd=3,col="slateblue4")
  
  
  for(i in 0:2)
  {
    current.line <- (max.expr- min.expr)/10
    
    polygon(x = c(24*i, 24*i+8, 24*i+8, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="royalblue3")
    polygon(x = c(24*i+8,24*(i+1),24*(i+1),24*i+8),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="royalblue3",col="royalblue3")
  }
  
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(72, 96, 96, 72),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="royalblue3", col="slateblue4")
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(96, 120, 120, 96),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="royalblue3", col="slateblue4")
  
  return(0)  
}

plot.sd.dd.wave <- function(gene.id, gene.name, genes.expression, period=1, wave_parameters)
{
  current.gene.expression.sd <- 0
  
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  sd.zt.i <- paste(rep(sd.zt,3), rep(1:3, each=6), sep="_")
  dd.sd.zt <- paste("dd_sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  dd.sd.zt.i <- paste(rep(dd.sd.zt,2), rep(1:2, each=6), sep="_")
  
  current.gene.expression.sd <- genes.expression[gene.id,c(sd.zt.i)]
  current.gene.expression.dd <- genes.expression[gene.id,c(dd.sd.zt.i)]
  
  max.expr <- max(c(current.gene.expression.sd, current.gene.expression.dd))
  min.expr <- min(c(current.gene.expression.sd, current.gene.expression.dd))
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="",
       ylim=c((min.expr-(max.expr- min.expr)/5), max.expr),xlim=c(0,120),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  #abline(v=72, lwd=2, col= "green")
  axis(side=2,lwd=3)
  zt.points<-paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4)),sep="_")
  zt.points.rep<-c(paste(zt.points,1,sep="_"),paste(zt.points,2,sep="_"),paste(zt.points,3,sep="_"), paste(zt.points,4,sep="_"), paste(zt.points,5,sep="_"))
  axis(side = 1,lwd = 3, at=seq(from=0,to=116,by=4),labels = zt.points.rep,las=2)
  
  #lines(x = seq(from=0,by=4,to=116),current.gene.expression.sd,type="o",lwd=3,col="blue")
  points(x = seq(from=0,by=4,to=68),current.gene.expression.sd,lwd=3,col="firebrick2", pch=16)
  points(x = seq(from=72,by=4,to=116),current.gene.expression.dd,lwd=3,col="red4", pch=16)
  
  wave.sd <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["sd.mesor.estimate"]+wave_parameters["sd.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["sd.peak.time.hours"]*pi/12)))
  lines(x = seq(from=0,by=0.1,to=68),wave.sd,lwd=3,col="firebrick2")
  
  wave.dd <- sapply(seq(from=72,by=0.1,to=116), function(x) 
    wave_parameters["dd.mesor.estimate"]+wave_parameters["dd.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["dd.peak.time.hours"]*pi/12)))
  lines(x = seq(from=72,by=0.1,to=116),wave.dd,lwd=3,col="red4")
  
  
  for(i in 0:2)
  {
    current.line <- (max.expr- min.expr)/10
    
    polygon(x = c(24*i, 24*i+8, 24*i+8, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="firebrick2")
    polygon(x = c(24*i+8,24*(i+1),24*(i+1),24*i+8),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="firebrick2",col="firebrick2")
  }
  
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(72, 96, 96, 72),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="firebrick2", col="red4")
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(96, 120, 120, 96),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="firebrick2", col="red4")
  
  return(0)  
}

create_motif_table <- function(non.redundant.motif.file)
{
  motifs_lines <- readLines(non.redundant.motif.file)
  motifs_clean <- motifs_lines[grep(">", motifs_lines)]
  motifs_desc <- sapply(strsplit(motifs_clean, "\t"), function(x) x[2])
  motifs_desc <- sapply(strsplit(motifs_desc, "BestGuess:"), function(x) x[2])
  motifs_data <- sapply(strsplit(motifs_clean, "\t"), function(x) x[6])
  motifs_p_value <- as.numeric(sapply(strsplit(motifs_data, "P:"), function(x) x[2]))
  motifs_target <- sapply(strsplit(motifs_data, "%"), function(x) x[1])
  motifs_target <- as.numeric(sapply(strsplit(motifs_target, "[(]"), function(x) x[2]))
  motifs_bg <- sapply(strsplit(motifs_data, "%"), function(x) x[2])
  motifs_bg <- as.numeric(sapply(strsplit(motifs_bg, "[(]"), function(x) x[2]))
  motifs_motif_table <- data.frame(motif=motifs_desc, target=motifs_target, 
                                   bg=motifs_bg, p_value=motifs_p_value)
  return(motifs_motif_table)
}

heatmap.gene <- function(gene.id)
{
  current.gene.mean.sd.gene.expression <- (scaled.sd.gene.expression[gene.id,seq(from=1,to=16,by=3)] + scaled.sd.gene.expression[gene.id,seq(from=2,to=17,by=3)] + scaled.sd.gene.expression[gene.id,seq(from=3,to=18,by=3)])/3
  current.gene.mean.sd.gene.expression <- c(current.gene.mean.sd.gene.expression, current.gene.mean.sd.gene.expression[1])
  spline.current.gene.mean.sd.gene.expression <- spline(x=seq(from=0,to=24,by=4),y=current.gene.mean.sd.gene.expression,n = 100)
  col_fun_gene = colorRamp2(c(min(spline.current.gene.mean.sd.gene.expression$y),
                              (min(spline.current.gene.mean.sd.gene.expression$y) + max(spline.current.gene.mean.sd.gene.expression$y)) /2,
                              max(spline.current.gene.mean.sd.gene.expression$y)), 
                            c("black","blue","yellow"))
  cols.gene <- col_fun_gene(x = spline.current.gene.mean.sd.gene.expression$y)
  
  par(mar=c(0,0,0,0))
  plot(x=1,y=1,col="white",xlim=c(-1,25),ylim=c(0,14),axes=F,xlab="",ylab="")
  
  polygon(x = c(0,8,8,0),y=c(2, 2, 0, 0),lwd=2,border="red")
  polygon(x = c(8,24,24,8),y=c(2, 2, 0, 0),lwd=2,border="red",col="red") 
  
  for(i in 1:99)
  {
    polygon(x = c(spline.current.gene.mean.sd.gene.expression$x[i],spline.current.gene.mean.sd.gene.expression$x[i+1],
                  spline.current.gene.mean.sd.gene.expression$x[i+1],spline.current.gene.mean.sd.gene.expression$x[i]),
            y=c(2.5, 2.5, 5, 5),lwd=2,border=cols.gene[i],col=cols.gene[i]) 
  }
  
  current.gene.mean.ld.gene.expression <- (scaled.ld.gene.expression[gene.id,seq(from=1,to=16,by=3)] + scaled.ld.gene.expression[gene.id,seq(from=2,to=17,by=3)] + scaled.ld.gene.expression[gene.id,seq(from=3,to=18,by=3)])/3
  current.gene.mean.ld.gene.expression <- c(current.gene.mean.ld.gene.expression, current.gene.mean.ld.gene.expression[1])
  spline.current.gene.mean.ld.gene.expression <- spline(x=seq(from=0,to=24,by=4),y=current.gene.mean.ld.gene.expression,n = 100)
  col_fun_gene = colorRamp2(c(min(spline.current.gene.mean.ld.gene.expression$y),
                              (min(spline.current.gene.mean.ld.gene.expression$y) + max(spline.current.gene.mean.ld.gene.expression$y)) /2,
                              max(spline.current.gene.mean.ld.gene.expression$y)), 
                            c("black","blue","yellow"))
  cols.gene <- col_fun_gene(x = spline.current.gene.mean.ld.gene.expression$y)
  
  
  polygon(x = c(0,16,16,0),y=c(8.5, 8.5, 6.5, 6.5),lwd=2,border="blue")
  polygon(x = c(16,24,24,16),y=c(8.5, 8.5, 6.5, 6.5),lwd=2,border="blue",col="blue") 
  
  for(i in 1:99)
  {
    polygon(x = c(spline.current.gene.mean.ld.gene.expression$x[i],spline.current.gene.mean.ld.gene.expression$x[i+1],
                  spline.current.gene.mean.ld.gene.expression$x[i+1],spline.current.gene.mean.ld.gene.expression$x[i]),
            y=c(9, 9, 11.5, 11.5),lwd=2,border=cols.gene[i],col=cols.gene[i]) 
  }
  
}


# Plot transitions with errors
# LD LL
plot.ld.ll.wave.error <- function(gene.id, gene.name="", genes.expression, period=1, 
                                  wave_parameters.df, data.points = F, 
                                  color.1 = "royalblue3", color.2="skyblue")
{
  wave_parameters <- wave_parameters.df[gene.id,]
  
  wave.ld <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["ld.mesor.estimate"]+wave_parameters["ld.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["ld.peak.time.hours"]*pi/12)))
  
  wave.ld.ggplot <- data.frame(x = seq(from=0,by=0.1,to=68), y = unlist(wave.ld),
                               min.bound = unlist(wave.ld) - wave_parameters$error_amp_1,
                               max.bound = unlist(wave.ld) + wave_parameters$error_amp_1)
  points_ld <- data.frame(x=seq(from=0,by=4,to=68), y=genes.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"), paste(ld.zt,3,sep="_"))])
  
  wave.ll <- sapply(seq(from=72,by=0.1,to=116), function(x) 
    wave_parameters["ll.mesor.estimate"]+wave_parameters["ll.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["ll.peak.time.hours"]*pi/12)))
  
  wave.ll.ggplot <- data.frame(x = seq(from=72,by=0.1,to=116), y = unlist(wave.ll),
                               min.bound = unlist(wave.ll) - wave_parameters$error_amp_2,
                               max.bound = unlist(wave.ll) + wave_parameters$error_amp_2)
  points_ll <- data.frame(x=seq(from=72,by=4,to=116), y=genes.expression[gene.id,c(paste(ll.ld.zt,1,sep="_"),
                                                                                   paste(ll.ld.zt,2,sep="_"))])
  
  min.limit <- min(c(points_ld$y), points_ll$y, wave.ld.ggplot$min.bound, wave.ll.ggplot$min.bound)
  max.limit <- max(c(points_ld$y), points_ll$y, wave.ld.ggplot$max.bound, wave.ll.ggplot$max.bound)
  
  df_polygon_ld <- data.frame(
    x = c(0,16,16,0,16,24,24,16,24,40,40,24,40,48,48,40,48,64,64,48,64,72,72,64,
          72,88,88,72,88,96,96,88,96,112,112,96,112,120,120,112),
    
    y = rep(c(min.limit - (max.limit - min.limit)/25, min.limit - (max.limit - min.limit)/25, 
              min.limit - (max.limit - min.limit)/50, min.limit - (max.limit - min.limit)/50), 10),
    group = factor(rep(1:10, each = 4))
  )
  
  if(data.points)
  {
    ggplot() +
      geom_ribbon(data = wave.ld.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.1, alpha = 0.2) +
      geom_line(data = wave.ld.ggplot, aes(x = x, y = y), color = color.1, size = 1) +
      geom_point(data = points_ld, aes(x = x, y = y), color = color.1, size = 1, alpha=0.3) +
      
      geom_ribbon(data = wave.ll.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.2, alpha = 0.2) +
      geom_line(data = wave.ll.ggplot, aes(x = x, y = y), color = color.2, size = 1) +
      geom_point(data = points_ll, aes(x = x, y = y), color = color.2, size = 1, alpha=0.3) +
      scale_x_continuous(
        breaks = seq(from=0,to=116,by=4),                
        labels = c(paste0("ZT", rep(seq(0,20,4),3)), paste0("CT", rep(seq(0,20,4),2)))) +
      scale_y_continuous(limits = c(min.limit - (max.limit - min.limit)/25, 
                                    max.limit)) +
      labs(title = gene.id,
           y = "FPKM",
           x = "") +
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
          size = 10, angle = 90, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.text.y = element_text(
          size = 12, angle = 0, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.title.y = element_text(
          size = 12, angle = 90, vjust = 2.5, hjust = 0.5,
          face = "bold", color = "black"
        )) +
      geom_polygon(
        data = df_polygon_ld,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.1
      ) +
      scale_fill_manual(values = c(rep(c("white", color.1), 3), rep(c("white", color.2), 2)))
  }
  else
  {
    ggplot() +
      geom_ribbon(data = wave.ld.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.1, alpha = 0.2) +
      geom_line(data = wave.ld.ggplot, aes(x = x, y = y), color = color.1, size = 1) +
      
      geom_ribbon(data = wave.ll.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.2, alpha = 0.2) +
      geom_line(data = wave.ll.ggplot, aes(x = x, y = y), color = color.2, size = 1) +
      scale_x_continuous(
        breaks = seq(from=0,to=116,by=4),                
        labels = c(paste0("ZT", rep(seq(0,20,4),3)), paste0("CT", rep(seq(0,20,4),2)))) +
      scale_y_continuous(limits = c(min.limit - (max.limit - min.limit)/25, 
                                    max.limit)) +
      labs(title = gene.id,
           y = "FPKM",
           x = "") +
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
          size = 10, angle = 90, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.text.y = element_text(
          size = 12, angle = 0, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.title.y = element_text(
          size = 12, angle = 90, vjust = 2.5, hjust = 0.5,
          face = "bold", color = "black"
        )) +
      geom_polygon(
        data = df_polygon_ld,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.1
      ) +
      scale_fill_manual(values = c(rep(c("white", color.1), 3), rep(c("white", color.2), 2)))
  }
  
}

# LD DD
plot.ld.dd.wave.error <- function(gene.id, gene.name="", genes.expression, period=1, 
                                  wave_parameters.df, data.points = F, 
                                  color.1 = "royalblue3", color.2="darkblue")
{
  wave_parameters <- wave_parameters.df[gene.id,]
  
  wave.ld <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["ld.mesor.estimate"]+wave_parameters["ld.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["ld.peak.time.hours"]*pi/12)))
  
  wave.ld.ggplot <- data.frame(x = seq(from=0,by=0.1,to=68), y = unlist(wave.ld),
                               min.bound = unlist(wave.ld) - wave_parameters$error_amp_2,
                               max.bound = unlist(wave.ld) + wave_parameters$error_amp_2)
  points_ld <- data.frame(x=seq(from=0,by=4,to=68), y=genes.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"), paste(ld.zt,3,sep="_"))])
  
  wave.dd <- sapply(seq(from=72,by=0.1,to=116), function(x) 
    wave_parameters["dd.mesor.estimate"]+wave_parameters["dd.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["dd.peak.time.hours"]*pi/12)))
  
  wave.dd.ggplot <- data.frame(x = seq(from=72,by=0.1,to=116), y = unlist(wave.dd),
                               min.bound = unlist(wave.dd) - wave_parameters$error_amp_1,
                               max.bound = unlist(wave.dd) + wave_parameters$error_amp_1)
  points_dd <- data.frame(x=seq(from=72,by=4,to=116), y=genes.expression[gene.id,c(paste(dd.ld.zt,1,sep="_"),
                                                                                   paste(dd.ld.zt,2,sep="_"))])
  
  min.limit <- min(c(points_ld$y), points_dd$y, wave.ld.ggplot$min.bound, wave.dd.ggplot$min.bound)
  max.limit <- max(c(points_ld$y), points_dd$y, wave.ld.ggplot$max.bound, wave.dd.ggplot$max.bound)
  
  df_polygon_ld <- data.frame(
    x = c(0,16,16,0,16,24,24,16,24,40,40,24,40,48,48,40,48,64,64,48,64,72,72,64,
          72,88,88,72,88,96,96,88,96,112,112,96,112,120,120,112),
    
    y = rep(c(min.limit - (max.limit - min.limit)/25, min.limit - (max.limit - min.limit)/25, 
              min.limit - (max.limit - min.limit)/50, min.limit - (max.limit - min.limit)/50), 10),
    group = factor(rep(1:10, each = 4))
  )
  
  if(data.points)
  {
    ggplot() +
      geom_ribbon(data = wave.ld.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.1, alpha = 0.2) +
      geom_line(data = wave.ld.ggplot, aes(x = x, y = y), color = color.1, size = 1) +
      geom_point(data = points_ld, aes(x = x, y = y), color = color.1, size = 1, alpha=0.3) +
      
      geom_ribbon(data = wave.dd.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.2, alpha = 0.2) +
      geom_line(data = wave.dd.ggplot, aes(x = x, y = y), color = color.2, size = 1) +
      geom_point(data = points_dd, aes(x = x, y = y), color = color.2, size = 1, alpha=0.3) +
      scale_x_continuous(
        breaks = seq(from=0,to=116,by=4),                
        labels = c(paste0("ZT", rep(seq(0,20,4),3)), paste0("CT", rep(seq(0,20,4),2)))) +
      scale_y_continuous(limits = c(min.limit - (max.limit - min.limit)/25, 
                                    max.limit)) +
      labs(title = gene.id,
           y = "FPKM",
           x = "") +
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
          size = 10, angle = 90, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.text.y = element_text(
          size = 12, angle = 0, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.title.y = element_text(
          size = 12, angle = 90, vjust = 2.5, hjust = 0.5,
          face = "bold", color = "black"
        )) +
      geom_polygon(
        data = df_polygon_ld,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.1
      ) +
      scale_fill_manual(values = c(rep(c("white", color.1), 3), rep(c(color.2, color.1), 2)))
  }
  else
  {
    ggplot() +
      geom_ribbon(data = wave.ld.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.1, alpha = 0.2) +
      geom_line(data = wave.ld.ggplot, aes(x = x, y = y), color = color.1, size = 1) +
      
      geom_ribbon(data = wave.dd.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.2, alpha = 0.2) +
      geom_line(data = wave.dd.ggplot, aes(x = x, y = y), color = color.2, size = 1) +
      scale_x_continuous(
        breaks = seq(from=0,to=116,by=4),                
        labels = c(paste0("ZT", rep(seq(0,20,4),3)), paste0("CT", rep(seq(0,20,4),2)))) +
      scale_y_continuous(limits = c(min.limit - (max.limit - min.limit)/25, 
                                    max.limit)) +
      labs(title = gene.id,
           y = "FPKM",
           x = "") +
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
          size = 10, angle = 90, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.text.y = element_text(
          size = 12, angle = 0, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.title.y = element_text(
          size = 12, angle = 90, vjust = 2.5, hjust = 0.5,
          face = "bold", color = "black"
        )) +
      geom_polygon(
        data = df_polygon_ld,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.1
      ) +
      scale_fill_manual(values = c(rep(c("white", color.1), 3), rep(c(color.2, color.1), 2)))
  }
  
}

# SD LL
plot.sd.ll.wave.error <- function(gene.id, gene.name="", genes.expression, period=1, 
                                  wave_parameters.df, data.points = F, 
                                  color.1 = "tomato3", color.2="salmon")
{
  wave_parameters <- wave_parameters.df[gene.id,]
  
  wave.sd <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["sd.mesor.estimate"]+wave_parameters["sd.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["sd.peak.time.hours"]*pi/12)))
  
  wave.sd.ggplot <- data.frame(x = seq(from=0,by=0.1,to=68), y = unlist(wave.sd),
                               min.bound = unlist(wave.sd) - wave_parameters$error_amp_2,
                               max.bound = unlist(wave.sd) + wave_parameters$error_amp_2)
  points_sd <- data.frame(x=seq(from=0,by=4,to=68), y=genes.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"), paste(sd.zt,3,sep="_"))])
  
  wave.ll <- sapply(seq(from=72,by=0.1,to=116), function(x) 
    wave_parameters["ll.mesor.estimate"]+wave_parameters["ll.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["ll.peak.time.hours"]*pi/12)))
  
  wave.ll.ggplot <- data.frame(x = seq(from=72,by=0.1,to=116), y = unlist(wave.ll),
                               min.bound = unlist(wave.ll) - wave_parameters$error_amp_1,
                               max.bound = unlist(wave.ll) + wave_parameters$error_amp_1)
  points_ll <- data.frame(x=seq(from=72,by=4,to=116), y=genes.expression[gene.id,c(paste(ll.sd.zt,1,sep="_"),
                                                                                   paste(ll.sd.zt,2,sep="_"))])
  
  min.limit <- min(c(points_sd$y), points_ll$y, wave.sd.ggplot$min.bound, wave.ll.ggplot$min.bound)
  max.limit <- max(c(points_sd$y), points_ll$y, wave.sd.ggplot$max.bound, wave.ll.ggplot$max.bound)
  
  df_polygon_sd <- data.frame(
    x = c(0,8,8,0,8,24,24,8,24,32,32,24,32,48,48,32,48,56,56,48,56,72,72,56,
          72,80,80,72,80,96,96,80,96,104,104,96,104,120,120,104),
    
    y = rep(c(min.limit - (max.limit - min.limit)/25, min.limit - (max.limit - min.limit)/25, 
              min.limit - (max.limit - min.limit)/50, min.limit - (max.limit - min.limit)/50), 10),
    group = factor(rep(1:10, each = 4))
  )
  
  if(data.points)
  {
    ggplot() +
      geom_ribbon(data = wave.sd.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.1, alpha = 0.2) +
      geom_line(data = wave.sd.ggplot, aes(x = x, y = y), color = color.1, size = 1) +
      geom_point(data = points_sd, aes(x = x, y = y), color = color.1, size = 1, alpha=0.3) +
      
      geom_ribbon(data = wave.ll.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.2, alpha = 0.2) +
      geom_line(data = wave.ll.ggplot, aes(x = x, y = y), color = color.2, size = 1) +
      geom_point(data = points_ll, aes(x = x, y = y), color = color.2, size = 1, alpha=0.3) +
      scale_x_continuous(
        breaks = seq(from=0,to=116,by=4),                
        labels = c(paste0("ZT", rep(seq(0,20,4),3)), paste0("CT", rep(seq(0,20,4),2)))) +
      scale_y_continuous(limits = c(min.limit - (max.limit - min.limit)/25, 
                                    max.limit)) +
      labs(title = gene.id,
           y = "FPKM",
           x = "") +
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
          size = 10, angle = 90, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.text.y = element_text(
          size = 12, angle = 0, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.title.y = element_text(
          size = 12, angle = 90, vjust = 2.5, hjust = 0.5,
          face = "bold", color = "black"
        )) +
      geom_polygon(
        data = df_polygon_sd,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.1
      ) +
      scale_fill_manual(values = c(rep(c("white", color.1), 3), rep(c("white", color.2), 2)))
  }
  else
  {
    ggplot() +
      geom_ribbon(data = wave.sd.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.1, alpha = 0.2) +
      geom_line(data = wave.sd.ggplot, aes(x = x, y = y), color = color.1, size = 1) +
      
      geom_ribbon(data = wave.ll.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.2, alpha = 0.2) +
      geom_line(data = wave.ll.ggplot, aes(x = x, y = y), color = color.2, size = 1) +
      scale_x_continuous(
        breaks = seq(from=0,to=116,by=4),                
        labels = c(paste0("ZT", rep(seq(0,20,4),3)), paste0("CT", rep(seq(0,20,4),2)))) +
      scale_y_continuous(limits = c(min.limit - (max.limit - min.limit)/25, 
                                    max.limit)) +
      labs(title = gene.id,
           y = "FPKM",
           x = "") +
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
          size = 10, angle = 90, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.text.y = element_text(
          size = 12, angle = 0, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.title.y = element_text(
          size = 12, angle = 90, vjust = 2.5, hjust = 0.5,
          face = "bold", color = "black"
        )) +
      geom_polygon(
        data = df_polygon_sd,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.1
      ) +
      scale_fill_manual(values = c(rep(c("white", color.1), 3), rep(c("white", color.2), 2)))
  }
  
}

# SD DD
plot.sd.dd.wave.error <- function(gene.id, gene.name="", genes.expression, period=1, 
                                  wave_parameters.df, data.points = F, 
                                  color.1 = "tomato3", color.2="red4")
{
  wave_parameters <- wave_parameters.df[gene.id,]
  
  wave.sd <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["sd.mesor.estimate"]+wave_parameters["sd.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["sd.peak.time.hours"]*pi/12)))
  
  wave.sd.ggplot <- data.frame(x = seq(from=0,by=0.1,to=68), y = unlist(wave.sd),
                               min.bound = unlist(wave.sd) - wave_parameters$error_amp_2,
                               max.bound = unlist(wave.sd) + wave_parameters$error_amp_2)
  points_sd <- data.frame(x=seq(from=0,by=4,to=68), y=genes.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"), paste(sd.zt,3,sep="_"))])
  
  wave.dd <- sapply(seq(from=72,by=0.1,to=116), function(x) 
    wave_parameters["dd.mesor.estimate"]+wave_parameters["dd.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["dd.peak.time.hours"]*pi/12)))
  
  wave.dd.ggplot <- data.frame(x = seq(from=72,by=0.1,to=116), y = unlist(wave.dd),
                               min.bound = unlist(wave.dd) - wave_parameters$error_amp_1,
                               max.bound = unlist(wave.dd) + wave_parameters$error_amp_1)
  points_dd <- data.frame(x=seq(from=72,by=4,to=116), y=genes.expression[gene.id,c(paste(dd.sd.zt,1,sep="_"),
                                                                                   paste(dd.sd.zt,2,sep="_"))])
  
  min.limit <- min(c(points_sd$y), points_dd$y, wave.sd.ggplot$min.bound, wave.dd.ggplot$min.bound)
  max.limit <- max(c(points_sd$y), points_dd$y, wave.sd.ggplot$max.bound, wave.dd.ggplot$max.bound)
  
  df_polygon_sd <- data.frame(
    x = c(0,8,8,0,8,24,24,8,24,32,32,24,32,48,48,32,48,56,56,48,56,72,72,56,
          72,80,80,72,80,96,96,80,96,104,104,96,104,120,120,104),
    
    y = rep(c(min.limit - (max.limit - min.limit)/25, min.limit - (max.limit - min.limit)/25, 
              min.limit - (max.limit - min.limit)/50, min.limit - (max.limit - min.limit)/50), 10),
    group = factor(rep(1:10, each = 4))
  )
  
  if(data.points)
  {
    ggplot() +
      geom_ribbon(data = wave.sd.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.1, alpha = 0.2) +
      geom_line(data = wave.sd.ggplot, aes(x = x, y = y), color = color.1, size = 1) +
      geom_point(data = points_sd, aes(x = x, y = y), color = color.1, size = 1, alpha=0.3) +
      
      geom_ribbon(data = wave.dd.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.2, alpha = 0.2) +
      geom_line(data = wave.dd.ggplot, aes(x = x, y = y), color = color.2, size = 1) +
      geom_point(data = points_dd, aes(x = x, y = y), color = color.2, size = 1, alpha=0.3) +
      scale_x_continuous(
        breaks = seq(from=0,to=116,by=4),                
        labels = c(paste0("ZT", rep(seq(0,20,4),3)), paste0("CT", rep(seq(0,20,4),2)))) +
      scale_y_continuous(limits = c(min.limit - (max.limit - min.limit)/25, 
                                    max.limit)) +
      labs(title = gene.id,
           y = "FPKM",
           x = "") +
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
          size = 10, angle = 90, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.text.y = element_text(
          size = 12, angle = 0, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.title.y = element_text(
          size = 12, angle = 90, vjust = 2.5, hjust = 0.5,
          face = "bold", color = "black"
        )) +
      geom_polygon(
        data = df_polygon_sd,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.1
      ) +
      scale_fill_manual(values = c(rep(c("white", color.1), 3), rep(c(color.2,color.1), 2)))
  }
  else
  {
    ggplot() +
      geom_ribbon(data = wave.sd.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.1, alpha = 0.2) +
      geom_line(data = wave.sd.ggplot, aes(x = x, y = y), color = color.1, size = 1) +
      
      geom_ribbon(data = wave.dd.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.2, alpha = 0.2) +
      geom_line(data = wave.dd.ggplot, aes(x = x, y = y), color = color.2, size = 1) +
      scale_x_continuous(
        breaks = seq(from=0,to=116,by=4),                
        labels = c(paste0("ZT", rep(seq(0,20,4),3)), paste0("CT", rep(seq(0,20,4),2)))) +
      scale_y_continuous(limits = c(min.limit - (max.limit - min.limit)/25, 
                                    max.limit)) +
      labs(title = gene.id,
           y = "FPKM",
           x = "") +
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
          size = 10, angle = 90, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.text.y = element_text(
          size = 12, angle = 0, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.title.y = element_text(
          size = 12, angle = 90, vjust = 2.5, hjust = 0.5,
          face = "bold", color = "black"
        )) +
      geom_polygon(
        data = df_polygon_sd,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.1
      ) +
      scale_fill_manual(values = c(rep(c("white", color.1), 3), rep(c(color.2, color.1), 2)))
  }
  
}

# LD SD
plot.ld.sd.wave.error <- function(gene.id, gene.name="", genes.expression, period=1, 
                                  wave_parameters.df, data.points = F, 
                                  color.1 = "royalblue3", color.2="tomato3")
{
  wave_parameters <- wave_parameters.df[gene.id,]
  
  wave.ld <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["ld.mesor.estimate"]+wave_parameters["ld.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["ld.peak.time.hours"]*pi/12)))
  
  wave.ld.ggplot <- data.frame(x = seq(from=0,by=0.1,to=68), y = unlist(wave.ld),
                               min.bound = unlist(wave.ld) - wave_parameters$error_amp_1,
                               max.bound = unlist(wave.ld) + wave_parameters$error_amp_1)
  points_ld <- data.frame(x=seq(from=0,by=4,to=68), y=genes.expression[gene.id,c(paste(ld.zt,1,sep="_"),
                                                                                 paste(ld.zt,2,sep="_"),
                                                                                 paste(ld.zt,3,sep="_"))])
  
  wave.sd <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["sd.mesor.estimate"]+wave_parameters["sd.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["sd.peak.time.hours"]*pi/12)))
  
  wave.sd.ggplot <- data.frame(x = seq(from=0,by=0.1,to=68), y = unlist(wave.sd),
                               min.bound = unlist(wave.sd) - wave_parameters$error_amp_2,
                               max.bound = unlist(wave.sd) + wave_parameters$error_amp_2)
  points_sd <- data.frame(x=seq(from=0,by=4,to=68), y=genes.expression[gene.id,c(paste(sd.zt,1,sep="_"),
                                                                                 paste(sd.zt,2,sep="_"),
                                                                                 paste(sd.zt,3,sep="_"))])
  
  min.limit <- min(c(points_ld$y), points_sd$y, wave.ld.ggplot$min.bound, wave.sd.ggplot$min.bound)
  max.limit <- max(c(points_ld$y), points_sd$y, wave.ld.ggplot$max.bound, wave.sd.ggplot$max.bound)
  
  df_polygon_ld <- data.frame(
    x = c(0,16,16,0,16,24,24,16,24,40,40,24,40,48,48,40,48,64,64,48,64,72,72,64),
    
    y = rep(c(min.limit - (max.limit - min.limit)/25, min.limit - (max.limit - min.limit)/25, 
              min.limit - (max.limit - min.limit)/50, min.limit - (max.limit - min.limit)/50), 6),
    group = factor(rep(1:6, each = 4))
  )
  
  df_polygon_sd <- data.frame(
    x = c(0,8,8,0,8,24,24,8,24,32,32,24,32,48,48,32,48,56,56,48,56,72,72,56),
    
    y = rep(c(min.limit - 3*(max.limit - min.limit)/50, min.limit - 3*(max.limit - min.limit)/50,
              min.limit - (max.limit - min.limit)/25, min.limit - (max.limit - min.limit)/25), 6),
    group = factor(rep(7:12, each = 4))
  )
  
  
  if(data.points)
  {
    ggplot() +
      geom_ribbon(data = wave.ld.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.1, alpha = 0.2) +
      geom_line(data = wave.ld.ggplot, aes(x = x, y = y), color = color.1, size = 1) +
      geom_point(data = points_ld, aes(x = x, y = y), color = color.1, size = 1, alpha=0.3) +
      
      geom_ribbon(data = wave.sd.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.2, alpha = 0.2) +
      geom_line(data = wave.sd.ggplot, aes(x = x, y = y), color = color.2, size = 1) +
      geom_point(data = points_sd, aes(x = x, y = y), color = color.2, size = 1, alpha=0.3) +
      scale_x_continuous(
        breaks = seq(from=0,to=68,by=4),                
        labels = paste0("ZT", rep(seq(0,20,4),3))) +
      scale_y_continuous(limits = c(min.limit - (max.limit - min.limit)/12.5, 
                                    max.limit)) +
      labs(title = gene.id,
           y = "FPKM",
           x = "") +
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
          size = 10, angle = 90, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.text.y = element_text(
          size = 12, angle = 0, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.title.y = element_text(
          size = 12, angle = 90, vjust = 2.5, hjust = 0.5,
          face = "bold", color = "black"
        )) +
      geom_polygon(
        data = df_polygon_ld,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.1
      ) +
      geom_polygon(
        data = df_polygon_sd,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.2
      ) +
      scale_fill_manual(values = c(rep(c("white", color.1), 3), rep(c("white", color.2),3)))
  }
  else
  {
    ggplot() +
      geom_ribbon(data = wave.ld.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.1, alpha = 0.2) +
      geom_line(data = wave.ld.ggplot, aes(x = x, y = y), color = color.1, size = 1) +
      
      geom_ribbon(data = wave.sd.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.2, alpha = 0.2) +
      geom_line(data = wave.sd.ggplot, aes(x = x, y = y), color = color.2, size = 1) +
      scale_x_continuous(
        breaks = seq(from=0,to=116,by=4),                
        labels = c(paste0("ZT", rep(seq(0,20,4),3)), paste0("ZT", rep(seq(0,20,4),2)))) +
      scale_y_continuous(limits = c(min.limit - (max.limit - min.limit)/12.5, 
                                    max.limit)) +
      labs(title = gene.id,
           y = "FPKM",
           x = "") +
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
          size = 10, angle = 90, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.text.y = element_text(
          size = 12, angle = 0, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.title.y = element_text(
          size = 12, angle = 90, vjust = 2.5, hjust = 0.5,
          face = "bold", color = "black"
        )) +
      geom_polygon(
        data = df_polygon_ld,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.1
      ) +
      geom_polygon(
        data = df_polygon_sd,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.2
      ) +
      scale_fill_manual(values = c(rep(c("white", color.1), 3), rep(c("white", color.2),3)))
  }
  
}

# LD SD with one or more conditions exhibiting a 12 h period
plot.ld.sd.wave.error.12 <- function(gene.id, gene.name="", genes.expression, period1=1, period2=1,
                                  wave_parameters.df1, wave_parameters.df2, data.points = F, 
                                  color.1 = "royalblue3", color.2="tomato3")
{
  wave_parameters1 <- wave_parameters.df1[gene.id,]
  wave_parameters2 <- wave_parameters.df2[gene.id,]
  
  wave.ld <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters1["mesor"]+wave_parameters1["amplitude"]*cos(period1*(x*pi/12-wave_parameters1["peak_time_hours"]*pi/12)))
  
  wave.ld.ggplot <- data.frame(x = seq(from=0,by=0.1,to=68), y = unlist(wave.ld),
                               min.bound = unlist(wave.ld) - wave_parameters1$error_amp,
                               max.bound = unlist(wave.ld) + wave_parameters1$error_amp)
  points_ld <- data.frame(x=seq(from=0,by=4,to=68), y=genes.expression[gene.id,c(paste(ld.zt,1,sep="_"),
                                                                                 paste(ld.zt,2,sep="_"),
                                                                                 paste(ld.zt,3,sep="_"))])
  
  wave.sd <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters2["mesor"]+wave_parameters2["amplitude"]*cos(period2*(x*pi/12-wave_parameters2["peak_time_hours"]*pi/12)))
  
  wave.sd.ggplot <- data.frame(x = seq(from=0,by=0.1,to=68), y = unlist(wave.sd),
                               min.bound = unlist(wave.sd) - wave_parameters2$error_amp,
                               max.bound = unlist(wave.sd) + wave_parameters2$error_amp)
  points_sd <- data.frame(x=seq(from=0,by=4,to=68), y=genes.expression[gene.id,c(paste(sd.zt,1,sep="_"),
                                                                                 paste(sd.zt,2,sep="_"),
                                                                                 paste(sd.zt,3,sep="_"))])
  
  min.limit <- min(c(points_ld$y), points_sd$y, wave.ld.ggplot$min.bound, wave.sd.ggplot$min.bound)
  max.limit <- max(c(points_ld$y), points_sd$y, wave.ld.ggplot$max.bound, wave.sd.ggplot$max.bound)
  
  df_polygon_ld <- data.frame(
    x = c(0,16,16,0,16,24,24,16,24,40,40,24,40,48,48,40,48,64,64,48,64,72,72,64),
    
    y = rep(c(min.limit - (max.limit - min.limit)/25, min.limit - (max.limit - min.limit)/25, 
              min.limit - (max.limit - min.limit)/50, min.limit - (max.limit - min.limit)/50), 6),
    group = factor(rep(1:6, each = 4))
  )
  
  df_polygon_sd <- data.frame(
    x = c(0,8,8,0,8,24,24,8,24,32,32,24,32,48,48,32,48,56,56,48,56,72,72,56),
    
    y = rep(c(min.limit - 3*(max.limit - min.limit)/50, min.limit - 3*(max.limit - min.limit)/50,
              min.limit - (max.limit - min.limit)/25, min.limit - (max.limit - min.limit)/25), 6),
    group = factor(rep(7:12, each = 4))
  )
  
  
  if(data.points)
  {
    ggplot() +
      geom_ribbon(data = wave.ld.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.1, alpha = 0.2) +
      geom_line(data = wave.ld.ggplot, aes(x = x, y = y), color = color.1, size = 1) +
      geom_point(data = points_ld, aes(x = x, y = y), color = color.1, size = 1, alpha=0.3) +
      
      geom_ribbon(data = wave.sd.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.2, alpha = 0.2) +
      geom_line(data = wave.sd.ggplot, aes(x = x, y = y), color = color.2, size = 1) +
      geom_point(data = points_sd, aes(x = x, y = y), color = color.2, size = 1, alpha=0.3) +
      scale_x_continuous(
        breaks = seq(from=0,to=68,by=4),                
        labels = paste0("ZT", rep(seq(0,20,4),3))) +
      scale_y_continuous(limits = c(min.limit - (max.limit - min.limit)/12.5, 
                                    max.limit)) +
      labs(title = gene.id,
           y = "FPKM",
           x = "") +
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
          size = 10, angle = 90, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.text.y = element_text(
          size = 12, angle = 0, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.title.y = element_text(
          size = 12, angle = 90, vjust = 2.5, hjust = 0.5,
          face = "bold", color = "black"
        )) +
      geom_polygon(
        data = df_polygon_ld,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.1
      ) +
      geom_polygon(
        data = df_polygon_sd,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.2
      ) +
      scale_fill_manual(values = c(rep(c("white", color.1), 3), rep(c("white", color.2),3)))
  }
  else
  {
    ggplot() +
      geom_ribbon(data = wave.ld.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.1, alpha = 0.2) +
      geom_line(data = wave.ld.ggplot, aes(x = x, y = y), color = color.1, size = 1) +
      
      geom_ribbon(data = wave.sd.ggplot, aes(x = x, ymin = min.bound, ymax = max.bound),
                  fill = color.2, alpha = 0.2) +
      geom_line(data = wave.sd.ggplot, aes(x = x, y = y), color = color.2, size = 1) +
      scale_x_continuous(
        breaks = seq(from=0,to=116,by=4),                
        labels = c(paste0("ZT", rep(seq(0,20,4),3)), paste0("ZT", rep(seq(0,20,4),2)))) +
      scale_y_continuous(limits = c(min.limit - (max.limit - min.limit)/12.5, 
                                    max.limit)) +
      labs(title = gene.id,
           y = "FPKM",
           x = "") +
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
          size = 10, angle = 90, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.text.y = element_text(
          size = 12, angle = 0, vjust = 0.5, hjust = 1,
          face = "bold", color = "black"),
        axis.title.y = element_text(
          size = 12, angle = 90, vjust = 2.5, hjust = 0.5,
          face = "bold", color = "black"
        )) +
      geom_polygon(
        data = df_polygon_ld,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.1
      ) +
      geom_polygon(
        data = df_polygon_sd,
        aes(x, y, fill = group, group = group),
        inherit.aes = FALSE,
        color = color.2
      ) +
      scale_fill_manual(values = c(rep(c("white", color.1), 3), rep(c("white", color.2),3)))
  }
  
}
