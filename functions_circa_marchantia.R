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
  points(x = seq(from=0,by=4,to=68),current.gene.expression.sd,lwd=3,col="red", pch=16)
  points(x = seq(from=72,by=4,to=116),current.gene.expression.ll,lwd=3,col="salmon", pch=16)
  
  wave.sd <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["sd.mesor.estimate"]+wave_parameters["sd.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["sd.peak.time.hours"]*pi/12)))
  lines(x = seq(from=0,by=0.1,to=68),wave.sd,lwd=3,col="red")
  
  wave.ll <- sapply(seq(from=72,by=0.1,to=116), function(x) 
    wave_parameters["ll.mesor.estimate"]+wave_parameters["ll.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["ll.peak.time.hours"]*pi/12)))
  lines(x = seq(from=72,by=0.1,to=116),wave.ll,lwd=3,col="salmon")
  
  
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
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red", col="lightpink")
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(96, 120, 120, 96),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="red", col="lightpink")
  
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
  points(x = seq(from=0,by=4,to=68),current.gene.expression.sd,lwd=3,col="tomato3", pch=16)
  points(x = seq(from=72,by=4,to=116),current.gene.expression.dd,lwd=3,col="red4", pch=16)
  
  wave.sd <- sapply(seq(from=0,by=0.1,to=68), function(x) 
    wave_parameters["sd.mesor.estimate"]+wave_parameters["sd.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["sd.peak.time.hours"]*pi/12)))
  lines(x = seq(from=0,by=0.1,to=68),wave.sd,lwd=3,col="tomato3")
  
  wave.dd <- sapply(seq(from=72,by=0.1,to=116), function(x) 
    wave_parameters["dd.mesor.estimate"]+wave_parameters["dd.amplitude.estimate"]*cos(period*(x*pi/12-wave_parameters["dd.peak.time.hours"]*pi/12)))
  lines(x = seq(from=72,by=0.1,to=116),wave.dd,lwd=3,col="red4")
  
  
  for(i in 0:2)
  {
    current.line <- (max.expr- min.expr)/10
    
    polygon(x = c(24*i, 24*i+8, 24*i+8, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="tomato3")
    polygon(x = c(24*i+8,24*(i+1),24*(i+1),24*i+8),
            y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
                min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="tomato3",col="tomato3")
  }
  
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(72, 96, 96, 72),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="tomato3", col="red4")
  
  current.line <- (max.expr- min.expr)/10
  
  polygon(x = c(96, 120, 120, 96),
          y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+(max.expr- min.expr)/40),
              min.expr-(current.line+(max.expr- min.expr)/40)), lwd=2,border="tomato3", col="red4")
  
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


