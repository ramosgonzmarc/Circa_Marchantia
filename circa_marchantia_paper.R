.libPaths("/home/marcos/R/x86_64-pc-linux-gnu-library/4.1")

############################ ANALYSIS ###################################
# Read data
genes.expression <- read.table(file ="gene_expression_marchantia.tsv",sep="\t",header=T,as.is=T)
dim(genes.expression)
head(genes.expression)
genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)

# Number of non-expressed genes and low-expressed genes
length(which(apply(X = genes.expression,MARGIN = 1,FUN = max) == 0))
length(which(apply(X = genes.expression,MARGIN = 1,FUN = max) < 1))
n_exp_genes <- number.genes - length(which(apply(X = genes.expression,MARGIN = 1,FUN = max) < 2))
length(which(apply(X = genes.expression,MARGIN = 1,FUN = max) < 5))

# GO enrichment for low-expressed genes
library("clusterProfiler")

#Install our annotation package for Marchantia polymorpha
#install.packages("./org.Mpolymorpha.eg.db/", repos=NULL, type="source")
library("org.Mpolymorpha71.eg.db")
columns(org.Mpolymorpha71.eg.db)

low.exp.genes <- names(which(apply(X = genes.expression,MARGIN = 1,FUN = max) < 2))
ego.low.exp.genes <- enrichGO(gene = low.exp.genes,
                        OrgDb         = org.Mpolymorpha71.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = F,
                        keyType = "GID")

as.data.frame(ego.low.exp.genes)$Description
barplot(ego.low.exp.genes, showCategory = 20)

# Gene expression matrix for samples under light-dark cycles
ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ld.zt.i <- sapply(X = ld.zt,FUN = function(x){ paste(x,1:3,sep="_")})
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt.i <- sapply(X = sd.zt,FUN = function(x){ paste(x,1:3,sep="_")})

# Free running conditions
dd.ld.zt <- paste("dd_ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
dd.ld.zt.i <- sapply(X = dd.ld.zt,FUN = function(x){ paste(x,1:2,sep="_")})
dd.sd.zt <- paste("dd_sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
dd.sd.zt.i <- sapply(X = dd.sd.zt,FUN = function(x){ paste(x,1:2,sep="_")})

ll.ld.zt <- paste("ll_ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ll.ld.zt.i <- sapply(X = ll.ld.zt,FUN = function(x){ paste(x,1:2,sep="_")})
ll.sd.zt <- paste("ll_sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ll.sd.zt.i <- sapply(X = ll.sd.zt,FUN = function(x){ paste(x,1:2,sep="_")})

# Define general variables
ld.gene.expression <- genes.expression[,ld.zt.i]
sd.gene.expression <- genes.expression[,sd.zt.i]
ld.sd.gene.expression <- genes.expression[,c(ld.zt.i,sd.zt.i)]

ld.dd.gene.expression <- genes.expression[,dd.ld.zt.i]
sd.dd.gene.expression <- genes.expression[,dd.sd.zt.i]
ld.ll.gene.expression <- genes.expression[,ll.ld.zt.i]
sd.ll.gene.expression <- genes.expression[,ll.sd.zt.i]


######### Statistical analysis of circadian genes with CircaCompare  ###########

# First of all, define new expressed-genes background set
bg_set <- names(which(apply(X = genes.expression,MARGIN = 1,FUN = max) >= 2))
bg_table <- genes.expression[bg_set,]
dim(bg_table)

# Determine rhythmic genes in LD 
# With 24 h period
library(circacompare)

results.ld <- matrix(nrow=length(bg_set),ncol=6)
rownames(results.ld) <- bg_set
dim(results.ld)

for(i in 1:length(bg_set))
{
  gene.i <- bg_set[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"), paste(ld.zt,3,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  ld.vector <- c(ld.expression.i)
  
  if ((length(unique(ld.vector[1:12])) > 2))
  {
    ld.df <- data.frame(time=time.points,
                        measure=ld.expression.i)
    
    out.i <- circa_single(x = ld.df, col_time = "time", col_outcome = "measure", period = 24,
                 alpha_threshold = 1)
    results.ld[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.ld[i,] <- NA
  }
}

colnames(results.ld) <- out.i[[2]][,1]
write.table(results.ld,file="circacompare_ld_24.tsv", sep="\t")

# With 12 h period
library(circacompare)

results.ld <- matrix(nrow=length(bg_set),ncol=6)
rownames(results.ld) <- bg_set
dim(results.ld)

for(i in 1:length(bg_set))
{
  gene.i <- bg_set[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"), paste(ld.zt,3,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  ld.vector <- c(ld.expression.i)
  
  if ((length(unique(ld.vector[1:12])) > 2))
  {
    ld.df <- data.frame(time=time.points,
                        measure=ld.expression.i)
    
    out.i <- circa_single(x = ld.df, col_time = "time", col_outcome = "measure", period = 12,
                          alpha_threshold = 1)
    results.ld[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.ld[i,] <- NA
  }
}

colnames(results.ld) <- out.i[[2]][,1]
write.table(results.ld,file="circacompare_ld_12.tsv", sep="\t")

# Determine rhythmic genes in SD 
# With 24 h period
library(circacompare)

results.sd <- matrix(nrow=length(bg_set),ncol=6)
rownames(results.sd) <- bg_set
dim(results.sd)

for(i in 1:length(bg_set))
{
  gene.i <- bg_set[i]
  
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"), paste(sd.zt,3,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  sd.vector <- c(sd.expression.i)
  
  if ((length(unique(sd.vector[1:12])) > 2))
  {
    sd.df <- data.frame(time=time.points,
                        measure=sd.expression.i)
    
    out.i <- circa_single(x = sd.df, col_time = "time", col_outcome = "measure", period = 24,
                          alpha_threshold = 1)
    results.sd[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.sd[i,] <- NA
  }
}

colnames(results.sd) <- out.i[[2]][,1]
write.table(results.sd,file="circacompare_sd_24.tsv", sep="\t")

# With 12 h period
library(circacompare)

results.sd <- matrix(nrow=length(bg_set),ncol=6)
rownames(results.sd) <- bg_set
dim(results.sd)

for(i in 1:length(bg_set))
{
  gene.i <- bg_set[i]
  
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"), paste(sd.zt,3,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  sd.vector <- c(sd.expression.i)
  
  if ((length(unique(sd.vector[1:12])) > 2))
  {
    sd.df <- data.frame(time=time.points,
                        measure=sd.expression.i)
    
    out.i <- circa_single(x = sd.df, col_time = "time", col_outcome = "measure", period = 12,
                          alpha_threshold = 1)
    results.sd[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.sd[i,] <- NA
  }
}

colnames(results.sd) <- out.i[[2]][,1]
write.table(results.sd,file="circacompare_sd_12.tsv", sep="\t")

# Rhythmicity under light/dark cycles
#LD
results.ld <- read.table("tables/circacompare_ld_24.tsv", header = T, sep = "\t")
sum(results.ld$rhythmic_p < 0.05, na.rm = T)/length(bg_set)
rhythmic.genes.ld <- rownames(subset(results.ld, rhythmic_p < 0.05))
length(rhythmic.genes.ld)

results.ld.12 <- read.table("tables/circacompare_ld_12.tsv", header = T, sep = "\t")
sum(results.ld.12$rhythmic_p < 0.05, na.rm = T)/length(bg_set)
rhythmic.genes.ld.12 <- rownames(subset(results.ld.12, rhythmic_p < 0.05))
length(rhythmic.genes.ld.12)

complete.ld.rhythmic.genes <- unique(c(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(complete.ld.rhythmic.genes)/length(bg_set)

non.rhythmic.ld.genes <- setdiff(mp.genes,complete.ld.rhythmic.genes)
length(non.rhythmic.ld.genes)

#SD
results.sd <- read.table("tables/circacompare_sd_24.tsv", header = T, sep = "\t")
sum(results.sd$rhythmic_p < 0.05, na.rm = T)/length(bg_set)
rhythmic.genes.sd <- rownames(subset(results.sd, rhythmic_p < 0.05))
length(rhythmic.genes.sd)

results.sd.12 <- read.table("tables/circacompare_sd_12.tsv", header = T, sep = "\t")
sum(results.sd.12$rhythmic_p < 0.05, na.rm = T)/length(bg_set)
rhythmic.genes.sd.12 <- rownames(subset(results.sd.12, rhythmic_p < 0.05))
length(rhythmic.genes.sd.12)

complete.sd.rhythmic.genes <- unique(c(rhythmic.genes.sd,rhythmic.genes.sd.12))
length(complete.sd.rhythmic.genes)/length(bg_set)

non.rhythmic.sd.genes <- setdiff(mp.genes,complete.sd.rhythmic.genes)
length(non.rhythmic.sd.genes)

#Rhythmic genes in SD conditions but not in LD
sd.specific.rhythmic.genes <- setdiff(complete.sd.rhythmic.genes,complete.ld.rhythmic.genes)
length(sd.specific.rhythmic.genes)

#Rhythmic genes in LD conditions but not in SD
ld.specific.rhytmic.genes <-setdiff(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(ld.specific.rhytmic.genes)

#Rhythmic genes in both LD and SD
complete.ld.sd.rhythmic.genes <- intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(complete.ld.sd.rhythmic.genes)

#Save lists of rhythmic genes
ld.rhytmic.table <- matrix(data = complete.ld.rhythmic.genes,ncol = 1)
write.table(ld.rhytmic.table,file="tables/ld.rhythmic.table.tsv")
sd.rhytmic.table <- matrix(data = complete.sd.rhythmic.genes,ncol = 1)
write.table(sd.rhytmic.table,file="tables/sd.rhythmic.table.tsv")

#### Plots of rhythmic genes

#Venn diagram of LD and SD rhythmic genes
library(VennDiagram)
png("ld_sd_venn.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(complete.ld.rhythmic.genes), area2 = length(complete.sd.rhythmic.genes),
                   cross.area = length(intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

#Barplots of rhythmic and non-rhythmic genes in each photoperiod
rhythmic.bars <- matrix(c(length(rhythmic.genes.ld.12),length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12)),
                          length(non.rhythmic.ld.genes),length(rhythmic.genes.sd.12),
                          length(setdiff(rhythmic.genes.sd,rhythmic.genes.sd.12)),length(non.rhythmic.sd.genes)),
                        nrow=3)

png("rhythm_by_period.png", width = 1500, height = 1700, res=300)
par(lwd=3)
barplot(rhythmic.bars, lwd=2, names.arg = c("LD","SD"),cex.names = 2,cex.axis = 1.5, width = 0.35, xlim = c(0,2),
        col=colors()[c(31,89,12)], border="black", space=0.25, font.axis=2, xlab="", ylim=c(0,20000), 
        legend.text = c("Rhythmic (12h)","Rhythmic (24h)","Non rhythmic"),
        main = "Rhythmic genes per photoperiod")
dev.off()


#### Expression level in rhythmic and non rhythmic genes
rhythmic.max.expression <- apply(X = genes.expression[complete.ld.sd.rhythmic.genes,
                                                      c(ld.zt.i,sd.zt.i)],MARGIN = 1,FUN = max)
non.rhythmic.max.expression <- apply(X = genes.expression[intersect(non.rhythmic.ld.genes,non.rhythmic.sd.genes),
                                                          c(ld.zt.i,sd.zt.i)], MARGIN = 1,FUN = max)

png("expression_vs_rhythmicity.png", width = 1500, height = 1500, res=300)
par(lwd=3)
boxplot(rhythmic.max.expression, non.rhythmic.max.expression, outline=F,col=c("chartreuse4","royalblue3"),
        names=c("Rhythmic","Non rhythmic"),ylab="FPKM",cex.lab=1.5,las=1, cex.axis=1.5, main= "Distribution of expression levels")
dev.off()

wilcox.test(rhythmic.max.expression,non.rhythmic.max.expression,alternative="greater")

##### Difference between LD and SD in a parametric manner, complete set
library(circacompare)

results.ld.sd <- matrix(nrow=length(bg_set),ncol=15)
rownames(results.ld.sd) <- bg_set
dim(results.ld.sd)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

for(i in 1:length(bg_set))
{
  gene.i <- bg_set[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"), paste(ld.zt,3,sep="_"))]
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"), paste(sd.zt,3,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  ld.sd.vector <- c(ld.expression.i, sd.expression.i)
  
  if ((length(unique(ld.sd.vector[1:12])) > 2) && (length(unique(ld.sd.vector[19:30])) > 2))
  {
    ld.sd.df <- data.frame(time=c(time.points,time.points),
                           measure=c(ld.expression.i, sd.expression.i),
                           group=c(rep("ld",18),rep("sd",18)))
    
    out.i <- circacompare(x = ld.sd.df, col_time = "time", col_group = "group",
                          col_outcome = "measure",alpha_threshold = 1)
    results.ld.sd[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.ld.sd[i,] <- NA
  }
}

colnames(results.ld.sd) <- out.i[[2]][,1]
write.table(results.ld.sd,file="circacompare_ld_sd_complete.tsv", sep="\t")

##### Difference between LD and SD in a parametric manner, rhythmic in both
library(circacompare)

intersect_set <- intersect(complete.ld.rhythmic.genes, complete.sd.rhythmic.genes)
intersect_table <- genes.expression[intersect_set,]
dim(intersect_table)

results.ld.sd.intersect <- matrix(nrow=length(intersect_set),ncol=15)
rownames(results.ld.sd.intersect) <- intersect_set
dim(results.ld.sd.intersect)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")


#Construct a table with the circadian data for all the genes rhythmic in both LD and SD conditions
for(i in 1:length(intersect_set))
{
  gene.i <- intersect_set[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"), paste(ld.zt,3,sep="_"))]
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"), paste(sd.zt,3,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  ld.sd.vector <- c(ld.expression.i, sd.expression.i)
  
  if ((length(unique(ld.sd.vector[1:12])) > 2) && (length(unique(ld.sd.vector[19:30])) > 2))
  {
    ld.sd.df <- data.frame(time=c(time.points,time.points),
                           measure=c(ld.expression.i, sd.expression.i),
                           group=c(rep("ld",18),rep("sd",18)))
    
    out.i <- circacompare(x = ld.sd.df, col_time = "time", col_group = "group",
                          col_outcome = "measure",alpha_threshold = 1)
    results.ld.sd.intersect[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.ld.sd.intersect[i,] <- NA
  }
}

colnames(results.ld.sd.intersect) <- out.i[[2]][,1]
write.table(results.ld.sd.intersect,file="circacompare_ld_sd_intersection.tsv", sep="\t")


################# Free-running ######################33
# LD - LL 24 h
# Only for rhythmic genes in LD
complete.ld.rhythmic.genes <- read.table("tables/ld.rhythmic.table.tsv", header = T)$V1
head(complete.ld.rhythmic.genes)
ld_set <- complete.ld.rhythmic.genes
ld_set_table <- genes.expression[ld_set,]
dim(ld_set_table)

results.ld.ll <- matrix(nrow=length(ld_set),ncol=6)
rownames(results.ld.ll) <- ld_set
dim(results.ld.ll)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ll.ld.zt <- paste("ll_ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

#Construct a table with the circadian data for all the genes rhythmic in both LD and SD conditions
for(i in 1:length(ld_set))
{
  gene.i <- ld_set[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,2,sep="_"))]
  ll.expression.i <- genes.expression[gene.i,c(paste(ll.ld.zt,1,sep="_"),paste(ll.ld.zt,2,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  ld.ll.vector <- c(ld.expression.i, ll.expression.i)
  
  if ((length(unique(ld.ll.vector[1:12])) > 6) & (length(unique(ll.expression.i)) > 3))
  {
    ld.ll.df <- data.frame(time=time.points,
                        measure=ld.ll.vector)
    
    out.i <- circa_single(x = ld.ll.df, col_time = "time", col_outcome = "measure", period = 24,
                          alpha_threshold = 1)
    results.ld.ll[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.ld.ll[i,] <- NA
  }
}
colnames(results.ld.ll) <- out.i[[2]][,1]
write.table(results.ld.ll,file="circacompare_ld_ll.tsv", sep="\t")

# LD - LL 12 h
# Only for rhythmic genes in LD
ld_set <- complete.ld.rhythmic.genes
ld_set_table <- genes.expression[ld_set,]
dim(ld_set_table)

results.ld.ll <- matrix(nrow=length(ld_set),ncol=6)
rownames(results.ld.ll) <- ld_set
dim(results.ld.ll)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ll.ld.zt <- paste("ll_ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

#Construct a table with the circadian data for all the genes rhythmic in both LD and SD conditions
for(i in 1:length(ld_set))
{
  gene.i <- ld_set[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,2,sep="_"))]
  ll.expression.i <- genes.expression[gene.i,c(paste(ll.ld.zt,1,sep="_"),paste(ll.ld.zt,2,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  ld.ll.vector <- c(ld.expression.i, ll.expression.i)
  
  if ((length(unique(ld.ll.vector[1:12])) > 6) & (length(unique(ll.expression.i)) > 3))
  {
    ld.ll.df <- data.frame(time=time.points,
                           measure=ld.ll.vector)
    
    out.i <- circa_single(x = ld.ll.df, col_time = "time", col_outcome = "measure", period = 12,
                          alpha_threshold = 1)
    results.ld.ll[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.ld.ll[i,] <- NA
  }
}
colnames(results.ld.ll) <- out.i[[2]][,1]
write.table(results.ld.ll,file="circacompare_ld_ll_12.tsv", sep="\t")

# LD - DD 24 h
# Only for rhythmic genes in LD
ld_set <- complete.ld.rhythmic.genes
ld_set_table <- genes.expression[ld_set,]
dim(ld_set_table)

results.ld.dd <- matrix(nrow=length(ld_set),ncol=6)
rownames(results.ld.dd) <- ld_set
dim(results.ld.dd)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
dd.ld.zt <- paste("dd_ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

for(i in 1:length(ld_set))
{
  gene.i <- ld_set[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,2,sep="_"))]
  dd.expression.i <- genes.expression[gene.i,c(paste(dd.ld.zt,1,sep="_"),paste(dd.ld.zt,2,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  ld.dd.vector <- c(ld.expression.i, dd.expression.i)
  
  if ((length(unique(ld.dd.vector[1:12])) > 6) & (length(unique(dd.expression.i)) > 3))
  {
    ld.dd.df <- data.frame(time=time.points,
                           measure=ld.dd.vector)
    
    out.i <- circa_single(x = ld.dd.df, col_time = "time", col_outcome = "measure", period = 24,
                          alpha_threshold = 1)
    results.ld.dd[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.ld.dd[i,] <- NA
  }
}
colnames(results.ld.dd) <- out.i[[2]][,1]
write.table(results.ld.dd,file="circacompare_ld_dd.tsv", sep="\t")

# LD - DD 12 h
# Only for rhythmic genes in LD
ld_set <- complete.ld.rhythmic.genes
ld_set_table <- genes.expression[ld_set,]
dim(ld_set_table)

results.ld.dd <- matrix(nrow=length(ld_set),ncol=6)
rownames(results.ld.dd) <- ld_set
dim(results.ld.dd)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
dd.ld.zt <- paste("dd_ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

for(i in 1:length(ld_set))
{
  gene.i <- ld_set[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,2,sep="_"))]
  dd.expression.i <- genes.expression[gene.i,c(paste(dd.ld.zt,1,sep="_"),paste(dd.ld.zt,2,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  ld.dd.vector <- c(ld.expression.i, dd.expression.i)
  
  if ((length(unique(ld.dd.vector[1:12])) > 6) & (length(unique(dd.expression.i)) > 3))
  {
    ld.dd.df <- data.frame(time=time.points,
                           measure=ld.dd.vector)
    
    out.i <- circa_single(x = ld.dd.df, col_time = "time", col_outcome = "measure", period = 12,
                          alpha_threshold = 1)
    results.ld.dd[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.ld.dd[i,] <- NA
  }
}
colnames(results.ld.dd) <- out.i[[2]][,1]
write.table(results.ld.dd,file="circacompare_ld_dd_12.tsv", sep="\t")

# SD - LL 24 h
# Only for rhythmic genes in SD
complete.sd.rhythmic.genes <- read.table("tables/sd.rhythmic.table.tsv", header = T)$V1
head(complete.sd.rhythmic.genes)
sd_set <- complete.sd.rhythmic.genes
sd_set_table <- genes.expression[sd_set,]
dim(sd_set_table)

results.sd.ll <- matrix(nrow=length(sd_set),ncol=6)
rownames(results.sd.ll) <- sd_set
dim(results.sd.ll)

sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ll.sd.zt <- paste("ll_sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

for(i in 1:length(sd_set))
{
  gene.i <- sd_set[i]
  
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,2,sep="_"))]
  ll.expression.i <- genes.expression[gene.i,c(paste(ll.sd.zt,1,sep="_"),paste(ll.sd.zt,2,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  sd.ll.vector <- c(sd.expression.i, ll.expression.i)
  
  if ((length(unique(sd.ll.vector[1:12])) > 6) & (length(unique(ll.expression.i)) > 3))
  {
    sd.ll.df <- data.frame(time=time.points,
                           measure=sd.ll.vector)
    
    out.i <- circa_single(x = sd.ll.df, col_time = "time", col_outcome = "measure", period = 24,
                          alpha_threshold = 1)
    results.sd.ll[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.sd.ll[i,] <- NA
  }
}
colnames(results.sd.ll) <- out.i[[2]][,1]
write.table(results.sd.ll,file="circacompare_sd_ll.tsv", sep="\t")

# SD - LL 12 h
# Only for rhythmic genes in SD
sd_set <- complete.sd.rhythmic.genes
sd_set_table <- genes.expression[sd_set,]
dim(sd_set_table)

results.sd.ll <- matrix(nrow=length(sd_set),ncol=6)
rownames(results.sd.ll) <- sd_set
dim(results.sd.ll)

sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ll.sd.zt <- paste("ll_sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

for(i in 1:length(sd_set))
{
  gene.i <- sd_set[i]
  
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,2,sep="_"))]
  ll.expression.i <- genes.expression[gene.i,c(paste(ll.sd.zt,1,sep="_"),paste(ll.sd.zt,2,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  sd.ll.vector <- c(sd.expression.i, ll.expression.i)
  
  if ((length(unique(sd.ll.vector[1:12])) > 6) & (length(unique(ll.expression.i)) > 3))
  {
    sd.ll.df <- data.frame(time=time.points,
                           measure=sd.ll.vector)
    
    out.i <- circa_single(x = sd.ll.df, col_time = "time", col_outcome = "measure", period = 12,
                          alpha_threshold = 1)
    results.sd.ll[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.sd.ll[i,] <- NA
  }
}
colnames(results.sd.ll) <- out.i[[2]][,1]
write.table(results.sd.ll,file="circacompare_sd_ll_12.tsv", sep="\t")

# SD - DD 24 h
# Only for rhythmic genes in SD
sd_set <- complete.sd.rhythmic.genes
sd_set_table <- genes.expression[sd_set,]
dim(sd_set_table)

results.sd.dd <- matrix(nrow=length(sd_set),ncol=6)
rownames(results.sd.dd) <- sd_set
dim(results.sd.dd)

sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
dd.sd.zt <- paste("dd_sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

for(i in 1:length(sd_set))
{
  gene.i <- sd_set[i]
  
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,2,sep="_"))]
  dd.expression.i <- genes.expression[gene.i,c(paste(dd.sd.zt,1,sep="_"),paste(dd.sd.zt,2,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  sd.dd.vector <- c(sd.expression.i, dd.expression.i)
  
  if ((length(unique(sd.dd.vector[1:12])) > 6) & (length(unique(dd.expression.i)) > 3))
  {
    sd.dd.df <- data.frame(time=time.points,
                           measure=sd.dd.vector)
    
    out.i <- circa_single(x = sd.dd.df, col_time = "time", col_outcome = "measure", period = 24,
                          alpha_threshold = 1)
    results.sd.dd[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.sd.dd[i,] <- NA
  }
}
colnames(results.sd.dd) <- out.i[[2]][,1]
write.table(results.sd.dd,file="circacompare_sd_dd.tsv", sep="\t")

# SD - DD 12 h
# Only for rhythmic genes in SD
sd_set <- complete.sd.rhythmic.genes
sd_set_table <- genes.expression[sd_set,]
dim(sd_set_table)

results.sd.dd <- matrix(nrow=length(sd_set),ncol=6)
rownames(results.sd.dd) <- sd_set
dim(results.sd.dd)

sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
dd.sd.zt <- paste("dd_sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

for(i in 1:length(sd_set))
{
  gene.i <- sd_set[i]
  
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,2,sep="_"))]
  dd.expression.i <- genes.expression[gene.i,c(paste(dd.sd.zt,1,sep="_"),paste(dd.sd.zt,2,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  sd.dd.vector <- c(sd.expression.i, dd.expression.i)
  
  if ((length(unique(sd.dd.vector[1:12])) > 6) & (length(unique(dd.expression.i)) > 3))
  {
    sd.dd.df <- data.frame(time=time.points,
                           measure=sd.dd.vector)
    
    out.i <- circa_single(x = sd.dd.df, col_time = "time", col_outcome = "measure", period = 12,
                          alpha_threshold = 1)
    results.sd.dd[i,] <- out.i[[2]][,2]
  }
  else
  {
    results.sd.dd[i,] <- NA
  }
}
colnames(results.sd.dd) <- out.i[[2]][,1]
write.table(results.sd.dd,file="circacompare_sd_dd_12.tsv", sep="\t")

# Exploration of free-running results
# LD - LL
results.ld.ll.24 <- read.table("tables/circacompare_ld_ll.tsv", header = T, sep = "\t")
sum(results.ld.ll.24$rhythmic_p < 0.05, na.rm = T)/length(bg_set)
rhythmic.genes.ld.ll.24 <- rownames(subset(results.ld.ll.24, rhythmic_p < 0.05))
length(rhythmic.genes.ld.ll.24)

results.ld.ll.12 <- read.table("tables/circacompare_ld_ll_12.tsv", header = T, sep = "\t")
sum(results.ld.ll.12$rhythmic_p < 0.05, na.rm = T)/length(bg_set)
rhythmic.genes.ld.ll.12 <- rownames(subset(results.ld.ll.12, rhythmic_p < 0.05))
length(rhythmic.genes.ld.ll.12)

complete.ld.ll.rhythmic.genes <- unique(c(rhythmic.genes.ld.ll.24,rhythmic.genes.ld.ll.12))
length(complete.ld.ll.rhythmic.genes)/length(bg_set)

write.table(complete.ld.ll.rhythmic.genes,file="tables/ld.ll.rhythmic.table.tsv", sep="\t")

# LD - DD
results.ld.dd.24 <- read.table("tables/circacompare_ld_dd.tsv", header = T, sep = "\t")
sum(results.ld.dd.24$rhythmic_p < 0.05, na.rm = T)/length(bg_set)
rhythmic.genes.ld.dd.24 <- rownames(subset(results.ld.dd.24, rhythmic_p < 0.05))
length(rhythmic.genes.ld.dd.24)

results.ld.dd.12 <- read.table("tables/circacompare_ld_dd_12.tsv", header = T, sep = "\t")
sum(results.ld.dd.12$rhythmic_p < 0.05, na.rm = T)/length(bg_set)
rhythmic.genes.ld.dd.12 <- rownames(subset(results.ld.dd.12, rhythmic_p < 0.05))
length(rhythmic.genes.ld.dd.12)

complete.ld.dd.rhythmic.genes <- unique(c(rhythmic.genes.ld.dd.24,rhythmic.genes.ld.dd.12))
length(complete.ld.dd.rhythmic.genes)/length(bg_set)

write.table(complete.ld.dd.rhythmic.genes,file="tables/ld.dd.rhythmic.table.tsv", sep="\t")

# SD - LL
results.sd.ll.24 <- read.table("tables/circacompare_sd_ll.tsv", header = T, sep = "\t")
sum(results.sd.ll.24$rhythmic_p < 0.05, na.rm = T)/length(bg_set)
rhythmic.genes.sd.ll.24 <- rownames(subset(results.sd.ll.24, rhythmic_p < 0.05))
length(rhythmic.genes.sd.ll.24)

results.sd.ll.12 <- read.table("tables/circacompare_sd_ll_12.tsv", header = T, sep = "\t")
sum(results.sd.ll.12$rhythmic_p < 0.05, na.rm = T)/length(bg_set)
rhythmic.genes.sd.ll.12 <- rownames(subset(results.sd.ll.12, rhythmic_p < 0.05))
length(rhythmic.genes.sd.ll.12)

complete.sd.ll.rhythmic.genes <- unique(c(rhythmic.genes.sd.ll.24,rhythmic.genes.sd.ll.12))
length(complete.sd.ll.rhythmic.genes)/length(bg_set)

write.table(complete.sd.ll.rhythmic.genes,file="tables/sd.ll.rhythmic.table.tsv", sep="\t")

# SD - DD
results.sd.dd.24 <- read.table("tables/circacompare_sd_dd.tsv", header = T, sep = "\t")
sum(results.sd.dd.24$rhythmic_p < 0.05, na.rm = T)/length(bg_set)
rhythmic.genes.sd.dd.24 <- rownames(subset(results.sd.dd.24, rhythmic_p < 0.05))
length(rhythmic.genes.sd.dd.24)

results.sd.dd.12 <- read.table("tables/circacompare_sd_dd_12.tsv", header = T, sep = "\t")
sum(results.sd.dd.12$rhythmic_p < 0.05, na.rm = T)/length(bg_set)
rhythmic.genes.sd.dd.12 <- rownames(subset(results.sd.dd.12, rhythmic_p < 0.05))
length(rhythmic.genes.sd.dd.12)

complete.sd.dd.rhythmic.genes <- unique(c(rhythmic.genes.sd.dd.24,rhythmic.genes.sd.dd.12))
length(complete.sd.dd.rhythmic.genes)/length(bg_set)

write.table(complete.sd.dd.rhythmic.genes,file="tables/sd.dd.rhythmic.table.tsv", sep="\t")

# Create barplots of each type of retained rhythmicity for LD transitions
bg_set <- names(which(apply(X = genes.expression,MARGIN = 1,FUN = max) >= 2))
intersect(names(which(apply(X = genes.expression[, ld.zt.i],MARGIN = 1,FUN = max) < 2 )), complete.ld.rhythmic.genes)

# We are considering the same non-expressed genes for LD and SD, as those that exhibit < 2 FPKM on all 84 samples
not.expressed.ld <- names(which(apply(X = genes.expression,MARGIN = 1,FUN = max) < 2))
non.rhythmic.ld.genes <- setdiff(mp.genes,complete.ld.rhythmic.genes)
non.rhythmic.ld.genes.but.expressed <- setdiff(non.rhythmic.ld.genes, not.expressed.ld)
rhythmic.ld.ll.dd <- intersect(complete.ld.rhythmic.genes, intersect(complete.ld.ll.rhythmic.genes,
                              complete.ld.dd.rhythmic.genes))
rhythmic.ld.ll.no.dd <- setdiff(complete.ld.ll.rhythmic.genes, rhythmic.ld.ll.dd)
rhythmic.ld.dd.no.ll <- setdiff(complete.ld.dd.rhythmic.genes, rhythmic.ld.ll.dd)
rhythmic.only.ld <- setdiff(setdiff(complete.ld.rhythmic.genes,complete.ld.ll.rhythmic.genes),
                           complete.ld.dd.rhythmic.genes)

data.ld.ll <- matrix(c( length(rhythmic.ld.ll.dd),
                        length(rhythmic.ld.dd.no.ll),
                        length(rhythmic.ld.ll.no.dd),
                        length(rhythmic.only.ld),
                        length(non.rhythmic.ld.genes.but.expressed),
                        length(not.expressed.ld)
),nrow=6)


png("stack_ld.png", width = 1500, height = 1500, res=300)
par(lwd=3)
barplot(data.ld.ll, lwd=3, names.arg = "",cex.names = 2,cex.axis = 1.5,
        col=c("blue", "darkblue", "lightskyblue", "lightcyan2", "white", "gray89"), 
        border="black", width = 0.2, xlim = c(0,1),
        font.axis=2, main = "LD rhythmic genes", 
        legend.text = c("Rhythmic LD/LL/DD", "Rhythmic LD/DD", "Rhythmic LD/LL", 
                        "Rhythmic only LD", "Non rhythmic LD", "Not expressed"),
        xlab="",ylim=c(0,20000))
dev.off()

# Version without non-expressed genes
png("stack_ld_only_rhythmic.png", width = 1500, height = 1500, res=300)
par(lwd=3)
barplot(as.matrix(data.ld.ll[1:5,]), lwd=3, names.arg = "",cex.names = 2,cex.axis = 1.5,
        col=c("blue", "darkblue","lightskyblue", "lightcyan2","white"), 
        border="black", width = 0.2, xlim = c(0,1),
        font.axis=2, main = "LD rhythmic genes", 
        legend.text = c("Rhythmic LD/LL/DD", "Rhythmic LD/DD", "Rhythmic LD/LL", 
                        "Rhythmic only LD", "Non rhythmic LD"),
        xlab="",ylim=c(0,13000))
dev.off()

# Venn for LD DD and LL
library(VennDiagram)
png("ld_dd_ll_venn.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(complete.ld.dd.rhythmic.genes), area2 = length(complete.ld.ll.rhythmic.genes),
                   cross.area = length(intersect(complete.ld.dd.rhythmic.genes,complete.ld.ll.rhythmic.genes)),
                   lwd = 3,category = c("LD-DD","LD-LL"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# Create barplots of each type of retained rhythmicity for SD transitions
not.expressed.sd <- names(which(apply(X = genes.expression,MARGIN = 1,FUN = max) < 2))
non.rhythmic.sd.genes <- setdiff(mp.genes,complete.sd.rhythmic.genes)
non.rhythmic.sd.genes.but.expressed <- setdiff(non.rhythmic.sd.genes, not.expressed.sd)
rhythmic.sd.ll.dd <- intersect(complete.sd.rhythmic.genes, intersect(complete.sd.ll.rhythmic.genes,
                                                                     complete.sd.dd.rhythmic.genes))
rhythmic.sd.ll.no.dd <- setdiff(complete.sd.ll.rhythmic.genes, rhythmic.sd.ll.dd)
rhythmic.sd.dd.no.ll <- setdiff(complete.sd.dd.rhythmic.genes, rhythmic.sd.ll.dd)
rhythmic.only.sd <- setdiff(setdiff(complete.sd.rhythmic.genes,complete.sd.ll.rhythmic.genes),
                            complete.sd.dd.rhythmic.genes)

data.sd.ll <- matrix(c( length(rhythmic.sd.ll.dd),
                        length(rhythmic.sd.dd.no.ll),
                        length(rhythmic.sd.ll.no.dd),
                        length(rhythmic.only.sd),
                        length(non.rhythmic.sd.genes.but.expressed),
                        length(not.expressed.sd)
),nrow=6)

png("stack_sd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
barplot(data.sd.ll, lwd=3, names.arg = "",cex.names = 2,cex.axis = 1.5,
        col=c("red", "darkred", "salmon", "rosybrown1", "white", "gray89"), 
        border="black", width = 0.2, xlim = c(0,1),
        font.axis=2, main = "SD-LL rhythmic genes", 
        legend.text = c("Rhythmic SD/LL/DD", "Rhythmic SD/DD", "Rhythmic SD/LL", 
                        "Rhythmic only SD", "Non rhythmic SD", "Not expressed"),
        xlab="",ylim=c(0,20000))
dev.off()

# Version without non-expressed genes
png("stack_sd_only_rhythmic.png", width = 1500, height = 1500, res=300)
par(lwd=3)
barplot(as.matrix(data.sd.ll[1:5,]), lwd=3, names.arg = "",cex.names = 2,cex.axis = 1.5,
        col=c("red", "darkred", "salmon", "rosybrown1", "white"), 
        border="black", width = 0.2, xlim = c(0,1),
        font.axis=2, main = "SD-LL rhythmic genes", 
        legend.text = c("Rhythmic SD/LL/DD", "Rhythmic SD/DD", "Rhythmic SD/LL", 
                        "Rhythmic only SD", "Non rhythmic SD"),
        xlab="",ylim=c(0,13000))
dev.off()

# Venn for SD DD and LL
png("sd_dd_ll_venn.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(complete.sd.dd.rhythmic.genes), area2 = length(complete.sd.ll.rhythmic.genes),
                   cross.area = length(intersect(complete.sd.dd.rhythmic.genes,complete.sd.ll.rhythmic.genes)),
                   lwd = 3,category = c("SD-DD","SD-LL"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# Venn for genes that retained rhythmicity in LL under LD and SD
png("ld_sd_ll.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(complete.ld.ll.rhythmic.genes), area2 = length(complete.sd.ll.rhythmic.genes),
                   cross.area = length(intersect(complete.ld.ll.rhythmic.genes,complete.sd.ll.rhythmic.genes)),
                   lwd = 3,category = c("LD_LL","SD_LL"),euler.d = T,col = "black",fill = c("blue","red"),
                   alpha = 0.65,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# Venn for genes that retained rhythmicity in DD under LD and SD
png("ld_sd_dd.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(complete.ld.dd.rhythmic.genes), area2 = length(complete.sd.dd.rhythmic.genes),
                   cross.area = length(intersect(complete.ld.dd.rhythmic.genes,complete.sd.dd.rhythmic.genes)),
                   lwd = 3,category = c("LD_DD","SD_DD"),euler.d = T,col = "black",fill = c("blue","red"),
                   alpha = 0.65,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()


# Bona fide circadian genes
png("bona_fide.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(rhythmic.ld.ll.dd), area2 = length(rhythmic.sd.ll.dd),
                   cross.area = length(intersect(rhythmic.ld.ll.dd,rhythmic.sd.ll.dd)),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = "black",fill = c("blue","red"),
                   alpha = 0.65,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()


#######Change from 1 to 2 peaks#########

only.rhythmic.genes.sd.12 <- setdiff(rhythmic.genes.sd.12, rhythmic.genes.sd)
only.rhythmic.genes.ld.24 <- setdiff(rhythmic.genes.ld, rhythmic.genes.ld.12)
two.peaks.change.genes <- intersect(only.rhythmic.genes.ld.24, only.rhythmic.genes.sd.12)
length(two.peaks.change.genes) # 722 genes changed from one peak in LD to two in SD

only.rhythmic.genes.ld.12 <- setdiff(rhythmic.genes.ld.12, rhythmic.genes.ld)
only.rhythmic.genes.sd.24 <- setdiff(rhythmic.genes.sd, rhythmic.genes.sd.12)
two.peaks.change.genes.reverse <- intersect(only.rhythmic.genes.sd.24, only.rhythmic.genes.ld.12)
length(two.peaks.change.genes.reverse) # 451 genes changed from one peak in SD to two in LD

#### Free running parameters
## Compare LD against LL for estimating parameters
# LD to LL
library(circacompare)

complete.ld.ll.rhythmic.genes <- read.table(file="tables/ld.ll.rhythmic.table.tsv", sep="\t")$x
circacompare.ld.ll <- matrix(nrow=length(complete.ld.ll.rhythmic.genes),ncol=15)
rownames(circacompare.ld.ll) <- complete.ld.ll.rhythmic.genes

for(i in 1:length(complete.ld.ll.rhythmic.genes))
{
  gene.i <- complete.ld.ll.rhythmic.genes[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
  ll.expression.i <- genes.expression[gene.i,c(paste(ll.ld.zt,1,sep="_"),paste(ll.ld.zt,2,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 18)
  
  if ((length(unique(sd.dd.vector[1:12])) > 6) & (length(unique(sd.dd.vector[1:12])) > 3))
  {
    ld.ll.df <- data.frame(time=c(time.points,time.points[1:12]),
                           measure=c(ld.expression.i, ll.expression.i),
                           group=c(rep("ld",18),rep("ll",12)))
    
    out.i <- circacompare(x = ld.ll.df, col_time = "time", col_group = "group",
                          col_outcome = "measure",alpha_threshold = 1)
    circacompare.ld.ll[i,] <- out.i[[2]][,2]
  }
  else
  {
    circacompare.ld.ll[i,] <- NA
  }
  
}

colnames(circacompare.ld.ll) <- out.i[[2]][,1]
write.table(circacompare.ld.ll,file="circacompare_ld_ll.tsv", sep="\t")


# LD to DD
library(circacompare)

complete.ld.dd.rhythmic.genes <- read.table(file="tables/ld.dd.rhythmic.table.tsv", sep="\t")$x
circacompare.ld.dd <- matrix(nrow=length(complete.ld.dd.rhythmic.genes),ncol=15)
rownames(circacompare.ld.dd) <- complete.ld.dd.rhythmic.genes

for(i in 1:length(complete.ld.dd.rhythmic.genes))
{
  gene.i <- complete.ld.dd.rhythmic.genes[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
  dd.expression.i <- genes.expression[gene.i,c(paste(dd.ld.zt,1,sep="_"),paste(dd.ld.zt,2,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 18)
  
  if ((length(unique(sd.dd.vector[1:12])) > 6) & (length(unique(sd.dd.vector[1:12])) > 3))
  {
    ld.dd.df <- data.frame(time=c(time.points,time.points[1:12]),
                           measure=c(ld.expression.i, dd.expression.i),
                           group=c(rep("ld",18),rep("dd",12)))
    
    out.i <- circacompare(x = ld.dd.df, col_time = "time", col_group = "group",
                          col_outcome = "measure",alpha_threshold = 1)
    circacompare.ld.dd[i,] <- out.i[[2]][,2]
  }
  else
  {
    circacompare.ld.dd[i,] <- NA
  }
  
}

colnames(circacompare.ld.dd) <- out.i[[2]][,1]
write.table(circacompare.ld.dd,file="circacompare_ld_dd.tsv", sep="\t")


# SD to LL
library(circacompare)

complete.sd.ll.rhythmic.genes <- read.table(file="tables/sd.ll.rhythmic.table.tsv", sep="\t")$x
circacompare.sd.ll <- matrix(nrow=length(complete.sd.ll.rhythmic.genes),ncol=15)
rownames(circacompare.sd.ll) <- complete.sd.ll.rhythmic.genes

for(i in 1:length(complete.sd.ll.rhythmic.genes))
{
  gene.i <- complete.sd.ll.rhythmic.genes[i]
  
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
  ll.expression.i <- genes.expression[gene.i,c(paste(ll.sd.zt,1,sep="_"),paste(ll.sd.zt,2,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 18)
  
  if ((length(unique(sd.dd.vector[1:12])) > 6) & (length(unique(sd.dd.vector[1:12])) > 3))
  {
    sd.ll.df <- data.frame(time=c(time.points,time.points[1:12]),
                           measure=c(sd.expression.i, ll.expression.i),
                           group=c(rep("sd",18),rep("ll",12)))
    
    out.i <- circacompare(x = sd.ll.df, col_time = "time", col_group = "group",
                          col_outcome = "measure",alpha_threshold = 1)
    circacompare.sd.ll[i,] <- out.i[[2]][,2]
  }
  else
  {
    circacompare.sd.ll[i,] <- NA
  }
  
}

colnames(circacompare.sd.ll) <- out.i[[2]][,1]
write.table(circacompare.sd.ll,file="circacompare_sd_ll.tsv", sep="\t")


# SD to DD
library(circacompare)

complete.sd.dd.rhythmic.genes <- read.table(file="tables/sd.dd.rhythmic.table.tsv", sep="\t")$x
circacompare.sd.dd <- matrix(nrow=length(complete.sd.dd.rhythmic.genes),ncol=15)
rownames(circacompare.sd.dd) <- complete.sd.dd.rhythmic.genes

for(i in 1:length(complete.sd.dd.rhythmic.genes))
{
  gene.i <- complete.sd.dd.rhythmic.genes[i]
  
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
  dd.expression.i <- genes.expression[gene.i,c(paste(dd.sd.zt,1,sep="_"),paste(dd.sd.zt,2,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 18)
  
  if ((length(unique(sd.dd.vector[1:12])) > 6) & (length(unique(sd.dd.vector[1:12])) > 3))
  {
    sd.dd.df <- data.frame(time=c(time.points,time.points[1:12]),
                           measure=c(sd.expression.i, dd.expression.i),
                           group=c(rep("sd",18),rep("dd",12)))
    
    out.i <- circacompare(x = sd.dd.df, col_time = "time", col_group = "group",
                          col_outcome = "measure",alpha_threshold = 1)
    circacompare.sd.dd[i,] <- out.i[[2]][,2]
  }
  else
  {
    circacompare.sd.dd[i,] <- NA
  }
  
}

colnames(circacompare.sd.dd) <- out.i[[2]][,1]
write.table(circacompare.sd.dd,file="circacompare_sd_dd.tsv", sep="\t")

# Read and analyze circacompare data
# LD-LL
library(data.table)
circacompare.ld.ll <- as.data.frame(fread("tables/estimates_ld_ll.tsv", sep="\t"))
circacompare.ld.ll.names <- circacompare.ld.ll$V1
circacompare.ld.ll <- circacompare.ld.ll[,-1]
rownames(circacompare.ld.ll) <- circacompare.ld.ll.names
head(circacompare.ld.ll)

png("boxplot_amp_ld_ll.png", width = 1500, height = 1500, res=300)
par(lwd=3)
boxplot(circacompare.ld.ll[,"ld amplitude estimate" ],circacompare.ld.ll[,"ll amplitude estimate" ],outline=F,col=c("blue","lightblue"),names=c("LD","LL"),cex.axis=2)
dev.off()
sum(circacompare.ld.ll[,"Amplitude difference estimate" ] < 0) # greater in LD
sum(circacompare.ld.ll[,"Amplitude difference estimate" ] < 0) / length(complete.ld.ll.rhythmic.genes)
sum(circacompare.ld.ll[,"P-value for amplitude difference" ] < 0.05)
sum(circacompare.ld.ll[,"P-value for amplitude difference" ] < 0.05)/length(complete.ld.ll.rhythmic.genes)
wilcox.test(x = circacompare.ld.ll[,"ld amplitude estimate" ], y = circacompare.ld.ll[,"ll amplitude estimate" ],alternative = "greater")$p.value

png("boxplot_mesor_ld_ll.png", width = 1500, height = 1500, res=300)
par(lwd=3)
boxplot(circacompare.ld.ll[,"ld mesor estimate" ],circacompare.ld.ll[,"ll mesor estimate" ],outline=F,col=c("blue","lightblue"),names=c("LD","LL"),cex.axis=2)
dev.off()
sum(circacompare.ld.ll[,"Mesor difference estimate" ] < 0) /length(complete.ld.ll.rhythmic.genes)
sum(circacompare.ld.ll[,"Mesor difference estimate" ] > 0) /length(complete.ld.ll.rhythmic.genes)
sum(circacompare.ld.ll[,"P-value for mesor difference" ] < 0.05)/length(complete.ld.ll.rhythmic.genes)
wilcox.test(x = circacompare.ld.ll[,"ld mesor estimate" ], y = circacompare.ld.ll[,"ll mesor estimate" ],alternative = "greater")$p.value

png("hist_ld_ll.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.ll[,"Phase difference estimate"],xlab="Phase Difference (h)",ylab="Number of Genes", col="lightblue",border="blue",main="",lwd=2,cex.axis=1.5,
     cex.lab=1.8)
abline(v=0, lty=2)
dev.off()
sum(circacompare.ld.ll[,"Phase difference estimate"] < 0)
sum(circacompare.ld.ll[,"Phase difference estimate"] > 0)/length(complete.ld.ll.rhythmic.genes)
sum(circacompare.ld.ll[,"P-value for difference in phase" ] < 0.05)
sum(circacompare.ld.ll[,"P-value for difference in phase" ] < 0.05)/length(complete.ld.ll.rhythmic.genes)

# LD-DD
library(data.table)
circacompare.ld.dd <- as.data.frame(fread("tables/estimates_ld_dd.tsv", sep="\t"))
circacompare.ld.dd.names <- circacompare.ld.dd$V1
circacompare.ld.dd <- circacompare.ld.dd[,-1]
rownames(circacompare.ld.dd) <- circacompare.ld.dd.names
head(circacompare.ld.dd)

png("boxplot_amp_ld_dd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
boxplot(circacompare.ld.dd[,"ld amplitude estimate" ],circacompare.ld.dd[,"dd amplitude estimate" ],outline=F,col=c("blue","darkblue"),names=c("LD","DD"),cex.axis=2)
dev.off()
sum(circacompare.ld.dd[,"Amplitude difference estimate" ] > 0) # lesser in DD
sum(circacompare.ld.dd[,"Amplitude difference estimate" ] > 0) / length(complete.ld.dd.rhythmic.genes)
sum(circacompare.ld.dd[,"P-value for amplitude difference" ] < 0.05)
sum(circacompare.ld.dd[,"P-value for amplitude difference" ] < 0.05)/length(complete.ld.dd.rhythmic.genes)
wilcox.test(x = circacompare.ld.dd[,"ld amplitude estimate" ], y = circacompare.ld.dd[,"dd amplitude estimate" ],alternative = "greater")$p.value

png("boxplot_mesor_ld_dd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
boxplot(circacompare.ld.dd[,"ld mesor estimate" ],circacompare.ld.dd[,"dd mesor estimate" ],outline=F,col=c("blue","darkblue"),names=c("LD","DD"),cex.axis=2)
dev.off()
sum(circacompare.ld.dd[,"Mesor difference estimate" ] < 0) /length(complete.ld.dd.rhythmic.genes)
sum(circacompare.ld.dd[,"Mesor difference estimate" ] > 0) /length(complete.ld.dd.rhythmic.genes)
sum(circacompare.ld.dd[,"P-value for mesor difference" ] < 0.05)/length(complete.ld.dd.rhythmic.genes)
wilcox.test(x = circacompare.ld.dd[,"ld mesor estimate" ], y = circacompare.ld.dd[,"dd mesor estimate" ],alternative = "l")$p.value

png("hist_ld_dd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(-circacompare.ld.dd[,"Phase difference estimate"],xlab="Phase Difference (h)",ylab="Number of Genes", col="darkblue",border="blue",main="",lwd=2,cex.axis=1.5,
     cex.lab=1.8)
abline(v=0, lty=2)
dev.off()
sum(circacompare.ld.dd[,"Phase difference estimate"] < 0)
sum(circacompare.ld.dd[,"Phase difference estimate"] > 0)/length(complete.ld.dd.rhythmic.genes) # Reduced in DD
sum(circacompare.ld.dd[,"P-value for difference in phase" ] < 0.05)
sum(circacompare.ld.dd[,"P-value for difference in phase" ] < 0.05)/length(complete.ld.dd.rhythmic.genes)

# Boxplots for amplitude and mesor of LD LL DD
png("boxplot_amp_ld_ll_dd_4.png", width = 1800, height = 1800, res=300)
par(lwd=3)
boxplot(circacompare.ld.ll[,"ld amplitude estimate" ],
        circacompare.ld.ll[,"ll amplitude estimate" ],
        circacompare.ld.dd[,"ld amplitude estimate" ],
        circacompare.ld.dd[,"dd amplitude estimate" ],
        outline=F,col=c("blue","lightblue", "blue", "darkblue"),names=c("LD-LL","LL","LD-DD","DD"),cex.axis=2)
dev.off()
wilcox.test(x = circacompare.ld.ll[,"ld amplitude estimate" ],y = circacompare.ld.ll[,"ll amplitude estimate" ], alternative = "greater")$p.value
wilcox.test(x = circacompare.ld.ll[,"ld amplitude estimate" ],y = circacompare.ld.dd[,"dd amplitude estimate" ], alternative = "greater")$p.value
wilcox.test(x = circacompare.ld.ll[,"ll amplitude estimate" ],y = circacompare.ld.dd[,"dd amplitude estimate" ], alternative = "g")$p.value

png("boxplot_mesor_ld_ll_dd_4.png", width = 1800, height = 1800, res=300)
par(lwd=3)
boxplot(circacompare.ld.ll[,"ld mesor estimate" ],
        circacompare.ld.ll[,"ll mesor estimate" ],
        circacompare.ld.dd[,"ld mesor estimate" ],
        circacompare.ld.dd[,"dd mesor estimate" ],
        outline=F,col=c("blue","lightblue", "blue","darkblue"),names=c("LD-LL","LL","LD-DD","DD"),cex.axis=2)
dev.off()

wilcox.test(x = circacompare.ld.ll[,"ld mesor estimate" ],y = circacompare.ld.ll[,"ll mesor estimate" ])$p.value
wilcox.test(x = circacompare.ld.ll[,"ld mesor estimate" ],y = circacompare.ld.dd[,"dd mesor estimate" ])$p.value
wilcox.test(x = circacompare.ld.ll[,"ll mesor estimate" ],y = circacompare.ld.dd[,"dd mesor estimate" ])$p.value


# SD-LL
library(data.table)
circacompare.sd.ll <- as.data.frame(fread("tables/estimates_sd_ll.tsv", sep="\t"))
circacompare.sd.ll.names <- circacompare.sd.ll$V1
circacompare.sd.ll <- circacompare.sd.ll[,-1]
rownames(circacompare.sd.ll) <- circacompare.sd.ll.names
head(circacompare.sd.ll)

png("boxplot_amp_sd_ll.png", width = 1500, height = 1500, res=300)
par(lwd=3)
boxplot(circacompare.sd.ll[,"sd amplitude estimate" ],circacompare.sd.ll[,"ll amplitude estimate" ],outline=F,col=c("red","lightsalmon"),names=c("SD","LL"),cex.axis=2)
dev.off()
sum(circacompare.sd.ll[,"Amplitude difference estimate" ] > 0) # lesser in LL
sum(circacompare.sd.ll[,"Amplitude difference estimate" ] > 0) / length(complete.sd.ll.rhythmic.genes)
sum(circacompare.sd.ll[,"P-value for amplitude difference" ] < 0.05)
sum(circacompare.sd.ll[,"P-value for amplitude difference" ] < 0.05)/length(complete.sd.ll.rhythmic.genes)
wilcox.test(x = circacompare.sd.ll[,"sd amplitude estimate" ],y = circacompare.sd.ll[,"ll amplitude estimate" ], alternative = "greater" )$p.value

png("boxplot_mesor_sd_ll.png", width = 1500, height = 1500, res=300)
par(lwd=3)
boxplot(circacompare.sd.ll[,"sd mesor estimate" ],circacompare.sd.ll[,"ll mesor estimate" ],outline=F,col=c("red","lightsalmon"),names=c("SD","LL"),cex.axis=2)
dev.off()
sum(circacompare.sd.ll[,"Mesor difference estimate" ] < 0) /length(complete.sd.ll.rhythmic.genes)
sum(circacompare.sd.ll[,"Mesor difference estimate" ] > 0) /length(complete.sd.ll.rhythmic.genes)
sum(circacompare.sd.ll[,"P-value for mesor difference" ] < 0.05)/length(complete.sd.ll.rhythmic.genes)
wilcox.test(x = circacompare.sd.ll[,"sd mesor estimate" ],y = circacompare.sd.ll[,"ll mesor estimate" ], alternative = "greater" )$p.value

png("hist_sd_ll.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(-circacompare.sd.ll[,"Phase difference estimate"],xlab="Phase Difference (h)",ylab="Number of Genes", col="lightsalmon",border="red",main="",lwd=2,cex.axis=1.5,
     cex.lab=1.8)
abline(v=0, lty=2)
dev.off()
sum(-circacompare.sd.ll[,"Phase difference estimate"] > 0)
sum(-circacompare.sd.ll[,"Phase difference estimate"] > 0)/length(complete.sd.ll.rhythmic.genes)
sum(circacompare.sd.ll[,"P-value for difference in phase" ] < 0.05)
sum(circacompare.sd.ll[,"P-value for difference in phase" ] < 0.05)/length(complete.sd.ll.rhythmic.genes)

# SD-DD
library(data.table)
circacompare.sd.dd <- as.data.frame(fread("tables/estimates_sd_dd.tsv", sep="\t"))
circacompare.sd.dd.names <- circacompare.sd.dd$V1
circacompare.sd.dd <- circacompare.sd.dd[,-1]
rownames(circacompare.sd.dd) <- circacompare.sd.dd.names
head(circacompare.sd.dd)

png("boxplot_amp_sd_dd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
boxplot(circacompare.sd.dd[,"sd amplitude estimate" ],circacompare.sd.dd[,"dd amplitude estimate" ],outline=F,col=c("red","firebrick4"),names=c("SD","DD"),cex.axis=2)
dev.off()
sum(circacompare.sd.dd[,"Amplitude difference estimate" ] > 0)
sum(circacompare.sd.dd[,"Amplitude difference estimate" ] > 0) / length(complete.sd.dd.rhythmic.genes)
sum(circacompare.sd.dd[,"P-value for amplitude difference" ] < 0.05)
sum(circacompare.sd.dd[,"P-value for amplitude difference" ] < 0.05)/length(complete.sd.dd.rhythmic.genes)
wilcox.test(x = circacompare.sd.dd[,"sd amplitude estimate" ],y = circacompare.sd.dd[,"dd amplitude estimate" ], alternative = "greater" )$p.value


png("boxplot_mesor_sd_dd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
boxplot(circacompare.sd.dd[,"sd mesor estimate" ],circacompare.sd.dd[,"dd mesor estimate" ],outline=F,col=c("red","firebrick4"),names=c("SD","DD"),cex.axis=2)
dev.off()
sum(circacompare.sd.dd[,"Mesor difference estimate" ] < 0) /length(complete.sd.dd.rhythmic.genes)
sum(circacompare.sd.dd[,"Mesor difference estimate" ] > 0) /length(complete.sd.dd.rhythmic.genes)
sum(circacompare.sd.dd[,"P-value for mesor difference" ] < 0.05)/length(complete.sd.dd.rhythmic.genes)
wilcox.test(x = circacompare.sd.dd[,"sd mesor estimate" ],y = circacompare.sd.dd[,"dd mesor estimate" ], alternative = "l" )$p.value

png("hist_sd_dd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(-circacompare.sd.dd[,"Phase difference estimate"],xlab="Phase Difference estimate (h)",ylab="Number of Genes", col="firebrick4",border="red",main="",lwd=2, cex.axis=1.5, cex.lab=1.8)
abline(v=0, lty=2)
dev.off()
sum(circacompare.sd.dd[,"Phase difference estimate"] > 0)
sum(circacompare.sd.dd[,"Phase difference estimate"] > 0)/length(complete.sd.dd.rhythmic.genes)
sum(circacompare.sd.dd[,"P-value for difference in phase" ] < 0.05)
sum(circacompare.sd.dd[,"P-value for difference in phase" ] < 0.05)/length(complete.sd.dd.rhythmic.genes)

# Boxplots for amplitude and mesor of SD LL DD
png("boxplot_amp_sd_ll_dd_4.png", width = 1800, height = 1800, res=300)
par(lwd=3)
boxplot(circacompare.sd.ll[,"sd amplitude estimate" ],
        circacompare.sd.ll[,"ll amplitude estimate" ],
        circacompare.sd.dd[,"sd amplitude estimate" ],
        circacompare.sd.dd[,"dd amplitude estimate" ],
        outline=F,col=c("red","lightsalmon","red","firebrick4"),names=c("SD-LL","LL","SD-DD","DD"),cex.axis=2)
dev.off()

wilcox.test(x = circacompare.sd.ll[,"sd amplitude estimate" ],y = circacompare.sd.ll[,"ll amplitude estimate" ], alternative = "greater")$p.value
wilcox.test(x = circacompare.sd.ll[,"sd amplitude estimate" ],y = circacompare.sd.dd[,"dd amplitude estimate" ], alternative = "g")$p.value
wilcox.test(x = circacompare.sd.ll[,"ll amplitude estimate" ],y = circacompare.sd.dd[,"dd amplitude estimate" ], alternative = "g")$p.value

png("boxplot_mesor_sd_ll_dd_4.png", width = 1800, height = 1800, res=300)
par(lwd=3)
boxplot(circacompare.sd.ll[,"sd mesor estimate" ],
        circacompare.sd.ll[,"ll mesor estimate" ],
        circacompare.sd.dd[,"sd mesor estimate" ],
        circacompare.sd.dd[,"dd mesor estimate" ],
        outline=F,col=c("red","lightsalmon","red","firebrick4"),names=c("SD-LL","LL", "SD-DD","DD"),cex.axis=2)
dev.off()


# Complete mean per condition
# Gene expression matrix for samples under light-dark cycles
ld.genes.exp.in.ld <- genes.expression[complete.ld.rhythmic.genes,ld.zt.i]
ld.genes.exp.in.ll <- genes.expression[complete.ld.rhythmic.genes,ll.ld.zt.i]
ld.genes.exp.in.dd <- genes.expression[complete.ld.rhythmic.genes,dd.ld.zt.i]

sd.genes.exp.in.sd <- genes.expression[complete.sd.rhythmic.genes,sd.zt.i]
sd.genes.exp.in.ll <- genes.expression[complete.sd.rhythmic.genes,ll.sd.zt.i]
sd.genes.exp.in.dd <- genes.expression[complete.sd.rhythmic.genes,dd.sd.zt.i]

mean(ld.genes.exp.in.ld)
mean(ld.genes.exp.in.ll)
mean(ld.genes.exp.in.dd)
wilcox.test(x = ld.genes.exp.in.ld, y = ld.genes.exp.in.ll, alternative = "l")$p.value

mean(sd.genes.exp.in.sd)
mean(sd.genes.exp.in.ll)
mean(sd.genes.exp.in.dd)
wilcox.test(x = sd.genes.exp.in.sd,y = sd.genes.exp.in.dd, alternative = "l")$p.value

### Effects of light on transitions
# Effects of light on LD transitions
complete.ld.rhythmic.genes <- read.table("tables/ld.rhythmic.table.tsv", header = T)$V1
complete.ld.ll.rhythmic.genes <- read.table(file="tables/ld.ll.rhythmic.table.tsv", sep="\t")$x
complete.ld.dd.rhythmic.genes <- read.table(file="tables/ld.dd.rhythmic.table.tsv", sep="\t")$x
length(complete.ld.rhythmic.genes)
length(complete.ld.ll.rhythmic.genes)
length(complete.ld.dd.rhythmic.genes)

not.expressed.ld <- names(which(apply(X = genes.expression,MARGIN = 1,FUN = max) < 2))
non.rhythmic.ld.genes <- setdiff(mp.genes,complete.ld.rhythmic.genes)
non.rhythmic.ld.genes.but.expressed <- setdiff(non.rhythmic.ld.genes, not.expressed.ld)
rhythmic.ld.ll.dd <- intersect(complete.ld.rhythmic.genes, intersect(complete.ld.ll.rhythmic.genes,
                                                                     complete.ld.dd.rhythmic.genes))
rhythmic.ld.ll.no.dd <- setdiff(complete.ld.ll.rhythmic.genes, rhythmic.ld.ll.dd)
rhythmic.ld.dd.no.ll <- setdiff(complete.ld.dd.rhythmic.genes, rhythmic.ld.ll.dd)
rhythmic.only.ld <- setdiff(setdiff(complete.ld.rhythmic.genes,complete.ld.ll.rhythmic.genes),
                            complete.ld.dd.rhythmic.genes)

ld.rhythmic.only.dd <- rhythmic.ld.dd.no.ll
ld.rhythmic.only.ll <- rhythmic.ld.ll.no.dd

# Dark regulation of LD
potential.light.dark.regulated <- setdiff(complete.ld.rhythmic.genes,unique(c(rhythmic.ld.ll.dd,ld.rhythmic.only.dd,ld.rhythmic.only.ll)))
length(potential.light.dark.regulated)

non.rhythmic.ld.dd.mean.expression.cycle <- apply(X = genes.expression[potential.light.dark.regulated,ld.zt.i],MARGIN = 1,FUN = mean)
non.rhythmic.ld.dd.mean.expression.dd <- apply(X = genes.expression[potential.light.dark.regulated,dd.ld.zt.i],MARGIN = 1,FUN = mean)

ld.potential.dark.repressed <- names(which(non.rhythmic.ld.dd.mean.expression.cycle > non.rhythmic.ld.dd.mean.expression.dd))
length(ld.potential.dark.repressed)
ld.potential.dark.repressed.expression.1 <- apply(X = genes.expression[ld.potential.dark.repressed,dd.ld.zt.i[1,]],MARGIN = 1,FUN = mean)
ld.potential.dark.repressed.expression.2 <- apply(X = genes.expression[ld.potential.dark.repressed,dd.ld.zt.i[2,]],MARGIN = 1,FUN = mean)
ld.dark.repressed <- names(which(ld.potential.dark.repressed.expression.1 >= ld.potential.dark.repressed.expression.2))
length(ld.dark.repressed)

ld.potential.dark.activated <- names(which(non.rhythmic.ld.dd.mean.expression.cycle < non.rhythmic.ld.dd.mean.expression.dd))
length(ld.potential.dark.activated)
ld.potential.dark.activated.expression.1 <- apply(X = genes.expression[ld.potential.dark.activated,dd.ld.zt.i[1,]],MARGIN = 1,FUN = mean)
ld.potential.dark.activated.expression.2 <- apply(X = genes.expression[ld.potential.dark.activated,dd.ld.zt.i[2,]],MARGIN = 1,FUN = mean)
ld.dark.activated <- names(which(ld.potential.dark.activated.expression.1 <= ld.potential.dark.activated.expression.2))
length(ld.dark.activated)

ld.no.clear.effect.dark <- setdiff(potential.light.dark.regulated, 
                                    c(ld.dark.activated,ld.dark.repressed))
length(ld.no.clear.effect.dark)

# Light regulation of LD
non.rhythmic.ld.ll.mean.expression.cycle <- apply(X = genes.expression[potential.light.dark.regulated,ld.zt.i],MARGIN = 1,FUN = mean)

non.rhythmic.ld.ll.mean.expression.ll <- apply(X = genes.expression[potential.light.dark.regulated,ll.ld.zt.i],MARGIN = 1,FUN = mean)

ld.potential.light.repressed <- names(which(non.rhythmic.ld.ll.mean.expression.cycle > non.rhythmic.ld.ll.mean.expression.ll))
length(ld.potential.light.repressed)
ld.potential.light.repressed.expression.1 <- apply(X = genes.expression[ld.potential.light.repressed,ll.ld.zt.i[1,]],MARGIN = 1,FUN = mean)
ld.potential.light.repressed.expression.2 <- apply(X = genes.expression[ld.potential.light.repressed,ll.ld.zt.i[2,]],MARGIN = 1,FUN = mean)
ld.light.repressed <- names(which(ld.potential.light.repressed.expression.1 >= ld.potential.light.repressed.expression.2))
length(ld.light.repressed)

ld.potential.light.activated <- names(which(non.rhythmic.ld.ll.mean.expression.cycle < non.rhythmic.ld.ll.mean.expression.ll))
length(ld.potential.light.activated)
ld.potential.light.activated.expression.1 <- apply(X = genes.expression[ld.potential.light.activated,ll.ld.zt.i[1,]],MARGIN = 1,FUN = mean)
ld.potential.light.activated.expression.2 <- apply(X = genes.expression[ld.potential.light.activated,ll.ld.zt.i[2,]],MARGIN = 1,FUN = mean)
ld.light.activated <- names(which(ld.potential.light.activated.expression.1 <= ld.potential.light.activated.expression.2))
length(ld.light.activated)

ld.no.clear.effect.light <- setdiff(potential.light.dark.regulated, 
                                    c(ld.light.activated,ld.light.repressed))
length(ld.no.clear.effect.light)

# Genes that are activated by light and repressed by dark
length(ld.light.activated)
length(ld.dark.repressed)

ld.light.activated.dark.repressed <- intersect(ld.light.activated,ld.dark.repressed)
length(ld.light.activated.dark.repressed)

# Genes that are repressed by light and activated by dark
length(ld.light.repressed)
length(ld.dark.activated)

ld.light.repressed.dark.activated <- intersect(ld.light.repressed,ld.dark.activated)
length(ld.light.repressed.dark.activated)

# Genes that are regulated in opposites directions by light and dark, total
ld.light.dark.regulated <- c(ld.light.activated.dark.repressed,ld.light.repressed.dark.activated)
length(ld.light.dark.regulated)

length(potential.light.dark.regulated)
ld.no.clear.effect.light.dark.regulation <- setdiff(potential.light.dark.regulated,ld.light.dark.regulated)
length(ld.no.clear.effect.light.dark.regulation)

write.table(x = ld.light.activated.dark.repressed, file="tables/LD_light_activated_dark_repressed.tsv",row.names = F)

write.table(x = ld.light.repressed.dark.activated, file="tables/LD_light_repressed_dark_activated.tsv",row.names = F)

write.table(x = ld.no.clear.effect.light.dark.regulation, file="tables/LD_no_clear_effect_light_dark_regulation.tsv",row.names = F)

# Effects of light on SD transitions
complete.sd.rhythmic.genes <- read.table("tables/sd.rhythmic.table.tsv", header = T)$V1
complete.sd.ll.rhythmic.genes <- read.table(file="tables/sd.ll.rhythmic.table.tsv", sep="\t")$x
complete.sd.dd.rhythmic.genes <- read.table(file="tables/sd.dd.rhythmic.table.tsv", sep="\t")$x
length(complete.sd.rhythmic.genes)
length(complete.sd.ll.rhythmic.genes)
length(complete.sd.dd.rhythmic.genes)

not.expressed.sd <- names(which(apply(X = genes.expression,MARGIN = 1,FUN = max) < 2))
non.rhythmic.sd.genes <- setdiff(mp.genes,complete.sd.rhythmic.genes)
non.rhythmic.sd.genes.but.expressed <- setdiff(non.rhythmic.sd.genes, not.expressed.sd)
rhythmic.sd.ll.dd <- intersect(complete.sd.rhythmic.genes, intersect(complete.sd.ll.rhythmic.genes,
                                                                     complete.sd.dd.rhythmic.genes))
rhythmic.sd.ll.no.dd <- setdiff(complete.sd.ll.rhythmic.genes, rhythmic.sd.ll.dd)
rhythmic.sd.dd.no.ll <- setdiff(complete.sd.dd.rhythmic.genes, rhythmic.sd.ll.dd)
rhythmic.only.sd <- setdiff(setdiff(complete.sd.rhythmic.genes,complete.sd.ll.rhythmic.genes),
                            complete.sd.dd.rhythmic.genes)

sd.rhythmic.only.dd <- rhythmic.sd.dd.no.ll
sd.rhythmic.only.ll <- rhythmic.sd.ll.no.dd

# Dark regulation of SD
length(complete.sd.rhythmic.genes)
potential.light.dark.regulated <- setdiff(complete.sd.rhythmic.genes,unique(c(rhythmic.sd.ll.dd,sd.rhythmic.only.dd,sd.rhythmic.only.ll)))
length(potential.light.dark.regulated)

non.rhythmic.sd.dd.mean.expression.cycle <- apply(X = genes.expression[potential.light.dark.regulated,sd.zt.i],MARGIN = 1,FUN = mean)
non.rhythmic.sd.dd.mean.expression.dd <- apply(X = genes.expression[potential.light.dark.regulated,dd.sd.zt.i],MARGIN = 1,FUN = mean)

sd.potential.dark.repressed <- names(which(non.rhythmic.sd.dd.mean.expression.cycle > non.rhythmic.sd.dd.mean.expression.dd))
length(sd.potential.dark.repressed)
sd.potential.dark.repressed.expression.1 <- apply(X = genes.expression[sd.potential.dark.repressed,dd.sd.zt.i[1,]],MARGIN = 1,FUN = mean)
sd.potential.dark.repressed.expression.2 <- apply(X = genes.expression[sd.potential.dark.repressed,dd.sd.zt.i[2,]],MARGIN = 1,FUN = mean)
sd.dark.repressed <- names(which(sd.potential.dark.repressed.expression.1 >= sd.potential.dark.repressed.expression.2))
length(sd.dark.repressed)

sd.potential.dark.activated <- names(which(non.rhythmic.sd.dd.mean.expression.cycle < non.rhythmic.sd.dd.mean.expression.dd))
length(sd.potential.dark.activated)
sd.potential.dark.activated.expression.1 <- apply(X = genes.expression[sd.potential.dark.activated,dd.sd.zt.i[1,]],MARGIN = 1,FUN = mean)
sd.potential.dark.activated.expression.2 <- apply(X = genes.expression[sd.potential.dark.activated,dd.sd.zt.i[2,]],MARGIN = 1,FUN = mean)
sd.dark.activated <- names(which(sd.potential.dark.activated.expression.1 <= sd.potential.dark.activated.expression.2))
length(sd.dark.activated)

sd.no.clear.effect.dark <- setdiff(potential.light.dark.regulated, 
                                   c(ld.dark.activated,ld.dark.repressed))
length(sd.no.clear.effect.dark)

# Light regulation of SD
non.rhythmic.sd.ll.mean.expression.cycle <- apply(X = genes.expression[potential.light.dark.regulated,sd.zt.i],MARGIN = 1,FUN = mean)

non.rhythmic.sd.ll.mean.expression.ll <- apply(X = genes.expression[potential.light.dark.regulated,ll.sd.zt.i],MARGIN = 1,FUN = mean)

sd.potential.light.repressed <- names(which(non.rhythmic.sd.ll.mean.expression.cycle > non.rhythmic.sd.ll.mean.expression.ll))
length(sd.potential.light.repressed)
sd.potential.light.repressed.expression.1 <- apply(X = genes.expression[sd.potential.light.repressed,ll.sd.zt.i[1,]],MARGIN = 1,FUN = mean)
sd.potential.light.repressed.expression.2 <- apply(X = genes.expression[sd.potential.light.repressed,ll.sd.zt.i[2,]],MARGIN = 1,FUN = mean)
sd.light.repressed <- names(which(sd.potential.light.repressed.expression.1 >= sd.potential.light.repressed.expression.2))
length(sd.light.repressed)

sd.potential.light.activated <- names(which(non.rhythmic.sd.ll.mean.expression.cycle < non.rhythmic.sd.ll.mean.expression.ll))
length(sd.potential.light.activated)
sd.potential.light.activated.expression.1 <- apply(X = genes.expression[sd.potential.light.activated,ll.sd.zt.i[1,]],MARGIN = 1,FUN = mean)
sd.potential.light.activated.expression.2 <- apply(X = genes.expression[sd.potential.light.activated,ll.sd.zt.i[2,]],MARGIN = 1,FUN = mean)
sd.light.activated <- names(which(sd.potential.light.activated.expression.1 <= sd.potential.light.activated.expression.2))
length(sd.light.activated)

sd.no.clear.effect.light <- setdiff(potential.light.dark.regulated, 
                                    c(sd.light.activated,sd.light.repressed))
length(sd.no.clear.effect.light)

length(sd.light.activated)
length(sd.dark.repressed)

sd.light.activated.dark.repressed <- intersect(sd.light.activated,sd.dark.repressed)
length(sd.light.activated.dark.repressed)

length(intersect(sd.light.activated.dark.repressed,ld.light.activated.dark.repressed))

length(sd.light.repressed)
length(sd.dark.activated)

sd.light.repressed.dark.activated <- intersect(sd.light.repressed,sd.dark.activated)
length(sd.light.repressed.dark.activated)

length(intersect(sd.light.repressed.dark.activated,ld.light.repressed.dark.activated))

sd.light.dark.regulated <- c(sd.light.activated.dark.repressed,sd.light.repressed.dark.activated)
length(sd.light.dark.regulated)

length(potential.light.dark.regulated)
sd.no.clear.effect.light.dark.regulation <- setdiff(potential.light.dark.regulated,sd.light.dark.regulated)
length(sd.no.clear.effect.light.dark.regulation)

write.table(x = sd.light.activated.dark.repressed, file="tables/SD_light_activated_dark_repressed.tsv",row.names = F)

write.table(x = sd.light.repressed.dark.activated, file="tables/SD_light_repressed_dark_activated.tsv",row.names = F)

write.table(x = sd.no.clear.effect.light.dark.regulation, file="tables/SD_no_clear_effect_light_dark_regulation.tsv",row.names = F)

# Circacompare for LD-SD differences
# Rhythmic genes in LD and SD 
complete.sd.rhythmic.genes <- read.table("tables/sd.rhythmic.table.tsv", header = T)$V1
complete.ld.rhythmic.genes <- read.table("tables/ld.rhythmic.table.tsv", header = T)$V1
both.ld.sd.rhythmic.genes <- intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(both.ld.sd.rhythmic.genes)

circacompare.ld.sd <- matrix(nrow=length(both.ld.sd.rhythmic.genes),ncol=15)
rownames(circacompare.ld.sd) <- both.ld.sd.rhythmic.genes

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

for(i in 1:length(both.ld.sd.rhythmic.genes))
{
  gene.i <- both.ld.sd.rhythmic.genes[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"), paste(ld.zt,3,sep="_"))]
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"), paste(sd.zt,3,sep="_"))]
  time.points <- seq(from=0,by=4,length.out = 18)
  
  if ((length(unique(ld.expression.i[1:12])) > 2) & (length(unique(sd.expression.i[1:12])) > 2))
  {
    ld.sd.df <- data.frame(time=c(time.points,time.points),
                           measure=c(ld.expression.i, sd.expression.i),
                           group=c(rep("ld",18),rep("sd",18)))
    
    out.i <- circacompare(x = ld.sd.df, col_time = "time", col_group = "group",
                          col_outcome = "measure",alpha_threshold = 1 )
    circacompare.ld.sd[i,] <- out.i[[2]][,2]
  }
  else
  {
    circacompare.ld.sd[i,] <- NA
  }
}

colnames(circacompare.ld.sd) <- out.i[[2]][,1]
write.table(circacompare.ld.sd,file="tables/circacompare_intersection.tsv", sep="\t")

# Read circacompare data
library(data.table)
circacompare.ld.sd <- as.data.frame(fread("tables/circacompare_intersection.tsv", sep="\t"))
circacompare.ld.sd.names <- circacompare.ld.sd$V1
circacompare.ld.sd <- circacompare.ld.sd[,-1]
rownames(circacompare.ld.sd) <- circacompare.ld.sd.names

#Amplitude
png("boxplot_amp_ld_sd_intersection.png", width = 1500, height = 1500, res=300)
par(lwd=3)
boxplot(circacompare.ld.sd[,"ld amplitude estimate" ],
        circacompare.ld.sd[,"sd amplitude estimate" ],outline=F,
        col=c("royalblue3","tomato3"),names=c("LD","SD"),cex.axis=1.5, main="Distribution of amplitudes")
dev.off()
sum(circacompare.ld.sd[,"Amplitude difference estimate" ] < 0)
sum(circacompare.ld.sd[,"Amplitude difference estimate" ] < 0) / length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"P-value for amplitude difference" ] < 0.05)
sum(circacompare.ld.sd[,"P-value for amplitude difference" ] < 0.05)/length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[circacompare.ld.sd[,"Amplitude difference estimate" ] > 0,"P-value for amplitude difference" ] < 0.05)
sum(circacompare.ld.sd[circacompare.ld.sd[,"Amplitude difference estimate" ] < 0,"P-value for amplitude difference" ] < 0.05)

wilcox.test(x = circacompare.ld.sd[,"ld amplitude estimate" ], y = circacompare.ld.sd[,"sd amplitude estimate" ], alternative = "g")$p.value

png("boxplot_mesor_ld_sd_intersection.png", width = 1500, height = 1500, res=300)
par(lwd=3)
boxplot(circacompare.ld.sd[,"ld mesor estimate" ],circacompare.ld.sd[,"sd mesor estimate" ],
        outline=F,col=c("royalblue3","tomato3"),names=c("LD","SD"),cex.axis=1.5, main="Distribution of mesor")
dev.off()
sum(circacompare.ld.sd[,"Mesor difference estimate" ] < 0) /length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"Mesor difference estimate" ] > 0) /length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"P-value for mesor difference" ] < 0.05)/length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[circacompare.ld.sd[,"Mesor difference estimate" ] > 0,"P-value for mesor difference" ] < 0.05)
sum(circacompare.ld.sd[circacompare.ld.sd[,"Mesor difference estimate" ] < 0,"P-value for mesor difference" ] < 0.05)

wilcox.test(x = circacompare.ld.sd[,"ld mesor estimate" ], y = circacompare.ld.sd[,"sd mesor estimate" ], alternative = "l")$p.value


#Phase
png("boxplot_phase_ld_sd_intersection.png", width = 1500, height = 1500, res=300)
par(lwd=3)
boxplot(circacompare.ld.sd[,"ld peak time hours" ],circacompare.ld.sd[,"sd peak time hours" ],
        outline=F,col=c("royalblue3","tomato3"),names=c("LD","SD"),cex.axis=1.5, main="Distribution of phase")
dev.off()

sum(circacompare.ld.sd[,"Phase difference estimate" ] < 0) /length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"Phase difference estimate" ] > 0) /length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"P-value for difference in phase" ] < 0.05)/length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"P-value for difference in phase" ] < 0.05)
mean(circacompare.ld.sd[,"Phase difference estimate" ])
mean(circacompare.ld.sd[,"ld peak time hours" ])-mean(circacompare.ld.sd[,"sd peak time hours" ])
median(circacompare.ld.sd[,"ld peak time hours" ])-median(circacompare.ld.sd[,"sd peak time hours" ])

# This should be circular
source("../../../PSEA_R/GeneSCARAB/PSEA_functions_package.R")
library(circular)
custom_set <- list(circacompare.ld.sd[,"ld peak time hours" ],
                   circacompare.ld.sd[,"sd peak time hours" ]) # Load vector of phases
ind_table <- sapply(custom_set, function(x) create_circular_table(circa_vector = x, 
                                                                  hr.on.large.sets = F , hr.on.large.sets.th = 100, 
                                                                  force.hr = F, force.hr.th = 0.05, 
                                                                  rao.on.large.sets = F, rao.on.large.sets.th = 100, 
                                                                  force.rao = F, force.rao.th = 0.05,
                                                                  iter.rao = 100, iter.hr = 100))

ind_table2 <- as.data.frame(t(ind_table))
ind_table2 <- matrix(c(as.numeric(ind_table[,2]), as.numeric(ind_table[,1])), ncol = 11, byrow =T)
colnames(ind_table2) <- rownames(ind_table)
rownames(ind_table2) <- c("SD", "LD")
ind_table2 <- as.data.frame(ind_table2)

# Transform negative means and quantiles to positive
pos_res <- apply(ind_table2[,2:5], MARGIN=2, function(x) ifelse(as.numeric(x) < 0, 24 + as.numeric(x) , as.numeric(x)))
ind_table2[,2:5] <- pos_res
ind_table2

library(circlize)
library(MetBrewer)
library(ComplexHeatmap)

circular_boxplot(res.table = ind_table2, "ld_sd_phases.png", tr.height = 0.14, color.palette = "Egypt")

# Comparison of lineal vs circular medians
median(circacompare.ld.sd[,"ld peak time hours" ])-median(circacompare.ld.sd[,"sd peak time hours" ]) # 3.18h linear
ind_table2$median[2] - ind_table2$median[1] # 5.75h circular
ind_table2$mean[2] - ind_table2$mean[1] # 5.32h circular
median(circacompare.ld.sd[,"Phase difference estimate" ]) # 3.54 median of gene-by-gene difference

# MANOVA to compare circular distributions
summary(manova(cbind(sin(c(custom_set[[1]]*pi/12, custom_set[[2]]*pi/12)),
                     cos(c(custom_set[[1]]*pi/12, custom_set[[2]]*pi/12))) 
               ~ c(rep("ld", length(custom_set[[1]])), rep("sd", length(custom_set[[2]]))))) # Significant

# Phase histograms
png("hist_phase_ld_intersection.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[,"ld peak time hours" ],col="blue",
     ylim=c(0,800),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

png("hist_phase_sd_intersection.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[,"sd peak time hours" ],col="red",
     ylim=c(0,1000),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=8, col="black", lwd=5, lty = 2)
dev.off()

###### Functional analysis
# Using GeneSCARAB
# First, load annotation and phases information for LD and SD
circa_table_genescarab <- read.table("tables/circacompare_intersection.tsv", header = T, sep = "\t")
total_phases_table_ld <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["ld.peak.time.hours"]]))
total_phases_table_sd <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["sd.peak.time.hours"]]))


### GO Analysis
# Select all GOs that appear in the annotation package
library(GO.db)
library(org.Mpolymorpha71.eg.db)
functional_data <- select(org.Mpolymorpha71.eg.db, 
                          keys = keys(org.Mpolymorpha71.eg.db,keytype="GID"), 
                          columns =  c("GID", "GO")) # load annotation data
complete_gos <- unique(functional_data$GO)
complete_gos <- complete_gos[!is.na(complete_gos)]

# Generate GO terms list with genes per GO term (including ancestors of selected GO terms)
go.list.test <- create_gene_list_go(go_vector = complete_gos, org.package = "org.Mpolymorpha71.eg.db",
                                    go_column = "GO", id_column = "GID")

# Generate GO terms list with phases per GO term for LD and SD
phases.list.ld <- gene.list.to.phases(go.list.test, total_phases_table_ld)
phases.list.sd <- gene.list.to.phases(go.list.test, total_phases_table_sd)

# Remove GOs with no associated rhytmic genes
phases.list.ld.clean <- phases.list.ld[which(sapply(phases.list.ld, nrow) != 0)]
phases.list.sd.clean <- phases.list.sd[which(sapply(phases.list.sd, nrow) != 0)]

# Test preferential concentration of GO terms during the diel cycle
library(circular)

go.circa.table.ld <- complete_circular_table(phase.list = phases.list.ld.clean, 
                      hr.on.large.sets = F , hr.on.large.sets.th = 200, 
                      force.hr = F, force.hr.th = 0.05, 
                      rao.on.large.sets = F, rao.on.large.sets.th = 200, 
                      force.rao = F, force.rao.th = 0.05,
                      iter.rao = 999, iter.hr = 999)

write.table(go.circa.table.ld,file="tables/go_circa_table_ld.tsv", sep="\t")

go.circa.table.sd <- complete_circular_table(phase.list = phases.list.sd.clean, 
                                             hr.on.large.sets = F , hr.on.large.sets.th = 200, 
                                             force.hr = F, force.hr.th = 0.05, 
                                             rao.on.large.sets = F, rao.on.large.sets.th = 200, 
                                             force.rao = F, force.rao.th = 0.05,
                                             iter.rao = 999, iter.hr = 999)
write.table(go.circa.table.sd,file="tables/go_circa_table_sd.tsv", sep="\t")

# Read and explore results
go.circa.table.ld <- read.table("tables/go_circa_table_ld.tsv", sep="\t", header = T)
go.circa.table.sd <- read.table("tables/go_circa_table_sd.tsv", sep="\t", header = T)

sum(go.circa.table.ld$rayleigh_p_value_adj < 0.05 | go.circa.table.ld$kuiper_p_value_adj < 0.05) # 1261
sum(go.circa.table.sd$rayleigh_p_value_adj < 0.05 | go.circa.table.sd$kuiper_p_value_adj < 0.05) # 1191
sum(go.circa.table.ld$hr_p_value_adj < 0.05, na.rm = T) # 1920
sum(go.circa.table.sd$hr_p_value_adj < 0.05, na.rm = T) # 1923
sum(go.circa.table.ld$rao_p_value_adj < 0.05, na.rm = T) # 0
sum(go.circa.table.sd$rao_p_value_adj < 0.05, na.rm = T) # 1

rownames(go.circa.table.sd)[which(go.circa.table.sd$rao_p_value_adj < 0.05)]
plot(circular(phases.list.sd$`GO:0051253`$phase)) # avoid night


unimod_ld <- rownames(go.circa.table.ld)[which(go.circa.table.ld$rayleigh_p_value_adj < 0.05 | go.circa.table.ld$kuiper_p_value_adj < 0.05)]
unimod_sd <- rownames(go.circa.table.sd)[which(go.circa.table.sd$rayleigh_p_value_adj < 0.05 | go.circa.table.sd$kuiper_p_value_adj < 0.05)]

multimod_ld <- rownames(go.circa.table.ld)[which(go.circa.table.ld$hr_p_value_adj < 0.05 | go.circa.table.ld$rao_p_value_adj < 0.05)]
multimod_sd <- rownames(go.circa.table.sd)[which(go.circa.table.sd$hr_p_value_adj < 0.05 | go.circa.table.sd$rao_p_value_adj < 0.05)]

library(VennDiagram)
png("ld_sd_unimod_go.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(unimod_ld), area2 = length(unimod_sd),
                   cross.area = length(intersect(unimod_ld,unimod_sd)),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

png("ld_sd_multimod_go.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(multimod_ld), area2 = length(multimod_sd),
                   cross.area = length(intersect(multimod_ld,multimod_sd)),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# For the GO terms not following an uniform distribution, test them against the general distribution of phases
# per photoperiod
non_uniform_gos_ld <- unique(c(unimod_ld, multimod_ld))
phases.list.ld.non.uni <- phases.list.ld.clean[non_uniform_gos_ld]
gen.dist.ld <- test_against_gen_dist(phase.list = phases.list.ld.non.uni, total.phase.table = total_phases_table_ld)
write.table(gen.dist.ld,file="tables/go_general_dist_ld.tsv", sep="\t")

non_uniform_gos_sd <- unique(c(unimod_sd, multimod_sd))
phases.list.sd.non.uni <- phases.list.sd.clean[non_uniform_gos_sd]
gen.dist.sd <- test_against_gen_dist(phase.list = phases.list.sd.non.uni, total.phase.table = total_phases_table_sd)
write.table(gen.dist.sd,file="tables/go_general_dist_sd.tsv", sep="\t")

gen.dist.ld <- read.table("tables/go_general_dist_ld.tsv", sep="\t", header = T)
gen.dist.sd <- read.table("tables/go_general_dist_sd.tsv", sep="\t", header = T)

sum(gen.dist.ld$gen_dist_p_value_adj < 0.05) # 357 non-uniform GO terms present a different distribution than that of global
# mRNA accumulation in LD
sum(gen.dist.sd$gen_dist_p_value_adj < 0.05) # 197 non-uniform GO terms present a different distribution than that of global
# mRNA accumulation in SD

go_dist_gen_ld <- rownames(gen.dist.ld)[which(gen.dist.ld$gen_dist_p_value_adj < 0.05)]
go_dist_gen_sd <- rownames(gen.dist.sd)[which(gen.dist.sd$gen_dist_p_value_adj < 0.05)]

png("ld_sd_general_dist.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(go_dist_gen_ld), area2 = length(go_dist_gen_sd),
                   cross.area = length(intersect(go_dist_gen_ld,go_dist_gen_sd)),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# Test shared GO terms not following an uniform distribution for differences by photoperiod
shared_non_uniform <- intersect(non_uniform_gos_ld, non_uniform_gos_sd)
phases.list.shared.ld.non.uni <- phases.list.ld.non.uni[shared_non_uniform]
phases.list.shared.sd.non.uni <- phases.list.sd.non.uni[shared_non_uniform]

go_diff_ld_sd <- mapply(function(x, y) summary(manova(cbind(sin(c(x$phase*pi/12, y$phase*pi/12)),cos(c(x$phase*pi/12, y$phase*pi/12))) 
                                     ~ c(rep("ld", nrow(x)), rep("sd", nrow(y)))))$stats[1,6],
                        phases.list.shared.ld.non.uni, phases.list.shared.sd.non.uni)

# Remove NULL results (due to a single phase in the set of one of them)
go_diff_ld_sd_clean <- go_diff_ld_sd[names(which(sapply(go_diff_ld_sd, length) != 0))]
go_diff_ld_sd_clean_bh <- p.adjust(unlist(go_diff_ld_sd_clean), method = "BH")

go_diff_ld_sd_table <- data.frame(gen_dist_p_value=unlist(go_diff_ld_sd_clean), 
                             gen_dist_p_value_adj=go_diff_ld_sd_clean_bh)

write.table(go_diff_ld_sd_table,file="tables/go_diff_ld_sd_table", sep="\t")

go_diff_ld_sd_table <- read.table("tables/go_diff_ld_sd_table", sep="\t", header = T)

sum(go_diff_ld_sd_table$gen_dist_p_value_adj < 0.05) # 573 change their distribution between LD and SD
change_gos <- rownames(go_diff_ld_sd_table)[which(go_diff_ld_sd_table$gen_dist_p_value_adj < 0.05)]

# Difference in medians for changing GOs
diff_gos_medians <- go.circa.table.ld[change_gos,"median"] - go.circa.table.sd[change_gos,"median"]

sd_anticipated <- change_gos[which((0 < diff_gos_medians & diff_gos_medians < 12) | (diff_gos_medians < -12))]
# 572 are anticipated in SD, only 1 in LD
length(sd_anticipated)
ld_anticipated <- setdiff(change_gos, sd_anticipated)

library(circlize)
library(ComplexHeatmap)
library(MetBrewer)

circ_ex_ld <- go.circa.table.ld[c("GO:0005978", "GO:0051015", "GO:0015693", "GO:0015095", "GO:0000079"),]
circular_boxplot(circ_ex_ld, "prueba_gos_ld.png", color.palette = "Tam",
                             tr.height = 0.1)
circ_ex_sd <- go.circa.table.sd[c("GO:0005978", "GO:0051015", "GO:0015693", "GO:0015095", "GO:0000079"),]
circular_boxplot(circ_ex_sd, "prueba_gos_sd.png", color.palette = "Tam",
                 tr.height = 0.1)


circ_ex_ld_big <- go.circa.table.ld[c("GO:0006139", "GO:0006725", "GO:0006807", "GO:0008150",
                                      "GO:0008152", "GO:0009987", "GO:0016070", "GO:0016071"),]
circular_boxplot(circ_ex_ld_big, "prueba_gos_ld_big.png", color.palette = "Tam",
                 tr.height = 0.1)
circ_ex_sd_big <- go.circa.table.sd[c("GO:0006139", "GO:0006725", "GO:0006807", "GO:0008150",
                                      "GO:0008152", "GO:0009987", "GO:0016070", "GO:0016071", "GO:1901565", "GO:0008233",
                                      "GO:0008152", "GO:0009987", "GO:0016070", "GO:0016071"),]
circular_boxplot(circ_ex_sd_big, "prueba_gos_sd_big3.png", color.palette = "Tam",
                 tr.height = 0.05)

# Description of changing GO terms
GO_complete_list <- as.list(GOTERM)
change_gos_description <- sapply(change_gos, function(x) GO_complete_list[[x]]@Term)
  
# Determine which ones correspond to BP, MF and CC
change_gos_ontology <- sapply(change_gos, function(x) GO_complete_list[[x]]@Ontology)
change_gos_description_bp <- change_gos_description[change_gos_ontology == "BP"]
change_gos_description_mf <- change_gos_description[change_gos_ontology == "MF"]
change_gos_description_cc <- change_gos_description[change_gos_ontology == "CC"]

# Add phase medians in each photoperiod

change_gos_bp_table <- data.frame(go=names(change_gos_description_bp), description=change_gos_description_bp,
                                  ld_phase=go.circa.table.ld[names(change_gos_description_bp), "median"],
                                  sd_phase=go.circa.table.sd[names(change_gos_description_bp), "median"])

change_gos_mf_table <- data.frame(go=names(change_gos_description_mf), description=change_gos_description_mf,
                                  ld_phase=go.circa.table.ld[names(change_gos_description_mf), "median"],
                                  sd_phase=go.circa.table.sd[names(change_gos_description_mf), "median"])

change_gos_cc_table <- data.frame(go=names(change_gos_description_cc), description=change_gos_description_cc,
                                  ld_phase=go.circa.table.ld[names(change_gos_description_cc), "median"],
                                  sd_phase=go.circa.table.sd[names(change_gos_description_cc), "median"])

write.table(change_gos_bp_table,file="tables/anticipated_gos_sd_bp.tsv", sep="\t", col.names = T)
write.table(change_gos_mf_table,file="tables/anticipated_gos_sd_mf.tsv", sep="\t", col.names = T)
write.table(change_gos_cc_table,file="tables/anticipated_gos_sd_cc.tsv", sep="\t", col.names = T)


### kegg Analysis
# Select all keggs that appear in the annotation package
library(KEGGREST)
library(org.Mpolymorpha71.eg.db)
library(stringr)

circa_table_genescarab <- read.table("tables/circacompare_intersection.tsv", header = T, sep = "\t")
total_phases_table_ld <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["ld.peak.time.hours"]]))
total_phases_table_sd <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["sd.peak.time.hours"]]))

functional_data <- select(org.Mpolymorpha71.eg.db, 
                          keys = keys(org.Mpolymorpha71.eg.db,keytype="GID"), 
                          columns =  c("GID", "KEGG")) # load annotation data
complete_kos <- unique(functional_data$KEGG)

complete_kos <- complete_kos[!is.na(complete_kos)]

# Generate kegg pathways list with genes per KO term
ko_vector = complete_kos
org.package = "org.Mpolymorpha71.eg.db"
ko_column = "KEGG"
id_column = "GID"

{
  functional_data <- select(get(org.package), 
                            keys = keys(get(org.package),keytype=id_column), 
                            columns =  c(id_column, ko_column))
  
  # Get pathways for selected KOs 
  
  present_pathways1 <- keggLink("pathway", paste0(ko_vector[1:500], collapse = "+"))
  present_pathways2 <- keggLink("pathway", paste0(ko_vector[501:1000], collapse = "+"))
  present_pathways3 <- keggLink("pathway", paste0(ko_vector[1001:1500], collapse = "+"))
  present_pathways4 <- keggLink("pathway", paste0(ko_vector[1501:2000], collapse = "+"))
  present_pathways5 <- keggLink("pathway", paste0(ko_vector[2001:2500], collapse = "+"))
  present_pathways6 <- keggLink("pathway", paste0(ko_vector[2501:3000], collapse = "+"))
  present_pathways7 <- keggLink("pathway", paste0(ko_vector[3001:3500], collapse = "+"))
  present_pathways8 <- keggLink("pathway", paste0(ko_vector[3501:4000], collapse = "+"))
  present_pathways9 <- keggLink("pathway", paste0(ko_vector[4001:4067], collapse = "+"))
  
  unique_pathways <- unique(unlist(mget(paste0("present_pathways", 1:9))))
  map_pathways <- unique_pathways[grep("map", unique_pathways)] # maybe change this for ko in other versions
  list_pathways <- sapply(strsplit(map_pathways, ":"), function(x) x[[2]])
  
  kos_pathways <- lapply(list_pathways, function(x) keggLink("ko",x))
  #kos_pathways <- keggLink("ko", paste0(list_pathways, collapse = "+"))
  names(kos_pathways) <- list_pathways
  kos_pathways_clean <- lapply(kos_pathways, function(x) unique(x))
  kos_pathways_clean <- lapply(kos_pathways_clean, function(x) sapply(strsplit(x, ":"), function(y) y[[2]]))
  #ko_list <- lapply(kos_pathways_clean, function(x) subset(functional_data, grepl(paste0(x, collapse = "|"), functional_data[["KO"]]))[["GID"]])
  ko_list <- lapply(kos_pathways_clean, function(x) subset(functional_data, str_detect(functional_data[["KEGG"]], paste0(x, collapse = "|")))[["GID"]])
  ko_list <- lapply(ko_list, unique)
}

# Filter out pathways not present in plants
kegg_plants <- read.csv("../../../PSEA_R/GeneSCARAB/pathways_plants.txt", sep = "\t", header = T)$x
kegg_plants <- paste0("map", sapply(strsplit(kegg_plants, "ko"), function(x) x[[2]]))

ko_list <- ko_list[names(ko_list) %in% kegg_plants]
length(ko_list)

# Generate kegg terms list with phases per kegg term for LD and SD
phases.list.ld <- gene.list.to.phases(ko_list, total_phases_table_ld)
phases.list.sd <- gene.list.to.phases(ko_list, total_phases_table_sd)

# Remove keggs with no associated rhytmic genes
phases.list.ld.clean <- phases.list.ld[which(sapply(phases.list.ld, nrow) != 0)]
phases.list.sd.clean <- phases.list.sd[which(sapply(phases.list.sd, nrow) != 0)]

# Test preferential concentration of kegg terms during the diel cycle
library(circular)

kegg.circa.table.ld <- complete_circular_table(phases.list.ld.clean, hr.on.large.sets = F, hr.on.large.sets.th=200,
                                             iter = 3999, force.hr=F, force.hr.th=0.05)

kegg.circa.table.ld <- complete_circular_table(phase.list = phases.list.ld.clean,
                                               hr.on.large.sets = F , hr.on.large.sets.th = 200,
                                               force.hr = F, force.hr.th = 0.05, 
                                               rao.on.large.sets = F, rao.on.large.sets.th = 200, 
                                               force.rao = F, force.rao.th = 0.05,
                                               iter.rao = 3999, iter.hr = 3999)

write.table(kegg.circa.table.ld,file="tables/kegg_circa_table_ld.tsv", sep="\t")


kegg.circa.table.sd <- complete_circular_table(phase.list = phases.list.sd.clean, 
                                               hr.on.large.sets = F , hr.on.large.sets.th = 200, 
                                               force.hr = F, force.hr.th = 0.05, 
                                               rao.on.large.sets = F, rao.on.large.sets.th = 200, 
                                               force.rao = F, force.rao.th = 0.05,
                                               iter.rao = 3999, iter.hr = 3999)

write.table(kegg.circa.table.sd,file="tables/kegg_circa_table_sd.tsv", sep="\t")

# Read and explore results
kegg.circa.table.ld <- read.table("tables/kegg_circa_table_ld.tsv", sep="\t", header = T)
kegg.circa.table.sd <- read.table("tables/kegg_circa_table_sd.tsv", sep="\t", header = T)

sum(kegg.circa.table.ld$rayleigh_p_value_adj < 0.05 | kegg.circa.table.ld$kuiper_p_value_adj < 0.05) # 73
sum(kegg.circa.table.sd$rayleigh_p_value_adj < 0.05 | kegg.circa.table.sd$kuiper_p_value_adj < 0.05) # 69
sum(kegg.circa.table.ld$hr_p_value_adj < 0.05, na.rm = T) # 2
sum(kegg.circa.table.sd$hr_p_value_adj < 0.05, na.rm = T) # 2
sum(kegg.circa.table.ld$rao_p_value_adj < 0.05, na.rm = T) # 0
sum(kegg.circa.table.sd$rao_p_value_adj < 0.05, na.rm = T) # 0

unimod_ld <- rownames(kegg.circa.table.ld)[which(kegg.circa.table.ld$rayleigh_p_value_adj < 0.05 | kegg.circa.table.ld$kuiper_p_value_adj < 0.05)]
unimod_sd <- rownames(kegg.circa.table.sd)[which(kegg.circa.table.sd$rayleigh_p_value_adj < 0.05 | kegg.circa.table.sd$kuiper_p_value_adj < 0.05)]

multimod_ld <- rownames(kegg.circa.table.ld)[which(kegg.circa.table.ld$hr_p_value_adj < 0.05 | kegg.circa.table.ld$rao_p_value_adj < 0.05)]
multimod_sd <- rownames(kegg.circa.table.sd)[which(kegg.circa.table.sd$hr_p_value_adj < 0.05 | kegg.circa.table.sd$rao_p_value_adj < 0.05)]

library(VennDiagram)
png("ld_sd_unimod_kegg.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(unimod_ld), area2 = length(unimod_sd),
                   cross.area = length(intersect(unimod_ld,unimod_sd)),
                   lwd = 3,catekeggry = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

png("ld_sd_multimod_kegg.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(multimod_ld), area2 = length(multimod_sd),
                   cross.area = length(intersect(multimod_ld,multimod_sd)),
                   lwd = 3,catekeggry = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# For the kegg terms not following an uniform distribution, test them against the general distribution of phases
# per photoperiod
non_uniform_keggs_ld <- unique(c(unimod_ld, multimod_ld))
phases.list.ld.non.uni <- phases.list.ld.clean[non_uniform_keggs_ld]
gen.dist.ld <- test_against_gen_dist(phase.list = phases.list.ld.non.uni, total.phase.table = total_phases_table_ld)
write.table(gen.dist.ld,file="tables/kegg_general_dist_ld.tsv", sep="\t")

non_uniform_keggs_sd <- unique(c(unimod_sd, multimod_sd))
phases.list.sd.non.uni <- phases.list.sd.clean[non_uniform_keggs_sd]
gen.dist.sd <- test_against_gen_dist(phase.list = phases.list.sd.non.uni, total.phase.table = total_phases_table_sd)
write.table(gen.dist.sd,file="tables/kegg_general_dist_sd.tsv", sep="\t")

gen.dist.ld <- read.table("tables/kegg_general_dist_ld.tsv", sep="\t", header = T)
gen.dist.sd <- read.table("tables/kegg_general_dist_sd.tsv", sep="\t", header = T)

sum(gen.dist.ld$gen_dist_p_value_adj < 0.05) # 25 non-uniform pathways present a different distribution than that of global
# mRNA accumulation in LD
sum(gen.dist.sd$gen_dist_p_value_adj < 0.05) # 21 non-uniform pathways present a different distribution than that of global
# mRNA accumulation in SD

kegg_dist_gen_ld <- rownames(gen.dist.ld)[which(gen.dist.ld$gen_dist_p_value_adj < 0.05)]
kegg_dist_gen_sd <- rownames(gen.dist.sd)[which(gen.dist.sd$gen_dist_p_value_adj < 0.05)]

png("ld_sd_general_dist_kegg.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(kegg_dist_gen_ld), area2 = length(kegg_dist_gen_sd),
                   cross.area = length(intersect(kegg_dist_gen_ld,kegg_dist_gen_sd)),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# Test shared kegg terms not following an uniform distribution for differences by photoperiod
shared_non_uniform <- intersect(non_uniform_keggs_ld, non_uniform_keggs_sd)
phases.list.shared.ld.non.uni <- phases.list.ld.non.uni[shared_non_uniform]
phases.list.shared.sd.non.uni <- phases.list.sd.non.uni[shared_non_uniform]

kegg_diff_ld_sd_table <- test_two_dist(phases.list.shared.ld.non.uni, phases.list.shared.sd.non.uni)

write.table(kegg_diff_ld_sd_table,file="tables/kegg_diff_ld_sd_table", sep="\t")

kegg_diff_ld_sd_table <- read.table("tables/kegg_diff_ld_sd_table", sep="\t", header = T)

sum(kegg_diff_ld_sd_table$gen_dist_p_value_adj < 0.05) # 27 change their distribution between LD and SD
change_keggs <- rownames(kegg_diff_ld_sd_table)[which(kegg_diff_ld_sd_table$gen_dist_p_value_adj < 0.05)]

# Difference in medians for changing keggs
diff_keggs_medians <- kegg.circa.table.ld[change_keggs,"median"] - kegg.circa.table.sd[change_keggs,"median"]

sd_anticipated <- change_keggs[which((0 < diff_keggs_medians & diff_keggs_medians < 12) | (diff_keggs_medians < -12))]
length(sd_anticipated)
# All changing keggs are anticipated in SD

library(circlize)
library(ComplexHeatmap)
library(MetBrewer)

circ_ex_ld_kegg1 <- kegg.circa.table.ld[change_keggs[1:14],]
circ_ex_ld_kegg2 <- kegg.circa.table.ld[change_keggs[15:27],]
circular_boxplot(circ_ex_ld_kegg1, "prueba_keggs1_ld.png", color.palette = "Tam",
                 tr.height = 0.055)
circular_boxplot(circ_ex_ld_kegg2, "prueba_keggs2_ld.png", color.palette = "Tam",
                 tr.height = 0.055)

circ_ex_sd_kegg1 <- kegg.circa.table.sd[change_keggs[1:14],]
circ_ex_sd_kegg2 <- kegg.circa.table.sd[change_keggs[15:27],]
circular_boxplot(circ_ex_sd_kegg1, "prueba_keggs1_sd.png", color.palette = "Tam",
                 tr.height = 0.055)
circular_boxplot(circ_ex_sd_kegg2, "prueba_keggs2_sd.png", color.palette = "Tam",
                 tr.height = 0.055)

# Description of changing kegg terms
change_keggs_description <- mapply(function(x) {keggFind("pathway", x)}, change_keggs, USE.NAMES = F)
names(change_keggs_description) <- sapply(strsplit(names(change_keggs_description), ":"), function(x) x[2])

# Add phase medians in each photoperiod
change_keggs_table <- data.frame(kegg=names(change_keggs_description), description=change_keggs_description,
                                  ld_phase=kegg.circa.table.ld[names(change_keggs_description), "median"],
                                  sd_phase=kegg.circa.table.sd[names(change_keggs_description), "median"])
write.table(change_keggs_table,file="tables/anticipated_ld_sd_kegg.tsv", sep="\t", col.names = T)

# Multimodality is not real but a consequence of only one gene in set
library(OptCirClust)
multi_ld_path1 <- multimodal_analysis(phases.list.ld[[multimod_ld[1]]])
multi_ld_path2 <- multimodal_analysis(phases.list.ld[[multimod_ld[2]]])

circular_dotplot(phases.list.ld[multimod_ld[1]], "multimodal_ld1.png", 
                 color.palette = "Tam", tr.height = 0.15)
circular_dotplot(phases.list.sd[multimod_sd[2]], "multimodal_ld2.png", 
                 color.palette = "Tam", tr.height = 0.15)

# Paths that change between photoperiods and are not experimentally distributed
change_keggs[change_keggs %in% kegg_dist_gen_ld]
change_keggs[change_keggs %in% kegg_dist_gen_sd]

# Paths that change between photoperiods and are experimentally distributed
change_keggs[!(change_keggs %in% kegg_dist_gen_ld)]
change_keggs[!(change_keggs %in% kegg_dist_gen_sd)]

# Paths that do not change between photoperiods and are not experimentally distributed
kegg_dist_gen_ld[!(kegg_dist_gen_ld %in% change_keggs)]
kegg_dist_gen_sd[!(kegg_dist_gen_sd %in% change_keggs)]



## Functional enrichment of free-running conditions
complete.ld.ll.rhythmic.genes <- read.table(file="tables/ld.ll.rhythmic.table.tsv", sep="\t",header = T)$x
complete.ld.dd.rhythmic.genes <- read.table(file="tables/ld.dd.rhythmic.table.tsv", sep="\t",header = T)$x
complete.sd.ll.rhythmic.genes <- read.table(file="tables/sd.ll.rhythmic.table.tsv", sep="\t",header = T)$x
complete.sd.dd.rhythmic.genes <- read.table(file="tables/sd.dd.rhythmic.table.tsv", sep="\t",header = T)$x

inter_ll <- intersect(complete.ld.ll.rhythmic.genes, complete.sd.ll.rhythmic.genes)
inter_dd <- intersect(complete.ld.dd.rhythmic.genes, complete.sd.dd.rhythmic.genes)
union_ll <- unique(c(complete.ld.ll.rhythmic.genes, complete.sd.ll.rhythmic.genes))
union_dd <- unique(c(complete.ld.dd.rhythmic.genes, complete.sd.dd.rhythmic.genes))

complete.ld.rhythmic.genes <- read.table(file="tables/ld.rhythmic.table.tsv",header = T)$V1
complete.sd.rhythmic.genes <- read.table(file="tables/sd.rhythmic.table.tsv",header = T)$V1

# For each one
library(clusterProfiler)
library("org.Mpolymorpha71.eg.db")
universe_set <- unique(c(complete.ld.rhythmic.genes, complete.sd.rhythmic.genes))

ego.inter.ll <- enrichGO(gene = inter_ll,
                         OrgDb = org.Mpolymorpha71.eg.db,
                         universe = universe_set,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = F,
                         keyType = "GID")

barplot(ego.inter.ll, showCategory = 20)
ego.inter.ll.df <- as.data.frame(ego.inter.ll)

ego.inter.dd <- enrichGO(gene = inter_dd,
                         OrgDb         = org.Mpolymorpha71.eg.db,
                         universe = universe_set,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = F,
                         keyType = "GID")

ego.union.ll <- enrichGO(gene = union_ll,
                         OrgDb         = org.Mpolymorpha71.eg.db,
                         universe = universe_set,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = F,
                         keyType = "GID")

barplot(ego.union.ll, showCategory = 20)
ego.union.ll.df <- as.data.frame(ego.union.ll)

ego.union.dd <- enrichGO(gene = union_dd,
                         OrgDb         = org.Mpolymorpha71.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = F,
                         keyType = "GID")
barplot(ego.union.dd, showCategory = 20)
ego.inter.ll.df <- as.data.frame(ego.inter.ll)


# Chloroplasts photorelocation plots
source("functions_circa_marchantia.R")

# Mp5g03810 PHOT
# Mp7g03540	CHUP1
# Mp3g03760	KAC
# Mp8g10870 BLUS1
# Mp1g07590	THRUMIN1
# Mp3g04860	THRUMIN1
# Mp6g17960 Mp6g06160 Mp3g09960 JAC1
# Mp5g14210 THRUMIN1 no glutaredoxin 
# Mp2g25860 PMI1 related
# "Mp5g06760" "Mp8g10000" "Mp8g18590" "Mp8g18700" PMI related
# No WEB1

# Photoreceptor PHOT and its substrate BLUS1
"Mp5g03810" %in% complete.ld.rhythmic.genes
"Mp5g03810" %in% complete.sd.rhythmic.genes

"Mp8g10870" %in% complete.ld.rhythmic.genes
"Mp8g10870" %in% complete.sd.rhythmic.genes

plot.ld.sd("Mp5g03810", "PHOT", genes.expression)

# Motility CHUP1, KAC, THRUM1 and PMI1
c("Mp7g03540", "Mp3g03760", "Mp1g07590", "Mp3g04860", "Mp2g25860") %in% complete.ld.rhythmic.genes
c("Mp7g03540", "Mp3g03760", "Mp1g07590", "Mp3g04860", "Mp2g25860") %in% complete.sd.rhythmic.genes

plot.ld.sd("Mp3g03760", "KAC", genes.expression)
plot.ld.sd("Mp7g03540", "CHUP1", genes.expression)

# Signaling JAC1, WEB1, and PMI2
c("Mp6g17960", "Mp6g06160", "Mp3g09960", "Mp5g06760", "Mp8g10000", "Mp8g18590", "Mp8g18700") %in% complete.ld.rhythmic.genes
c("Mp6g17960", "Mp6g06160", "Mp3g09960", "Mp5g06760", "Mp8g10000", "Mp8g18590", "Mp8g18700") %in% complete.sd.rhythmic.genes

plot.ld.sd("Mp5g06760", "PMIr", genes.expression)
png("phot.png")
plot.ld.sd.wave("Mp5g03810", "PHOT",  genes.expression, wave_parameters =  circa_table_genescarab["Mp5g03810",])
dev.off()

plot.ld.sd.wave("Mp2g14310", "RVE1", genes.expression, wave_parameters =  circa_table_genescarab["Mp2g14310",])
plot.ld.ll("Mp2g14310", "RVE1", genes.expression)
plot.sd.ll("Mp2g14310", "RVE1", genes.expression)

# Create rhythmic plots
circa_table_genescarab <- read.table("tables/circacompare_intersection.tsv", header = T, sep = "\t")
head(circa_table_genescarab)
genes.expression
ld.sd.gene.expression <- genes.expression[,c(ld.zt.i,sd.zt.i)]
plot.ld.sd("Mp1g00120", "", genes.expression)

parameters_example <- circa_table_genescarab["Mp1g00120",]
data_example <- ld.sd.gene.expression["Mp1g00120",]


wave.ld <- sapply(seq(from=0,by=0.1,to=24), function(x) (12 + pi*cos(24*(x-3))))
wave.ld <- sapply(seq(from=0,by=0.1,to=68), function(x) 12+3*cos(1*(x*pi/12-3*pi/12)))

plot(x = seq(from=0,by=0.1,to=68),wave.ld,lwd=3,col="blue", type="l")

plot.ld.sd.wave(gene.id = "Mp1g00120", gene.name = "", genes.expression)


plot.ld.sd.wave(gene.id = "Mp1g00120", gene.name = "", genes.expression, wave_parameters = parameters_example)

# LD SD
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

# SD LL
circa_table_sd_ll <- read.table("tables/circacompare_sd_ll.tsv", header = T, sep = "\t")
head(circa_table_sd_ll)
plot.sd.ll.wave(gene.id = "Mp1g00120", gene.name = "", genes.expression, wave_parameters = circa_table_sd_ll["Mp1g00120",])
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


# ld LL
circa_table_ld_ll <- read.table("tables/circacompare_ld_ll.tsv", header = T, sep = "\t")
head(circa_table_ld_ll)

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

plot.ld.ll.wave(gene.id = "Mp1g00160", gene.name = "", genes.expression, wave_parameters = circa_table_ld_ll["Mp1g00160",])

# ld dd
circa_table_ld_dd <- read.table("tables/circacompare_ld_dd.tsv", header = T, sep = "\t")
head(circa_table_ld_dd)

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

plot.ld.dd.wave(gene.id = "Mp1g00190", gene.name = "", genes.expression, wave_parameters = circa_table_ld_dd["Mp1g00190",])


# sd dd
circa_table_sd_dd <- read.table("tables/circacompare_sd_dd.tsv", header = T, sep = "\t")
head(circa_table_sd_dd)

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

plot.sd.dd.wave(gene.id = "Mp1g02630", gene.name = "", genes.expression, wave_parameters = circa_table_sd_dd["Mp1g02630",])


# Study of TFs

# Load names of Marchantia's TFs
tf_list <- unique(c("Mp1g00510",
"Mp1g01010",
"Mp1g01110",
"Mp1g01160",
"Mp1g01900",
"Mp1g01960",
"Mp1g02320",
"Mp1g02860",
"Mp1g02890",
"Mp1g03040",
  "Mp1g03100",
"Mp1g03200",
"Mp1g03580",
"Mp1g03950",
"Mp1g04480",
"Mp1g04550",
"Mp1g05310",
"Mp1g05400",
"Mp1g05630",
"Mp1g06180",
"Mp1g07040",
"Mp1g07070",
"Mp1g08640",
"Mp1g08660",
"Mp1g08890",
"Mp1g08960",
"Mp1g09040",
"Mp1g09420",
"Mp1g09500",
"Mp1g10020",
"Mp1g10030",
"Mp1g10300",
"Mp1g10310",
"Mp1g10330",
"Mp1g10440",
"Mp1g11660",
"Mp1g11830",
"Mp1g12700",
"Mp1g12750",
"Mp1g13010",
"Mp1g13260",
"Mp1g13640",
"Mp1g13740",
"Mp1g14500",
"Mp1g15660",
"Mp1g16790",
"Mp1g16800",
"Mp1g16880",
"Mp1g17210",
"Mp1g18540",
"Mp1g19210",
"Mp1g19590",
"Mp1g19730",
"Mp1g20040",
"Mp1g20490",
"Mp1g21450",
"Mp1g22880",
"Mp1g23310",
"Mp1g24140",
"Mp1g24520",
"Mp1g24950",
"Mp1g25020",
"Mp1g25140",
"Mp1g25720",
"Mp1g26740",
"Mp1g27060",
"Mp1g28480",
"Mp1g28490",
"Mp1g29240",
"Mp2g00500",
"Mp2g00510",
"Mp2g00890",
"Mp2g00910",
"Mp2g00930",
"Mp2g02150",
"Mp2g02230",
"Mp2g02280",
"Mp2g02990",
"Mp2g03850",
"Mp2g04160",
"Mp2g04180",
"Mp2g04190",
"Mp2g04200",
"Mp2g05060",
"Mp2g05070",
"Mp2g05160",
"Mp2g06450",
"Mp2g06460",
"Mp2g07110",
"Mp2g07130",
"Mp2g07150",
"Mp2g07170",
"Mp2g07720",
"Mp2g08260",
"Mp2g08360",
"Mp2g08790",
"Mp2g09320",
"Mp2g10620",
"Mp2g11080",
"Mp2g11220",
"Mp2g11330",
"Mp2g11520",
"Mp2g12750",
  "Mp2g13020",
"Mp2g13330",
"Mp2g13910",
"Mp2g14310",
"Mp2g14450",
"Mp2g16340",
"Mp2g16350",
"Mp2g19970",
  "Mp2g20695",
"Mp2g20790",
"Mp2g20960",
"Mp2g20990",
"Mp2g22110",
"Mp2g22250",
"Mp2g22560",
"Mp2g22820",
"Mp2g23000",
"Mp2g23600",
"Mp2g24110",
"Mp2g24200",
"Mp2g25260",
"Mp2g25420",
"Mp2g25660",
"Mp2g25840",
"Mp2g26400",
"Mp3g00010",
"Mp3g00360",
"Mp3g00600",
"Mp3g01640",
"Mp3g02320",
"Mp3g02770",
"Mp3g02830",
"Mp3g04030",
"Mp3g04360",
"Mp3g04970",
"Mp3g05000",
"Mp3g05910",
"Mp3g06860",
"Mp3g07510",
"Mp3g07610",
"Mp3g09110",
"Mp3g09350",
"Mp3g09590",
"Mp3g10500",
"Mp3g10830",
"Mp3g10970",
"Mp3g11570",
"Mp3g11900",
"Mp3g12060",
"Mp3g12980",
"Mp3g16240",
"Mp3g17260",
"Mp3g17350",
"Mp3g17660",
"Mp3g17930",
"Mp3g18590",
"Mp3g19160",
"Mp3g19510",
"Mp3g19580",
"Mp3g19670",
"Mp3g20870",
"Mp3g20900",
"Mp3g20980",
"Mp3g20990",
"Mp3g21340",
"Mp3g21490",
"Mp3g23170",
"Mp3g23300",
"Mp3g23720",
"Mp3g24320",
"Mp3g24780",
"Mp3g25200",
"Mp4g00180",
"Mp4g00200",
"Mp4g00280",
"Mp4g00380",
"Mp4g01170",
"Mp4g01540",
"Mp4g01560",
"Mp4g02060",
"Mp4g02620",
"Mp4g04750",
"Mp4g04760",
"Mp4g04910",
"Mp4g04920",
"Mp4g06170",
"Mp4g06490",
"Mp4g06740",
"Mp4g07330",
"Mp4g07910",
"Mp4g08690",
"Mp4g08700",
"Mp4g09650",
"Mp4g09730",
"Mp4g09860",
"Mp4g09880",
"Mp4g10240",
"Mp4g10500",
"Mp4g10710",
"Mp4g11030",
"Mp4g11380",
"Mp4g11760",
"Mp4g11820",
"Mp4g11860",
"Mp4g11910",
"Mp4g12230",
"Mp4g12420",
"Mp4g12490",
"Mp4g13360",
"Mp4g13860",
"Mp4g13880",
"Mp4g14520",
"Mp4g14530",
"Mp4g15170",
"Mp4g15240",
"Mp4g15460",
"Mp4g15560",
"Mp4g17430",
"Mp4g18230",
"Mp4g19650",
"Mp4g20600",
"Mp4g21220",
"Mp4g21300",
"Mp4g21790",
"Mp4g22280",
"Mp4g22290",
"Mp4g22310",
"Mp4g22430",
"Mp4g22890",
"Mp5g00110",
"Mp5g01050",
"Mp5g01530",
"Mp5g01600",
"Mp5g02590",
"Mp5g03020",
"Mp5g05560",
"Mp5g05940",
"Mp5g06970",
"Mp5g07255",
"Mp5g07890",
"Mp5g08310",
"Mp5g08410",
"Mp5g08810",
"Mp5g09130",
"Mp5g09220",
"Mp5g09660",
"Mp5g09710",
"Mp5g10280",
"Mp5g10960",
"Mp5g11060",
"Mp5g11100",
"Mp5g11460",
"Mp5g11830",
"Mp5g12220",
"Mp5g12480",
"Mp5g13030",
"Mp5g13040",
"Mp5g13060",
"Mp5g13130",
"Mp5g13160",
"Mp5g14610",
"Mp5g14820",
"Mp5g17820",
"Mp5g18910",
"Mp5g19050",
"Mp5g20280",
"Mp5g20530",
"Mp5g20660",
"Mp5g21080",
"Mp5g21990",
"Mp5g22160",
"Mp5g22280",
"Mp5g23860",
"Mp6g00540",
"Mp6g02620",
"Mp6g02670",
"Mp6g02950",
"Mp6g03670",
"Mp6g03920",
"Mp6g04430",
"Mp6g04650",
"Mp6g04830",
"Mp6g04880",
"Mp6g05720",
"Mp6g05810",
"Mp6g06850",
"Mp6g06870",
"Mp6g07410",
"Mp6g07440",
"Mp6g08290",
"Mp6g08690",
"Mp6g09270",
"Mp6g09280",
"Mp6g10210",
"Mp6g10570",
"Mp6g10690",
"Mp6g11770",
"Mp6g13820",
"Mp6g15470",
"Mp6g16170",
"Mp6g16760",
"Mp6g16800",
"Mp6g17420",
"Mp6g18480",
"Mp6g18650",
"Mp6g18860",
"Mp6g20920",
"Mp6g21430",
"Mp6g21470",
"Mp7g00210",
"Mp7g00650",
"Mp7g00860",
"Mp7g01800",
"Mp7g02220",
"Mp7g02640",
"Mp7g03150",
"Mp7g03490",
"Mp7g03700",
"Mp7g06280",
"Mp7g06290",
"Mp7g06550",
"Mp7g07550",
"Mp7g09260",
"Mp7g09350",
"Mp7g09490",
"Mp7g09710",
"Mp7g09740",
"Mp7g11330",
"Mp7g11500",
"Mp7g12630",
"Mp7g12800",
"Mp7g13270",
"Mp7g13760",
"Mp7g14310",
"Mp7g14810",
"Mp7g15840",
"Mp7g16300",
"Mp7g16910",
"Mp7g17020",
"Mp7g17200",
"Mp7g17250",
"Mp7g18530",
"Mp7g18640",
"Mp7g18690",
"Mp7g19570",
"Mp8g01190",
"Mp8g01770",
"Mp8g02290",
"Mp8g03820",
"Mp8g03980",
"Mp8g04010",
"Mp8g04130",
"Mp8g04330",
"Mp8g04400",
"Mp8g06180",
"Mp8g06230",
"Mp8g06450",
"Mp8g07390",
"Mp8g08100",
"Mp8g08190",
"Mp8g09250",
"Mp8g09500",
"Mp8g10640",
"Mp8g11450",
"Mp8g11560",
"Mp8g11850",
"Mp8g11870",
"Mp8g12000",
"Mp8g12100",
"Mp8g12510",
"Mp8g12800",
"Mp8g13205",
"Mp8g13470",
"Mp8g14220",
"Mp8g15500",
"Mp8g16170",
"Mp8g16900",
"Mp8g18080",
"Mp8g18310",
"Mp8g18820",
"MpUg00060",
"MpUg00230",
"MpUg00270",
"MpUg00290",
"MpUg00370",
"MpVg00300",
"MpVg00340",
"MpVg00350",
"MpVg00470"))

# Rhythmic TFs
circa_table_tf_inter <- read.table("tables/circacompare_intersection.tsv", header = T, sep = "\t")
cycle_inter_tfs_table <- circa_table_tf_inter[rownames(circa_table_tf_inter) %in% tf_list,]
write.table(cycle_inter_tfs_table, "tfs_ld_sd_circacompare_inter.tsv", col.names = T, sep="\t")

# Plot TFs
sorted_cycle <- sort_by(cycle_inter_tfs_table, cycle_inter_tfs_table$ld.peak.time.hours)
k_names <- c(paste0(rownames(sorted_cycle),"_ld"),
             paste0(rownames(sorted_cycle),"_sd"))

k_values <- c(sorted_cycle$ld.peak.time.hours, 
              sorted_cycle$sd.peak.time.hours)

k_table <- data.frame(names = k_names, values = k_values,
                      cat = factor(rownames(sorted_cycle), levels = rownames(sorted_cycle)), 
                      photo = rep(c("ld", "sd"), each= nrow(sorted_cycle)))
library(ggplot2)

png("ordered_tfs.png", width = 480, height = 1440, res=100)
k_plot <- ggplot(k_table, aes(x = values, y = cat)) +
  #geom_line(aes(color = allowed)) +
  geom_line() +
  geom_point(aes(color = photo), size = 3) +
  scale_color_brewer(palette = "Set1", direction = -1) +
  theme_bw() +
  theme(legend.position = "bottom")
plot(k_plot)
dev.off()

# Clock genes
plot.ld.sd("Mp3g09460", "TOC1", genes.expression)
plot.ld.ll("Mp3g09460", "TOC1", genes.expression)
plot.sd.ll("Mp3g09460", "TOC1", genes.expression)
plot.ld.sd("Mp2g16560", "PRR", genes.expression)
plot.ld.ll("Mp2g16560", "PRR", genes.expression)
plot.sd.ll("Mp2g16560", "PRR", genes.expression)
plot.ld.sd("Mp8g00970", "DET1", genes.expression)
plot.ld.ll("Mp8g00970", "DET1", genes.expression)
plot.sd.ll("Mp8g00970", "DET1", genes.expression)
plot.ld.dd("Mp8g00970", "DET1", genes.expression)
plot.ld.sd("Mp2g16350", "COL2", genes.expression)
plot.ld.sd("Mp2g14310", "RVE", genes.expression)
plot.ld.sd("Mp1g10620", "LNK", genes.expression)
plot.ld.sd("Mp1g13750", "GI", genes.expression)
plot.ld.sd("Mp5g05940", "COL6", genes.expression)
plot.ld.sd("Mp1g10620", "LNK", genes.expression)
plot.ld.sd("Mp1g13750", "GI", genes.expression)
plot.ld.sd("Mp5g01830", "WDR2", genes.expression)
plot.ld.sd("Mp1g18540", "CDF1", genes.expression)


plot.ld.sd.wave("Mp2g24110", "IDD3",  genes.expression, wave_parameters =  circa_table_genescarab["Mp2g24110",])
plot.ld.ll("Mp2g24110", "", genes.expression)
plot.sd.ll("Mp2g24110", "", genes.expression)
plot.ld.ll("Mp1g10620", "LNK", genes.expression)
plot.sd.ll("Mp1g10620", "LNK", genes.expression)

# Potential regulation of DET1
plot.ld.sd("Mp1g15660", "SGR5", genes.expression) # not rhythmic

plot.ld.sd("Mp8g12510", "IDD5", genes.expression) # change between photoperiods
plot.ld.sd.wave("Mp8g12510", "IDD5",  genes.expression, wave_parameters =  circa_table_genescarab["Mp8g12510",])
plot.ld.ll("Mp8g12510", "IDD5", genes.expression)
plot.sd.ll("Mp8g12510", "IDD5", genes.expression)

plot.ld.sd("Mp5g20530", "IDD3",  genes.expression) # at1g14580, MpIDDL3
plot.ld.sd.wave("Mp5g20530", "IDD3",  genes.expression, wave_parameters =  circa_table_genescarab["Mp5g20530",]) # good candidate
plot.ld.ll("Mp5g20530", "", genes.expression)
plot.sd.ll("Mp5g20530", "", genes.expression)

plot.ld.sd("Mp2g24110", "C2H2-10",  genes.expression) # at1g14580, MpC2H2-10
plot.ld.sd.wave("Mp2g24110", "C2H2-10",  genes.expression, wave_parameters =  circa_table_genescarab["Mp2g24110",]) # good candidate
plot.ld.ll("Mp2g24110", "", genes.expression)
plot.sd.ll("Mp2g24110", "", genes.expression) # except for second day in SD-LL

plot.ld.sd("Mp4g15170", "MpGARP6",  genes.expression) # AT1G68670, MpGARP6
plot.ld.sd.wave("Mp4g15170", "MpGARP6",  genes.expression, wave_parameters =  circa_table_genescarab["Mp4g15170",]) # incorrect phase
plot.ld.ll("Mp4g15170", "", genes.expression)
plot.sd.ll("Mp4g15170", "", genes.expression) 

plot.ld.sd("Mp4g01560", "MpGARP3",  genes.expression) # AT3G12730, MpGARP3
plot.ld.sd.wave("Mp4g01560", "MpGARP3",  genes.expression, wave_parameters =  circa_table_genescarab["Mp4g01560",]) # incorrect phase
plot.ld.ll("Mp4g01560", "", genes.expression)
plot.sd.ll("Mp4g01560", "", genes.expression) 

plot.ld.sd("Mp1g15660", "IDD2",  genes.expression) # MpIDDL2, not rhythmic
plot.ld.sd.wave("Mp1g15660", "IDD2",  genes.expression, wave_parameters =  circa_table_genescarab["Mp1g15660",])
plot.ld.ll("Mp1g15660", "", genes.expression)
plot.sd.ll("Mp1g15660", "", genes.expression) 

plot.ld.sd("MpVg00300", "ZML2",  genes.expression) # ZML2, MpGATA6
plot.ld.sd.wave("MpVg00300", "ZML2",  genes.expression, wave_parameters =  circa_table_genescarab["MpVg00300",]) # incorrect phase
plot.ld.ll("MpVg00300", "", genes.expression)
plot.sd.ll("MpVg00300", "", genes.expression) 

plot.ld.sd("Mp4g12420", "BIM1",  genes.expression) # BIM1, MpBHLH44
plot.ld.sd.wave("Mp4g12420", "BIM1",  genes.expression, wave_parameters =  circa_table_genescarab["Mp4g12420",]) # incorrect phase
plot.ld.ll("Mp4g12420", "", genes.expression)
plot.sd.ll("Mp4g12420", "", genes.expression)

plot.ld.sd("MpVg00350", "BPC1",  genes.expression) # BPC1, MpBPCV, not rhythmic
plot.ld.sd.wave("MpVg00350", "BPC1",  genes.expression, wave_parameters =  circa_table_genescarab["MpVg00350",])
plot.ld.ll("MpVg00350", "", genes.expression)
plot.sd.ll("MpVg00350", "", genes.expression)

plot.ld.sd("MpVg00470", "MpHSF3",  genes.expression) # HSFA1E, MpHSF3
plot.ld.sd.wave("MpVg00470", "MpHSF3",  genes.expression, wave_parameters =  circa_table_genescarab["MpVg00470",]) # incorrect phase
plot.ld.ll("MpVg00470", "", genes.expression)
plot.sd.ll("MpVg00470", "", genes.expression)

plot.ld.sd("Mp4g12230", "MpHSF1",  genes.expression) # HSFB3, MpHSF1
plot.ld.sd.wave("Mp4g12230", "MpHSF1",  genes.expression, wave_parameters =  circa_table_genescarab["Mp4g12230",]) # incorrect phase
plot.ld.ll("Mp4g12230", "", genes.expression)
plot.sd.ll("Mp4g12230", "", genes.expression)

# HY5
plot.ld.sd("Mp1g16800", "MpHY5",  genes.expression)
plot.ld.sd.wave("Mp1g16800", "MpHY5",  genes.expression, wave_parameters =  circa_table_genescarab["Mp1g16800",]) # incorrect phase
plot.ld.ll("Mp1g16800", "MpHY5", genes.expression)
plot.sd.ll("Mp1g16800", "MpHY5", genes.expression)



# Circadian TFs
length(tf_list)
circa_table_sd_dd <- read.table("tables/circacompare_sd_dd.tsv", header = T, sep = "\t")
circa_table_ld_dd <- read.table("tables/circacompare_ld_dd.tsv", header = T, sep = "\t")
circa_table_sd_ll <- read.table("tables/circacompare_sd_ll.tsv", header = T, sep = "\t")
circa_table_ld_ll <- read.table("tables/circacompare_ld_ll.tsv", header = T, sep = "\t")

tfs_ld_ll <- tf_list[tf_list %in% rownames(circa_table_ld_ll)]
tfs_sd_ll <- tf_list[tf_list %in% rownames(circa_table_sd_ll)]
tfs_ld_dd <- tf_list[tf_list %in% rownames(circa_table_ld_dd)]
tfs_sd_dd <- tf_list[tf_list %in% rownames(circa_table_sd_dd)]

length(tfs_sd_dd) # 13 ld ll, 37 sd ll, 24 ld dd, 5 sd dd

intersect(tfs_ld_ll, tfs_sd_ll)
intersect(tfs_ld_dd, tfs_sd_dd)
# SD DD bzip, nac4, c3h, c3h, c3h -> C3H se une a ZF, aparece
# no usar Known results.

c("Mp2g19970", "Mp3g18590", "Mp5g05940") %in% tfs_ld_dd 
c("Mp2g19970", "Mp3g18590", "Mp5g05940") %in% tfs_sd_dd

# "Mp2g19970" and Mp3g18590" cycle in LL in both photoperiods, but not in DD
# They are MpBHLH45 and MpTRIHELIX32

intersect(tfs_ld_ll, tfs_ld_dd) # trihelix31, bbx1
intersect(tfs_sd_ll, tfs_sd_dd) # no one
c("Mp4g09730", "Mp5g05940") %in% tfs_sd_ll # bbx1 also in SD-LL

# Enriched motifs
# Enriched in LL with respect to both cycle and free running conditions
ll_vs_complete_motif_table <- create_motif_table("promoters/homer_results/ll_inter_vs_complete/nonRedundant.motifs")
ll_vs_free_running <- create_motif_table("promoters/homer_results/ll_inter_vs_free_running/nonRedundant.motifs")

intersect(ll_vs_complete_motif_table$motif, ll_vs_free_running$motif)

# IDD11 motif is 
c("Mp2g24110", "Mp5g20530", "Mp1g15660","Mpzg00170") %in% tfs_ld_ll
c("Mp2g24110", "Mp5g20530", "Mp1g15660","Mpzg00170") %in% tfs_sd_ll
c("Mp2g24110", "Mp5g20530", "Mp1g15660","Mpzg00170") %in% tfs_ld_dd
c("Mp2g24110", "Mp5g20530", "Mp1g15660","Mpzg00170") %in% tfs_sd_dd

# Look for TFs that cycle in LL but not in DD
setdiff(intersect(tfs_sd_ll, tfs_ld_ll), union(tfs_sd_dd, tfs_ld_dd)) # bhlh45, trihelix32

setdiff(union(tfs_sd_dd, tfs_ld_dd), union(tfs_sd_ll, tfs_ld_ll)) # "Mp2g05070" is bhlh

"Mp1g24140" %in% union(tfs_sd_ll, tfs_ld_ll) # puede serlo, homeobox ATHB
# Phytochrome-interacting factor 4 and 5 (PIF4 and PIF5) activate the homeobox ATHB2 and auxin-inducible 
# IAA29 genes in the coincidence mechanism underlying photoperiodic control of plant growth of Arabidopsis thaliana
plot.sd.ll("Mp3g09460", "", genes.expression)

"Mp3g09460" %in% union(tfs_sd_ll, tfs_ld_ll)

# DD
dd_vs_complete_motif_table <- create_motif_table("promoters/homer_results/dd_inter_vs_complete/nonRedundant.motifs")
dd_vs_free_running <- create_motif_table("promoters/homer_results/dd_inter_vs_free_running/nonRedundant.motifs")
intersect(dd_vs_complete_motif_table$motif, dd_vs_free_running$motif)


# Circacompare for a less stringent test for rhythmicity in TFs
# LD to LL
library(circacompare)

circacompare.ld.ll.tf <- matrix(nrow=length(tf_list),ncol=15)
rownames(circacompare.ld.ll.tf) <- tf_list

# For RVE1
circacompare.ll.rve1 <- matrix(nrow=2,ncol=15)
gene.i <- "Mp2g14310"

sd.expression.i <- rep(genes.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"), paste(sd.zt,3,sep="_"))])
ll.expression.i <- genes.expression[gene.i,c(paste(sd.zt,2,sep="_"),paste(ll.sd.zt,1,sep="_"),paste(ll.sd.zt,2,sep="_"))]
time.points <- seq(from=0,by=4,length.out = 18)
sd.ll.df <- data.frame(time=c(time.points,time.points),
                       measure=c(sd.expression.i, ll.expression.i),
                       group=c(rep("cycle",18),rep("ll",18)))
out.i <- circacompare(x = sd.ll.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
circacompare.ll.rve1[1,] <- out.i[[2]][,2]
out.i[[1]]
plot.sd.ll("Mp2g14310", "RVE1", genes.expression)

ld.expression.i <- rep(genes.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"), paste(ld.zt,3,sep="_"))])
ll.expression.i <- genes.expression[gene.i,c(paste(ld.zt,2,sep="_"),paste(ll.ld.zt,1,sep="_"),paste(ll.ld.zt,2,sep="_"))]
time.points <- seq(from=0,by=4,length.out = 18)
ld.ll.df <- data.frame(time=c(time.points,time.points),
                       measure=c(ld.expression.i, ll.expression.i),
                       group=c(rep("cycle",18),rep("ll",18)))
out.i <- circacompare(x = ld.ll.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 0.05)
circacompare.ll.rve1[2,] <- out.i[[2]][,2]
colnames(circacompare.ll.rve1) <- out.i[[2]][,1]
out.i[[1]]
rownames(circacompare.ll.rve1) <- c("SD", "LD")
plot.ld.ll("Mp2g14310", "RVE1", genes.expression)







for(i in 1:length(tf_list))
{
  gene.i <- tf_list[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"))]
  ll.expression.i <- genes.expression[gene.i,c(paste(ll.ld.zt,1,sep="_"),paste(ll.ld.zt,2,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 12)
  
  ld.ll.df <- data.frame(time=c(time.points,time.points),
                         measure=c(ld.expression.i, ll.expression.i),
                         group=c(rep("ld",12),rep("ll",12)))
  
  out.i <- circacompare(x = ld.ll.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  circacompare.ld.ll.tf[i,] <- out.i[[2]][,2]
}

colnames(circacompare.ld.ll) <- out.i[[2]][,1]
write.table(circacompare.ld.ll,file="circacompare_ld_ll.tsv", sep="\t")



# Go enrichment for free running
circa_table_sd_dd <- read.table("tables/circacompare_sd_dd.tsv", header = T, sep = "\t")
circa_table_ld_dd <- read.table("tables/circacompare_ld_dd.tsv", header = T, sep = "\t")
circa_table_sd_ll <- read.table("tables/circacompare_sd_ll.tsv", header = T, sep = "\t")
circa_table_ld_ll <- read.table("tables/circacompare_ld_ll.tsv", header = T, sep = "\t")
circa_table_ld_sd_union <- read.table("tables/circacompare_union.tsv", header = T, sep = "\t")
ll_union <- union(rownames(circa_table_sd_ll), rownames(circa_table_ld_ll))
ll_inter <- intersect(rownames(circa_table_sd_ll), rownames(circa_table_ld_ll))
dd_union <- union(rownames(circa_table_sd_dd), rownames(circa_table_ld_dd))
dd_inter <- intersect(rownames(circa_table_sd_dd), rownames(circa_table_ld_dd))
bk_cycle <- rownames(circa_table_ld_sd_union)

# GO enrichment for low-expressed genes
library("clusterProfiler")

#Install our annotation package for Marchantia polymorpha
#install.packages("./org.Mpolymorpha.eg.db/", repos=NULL, type="source")
library("org.Mpolymorpha.eg.db")
columns(org.Mpolymorpha.eg.db)

ego.results <- lapply(list(ll_union=ll_union, ll_inter=ll_inter, 
                           dd_union=dd_union, dd_inter=dd_inter),
       function(x) enrichGO(gene = x, #universe = bk_cycle,
                OrgDb         = org.Mpolymorpha.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE,
                keyType = "GID")
) 


barplot(ego.results$ll_union)
lapply(ego.results, function(x) as.data.frame(x))
































################################################## 
library(circacompare)
rain.ld.rhythmic.table <- read.table("ld.rhytmic.table.tsv",header=T)
rain.ld.rhythmic.genes <- rain.ld.rhythmic.table$V1
length(rain.ld.rhythmic.genes)
rain.sd.rhythmic.table <- read.table("sd.rhytmic.table.tsv",header=T)
rain.sd.rhythmic.genes <- rain.sd.rhythmic.table$V1
length(rain.sd.rhythmic.genes)

genes.expression <- read.table(file ="marchantia_gene_expression.tsv",sep="\t",header=T,as.is=T)
head(genes.expression)

#### Rhythmic genes in LD and SD (as determined by RAIN)
both.ld.sd.rhythmic.genes <- intersect(rain.ld.rhythmic.genes,rain.sd.rhythmic.genes)
length(both.ld.sd.rhythmic.genes)

circacompare.ld.sd <- matrix(nrow=length(both.ld.sd.rhythmic.genes),ncol=15)
rownames(circacompare.ld.sd) <- both.ld.sd.rhythmic.genes

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

#Construct a table with the circadian data for all the genes rhythmic in both LD and SD conditions
for(i in 1:length(both.ld.sd.rhythmic.genes))
{
  gene.i <- both.ld.sd.rhythmic.genes[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"), paste(ld.zt,3,sep="_"))]
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"), paste(ld.zt,3,sep="_"))]
  ld.sd.expression.i <- cbind(ld.expression.i,sd.expression.i)
  ld.sd.vector <- c()
  for (j in 1:ncol(ld.sd.expression.i))
  {
    ld.sd.vector <- c(ld.sd.vector,ld.sd.expression.i[1,j])
  }
  
  time.points <- seq(from=0,by=4,length.out = 18)
  
  ld.sd.df <- data.frame(time=c(time.points,time.points),
                         measure=ld.sd.vector,
                         group=c(rep("ld",18),rep("sd",18)))
  
  out.i <- circacompare(x = ld.sd.df, col_time = "time", col_group = "group",
                        col_outcome = "measure",alpha_threshold = 1 )
  circacompare.ld.sd[i,] <- out.i[[2]][,2]
}

colnames(circacompare.ld.sd) <- out.i[[2]][,1]
write.table(circacompare.ld.sd,file="circacompare_intersection.tsv", sep="\t")

#Amplitude
par(lwd=3)
boxplot(circacompare.ld.sd[,"ld amplitude estimate" ],
        circacompare.ld.sd[,"sd amplitude estimate" ],outline=F,
        col=c("royalblue3","tomato3"),names=c("LD","SD"),cex.axis=1.5, main="Distribution of amplitudes")

sum(circacompare.ld.sd[,"Amplitude difference estimate" ] < 0)
sum(circacompare.ld.sd[,"Amplitude difference estimate" ] < 0) / length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"P-value for amplitude difference" ] < 0.05)
sum(circacompare.ld.sd[,"P-value for amplitude difference" ] < 0.05)/length(both.ld.sd.rhythmic.genes)

#Mesor
boxplot(circacompare.ld.sd[,"ld mesor estimate" ],circacompare.ld.sd[,"sd mesor estimate" ],
        outline=F,col=c("royalblue3","tomato3"),names=c("LD","SD"),cex.axis=1.5, main="Distribution of mesor")
sum(circacompare.ld.sd[,"Mesor difference estimate" ] < 0) /length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"Mesor difference estimate" ] > 0) /length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"P-value for mesor difference" ] < 0.05)/length(both.ld.sd.rhythmic.genes)


#Phase
boxplot(circacompare.ld.sd[,"ld peak time hours" ],circacompare.ld.sd[,"sd peak time hours" ],
        outline=F,col=c("royalblue3","tomato3"),names=c("LD","SD"),cex.axis=1.5, main="Distribution of phase")
median(circacompare.ld.sd[,"ld peak time hours" ])-median(circacompare.ld.sd[,"sd peak time hours" ])
sum(circacompare.ld.sd[,"Phase difference estimate" ] < 0) /length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"Phase difference estimate" ] > 0) /length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"P-value for difference in phase" ] < 0.05)/length(both.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"P-value for difference in phase" ] < 0.05)

# This should be circular
custom_set <- list(circacompare.ld.sd[,"ld peak time hours" ],
                   circacompare.ld.sd[,"sd peak time hours" ]) # Load vector of phases
ind_table <- sapply(custom_set, function(x) create_circular_table(x, iter=10))
ind_table2 <- as.data.frame(t(ind_table))
ind_table2 <- matrix(c(as.numeric(ind_table[,2]), as.numeric(ind_table[,1])), ncol = 10, byrow =T)
colnames(ind_table2) <- rownames(ind_table)
rownames(ind_table2) <- c("SD", "LD")
ind_table2 <- as.data.frame(ind_table2)
nn_ind_table <- non_negative_table(ind_table2)
circular_boxplot(res.table = nn_ind_table, "ld_sd_phases.png", tr.height = 0.14, color.palette = "Egypt")

nn_ind_table$median[2] - nn_ind_table$median[1]

 ##### Rhythmic genes in LD or SD (as determined by RAIN)
total.rhythmic.genes <- unique(c(rain.ld.rhythmic.genes,rain.sd.rhythmic.genes))
length(total.rhythmic.genes)

circacompare.ld.sd.total <- matrix(nrow=length(total.rhythmic.genes),ncol=15)
rownames(circacompare.ld.sd.total) <- total.rhythmic.genes

#Construct a table with the circadian data for all the genes rhythmic in LD, SD or both conditions
#Genes must not have an expression value of 0 at the six temporal points of one or both conditions 
#in order extract the circadian data using CircaCompare
for(i in 1:length(total.rhythmic.genes))
{
  gene.i <- total.rhythmic.genes[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"), c(paste(ld.zt,3,sep="_")))]
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"), c(paste(sd.zt,3,sep="_")))]
  ld.sd.expression.i <- cbind(ld.expression.i,sd.expression.i)
  ld.sd.vector <- c()
  for (j in 1:ncol(ld.sd.expression.i))
  {
    ld.sd.vector <- c(ld.sd.vector,ld.sd.expression.i[1,j])
  }
  
  if ((sum(ld.sd.vector[1:18] == 0) != 18) && (sum(ld.sd.vector[19:36] == 0)) != 18)
  {
    time.points <- seq(from=0,by=4,length.out = 18)
    
    ld.sd.df <- data.frame(time=c(time.points,time.points),
                           measure=ld.sd.vector,
                           group=c(rep("ld",18),rep("sd",18)))
    
    out.i <- circacompare(x = ld.sd.df, col_time = "time", col_group = "group",
                          col_outcome = "measure",alpha_threshold = 1 )
    circacompare.ld.sd.total[i,] <- out.i[[2]][,2]
  }
  else
  {
    circacompare.ld.sd.total[i,] <- NA
  }
}


colnames(circacompare.ld.sd.total) <- out.i[[2]][,1]
out.i[[1]]
circacompare.ld.sd.total
dim(circacompare.ld.sd.total)
circacompare.ld.sd.total <- circacompare.ld.sd.total[!is.na(circacompare.ld.sd.total[,1]),]
write.table(circacompare.ld.sd.total, file = "circacompare_complete.tsv", sep="\t",row.names = T)

#Amplitude
par(lwd=3)
boxplot(circacompare.ld.sd.total[,"ld amplitude estimate" ],
        circacompare.ld.sd.total[,"sd amplitude estimate" ],outline=F,
        col=c("royalblue3","tomato3"),names=c("LD","SD"),cex.axis=1.5, main="Distribution of amplitude")
sum(circacompare.ld.sd.total[,"Amplitude difference estimate" ] < 0)
sum(circacompare.ld.sd.total[,"Amplitude difference estimate" ] < 0) / length(total.rhythmic.genes)
sum(circacompare.ld.sd.total[,"P-value for amplitude difference" ] < 0.05)
sum(circacompare.ld.sd.total[,"P-value for amplitude difference" ] < 0.05)/length(total.rhythmic.genes)

#Mesor
boxplot(circacompare.ld.sd.total[,"ld mesor estimate" ],
        circacompare.ld.sd.total[,"sd mesor estimate" ],outline=F,
        col=c("royalblue3","tomato3"),names=c("LD","SD"),cex.axis=1.5, main="Distribution of mesor")
sum(circacompare.ld.sd.total[,"Mesor difference estimate" ] < 0)
sum(circacompare.ld.sd.total[,"Mesor difference estimate" ] < 0) / length(total.rhythmic.genes)
sum(circacompare.ld.sd.total[,"P-value for mesor difference" ] < 0.05)
sum(circacompare.ld.sd.total[,"P-value for mesor difference" ] < 0.05)/length(total.rhythmic.genes)

#Phase
boxplot(circacompare.ld.sd.total[,"ld peak time hours" ],circacompare.ld.sd.total[,"sd peak time hours" ],
        outline=F,col=c("royalblue3","tomato3"),names=c("LD","SD"),cex.axis=1.5, main="Distribution of phase")
median(circacompare.ld.sd.total[,"ld peak time hours" ])-median(circacompare.ld.sd.total[,"sd peak time hours" ])
sum(circacompare.ld.sd.total[,"Phase difference estimate" ] < 0) /length(total.rhythmic.genes)
sum(circacompare.ld.sd.total[,"Phase difference estimate" ] > 0) /length(total.rhythmic.genes)
sum(circacompare.ld.sd.total[,"P-value for difference in phase" ] < 0.05)/length(total.rhythmic.genes)
sum(circacompare.ld.sd.total[,"P-value for difference in phase" ] < 0.05)

#### LL-LD conditions
circacompare.ld.ll <- matrix(nrow=length(complete.ld.ll.rhythmic.genes),ncol=15)
rownames(circacompare.ld.ll) <- complete.ld.ll.rhythmic.genes

for(i in 1:length(complete.ld.ll.rhythmic.genes))
{
  gene.i <- complete.ld.ll.rhythmic.genes[i]
  ld.expression.i <- ld.ll.gene.expression[gene.i,c(paste(ld.zt,2,sep="_"))]
  ll.expression.i <- ld.ll.gene.expression[gene.i,c(paste(ld.zt,3,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 6)
  
  ld.ll.df <- data.frame(time=c(time.points,time.points),
                         measure=c(ld.expression.i, ll.expression.i),
                         group=c(rep("ld",6),rep("ll",6)))
  
  out.i <- circacompare(x = ld.ll.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  circacompare.ld.ll[i,] <- out.i[[2]][,2]
}

colnames(circacompare.ld.ll) <- out.i[[2]][,1]
write.table(circacompare.ld.ll, file = "circacompare_ld_ll.tsv", sep="\t",row.names = T)

#Amplitude
boxplot(circacompare.ld.ll[,"ld amplitude estimate" ],
        circacompare.ld.ll[,"ll amplitude estimate" ],outline=F,
        col=c("royalblue3","lightblue"),names=c("LD","LL"),cex.axis=1.5, main="Distribution of amplitude")
sum(circacompare.ld.ll[,"Amplitude difference estimate" ] < 0)
sum(circacompare.ld.ll[,"Amplitude difference estimate" ] < 0) / length(complete.ld.ll.rhythmic.genes)
sum(circacompare.ld.ll[,"P-value for amplitude difference" ] < 0.05)
sum(circacompare.ld.ll[,"P-value for amplitude difference" ] < 0.05)/length(complete.ld.ll.rhythmic.genes)

#Mesor
boxplot(circacompare.ld.ll[,"ld mesor estimate" ],
        circacompare.ld.ll[,"ll mesor estimate" ],outline=F,
        col=c("royalblue3","lightblue"),names=c("LD","LL"),cex.axis=1.5, main="Distribution of mesor")
sum(circacompare.ld.ll[,"Mesor difference estimate" ] < 0)
sum(circacompare.ld.ll[,"Mesor difference estimate" ] < 0) / length(complete.ld.ll.rhythmic.genes)
sum(circacompare.ld.ll[,"P-value for mesor difference"] < 0.05)
sum(circacompare.ld.ll[,"P-value for mesor difference" ] < 0.05)/length(complete.ld.ll.rhythmic.genes)

#Phase
boxplot(circacompare.ld.ll[,"ld peak time" ],circacompare.ld.ll[,"ll peak time" ],
        outline=F,col=c("royalblue3","lightblue"),names=c("LD","LL"),cex.axis=1.5, main="Distribution of phase")
median(circacompare.ld.ll[,"ld peak time" ])-median(circacompare.ld.ll[,"ll peak time" ])
sum(circacompare.ld.ll[,"Phase difference estimate" ] < 0) /length(complete.ld.ll.rhythmic.genes)
sum(circacompare.ld.ll[,"Phase difference estimate" ] > 0) /length(complete.ld.ll.rhythmic.genes)
sum(circacompare.ld.ll[,"P-value for difference in phase" ] < 0.05)/length(complete.ld.ll.rhythmic.genes)
sum(circacompare.ld.ll[,"P-value for difference in phase" ] < 0.05)


##### SD-LL
circacompare.sd.ll <- matrix(nrow=length(complete.sd.ll.rhythmic.genes),ncol=15)
rownames(circacompare.sd.ll) <- complete.sd.ll.rhythmic.genes

for(i in 1:length(complete.sd.ll.rhythmic.genes))
{
  gene.i <- complete.sd.ll.rhythmic.genes[i]
  
  sd.expression.i <- sd.ll.gene.expression[gene.i,c(paste(sd.zt,2,sep="_"))]
  ll.expression.i <- sd.ll.gene.expression[gene.i,c(paste(sd.zt,3,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 6)
  
  sd.ll.df <- data.frame(time=c(time.points,time.points),
                         measure=c(sd.expression.i, ll.expression.i),
                         group=c(rep("sd",6),rep("ll",6)))
  
  out.i <- circacompare(x = sd.ll.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  circacompare.sd.ll[i,] <- out.i[[2]][,2]
}

colnames(circacompare.sd.ll) <- out.i[[2]][,1]
write.table(circacompare.sd.ll, file = "circacompare_sd_ll.tsv", sep="\t",row.names = T)

#Amplitude
boxplot(circacompare.sd.ll[,"sd amplitude estimate" ],
        circacompare.sd.ll[,"ll amplitude estimate" ],outline=F,
        col=c("tomato3","lightpink"),names=c("SD","LL"),cex.axis=1.5, main="Distribution of amplitude")

sum(circacompare.sd.ll[,"Amplitude difference estimate" ] < 0)
sum(circacompare.sd.ll[,"Amplitude difference estimate" ] < 0) / length(complete.sd.ll.rhythmic.genes)
sum(circacompare.sd.ll[,"P-value for amplitude difference" ] < 0.05)
sum(circacompare.sd.ll[,"P-value for amplitude difference" ] < 0.05)/length(complete.sd.ll.rhythmic.genes)

#Mesor
boxplot(circacompare.sd.ll[,"sd mesor estimate" ],
        circacompare.sd.ll[,"ll mesor estimate" ],outline=F,
        col=c("tomato3","lightpink"),names=c("SD","LL"),cex.axis=1.5, main="Distribution of mesor")
sum(circacompare.sd.ll[,"Mesor difference estimate" ] < 0)
sum(circacompare.sd.ll[,"Mesor difference estimate" ] < 0) / length(complete.sd.ll.rhythmic.genes)
sum(circacompare.sd.ll[,"P-value for mesor difference" ] < 0.05)
sum(circacompare.sd.ll[,"P-value for mesor difference" ] < 0.05)/length(complete.sd.ll.rhythmic.genes)

#Phase
boxplot(circacompare.sd.ll[,"sd peak time" ],circacompare.sd.ll[,"ll peak time" ],
        outline=F,col=c("tomato3","lightpink"),names=c("SD","LL"),cex.axis=1.5, main="Distribution of phase")
median(circacompare.sd.ll[,"sd peak time" ])-median(circacompare.sd.ll[,"ll peak time" ])
sum(circacompare.sd.ll[,"Phase difference estimate" ] < 0) /length(complete.sd.ll.rhythmic.genes)
sum(circacompare.sd.ll[,"Phase difference estimate" ] > 0) /length(complete.sd.ll.rhythmic.genes)
sum(circacompare.sd.ll[,"P-value for difference in phase" ] < 0.05)/length(complete.sd.ll.rhythmic.genes)
sum(circacompare.sd.ll[,"P-value for difference in phase" ] < 0.05)


#Representative phase-changing genes between LD and SD conditions
prueba1.df <- data.frame(time=c(seq(from=0,by=4,length.out = 12),seq(from=0,by=4,length.out = 12)),
                         measure=c(genes.expression["Mp1g00120",c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"))],
                                   genes.expression["Mp1g00120",c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"))]),
                         group=c(rep("ld",12),rep("sd",12)))
prueba1 <- circacompare(x = prueba1.df, col_time = "time", col_group = "group", 
                        col_outcome = "measure",alpha_threshold = 0.05)
prueba1[[1]]

prueba2.df <- data.frame(time=c(seq(from=0,by=4,length.out = 12),seq(from=0,by=4,length.out = 12)),
                         measure=c(genes.expression["Mp1g01050",c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"))],
                                   genes.expression["Mp1g01050",c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"))]),
                         group=c(rep("ld",12),rep("sd",12)))
prueba2 <- circacompare(x = prueba2.df, col_time = "time", col_group = "group", 
                        col_outcome = "measure",alpha_threshold = 0.05)
prueba2[[1]]

#Representative non phase-changing genes between LD and SD conditions
prueba3.df <- data.frame(time=c(seq(from=0,by=4,length.out = 12),seq(from=0,by=4,length.out = 12)),
                         measure=c(genes.expression["Mp2g01260",c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"))],
                                   genes.expression["Mp2g01260",c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"))]),
                         group=c(rep("ld",12),rep("sd",12)))
prueba3 <- circacompare(x = prueba3.df, col_time = "time", col_group = "group", 
                        col_outcome = "measure",alpha_threshold = 0.05)
prueba3[[1]]


####Circadian, light-responsive and darkness-responsive in LD
library(circacompare)

circacompare.ld.ll.resp <- matrix(nrow=length(ld.responsive),ncol=15)
rownames(circacompare.ld.ll.resp) <- ld.responsive

for(i in 1:length(ld.responsive))
{
  gene.i <- ld.responsive[i]
  
  ld.expression.i.resp <- ld.ll.gene.expression[gene.i,c(paste(ld.zt,2,sep="_"))]
  ll.expression.i.resp <- ld.ll.gene.expression[gene.i,c(paste(ld.zt,3,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 6)
  
  ld.ll.df.resp <- data.frame(time=c(time.points,time.points),
                         measure=c(ld.expression.i.resp, ll.expression.i.resp),
                         group=c(rep("ld",6),rep("ll",6)))
  
  out.i <- circacompare(x = ld.ll.df.resp, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  circacompare.ld.ll.resp[i,] <- out.i[[2]][,2]
}
colnames(circacompare.ld.ll.resp) <- out.i[[2]][,1]
nrow(circacompare.ld.ll.resp)
ld.light.resp <- subset(circacompare.ld.ll.resp, circacompare.ld.ll.resp[,"Mesor difference estimate" ] > 0)
nrow(ld.light.resp)
ld.dark.resp <- subset(circacompare.ld.ll.resp, circacompare.ld.ll.resp[,"Mesor difference estimate" ] < 0)
nrow(ld.dark.resp)

data.ld.ll.resp <- matrix(c(length(complete.ld.ll.rhythmic.genes),nrow(ld.light.resp),nrow(ld.dark.resp),length(non.rhythmic.ld.genes)),nrow=4)

barplot(data.ld.ll.resp, lwd=3, names.arg = "",cex.names = 2,cex.axis = 1.5,
        col=c("blue","purple","lightblue","white"), 
        border="black", width = 0.2, xlim = c(0,1),
        font.axis=2, main = "LD-LL rhythmic genes", legend.text = c("Circadian", "Inducted by light", "Repressed by light", "Non rhythmic"),
        xlab="",ylim=c(0,25000))

non.expressed.genes.ld <- which(apply(X = genes.expression[, c(ld.zt.i)],MARGIN = 1,FUN = max) < 1)
length(non.expressed.genes.ld)
non.expressed.genes.ld.names <- names(non.expressed.genes.ld)

data.ld.ll.resp <- matrix(c(length(complete.ld.ll.rhythmic.genes),nrow(ld.light.resp),nrow(ld.dark.resp),
                            length(setdiff(non.rhythmic.ld.genes,non.expressed.genes.ld.names))),nrow=4)
par(lwd=3)
barplot(data.ld.ll.resp, lwd=3, names.arg = "",cex.names = 2,cex.axis = 1.5,
        col=c("blue","purple","lightblue","white"), 
        border="black", width = 0.2, xlim = c(0,1),
        font.axis=2, main = "LD rhythmic genes", legend.text = c("Circadian", "Inducted by light", "Repressed by light", "Non rhythmic"),
        xlab="",ylim=c(0,25000))


write(non.expressed.genes.ld.names, file = "gene_lists/no_exp_ld_genes.txt")
write(setdiff(non.rhythmic.ld.genes,non.expressed.genes.ld.names), file = "gene_lists/exp_but_no_rhythmic_ld_genes.txt")
write(rownames(ld.light.resp), file = "gene_lists/light_ind_ld_genes.txt")
write(rownames(ld.dark.resp), file = "gene_lists/light_rep_ld_genes.txt")
write(complete.ld.ll.rhythmic.genes, file = "gene_lists/circadian_ld_genes.txt")

####Circadian, light-responsive and darkness-responsive in SD
circacompare.sd.ll.resp <- matrix(nrow=length(sd.responsive),ncol=15)
rownames(circacompare.sd.ll.resp) <- sd.responsive

for(i in 1:length(sd.responsive))
{
  gene.i <- sd.responsive[i]
  
  sd.expression.i <- sd.ll.gene.expression[gene.i,c(paste(sd.zt,2,sep="_"))]
  ll.expression.i <- sd.ll.gene.expression[gene.i,c(paste(sd.zt,3,sep="_"))]
  sd.ll.vector <- c(sd.expression.i,ll.expression.i)
  
  if ((sum(sd.ll.vector[1:6] == 0) != 6) && (sum(sd.ll.vector[7:12] == 0)) != 6)
  {
    time.points <- seq(from=0,by=4,length.out = 6)
    
    sd.ll.df <- data.frame(time=c(time.points,time.points),
                           measure=sd.ll.vector,
                           group=c(rep("sd",6),rep("ll",6)))
    
    out.i <- circacompare(x = sd.ll.df, col_time = "time", col_group = "group",
                          col_outcome = "measure",alpha_threshold = 1 )
    circacompare.sd.ll.resp[i,] <- out.i[[2]][,2]
  }
  else
  {
    circacompare.sd.ll.resp[i,] <- NA
  }
}
colnames(circacompare.sd.ll.resp) <- out.i[[2]][,1]
nrow(circacompare.sd.ll.resp)
# <0 for response to light because the mesor difference in this case corresponds to sd-ll
sd.light.resp <- subset(circacompare.sd.ll.resp, circacompare.sd.ll.resp[,"Mesor difference estimate" ] < 0) 
nrow(sd.light.resp)
sd.dark.resp <- subset(circacompare.sd.ll.resp, circacompare.sd.ll.resp[,"Mesor difference estimate" ] > 0)
nrow(sd.dark.resp)
sd.na.resp <- subset(circacompare.sd.ll.resp, is.na(circacompare.sd.ll.resp[,"Mesor difference estimate" ]))
nrow(sd.na.resp)

data.sd.ll.resp <- matrix(c(length(complete.sd.ll.rhythmic.genes),nrow(sd.light.resp),nrow(sd.dark.resp),length(non.rhythmic.sd.genes)),nrow=4)

par(lwd=3)
barplot(data.sd.ll.resp, lwd=3, names.arg = "",cex.names = 2,cex.axis = 1.5,
        col=c("red","tomato","orange","white"), 
        border="black", width = 0.2, xlim = c(0,1),
        font.axis=2, main = "SD-LL rhythmic genes", legend.text = c("Circadian", "Inducted by light", "Repressed by light", "Non rhythmic"),
        xlab="",ylim=c(0,25000))

non.expressed.genes.sd <- which(apply(X = genes.expression[, c(sd.zt.i)],MARGIN = 1,FUN = max) < 1)
length(non.expressed.genes.sd)
non.expressed.genes.sd.names <- names(non.expressed.genes.sd)

data.sd.ll.resp <- matrix(c(length(complete.sd.ll.rhythmic.genes),nrow(sd.light.resp),nrow(sd.dark.resp),
                            length(setdiff(non.rhythmic.sd.genes,non.expressed.genes.sd.names))),nrow=4)
par(lwd=3)
barplot(data.sd.ll.resp, lwd=3, names.arg = "",cex.names = 2,cex.axis = 1.5,
        col=c("red","tomato","orange","white"), 
        border="black", width = 0.2, xlim = c(0,1),
        font.axis=2, main = "SD rhythmic genes", legend.text = c("Circadian", "Inducted by light", "Repressed by light", "Non rhythmic"),
        xlab="",ylim=c(0,25000))

write(non.expressed.genes.sd.names, file = "gene_lists/no_exp_sd_genes.txt")
write(setdiff(non.rhythmic.sd.genes,non.expressed.genes.sd.names), file = "gene_lists/exp_but_no_rhythmic_sd_genes.txt")
write(rownames(sd.light.resp), file = "gene_lists/light_ind_sd_genes.txt")
write(rownames(sd.dark.resp), file = "gene_lists/light_rep_sd_genes.txt")
write(complete.sd.ll.rhythmic.genes, file = "gene_lists/circadian_sd_genes.txt")

#### Clustering of rhythmic genes by the ZT at which they show their expression peak

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

#Clusters for rhythmic genes in LD with a 24h period
clusters(rain.results = results.ld,gene.set = rhythmic.genes.ld,photoperiod = "LD")

#Clusters for rhythmic genes in SD with a 24h period
clusters(rain.results = results.sd,gene.set = rhythmic.genes.sd,photoperiod = "SD")

#Clusters for rhythmic genes in LD with a 12h period
clusters(rain.results = results.ld.12,gene.set = rhythmic.genes.ld.12,photoperiod = "LD_double")

#Clusters for rhythmic genes in SD with a 12h period
clusters(rain.results = results.sd.12,gene.set = rhythmic.genes.sd.12,photoperiod = "SD_double")

## Clusters plots

#24h period rhythmic genes
cluster.LD.ZT00<-read.table("LDpeak_clusters/cluster_ZT0.txt")
vector.LD.ZT00<-as.vector(cluster.LD.ZT00$V1)
cluster.LD.ZT04<-read.table("LDpeak_clusters/cluster_ZT4.txt")
vector.LD.ZT04<-as.vector(cluster.LD.ZT04$V1)
cluster.LD.ZT08<-read.table("LDpeak_clusters/cluster_ZT8.txt")
vector.LD.ZT08<-as.vector(cluster.LD.ZT08$V1)
cluster.LD.ZT12<-read.table("LDpeak_clusters/cluster_ZT12.txt")
vector.LD.ZT12<-as.vector(cluster.LD.ZT12$V1)
cluster.LD.ZT16<-read.table("LDpeak_clusters/cluster_ZT16.txt")
vector.LD.ZT16<-as.vector(cluster.LD.ZT16$V1)
cluster.LD.ZT20<-read.table("LDpeak_clusters/cluster_ZT20.txt")
vector.LD.ZT20<-as.vector(cluster.LD.ZT20$V1)
cluster.SD.ZT00<-read.table("SDpeak_clusters/cluster_ZT0.txt")
vector.SD.ZT00<-as.vector(cluster.SD.ZT00$V1)
cluster.SD.ZT04<-read.table("SDpeak_clusters/cluster_ZT4.txt")
vector.SD.ZT04<-as.vector(cluster.SD.ZT04$V1)
cluster.SD.ZT08<-read.table("SDpeak_clusters/cluster_ZT8.txt")
vector.SD.ZT08<-as.vector(cluster.SD.ZT08$V1)
cluster.SD.ZT12<-read.table("SDpeak_clusters/cluster_ZT12.txt")
vector.SD.ZT12<-as.vector(cluster.SD.ZT12$V1)
cluster.SD.ZT16<-read.table("SDpeak_clusters/cluster_ZT16.txt")
vector.SD.ZT16<-as.vector(cluster.SD.ZT16$V1)
cluster.SD.ZT20<-read.table("SDpeak_clusters/cluster_ZT20.txt")
vector.SD.ZT20<-as.vector(cluster.SD.ZT20$V1)

clusters.bars <- matrix(c(length(vector.LD.ZT00),length(vector.LD.ZT04),length(vector.LD.ZT08),
                          length(vector.LD.ZT12),length(vector.LD.ZT16),length(vector.LD.ZT20), 
                          length(vector.SD.ZT00),length(vector.SD.ZT04),length(vector.SD.ZT08),
                          length(vector.SD.ZT12),length(vector.SD.ZT16),length(vector.SD.ZT20)),nrow=6)

barplot(clusters.bars, lwd=2, names.arg = c("LD","SD"),cex.names = 2,cex.axis = 1.5, width = 0.35, xlim = c(0,2),
        col=c("gray","darkolivegreen1","darkolivegreen2","darkolivegreen3","darkolivegreen4","darkolivegreen"),
        border="black", space=0.25, font.axis=2, xlab="", ylim=c(0,7000), legend.text = c("ZT00","ZT04","ZT08","ZT12","ZT16","ZT20"),
        main = "Number of genes per cluster")

#12h period rhythmic genes
cluster.LD.double.ZT00<-read.table("LD_doublepeak_clusters/cluster_ZT0.txt")
vector.LD.double.ZT00<-as.vector(cluster.LD.double.ZT00$V1)
cluster.LD.double.ZT04<-read.table("LD_doublepeak_clusters/cluster_ZT4.txt")
vector.LD.double.ZT04<-as.vector(cluster.LD.double.ZT04$V1)
cluster.LD.double.ZT08<-read.table("LD_doublepeak_clusters/cluster_ZT8.txt")
vector.LD.double.ZT08<-as.vector(cluster.LD.double.ZT08$V1)
cluster.SD.double.ZT00<-read.table("SD_doublepeak_clusters/cluster_ZT0.txt")
vector.SD.double.ZT00<-as.vector(cluster.SD.double.ZT00$V1)
cluster.SD.double.ZT04<-read.table("SD_doublepeak_clusters/cluster_ZT4.txt")
vector.SD.double.ZT04<-as.vector(cluster.SD.double.ZT04$V1)
cluster.SD.double.ZT08<-read.table("SD_doublepeak_clusters/cluster_ZT8.txt")
vector.SD.double.ZT08<-as.vector(cluster.SD.double.ZT08$V1)

clusters.double.bars <- matrix(c(length(vector.LD.double.ZT00),length(vector.LD.double.ZT04),
                                 length(vector.LD.double.ZT08),length(vector.SD.double.ZT00),
                                 length(vector.SD.double.ZT04),length(vector.SD.double.ZT08)),nrow=3)

barplot(clusters.double.bars, lwd=2, names.arg = c("LD","SD"),cex.names = 2,cex.axis = 1.5, width = 0.35, xlim = c(0,2),
        col=c("darkolivegreen1","darkolivegreen3","darkolivegreen"),
        border="black", space=0.25, font.axis=2, xlab="", ylim=c(0,3000), legend.text = c("ZT00","ZT04","ZT08"),
        main = "Number of genes per cluster")

##Circacompare histograms##
head(circacompare.ld.sd)
head(circacompare.ld.sd.total)

length(rain.ld.rhythmic.genes)
length(intersect(rain.ld.rhythmic.genes, rownames(circacompare.ld.sd.total)))
rhythmic.for.circa.ld <- intersect(rain.ld.rhythmic.genes, rownames(circacompare.ld.sd.total))

length(rain.sd.rhythmic.genes)
length(intersect(rain.sd.rhythmic.genes, rownames(circacompare.ld.sd.total)))
rhythmic.for.circa.sd <- intersect(rain.sd.rhythmic.genes, rownames(circacompare.ld.sd.total))

circacompare.ld.rhythmic.peaks <- circacompare.ld.sd.total[rhythmic.for.circa.ld,12]
head(circacompare.ld.rhythmic.peaks)

circacompare.sd.rhythmic.peaks <- circacompare.ld.sd.total[rhythmic.for.circa.sd,13]
head(circacompare.sd.rhythmic.peaks)

par(lwd=2)
hist(circacompare.ld.rhythmic.peaks, freq = T, col = "blue", xlim = c(0,24), lwd= 2, xaxt='n',
     main= "Pico de expresión de los genes rítmicos en LD", xlab= "Tiempo Zeitberg", ylab= "Número de genes", yaxt='n')
axis(side=1, at=seq(0,24,4), seq(0,24,4), lwd=2)
axis(side=2, at=seq(0,900,300), seq(0,900,300), lwd=2)
abline(v=16, col="black", lwd=5, lty = 2)

par(lwd=2)
hist(circacompare.sd.rhythmic.peaks, freq = T, col = "red", xlim = c(0,24), lwd= 2, xaxt='n',
     main= "Pico de expresión de los genes rítmicos en SD", xlab= "Tiempo Zeitberg", yaxt='n', ylab= "Número de genes")
axis(side=1, at=seq(0,24,4), seq(0,24,4), lwd=2)
axis(side=2, at=seq(0,1200,300), seq(0,1200,300), lwd=2)
abline(v=8, col="black", lwd=5, lty = 2)



#### FUNCTIONAL ANALYSIS ####

#### GO terms enrichment 
library("clusterProfiler")

#Install our annotation package for Marchantia polymorpha
#install.packages("./org.Mpolymorpha.eg.db/", repos=NULL, type="source")
library("org.Mpolymorpha.eg.db")
columns(org.Mpolymorpha.eg.db)


#mp.universe <- unique(select(org.Mpolymorpha.eg.db,columns = c("GO"),keys=keys(org.Mpolymorpha.eg.db,keytype = "GID"))[["GID"]])
#length(mp.universe)

#SD
#Cluster ZT00 SD
cluster.ZT00.SD<-read.table("SDpeak_clusters/cluster_ZT0.txt")
cluster.ZT00.SD<-as.vector(cluster.ZT00.SD$V1)

ego.ZT00.SD <- enrichGO(gene = unique(cluster.ZT00.SD),
                        OrgDb         = org.Mpolymorpha.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        keyType = "GID")

barplot(ego.ZT00.SD, title = "SD ZT00")

#Cluster ZT04 SD
cluster.ZT04.SD<-read.table("SDpeak_clusters/cluster_ZT4.txt")
cluster.ZT04.SD<-as.vector(cluster.ZT04.SD$V1)

ego.ZT04.SD <- enrichGO(gene = unique(cluster.ZT04.SD),
                        OrgDb         = org.Mpolymorpha.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        keyType = "GID")

barplot(ego.ZT04.SD, title = "SD ZT04")

#Cluster ZT08 SD
cluster.ZT08.SD<-read.table("SDpeak_clusters/cluster_ZT8.txt")
cluster.ZT08.SD<-as.vector(cluster.ZT08.SD$V1)

ego.ZT08.SD <- enrichGO(gene = unique(cluster.ZT08.SD),
                OrgDb         = org.Mpolymorpha.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE,
                keyType = "GID")

barplot(ego.ZT08.SD, title = "SD ZT08")

#Cluster ZT12 SD
cluster.ZT12.SD<-read.table("SDpeak_clusters/cluster_ZT12.txt")
cluster.ZT12.SD<-as.vector(cluster.ZT12.SD$V1)

ego.ZT12.SD <- enrichGO(gene = unique(cluster.ZT12.SD),
                        OrgDb         = org.Mpolymorpha.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        keyType = "GID")

barplot(ego.ZT12.SD, title = "SD ZT12")

#Cluster ZT16 SD
cluster.ZT16.SD<-read.table("SDpeak_clusters/cluster_ZT16.txt")
cluster.ZT16.SD<-as.vector(cluster.ZT16.SD$V1)

ego.ZT16.SD <- enrichGO(gene = unique(cluster.ZT16.SD),
                        OrgDb         = org.Mpolymorpha.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        keyType = "GID")

barplot(ego.ZT16.SD, title = "SD ZT16")

#Cluster ZT20 SD
cluster.ZT20.SD<-read.table("SDpeak_clusters/cluster_ZT20.txt")
cluster.ZT20.SD<-as.vector(cluster.ZT20.SD$V1)

ego.ZT20.SD <- enrichGO(gene = unique(cluster.ZT20.SD),
                        OrgDb         = org.Mpolymorpha.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        keyType = "GID")

barplot(ego.ZT20.SD, title = "SD ZT20")

#LD
#Cluster ZT00 LD
cluster.ZT00.LD<-read.table("LDpeak_clusters/cluster_ZT0.txt")
cluster.ZT00.LD<-as.vector(cluster.ZT00.LD$V1)

ego.ZT00.LD <- enrichGO(gene = unique(cluster.ZT00.LD),
                        OrgDb         = org.Mpolymorpha.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        keyType = "GID")

barplot(ego.ZT00.LD, title = "LD ZT00")

#Cluster ZT04 LD
cluster.ZT04.LD<-read.table("LDpeak_clusters/cluster_ZT4.txt")
cluster.ZT04.LD<-as.vector(cluster.ZT04.LD$V1)

ego.ZT04.LD <- enrichGO(gene = unique(cluster.ZT04.LD),
                        OrgDb         = org.Mpolymorpha.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        keyType = "GID")

barplot(ego.ZT04.LD, title = "LD ZT04")

#Cluster ZT08 LD
cluster.ZT08.LD<-read.table("LDpeak_clusters/cluster_ZT8.txt")
cluster.ZT08.LD<-as.vector(cluster.ZT08.LD$V1)

ego.ZT08.LD <- enrichGO(gene = unique(cluster.ZT08.LD),
                        OrgDb         = org.Mpolymorpha.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        keyType = "GID")

barplot(ego.ZT08.LD, title = "LD ZT08")

#Cluster ZT12 LD
cluster.ZT12.LD<-read.table("LDpeak_clusters/cluster_ZT12.txt")
cluster.ZT12.LD<-as.vector(cluster.ZT12.LD$V1)

ego.ZT12.LD <- enrichGO(gene = unique(cluster.ZT12.LD),
                        OrgDb         = org.Mpolymorpha.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        keyType = "GID")

barplot(ego.ZT12.LD,title = "LD ZT12")

#Cluster ZT16 LD
cluster.ZT16.LD<-read.table("LDpeak_clusters/cluster_ZT16.txt")
cluster.ZT16.LD<-as.vector(cluster.ZT16.LD$V1)

ego.ZT16.LD <- enrichGO(gene = unique(cluster.ZT16.LD),
                        OrgDb         = org.Mpolymorpha.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        keyType = "GID")

barplot(ego.ZT16.LD, title = "LD ZT16")

#Cluster ZT20 LD
cluster.ZT20.LD<-read.table("LDpeak_clusters/cluster_ZT20.txt")
cluster.ZT20.LD<-as.vector(cluster.ZT20.LD$V1)

ego.ZT20.LD <- enrichGO(gene = unique(cluster.ZT20.LD),
                        OrgDb         = org.Mpolymorpha.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        keyType = "GID")

barplot(ego.ZT20.LD, title = "LD ZT20")


#### GO terms enrichment based in orthologs of Arabidopsis thaliana
library(org.Mpolymorpha.eg.db)

ath_orth <- as.vector(select(org.Mpolymorpha.eg.db,columns = c("ATH_ORTH"),keys=keys(org.Mpolymorpha.eg.db,keytype = "GID"))$ATH_ORTH)
names(ath_orth) <- as.vector(select(org.Mpolymorpha.eg.db,columns = c("ATH_ORTH"),keys=keys(org.Mpolymorpha.eg.db,keytype = "GID"))$GID)
head(ath_orth)

#LD
#Cluster ZT00 LD

cluster.ZT00.LD<-read.table("LDpeak_clusters/cluster_ZT0.txt")
cluster.ZT00.LD<-as.vector(cluster.ZT00.LD$V1)
head (cluster.ZT00.LD)
ZT00.arabi <-ath_orth[cluster.ZT00.LD]
ZT00.arabi <- ZT00.arabi[!(is.na(ZT00.arabi))]
names(ZT00.arabi) <- NULL
ZT00.arabi
ZT00.arabi <- ZT00.arabi[ZT00.arabi != ""]
ZT00.arabi

library("org.At.tair.db")
ego.ZT00.LD.ath <- enrichGO(gene = unique(ZT00.arabi),
                        OrgDb         = org.At.tair.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = FALSE,
                        keyType = "TAIR")
barplot(ego.ZT00.LD.ath, showCategory = 30)

#Cluster ZT04 LD
ZT04.arabi <-ath_orth[cluster.ZT04.LD]
ZT04.arabi <- ZT04.arabi[!(is.na(ZT04.arabi))]
names(ZT04.arabi) <- NULL
ZT04.arabi
ZT04.arabi <- ZT04.arabi[ZT04.arabi != ""]
ZT04.arabi
ego.ZT04.LD.ath <- enrichGO(gene = unique(ZT04.arabi),
                            OrgDb         = org.At.tair.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = FALSE,
                            keyType = "TAIR")
barplot(ego.ZT04.LD.ath)

#Cluster ZT08 LD
ZT08.arabi <-ath_orth[cluster.ZT08.LD]
ZT08.arabi <- ZT08.arabi[!(is.na(ZT08.arabi))]
names(ZT08.arabi) <- NULL
ZT08.arabi
ZT08.arabi <- ZT08.arabi[ZT08.arabi != ""]
ZT08.arabi
ego.ZT08.LD.ath <- enrichGO(gene = unique(ZT08.arabi),
                            OrgDb         = org.At.tair.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = FALSE,
                            keyType = "TAIR")
barplot(ego.ZT08.LD.ath)

#Cluster ZT12 LD
ZT12.arabi <-ath_orth[cluster.ZT12.LD]
ZT12.arabi <- ZT12.arabi[!(is.na(ZT12.arabi))]
names(ZT12.arabi) <- NULL
ZT12.arabi
ZT12.arabi <- ZT12.arabi[ZT12.arabi != ""]
ZT12.arabi
ego.ZT12.LD.ath <- enrichGO(gene = unique(ZT12.arabi),
                            OrgDb         = org.At.tair.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = FALSE,
                            keyType = "TAIR")
barplot(ego.ZT12.LD.ath)

#Cluster ZT16 LD
ZT16.arabi <-ath_orth[cluster.ZT16.LD]
ZT16.arabi <- ZT16.arabi[!(is.na(ZT16.arabi))]
names(ZT16.arabi) <- NULL
ZT16.arabi
ZT16.arabi <- ZT16.arabi[ZT16.arabi != ""]
ZT16.arabi
ego.ZT16.LD.ath <- enrichGO(gene = unique(ZT16.arabi),
                            OrgDb         = org.At.tair.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = FALSE,
                            keyType = "TAIR")
barplot(ego.ZT16.LD.ath)

#Cluster ZT20 LD
ZT20.arabi <-ath_orth[cluster.ZT20.LD]
ZT20.arabi <- ZT20.arabi[!(is.na(ZT20.arabi))]
names(ZT20.arabi) <- NULL
ZT20.arabi
ZT20.arabi <- ZT20.arabi[ZT20.arabi != ""]
ZT20.arabi
ego.ZT20.LD.ath <- enrichGO(gene = unique(ZT20.arabi),
                            OrgDb         = org.At.tair.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = FALSE,
                            keyType = "TAIR")
barplot(ego.ZT20.LD.ath)

#SD
#Cluster ZT00 SD
ZT00.arabi.SD <-ath_orth[cluster.ZT00.SD]
ZT00.arabi.SD <- ZT00.arabi.SD[!(is.na(ZT00.arabi.SD))]
names(ZT00.arabi.SD) <- NULL
ZT00.arabi.SD
ZT00.arabi.SD <- ZT00.arabi.SD[ZT00.arabi.SD != ""]
ZT00.arabi.SD
ego.ZT00.SD.ath <- enrichGO(gene = unique(ZT00.arabi.SD),
                            OrgDb         = org.At.tair.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = FALSE,
                            keyType = "TAIR")
barplot(ego.ZT00.SD.ath)

#Cluster ZT04 SD
ZT04.arabi.SD <-ath_orth[cluster.ZT04.SD]
ZT04.arabi.SD <- ZT04.arabi.SD[!(is.na(ZT04.arabi.SD))]
names(ZT04.arabi.SD) <- NULL
ZT04.arabi.SD
ZT04.arabi.SD <- ZT04.arabi.SD[ZT04.arabi.SD != ""]
ZT04.arabi.SD
ego.ZT04.SD.ath <- enrichGO(gene = unique(ZT04.arabi.SD),
                            OrgDb         = org.At.tair.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = FALSE,
                            keyType = "TAIR")
barplot(ego.ZT04.SD.ath)

#Cluster ZT08 SD
ZT08.arabi.SD <-ath_orth[cluster.ZT08.SD]
ZT08.arabi.SD <- ZT08.arabi.SD[!(is.na(ZT08.arabi.SD))]
names(ZT08.arabi.SD) <- NULL
ZT08.arabi.SD
ZT08.arabi.SD <- ZT08.arabi.SD[ZT08.arabi.SD != ""]
ZT08.arabi.SD
ego.ZT08.SD.ath <- enrichGO(gene = unique(ZT08.arabi.SD),
                            OrgDb         = org.At.tair.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = FALSE,
                            keyType = "TAIR")
barplot(ego.ZT08.SD.ath)

#Cluster ZT12 SD
ZT12.arabi.SD <-ath_orth[cluster.ZT12.SD]
ZT12.arabi.SD <- ZT12.arabi.SD[!(is.na(ZT12.arabi.SD))]
names(ZT12.arabi.SD) <- NULL
ZT12.arabi.SD
ZT12.arabi.SD <- ZT12.arabi.SD[ZT12.arabi.SD != ""]
ZT12.arabi.SD
ego.ZT12.SD.ath <- enrichGO(gene = unique(ZT12.arabi.SD),
                            OrgDb         = org.At.tair.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = FALSE,
                            keyType = "TAIR")
barplot(ego.ZT12.SD.ath)

#Cluster ZT16 SD
ZT16.arabi.SD <-ath_orth[cluster.ZT16.SD]
ZT16.arabi.SD <- ZT16.arabi.SD[!(is.na(ZT16.arabi.SD))]
names(ZT16.arabi.SD) <- NULL
ZT16.arabi.SD
ZT16.arabi.SD <- ZT16.arabi.SD[ZT16.arabi.SD != ""]
ZT16.arabi.SD
ego.ZT16.SD.ath <- enrichGO(gene = unique(ZT16.arabi.SD),
                            OrgDb         = org.At.tair.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = FALSE,
                            keyType = "TAIR")
barplot(ego.ZT16.SD.ath)

#Cluster ZT20 SD
ZT20.arabi.SD <-ath_orth[cluster.ZT20.SD]
ZT20.arabi.SD <- ZT20.arabi.SD[!(is.na(ZT20.arabi.SD))]
names(ZT20.arabi.SD) <- NULL
ZT20.arabi.SD
ZT20.arabi.SD <- ZT20.arabi.SD[ZT20.arabi.SD != ""]
ZT20.arabi.SD
ego.ZT20.SD.ath <- enrichGO(gene = unique(ZT20.arabi.SD),
                            OrgDb         = org.At.tair.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = FALSE,
                            keyType = "TAIR")
barplot(ego.ZT20.SD.ath)

#Save GO terms enrichment for M.poly genes as data frames
data.go.sd.00 <- as.data.frame(ego.ZT00.SD)
write.table(data.go.sd.00,file = "data.go.sd.00.tsv")

data.go.sd.04 <- as.data.frame(ego.ZT04.SD)
write.table(data.go.sd.04,file = "data.go.sd.04.tsv")

data.go.sd.08 <- as.data.frame(ego.ZT08.SD)
write.table(data.go.sd.08,file = "data.go.sd.08.tsv")

data.go.sd.12 <- as.data.frame(ego.ZT12.SD)
write.table(data.go.sd.12,file = "data.go.sd.12.tsv")

data.go.sd.16 <- as.data.frame(ego.ZT16.SD)
write.table(data.go.sd.16,file = "data.go.sd.16.tsv")

data.go.sd.20 <- as.data.frame(ego.ZT20.SD)
write.table(data.go.sd.20,file = "data.go.sd.20.tsv")

data.go.ld.00 <- as.data.frame(ego.ZT00.LD)
write.table(data.go.ld.00,file = "data.go.ld.00.tsv")

data.go.ld.04 <- as.data.frame(ego.ZT04.LD)
write.table(data.go.ld.04,file = "data.go.ld.04.tsv")

data.go.ld.08 <- as.data.frame(ego.ZT08.LD)
write.table(data.go.ld.08,file = "data.go.ld.08.tsv")

data.go.ld.12 <- as.data.frame(ego.ZT12.LD)
write.table(data.go.ld.12,file = "data.go.ld.12.tsv")

data.go.ld.16 <- as.data.frame(ego.ZT16.LD)
write.table(data.go.ld.16,file = "data.go.ld.16.tsv")

data.go.ld.20 <- as.data.frame(ego.ZT20.LD)
write.table(data.go.ld.20,file = "data.go.ld.20.tsv")


#Save GO terms enrichment for A.thaliana orthologs as data frames
data.go.ld.00.ath <- as.data.frame(ego.ZT00.LD.ath)
write.table(data.go.ld.00.ath,file = "data.go.ld.00.ath.tsv")

data.go.ld.04.ath <- as.data.frame(ego.ZT04.LD.ath)
write.table(data.go.ld.04.ath,file = "data.go.ld.04.ath.tsv")

data.go.ld.08.ath <- as.data.frame(ego.ZT08.LD.ath)
write.table(data.go.ld.08.ath,file = "data.go.ld.08.ath.tsv")

data.go.ld.12.ath <- as.data.frame(ego.ZT12.LD.ath)
write.table(data.go.ld.12.ath,file = "data.go.ld.12.ath.tsv")

data.go.ld.16.ath <- as.data.frame(ego.ZT16.LD.ath)
write.table(data.go.ld.16.ath,file = "data.go.ld.16.ath.tsv")

data.go.ld.20.ath <- as.data.frame(ego.ZT20.LD.ath)
write.table(data.go.ld.20.ath,file = "data.go.ld.20.ath.tsv")

data.go.sd.00.ath <- as.data.frame(ego.ZT00.SD.ath)
write.table(data.go.sd.00.ath,file = "data.go.sd.00.ath.tsv")

data.go.sd.04.ath <- as.data.frame(ego.ZT04.SD.ath)
write.table(data.go.sd.04.ath,file = "data.go.sd.04.ath.tsv")

data.go.sd.08.ath <- as.data.frame(ego.ZT08.SD.ath)
write.table(data.go.sd.08.ath,file = "data.go.sd.08.ath.tsv")

data.go.sd.12.ath <- as.data.frame(ego.ZT12.SD.ath)
write.table(data.go.sd.12.ath,file = "data.go.sd.12.ath.tsv")

data.go.sd.16.ath <- as.data.frame(ego.ZT16.SD.ath)
write.table(data.go.sd.16.ath,file = "data.go.sd.16.ath.tsv")

data.go.sd.20.ath <- as.data.frame(ego.ZT20.SD.ath)
write.table(data.go.sd.20.ath,file = "data.go.sd.20.ath.tsv")


#### KEGG pathways enrichment

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("pathview")

library(clusterProfiler)
library("pathview")

## KO enrichment using M.poly annotation
annotation.table <- read.table("Marchantia_polymorpha_v5.1_final_annotation.tsv", header = T, as.is = T)
head(annotation.table)
dim(annotation.table)
ko.table <- annotation.table[,c("ID5","KO")]
head(ko.table)
dim(ko.table)
ko.final.table <- subset(ko.table, ko.table$KO != "")
dim(ko.final.table)
head(ko.final.table)

ko.universe <- ko.final.table$KO

for (photoperiod in c("LD","SD"))
{
  for (i in seq(from=0,to=20,by=4))
  {
    target.cluster <- paste0(photoperiod,"peak_clusters/cluster_ZT",i,".txt")
    target.names.table <- read.table(target.cluster)
    target.ids <- target.names.table[,1]
    target.ko <- subset(ko.final.table, ko.final.table$ID5 %in% target.ids)$KO
    target.ko <- target.ko[!is.na(target.ko)]
    pathway.enrich <- as.data.frame(enrichKEGG(gene=target.ko, organism="ko", universe=ko.universe, qvalueCutoff=0.05))
    write.table(file = paste0("KO_enrichment/KO_ZT",i,"_",photoperiod,".tsv"),pathway.enrich)
  }
}

for (photoperiod in c("LD","SD"))
{
  for (i in seq(from=0,to=20,by=4))
  {
    target.cluster <- paste0(photoperiod,"peak_clusters/cluster_ZT",i,".txt")
    target.names.table <- read.table(target.cluster)
    target.ids <- target.names.table[,1]
    target.ko <- subset(ko.final.table, ko.final.table$ID5 %in% target.ids)$KO
    gene.pathway <- rep(0, length(ko.universe))
    names(gene.pathway) <-  ko.universe
    gene.pathway[target.ko] <-1
    ko.cluster.name <- paste0("KO_enrichment/KO_ZT",i,"_",photoperiod,".tsv")
    ko.cluster <- read.table(ko.cluster.name,header = T)
    for (j in 1:nrow(ko.cluster))
      {
        pathway.current.id <- ko.cluster[j,1]
        pathview(gene.data = sort(gene.pathway, decreasing = T), pathway.id = pathway.current.id,species = "ko",
                 limit=list(gene=max(abs(gene.pathway)),cpd=1),gene.idtype = "kegg",
                 kegg.dir = paste0("kegg_images/KO_ZT",i,"_",photoperiod))
      }
  }
}  


#Plots of genes involved in pathways of interest in each temporal point and photoperiod
goi.mpoly <- read.table("../pathways_gene_correspondece.tsv",header=T, sep="\t")
head(goi.mpoly)
goi.mpoly.sd.ids<-goi.mpoly[,3]
goi.mpoly.sd.ids

for (i in 1:length(goi.mpoly.sd.ids))
{
  genes.vector <- strsplit(goi.mpoly.sd.ids[i],split = "/")[[1]]
  if (i %in% 1:19)
  {
    pathway.photoperiod <- "LD"
  }
  else
  {
    pathway.photoperiod <- "SD"
  }
  if (i %in% c(1,2,3,20,21,22))
  {
    pathway.point <- "ZT00" 
  }
  else if (i %in% c(4,5,6,7,23,24,25,26))
  {
    pathway.point <- "ZT04" 
  }
  else if (i %in% c(8,9,10,27,28,29))
  {
    pathway.point <- "ZT08" 
  }
  else if (i %in% c(11,12,13,30,31,32))
  {
    pathway.point <- "ZT12" 
  }
  else if (i %in% c(14,15,16,33,34,35))
  {
    pathway.point <- "ZT16" 
  }
  else if (i %in% c(17,18,19,36,37,38))
  {
    pathway.point <- "ZT20" 
  }
  gene.temporal.point <- paste(pathway.photoperiod, pathway.point, sep = "_")
  for (j in 1:length(genes.vector))
  {
    plot.ld.sd(gene.id = genes.vector[j],gene.name = gene.temporal.point, 
               ld.gene.expression,sd.gene.expression)
  }
}

## KO enrichment using A.thaliana orthologs
library(clusterProfiler)
library(pathview)
library(org.At.tair.db)

for (photoperiod in c("LD","SD"))
{
  for (i in seq(from=0,to=20,by=4))
  {
    target.cluster <- paste0(photoperiod,"peak_clusters/cluster_ZT",i,".txt")
    target.names.table <- read.table(target.cluster)
    target.vector <- as.vector(target.names.table$V1)
    target.arabi <- ath_orth[target.vector]
    target.arabi <- target.arabi[!(is.na(target.arabi))]
    names(target.arabi) <- NULL
    target.arabi <- target.arabi[target.arabi != ""]
    pathway.enrich <- as.data.frame(enrichKEGG(gene = unique(target.arabi),
                                   organism      = "arabidopsis",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   keyType = "kegg"))
    
    write.table(file = paste0("KO_arabi/KO_ara_ZT",i,"_",photoperiod,".tsv"),pathway.enrich)
    for (j in 1:nrow(pathway.enrich))
    {
      pathway.current.id <- pathway.enrich$ID[j]
      pathview(gene.data = target.arabi, pathway.id = pathway.current.id,species = "ath",
               gene.idtype = "TAIR", kegg.dir = paste0("kegg_images_ara/KO_ara_ZT",i,"_",photoperiod))
    }
  }
}


##CHLOROPLAST DIVISION AND ORGANIZATION##

#Identifying genes associated with chloroplast division and organization based on Arabidopsis orthologs

library(AnnotationForge)
library(org.At.tair.db)
library(org.Mpolymorpha.eg.db)
columns(org.At.tair.db)
columns(org.Mpolymorpha.eg.db)

#GO terms for this biological processes
chloro_go <- c("GO:0010020", "GO:0009658")

head(keys(org.At.tair.db,keytype="GO"))

#Selection of proper genes in Arabidopsis
chloro_ath <- select(org.At.tair.db, keys=chloro_go, columns=c("TAIR"), keytype="GO")
chloro_ath_vector <- chloro_ath$TAIR
length(chloro_ath_vector)

#Selection of Marchantia orthologs
columns(org.Mpolymorpha.eg.db)
chloro_mp <- select(org.Mpolymorpha.eg.db, keys=chloro_ath_vector, columns=c("GID"), keytype="ATH_ORTH")
chloro_mp_vector <- chloro_mp$GID
chloro_mp_vector <- unique(chloro_mp_vector[!is.na(chloro_mp_vector)])
length(unique(chloro_mp_vector))
chloro_expr <- genes.expression[chloro_mp_vector,]
chloro_expr

###Circular heatmaps
mean.chloro.ld.00 <- (chloro_expr[,1]+chloro_expr[,2])/2
mean.chloro.ld.04 <- (chloro_expr[,4]+chloro_expr[,5])/2
mean.chloro.ld.08 <- (chloro_expr[,7]+chloro_expr[,8])/2
mean.chloro.ld.12 <- (chloro_expr[,10]+chloro_expr[,11])/2
mean.chloro.ld.16 <- (chloro_expr[,13]+chloro_expr[,14])/2
mean.chloro.ld.20 <- (chloro_expr[,16]+chloro_expr[,17])/2
mean.chloro.sd.00 <- (chloro_expr[,19]+chloro_expr[,20])/2
mean.chloro.sd.04 <- (chloro_expr[,22]+chloro_expr[,23])/2
mean.chloro.sd.08 <- (chloro_expr[,25]+chloro_expr[,26])/2
mean.chloro.sd.12 <- (chloro_expr[,28]+chloro_expr[,29])/2
mean.chloro.sd.16 <- (chloro_expr[,31]+chloro_expr[,32])/2
mean.chloro.sd.20 <- (chloro_expr[,34]+chloro_expr[,35])/2
mean.ld.expression <- matrix(c(mean.chloro.ld.00, mean.chloro.ld.04, mean.chloro.ld.08, mean.chloro.ld.12, mean.chloro.ld.16, 
                               mean.chloro.ld.20), ncol = 6, byrow = F)
colnames(mean.ld.expression) <- c("zt00","zt04","zt08","zt12","zt16","zt20")
rownames(mean.ld.expression) <- rownames(chloro_expr)
mean.ld.expression 

mean.sd.expression <- matrix(c(mean.chloro.sd.00, mean.chloro.sd.04, mean.chloro.sd.08, mean.chloro.sd.12, mean.chloro.sd.16, 
                               mean.chloro.sd.20), ncol = 6, byrow = F)
colnames(mean.sd.expression) <- c("zt00","zt04","zt08","zt12","zt16","zt20")
rownames(mean.sd.expression) <- rownames(chloro_expr)

#SD
sd.gene.classification <- classify.genes(rain.results = results.sd, gene.set = chloro_mp_vector)

zt0.sd <- subset(sd.gene.classification, peak == 0)$gene
zt4.sd <- subset(sd.gene.classification, peak == 4)$gene
zt8.sd <- subset(sd.gene.classification, peak == 8)$gene
zt12.sd <- subset(sd.gene.classification, peak == 12)$gene
zt16.sd <- subset(sd.gene.classification, peak == 16)$gene
zt20.sd <- subset(sd.gene.classification, peak == 20)$gene

sd.mean.expression.sorted <- mean.sd.expression[c(zt0.sd,zt4.sd,zt8.sd,
                                                       zt12.sd,zt16.sd,zt20.sd),]

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


scaled.sd.mean.expression.sorted <- scale.and.interpolate.mean.expression(mean.expression = sd.mean.expression.sorted)
nrow(scaled.sd.mean.expression.sorted)
no.na.scaled.sd.mean.expression.sorted <- na.omit(scaled.sd.mean.expression.sorted)
nrow(no.na.scaled.sd.mean.expression.sorted)

library(gplots)
colfunc <- colorRampPalette(c("black","tomato","yellow"))

png(filename = "heatmap_SD.png",width = 400)
heatmap.2(no.na.scaled.sd.mean.expression.sorted,Colv =FALSE,trace=c("none"),dendrogram = "none",col=colfunc(100),density.info=c("none"),key=FALSE,margins = c(2,8),cexRow = 1.2,Rowv=FALSE,labCol = c(""),labRow = c(""))
dev.off()

ind.1 <- seq(from=1,by=3,to=nrow(no.na.scaled.sd.mean.expression.sorted))
length(ind.1)
ind.2 <- seq(from=2,by=3,to=nrow(no.na.scaled.sd.mean.expression.sorted))
length(ind.2)
ind.3 <- seq(from=3,by=3,to=nrow(no.na.scaled.sd.mean.expression.sorted))
length(ind.3)

sd.cir.heatmap <- (no.na.scaled.sd.mean.expression.sorted[ind.1[1:42],] +
                     no.na.scaled.sd.mean.expression.sorted[ind.2,] +
                     no.na.scaled.sd.mean.expression.sorted[ind.3,])/3 #the 42 shows the number of elements of ind.1 to be used

#install.packages("circlize")
library(circlize)
library(ComplexHeatmap)

circos.par("track.height" = 0.3, "circle.margin" = c(0.2))
col_fun_g1 = colorRamp2(c(0,0.5,1), c("darkred", "red", "greenyellow"))
circos.heatmap(t(sd.cir.heatmap), #data in a matrix
               col = col_fun_g1, #color function defined previosuly
               show.sector.labels = F, #labels of track (col names) or levels
               cluster = F) #F to avoid cluster organization 

draw.sector(start.degree = 0, end.degree = 360, clock.wise = T, rou1 = 1.1, rou2 = 1.2, center = c(0,0), col = "white", border = "red", lwd = 2)
draw.sector(start.degree = 240, end.degree = 360, clock.wise = T, rou1 = 1.1, rou2 = 1.2, center = c(0,0), col = "red", border = "red", lwd = 2)
circos.clear()


#LD
ld.gene.classification <- classify.genes(rain.results = results.ld, gene.set = chloro_mp_vector)

zt0.ld <- subset(ld.gene.classification, peak == 0)$gene
zt4.ld <- subset(ld.gene.classification, peak == 4)$gene
zt8.ld <- subset(ld.gene.classification, peak == 8)$gene
zt12.ld <- subset(ld.gene.classification, peak == 12)$gene
zt16.ld <- subset(ld.gene.classification, peak == 16)$gene
zt20.ld <- subset(ld.gene.classification, peak == 20)$gene

ld.mean.expression.sorted <- mean.ld.expression[c(zt0.ld,zt4.ld,zt8.ld,
                                                  zt12.ld,zt16.ld,zt20.ld),]

scaled.ld.mean.expression.sorted <- scale.and.interpolate.mean.expression(mean.expression = ld.mean.expression.sorted)
nrow(scaled.ld.mean.expression.sorted)
no.na.scaled.ld.mean.expression.sorted <- na.omit(scaled.ld.mean.expression.sorted)
nrow(no.na.scaled.ld.mean.expression.sorted)

colfunc <- colorRampPalette(c("black","blue","yellow"))

png(filename = "heatmap_LD.png",width = 400)
heatmap.2(no.na.scaled.ld.mean.expression.sorted,Colv =FALSE,trace=c("none"),dendrogram = "none",col=colfunc(100),density.info=c("none"),key=FALSE,margins = c(2,8),cexRow = 1.2,Rowv=FALSE,labCol = c(""),labRow = c(""))
dev.off()

ld.cir.heatmap <- (no.na.scaled.ld.mean.expression.sorted[ind.1[1:42],] +
                     no.na.scaled.ld.mean.expression.sorted[ind.2,] +
                     no.na.scaled.ld.mean.expression.sorted[ind.3,])/3

circos.par("track.height" = 0.3, "circle.margin" = c(0.2))
col_fun_g1 = colorRamp2(c(0,0.5,1), c("black", "blue", "yellow"))
circos.heatmap(t(ld.cir.heatmap), #data in a matrix
               col = col_fun_g1, #color function defined previosly
               show.sector.labels = F, #labels of track (col names) or levels
               cluster = F) #F to avoid cluster organization 

draw.sector(start.degree = 0, end.degree = 360, clock.wise = T, rou1 = 1.1, rou2 = 1.2, center = c(0,0), col = "white", border = "blue", lwd = 2)
draw.sector(start.degree = 120, end.degree = 360, clock.wise = T, rou1 = 1.1, rou2 = 1.2, center = c(0,0), col = "blue", border = "blue", lwd = 2)
circos.clear()


##Division and organization apart##

division_go <- c("GO:0010020")
org_go <- c("GO:0009658")

#Selection of proper genes in Arabidopsis
div_ath <- select(org.At.tair.db, keys=division_go, columns=c("TAIR"), keytype="GO")
div_ath_vector <- div_ath$TAIR
length(div_ath_vector)

org_ath <- select(org.At.tair.db, keys=org_go, columns=c("TAIR"), keytype="GO")
org_ath_vector <- org_ath$TAIR
length(org_ath_vector)

#Selection of Marchantia orthologs
columns(org.Mpolymorpha.eg.db)
div_mp <- select(org.Mpolymorpha.eg.db, keys=div_ath_vector, columns=c("GID"), keytype="ATH_ORTH")
div_mp_vector <- div_mp$GID
div_mp_vector <- unique(div_mp_vector[!is.na(div_mp_vector)])
length(unique(div_mp_vector))
div_expr <- genes.expression[div_mp_vector,]
div_expr

org_mp <- select(org.Mpolymorpha.eg.db, keys=org_ath_vector, columns=c("GID"), keytype="ATH_ORTH")
org_mp_vector <- org_mp$GID
org_mp_vector <- unique(org_mp_vector[!is.na(org_mp_vector)])
length(unique(org_mp_vector))
org_expr <- genes.expression[org_mp_vector,]
org_expr

###Circular heatmaps
mean.div.ld.00 <- (div_expr[,1]+div_expr[,2])/2
mean.div.ld.04 <- (div_expr[,4]+div_expr[,5])/2
mean.div.ld.08 <- (div_expr[,7]+div_expr[,8])/2
mean.div.ld.12 <- (div_expr[,10]+div_expr[,11])/2
mean.div.ld.16 <- (div_expr[,13]+div_expr[,14])/2
mean.div.ld.20 <- (div_expr[,16]+div_expr[,17])/2
mean.div.sd.00 <- (div_expr[,19]+div_expr[,20])/2
mean.div.sd.04 <- (div_expr[,22]+div_expr[,23])/2
mean.div.sd.08 <- (div_expr[,25]+div_expr[,26])/2
mean.div.sd.12 <- (div_expr[,28]+div_expr[,29])/2
mean.div.sd.16 <- (div_expr[,31]+div_expr[,32])/2
mean.div.sd.20 <- (div_expr[,34]+div_expr[,35])/2
mean.ld.div <- matrix(c(mean.div.ld.00, mean.div.ld.04, mean.div.ld.08, mean.div.ld.12, mean.div.ld.16, 
                               mean.div.ld.20), ncol = 6, byrow = F)
colnames(mean.ld.div) <- c("zt00","zt04","zt08","zt12","zt16","zt20")
rownames(mean.ld.div) <- rownames(div_expr)
mean.ld.div 

mean.sd.div <- matrix(c(mean.div.sd.00, mean.div.sd.04, mean.div.sd.08, mean.div.sd.12, mean.div.sd.16, 
                               mean.div.sd.20), ncol = 6, byrow = F)
colnames(mean.sd.div) <- c("zt00","zt04","zt08","zt12","zt16","zt20")
rownames(mean.sd.div) <- rownames(div_expr)
mean.sd.div

mean.org.ld.00 <- (org_expr[,1]+org_expr[,2])/2
mean.org.ld.04 <- (org_expr[,4]+org_expr[,5])/2
mean.org.ld.08 <- (org_expr[,7]+org_expr[,8])/2
mean.org.ld.12 <- (org_expr[,10]+org_expr[,11])/2
mean.org.ld.16 <- (org_expr[,13]+org_expr[,14])/2
mean.org.ld.20 <- (org_expr[,16]+org_expr[,17])/2
mean.org.sd.00 <- (org_expr[,19]+org_expr[,20])/2
mean.org.sd.04 <- (org_expr[,22]+org_expr[,23])/2
mean.org.sd.08 <- (org_expr[,25]+org_expr[,26])/2
mean.org.sd.12 <- (org_expr[,28]+org_expr[,29])/2
mean.org.sd.16 <- (org_expr[,31]+org_expr[,32])/2
mean.org.sd.20 <- (org_expr[,34]+org_expr[,35])/2
mean.ld.org <- matrix(c(mean.org.ld.00, mean.org.ld.04, mean.org.ld.08, mean.org.ld.12, mean.org.ld.16, 
                        mean.org.ld.20), ncol = 6, byrow = F)
colnames(mean.ld.org) <- c("zt00","zt04","zt08","zt12","zt16","zt20")
rownames(mean.ld.org) <- rownames(org_expr)
mean.ld.org

mean.sd.org <- matrix(c(mean.org.sd.00, mean.org.sd.04, mean.org.sd.08, mean.org.sd.12, mean.org.sd.16, 
                        mean.org.sd.20), ncol = 6, byrow = F)
colnames(mean.sd.org) <- c("zt00","zt04","zt08","zt12","zt16","zt20")
rownames(mean.sd.org) <- rownames(org_expr)
mean.sd.org

#SD
sd.div.classification <- classify.genes(rain.results = results.sd, gene.set = div_mp_vector)

zt0.sd.div <- subset(sd.div.classification, peak == 0)$gene
zt4.sd.div <- subset(sd.div.classification, peak == 4)$gene
zt8.sd.div <- subset(sd.div.classification, peak == 8)$gene
zt12.sd.div <- subset(sd.div.classification, peak == 12)$gene
zt16.sd.div <- subset(sd.div.classification, peak == 16)$gene
zt20.sd.div <- subset(sd.div.classification, peak == 20)$gene

sd.mean.div.sorted <- mean.sd.div[c(zt0.sd.div,zt4.sd.div,zt8.sd.div,
                                                  zt12.sd.div,zt16.sd.div,zt20.sd.div),]

scaled.sd.mean.div.sorted <- scale.and.interpolate.mean.expression(mean.expression = sd.mean.div.sorted)
nrow(scaled.sd.mean.div.sorted)
no.na.scaled.sd.mean.div.sorted <- na.omit(scaled.sd.mean.div.sorted)
nrow(no.na.scaled.sd.mean.div.sorted)



sd.org.classification <- classify.genes(rain.results = results.sd, gene.set = org_mp_vector)

zt0.sd.org <- subset(sd.org.classification, peak == 0)$gene
zt4.sd.org <- subset(sd.org.classification, peak == 4)$gene
zt8.sd.org <- subset(sd.org.classification, peak == 8)$gene
zt12.sd.org <- subset(sd.org.classification, peak == 12)$gene
zt16.sd.org <- subset(sd.org.classification, peak == 16)$gene
zt20.sd.org <- subset(sd.org.classification, peak == 20)$gene

sd.mean.org.sorted <- mean.sd.org[c(zt0.sd.org,zt4.sd.org,zt8.sd.org,
                                    zt12.sd.org,zt16.sd.org,zt20.sd.org),]

scaled.sd.mean.org.sorted <- scale.and.interpolate.mean.expression(mean.expression = sd.mean.org.sorted)
nrow(scaled.sd.mean.org.sorted)
no.na.scaled.sd.mean.org.sorted <- na.omit(scaled.sd.mean.org.sorted)
nrow(no.na.scaled.sd.mean.org.sorted)


#ind.1.div <- seq(from=1,by=3,to=nrow(no.na.scaled.sd.mean.div.sorted))
#length(ind.1.div)
#ind.2.div <- seq(from=2,by=3,to=nrow(no.na.scaled.sd.mean.div.sorted))
#length(ind.2.div)
#ind.3.div <- seq(from=3,by=3,to=nrow(no.na.scaled.sd.mean.div.sorted))
#length(ind.3.div)

#sd.cir.heatmap.div <- (no.na.scaled.sd.mean.div.sorted[ind.1.div[1:6],] +
#                     no.na.scaled.sd.mean.div.sorted[ind.2.div,] +
#                     no.na.scaled.sd.mean.div.sorted[ind.3.div,])/3 #the 6 shows the number of elements of ind.1 to be used


#ind.1.org <- seq(from=1,by=3,to=nrow(no.na.scaled.sd.mean.org.sorted))
#length(ind.1.org)
#ind.2.org <- seq(from=2,by=3,to=nrow(no.na.scaled.sd.mean.org.sorted))
#length(ind.2.org)
#ind.3.org <- seq(from=3,by=3,to=nrow(no.na.scaled.sd.mean.org.sorted))
#length(ind.3.org)

#sd.cir.heatmap.org <- (no.na.scaled.sd.mean.org.sorted[ind.1.org[1:38],] +
#                         no.na.scaled.sd.mean.org.sorted[ind.2.org,] +
#                         no.na.scaled.sd.mean.org.sorted[ind.3.org,])/3


circos.par("track.height" = 0.3, "circle.margin" = c(0.2))
col_fun_g1 = colorRamp2(c(0,0.5,1), c("darkred", "red", "greenyellow"))
circos.heatmap(t(no.na.scaled.sd.mean.org.sorted), #data in a matrix
               col = col_fun_g1, #color function defined previosuly
               show.sector.labels = F, #labels of track (col names) or levels
               cluster = F) #F to avoid cluster organization 
draw.sector(start.degree = 0, end.degree = 360, clock.wise = T, rou1 = 0.69, rou2 = 0.7, center = c(0,0), col = "black", border = "black", lwd = 2)
circos.par("track.height" = 0.2)
circos.heatmap(t(no.na.scaled.sd.mean.div.sorted), #data in a matrix
               col = col_fun_g1, #color function defined previosuly
               show.sector.labels = F, #labels of track (col names) or levels
               cluster = F) #F to avoid cluster organization

draw.sector(start.degree = 0, end.degree = 360, clock.wise = T, rou1 = 1.1, rou2 = 1.2, center = c(0,0), col = "white", border = "red", lwd = 2)
draw.sector(start.degree = 240, end.degree = 360, clock.wise = T, rou1 = 1.1, rou2 = 1.2, center = c(0,0), col = "red", border = "red", lwd = 2)
circos.clear()

lgd = Legend(title = "", col_fun = col_fun_g1, direction = "vertical", grid_width = unit(10, "mm"), grid_height = unit(12,"mm"))
grid.draw(lgd)
grid.newpage()

#LD
ld.div.classification <- classify.genes(rain.results = results.ld, gene.set = div_mp_vector)

zt0.ld.div <- subset(ld.div.classification, peak == 0)$gene
zt4.ld.div <- subset(ld.div.classification, peak == 4)$gene
zt8.ld.div <- subset(ld.div.classification, peak == 8)$gene
zt12.ld.div <- subset(ld.div.classification, peak == 12)$gene
zt16.ld.div <- subset(ld.div.classification, peak == 16)$gene
zt20.ld.div <- subset(ld.div.classification, peak == 20)$gene

ld.mean.div.sorted <- mean.ld.div[c(zt0.ld.div,zt4.ld.div,zt8.ld.div,
                                    zt12.ld.div,zt16.ld.div,zt20.ld.div),]

scaled.ld.mean.div.sorted <- scale.and.interpolate.mean.expression(mean.expression = ld.mean.div.sorted)
nrow(scaled.ld.mean.div.sorted)
no.na.scaled.ld.mean.div.sorted <- na.omit(scaled.ld.mean.div.sorted)
nrow(no.na.scaled.ld.mean.div.sorted)



ld.org.classification <- classify.genes(rain.results = results.ld, gene.set = org_mp_vector)

zt0.ld.org <- subset(ld.org.classification, peak == 0)$gene
zt4.ld.org <- subset(ld.org.classification, peak == 4)$gene
zt8.ld.org <- subset(ld.org.classification, peak == 8)$gene
zt12.ld.org <- subset(ld.org.classification, peak == 12)$gene
zt16.ld.org <- subset(ld.org.classification, peak == 16)$gene
zt20.ld.org <- subset(ld.org.classification, peak == 20)$gene

ld.mean.org.sorted <- mean.ld.org[c(zt0.ld.org,zt4.ld.org,zt8.ld.org,
                                    zt12.ld.org,zt16.ld.org,zt20.ld.org),]

scaled.ld.mean.org.sorted <- scale.and.interpolate.mean.expression(mean.expression = ld.mean.org.sorted)
nrow(scaled.ld.mean.org.sorted)
no.na.scaled.ld.mean.org.sorted <- na.omit(scaled.ld.mean.org.sorted)
nrow(no.na.scaled.ld.mean.org.sorted)


#ld.cir.heatmap.div <- (no.na.scaled.ld.mean.div.sorted[ind.1.div[1:6],] +
#                         no.na.scaled.ld.mean.div.sorted[ind.2.div,] +
#                         no.na.scaled.ld.mean.div.sorted[ind.3.div,])/3 #the 6 shows the number of elements of ind.1 to be used


#ld.cir.heatmap.org <- (no.na.scaled.ld.mean.org.sorted[ind.1.org[1:38],] +
#                         no.na.scaled.ld.mean.org.sorted[ind.2.org,] +
#                         no.na.scaled.ld.mean.org.sorted[ind.3.org,])/3


circos.par("track.height" = 0.3, "circle.margin" = c(0.2))
col_fun_g1 = colorRamp2(c(0,0.5,1), c("black", "blue", "yellow"))
circos.heatmap(t(no.na.scaled.ld.mean.org.sorted), #data in a matrix
               col = col_fun_g1, #color function defined previosuly
               show.sector.labels = F, #labels of track (col names) or levels
               cluster = F) #F to avoid cluster organization 
draw.sector(start.degree = 0, end.degree = 360, clock.wise = T, rou1 = 0.69, rou2 = 0.7, center = c(0,0), col = "navy", border = "navy", lwd = 2)
circos.heatmap(t(no.na.scaled.ld.mean.div.sorted), #data in a matrix
               col = col_fun_g1, #color function defined previosuly
               show.sector.labels = F, #labels of track (col names) or levels
               cluster = F) #F to avoid cluster organization

draw.sector(start.degree = 0, end.degree = 360, clock.wise = T, rou1 = 1.1, rou2 = 1.2, center = c(0,0), col = "white", border = "blue", lwd = 2)
draw.sector(start.degree = 120, end.degree = 360, clock.wise = T, rou1 = 1.1, rou2 = 1.2, center = c(0,0), col = "blue", border = "blue", lwd = 2)
circos.clear()

lgd = Legend(title = "", col_fun = col_fun_g1, direction = "vertical", grid_width = unit(10, "mm"), grid_height = unit(12,"mm"))
grid.draw(lgd)
grid.newpage()

library(ComplexHeatmap)


##Electron flow##

plot.ld.sd(gene.id = "Mp7g15130",gene.name = "STN7", ld.gene.expression,sd.gene.expression)
plot.ld.sd(gene.id = "Mp7g16020",gene.name = "STN8", ld.gene.expression,sd.gene.expression)
plot.ld.sd(gene.id = "Mp6g19810",gene.name = "PPH1", ld.gene.expression,sd.gene.expression)
plot.ld.sd(gene.id = "Mp3g23840",gene.name = "PGR5", ld.gene.expression,sd.gene.expression)
plot.ld.sd(gene.id = "Mp5g09890",gene.name = "PGRL1", ld.gene.expression,sd.gene.expression)

plot.ld.sd.ll(gene.id = "Mp7g15130",gene.name = "STN7", genes.expression)
plot.ld.sd.ll(gene.id = "Mp7g16020",gene.name = "STN8", genes.expression)
plot.ld.sd.ll(gene.id = "Mp6g19810",gene.name = "PPH1", genes.expression)
plot.ld.sd.ll(gene.id = "Mp3g23840",gene.name = "PGR5", genes.expression)
plot.ld.sd.ll(gene.id = "Mp5g09890",gene.name = "PGRL1", genes.expression)

flow_vector <- c("Mp7g15130","Mp7g16020","Mp6g19810","Mp3g23840","Mp5g09890")
flow_expr <- genes.expression[flow_vector,]

mean.flow.ld.00 <- (flow_expr[,1]+flow_expr[,2])/2
mean.flow.ld.04 <- (flow_expr[,4]+flow_expr[,5])/2
mean.flow.ld.08 <- (flow_expr[,7]+flow_expr[,8])/2
mean.flow.ld.12 <- (flow_expr[,10]+flow_expr[,11])/2
mean.flow.ld.16 <- (flow_expr[,13]+flow_expr[,14])/2
mean.flow.ld.20 <- (flow_expr[,16]+flow_expr[,17])/2
mean.flow.sd.00 <- (flow_expr[,19]+flow_expr[,20])/2
mean.flow.sd.04 <- (flow_expr[,22]+flow_expr[,23])/2
mean.flow.sd.08 <- (flow_expr[,25]+flow_expr[,26])/2
mean.flow.sd.12 <- (flow_expr[,28]+flow_expr[,29])/2
mean.flow.sd.16 <- (flow_expr[,31]+flow_expr[,32])/2
mean.flow.sd.20 <- (flow_expr[,34]+flow_expr[,35])/2
mean.ld.flow <- matrix(c(mean.flow.ld.00, mean.flow.ld.04, mean.flow.ld.08, mean.flow.ld.12, mean.flow.ld.16, 
                        mean.flow.ld.20), ncol = 6, byrow = F)
colnames(mean.ld.flow) <- c("zt00","zt04","zt08","zt12","zt16","zt20")
rownames(mean.ld.flow) <- rownames(flow_expr)
mean.ld.flow

mean.sd.flow <- matrix(c(mean.flow.sd.00, mean.flow.sd.04, mean.flow.sd.08, mean.flow.sd.12, mean.flow.sd.16, 
                        mean.flow.sd.20), ncol = 6, byrow = F)
colnames(mean.sd.flow) <- c("zt00","zt04","zt08","zt12","zt16","zt20")
rownames(mean.sd.flow) <- rownames(flow_expr)
mean.sd.flow

#SD
sd.flow.classification <- classify.genes(rain.results = results.sd, gene.set = flow_vector)

zt0.sd.flow <- subset(sd.flow.classification, peak == 0)$gene
zt4.sd.flow <- subset(sd.flow.classification, peak == 4)$gene
zt8.sd.flow <- subset(sd.flow.classification, peak == 8)$gene
zt12.sd.flow <- subset(sd.flow.classification, peak == 12)$gene
zt16.sd.flow <- subset(sd.flow.classification, peak == 16)$gene
zt20.sd.flow <- subset(sd.flow.classification, peak == 20)$gene

sd.mean.flow.sorted <- mean.sd.flow[c(zt0.sd.flow,zt4.sd.flow,zt8.sd.flow,
                                    zt12.sd.flow,zt16.sd.flow,zt20.sd.flow),]

scaled.sd.mean.flow.sorted <- scale.and.interpolate.mean.expression(mean.expression = sd.mean.flow.sorted)
nrow(scaled.sd.mean.flow.sorted)
no.na.scaled.sd.mean.flow.sorted <- na.omit(scaled.sd.mean.flow.sorted)
nrow(no.na.scaled.sd.mean.flow.sorted)

library(circlize)
library(ComplexHeatmap)
circos.par("track.height" = 0.3, "circle.margin" = c(0.2))
col_fun_g1 = colorRamp2(c(0,0.5,1), c("darkred", "red", "greenyellow"))
circos.heatmap(t(no.na.scaled.sd.mean.flow.sorted), #data in a matrix
               col = col_fun_g1, #color function defined previosuly
               show.sector.labels = F, #labels of track (col names) or levels
               cluster = F) #F to avoid cluster organization 

draw.sector(start.degree = 0, end.degree = 360, clock.wise = T, rou1 = 1.1, rou2 = 1.2, center = c(0,0), col = "white", border = "red", lwd = 2)
draw.sector(start.degree = 240, end.degree = 360, clock.wise = T, rou1 = 1.1, rou2 = 1.2, center = c(0,0), col = "red", border = "red", lwd = 2)
circos.clear()
grid.newpage()

#LD
ld.flow.classification <- classify.genes(rain.results = results.ld, gene.set = flow_vector)

zt0.ld.flow <- subset(ld.flow.classification, peak == 0)$gene
zt4.ld.flow <- subset(ld.flow.classification, peak == 4)$gene
zt8.ld.flow <- subset(ld.flow.classification, peak == 8)$gene
zt12.ld.flow <- subset(ld.flow.classification, peak == 12)$gene
zt16.ld.flow <- subset(ld.flow.classification, peak == 16)$gene
zt20.ld.flow <- subset(ld.flow.classification, peak == 20)$gene

ld.mean.flow.sorted <- mean.ld.flow[c(zt0.ld.flow,zt4.ld.flow,zt8.ld.flow,
                                      zt12.ld.flow,zt16.ld.flow,zt20.ld.flow),]

scaled.ld.mean.flow.sorted <- scale.and.interpolate.mean.expression(mean.expression = ld.mean.flow.sorted)
nrow(scaled.ld.mean.flow.sorted)
no.na.scaled.ld.mean.flow.sorted <- na.omit(scaled.ld.mean.flow.sorted)
nrow(no.na.scaled.ld.mean.flow.sorted)

circos.par("track.height" = 0.3, "circle.margin" = c(0.2))
col_fun_g1 = colorRamp2(c(0,0.5,1), c("black", "blue", "yellow"))
circos.heatmap(t(no.na.scaled.ld.mean.flow.sorted), #data in a matrix
               col = col_fun_g1, #color function defined previosuly
               show.sector.labels = F, #labels of track (col names) or levels
               cluster = F) #F to avoid cluster organization 

draw.sector(start.degree = 0, end.degree = 360, clock.wise = T, rou1 = 1.1, rou2 = 1.2, center = c(0,0), col = "white", border = "blue", lwd = 2)
draw.sector(start.degree = 120, end.degree = 360, clock.wise = T, rou1 = 1.1, rou2 = 1.2, center = c(0,0), col = "blue", border = "blue", lwd = 2)
circos.clear()
grid.newpage()

plot.ld.sd(gene.id = "Mp8g05910", gene.name = "PsbS", ld.gene.expression, sd.gene.expression)



### Genes involved in photoreolcation of chloroplasts

# MpPHOT (Mp5g03810) expression trough LD and SD
png("../../../Figs_for_Windows/Figcloro/phot.png", width = 700, height = 400)
plot.ld.sd("Mp5g03810", "MpPHOT", ld.gene.expression, sd.gene.expression)
dev.off()
"Mp5g03810" %in% rhythmic.genes.ld
"Mp5g03810" %in% rhythmic.genes.sd
plot.ld.sd.ll("Mp5g03810", "MpPHOT", genes.expression)
# MpPHOT shows a rhythmic expression pattern in both photoperiods, with marked reductions at 
# dusk, which could explain the relocation of chloroplasts from the periclinal to the anticlinal
# cell walls during the dark periods
# 

# CHUP 1 (Mp7g03540) expression trough LD and SD
plot.ld.sd("Mp7g03540", "CHUP1", ld.gene.expression, sd.gene.expression)
"Mp7g03540" %in% rhythmic.genes.ld
"Mp7g03540" %in% rhythmic.genes.sd
plot.ld.sd.ll("Mp7g03540", "CHUP1", genes.expression)

# KAC1 (Mp3g03760) expression trough LD and SD
plot.ld.sd("Mp3g03760", "KAC1", ld.gene.expression, sd.gene.expression)
"Mp3g03760" %in% rhythmic.genes.ld
"Mp3g03760" %in% rhythmic.genes.sd
plot.ld.sd.ll("Mp3g03760", "KAC1", genes.expression)
# Shows rhythmicity in SD but not in LD

# BLUS1 (Mp8g10870) expression trough LD and SD
plot.ld.sd("Mp8g10870", "BLUS1", ld.gene.expression, sd.gene.expression)
"Mp8g10870" %in% rhythmic.genes.ld
"Mp8g10870" %in% rhythmic.genes.sd
plot.ld.sd.ll("Mp8g10870", "BLUS1", genes.expression)
# Shows rhythmicity in LD but not in SD

# THRUMIN1 (Mp1g07590 and Mp3g04860) expression trough LD and SD
plot.ld.sd("Mp1g07590", "THRUMIN1", ld.gene.expression, sd.gene.expression)
"Mp1g07590" %in% rhythmic.genes.ld
"Mp1g07590" %in% rhythmic.genes.sd
plot.ld.sd.ll("Mp1g07590", "THRUMIN1", genes.expression)
# Shows rhythmicity in SD but not in LD

plot.ld.sd("Mp3g04860", "THRUMIN1", ld.gene.expression, sd.gene.expression)
"Mp3g04860" %in% rhythmic.genes.ld
"Mp3g04860" %in% rhythmic.genes.sd
plot.ld.sd.ll("Mp3g04860", "THRUMIN1", genes.expression)
# Shows rhythmicity in SD but not in LD

# JAC1 (Mp6g17960) expression trough LD and SD
plot.ld.sd("Mp6g17960", "JAC1", ld.gene.expression, sd.gene.expression)
"Mp6g17960" %in% rhythmic.genes.ld
"Mp6g17960" %in% rhythmic.genes.sd
plot.ld.sd.ll("Mp6g17960", "JAC1", genes.expression)
# Does not show a rhythmic behaviour beacause of a lack of expression in these
# conditions

# Homologs for OEP7, WEB1, PMI1 and PMI2 were not found in Marchantiaa's genome
