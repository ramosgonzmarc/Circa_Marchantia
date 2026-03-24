
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
library(enrichplot)
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

barplot(ego.low.exp.genes, showCategory = 20)

png("ego_low_exp_genes.png", res = 230, width = 1700, height = 1700)
emapplot(pairwise_termsim(ego.low.exp.genes))
dev.off()


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


# Global exploration of gene expression in each condition
mean(ld.gene.expression)
mean(ld.ll.gene.expression)
mean(ld.dd.gene.expression)
mean(sd.gene.expression)
mean(sd.ll.gene.expression)
mean(sd.dd.gene.expression)

median(ld.gene.expression)
median(ld.ll.gene.expression)
median(ld.dd.gene.expression)
median(sd.gene.expression)
median(sd.ll.gene.expression)
median(sd.dd.gene.expression)

# Plot expression
boxplot(log2(as.numeric(apply(ld.gene.expression, MARGIN = 1, mean))+1),
        log2(as.numeric(apply(ld.ll.gene.expression, MARGIN = 1, mean))+1),
        log2(as.numeric(apply(ld.dd.gene.expression, MARGIN = 1, mean))+1),
        log2(as.numeric(apply(sd.gene.expression, MARGIN = 1, mean))+1),
        log2(as.numeric(apply(sd.ll.gene.expression, MARGIN = 1, mean))+1),
        log2(as.numeric(apply(sd.dd.gene.expression, MARGIN = 1, mean))+1))




######### Statistical analysis of circadian genes with CircaCompare  ###########

# First of all, define new expressed-genes background set
bg_set <- names(which(apply(X = genes.expression,MARGIN = 1,FUN = max) >= 2))
bg_table <- genes.expression[bg_set,]
dim(bg_table)

# Determine rhythmic genes in LD 
# With 24 h period
library(circacompare)

results.ld <- matrix(nrow=length(bg_set),ncol=9)
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
    #results.ld[i,] <- out.i[[2]][,2]
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    results.ld[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    results.ld[i,] <- NA
  }
}

colnames(results.ld) <- c(out.i[[2]][,1], "error_mesor", "error_amp", "error_phase")
write.table(results.ld,file="circacompare_ld_24_with_errors.tsv", sep="\t")

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

results.sd <- matrix(nrow=length(bg_set),ncol=9)
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
    #results.sd[i,] <- out.i[[2]][,2]
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    results.sd[i,] <- c(out.i[[2]][,2], error.i)
    
  }
  else
  {
    results.sd[i,] <- NA
  }
}

colnames(results.sd) <- c(out.i[[2]][,1], "error_mesor", "error_amp", "error_phase")
write.table(results.sd,file="circacompare_sd_12_with_errors.tsv", sep="\t")

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

png("expression_vs_rhythmicity.png", width = 1300, height = 1500, res=300)
par(lwd=3)
boxplot(rhythmic.max.expression, non.rhythmic.max.expression, outline=F,col=c("chartreuse4","royalblue3"),
        names=c("Rhythmic","Arrhythmic"),ylab="FPKM",cex.lab=1.5,las=1, cex.axis=1.5, main= "Distribution of expression levels")
dev.off()

wilcox.test(rhythmic.max.expression,non.rhythmic.max.expression,alternative="greater")$p.value



##### Difference between LD and SD in a parametric manner, complete set
library(circacompare)

results.ld.sd <- matrix(nrow=length(bg_set),ncol=21)
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
    
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    results.ld.sd[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    results.ld.sd[i,] <- NA
  }
}

# k, k1, alpha, alpha1, phi, phi1
colnames(results.ld.sd) <- c(out.i[[2]][,1], "error_mesor_1", "error_mesor_2",
                             "error_amp_1", "error_amp_2", "error_phase_1", 
                             "error_phase_2")
write.table(results.ld.sd,file="circacompare_ld_sd_complete_with_errors.tsv", sep="\t")

##### Difference between LD and SD in a parametric manner, rhythmic in both
library(circacompare)

complete.ld.rhythmic.genes <- read.table("tables/ld.rhythmic.table.tsv", header = T)$V1
complete.sd.rhythmic.genes <- read.table("tables/sd.rhythmic.table.tsv", header = T)$V1

intersect_set <- intersect(complete.ld.rhythmic.genes, complete.sd.rhythmic.genes)
intersect_table <- genes.expression[intersect_set,]
dim(intersect_table)

results.ld.sd.intersect <- matrix(nrow=length(intersect_set),ncol=21)
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
    
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    
    results.ld.sd.intersect[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    results.ld.sd.intersect[i,] <- NA
  }
}

colnames(results.ld.sd.intersect) <- c(out.i[[2]][,1], "error_mesor_1", "error_mesor_2",
                                       "error_amp_1", "error_amp_2", "error_phase_1", 
                                       "error_phase_2")
write.table(results.ld.sd.intersect,file="circacompare_ld_sd_intersection_with_errors.tsv", sep="\t")


################# Free-running ######################33
# LD - LL 24 h
# Only for rhythmic genes in LD
complete.ld.rhythmic.genes <- read.table("tables/ld.rhythmic.table.tsv", header = T)$V1
head(complete.ld.rhythmic.genes)
ld_set <- complete.ld.rhythmic.genes
ld_set_table <- genes.expression[ld_set,]
dim(ld_set_table)

results.ld.ll <- matrix(nrow=length(ld_set),ncol=9)
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
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    
    results.ld.ll[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    results.ld.ll[i,] <- NA
  }
}
colnames(results.ld.ll) <- c(out.i[[2]][,1], "error_mesor", "error_amp", "error_phase")
write.table(results.ld.ll,file="circacompare_ld_ll_with_errors.tsv", sep="\t")

# LD - LL 12 h
# Only for rhythmic genes in LD
ld_set <- complete.ld.rhythmic.genes
ld_set_table <- genes.expression[ld_set,]
dim(ld_set_table)

results.ld.ll <- matrix(nrow=length(ld_set),ncol=9)
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
    
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    
    results.ld.ll[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    results.ld.ll[i,] <- NA
  }
}
colnames(results.ld.ll) <- c(out.i[[2]][,1], "error_mesor", "error_amp", "error_phase")
write.table(results.ld.ll,file="circacompare_ld_ll_12_with_errors.tsv", sep="\t")

# LD - DD 24 h
# Only for rhythmic genes in LD
ld_set <- complete.ld.rhythmic.genes
ld_set_table <- genes.expression[ld_set,]
dim(ld_set_table)

results.ld.dd <- matrix(nrow=length(ld_set),ncol=9)
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
    
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    
    results.ld.dd[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    results.ld.dd[i,] <- NA
  }
}
colnames(results.ld.dd) <- c(out.i[[2]][,1], "error_mesor", "error_amp", "error_phase")
write.table(results.ld.dd,file="circacompare_ld_dd_with_errors.tsv", sep="\t")

# LD - DD 12 h
# Only for rhythmic genes in LD
ld_set <- complete.ld.rhythmic.genes
ld_set_table <- genes.expression[ld_set,]
dim(ld_set_table)

results.ld.dd <- matrix(nrow=length(ld_set),ncol=9)
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
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    
    results.ld.dd[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    results.ld.dd[i,] <- NA
  }
}
colnames(results.ld.dd) <- c(out.i[[2]][,1], "error_mesor", "error_amp", "error_phase")
write.table(results.ld.dd,file="circacompare_ld_dd_12_with_errors.tsv", sep="\t")

# SD - LL 24 h
# Only for rhythmic genes in SD
complete.sd.rhythmic.genes <- read.table("tables/sd.rhythmic.table.tsv", header = T)$V1
head(complete.sd.rhythmic.genes)
sd_set <- complete.sd.rhythmic.genes
sd_set_table <- genes.expression[sd_set,]
dim(sd_set_table)

results.sd.ll <- matrix(nrow=length(sd_set),ncol=9)
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
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    
    results.sd.ll[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    results.sd.ll[i,] <- NA
  }
}
colnames(results.sd.ll) <- c(out.i[[2]][,1], "error_mesor", "error_amp", "error_phase")
write.table(results.sd.ll,file="circacompare_sd_ll_with_errors.tsv", sep="\t")

# SD - LL 12 h
# Only for rhythmic genes in SD
sd_set <- complete.sd.rhythmic.genes
sd_set_table <- genes.expression[sd_set,]
dim(sd_set_table)

results.sd.ll <- matrix(nrow=length(sd_set),ncol=9)
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
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    
    results.sd.ll[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    results.sd.ll[i,] <- NA
  }
}
colnames(results.sd.ll) <- c(out.i[[2]][,1], "error_mesor", "error_amp", "error_phase")
write.table(results.sd.ll,file="circacompare_sd_ll_12_with_errors.tsv", sep="\t")

# SD - DD 24 h
# Only for rhythmic genes in SD
sd_set <- complete.sd.rhythmic.genes
sd_set_table <- genes.expression[sd_set,]
dim(sd_set_table)

results.sd.dd <- matrix(nrow=length(sd_set),ncol=9)
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
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    
    results.sd.dd[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    results.sd.dd[i,] <- NA
  }
}
colnames(results.sd.dd) <- c(out.i[[2]][,1], "error_mesor", "error_amp", "error_phase")
write.table(results.sd.dd,file="circacompare_sd_dd_with_errors.tsv", sep="\t")

# SD - DD 12 h
# Only for rhythmic genes in SD
sd_set <- complete.sd.rhythmic.genes
sd_set_table <- genes.expression[sd_set,]
dim(sd_set_table)

results.sd.dd <- matrix(nrow=length(sd_set),ncol=9)
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
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    
    results.sd.dd[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    results.sd.dd[i,] <- NA
  }
}
colnames(results.sd.dd) <- c(out.i[[2]][,1], "error_mesor", "error_amp", "error_phase")
write.table(results.sd.dd,file="circacompare_sd_dd_12_with_errors.tsv", sep="\t")



############# Exploration of free-running results
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

png("hist_phase_ld_ll.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(results.ld.ll.24[unique(c(rhythmic.genes.ld.ll.24, 
                               rhythmic.genes.ld.ll.12)),"peak_time_hours" ],
     col="lightblue", ylim=c(0,400),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
dev.off()


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

png("hist_phase_ld_dd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(results.ld.dd.24[unique(c(rhythmic.genes.ld.dd.24, 
                               rhythmic.genes.ld.dd.12)),"peak_time_hours" ],
     col="darkblue", ylim=c(0,400),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
dev.off()

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

png("hist_phase_sd_ll.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(results.sd.ll.24[unique(c(rhythmic.genes.sd.ll.24, 
                               rhythmic.genes.sd.ll.12)),"peak_time_hours" ],
     col="lightsalmon", ylim=c(0,400),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
dev.off()

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

png("hist_phase_sd_dd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(results.sd.dd.24[unique(c(rhythmic.genes.sd.dd.24, 
                               rhythmic.genes.sd.dd.12)),"peak_time_hours" ],
     col="firebrick4", ylim=c(0,400),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
dev.off()


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
length(rhythmic.genes.sd.12)
length(rhythmic.genes.ld.12)

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
circacompare.ld.ll <- matrix(nrow=length(complete.ld.ll.rhythmic.genes),ncol=21)
rownames(circacompare.ld.ll) <- complete.ld.ll.rhythmic.genes

for(i in 1:length(complete.ld.ll.rhythmic.genes))
{
  gene.i <- complete.ld.ll.rhythmic.genes[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
  ll.expression.i <- genes.expression[gene.i,c(paste(ll.ld.zt,1,sep="_"),paste(ll.ld.zt,2,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 18)
  
  if ((length(unique(ld.expression.i)) > 6) & (length(unique(ll.expression.i)) > 3))
  {
    ld.ll.df <- data.frame(time=c(time.points,time.points[1:12]),
                           measure=c(ld.expression.i, ll.expression.i),
                           group=c(rep("ld",18),rep("ll",12)))
    
    out.i <- circacompare(x = ld.ll.df, col_time = "time", col_group = "group",
                          col_outcome = "measure",alpha_threshold = 1)
    #circacompare.ld.ll[i,] <- out.i[[2]][,2]
    
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    circacompare.ld.ll[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    circacompare.ld.ll[i,] <- NA
  }
  
}

colnames(circacompare.ld.ll) <- c(out.i[[2]][,1], "error_mesor_1", "error_mesor_2",
                             "error_amp_1", "error_amp_2", "error_phase_1", 
                             "error_phase_2")
head(circacompare.ld.ll)
write.table(circacompare.ld.ll,file="tables/estimates_ld_ll_with_errors.tsv", sep="\t")

# LD to DD
library(circacompare)

complete.ld.dd.rhythmic.genes <- read.table(file="tables/ld.dd.rhythmic.table.tsv", sep="\t")$x
circacompare.ld.dd <- matrix(nrow=length(complete.ld.dd.rhythmic.genes),ncol=21)
rownames(circacompare.ld.dd) <- complete.ld.dd.rhythmic.genes

for(i in 1:length(complete.ld.dd.rhythmic.genes))
{
  gene.i <- complete.ld.dd.rhythmic.genes[i]
  
  ld.expression.i <- genes.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
  dd.expression.i <- genes.expression[gene.i,c(paste(dd.ld.zt,1,sep="_"),paste(dd.ld.zt,2,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 18)
  
  if ((length(unique(ld.expression.i)) > 6) & (length(unique(dd.expression.i)) > 3))
  {
    ld.dd.df <- data.frame(time=c(time.points,time.points[1:12]),
                           measure=c(ld.expression.i, dd.expression.i),
                           group=c(rep("ld",18),rep("dd",12)))
    
    out.i <- circacompare(x = ld.dd.df, col_time = "time", col_group = "group",
                          col_outcome = "measure",alpha_threshold = 1)
    #circacompare.ld.dd[i,] <- out.i[[2]][,2]
    
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    circacompare.ld.dd[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    circacompare.ld.dd[i,] <- NA
  }
  
}

colnames(circacompare.ld.dd) <- c(out.i[[2]][,1], "error_mesor_1", "error_mesor_2",
                                  "error_amp_1", "error_amp_2", "error_phase_1", 
                                  "error_phase_2")
write.table(circacompare.ld.dd,file="tables/estimates_ld_dd_with_errors.tsv", sep="\t")


# SD to LL
library(circacompare)

complete.sd.ll.rhythmic.genes <- read.table(file="tables/sd.ll.rhythmic.table.tsv", sep="\t")$x
circacompare.sd.ll <- matrix(nrow=length(complete.sd.ll.rhythmic.genes),ncol=21)
rownames(circacompare.sd.ll) <- complete.sd.ll.rhythmic.genes

for(i in 1:length(complete.sd.ll.rhythmic.genes))
{
  gene.i <- complete.sd.ll.rhythmic.genes[i]
  
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
  ll.expression.i <- genes.expression[gene.i,c(paste(ll.sd.zt,1,sep="_"),paste(ll.sd.zt,2,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 18)
  
  if ((length(unique(sd.expression.i)) > 6) & (length(unique(ll.expression.i)) > 3))
  {
    sd.ll.df <- data.frame(time=c(time.points,time.points[1:12]),
                           measure=c(sd.expression.i, ll.expression.i),
                           group=c(rep("sd",18),rep("ll",12)))
    
    out.i <- circacompare(x = sd.ll.df, col_time = "time", col_group = "group",
                          col_outcome = "measure",alpha_threshold = 1)
    #circacompare.sd.ll[i,] <- out.i[[2]][,2]
    
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    circacompare.sd.ll[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    circacompare.sd.ll[i,] <- NA
  }
  
}

colnames(circacompare.sd.ll) <- c(out.i[[2]][,1], "error_mesor_1", "error_mesor_2",
                                  "error_amp_1", "error_amp_2", "error_phase_1", 
                                  "error_phase_2")
head(circacompare.sd.ll)
write.table(circacompare.sd.ll,file="tables/estimates_sd_ll_with_errors.tsv", sep="\t")

# SD to DD
library(circacompare)

complete.sd.dd.rhythmic.genes <- read.table(file="tables/sd.dd.rhythmic.table.tsv", sep="\t")$x
circacompare.sd.dd <- matrix(nrow=length(complete.sd.dd.rhythmic.genes),ncol=21)
rownames(circacompare.sd.dd) <- complete.sd.dd.rhythmic.genes

for(i in 1:length(complete.sd.dd.rhythmic.genes))
{
  gene.i <- complete.sd.dd.rhythmic.genes[i]
  
  sd.expression.i <- genes.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
  dd.expression.i <- genes.expression[gene.i,c(paste(dd.sd.zt,1,sep="_"),paste(dd.sd.zt,2,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 18)
  
  if ((length(unique(sd.expression.i)) > 6) & (length(unique(dd.expression.i)) > 3))
  {
    sd.dd.df <- data.frame(time=c(time.points,time.points[1:12]),
                           measure=c(sd.expression.i, dd.expression.i),
                           group=c(rep("sd",18),rep("dd",12)))
    
    out.i <- circacompare(x = sd.dd.df, col_time = "time", col_group = "group",
                          col_outcome = "measure",alpha_threshold = 1)
    #circacompare.sd.dd[i,] <- out.i[[2]][,2]
    error.i <- as.data.frame(circacompare:::extract_model_coefs(out.i$fit))[,2]
    circacompare.sd.dd[i,] <- c(out.i[[2]][,2], error.i)
  }
  else
  {
    circacompare.sd.dd[i,] <- NA
  }
  
}

colnames(circacompare.sd.dd) <- c(out.i[[2]][,1], "error_mesor_1", "error_mesor_2",
                                  "error_amp_1", "error_amp_2", "error_phase_1", 
                                  "error_phase_2")
head(circacompare.sd.dd)
write.table(circacompare.sd.dd,file="tables/estimates_sd_dd_with_errors.tsv", sep="\t")

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

wilcox.test(x = circacompare.sd.ll[,"sd mesor estimate" ],y = circacompare.sd.ll[,"ll mesor estimate" ], alternative = "greater")$p.value
wilcox.test(x = circacompare.sd.ll[,"sd mesor estimate" ],y = circacompare.sd.dd[,"dd mesor estimate" ], alternative = "g")$p.value
wilcox.test(x = circacompare.sd.ll[,"ll mesor estimate" ],y = circacompare.sd.dd[,"dd mesor estimate" ], alternative = "g")$p.value


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
wilcox.test(x = rowMeans(ld.genes.exp.in.ld), y = rowMeans(ld.genes.exp.in.ll), alternative = "g")$p.value # LL lesser
wilcox.test(x = rowMeans(ld.genes.exp.in.ld),y = rowMeans(ld.genes.exp.in.dd), alternative = "l")$p.value ## DD greater

mean(sd.genes.exp.in.sd)
mean(sd.genes.exp.in.ll)
mean(sd.genes.exp.in.dd)
wilcox.test(x = rowMeans(sd.genes.exp.in.sd), y = rowMeans(sd.genes.exp.in.ll), alternative = "g")$p.value # LL lesser
wilcox.test(x = rowMeans(sd.genes.exp.in.sd),y = rowMeans(sd.genes.exp.in.dd), alternative = "l")$p.value ## DD greater

png("boxplot_mesor_sd_ll_dd_4.png", width = 1800, height = 1800, res=300)
par(lwd=3)
boxplot(rowMeans(ld.genes.exp.in.ld), rowMeans(ld.genes.exp.in.ll),
        rowMeans(ld.genes.exp.in.dd), rowMeans(sd.genes.exp.in.sd),
        rowMeans(sd.genes.exp.in.ll), rowMeans(sd.genes.exp.in.dd),
        outline=F,col=c("blue", "lightblue", "darkblue","red","lightsalmon","firebrick4"),
        names=c("LD", "LD -> LL","LD -> DD", "SD","SD -> LL","SD -> DD"),cex.axis=2, las = 2)
dev.off()

# Mean for rhythmic non-circadian genes in each condition
complete.ld.ll.rhythmic.genes <- read.table(file="tables/ld.ll.rhythmic.table.tsv", sep="\t")$x
complete.ld.dd.rhythmic.genes <- read.table(file="tables/ld.dd.rhythmic.table.tsv", sep="\t")$x
complete.sd.ll.rhythmic.genes <- read.table(file="tables/sd.ll.rhythmic.table.tsv", sep="\t")$x
complete.sd.dd.rhythmic.genes <- read.table(file="tables/sd.dd.rhythmic.table.tsv", sep="\t")$x

non.circadian.rhythmic.genes.ld <- setdiff(complete.ld.rhythmic.genes, 
                                          unique(c(complete.ld.ll.rhythmic.genes,complete.ld.dd.rhythmic.genes)))

non.circadian.rhythmic.genes.sd <- setdiff(complete.sd.rhythmic.genes, 
                                             unique(c(complete.sd.ll.rhythmic.genes,complete.sd.dd.rhythmic.genes)))

ld.genes.exp.in.ld.non.circ <- genes.expression[non.circadian.rhythmic.genes.ld, ld.zt.i]
ld.genes.exp.in.ll.non.circ <- genes.expression[non.circadian.rhythmic.genes.ld, ll.ld.zt.i]
ld.genes.exp.in.dd.non.circ <- genes.expression[non.circadian.rhythmic.genes.ld, dd.ld.zt.i]

wilcox.test(x = rowMeans(ld.genes.exp.in.ld.non.circ), 
            y = rowMeans(ld.genes.exp.in.ll.non.circ), alternative = "g")$p.value # sig
wilcox.test(x = rowMeans(ld.genes.exp.in.ld.non.circ),
            y = rowMeans(ld.genes.exp.in.dd.non.circ), alternative = "l")$p.value # sig

sd.genes.exp.in.sd.non.circ <- genes.expression[non.circadian.rhythmic.genes.sd, sd.zt.i]
sd.genes.exp.in.ll.non.circ <- genes.expression[non.circadian.rhythmic.genes.sd, ll.sd.zt.i]
sd.genes.exp.in.dd.non.circ <- genes.expression[non.circadian.rhythmic.genes.sd, dd.sd.zt.i]

wilcox.test(x = rowMeans(sd.genes.exp.in.sd.non.circ), 
            y = rowMeans(sd.genes.exp.in.ll.non.circ), alternative = "g")$p.value # no sig
wilcox.test(x = rowMeans(sd.genes.exp.in.sd.non.circ),
            y = rowMeans(sd.genes.exp.in.dd.non.circ), alternative = "l")$p.value # DD greater


png("boxplot_expression_non_circadian_genes_by_cond.png", width = 1800, height = 1800, res=200)
par(lwd=3)
boxplot(rowMeans(ld.genes.exp.in.ld.non.circ), rowMeans(ld.genes.exp.in.ll.non.circ),
        rowMeans(ld.genes.exp.in.dd.non.circ), rowMeans(sd.genes.exp.in.sd.non.circ),
        rowMeans(sd.genes.exp.in.ll.non.circ), rowMeans(sd.genes.exp.in.dd.non.circ),
        outline=F,col=c("blue", "lightblue", "darkblue","red","lightsalmon","firebrick4"),
        names=c("LD", "LD -> LL","LD -> DD", "SD","SD -> LL","SD -> DD"),cex.axis=2, las = 2)
dev.off()


################# Effects of light on transitions
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


# Venn diagrams
library(VennDiagram)

png("response_to_light_ld_sd.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(ld.light.activated.dark.repressed), area2 = length(sd.light.activated.dark.repressed),
                   cross.area = length(intersect(ld.light.activated.dark.repressed,sd.light.activated.dark.repressed)),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()


png("response_to_dark_ld_sd.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(ld.light.repressed.dark.activated), area2 = length(sd.light.repressed.dark.activated),
                   cross.area = length(intersect(ld.light.repressed.dark.activated,sd.light.repressed.dark.activated)),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

png("response_to_dark_ld_response_to_light_sd.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(ld.light.repressed.dark.activated), area2 = length(sd.light.activated.dark.repressed),
                   cross.area = length(intersect(ld.light.repressed.dark.activated,sd.light.activated.dark.repressed)),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

png("response_to_dark_sd_response_to_light_ld.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(ld.light.activated.dark.repressed), area2 = length(sd.light.repressed.dark.activated),
                   cross.area = length(intersect(ld.light.activated.dark.repressed,sd.light.repressed.dark.activated)),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()



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
source("../../../GeneSCARAB/PSEA_functions_package.R")
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



# Histograms for non-circadian genes (only the ones in intersect, to explain the diel histograms)
library(data.table)
circacompare.ld.sd <- as.data.frame(fread("tables/circacompare_intersection.tsv", sep="\t"))
circacompare.ld.sd.names <- circacompare.ld.sd$V1
circacompare.ld.sd <- circacompare.ld.sd[,-1]
rownames(circacompare.ld.sd) <- circacompare.ld.sd.names

rhythmic.only.ld
rhythmic.only.sd

png("hist_phase_ld_intersection_non_circadian.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[rhythmic.only.ld,"ld peak time hours" ],col="blue",
     ylim=c(0,800),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

png("hist_phase_sd_intersection_non_circadian.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[rhythmic.only.sd,"sd peak time hours" ],col="red",
     ylim=c(0,1000),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

# Histograms for circadian genes in light/dark cycles (only the ones in intersect, to explain the diel histograms)
png("hist_phase_ld_ll_in_cycle.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[unique(c(rhythmic.genes.ld.ll.24, 
                                 rhythmic.genes.ld.ll.12)),"ld peak time hours" ],col="blue",
     ylim=c(0,300),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

png("hist_phase_ld_dd_in_cycle.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[unique(c(rhythmic.genes.ld.dd.24, 
                                 rhythmic.genes.ld.dd.12)),"ld peak time hours" ],col="blue",
     ylim=c(0,300),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

png("hist_phase_sd_ll_in_cycle.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[unique(c(rhythmic.genes.sd.ll.24, 
                                 rhythmic.genes.sd.ll.12)),"sd peak time hours" ],col="red",
     ylim=c(0,300),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

png("hist_phase_sd_dd_in_cycle.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[unique(c(rhythmic.genes.sd.dd.24, 
                                 rhythmic.genes.sd.dd.12)),"ld peak time hours" ],col="red",
     ylim=c(0,300),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()



# Histograms for bona fide genes
bona_fide <- intersect(rhythmic.ld.ll.dd,rhythmic.sd.ll.dd)

png("hist_bona_fide_ld.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[bona_fide,"ld peak time hours" ],col="blue",
     ylim=c(0,20),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="", breaks = seq(0,24,2))
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

png("hist_bona_fide_ld_ll.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(results.ld.ll.24[bona_fide,"peak_time_hours" ],col="lightblue",
     ylim=c(0,20),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="", breaks = seq(0,24,2))
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

png("hist_bona_fide_ld_dd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(results.ld.dd.24[bona_fide,"peak_time_hours" ],col="darkblue",
     ylim=c(0,20),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="", breaks = seq(0,24,2))
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

png("hist_bona_fide_sd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[bona_fide,"sd peak time hours" ],col="red",
     ylim=c(0,20),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="", breaks = seq(0,24,2))
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

png("hist_bona_fide_sd_ll.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(results.sd.ll.24[bona_fide,"peak_time_hours" ],col="lightsalmon",
     ylim=c(0,20),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="", breaks = seq(0,24,2))
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

png("hist_bona_fide_sd_dd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(results.sd.dd.24[bona_fide,"peak_time_hours" ],col="firebrick4",
     ylim=c(0,20),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="", breaks = seq(0,24,2))
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()


# Histograms for light and dark-regulated genes
ld.light.activated.dark.repressed <- read.table("tables/LD_light_activated_dark_repressed.tsv", sep="\t", 
                                                header = T)$x
ld.light.repressed.dark.activated <- read.table("tables/LD_light_repressed_dark_activated.tsv", sep="\t", 
                                                header = T)$x
sd.light.activated.dark.repressed <- read.table("tables/SD_light_activated_dark_repressed.tsv", sep="\t", 
                                                header = T)$x
sd.light.repressed.dark.activated <- read.table("tables/SD_light_repressed_dark_activated.tsv", sep="\t", 
                                                header = T)$x



png("hist_light_act_dark_rep_ld.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[ld.light.activated.dark.repressed,"ld peak time hours" ],col="blue",
     ylim=c(0,250),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="", breaks = seq(0,24,2))
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

png("hist_light_rep_dark_act_ld.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[ld.light.repressed.dark.activated,"ld peak time hours" ],col="blue",
     ylim=c(0,250),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="", breaks = seq(0,24,2))
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

png("hist_light_act_dark_rep_sd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[sd.light.activated.dark.repressed,"sd peak time hours" ],col="red",
     ylim=c(0,250),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="", breaks = seq(0,24,2))
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()

png("hist_light_rep_dark_act_sd.png", width = 1500, height = 1500, res=300)
par(lwd=3)
hist(circacompare.ld.sd[sd.light.repressed.dark.activated,"sd peak time hours" ],col="red",
     ylim=c(0,250),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="", breaks = seq(0,24,2))
axis(side = 2,lwd=3,cex.axis=1.3,las=2)
#abline(v=16, col="black", lwd=5, lty = 2)
dev.off()





###### Functional analysis
# Using GeneSCARAB
# First, load annotation and phases information for LD and SD
source("../../../GeneSCARAB/PSEA_functions_package.R")
library(circular)
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

# HR positive results may be due to low number of associated genes
sum(go.circa.table.ld[go.circa.table.ld$n > 3,]$hr_p_value_adj < 0.05 | go.circa.table.ld[go.circa.table.ld$n > 3,]$rao_p_value_adj < 0.05, na.rm = T) #18
sum(go.circa.table.sd[go.circa.table.sd$n > 3,]$hr_p_value_adj < 0.05 | go.circa.table.sd[go.circa.table.sd$n > 3,]$rao_p_value_adj < 0.05, na.rm = T) #22

go.circa.table.sd[which(go.circa.table.sd$rao_p_value_adj < 0.05),]
plot(circular(phases.list.sd$`GO:0051253`$phase)) # avoid night


unimod_ld <- rownames(go.circa.table.ld)[which(go.circa.table.ld$rayleigh_p_value_adj < 0.05 | go.circa.table.ld$kuiper_p_value_adj < 0.05)]
unimod_sd <- rownames(go.circa.table.sd)[which(go.circa.table.sd$rayleigh_p_value_adj < 0.05 | go.circa.table.sd$kuiper_p_value_adj < 0.05)]

multimod_ld <- rownames(go.circa.table.ld[go.circa.table.ld$n > 3,])[which(go.circa.table.ld[go.circa.table.ld$n > 3,]$hr_p_value_adj < 0.05 | go.circa.table.ld[go.circa.table.ld$n > 3,]$rao_p_value_adj < 0.05)]
multimod_sd <- rownames(go.circa.table.sd[go.circa.table.sd$n > 3,])[which(go.circa.table.sd[go.circa.table.sd$n > 3,]$hr_p_value_adj < 0.05 | go.circa.table.sd[go.circa.table.sd$n > 3,]$rao_p_value_adj < 0.05)]


write.table(unimod_ld, "tables/gos_ld_unimodal.tsv", sep = "\t", row.names = F)
write.table(unimod_sd, "tables/gos_sd_unimodal.tsv", sep = "\t", row.names = F)
write.table(multimod_ld, "tables/gos_ld_multimodal.tsv", sep = "\t", row.names = F)
write.table(multimod_sd, "tables/gos_sd_multimodal.tsv", sep = "\t", row.names = F)


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

# Some examples
# DNA replication
plot(circular(phases.list.sd$`GO:0006261`$phase))
plot(circular(phases.list.ld$`GO:0006261`$phase))


# Check if there exist modality changes between GOs
intersect(unimod_ld, multimod_sd)
intersect(unimod_sd, multimod_ld)
plot(circular(phases.list.ld$`GO:0005247`$phase))

# For the GO terms not following an uniform distribution, test them against the general distribution of phases
# per photoperiod
non_uniform_gos_ld <- unique(c(unimod_ld, multimod_ld))
non_uniform_gos_ld <- non_uniform_gos_ld[non_uniform_gos_ld %in% names(phases.list.ld.clean)]
phases.list.ld.non.uni <- phases.list.ld.clean[non_uniform_gos_ld]

gen.dist.ld <- test_against_gen_dist(phase.list = phases.list.ld.non.uni, total.phase.table = total_phases_table_ld)
write.table(gen.dist.ld,file="tables/go_general_dist_ld_filtered.tsv", sep="\t")

non_uniform_gos_sd <- unique(c(unimod_sd, multimod_sd))
non_uniform_gos_sd <- non_uniform_gos_sd[non_uniform_gos_sd %in% names(phases.list.sd.clean)]
phases.list.sd.non.uni <- phases.list.sd.clean[non_uniform_gos_sd]
gen.dist.sd <- test_against_gen_dist(phase.list = phases.list.sd.non.uni, total.phase.table = total_phases_table_sd)
write.table(gen.dist.sd,file="tables/go_general_dist_sd_filtered.tsv", sep="\t")

gen.dist.ld <- read.table("tables/go_general_dist_ld_filtered.tsv", sep="\t", header = T)
gen.dist.sd <- read.table("tables/go_general_dist_sd_filtered.tsv", sep="\t", header = T)

sum(gen.dist.ld$gen_dist_p_value_adj < 0.05) # 451 non-uniform GO terms present a different distribution than that of global
# mRNA accumulation in LD
sum(gen.dist.sd$gen_dist_p_value_adj < 0.05) # 256 non-uniform GO terms present a different distribution than that of global
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

write.table(go_dist_gen_ld, "tables/go_diff_general_ld.tsv", sep = "\t", row.names = F)
write.table(go_dist_gen_sd, "tables/go_diff_general_sd.tsv", sep = "\t", row.names = F)

# Add description of non-experimental distributed GOs
GO_complete_list <- as.list(GOTERM)

go_dist_gen_ld <- go_dist_gen_ld[go_dist_gen_ld %in% names(GO_complete_list)]
go_dist_gen_ld_description <- sapply(go_dist_gen_ld, function(x) GO_complete_list[[x]]@Term)
go_dist_gen_ld_df <- data.frame(go=names(go_dist_gen_ld_description), description=go_dist_gen_ld_description)
write.table(go_dist_gen_ld_df,file="tables/go_general_dist_ld_filt_description.tsv", sep="\t", col.names = T)

go_dist_gen_sd <- go_dist_gen_sd[go_dist_gen_sd %in% names(GO_complete_list)]
go_dist_gen_sd_description <- sapply(go_dist_gen_sd, function(x) GO_complete_list[[x]]@Term)
go_dist_gen_sd_df <- data.frame(go=names(go_dist_gen_sd_description), description=go_dist_gen_sd_description)
write.table(go_dist_gen_sd_df,file="tables/go_general_dist_sd_filt_description.tsv", sep="\t", col.names = T)

shared_non_exp_dist <- intersect(go_dist_gen_ld, go_dist_gen_sd)
go_dist_gen_ld_sd_description <- sapply(shared_non_exp_dist, function(x) GO_complete_list[[x]]@Term)
go_dist_gen_sd_df <- data.frame(go=names(go_dist_gen_ld_sd_description), description=go_dist_gen_ld_sd_description)
write.table(go_dist_gen_sd_df,file="tables/go_general_dist_ld_sd_filt_description.tsv", sep="\t", col.names = T, row.names = F)


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

write.table(go_diff_ld_sd_table,file="tables/go_diff_ld_sd_table_filtered", sep="\t")

go_diff_ld_sd_table <- read.table("tables/go_diff_ld_sd_table_filtered", sep="\t", header = T)

sum(go_diff_ld_sd_table$gen_dist_p_value_adj < 0.05) # 525 change their distribution between LD and SD
change_gos <- rownames(go_diff_ld_sd_table)[which(go_diff_ld_sd_table$gen_dist_p_value_adj < 0.05)]

# Difference in medians for changing GOs
diff_gos_medians <- go.circa.table.ld[change_gos,"median"] - go.circa.table.sd[change_gos,"median"]

sd_anticipated <- change_gos[which((0 < diff_gos_medians & diff_gos_medians < 12) | (diff_gos_medians < -12))]
# 524 are anticipated in SD, only 1 in LD
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

change_gos <- change_gos[change_gos %in% names(GO_complete_list)]
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

write.table(change_gos_bp_table,file="tables/anticipated_gos_sd_bp_filt.tsv", sep="\t", col.names = T)
write.table(change_gos_mf_table,file="tables/anticipated_gos_sd_mf_filt.tsv", sep="\t", col.names = T)
write.table(change_gos_cc_table,file="tables/anticipated_gos_sd_cc_filt.tsv", sep="\t", col.names = T)


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
kegg_plants <- read.csv("../../../GeneSCARAB/pathways_plants.txt", sep = "\t", header = T)$x
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

# HR positive results may be due to low number of associated genes
sum(kegg.circa.table.ld[kegg.circa.table.ld$n > 2,]$hr_p_value_adj < 0.05, na.rm = T) #0
sum(kegg.circa.table.sd[kegg.circa.table.sd$n > 2,]$hr_p_value_adj < 0.05, na.rm = T) #0

unimod_ld <- rownames(kegg.circa.table.ld)[which(kegg.circa.table.ld$rayleigh_p_value_adj < 0.05 | kegg.circa.table.ld$kuiper_p_value_adj < 0.05)]
unimod_sd <- rownames(kegg.circa.table.sd)[which(kegg.circa.table.sd$rayleigh_p_value_adj < 0.05 | kegg.circa.table.sd$kuiper_p_value_adj < 0.05)]

multimod_ld <- rownames(kegg.circa.table.ld)[which(kegg.circa.table.ld$hr_p_value_adj < 0.05 | kegg.circa.table.ld$rao_p_value_adj < 0.05)]
multimod_sd <- rownames(kegg.circa.table.sd)[which(kegg.circa.table.sd$hr_p_value_adj < 0.05 | kegg.circa.table.sd$rao_p_value_adj < 0.05)]

kegg.circa.table.ld[multimod_ld,] # only one data point
kegg.circa.table.sd[multimod_sd,] # only one data point
# No significant multimodal pathways

write.table(unimod_ld, "tables/kegg_ld_unimodal.tsv", sep = "\t", row.names = F)
write.table(unimod_sd, "tables/kegg_sd_unimodal.tsv", sep = "\t", row.names = F)

library(VennDiagram)
png("ld_sd_unimod_kegg.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(unimod_ld), area2 = length(unimod_sd),
                   cross.area = length(intersect(unimod_ld,unimod_sd)),
                   lwd = 3,catekeggry = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# png("ld_sd_multimod_kegg.png", res=300, width = 1500, height = 1500)
# grid.newpage()
# draw.pairwise.venn(area1 = length(multimod_ld), area2 = length(multimod_sd),
#                    cross.area = length(intersect(multimod_ld,multimod_sd)),
#                    lwd = 3,catekeggry = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
#                    alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
# dev.off()

# For the kegg terms not following an uniform distribution, test them against the general distribution of phases
# per photoperiod
non_uniform_keggs_ld <- unique(c(unimod_ld)) #no significant
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

write.table(kegg_dist_gen_ld, "tables/kegg_general_dist_ld.tsv", sep = "\t", row.names = F)
write.table(kegg_dist_gen_sd, "tables/kegg_general_dist_sd.tsv", sep = "\t", row.names = F)

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

write.table(kegg_diff_ld_sd_table,file="tables/kegg_diff_ld_sd_table.tsv", sep="\t")

kegg_diff_ld_sd_table <- read.table("tables/kegg_diff_ld_sd_table.tsv", sep="\t", header = T)

sum(kegg_diff_ld_sd_table$gen_dist_p_value_adj < 0.05) # 27 change their distribution between LD and SD
change_keggs <- rownames(kegg_diff_ld_sd_table)[which(kegg_diff_ld_sd_table$gen_dist_p_value_adj < 0.05)]

# Difference in medians for changing keggs
diff_keggs_medians <- kegg.circa.table.ld[change_keggs,"median"] - kegg.circa.table.sd[change_keggs,"median"]

sd_anticipated <- change_keggs[which((0 < diff_keggs_medians & diff_keggs_medians < 12) | (diff_keggs_medians < -12))]
length(sd_anticipated)
# All 27 changing keggs are anticipated in SD

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


# Select some pathways to plot, non-uniform and SD-anticipated
# Photosynthesis, Photosynthesis - antenna proteins, Carbon fixation by Calvin cycle,
# Spliceosome, RNA degradation, Motor proteins, Endocytosis, Ribosome,
# Protein processing in endoplasmic reticulum, Proteasome, Plant hormone signal transduction

# Same order:
# map00195, map00196, map00710, map03040, map03018, map04814, map04144, map03010, 
# map04141, map03050, map04075
library(circlize)
library(ComplexHeatmap)
library(MetBrewer)

selected_kegg_ld <- kegg.circa.table.ld[c("map00195", "map00196", "map00710", "map03040",
                                          "map03018", "map04814", "map04144", "map03010", 
                                          "map04141", "map03050", "map04075"),]

selected_kegg_sd <- kegg.circa.table.sd[c("map00195", "map00196", "map00710", "map03040",
                                          "map03018", "map04814", "map04144", "map03010", 
                                          "map04141", "map03050", "map04075"),]

circular_boxplot(selected_kegg_ld, "selected_keggs_ld.png", color.palette = "Tam",
                 tr.height = 0.075)
circular_boxplot(selected_kegg_sd, "selected_keggs_sd.png", color.palette = "Tam",
                 tr.height = 0.075)




##################  Functional enrichment of free-running conditions  ##########################
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
library(enrichplot)
library("org.Mpolymorpha71.eg.db")
universe_set <- unique(c(complete.ld.rhythmic.genes, complete.sd.rhythmic.genes))

ego.inter.ll <- enrichGO(gene = inter_ll,
                         OrgDb = org.Mpolymorpha71.eg.db,
                         #universe = universe_set,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = F,
                         keyType = "GID")

barplot(ego.inter.ll, showCategory = 30)

png("ego_inter_ll_no_universe.png", res=200, width = 1700, height = 1700)
emapplot(pairwise_termsim(ego.inter.ll))
dev.off()

treeplot(pairwise_termsim(ego.inter.ll))

ego.inter.ll.df <- as.data.frame(ego.inter.ll)

ego.inter.dd <- enrichGO(gene = inter_dd,
                         OrgDb         = org.Mpolymorpha71.eg.db,
                         #universe = universe_set,
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = F,
                         keyType = "GID")

ego.union.ll <- enrichGO(gene = union_ll,
                         OrgDb         = org.Mpolymorpha71.eg.db,
                         #universe = universe_set,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = F,
                         keyType = "GID")

barplot(ego.union.ll, showCategory = 20)
png("ego_union_ll.png", res=200, width = 1700, height = 1700)
emapplot(pairwise_termsim(ego.union.ll))
dev.off()
treeplot(pairwise_termsim(ego.union.ll))
ego.union.ll.df <- as.data.frame(ego.union.ll)

ego.union.dd <- enrichGO(gene = union_dd,
                         OrgDb         = org.Mpolymorpha71.eg.db,
                         #universe = universe_set,
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = F,
                         keyType = "GID")

png("ego_union_dd_no_universe.png", res=200, width = 1800, height = 1200)
barplot(ego.union.dd, showCategory = 20)
dev.off()

# Functional annotation of bona_fide genes
complete.sd.dd.rhythmic.genes <- read.table("tables/sd.dd.rhythmic.table.tsv", sep="\t")$x
complete.ld.dd.rhythmic.genes <- read.table("tables/ld.dd.rhythmic.table.tsv", sep="\t")$x
complete.sd.ll.rhythmic.genes <- read.table("tables/sd.ll.rhythmic.table.tsv", sep="\t")$x
complete.ld.ll.rhythmic.genes <- read.table("tables/ld.ll.rhythmic.table.tsv", sep="\t")$x

complete.sd.ll.dd <- intersect(complete.sd.dd.rhythmic.genes, complete.sd.ll.rhythmic.genes)
complete.ld.ll.dd <- intersect(complete.ld.dd.rhythmic.genes, complete.ld.ll.rhythmic.genes)
bona_fide <- intersect(complete.ld.ll.dd, complete.sd.ll.dd)

ego.bona.fide <- enrichGO(gene = bona_fide,
                         OrgDb = org.Mpolymorpha71.eg.db,
                         #universe = universe_set,
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = F,
                         keyType = "GID")

barplot(ego.bona.fide, showCategory = 20)

png("ego_bona_fide.png", res=300, width = 1700, height = 1700)
emapplot(pairwise_termsim(ego.bona.fide))
dev.off()

png("ego_bona_fide_tree.png", res=300, width = 3500, height = 1700)
treeplot(pairwise_termsim(ego.bona.fide), cluster.params = list(n = 2))
dev.off()

# LL LD-specific, not bona fide
ll.ld.specific <- setdiff(complete.ld.ll.rhythmic.genes, complete.sd.ll.rhythmic.genes)
ego.ll.ld.specific <- enrichGO(gene = ll.ld.specific,
                                 OrgDb = org.Mpolymorpha71.eg.db,
                                 universe = universe_set,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05,
                                 readable      = F,
                                 keyType = "GID")

barplot(ego.ll.ld.specific, showCategory = 20)

png("ego_ll_ld_specific.png", res=200, width = 1700, height = 1700)
emapplot(pairwise_termsim(ego.ll.ld.specific))
dev.off()

write.table(as.data.frame(ego.ll.ld.specific), "tables/ld_ll_specific_gos.tsv", sep = "\t", row.names = F)

# LL SD-specific, not bona fide
ll.sd.specific <- setdiff(complete.sd.ll.rhythmic.genes, complete.ld.ll.rhythmic.genes)
ego.ll.sd.specific <- enrichGO(gene = ll.sd.specific,
                               OrgDb = org.Mpolymorpha71.eg.db,
                               universe = universe_set,
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05,
                               readable      = F,
                               keyType = "GID")

barplot(ego.ll.sd.specific, showCategory = 20)

png("ego_ll_sd_specific.png", res=200, width = 1700, height = 1700)
emapplot(pairwise_termsim(ego.ll.sd.specific))
dev.off()

write.table(as.data.frame(ego.ll.sd.specific), "tables/sd_ll_specific_gos.tsv", sep = "\t", row.names = F)

# DD LD-specific, not bona fide
dd.ld.specific <- setdiff(complete.ld.dd.rhythmic.genes, complete.sd.dd.rhythmic.genes)
ego.dd.ld.specific <- enrichGO(gene = dd.ld.specific,
                               OrgDb = org.Mpolymorpha71.eg.db,
                               universe = universe_set,
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05,
                               readable      = F,
                               keyType = "GID")

barplot(ego.dd.ld.specific, showCategory = 20)

# DD SD-specific, not bona fide
dd.sd.specific <- setdiff(complete.sd.dd.rhythmic.genes, complete.ld.dd.rhythmic.genes)
ego.dd.sd.specific <- enrichGO(gene = dd.sd.specific,
                               OrgDb = org.Mpolymorpha71.eg.db,
                               universe = universe_set,
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05,
                               readable      = F,
                               keyType = "GID")

barplot(ego.dd.sd.specific, showCategory = 20)



# LL-DD LD-specific, not bona fide
free.ld.specific <- setdiff(complete.ld.ll.dd, complete.sd.ll.dd)
ego.free.ld.specific <- enrichGO(gene = free.ld.specific,
                          OrgDb = org.Mpolymorpha71.eg.db,
                          universe = universe_set,
                          ont           = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = F,
                          keyType = "GID")

barplot(ego.free.ld.specific, showCategory = 20)

# LL-DD SD-specific, not bona fide
free.sd.specific <- setdiff(complete.sd.ll.dd, complete.ld.ll.dd)
ego.free.sd.specific <- enrichGO(gene = free.sd.specific,
                                 OrgDb = org.Mpolymorpha71.eg.db,
                                 universe = universe_set,
                                 ont           = "ALL",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05,
                                 readable      = F,
                                 keyType = "GID")

barplot(ego.free.sd.specific, showCategory = 20)





##########################     TFs ANALYSIS    #############################

tf_list <- unique(c("Mp1g00510","Mp1g01010","Mp1g01110","Mp1g01160","Mp1g01900","Mp1g01960","Mp1g02320","Mp1g02860",
         "Mp1g02890","Mp1g03040","Mp1g03100","Mp1g03200","Mp1g03580","Mp1g03950","Mp1g04480","Mp1g04550",
         "Mp1g05310","Mp1g05400","Mp1g05630","Mp1g06180","Mp1g07040","Mp1g07070","Mp1g08640","Mp1g08660",
         "Mp1g08890","Mp1g08960","Mp1g09040","Mp1g09420","Mp1g09500","Mp1g10020","Mp1g10030","Mp1g10300",
         "Mp1g10310","Mp1g10330","Mp1g10440","Mp1g11660","Mp1g11830","Mp1g12700","Mp1g12750","Mp1g13010",
         "Mp1g13260","Mp1g13640","Mp1g13740","Mp1g14500","Mp1g15660","Mp1g16790","Mp1g16800","Mp1g16880",
         "Mp1g17210","Mp1g18540","Mp1g19210","Mp1g19590","Mp1g19730","Mp1g20040","Mp1g20490","Mp1g21450",
         "Mp1g22880","Mp1g23310","Mp1g24140","Mp1g24520","Mp1g24950","Mp1g25020","Mp1g25140","Mp1g25720",
         "Mp1g26740","Mp1g27060","Mp1g28480","Mp1g28490","Mp1g29240","Mp2g00500","Mp2g00510","Mp2g00890",
         "Mp2g00910","Mp2g00930","Mp2g02150","Mp2g02230","Mp2g02280","Mp2g02990","Mp2g03850","Mp2g04160",
         "Mp2g04180","Mp2g04190","Mp2g04200","Mp2g05060","Mp2g05070","Mp2g05160","Mp2g06450","Mp2g06460",
         "Mp2g07110","Mp2g07130","Mp2g07150","Mp2g07170","Mp2g07720","Mp2g08260","Mp2g08360","Mp2g08790",
         "Mp2g09320","Mp2g10620","Mp2g11080","Mp2g11220","Mp2g11330","Mp2g11520","Mp2g12750","Mp2g13020",
         "Mp2g13330","Mp2g13910","Mp2g14310","Mp2g14450","Mp2g16340","Mp2g16350","Mp2g19970","Mp2g20695",
         "Mp2g20790","Mp2g20960","Mp2g20990","Mp2g22110","Mp2g22250","Mp2g22560","Mp2g22820","Mp2g23000",
         "Mp2g23600","Mp2g24110","Mp2g24200","Mp2g25260","Mp2g25420","Mp2g25660","Mp2g25840","Mp2g26400",
         "Mp3g00010","Mp3g00360","Mp3g00600","Mp3g01640","Mp3g02320","Mp3g02770","Mp3g02830","Mp3g04030",
         "Mp3g04360","Mp3g04970","Mp3g05000","Mp3g05910","Mp3g06860","Mp3g07510","Mp3g07610","Mp3g09110",
         "Mp3g09350","Mp3g09590","Mp3g10500","Mp3g10830","Mp3g10970","Mp3g11570","Mp3g11900","Mp3g12060",
         "Mp3g12980","Mp3g16240","Mp3g17260","Mp3g17350","Mp3g17660","Mp3g17930","Mp3g18590","Mp3g19160",
         "Mp3g19510","Mp3g19580","Mp3g19670","Mp3g20870","Mp3g20900","Mp3g20980","Mp3g20990","Mp3g21340",
         "Mp3g21490","Mp3g23170","Mp3g23300","Mp3g23720","Mp3g24320","Mp3g24780","Mp3g25200","Mp4g00180",
         "Mp4g00200","Mp4g00280","Mp4g00380","Mp4g01170","Mp4g01540","Mp4g01560","Mp4g02060","Mp4g02620",
         "Mp4g04750","Mp4g04760","Mp4g04910","Mp4g04920","Mp4g06170","Mp4g06490","Mp4g06740","Mp4g07330",
         "Mp4g07910","Mp4g08690","Mp4g08700","Mp4g09650","Mp4g09730","Mp4g09860","Mp4g09880","Mp4g10240",
         "Mp4g10500","Mp4g10710","Mp4g11030","Mp4g11380","Mp4g11760","Mp4g11820","Mp4g11860","Mp4g11910",
         "Mp4g12230","Mp4g12420","Mp4g12490","Mp4g13360","Mp4g13860","Mp4g13880","Mp4g14520","Mp4g14530",
         "Mp4g15170","Mp4g15240","Mp4g15460","Mp4g15560","Mp4g17430","Mp4g18230","Mp4g19650","Mp4g20600",
         "Mp4g21220","Mp4g21300","Mp4g21790","Mp4g22280","Mp4g22290","Mp4g22310","Mp4g22430","Mp4g22890",
         "Mp5g00110","Mp5g01050","Mp5g01530","Mp5g01600","Mp5g02590","Mp5g03020","Mp5g05560","Mp5g05940",
         "Mp5g06970","Mp5g07255","Mp5g07890","Mp5g08310","Mp5g08410","Mp5g08810","Mp5g09130","Mp5g09220",
         "Mp5g09660","Mp5g09710","Mp5g10280","Mp5g10960","Mp5g11060","Mp5g11100","Mp5g11460","Mp5g11830",
         "Mp5g12220","Mp5g12480","Mp5g13030","Mp5g13040","Mp5g13060","Mp5g13130","Mp5g13160","Mp5g14610",
         "Mp5g14820","Mp5g17820","Mp5g18910","Mp5g19050","Mp5g20280","Mp5g20530","Mp5g20660","Mp5g21080",
         "Mp5g21990","Mp5g22160","Mp5g22280","Mp5g23860","Mp6g00540","Mp6g02620","Mp6g02670","Mp6g02950",
         "Mp6g03670","Mp6g03920","Mp6g04430","Mp6g04650","Mp6g04830","Mp6g04880","Mp6g05720","Mp6g05810",
         "Mp6g06850","Mp6g06870","Mp6g07410","Mp6g07440","Mp6g08290","Mp6g08690","Mp6g09270","Mp6g09280",
         "Mp6g10210","Mp6g10570","Mp6g10690","Mp6g11770","Mp6g13820","Mp6g15470","Mp6g16170","Mp6g16760",
         "Mp6g16800","Mp6g17420","Mp6g18480","Mp6g18650","Mp6g18860","Mp6g20920","Mp6g21430","Mp6g21470",
         "Mp7g00210","Mp7g00650","Mp7g00860","Mp7g01800","Mp7g02220","Mp7g02640","Mp7g03150","Mp7g03490",
         "Mp7g03700","Mp7g06280","Mp7g06290","Mp7g06550","Mp7g07550","Mp7g09260","Mp7g09350","Mp7g09490",
         "Mp7g09710","Mp7g09740","Mp7g11330","Mp7g11500","Mp7g12630","Mp7g12800","Mp7g13270","Mp7g13760",
         "Mp7g14310","Mp7g14810","Mp7g15840","Mp7g16300","Mp7g16910","Mp7g17020","Mp7g17200","Mp7g17250",
         "Mp7g18530","Mp7g18640","Mp7g18690","Mp7g19570","Mp8g01190","Mp8g01770","Mp8g02290","Mp8g03820",
         "Mp8g03980","Mp8g04010","Mp8g04130","Mp8g04330","Mp8g04400","Mp8g06180","Mp8g06230","Mp8g06450",
         "Mp8g07390","Mp8g08100","Mp8g08190","Mp8g09250","Mp8g09500","Mp8g10640","Mp8g11450","Mp8g11560",
         "Mp8g11850","Mp8g11870","Mp8g12000","Mp8g12100","Mp8g12510","Mp8g12800","Mp8g13205","Mp8g13470",
         "Mp8g14220","Mp8g15500","Mp8g16170","Mp8g16900","Mp8g18080","Mp8g18310","Mp8g18820","MpUg00060",
         "MpUg00230","MpUg00270","MpUg00290","MpUg00370","MpVg00300","MpVg00340","MpVg00350","MpVg00470"))


# Rhythmic TFs in LD
complete.ld.rhythmic.genes <- read.table("tables/ld.rhythmic.table.tsv", header = T)$V1
complete.sd.rhythmic.genes <- read.table("tables/sd.rhythmic.table.tsv", header = T)$V1

ld_tfs <- tf_list[tf_list %in% complete.ld.rhythmic.genes]
sd_tfs <- tf_list[tf_list %in% complete.sd.rhythmic.genes]

complete.sd.dd.rhythmic.genes <- read.table("tables/sd.dd.rhythmic.table.tsv", sep="\t")$x
complete.ld.dd.rhythmic.genes <- read.table("tables/ld.dd.rhythmic.table.tsv", sep="\t")$x
complete.sd.ll.rhythmic.genes <- read.table("tables/sd.ll.rhythmic.table.tsv", sep="\t")$x
complete.ld.ll.rhythmic.genes <- read.table("tables/ld.ll.rhythmic.table.tsv", sep="\t")$x


ld_ll_tfs <- tf_list[tf_list %in% complete.ld.ll.rhythmic.genes]
ld_dd_tfs <- tf_list[tf_list %in% complete.ld.dd.rhythmic.genes]
sd_ll_tfs <- tf_list[tf_list %in% complete.sd.ll.rhythmic.genes]
sd_dd_tfs <- tf_list[tf_list %in% complete.sd.dd.rhythmic.genes]

ll_inter_tfs <- intersect(ld_ll_tfs, sd_ll_tfs)
dd_inter_tfs <- intersect(ld_dd_tfs, sd_dd_tfs)
bona_fide_tfs <- intersect(ll_inter_tfs, dd_inter_tfs)


# Venn diagrams for TFs
library(VennDiagram)
png("tfs_ll_inter.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(ld_ll_tfs), area2 = length(sd_ll_tfs),
                   cross.area = length(intersect(ld_ll_tfs,sd_ll_tfs)),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

png("tfs_dd_inter.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(ld_dd_tfs), area2 = length(sd_dd_tfs),
                   cross.area = length(intersect(ld_dd_tfs,sd_dd_tfs)),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# Correlation of free running-rhythmic TFs with free running-rhythhic genes
source("functions_circa_marchantia.R")

circa_sd_ll <- read.table("tables/estimates_sd_ll.tsv", sep="\t")
circa_sd_dd <- read.table("tables/estimates_sd_dd.tsv", sep="\t")
circa_ld_ll <- read.table("tables/estimates_ld_ll.tsv", sep="\t")
circa_ld_dd <- read.table("tables/estimates_ld_dd.tsv", sep="\t")


png("hy5.png", width = 3100, height = 2500, res = 300)
par(mfrow = c(2, 2))
plot.ld.ll.wave("Mp1g16800", "MpHY5", genes.expression, period=1, 
                wave_parameters=circa_ld_ll["Mp1g16800",])
plot.ld.dd.wave("Mp1g16800", "MpHY5", genes.expression, period=1, 
                wave_parameters=circa_ld_dd["Mp1g16800",])
plot.sd.ll.wave("Mp1g16800", "MpHY5", genes.expression, period=1, 
                wave_parameters=circa_sd_ll["Mp1g16800",])
plot.sd.dd.wave("Mp1g16800", "MpHY5", genes.expression, period=1, 
                wave_parameters=circa_sd_dd["Mp1g16800",])
dev.off()

png("trihelix11.png", width = 3100, height = 2500, res = 300)
par(mfrow = c(2, 2))
plot.ld.ll.wave("Mp2g08360", "MpTRIHELIX11", genes.expression, period=1, 
                wave_parameters=circa_ld_ll["Mp2g08360",])
plot.ld.dd.wave("Mp2g08360", "MpTRIHELIX11", genes.expression, period=1, 
                wave_parameters=circa_ld_dd["Mp2g08360",])
plot.sd.ll.wave("Mp2g08360", "MpTRIHELIX11", genes.expression, period=1, 
                wave_parameters=circa_sd_ll["Mp2g08360",])
plot.sd.dd.wave("Mp2g08360", "MpTRIHELIX11", genes.expression, period=1, 
                wave_parameters=circa_sd_dd["Mp2g08360",])
dev.off()

png("trihelix31.png", width = 3100, height = 2500, res = 300)
par(mfrow = c(2, 2))
plot.ld.ll.wave("Mp4g09730", "MpTRIHELIX31", genes.expression, period=1, 
                wave_parameters=circa_ld_ll["Mp4g09730",])
plot.ld.dd.wave("Mp4g09730", "MpTRIHELIX31", genes.expression, period=1, 
                wave_parameters=circa_ld_dd["Mp4g09730",])
plot.sd.ll.wave("Mp4g09730", "MpTRIHELIX31", genes.expression, period=1, 
                wave_parameters=circa_sd_ll["Mp4g09730",])
plot.sd.dd.wave("Mp4g09730", "MpTRIHELIX31", genes.expression, period=1, 
                wave_parameters=circa_sd_dd["Mp4g09730",])
dev.off()


circa_table_tf_inter <- read.table("tables/circacompare_intersection.tsv", header = T, sep = "\t")
cycle_inter_tfs_table <- circa_table_tf_inter[rownames(circa_table_tf_inter) %in% tf_list,]
write.table(cycle_inter_tfs_table, "tables/tfs_ld_sd_circacompare_inter.tsv", col.names = T, sep="\t")

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

png("ordered_tfs_complete.png", width = 480, height = 1440, res=100)
k_plot <- ggplot(k_table, aes(x = values, y = cat)) +
  #geom_line(aes(color = allowed)) +
  geom_line() +
  geom_point(aes(color = photo), size = 3) +
  scale_color_brewer(palette = "Set1", direction = -1) +
  theme_bw() +
  theme(legend.position = "bottom")
plot(k_plot)
dev.off()

# For TFs retaining rhythmicity under LL
circa_ld_ll_inter_tfs <- circa_ld_ll[ll_inter_tfs,]
circa_sd_ll_inter_tfs <- circa_sd_ll[ll_inter_tfs,]
rownames(circa_ld_ll_inter_tfs) == rownames(circa_sd_ll_inter_tfs)

k_names <- c(paste0(rownames(circa_ld_ll_inter_tfs),"_ld"),
             paste0(rownames(circa_ld_ll_inter_tfs),"_ld_ll"),
             paste0(rownames(circa_sd_ll_inter_tfs),"_sd"),
             paste0(rownames(circa_sd_ll_inter_tfs),"_sd_ll"))

k_values <- c(circa_ld_ll_inter_tfs$ld.peak.time.hours, 
              circa_ld_ll_inter_tfs$ll.peak.time.hours,
              circa_sd_ll_inter_tfs$sd.peak.time.hours, 
              circa_sd_ll_inter_tfs$ll.peak.time.hours)

k_table <- data.frame(names = k_names, values = k_values,
                      cat = factor(rownames(circa_ld_ll_inter_tfs), levels = rownames(circa_ld_ll_inter_tfs)), 
                      photo = rep(c("ld", "ll_ld", "sd", "sd_ll"), each= nrow(circa_ld_ll_inter_tfs)))
library(ggplot2)

png("ll_inter_tfs.png", width = 1440, height = 1400, res=300)
k_plot <- ggplot(k_table, aes(x = values, y = cat)) +
  #geom_line(aes(color = allowed)) +
  #geom_line() +
  geom_point(aes(color = photo), size = 3) +
  #scale_color_brewer(palette = "Set1", direction = -1) +
  scale_color_manual(values = c("royalblue3", "skyblue","tomato3", "salmon")) +
  theme_bw() +
  theme(legend.position = "bottom")
plot(k_plot)
dev.off()

# For TFs retaining rhythmicity under DD
circa_ld_dd_inter_tfs <- circa_ld_dd[dd_inter_tfs,]
circa_sd_dd_inter_tfs <- circa_sd_dd[dd_inter_tfs,]
rownames(circa_ld_dd_inter_tfs) == rownames(circa_sd_dd_inter_tfs)

k_names <- c(paste0(rownames(circa_ld_dd_inter_tfs),"_ld"),
             paste0(rownames(circa_ld_dd_inter_tfs),"_ld_dd"),
             paste0(rownames(circa_sd_dd_inter_tfs),"_sd"),
             paste0(rownames(circa_sd_dd_inter_tfs),"_sd_dd"))

k_values <- c(circa_ld_dd_inter_tfs$ld.peak.time.hours, 
              circa_ld_dd_inter_tfs$dd.peak.time.hours,
              circa_sd_dd_inter_tfs$sd.peak.time.hours, 
              circa_sd_dd_inter_tfs$dd.peak.time.hours)

k_table <- data.frame(names = k_names, values = k_values,
                      cat = factor(rownames(circa_ld_dd_inter_tfs), levels = rownames(circa_ld_dd_inter_tfs)), 
                      photo = rep(c("ld", "dd_ld", "sd", "sd_dd"), each= nrow(circa_ld_dd_inter_tfs)))
library(ggplot2)

png("dd_inter_tfs.png", width = 1440, height = 1400, res=300)
k_plot <- ggplot(k_table, aes(x = values, y = cat)) +
  #geom_line(aes(color = addowed)) +
  #geom_line() +
  geom_point(aes(color = photo), size = 3) +
  #scale_color_brewer(palette = "Set1", direction = -1) +
  scale_color_manual(values = c("slateblue4", "royalblue3", "tomato3", "red4")) +
  theme_bw() +
  theme(legend.position = "bottom")
plot(k_plot)
dev.off()

# Clock components
# TOC1
"Mp3g09460" %in% rownames(circa_ld_ll)
"Mp3g09460" %in% rownames(circa_ld_dd)
"Mp3g09460" %in% rownames(circa_sd_ll)
"Mp3g09460" %in% rownames(circa_sd_dd)

# PRR
"Mp2g16560" %in% rownames(circa_ld_ll)
"Mp2g16560" %in% rownames(circa_ld_dd)
"Mp2g16560" %in% rownames(circa_sd_ll)
"Mp2g16560" %in% rownames(circa_sd_dd)

# RVE
"Mp2g14310" %in% rownames(circa_ld_ll)
"Mp2g14310" %in% rownames(circa_ld_dd)
"Mp2g14310" %in% rownames(circa_sd_ll)
"Mp2g14310" %in% rownames(circa_sd_dd)

# DET1
"Mp8g00970" %in% rownames(circa_ld_ll)
"Mp8g00970" %in% rownames(circa_ld_dd)
"Mp8g00970" %in% rownames(circa_sd_ll)
"Mp8g00970" %in% rownames(circa_sd_dd)

png("det1.png", width = 3100, height = 2500, res = 300)
par(mfrow = c(3, 1))
plot.ld.ll.wave("Mp8g00970", "DET1", genes.expression, period=1, 
                wave_parameters=circa_ld_ll["Mp8g00970",])
plot.ld.dd.wave("Mp8g00970", "DET1", genes.expression, period=1, 
                wave_parameters=circa_ld_dd["Mp8g00970",])
plot.sd.ll.wave("Mp8g00970", "DET1", genes.expression, period=1, 
                wave_parameters=circa_sd_ll["Mp8g00970",])
dev.off()

plot.ld.ll.wave("Mp2g08360", "MpTRIHELIX11", genes.expression, period=1, 
                wave_parameters=circa_ld_ll["Mp2g08360",])

# GI
"Mp1g13750" %in% rownames(circa_ld_ll)
"Mp1g13750" %in% rownames(circa_ld_dd)
"Mp1g13750" %in% rownames(circa_sd_ll)
"Mp1g13750" %in% rownames(circa_sd_dd)

# COL2
"Mp2g16350" %in% rownames(circa_ld_ll)
"Mp2g16350" %in% rownames(circa_ld_dd)
"Mp2g16350" %in% rownames(circa_sd_ll)
"Mp2g16350" %in% rownames(circa_sd_dd)

# LNK
"Mp1g10620" %in% rownames(circa_ld_ll)
"Mp1g10620" %in% rownames(circa_ld_dd)
"Mp1g10620" %in% rownames(circa_sd_ll)
"Mp1g10620" %in% rownames(circa_sd_dd)

# COL6
"Mp5g05940" %in% rownames(circa_ld_ll)
"Mp5g05940" %in% rownames(circa_ld_dd)
"Mp5g05940" %in% rownames(circa_sd_ll)
"Mp5g05940" %in% rownames(circa_sd_dd)

# WDR2
"Mp5g01830" %in% rownames(circa_ld_ll)
"Mp5g01830" %in% rownames(circa_ld_dd)
"Mp5g01830" %in% rownames(circa_sd_ll)
"Mp5g01830" %in% rownames(circa_sd_dd)

# CDF1
"Mp1g18540" %in% rownames(circa_ld_ll)
"Mp1g18540" %in% rownames(circa_ld_dd)
"Mp1g18540" %in% rownames(circa_sd_ll)
"Mp1g18540" %in% rownames(circa_sd_dd)


# ELF3
"Mp1g10870" %in% rownames(circa_ld_ll)
"Mp1g10870" %in% rownames(circa_ld_dd)
"Mp1g10870" %in% rownames(circa_sd_ll)
"Mp1g10870" %in% rownames(circa_sd_dd)


# Create clock genes table in each conditions
clock.names <- c("TOC1", "PRR", "RVE", "DET1", "GI", "COL2", 
                 "LNK", "COL6", "WDR2", "CDF1", "ELF3")

clock.genes <- c("Mp3g09460", "Mp2g16560", "Mp2g14310", "Mp8g00970", "Mp1g13750", "Mp2g16350",
                 "Mp1g10620", "Mp5g05940", "Mp5g01830", "Mp1g18540", "Mp1g10870")

circacompare.sd.dd <- read.table("tables/circacompare_sd_dd.tsv", header = T, sep = "\t")
circacompare.sd.ll <- read.table("tables/circacompare_sd_ll.tsv", header = T, sep = "\t")
circacompare.ld.dd <- read.table("tables/circacompare_ld_dd.tsv", header = T, sep = "\t")
circacompare.ld.ll <- read.table("tables/circacompare_ld_ll.tsv", header = T, sep = "\t")
circacompare.sd.dd.12 <- read.table("tables/circacompare_sd_dd_12.tsv", header = T, sep = "\t")
circacompare.sd.ll.12 <- read.table("tables/circacompare_sd_ll_12.tsv", header = T, sep = "\t")
circacompare.ld.dd.12 <- read.table("tables/circacompare_ld_dd_12.tsv", header = T, sep = "\t")
circacompare.ld.ll.12 <- read.table("tables/circacompare_ld_ll_12.tsv", header = T, sep = "\t")


clock_df <- data.frame(Symbol = clock.names, ID = clock.genes, 
           `LD-LL p-value` = pmin(circacompare.ld.ll[clock.genes, "rhythmic_p"],
                                  circacompare.ld.ll.12[clock.genes, "rhythmic_p"]),
           `SD-LL p-value` = pmin(circacompare.sd.ll[clock.genes, "rhythmic_p"],
                                  circacompare.sd.ll.12[clock.genes, "rhythmic_p"]),
           `LD-DD p-value` = pmin(circacompare.ld.dd[clock.genes, "rhythmic_p"],
                                  circacompare.ld.dd.12[clock.genes, "rhythmic_p"]),
           `SD-DD p-value` = pmin(circacompare.sd.dd[clock.genes, "rhythmic_p"],
                                  circacompare.sd.dd.12[clock.genes, "rhythmic_p"]))

clock_df_ordered <- clock_df[order(clock_df$LD.LL.p.value),]

write.table(clock_df_ordered,file="clock_genes_df.tsv", sep="\t", row.names = F)





###########################  PHOTOSYNTHESIS ##########################
scaled.sd.gene.expression <- t(scale(t(sd.gene.expression)))
scaled.ld.gene.expression <- t(scale(t(ld.gene.expression)))

source("functions_circa_marchantia.R")
library(circlize)

heatmap.gene(gene.id = "Mp2g16560")

# PSII

# PsbA Mp5g10950
# PsbD Mp1g11630 Mp1g19900 Mp4g18690 Mp5g07440
# PsbP Mp1g27590 Mp2g21180 Mp5g12920 Mp1g26600 Mp2g17390 Mp5g20290 Mp5g20300
# Mp8g10040 Mp1g00680 Mp5g12900 Mp7g05810 Mp8g10770 Mp2g02710 Mp1g02960 Mp5g12920
# PsbQ Mp1g20230 Mp2g15670 Mp4g07510
# PsbR Mp1g17690 Mp3g05840
# PsbW Mp5g14190
# PsbX Mp6g21460
# PsbY Mp8g12730
# PsbP Mp1g27590 Mp2g21180 Mp5g12920 Mp1g26600 Mp2g17390 Mp5g20290 Mp5g20300
# Mp8g10040 Mp1g00680 Mp5g12900 Mp7g05810 Mp8g10770 Mp2g02710 Mp1g02960 Mp5g12920

# LHCB2 Mp1g29620
# LHCB3 Mp7g08940
# LHCB7 Mp8g12010
# LHCB5 Mp4g10900
# lhcb4 Mp1g16850
# LHCB6 

# Others
# Mp2g13460 Mp7g09340 
# Mp7g09180 Mp7g06780 Mp7g06760 Mp7g06790 Mp7g06770 Mp7g06750 Mp7g06720 
# Mp7g06730 Mp7g06740 Mp7g05530 Mp7g05980 Mp7g06710 Mp7g05880 Mp7g05890 

#https://www.sciencedirect.com/science/article/pii/S2211124724010477#sec8


psb_ids <- c("Mp1g00680", "Mp2g02710", "Mp1g27590", "Mp7g05810", "Mp8g10040",
              "Mp8g10770", "Mp5g12920", "Mp5g20290", "Mp5g20300", "Mp4g07510",
              "Mp1g17690", "Mp3g05840", "Mp6g21460", "Mp8g12730", "Mp5g14190",
              "Mp2g14630")

lhcb_ids <- c("Mp7g05890", "Mp7g06740", "Mp7g09340", "Mp7g05880", "Mp7g05980",
              "Mp7g06710", "Mp7g06720", "Mp7g06730", "Mp7g06750", "Mp7g06770",
              "Mp7g06780", "Mp7g06790", "Mp7g09180", "Mp1g29620", "Mp2g13460",
              "Mp7g08940", "Mp1g16850", "Mp4g10900", "Mp6g01650", "Mp8g12010")

lhca_ids <- c("Mp2g15420", "Mp8g13180", "Mp4g00930", "Mp3g07840")

psa_ids <- c("Mp5g04200", "Mp1g24020", "Mp2g02960", "Mp6g18440", "Mp7g13220",
             "Mp4g06880", "Mp7g18150", "Mp4g20410")

pet_ids <- c("Mp2g07880", "Mp4g02720", "Mp2g15250", "Mp4g16550", "Mp6g00090",
             "Mp1g07170")

atp_ids <- c("Mp1g25120", "Mp6g06680")


# PsbP Mp1g00680 Mp2g02710 Mp1g27590 Mp7g05810 Mp8g10040 Mp8g10770 Mp5g12920 Mp5g20290 Mp5g20300
# PsbQ Mp4g07510
# PsbR Mp1g17690 Mp3g05840
# PsbX Mp6g21460
# PsbY Mp8g12730
# PsbW Mp5g14190
# Psb28 Mp2g14630

# LHCB1 Mp7g05890 Mp7g06740 Mp7g09340
# LHCB2 Mp7g05880 Mp7g05980 Mp7g06710 Mp7g06720 Mp7g06730 Mp7g06750 Mp7g06770
# Mp7g06780 Mp7g06790 Mp7g09180 Mp1g29620 Mp2g13460 
# Estos son nuevos de PharaohFUN: Mp7g06760 Mp7g05530   

# LHCB3 <- "Mp7g08940"
# LHCB4 <- "Mp1g16850"
# LHCB5 <- "Mp4g10900"
# LHCB6 <- "Mp6g01650"
# LHCB7 <- "Mp8g12010"
# 
# LHCA1 <- "Mp2g15420"
# LHCA2 <- "Mp8g13180"
# LHCA3 <- "Mp4g00930"
# LHCA4 <- "Mp3g07840"
# 
# psaD <- "Mp5g04200"
# psaE <- "Mp1g24020"
# psaF <- "Mp2g02960"
# psaG <- "Mp6g18440"
# psaH <- "Mp7g13220"
# psaK <- "Mp4g06880"
# psaL <- "Mp7g18150"
# psaN <- "Mp4g20410"
# 
# petC <- "Mp2g07880"
# petE <- "Mp4g02720"
# petF <- "Mp2g15250"
# petF <- "Mp4g16550"
# petH <- "Mp6g00090"
# petM <- "Mp1g07170"
# atpF <- "Mp1g25120"
# atpH <- "Mp6g06680"
# 
# PsbP <- c("Mp1g00680", "Mp2g02710", "Mp1g27590", "Mp7g05810", "Mp8g10040",
#           "Mp8g10770", "Mp5g12920", "Mp5g20290", "Mp5g20300")
# PsbQ <- "Mp4g07510"
# PsbR <- c("Mp1g17690", "Mp3g05840")
# PsbX <- "Mp6g21460"
# PsbY <- "Mp8g12730"
# PsbW <- "Mp5g14190"
# Psb28<- "Mp2g14630"
# 
# LHCB1 <- c("Mp7g05890", "Mp7g06740", "Mp7g09340")
# LHCB2 <- c("Mp7g05880", "Mp7g05980", "Mp7g06710", "Mp7g06720", "Mp7g06730", "Mp7g06750", "Mp7g06770",
#            "Mp7g06780", "Mp7g06790", "Mp7g09180", "Mp1g29620", "Mp2g13460")
# LHCB3 <- "Mp7g08940"
# LHCB4 <- "Mp1g16850"
# LHCB5 <- "Mp4g10900"
# LHCB6 <- "Mp6g01650"
# LHCB7 <- "Mp8g12010"
# LHCA1 <- "Mp2g15420"
# LHCA2 <- "Mp8g13180"
# LHCA3 <- "Mp4g00930"
# LHCA4 <- "Mp3g07840"
# psaD <- "Mp5g04200"
# psaE <- "Mp1g24020"
# psaF <- "Mp2g02960"
# psaG <- "Mp6g18440"
# psaH <- "Mp7g13220"
# psaK <- "Mp4g06880"
# psaL <- "Mp7g18150"
# psaN <- "Mp4g20410"
# petC <- "Mp2g07880"
# petE <- "Mp4g02720"
# petF <- "Mp2g15250"
# petF <- "Mp4g16550"
# petH <- "Mp6g00090"
# petM <- "Mp1g07170"
# atpF <- "Mp1g25120"
# atpH <- "Mp6g06680"

# PsbP <- ("Mp1g00680 Mp2g02710 Mp1g27590 Mp7g05810 Mp8g10040 Mp8g10770 Mp5g12920 Mp5g20290 Mp5g20300
# PsbQ Mp4g07510
# PsbR Mp1g17690 Mp3g05840
# PsbX Mp6g21460
# PsbY Mp8g12730
# PsbW Mp5g14190
# Psb28 Mp2g14630

# LHCB1 Mp7g05890 Mp7g06740 Mp7g09340
# LHCB2 Mp7g05880 Mp7g05980 Mp7g06710 Mp7g06720 Mp7g06730 Mp7g06750 Mp7g06770
# Mp7g06780 Mp7g06790 Mp7g09180 Mp1g29620 Mp2g13460 
# Estos son nuevos de PharaohFUN: Mp7g06760 Mp7g05530 

LHCA1 <- "Mp2g15420"
LHCA2 <- "Mp8g13180"
LHCA3 <- "Mp4g00930"
LHCA4 <- "Mp3g07840"

psaD <- "Mp5g04200"
psaE <- "Mp1g24020"
psaF <- "Mp2g02960"
psaG <- "Mp6g18440"
psaH <- "Mp7g13220"
psaK <- "Mp4g06880"
psaL <- "Mp7g18150"
psaN <- "Mp4g20410"

petC <- "Mp2g07880"
petE <- "Mp4g02720"
petF <- c("Mp2g15250", "Mp4g16550")
petH <- "Mp6g00090"
petM <- "Mp1g07170"

atpF <- "Mp1g25120"
atpH <- "Mp6g06680"

PsbP <- c("Mp1g00680", "Mp2g02710", "Mp1g27590", "Mp7g05810", "Mp8g10040",
          "Mp8g10770", "Mp5g12920", "Mp5g20290", "Mp5g20300")
PsbQ <- "Mp4g07510"
PsbR <- c("Mp1g17690", "Mp3g05840")
PsbX <- "Mp6g21460"
PsbY <- "Mp8g12730"
PsbW <- "Mp5g14190"
Psb28<- "Mp2g14630"

LHCB1 <- c("Mp7g05890", "Mp7g06740", "Mp7g09340")
LHCB2 <- c("Mp7g05880", "Mp7g05980", "Mp7g06710", "Mp7g06720", "Mp7g06730", "Mp7g06750", "Mp7g06770",
           "Mp7g06780", "Mp7g06790", "Mp7g09180", "Mp1g29620", "Mp2g13460")

# Estos son nuevos de PharaohFUN, decidir si incluirlos en LHCB2: Mp7g06760 Mp7g05530 

LHCB3 <- "Mp7g08940"
LHCB4 <- "Mp1g16850"
LHCB5 <- "Mp4g10900"
LHCB6 <- "Mp6g01650"
LHCB7 <- "Mp8g12010"



for (i in 1:4)
{
  png(paste0("photo_heatmaps/LHCA/LHCA", i))
  heatmap.gene(gene.id = c(LHCA1, LHCA2, LHCA3, LHCA4)[i])
  dev.off()
}

for (i in 1:8)
{
  png(paste0("photo_heatmaps/psa/psa", c("D","E","F","G","H","K","L","N")[i]))
  heatmap.gene(gene.id = c(psaD,psaE,psaF,psaG,psaH,psaK,psaL,psaN)[i])
  dev.off()
}

for (i in 1:6)
{
  png(paste0("photo_heatmaps/pet/pet", c("C","E","F_1","F_2","H","M")[i]))
  heatmap.gene(gene.id = c(petC,petE,petF,petH,petM)[i])
  dev.off()
}

png(paste0("photo_heatmaps/atp/atpF"))
heatmap.gene(gene.id = "Mp1g25120")
dev.off()

png(paste0("photo_heatmaps/atp/atpH"))
heatmap.gene(gene.id = "Mp6g06680")
dev.off()

for (i in 1:16)
{
  png(paste0("photo_heatmaps/psb/psb", c(paste0(rep("P_", length(PsbP)), 1:length(PsbP)),
                                         paste0(rep("Q_", length(PsbQ)), 1:length(PsbQ)),
                                         paste0(rep("R_", length(PsbR)), 1:length(PsbR)),
                                         paste0(rep("X_", length(PsbX)), 1:length(PsbX)),
                                         paste0(rep("Y_", length(PsbY)), 1:length(PsbY)),
                                         paste0(rep("W_", length(PsbW)), 1:length(PsbW)),
                                         paste0(rep("28_", length(Psb28)), 1:length(Psb28)))[i]))
                                         
  heatmap.gene(gene.id = c(PsbP,PsbQ,PsbR,PsbX,PsbY,PsbW,Psb28)[i])
  dev.off()
}

for (i in 1:20)
{
  png(paste0("photo_heatmaps/LHCB/LHCB", c(paste0(rep("1_", length(LHCB1)), 1:length(LHCB1)),
                                         paste0(rep("2_", length(LHCB2)), 1:length(LHCB2)),
                                         "3", "4", "5", "6", "7")[i]))
  
  heatmap.gene(gene.id = c(LHCB1, LHCB2, LHCB3, LHCB4, LHCB5, LHCB6, LHCB7)[i])
  dev.off()
}

# PsbS
png("images/figure_photo/heatmap_psbs.png", width = 1400, height = 700)
heatmap.gene(gene.id = "Mp8g05910")
dev.off()

# SEP1
png("images/figure_photo/heatmap_sep1.png", width = 1400, height = 700)
heatmap.gene(gene.id = "Mp3g06770")
dev.off()

# VDE
png("images/figure_photo/heatmap_vde.png", width = 1400, height = 700)
heatmap.gene(gene.id = "Mp7g15960")
dev.off()


# Explore rhythmicity of PsbS and VDE
complete.ld.rhythmic.genes <- read.table("tables/ld.rhythmic.table.tsv", header = T)$V1
complete.sd.rhythmic.genes <- read.table("tables/sd.rhythmic.table.tsv", header = T)$V1
complete.sd.dd.rhythmic.genes <- read.table("tables/sd.dd.rhythmic.table.tsv", sep="\t")$x
complete.ld.dd.rhythmic.genes <- read.table("tables/ld.dd.rhythmic.table.tsv", sep="\t")$x
complete.sd.ll.rhythmic.genes <- read.table("tables/sd.ll.rhythmic.table.tsv", sep="\t")$x
complete.ld.ll.rhythmic.genes <- read.table("tables/ld.ll.rhythmic.table.tsv", sep="\t")$x

"Mp8g05910" %in% complete.ld.rhythmic.genes
"Mp8g05910" %in% complete.sd.rhythmic.genes
"Mp8g05910" %in% complete.ld.ll.rhythmic.genes
"Mp8g05910" %in% complete.sd.ll.rhythmic.genes
"Mp8g05910" %in% complete.ld.dd.rhythmic.genes
"Mp8g05910" %in% complete.sd.dd.rhythmic.genes # F

"Mp7g15960" %in% complete.ld.rhythmic.genes
"Mp7g15960" %in% complete.sd.rhythmic.genes
"Mp7g15960" %in% complete.ld.ll.rhythmic.genes
"Mp7g15960" %in% complete.sd.ll.rhythmic.genes # F
"Mp7g15960" %in% complete.ld.dd.rhythmic.genes # F
"Mp7g15960" %in% complete.sd.dd.rhythmic.genes # F



# Summary of distributions
# LHCB + psb
# SD
combined_lhcb_psb <- c(LHCB1, LHCB2, LHCB3, LHCB4, LHCB5, LHCB6, LHCB7,
                       PsbP,PsbQ,PsbR,PsbX,PsbY,PsbW,Psb28) 


scaled_sd_combined_lhcb_psb <- scaled.sd.gene.expression[combined_lhcb_psb,]
scaled_sd_combined_lhcb_psb_mean <- colMeans(scaled_sd_combined_lhcb_psb)

scaled_sd_combined_lhcb_psb_mean_zt <- tapply(scaled_sd_combined_lhcb_psb_mean, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled_sd_combined_lhcb_psb_mean_zt <- c(scaled_sd_combined_lhcb_psb_mean_zt, scaled_sd_combined_lhcb_psb_mean_zt[1])
spline.current.gene.mean.sd.gene.expression <- spline(x=seq(from=0,to=24,by=4),y=scaled_sd_combined_lhcb_psb_mean_zt,n = 100)

# Adjust min value to be 0
spline.current.gene.mean.sd.gene.expression$y <- spline.current.gene.mean.sd.gene.expression$y - 
  min(spline.current.gene.mean.sd.gene.expression$y)

data_sd_combined_lhcb_psb <- data.frame(x = spline.current.gene.mean.sd.gene.expression$x,
                                        y = spline.current.gene.mean.sd.gene.expression$y)

ggplot(data_sd_combined_lhcb_psb, aes(x=x, y=y)) +
  geom_area(fill="firebrick2", colour = "darkred", linewidth = 1.5) + theme(
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  )
ggsave("lhcb_psb_summary_sd.png", bg = "transparent", width = 12, height = 6)


# LD
scaled_ld_combined_lhcb_psb <- scaled.ld.gene.expression[combined_lhcb_psb,]
scaled_ld_combined_lhcb_psb_mean <- colMeans(scaled_ld_combined_lhcb_psb)

scaled_ld_combined_lhcb_psb_mean_zt <- tapply(scaled_ld_combined_lhcb_psb_mean, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled_ld_combined_lhcb_psb_mean_zt <- c(scaled_ld_combined_lhcb_psb_mean_zt, scaled_ld_combined_lhcb_psb_mean_zt[1])
spline.current.gene.mean.ld.gene.expression <- spline(x=seq(from=0,to=24,by=4),y=scaled_ld_combined_lhcb_psb_mean_zt,n = 100)

# Adjust min value to be 0
spline.current.gene.mean.ld.gene.expression$y <- spline.current.gene.mean.ld.gene.expression$y - 
  min(spline.current.gene.mean.ld.gene.expression$y)

data_ld_combined_lhcb_psb <- data.frame(x = spline.current.gene.mean.ld.gene.expression$x,
                                        y = spline.current.gene.mean.ld.gene.expression$y)

ggplot(data_ld_combined_lhcb_psb, aes(x=x, y=y)) +
  geom_area(fill="royalblue3", colour = "darkblue", linewidth = 1.5) + theme(
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  )
ggsave("lhcb_psb_summary_ld.png", bg = "transparent", width = 12, height = 6)

# pet
combined_pet <- c(petC,petE,petF,petH,petM) 


scaled_sd_combined_pet <- scaled.sd.gene.expression[combined_pet,]
scaled_sd_combined_pet_mean <- colMeans(scaled_sd_combined_pet)

scaled_sd_combined_pet_mean_zt <- tapply(scaled_sd_combined_pet_mean, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled_sd_combined_pet_mean_zt <- c(scaled_sd_combined_pet_mean_zt, scaled_sd_combined_pet_mean_zt[1])
spline.current.gene.mean.sd.gene.expression <- spline(x=seq(from=0,to=24,by=4),y=scaled_sd_combined_pet_mean_zt,n = 100)

# Adjust min value to be 0
spline.current.gene.mean.sd.gene.expression$y <- spline.current.gene.mean.sd.gene.expression$y - 
  min(spline.current.gene.mean.sd.gene.expression$y)

data_sd_combined_pet <- data.frame(x = spline.current.gene.mean.sd.gene.expression$x,
                                        y = spline.current.gene.mean.sd.gene.expression$y)

ggplot(data_sd_combined_pet, aes(x=x, y=y)) +
  geom_area(fill="firebrick2", colour = "darkred", linewidth = 1.5) + theme(
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  )
ggsave("pet_summary_sd.png", bg = "transparent", width = 12, height = 6)

scaled_ld_combined_pet <- scaled.ld.gene.expression[combined_pet,]
scaled_ld_combined_pet_mean <- colMeans(scaled_ld_combined_pet)

scaled_ld_combined_pet_mean_zt <- tapply(scaled_ld_combined_pet_mean, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled_ld_combined_pet_mean_zt <- c(scaled_ld_combined_pet_mean_zt, scaled_ld_combined_pet_mean_zt[1])
spline.current.gene.mean.ld.gene.expression <- spline(x=seq(from=0,to=24,by=4),y=scaled_ld_combined_pet_mean_zt,n = 100)

# Adjust min value to be 0
spline.current.gene.mean.ld.gene.expression$y <- spline.current.gene.mean.ld.gene.expression$y - 
  min(spline.current.gene.mean.ld.gene.expression$y)

data_ld_combined_pet <- data.frame(x = spline.current.gene.mean.ld.gene.expression$x,
                                        y = spline.current.gene.mean.ld.gene.expression$y)

ggplot(data_ld_combined_pet, aes(x=x, y=y)) +
  geom_area(fill="royalblue3", colour = "darkblue", linewidth = 1.5) + theme(
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  )
ggsave("pet_summary_ld.png", bg = "transparent", width = 12, height = 6)

# lhca + psa
combined_lhca_psa <- c(LHCA1, LHCA2, LHCA3, LHCA4, 
                       psaD,psaE,psaF,psaG,psaH,psaK,psaL,psaN)


scaled_sd_combined_lhca_psa <- scaled.sd.gene.expression[combined_lhca_psa,]
scaled_sd_combined_lhca_psa_mean <- colMeans(scaled_sd_combined_lhca_psa)

scaled_sd_combined_lhca_psa_mean_zt <- tapply(scaled_sd_combined_lhca_psa_mean, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled_sd_combined_lhca_psa_mean_zt <- c(scaled_sd_combined_lhca_psa_mean_zt, scaled_sd_combined_lhca_psa_mean_zt[1])
spline.current.gene.mean.sd.gene.expression <- spline(x=seq(from=0,to=24,by=4),y=scaled_sd_combined_lhca_psa_mean_zt,n = 100)

# Adjust min value to be 0
spline.current.gene.mean.sd.gene.expression$y <- spline.current.gene.mean.sd.gene.expression$y - 
  min(spline.current.gene.mean.sd.gene.expression$y)

data_sd_combined_lhca_psa <- data.frame(x = spline.current.gene.mean.sd.gene.expression$x,
                                   y = spline.current.gene.mean.sd.gene.expression$y)

ggplot(data_sd_combined_lhca_psa, aes(x=x, y=y)) +
  geom_area(fill="firebrick2", colour = "darkred", linewidth = 1.5) + theme(
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  )
ggsave("lhca_psa_summary_sd.png", bg = "transparent", width = 12, height = 6)

scaled_ld_combined_lhca_psa <- scaled.ld.gene.expression[combined_lhca_psa,]
scaled_ld_combined_lhca_psa_mean <- colMeans(scaled_ld_combined_lhca_psa)

scaled_ld_combined_lhca_psa_mean_zt <- tapply(scaled_ld_combined_lhca_psa_mean, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled_ld_combined_lhca_psa_mean_zt <- c(scaled_ld_combined_lhca_psa_mean_zt, scaled_ld_combined_lhca_psa_mean_zt[1])
spline.current.gene.mean.ld.gene.expression <- spline(x=seq(from=0,to=24,by=4),y=scaled_ld_combined_lhca_psa_mean_zt,n = 100)

# Adjust min value to be 0
spline.current.gene.mean.ld.gene.expression$y <- spline.current.gene.mean.ld.gene.expression$y - 
  min(spline.current.gene.mean.ld.gene.expression$y)

data_ld_combined_lhca_psa <- data.frame(x = spline.current.gene.mean.ld.gene.expression$x,
                                   y = spline.current.gene.mean.ld.gene.expression$y)

ggplot(data_ld_combined_lhca_psa, aes(x=x, y=y)) +
  geom_area(fill="royalblue3", colour = "darkblue", linewidth = 1.5) + theme(
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  )
ggsave("lhca_psa_summary_ld.png", bg = "transparent", width = 12, height = 6)


# atp
combined_atp <- c(atpF, atpH)


scaled_sd_combined_atp <- scaled.sd.gene.expression[combined_atp,]
scaled_sd_combined_atp_mean <- colMeans(scaled_sd_combined_atp)

scaled_sd_combined_atp_mean_zt <- tapply(scaled_sd_combined_atp_mean, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled_sd_combined_atp_mean_zt <- c(scaled_sd_combined_atp_mean_zt, scaled_sd_combined_atp_mean_zt[1])
spline.current.gene.mean.sd.gene.expression <- spline(x=seq(from=0,to=24,by=4),y=scaled_sd_combined_atp_mean_zt,n = 100)

# Adjust min value to be 0
spline.current.gene.mean.sd.gene.expression$y <- spline.current.gene.mean.sd.gene.expression$y - 
  min(spline.current.gene.mean.sd.gene.expression$y)

data_sd_combined_atp <- data.frame(x = spline.current.gene.mean.sd.gene.expression$x,
                                        y = spline.current.gene.mean.sd.gene.expression$y)

ggplot(data_sd_combined_atp, aes(x=x, y=y)) +
  geom_area(fill="firebrick2", colour = "darkred", linewidth = 1.5) + theme(
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  )
ggsave("atp_summary_sd.png", bg = "transparent", width = 12, height = 6)

scaled_ld_combined_atp <- scaled.ld.gene.expression[combined_atp,]
scaled_ld_combined_atp_mean <- colMeans(scaled_ld_combined_atp)

scaled_ld_combined_atp_mean_zt <- tapply(scaled_ld_combined_atp_mean, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled_ld_combined_atp_mean_zt <- c(scaled_ld_combined_atp_mean_zt, scaled_ld_combined_atp_mean_zt[1])
spline.current.gene.mean.ld.gene.expression <- spline(x=seq(from=0,to=24,by=4),y=scaled_ld_combined_atp_mean_zt,n = 100)

# Adjust min value to be 0
spline.current.gene.mean.ld.gene.expression$y <- spline.current.gene.mean.ld.gene.expression$y - 
  min(spline.current.gene.mean.ld.gene.expression$y)

data_ld_combined_atp <- data.frame(x = spline.current.gene.mean.ld.gene.expression$x,
                                        y = spline.current.gene.mean.ld.gene.expression$y)

ggplot(data_ld_combined_atp, aes(x=x, y=y)) +
  geom_area(fill="royalblue3", colour = "darkblue", linewidth = 1.5) + theme(
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  )
ggsave("atp_summary_ld.png", bg = "transparent", width = 12, height = 6)


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

# Heatmaps
png(paste0("chloro_heatmaps/phot.png"))
heatmap.gene(gene.id = "Mp5g03810")
dev.off()

png(paste0("chloro_heatmaps/chup1.png"))
heatmap.gene(gene.id = "Mp7g03540")
dev.off()

png(paste0("chloro_heatmaps/kac.png"))
heatmap.gene(gene.id = "Mp3g03760")
dev.off()

png(paste0("chloro_heatmaps/blus1.png"))
heatmap.gene(gene.id = "Mp8g10870")
dev.off()

png(paste0("chloro_heatmaps/thrumin1_1.png"))
heatmap.gene(gene.id = "Mp1g07590")
dev.off()

png(paste0("chloro_heatmaps/thrumin1_2.png"))
heatmap.gene(gene.id = "Mp3g04860")
dev.off()

png(paste0("chloro_heatmaps/jac1_1.png"))
heatmap.gene(gene.id = "Mp6g17960") # Not expressed
dev.off()

png(paste0("chloro_heatmaps/jac1_2.png"))
heatmap.gene(gene.id = "Mp6g06160") # Not expressed
dev.off()

png(paste0("chloro_heatmaps/jac1_3.png"))
heatmap.gene(gene.id = "Mp3g09960")
dev.off()

png(paste0("chloro_heatmaps/thrumin1_4.png"))
heatmap.gene(gene.id = "Mp5g14210")
dev.off()

png(paste0("chloro_heatmaps/pmi1_rel_1.png"))
heatmap.gene(gene.id = "Mp2g25860")
dev.off()

png(paste0("chloro_heatmaps/pmi1_rel_2.png"))
heatmap.gene(gene.id = "Mp5g06760")
dev.off()

png(paste0("chloro_heatmaps/pmi1_rel_3.png"))
heatmap.gene(gene.id = "Mp8g10000")
dev.off()

png(paste0("chloro_heatmaps/pmi1_rel_4.png"))
heatmap.gene(gene.id = "Mp8g18590")
dev.off()

png(paste0("chloro_heatmaps/pmi1_rel_5.png"))
heatmap.gene(gene.id = "Mp8g18700")
dev.off()

png(paste0("chloro_heatmaps/nch1.png"))
heatmap.gene(gene.id = "Mp5g07060")
dev.off()

png(paste0("chloro_heatmaps/rtn1.png"))
heatmap.gene(gene.id = "Mp6g03260")
dev.off()
# Exploration of rhythmicity among photorelocation genes
source("functions_circa_marchantia.R")
complete.ld.rhythmic.genes <- read.table("tables/ld.rhythmic.table.tsv", header = T)$V1
complete.sd.rhythmic.genes <- read.table("tables/sd.rhythmic.table.tsv", header = T)$V1

# Photoreceptor PHOT and its substrate BLUS1
"Mp5g03810" %in% complete.ld.rhythmic.genes
"Mp5g03810" %in% complete.sd.rhythmic.genes

"Mp8g10870" %in% complete.ld.rhythmic.genes
"Mp8g10870" %in% complete.sd.rhythmic.genes

plot.ld.sd("Mp5g03810", "PHOT", genes.expression)
plot.ld.sd("Mp8g10870", "BLUS", genes.expression)

# BLUS1 Not involved in this
# RTN1
"Mp6g03260" %in% complete.ld.rhythmic.genes
"Mp6g03260" %in% complete.sd.rhythmic.genes

#NCH1
"Mp5g07060" %in% complete.ld.rhythmic.genes
"Mp5g07060" %in% complete.sd.rhythmic.genes



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

# Motility CHUP1, KAC, THRUM1, THRUM1 and PMI1
c("Mp7g03540", "Mp3g03760", "Mp1g07590", "Mp3g04860", "Mp2g25860") %in% complete.ld.rhythmic.genes
c("Mp7g03540", "Mp3g03760", "Mp1g07590", "Mp3g04860", "Mp2g25860") %in% complete.sd.rhythmic.genes

plot.ld.sd("Mp3g03760", "KAC", genes.expression)
plot.ld.sd("Mp7g03540", "CHUP1", genes.expression)

# Signaling JAC1, JAC1, JAC1 and PMI2x4
c("Mp6g17960", "Mp6g06160", "Mp3g09960", "Mp5g06760", "Mp8g10000", "Mp8g18590", "Mp8g18700") %in% complete.ld.rhythmic.genes
c("Mp6g17960", "Mp6g06160", "Mp3g09960", "Mp5g06760", "Mp8g10000", "Mp8g18590", "Mp8g18700") %in% complete.sd.rhythmic.genes

plot.ld.sd("Mp5g06760", "PMIr", genes.expression)
png("phot.png")
plot.ld.sd.wave("Mp5g03810", "PHOT",  genes.expression, wave_parameters =  circa_table_genescarab["Mp5g03810",])
dev.off()


### Photoreceptors ###

complete.ld.rhythmic.genes <- read.table("tables/ld.rhythmic.table.tsv", header = T)$V1
complete.sd.rhythmic.genes <- read.table("tables/sd.rhythmic.table.tsv", header = T)$V1

source("functions_circa_marchantia.R")
# MpPHY, MpCRY, MpUVR8, MpPHOT

photoreceptors <- c("Mp2g16090", "Mp2g17590", "Mp2g11590", "Mp5g03810")
photoreceptors %in% complete.ld.rhythmic.genes
photoreceptors %in% complete.sd.rhythmic.genes

plot.ld.sd("Mp2g16090", "MpPHY", genes.expression = genes.expression)
plot.ld.sd("Mp2g17590", "MpCRY", genes.expression = genes.expression)
plot.ld.sd("Mp2g11590", "MpUVR8", genes.expression = genes.expression)
plot.ld.sd("Mp5g03810", "MpPHOT", genes.expression = genes.expression)


circa_table_sd_dd <- subset(read.table("tables/circacompare_sd_dd.tsv", header = T, sep = "\t"), rhythmic_p < 0.05)
circa_table_ld_dd <- subset(read.table("tables/circacompare_ld_dd.tsv", header = T, sep = "\t"), rhythmic_p < 0.05)
circa_table_sd_ll <- subset(read.table("tables/circacompare_sd_ll.tsv", header = T, sep = "\t"), rhythmic_p < 0.05)
circa_table_ld_ll <- subset(read.table("tables/circacompare_ld_ll.tsv", header = T, sep = "\t"), rhythmic_p < 0.05)

photoreceptors[photoreceptors %in% rownames(circa_table_ld_ll)]
photoreceptors[photoreceptors %in% rownames(circa_table_sd_ll)]
photoreceptors[photoreceptors %in% rownames(circa_table_ld_dd)]
photoreceptors[photoreceptors %in% rownames(circa_table_sd_dd)]

# MpPHOT is rhythmic in LD-LL and SD-LL


# Expression in each condition

mean(ld.gene.expression[photoreceptors,])
mean(ld.ll.gene.expression[photoreceptors,])
mean(ld.dd.gene.expression[photoreceptors,])
mean(sd.gene.expression[photoreceptors,])
mean(sd.ll.gene.expression[photoreceptors,])
mean(sd.dd.gene.expression[photoreceptors,])




# ESTO SE TIENE QUE IR DE AQUÍ EN CUANTO SE CONFIRME QUE ESTÁ EN EL DE FUNCIONES
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
