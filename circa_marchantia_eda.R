.libPaths("/home/marcos/R/x86_64-pc-linux-gnu-library/4.1")

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


# PCAs to visualize time point and harvesting day effects
library(FactoMineR)
library(factoextra)

# Cycle conditions
# Clustering of LD and SD samples
pca.gene.expression.ld.sd <- data.frame(colnames(ld.sd.gene.expression),t(ld.sd.gene.expression))
dim(pca.gene.expression.ld.sd)

head(pca.gene.expression.ld.sd[,1:9])
colnames(pca.gene.expression.ld.sd)[1] <- "Time point"
res.pca.ld.sd <- PCA(pca.gene.expression.ld.sd, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
res.hcpc.ld.sd <- HCPC(res.pca.ld.sd, graph=FALSE)
png("hcpc_ld_sd.png", width = 2000, height = 2000, res=300)
fviz_dend(res.hcpc.ld.sd,k=2,cex = 1,palette = "jco",rect = TRUE, rect_fill = TRUE,rect_border = "jco",type="rectangle",
          labels_track_height = 500)
dev.off()

# PCA LD
ld.gene.expression <- genes.expression[,ld.zt.i]
dim(ld.gene.expression)
pca.gene.expression.ld <- data.frame(colnames(ld.gene.expression),t(ld.gene.expression))
colnames(pca.gene.expression.ld)[1] <- "Time point"
pca.gene.expression.ld[1:12,1:5]
res.pca.ld <- PCA(pca.gene.expression.ld, graph = FALSE,scale.unit = F,quali.sup = 1 )
png("pca_ld.png", width = 1500, height = 1500, res=300)
fviz_pca_ind(res.pca.ld, col.ind = c(rep("ZT00", 3),rep("ZT04", 3), rep("ZT08", 3), 
                                     rep("ZT12", 3), rep("ZT16", 3), rep("ZT20", 3)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)
dev.off()

# PCA SD
sd.gene.expression <- genes.expression[,sd.zt.i]
pca.gene.expression.sd <- data.frame(colnames(sd.gene.expression),t(sd.gene.expression))
pca.gene.expression.sd[1:12,1:12]
colnames(pca.gene.expression.sd)[1] <- "Time point"
res.pca.sd <- PCA(pca.gene.expression.sd, graph = FALSE,scale.unit = F,quali.sup = 1 )
png("pca_sd.png", width = 1500, height = 1500, res=300)
fviz_pca_ind(res.pca.sd, col.ind = c(rep("ZT00", 3),rep("ZT04", 3), rep("ZT08", 3), 
                                     rep("ZT12", 3), rep("ZT16", 3), rep("ZT20", 3)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)
dev.off()

# Gene expression levels per sample
boxplot(log2(genes.expression + 1),col=rainbow(ncol(genes.expression)),ylab="log2(FPKM + 1)",cex.lab=1.5,las=2)
