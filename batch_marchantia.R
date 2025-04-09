.libPaths("/home/marcos/R/x86_64-pc-linux-gnu-library/4.1")

# Read data
genes.expression <- read.table(file ="marchantia_gene_expression.tsv",sep="\t",header=T,as.is=T)
dim(genes.expression)
head(genes.expression)
genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)


# Gene expression matrix for samples under light-dark cycles
ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ld.zt.i <- sapply(X = ld.zt,FUN = function(x){ paste(x,1:3,sep="_")})
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt.i <- sapply(X = sd.zt,FUN = function(x){ paste(x,1:3,sep="_")})

ld.sd.gene.expression <- genes.expression[,c(ld.zt.i,sd.zt.i)]
head(ld.sd.gene.expression)
dim(ld.sd.gene.expression)

# Free running conditions
dd.ld.zt <- paste("dd_ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
dd.ld.zt.i <- sapply(X = dd.ld.zt,FUN = function(x){ paste(x,1:2,sep="_")})
dd.sd.zt <- paste("dd_sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
dd.sd.zt.i <- sapply(X = dd.sd.zt,FUN = function(x){ paste(x,1:2,sep="_")})

ll.ld.zt <- paste("ll_ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ll.ld.zt.i <- sapply(X = ll.ld.zt,FUN = function(x){ paste(x,1:2,sep="_")})
ll.sd.zt <- paste("ll_sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ll.sd.zt.i <- sapply(X = ll.sd.zt,FUN = function(x){ paste(x,1:2,sep="_")})

dd.gene.expression <- genes.expression[,c(dd.ld.zt.i,dd.sd.zt.i)]
head(dd.gene.expression)
dim(dd.gene.expression)

ll.gene.expression <- genes.expression[,c(ll.ld.zt.i,ll.sd.zt.i)]
head(ll.gene.expression)
dim(ll.gene.expression)

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
fviz_dend(res.hcpc.ld.sd,k=2,cex = 1,palette = "jco",rect = TRUE, rect_fill = TRUE,rect_border = "jco",type="rectangle",
          labels_track_height = 500)


# PCA LD
ld.gene.expression <- genes.expression[,ld.zt.i]
dim(ld.gene.expression)
pca.gene.expression.ld <- data.frame(colnames(ld.gene.expression),t(ld.gene.expression))
colnames(pca.gene.expression.ld)[1] <- "Time point"
pca.gene.expression.ld[1:12,1:5]
res.pca.ld <- PCA(pca.gene.expression.ld, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_pca_ind(res.pca.ld, col.ind = c(rep("ZT00", 3),rep("ZT04", 3), rep("ZT08", 3), 
                                     rep("ZT12", 3), rep("ZT16", 3), rep("ZT20", 3)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)


# PCA SD
sd.gene.expression <- genes.expression[,sd.zt.i]
pca.gene.expression.sd <- data.frame(colnames(sd.gene.expression),t(sd.gene.expression))
pca.gene.expression.sd[1:12,1:12]
colnames(pca.gene.expression.sd)[1] <- "Time point"
res.pca.sd <- PCA(pca.gene.expression.sd, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_pca_ind(res.pca.sd, col.ind = c(rep("ZT00", 3),rep("ZT04", 3), rep("ZT08", 3), 
                                     rep("ZT12", 3), rep("ZT16", 3), rep("ZT20", 3)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)


# Free-running conditions DD
# Clustering of DD-running conditions samples
pca.gene.expression.dd <- data.frame(colnames(dd.gene.expression),t(dd.gene.expression))
dim(pca.gene.expression.dd)
colnames(pca.gene.expression.dd)[1] <- "Time point"
res.pca.dd<- PCA(pca.gene.expression.dd, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
res.hcpc.dd <- HCPC(res.pca.dd, graph=FALSE, nb.clust = 12)
fviz_dend(res.hcpc.dd,k=12,cex = 1,palette = "jco",rect = TRUE, rect_fill = TRUE,rect_border = "jco",type="rectangle",
          labels_track_height = 500)

# PCA DD LD
dd.ld.gene.expression <- genes.expression[,dd.ld.zt.i]
pca.gene.expression.dd.ld <- data.frame(colnames(dd.ld.gene.expression),t(dd.ld.gene.expression))
pca.gene.expression.dd.ld[1:12,1:12]
colnames(pca.gene.expression.dd.ld)[1] <- "Time point"
res.pca.dd.ld <- PCA(pca.gene.expression.dd.ld, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_pca_ind(res.pca.dd.ld, col.ind = c(rep("ZT00", 2),rep("ZT04", 2), rep("ZT08", 2), 
                                        rep("ZT12", 2), rep("ZT16", 2), rep("ZT20", 2)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)

# PCA DD SD
dd.sd.gene.expression <- genes.expression[,dd.sd.zt.i]
pca.gene.expression.dd.sd <- data.frame(colnames(dd.sd.gene.expression),t(dd.sd.gene.expression))
pca.gene.expression.dd.sd[1:12,1:12]
colnames(pca.gene.expression.dd.sd)[1] <- "Time point"
res.pca.dd.sd <- PCA(pca.gene.expression.dd.sd, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_pca_ind(res.pca.dd.sd, col.ind = c(rep("ZT00", 2),rep("ZT04", 2), rep("ZT08", 2), 
                                     rep("ZT12", 2), rep("ZT16", 2), rep("ZT20", 2)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)

# Free-running conditions LL
# Clustering of LL-running conditions samples
pca.gene.expression.ll <- data.frame(colnames(ll.gene.expression),t(ll.gene.expression))
dim(pca.gene.expression.ll)
colnames(pca.gene.expression.ll)[1] <- "Time point"
res.pca.ll<- PCA(pca.gene.expression.ll, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
res.hcpc.ll <- HCPC(res.pca.ll, graph=FALSE, nb.clust = 12)
fviz_dend(res.hcpc.ll,k=12,cex = 1,palette = "jco",rect = TRUE, rect_fill = TRUE,rect_border = "jco",type="rectangle",
          labels_track_height = 500)

# PCA LL LD
ll.ld.gene.expression <- genes.expression[,ll.ld.zt.i]
pca.gene.expression.ll.ld <- data.frame(colnames(ll.ld.gene.expression),t(ll.ld.gene.expression))
pca.gene.expression.ll.ld[1:12,1:12]
colnames(pca.gene.expression.ll.ld)[1] <- "Time point"
res.pca.ll.ld <- PCA(pca.gene.expression.ll.ld, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
res.pca.ll.ld <- res.pca.ll_ld
fviz_pca_ind(res.pca.ll.ld, col.ind = rep(c("ZT00", "ZT04", "ZT08",  "ZT12", "ZT16", "ZT20"), each = 2),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)



# PCA LL SD
ll.sd.gene.expression <- genes.expression[,ll.sd.zt.i]
pca.gene.expression.ll.sd <- data.frame(colnames(ll.sd.gene.expression),t(ll.sd.gene.expression))
pca.gene.expression.ll.sd[1:12,1:12]
colnames(pca.gene.expression.ll.sd)[1] <- "Time point"
res.pca.ll.sd <- PCA(pca.gene.expression.ll.sd, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_pca_ind(res.pca.ll.sd, col.ind = c(rep("ZT00", 2),rep("ZT04", 2), rep("ZT08", 2), 
                                        rep("ZT12", 2), rep("ZT16", 2), rep("ZT20", 2)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)


# PCA of LL LD samples
pca.gene.expression.ll_ld <- data.frame(colnames(ll.ld.gene.expression),t(ll.ld.gene.expression))
colnames(pca.gene.expression.ll_ld)[1] <- "Time point"
pca.gene.expression.ll_ld[1:12,1:5]
res.pca.ll_ld <- PCA(pca.gene.expression.ll_ld, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_pca_ind(res.pca.ll_ld, col.ind = rep(c("ZT00", "ZT04", "ZT08",  "ZT12", "ZT16", "ZT20"), each = 2),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)
fviz_pca_ind(res.pca.ll_ld, col.ind = rep(c("Day1", "Day2"),6),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)

# PCA of LL SD samples
pca.gene.expression.ll_sd <- data.frame(colnames(ll_sd.gene.expression),t(ll_sd.gene.expression))
colnames(pca.gene.expression.ll_sd)[1] <- "Time point"
pca.gene.expression.ll_sd[1:12,1:5]
res.pca.ll_sd <- PCA(pca.gene.expression.ll_sd, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_pca_ind(res.pca.ll_sd, col.ind = rep(c("ZT00", "ZT04", "ZT08",  "ZT12", "ZT16", "ZT20"), each = 2),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)
fviz_pca_ind(res.pca.ll_sd, col.ind = rep(c("Day1", "Day2"),6),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)


# PCA of DD LD samples
pca.gene.expression.dd_ld <- data.frame(colnames(dd.ld.gene.expression),t(dd.ld.gene.expression))
colnames(pca.gene.expression.dd_ld)[1] <- "Time point"
pca.gene.expression.dd_ld[1:12,1:5]
res.pca.dd_ld <- PCA(pca.gene.expression.dd_ld, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_pca_ind(res.pca.dd_ld, col.ind = rep(c("ZT00", "ZT04", "ZT08",  "ZT12", "ZT16", "ZT20"), each = 2),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)

fviz_pca_ind(res.pca.dd_ld, col.ind = rep(c("Day1", "Day2"),6),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)


# PCA of DD SD samples
pca.gene.expression.dd_sd <- data.frame(colnames(dd.sd.gene.expression),t(dd.sd.gene.expression))
colnames(pca.gene.expression.dd_sd)[1] <- "Time point"
pca.gene.expression.dd_sd[1:12,1:5]
res.pca.dd_sd <- PCA(pca.gene.expression.dd_sd, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_pca_ind(res.pca.dd_sd, col.ind = rep(c("ZT00", "ZT04", "ZT08",  "ZT12", "ZT16", "ZT20"), each = 2),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)

fviz_pca_ind(res.pca.dd_sd, col.ind = rep(c("Day1", "Day2"),6),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)


# Correct batch effect
browseVignettes("sva")


# Clustering of LD and SD samples only for rhytmic genes
sd.rhythmic <- read.table("sd.rhytmic.table.tsv", header = T)$V1
ld.rhythmic <- read.table("ld.rhytmic.table.tsv", header = T)$V1
ld.sd.gene.expression.rhythmic <- ld.sd.gene.expression[unique(c(ld.rhythmic,sd.rhythmic)),]
dim(ld.sd.gene.expression.rhythmic)
pca.gene.expression.ld.sd.rhythmic <- data.frame(colnames(ld.sd.gene.expression.rhythmic),t(ld.sd.gene.expression.rhythmic))
dim(pca.gene.expression.ld.sd.rhythmic)

colnames(pca.gene.expression.ld.sd.rhythmic)[1] <- "Time point"
res.pca.ld.sd.rhythmic <- PCA(pca.gene.expression.ld.sd.rhythmic, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
res.hcpc.ld.sd.rhythmic <- HCPC(res.pca.ld.sd.rhythmic, graph=FALSE)
fviz_dend(res.hcpc.ld.sd.rhythmic,k=2,cex = 1,palette = "jco",rect = TRUE, rect_fill = TRUE,rect_border = "jco",type="rectangle",
          labels_track_height = 500)

# PCA of rhythmic genes LD
ld.rhythmic <- read.table("ld.rhytmic.table.tsv", header = T)$V1
ld.gene.expression <- genes.expression[,ld.zt.i]
ld.gene.expression.rhythmic <- ld.gene.expression[ld.rhythmic,]
dim(ld.gene.expression.rhythmic)
pca.gene.expression.ld.rhythmic <- data.frame(colnames(ld.gene.expression.rhythmic),t(ld.gene.expression.rhythmic))
colnames(pca.gene.expression.ld.rhythmic)[1] <- "Time point"
res.pca.ld.rhythmic  <- PCA(pca.gene.expression.ld.rhythmic, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
png("pca_ld_rhythmic.png", width = 1500, height = 1500, res = 300)
fviz_pca_ind(res.pca.ld.rhythmic, col.ind = c(rep("ZT00", 3),rep("ZT04", 3), rep("ZT08", 3), 
                                     rep("ZT12", 3), rep("ZT16", 3), rep("ZT20", 3)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)
dev.off()

res.hcpc.ld.rhythmic <- HCPC(res.pca.ld.rhythmic, graph=FALSE)
png("hcpc_ld_rhythmic.png", width = 1500, height = 1500, res = 300)
fviz_dend(res.hcpc.ld.rhythmic,k=2,cex = 1,palette = "jco",rect = TRUE, rect_fill = TRUE,rect_border = "jco",type="rectangle",
          labels_track_height = 500)
dev.off()

# PCA of rhythmic genes SD
sd.rhythmic <- read.table("sd.rhytmic.table.tsv", header = T)$V1
sd.gene.expression <- genes.expression[,sd.zt.i]
sd.gene.expression.rhythmic <- sd.gene.expression[sd.rhythmic,]
dim(sd.gene.expression.rhythmic)
pca.gene.expression.sd.rhythmic <- data.frame(colnames(sd.gene.expression.rhythmic),t(sd.gene.expression.rhythmic))
colnames(pca.gene.expression.sd.rhythmic)[1] <- "Time point"
res.pca.sd.rhythmic  <- PCA(pca.gene.expression.sd.rhythmic, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
png("pca_sd_rhythmic.png", width = 1500, height = 1500, res = 300)
fviz_pca_ind(res.pca.sd.rhythmic, col.ind = c(rep("ZT00", 3),rep("ZT04", 3), rep("ZT08", 3), 
                                              rep("ZT12", 3), rep("ZT16", 3), rep("ZT20", 3)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)
dev.off()

# PCA of rhythmic genes LD 2 days
ld.rhythmic.2 <- read.table("ld.rhytmic.table.2.tsv", header = T)$V1
ld.zt.i <- sapply(X = ld.zt,FUN = function(x){ paste(x,1:2,sep="_")})
ld.gene.expression.2 <- genes.expression[,ld.zt.i]
ld.gene.expression.rhythmic.2 <- ld.gene.expression.2[ld.rhythmic.2,]
dim(ld.gene.expression.rhythmic.2)
pca.gene.expression.ld.rhythmic.2 <- data.frame(colnames(ld.gene.expression.rhythmic.2),t(ld.gene.expression.rhythmic.2))
colnames(pca.gene.expression.ld.rhythmic.2)[1] <- "Time point"
res.pca.ld.rhythmic.2  <- PCA(pca.gene.expression.ld.rhythmic.2, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
png("pca_ld_rhythmic_2_days.png", width = 1500, height = 1500, res = 300)
fviz_pca_ind(res.pca.ld.rhythmic.2, col.ind = c(rep("ZT00", 2),rep("ZT04", 2), rep("ZT08", 2), 
                                              rep("ZT12", 2), rep("ZT16", 2), rep("ZT20", 2)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)
dev.off()

# PCA of rhythmic genes SD 2 days
sd.rhythmic.2 <- read.table("sd.rhytmic.table.2.tsv", header = T)$V1
sd.zt.i <- sapply(X = sd.zt,FUN = function(x){ paste(x,1:2,sep="_")})
sd.gene.expression.2 <- genes.expression[,sd.zt.i]
sd.gene.expression.rhythmic.2 <- sd.gene.expression.2[sd.rhythmic.2,]
dim(sd.gene.expression.rhythmic.2)
pca.gene.expression.sd.rhythmic.2 <- data.frame(colnames(sd.gene.expression.rhythmic.2),t(sd.gene.expression.rhythmic.2))
colnames(pca.gene.expression.sd.rhythmic.2)[1] <- "Time point"
res.pca.sd.rhythmic.2  <- PCA(pca.gene.expression.sd.rhythmic.2, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
png("pca_sd_rhythmic_2_days.png", width = 1500, height = 1500, res = 300)
fviz_pca_ind(res.pca.sd.rhythmic.2, col.ind = c(rep("ZT00", 2),rep("ZT04", 2), rep("ZT08", 2), 
                                                rep("ZT12", 2), rep("ZT16", 2), rep("ZT20", 2)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)
dev.off()

res.hcpc.sd.rhythmic <- HCPC(res.pca.sd.rhythmic, graph=FALSE)
png("hcpc_sd_rhythmic.png", width = 1500, height = 1500, res = 300)
fviz_dend(res.hcpc.sd.rhythmic,k=2,cex = 1,palette = "jco",rect = TRUE, rect_fill = TRUE,rect_border = "jco",type="rectangle",
          labels_track_height = 500)
dev.off()

# Expression boxplot for all samples
head(ld.sd.gene.expression)
png("log2_boxplot_fpkm_per_sample.png", width = 2500, height = 1500, res=300)
boxplot(log2(ld.sd.gene.expression + 1),col=rainbow(ncol(ld.sd.gene.expression)),ylab="log2(FPKM + 1)",cex.lab=1.5,las=2)
dev.off()

# Normalyze data
library("NormalyzerDE")

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ld.zt.i <- sapply(X = ld.zt,FUN = function(x){ paste(x,1:3,sep="_")})
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt.i <- sapply(X = sd.zt,FUN = function(x){ paste(x,1:3,sep="_")})

ld.sd.gene.expression <- genes.expression[,c(ld.zt.i,sd.zt.i)]
groups_design <- paste0(c(rep("ld_",18), rep("sd_",18)), 
                        sapply(strsplit(colnames(ld.sd.gene.expression),"_"), function(x) x[[2]]))
design <- data.frame(sample=colnames(ld.sd.gene.expression),
                     group=groups_design)
ld.sd.gene.expression.normalyzer <- data.frame(gene=rownames(ld.sd.gene.expression), ld.sd.gene.expression)
head(ld.sd.gene.expression.normalyzer)

write.table(x = design,file = "normalyzer_design.tsv",quote = F,row.names = F,
            sep = "\t")

rows_above_one <- apply(ld.sd.gene.expression.normalyzer[, -1], 1, min) > 1 
df_above_one <- ld.sd.gene.expression.normalyzer[rows_above_one, ]
write.table(x = df_above_one,file = "normalyzer_expression.tsv",quote = F,row.names = F,
            sep = "\t")

normalyzer(jobName = "normalization_mp",designPath = "normalyzer_design.tsv",
           dataPath = "normalyzer_expression.tsv",outputDir = "." )


# Normalized gene expression
# VSN is the best one
vsn.normalized.data <- read.table(file="normalization_mp/Quantile-normalized.txt",
                                  header = T,as.is=T,sep="\t", row.names = 1)
head(vsn.normalized.data)
boxplot(vsn.normalized.data,outline=F,las=2,col=c(rep("red",8), rep("blue",6)),
        main="Datos normalizados", cex.main=1.5,ylab="Expresión génica")
boxplot(vsn.normalized.data,col=rainbow(ncol(vsn.normalized.data)),ylab="VSN Expression",cex.lab=1.5,las=2)

dim(vsn.normalized.data)

ld.gene.expression <- vsn.normalized.data[,ld.zt.i]
dim(ld.gene.expression)
pca.gene.expression.ld <- data.frame(colnames(ld.gene.expression),t(ld.gene.expression))
colnames(pca.gene.expression.ld)[1] <- "Time point"
pca.gene.expression.ld[1:12,1:5]
res.pca.ld <- PCA(pca.gene.expression.ld, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_pca_ind(res.pca.ld, col.ind = c(rep("ZT00", 3),rep("ZT04", 3), rep("ZT08", 3), 
                                     rep("ZT12", 3), rep("ZT16", 3), rep("ZT20", 3)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)

sd.gene.expression <- vsn.normalized.data[,sd.zt.i]
dim(sd.gene.expression)
pca.gene.expression.sd <- data.frame(colnames(sd.gene.expression),t(sd.gene.expression))
colnames(pca.gene.expression.sd)[1] <- "Time point"
pca.gene.expression.sd[1:12,1:5]
res.pca.sd <- PCA(pca.gene.expression.sd, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_pca_ind(res.pca.sd, col.ind = c(rep("ZT00", 3),rep("ZT04", 3), rep("ZT08", 3), 
                                     rep("ZT12", 3), rep("ZT16", 3), rep("ZT20", 3)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ld.zt.i <- sapply(X = ld.zt,FUN = function(x){ paste(x,1:3,sep="_")})
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt.i <- sapply(X = sd.zt,FUN = function(x){ paste(x,1:3,sep="_")})

# Add genes not represented in the normalization
gene.expression.comp <- read.table(file ="marchantia_gene_expression.tsv",sep="\t",header=T,as.is=T)
head(gene.expression.comp)

gene.expression.comp.ld.sd <- gene.expression.comp[,c(ld.zt.i, sd.zt.i)]

added.expression <- gene.expression.comp.ld.sd[which(!(rownames(gene.expression.comp.ld.sd) %in% rownames(vsn.normalized.data))),]

genes.expression <- rbind(vsn.normalized.data, added.expression)
dim(genes.expression)

#LD
library(rain)
ld.gene.expression <- genes.expression[,ld.zt.i]

results.ld <- rain(t(ld.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.ld$pVal < 0.05)/number.genes
rhythmic.genes.ld <- rownames(subset(results.ld, pVal < 0.05))
length(rhythmic.genes.ld)

results.ld.12 <- rain(t(ld.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.ld.12$pVal < 0.05)/number.genes
rhythmic.genes.ld.12 <- rownames(subset(results.ld.12, pVal < 0.05))
length(rhythmic.genes.ld.12)

complete.ld.rhythmic.genes <- unique(c(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(complete.ld.rhythmic.genes)/number.genes

non.rhythmic.ld.genes <- setdiff(mp.genes,complete.ld.rhythmic.genes)
length(non.rhythmic.ld.genes)
length(non.rhythmic.ld.genes)/number.genes
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))/number.genes

#SD
sd.gene.expression <- genes.expression[,sd.zt.i]
results.sd <- rain(t(sd.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.sd$pVal < 0.05)/number.genes
rhythmic.genes.sd <- rownames(subset(results.sd, pVal < 0.05))
length(rhythmic.genes.sd)

results.sd.12 <- rain(t(sd.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.sd.12$pVal < 0.05)/number.genes
rhythmic.genes.sd.12 <- rownames(subset(results.sd.12, pVal < 0.05))
length(rhythmic.genes.sd.12)

complete.sd.rhythmic.genes <- unique(c(rhythmic.genes.sd,rhythmic.genes.sd.12))
length(complete.sd.rhythmic.genes)
length(complete.sd.rhythmic.genes)/number.genes

non.rhythmic.sd.genes <- setdiff(mp.genes,complete.sd.rhythmic.genes)
length(non.rhythmic.sd.genes)
length(non.rhythmic.sd.genes)/number.genes

#Rhythmic genes in SD conditions but not in LD
sd.specific.rhythmic.genes <- setdiff(complete.sd.rhythmic.genes,complete.ld.rhythmic.genes)
length(sd.specific.rhythmic.genes)

#Rhythmic genes in LD conditions but not in SD
ld.specific.rhytmic.genes <-setdiff(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(ld.specific.rhytmic.genes)

#Rhythmic genes in both LD and SD
complete.ld.sd.rhythmic.genes<-intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(complete.ld.sd.rhythmic.genes)


#Save lists of rhythmic genes
ld.rhytmic.table<-matrix(data = complete.ld.rhythmic.genes,ncol = 1)
write.table(ld.rhytmic.table,file="ld.rhytmic.table.join.norm.quan.tsv")
sd.rhytmic.table<-matrix(data = complete.sd.rhythmic.genes,ncol = 1)
write.table(sd.rhytmic.table,file="sd.rhytmic.table.join.norm.quan.tsv")

#### Plots of rhythmic genes

#Venn diagram of LD and SD rhythmic genes
library(VennDiagram)
ld.rhythmic.join <- read.table("ld.rhytmic.table.join.norm.quan.tsv", header = T)$V1
sd.rhythmic.join <- read.table("sd.rhytmic.table.join.norm.quan.tsv", header = T)$V1
ld.rhythmic.2 <- read.table("ld.rhytmic.table.2.tsv", header = T)$V1
sd.rhythmic.2 <- read.table("sd.rhytmic.table.2.tsv", header = T)$V1

png("venn_2_vs_join_norm_ld.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(ld.rhythmic.join), area2 = length(ld.rhythmic.2),
                   cross.area = length(intersect(ld.rhythmic.join,ld.rhythmic.2)),
                   lwd = 3,category = c("LD_norm","LD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

png("venn_2_vs_join_norm_sd.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(sd.rhythmic.join), area2 = length(sd.rhythmic.2),
                   cross.area = length(intersect(sd.rhythmic.join,sd.rhythmic.2)),
                   lwd = 3,category = c("SD_norm","SD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()



# SVA
# https://support.bioconductor.org/p/9153697/

# Prepare data (cycle)
library(sva)
sva.gene.expression <- genes.expression[,c(ld.zt.i, sd.zt.i)]
pheno <- data.frame(sample=1:ncol(sva.gene.expression),
                    photoperiod=rep(c("LD","SD"), each=18),
                    zt=rep(rep(c("ZT00","ZT04","ZT08","ZT12","ZT16","ZT20"), each=3), 2),
                    batch=c(1,1,2,1,1,1,rep(c(1,1,2),10))
                    )
rownames(pheno) <- colnames(sva.gene.expression)
head(sva.gene.expression)

batch <- pheno$batch
edata <- sva.gene.expression

modcombat = model.matrix(~1, data=pheno)
modcombat2 = model.matrix(~ as.factor(photoperiod) + as.factor(zt), data=pheno)

combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
write.table(combat_edata,file = "combat_gene_expression.tsv",sep="\t") # Better without covariates

# RAIN
library(rain)

# Read data
genes.expression <- read.table(file ="combat_gene_expression.tsv",sep="\t",header=T,as.is=T)
dim(genes.expression)
head(genes.expression)
genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)

#LD
ld.gene.expression <- genes.expression[,ld.zt.i]

results.ld <- rain(t(ld.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.ld$pVal < 0.05)/number.genes
rhythmic.genes.ld <- rownames(subset(results.ld, pVal < 0.05))
length(rhythmic.genes.ld)

results.ld.12 <- rain(t(ld.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.ld.12$pVal < 0.05)/number.genes
rhythmic.genes.ld.12 <- rownames(subset(results.ld.12, pVal < 0.05))
length(rhythmic.genes.ld.12)

complete.ld.rhythmic.genes <- unique(c(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(complete.ld.rhythmic.genes)/number.genes

non.rhythmic.ld.genes <- setdiff(mp.genes,complete.ld.rhythmic.genes)
length(non.rhythmic.ld.genes)
length(non.rhythmic.ld.genes)/number.genes
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))/number.genes

#SD
sd.gene.expression <- genes.expression[,sd.zt.i]
results.sd <- rain(t(sd.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.sd$pVal < 0.05)/number.genes
rhythmic.genes.sd <- rownames(subset(results.sd, pVal < 0.05))
length(rhythmic.genes.sd)

results.sd.12 <- rain(t(sd.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.sd.12$pVal < 0.05)/number.genes
rhythmic.genes.sd.12 <- rownames(subset(results.sd.12, pVal < 0.05))
length(rhythmic.genes.sd.12)

complete.sd.rhythmic.genes <- unique(c(rhythmic.genes.sd,rhythmic.genes.sd.12))
length(complete.sd.rhythmic.genes)
length(complete.sd.rhythmic.genes)/number.genes

non.rhythmic.sd.genes <- setdiff(mp.genes,complete.sd.rhythmic.genes)
length(non.rhythmic.sd.genes)
length(non.rhythmic.sd.genes)/number.genes

#Rhythmic genes in SD conditions but not in LD
sd.specific.rhythmic.genes <- setdiff(complete.sd.rhythmic.genes,complete.ld.rhythmic.genes)
length(sd.specific.rhythmic.genes)

#Rhythmic genes in LD conditions but not in SD
ld.specific.rhytmic.genes <-setdiff(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(ld.specific.rhytmic.genes)

#Rhythmic genes in both LD and SD
complete.ld.sd.rhythmic.genes<-intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(complete.ld.sd.rhythmic.genes)


#Save lists of rhythmic genes
ld.rhytmic.table<-matrix(data = complete.ld.rhythmic.genes,ncol = 1)
write.table(ld.rhytmic.table,file="ld.rhytmic.table.combat.tsv")
sd.rhytmic.table<-matrix(data = complete.sd.rhythmic.genes,ncol = 1)
write.table(sd.rhytmic.table,file="sd.rhytmic.table.combat.tsv")

#### Plots of rhythmic genes

#Venn diagram of LD and SD rhythmic genes
library(VennDiagram)
ld.rhythmic.combat<- read.table("ld.rhytmic.table.combat.tsv", header = T)$V1
sd.rhythmic.combat <- read.table("sd.rhytmic.table.combat.tsv", header = T)$V1
ld.rhythmic.2 <- read.table("ld.rhytmic.table.2.tsv", header = T)$V1
sd.rhythmic.2 <- read.table("sd.rhytmic.table.2.tsv", header = T)$V1

png("venn_2_vs_combat_ld.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(ld.rhythmic.combat), area2 = length(ld.rhythmic.2),
                   cross.area = length(intersect(ld.rhythmic.combat,ld.rhythmic.2)),
                   lwd = 3,category = c("LD_combat","LD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

png("venn_2_vs_combat_sd.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(sd.rhythmic.combat), area2 = length(sd.rhythmic.2),
                   cross.area = length(intersect(sd.rhythmic.combat,sd.rhythmic.2)),
                   lwd = 3,category = c("SD_combat","SD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()


### Change only LD ZT08 Day 3
genes.expression <- read.table(file ="marchantia_gene_expression.tsv",sep="\t",header=T,as.is=T)
dim(genes.expression)
head(genes.expression)
genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ld.zt.i <- sapply(X = ld.zt,FUN = function(x){ paste(x,1:3,sep="_")})
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt.i <- sapply(X = sd.zt,FUN = function(x){ paste(x,1:3,sep="_")})

ld.gene.expression <- genes.expression[,ld.zt.i]
sd.gene.expression <- genes.expression[,sd.zt.i]

ld.sd.gene.expression <- genes.expression[,c(ld.zt.i,sd.zt.i)]
ld.gene.expression[,9] <- (ld.gene.expression[,7] + ld.gene.expression[,8])/2
ld.sd.gene.expression[,9] <- (ld.sd.gene.expression[,7] + ld.sd.gene.expression[,8])/2

#LD
library(rain)
results.ld <- rain(t(ld.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.ld$pVal < 0.05)/number.genes
rhythmic.genes.ld <- rownames(subset(results.ld, pVal < 0.05))
length(rhythmic.genes.ld)

results.ld.12 <- rain(t(ld.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.ld.12$pVal < 0.05)/number.genes
rhythmic.genes.ld.12 <- rownames(subset(results.ld.12, pVal < 0.05))
length(rhythmic.genes.ld.12)

complete.ld.rhythmic.genes <- unique(c(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(complete.ld.rhythmic.genes)/number.genes

non.rhythmic.ld.genes <- setdiff(mp.genes,complete.ld.rhythmic.genes)
length(non.rhythmic.ld.genes)
length(non.rhythmic.ld.genes)/number.genes
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))/number.genes

#Save lists of rhythmic genes
ld.rhytmic.table<-matrix(data = complete.ld.rhythmic.genes,ncol = 1)
write.table(ld.rhytmic.table,file="ld.rhytmic.table.ld8.day3.changed.tsv")

#Venn diagram of LD rhythmic genes
library(VennDiagram)
ld.rhythmic.ld8<- read.table("ld.rhytmic.table.ld8.day3.changed.tsv", header = T)$V1
ld.rhythmic.2 <- read.table("ld.rhytmic.table.2.tsv", header = T)$V1

png("venn_2_vs_ld8_changed.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(ld.rhythmic.ld8), area2 = length(ld.rhythmic.2),
                   cross.area = length(intersect(ld.rhythmic.ld8,ld.rhythmic.2)),
                   lwd = 3,category = c("LD_8_changed","LD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()


### # Intragroup normalization
# Separate all groups
genes.expression <- read.table(file ="marchantia_gene_expression.tsv",sep="\t",header=T,as.is=T)
dim(genes.expression)
head(genes.expression)
genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)

# Number of non-expressed genes and low-expressed genes
length(which(apply(X = genes.expression,MARGIN = 1,FUN = max) == 0))
length(which(apply(X = genes.expression,MARGIN = 1,FUN = max) < 1))

# Gene expression matrix for samples under light-dark cycles
ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ld.zt.i <- sapply(X = ld.zt,FUN = function(x){ paste(x,1:3,sep="_")})
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt.i <- sapply(X = sd.zt,FUN = function(x){ paste(x,1:3,sep="_")})

ld.gene.expression <- genes.expression[,ld.zt.i]
sd.gene.expression <- genes.expression[,sd.zt.i]

ld.sd.gene.expression <- genes.expression[,c(ld.zt.i,sd.zt.i)]
ld.sd.gene.expression.normalyzer <- data.frame(gene=rownames(ld.sd.gene.expression), ld.sd.gene.expression)
rows_above_one <- apply(ld.sd.gene.expression.normalyzer[, -1], 1, min) > 1 
df_above_one <- ld.sd.gene.expression.normalyzer[rows_above_one, ]


# Prepare files for iterative NormalyzerDE
for (i in 1:length(c(ld.zt, sd.zt)))
{
  df_above_one_iter <- df_above_one[,-1]
  design <- data.frame(sample=colnames(df_above_one_iter[,((i-1)*3+1):((i-1)*3+3)]),
                       group=1)
  expression.norm <- cbind(genes=df_above_one[,1] , df_above_one_iter[,((i-1)*3+1):((i-1)*3+3)])
  
  write.table(x = design,file = paste0(c(ld.zt, sd.zt)[i],"_normalyzer_design.tsv"),quote = F,row.names = F,
              sep = "\t")
  
  write.table(x = expression.norm,file = paste0(c(ld.zt, sd.zt)[i],"_normalyzer_expression.tsv"),quote = F,row.names = F,
              sep = "\t")
  
}

ld.sd.gene.expression <- genes.expression[,c(ld.zt.i,sd.zt.i)]
groups_design <- paste0(c(rep("ld_",18), rep("sd_",18)), 
                        sapply(strsplit(colnames(ld.sd.gene.expression),"_"), function(x) x[[2]]))
design <- data.frame(sample=colnames(ld.sd.gene.expression),
                     group=groups_design)
ld.sd.gene.expression.normalyzer <- data.frame(gene=rownames(ld.sd.gene.expression), ld.sd.gene.expression)
head(ld.sd.gene.expression.normalyzer)

write.table(x = design,file = "normalyzer_design.tsv",quote = F,row.names = F,
            sep = "\t")


write.table(x = df_above_one,file = "normalyzer_expression.tsv",quote = F,row.names = F,
            sep = "\t")

normalyzer(jobName = "normalization_mp",designPath = "normalyzer_design.tsv",
           dataPath = "normalyzer_expression.tsv",outputDir = "." )


# Read normalyze data
expression_folders <- list.dirs("normalyzed_per_groups", full.names = F, recursive = F)
norm_list <- lapply(expression_folders, function(x) read.table(paste0("normalyzed_per_groups/", x, "/VSN-normalized.txt"), header = T))
sapply(norm_list, dim)

norm_joined_table <- data.frame(genes=norm_list[[1]]$genes)

for (i in 1:length(norm_list))
{
  norm_joined_table <- cbind(norm_joined_table, norm_list[[i]][,-1])
}

head(norm_joined_table)  
write.table(norm_joined_table, "separated_groups_norm_expression.tsv", sep="\t")

# RAIN
library(rain)

# Read data
genes.expression <- read.table(file ="separated_groups_norm_expression.tsv",sep="\t",header=T,as.is=T)
rownames(genes.expression) <- genes.expression$genes
genes.expression <- genes.expression[,-1]
dim(genes.expression)
head(genes.expression)

genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ld.zt.i <- sapply(X = ld.zt,FUN = function(x){ paste(x,1:3,sep="_")})
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt.i <- sapply(X = sd.zt,FUN = function(x){ paste(x,1:3,sep="_")})

# Add genes not represented in the normalization
gene.expression.comp <- read.table(file ="marchantia_gene_expression.tsv",sep="\t",header=T,as.is=T)
head(gene.expression.comp)

gene.expression.comp.ld.sd <- gene.expression.comp[,c(ld.zt.i, sd.zt.i)]

added.expression <- gene.expression.comp.ld.sd[which(!(rownames(gene.expression.comp.ld.sd) %in% rownames(genes.expression))),]

genes.expression <- rbind(genes.expression, added.expression)
dim(genes.expression)

#LD
ld.gene.expression <- genes.expression[,ld.zt.i]

results.ld <- rain(t(ld.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.ld$pVal < 0.05)/number.genes
rhythmic.genes.ld <- rownames(subset(results.ld, pVal < 0.05))
length(rhythmic.genes.ld)

results.ld.12 <- rain(t(ld.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.ld.12$pVal < 0.05)/number.genes
rhythmic.genes.ld.12 <- rownames(subset(results.ld.12, pVal < 0.05))
length(rhythmic.genes.ld.12)

complete.ld.rhythmic.genes <- unique(c(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(complete.ld.rhythmic.genes)/number.genes

non.rhythmic.ld.genes <- setdiff(mp.genes,complete.ld.rhythmic.genes)
length(non.rhythmic.ld.genes)
length(non.rhythmic.ld.genes)/number.genes
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))/number.genes

#SD
sd.gene.expression <- genes.expression[,sd.zt.i]
results.sd <- rain(t(sd.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.sd$pVal < 0.05)/number.genes
rhythmic.genes.sd <- rownames(subset(results.sd, pVal < 0.05))
length(rhythmic.genes.sd)

results.sd.12 <- rain(t(sd.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.sd.12$pVal < 0.05)/number.genes
rhythmic.genes.sd.12 <- rownames(subset(results.sd.12, pVal < 0.05))
length(rhythmic.genes.sd.12)

complete.sd.rhythmic.genes <- unique(c(rhythmic.genes.sd,rhythmic.genes.sd.12))
length(complete.sd.rhythmic.genes)
length(complete.sd.rhythmic.genes)/number.genes

non.rhythmic.sd.genes <- setdiff(mp.genes,complete.sd.rhythmic.genes)
length(non.rhythmic.sd.genes)
length(non.rhythmic.sd.genes)/number.genes

#Rhythmic genes in SD conditions but not in LD
sd.specific.rhythmic.genes <- setdiff(complete.sd.rhythmic.genes,complete.ld.rhythmic.genes)
length(sd.specific.rhythmic.genes)

#Rhythmic genes in LD conditions but not in SD
ld.specific.rhytmic.genes <-setdiff(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(ld.specific.rhytmic.genes)

#Rhythmic genes in both LD and SD
complete.ld.sd.rhythmic.genes<-intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(complete.ld.sd.rhythmic.genes)


#Save lists of rhythmic genes
ld.rhytmic.table<-matrix(data = complete.ld.rhythmic.genes,ncol = 1)
write.table(ld.rhytmic.table,file="ld.rhytmic.table.sep.norm.tsv")
sd.rhytmic.table<-matrix(data = complete.sd.rhythmic.genes,ncol = 1)
write.table(sd.rhytmic.table,file="sd.rhytmic.table.sep.norm.tsv")

#### Plots of rhythmic genes

#Venn diagram of LD and SD rhythmic genes
library(VennDiagram)
ld.rhythmic.norm<- read.table("ld.rhytmic.table.sep.norm.tsv", header = T)$V1
sd.rhythmic.norm <- read.table("sd.rhytmic.table.sep.norm.tsv", header = T)$V1
ld.rhythmic.2 <- read.table("ld.rhytmic.table.2.tsv", header = T)$V1
sd.rhythmic.2 <- read.table("sd.rhytmic.table.2.tsv", header = T)$V1

png("venn_2_vs_sep_norm_ld.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(ld.rhythmic.norm), area2 = length(ld.rhythmic.2),
                   cross.area = length(intersect(ld.rhythmic.norm,ld.rhythmic.2)),
                   lwd = 3,category = c("LD_norm","LD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

png("venn_2_vs_sep_norm_sd.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(sd.rhythmic.norm), area2 = length(sd.rhythmic.2),
                   cross.area = length(intersect(sd.rhythmic.norm,sd.rhythmic.2)),
                   lwd = 3,category = c("SD_norm","SD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()


# PCA tras combat
# Clustering of LD and SD samples only for rhytmic genes
genes.expression <- read.table(file ="combat_gene_expression.tsv",sep="\t",header=T,as.is=T)
dim(genes.expression)
head(genes.expression)
genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ld.zt.i <- sapply(X = ld.zt,FUN = function(x){ paste(x,1:3,sep="_")})
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt.i <- sapply(X = sd.zt,FUN = function(x){ paste(x,1:3,sep="_")})

ld.sd.gene.expression <- genes.expression[,c(ld.zt.i,sd.zt.i)]
head(ld.sd.gene.expression)
dim(ld.sd.gene.expression)

sd.rhythmic <- read.table("sd.rhytmic.table.combat.tsv", header = T)$V1
ld.rhythmic <- read.table("ld.rhytmic.table.combat.tsv", header = T)$V1
ld.sd.gene.expression.rhythmic <- ld.sd.gene.expression[unique(c(ld.rhythmic,sd.rhythmic)),]
dim(ld.sd.gene.expression.rhythmic)
pca.gene.expression.ld.sd.rhythmic <- data.frame(colnames(ld.sd.gene.expression.rhythmic),t(ld.sd.gene.expression.rhythmic))
dim(pca.gene.expression.ld.sd.rhythmic)

library(factoextra)
library(FactoMineR)

colnames(pca.gene.expression.ld.sd.rhythmic)[1] <- "Time point"
res.pca.ld.sd.rhythmic <- PCA(pca.gene.expression.ld.sd.rhythmic, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
res.hcpc.ld.sd.rhythmic <- HCPC(res.pca.ld.sd.rhythmic, graph=FALSE)
fviz_dend(res.hcpc.ld.sd.rhythmic,k=2,cex = 1,palette = "jco",rect = TRUE, rect_fill = TRUE,rect_border = "jco",type="rectangle",
          labels_track_height = 500)

# PCA of rhythmic genes LD
ld.gene.expression <- genes.expression[,ld.zt.i]
ld.gene.expression.rhythmic <- ld.gene.expression[ld.rhythmic,]
dim(ld.gene.expression.rhythmic)
pca.gene.expression.ld.rhythmic <- data.frame(colnames(ld.gene.expression.rhythmic),t(ld.gene.expression.rhythmic))
colnames(pca.gene.expression.ld.rhythmic)[1] <- "Time point"
res.pca.ld.rhythmic  <- PCA(pca.gene.expression.ld.rhythmic, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
png("pca_ld_rhythmic_combat.png", width = 1500, height = 1500, res = 300)
fviz_pca_ind(res.pca.ld.rhythmic, col.ind = c(rep("ZT00", 3),rep("ZT04", 3), rep("ZT08", 3), 
                                              rep("ZT12", 3), rep("ZT16", 3), rep("ZT20", 3)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)
dev.off()

res.hcpc.ld.rhythmic <- HCPC(res.pca.ld.rhythmic, graph=FALSE)
png("hcpc_ld_rhythmic.png", width = 1500, height = 1500, res = 300)
fviz_dend(res.hcpc.ld.rhythmic,k=2,cex = 1,palette = "jco",rect = TRUE, rect_fill = TRUE,rect_border = "jco",type="rectangle",
          labels_track_height = 500)
dev.off()

# PCA of rhythmic genes SD
sd.gene.expression <- genes.expression[,sd.zt.i]
sd.gene.expression.rhythmic <- sd.gene.expression[sd.rhythmic,]
dim(sd.gene.expression.rhythmic)
pca.gene.expression.sd.rhythmic <- data.frame(colnames(sd.gene.expression.rhythmic),t(sd.gene.expression.rhythmic))
colnames(pca.gene.expression.sd.rhythmic)[1] <- "Time point"
res.pca.sd.rhythmic  <- PCA(pca.gene.expression.sd.rhythmic, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
png("pca_sd_rhythmic_combat.png", width = 1500, height = 1500, res = 300)
fviz_pca_ind(res.pca.sd.rhythmic, col.ind = c(rep("ZT00", 3),rep("ZT04", 3), rep("ZT08", 3), 
                                              rep("ZT12", 3), rep("ZT16", 3), rep("ZT20", 3)),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)
dev.off()

# Combat adjusted
# Prepare data (cycle)
library(sva)
sva.gene.expression <- genes.expression[,c(ld.zt.i, sd.zt.i)]
pheno <- data.frame(sample=1:ncol(sva.gene.expression),
                    photoperiod=rep(c("LD","SD"), each=18),
                    zt=rep(rep(c("ZT00","ZT04","ZT08","ZT12","ZT16","ZT20"), each=3), 2),
                    batch=c(1,1,2,1,1,1,rep(c(1,1,2),10))
)
rownames(pheno) <- colnames(sva.gene.expression)
head(sva.gene.expression)

batch <- pheno$batch
edata <- sva.gene.expression

modcombat = model.matrix(~1, data=pheno)
modcombat2 = model.matrix(~ as.factor(photoperiod) + as.factor(zt), data=pheno)

combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE, ref.batch = 1)
write.table(combat_edata,file = "combat_gene_expression_ref.tsv",sep="\t") # Better without covariates

# RAIN
library(rain)

# Read data
genes.expression <- read.table(file ="combat_gene_expression_ref.tsv",sep="\t",header=T,as.is=T)
dim(genes.expression)
head(genes.expression)
genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)

#LD
ld.gene.expression <- genes.expression[,ld.zt.i]

results.ld <- rain(t(ld.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.ld$pVal < 0.05)/number.genes
rhythmic.genes.ld <- rownames(subset(results.ld, pVal < 0.05))
length(rhythmic.genes.ld)

results.ld.12 <- rain(t(ld.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.ld.12$pVal < 0.05)/number.genes
rhythmic.genes.ld.12 <- rownames(subset(results.ld.12, pVal < 0.05))
length(rhythmic.genes.ld.12)

complete.ld.rhythmic.genes <- unique(c(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(complete.ld.rhythmic.genes)/number.genes

non.rhythmic.ld.genes <- setdiff(mp.genes,complete.ld.rhythmic.genes)
length(non.rhythmic.ld.genes)
length(non.rhythmic.ld.genes)/number.genes
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))/number.genes

#SD
sd.gene.expression <- genes.expression[,sd.zt.i]
results.sd <- rain(t(sd.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.sd$pVal < 0.05)/number.genes
rhythmic.genes.sd <- rownames(subset(results.sd, pVal < 0.05))
length(rhythmic.genes.sd)

results.sd.12 <- rain(t(sd.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.sd.12$pVal < 0.05)/number.genes
rhythmic.genes.sd.12 <- rownames(subset(results.sd.12, pVal < 0.05))
length(rhythmic.genes.sd.12)

complete.sd.rhythmic.genes <- unique(c(rhythmic.genes.sd,rhythmic.genes.sd.12))
length(complete.sd.rhythmic.genes)
length(complete.sd.rhythmic.genes)/number.genes

non.rhythmic.sd.genes <- setdiff(mp.genes,complete.sd.rhythmic.genes)
length(non.rhythmic.sd.genes)
length(non.rhythmic.sd.genes)/number.genes

#Rhythmic genes in SD conditions but not in LD
sd.specific.rhythmic.genes <- setdiff(complete.sd.rhythmic.genes,complete.ld.rhythmic.genes)
length(sd.specific.rhythmic.genes)

#Rhythmic genes in LD conditions but not in SD
ld.specific.rhytmic.genes <-setdiff(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(ld.specific.rhytmic.genes)

#Rhythmic genes in both LD and SD
complete.ld.sd.rhythmic.genes<-intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(complete.ld.sd.rhythmic.genes)


#Save lists of rhythmic genes
ld.rhytmic.table<-matrix(data = complete.ld.rhythmic.genes,ncol = 1)
write.table(ld.rhytmic.table,file="ld.rhytmic.table.combat.ref.tsv")
sd.rhytmic.table<-matrix(data = complete.sd.rhythmic.genes,ncol = 1)
write.table(sd.rhytmic.table,file="sd.rhytmic.table.combat.ref.tsv")

#### Plots of rhythmic genes

#Venn diagram of LD and SD rhythmic genes
library(VennDiagram)
ld.rhythmic.combat<- read.table("ld.rhytmic.table.combat.ref.tsv", header = T)$V1
sd.rhythmic.combat <- read.table("sd.rhytmic.table.combat.ref.tsv", header = T)$V1
ld.rhythmic.2 <- read.table("ld.rhytmic.table.2.tsv", header = T)$V1
sd.rhythmic.2 <- read.table("sd.rhytmic.table.2.tsv", header = T)$V1

png("venn_2_vs_combat_ref_ld.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(ld.rhythmic.combat), area2 = length(ld.rhythmic.2),
                   cross.area = length(intersect(ld.rhythmic.combat,ld.rhythmic.2)),
                   lwd = 3,category = c("LD_combat","LD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

png("venn_2_vs_combat_ref_sd.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(sd.rhythmic.combat), area2 = length(sd.rhythmic.2),
                   cross.area = length(intersect(sd.rhythmic.combat,sd.rhythmic.2)),
                   lwd = 3,category = c("SD_combat","SD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# Combat photoperiod as covariate
# Prepare data (cycle)
library(sva)
sva.gene.expression <- genes.expression[,c(ld.zt.i, sd.zt.i)]
pheno <- data.frame(sample=1:ncol(sva.gene.expression),
                    photoperiod=rep(c("LD","SD"), each=18),
                    zt=rep(rep(c("ZT00","ZT04","ZT08","ZT12","ZT16","ZT20"), each=3), 2),
                    batch=c(1,1,2,1,1,1,rep(c(1,1,2),10))
)
rownames(pheno) <- colnames(sva.gene.expression)
head(sva.gene.expression)

batch <- pheno$batch
edata <- sva.gene.expression

modcombat = model.matrix(~ as.factor(photoperiod), data=pheno)

combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE, ref.batch = 1)
write.table(combat_edata,file = "combat_gene_expression_covar.tsv",sep="\t") # Better without covariates

# RAIN
library(rain)

# Read data
genes.expression <- read.table(file ="combat_gene_expression_covar.tsv",sep="\t",header=T,as.is=T)
dim(genes.expression)
head(genes.expression)
genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)

#LD
ld.gene.expression <- genes.expression[,ld.zt.i]

results.ld <- rain(t(ld.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.ld$pVal < 0.05)/number.genes
rhythmic.genes.ld <- rownames(subset(results.ld, pVal < 0.05))
length(rhythmic.genes.ld)

results.ld.12 <- rain(t(ld.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.ld.12$pVal < 0.05)/number.genes
rhythmic.genes.ld.12 <- rownames(subset(results.ld.12, pVal < 0.05))
length(rhythmic.genes.ld.12)

complete.ld.rhythmic.genes <- unique(c(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(complete.ld.rhythmic.genes)/number.genes

non.rhythmic.ld.genes <- setdiff(mp.genes,complete.ld.rhythmic.genes)
length(non.rhythmic.ld.genes)
length(non.rhythmic.ld.genes)/number.genes
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))/number.genes

#SD
sd.gene.expression <- genes.expression[,sd.zt.i]
results.sd <- rain(t(sd.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.sd$pVal < 0.05)/number.genes
rhythmic.genes.sd <- rownames(subset(results.sd, pVal < 0.05))
length(rhythmic.genes.sd)

results.sd.12 <- rain(t(sd.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.sd.12$pVal < 0.05)/number.genes
rhythmic.genes.sd.12 <- rownames(subset(results.sd.12, pVal < 0.05))
length(rhythmic.genes.sd.12)

complete.sd.rhythmic.genes <- unique(c(rhythmic.genes.sd,rhythmic.genes.sd.12))
length(complete.sd.rhythmic.genes)
length(complete.sd.rhythmic.genes)/number.genes

non.rhythmic.sd.genes <- setdiff(mp.genes,complete.sd.rhythmic.genes)
length(non.rhythmic.sd.genes)
length(non.rhythmic.sd.genes)/number.genes

#Rhythmic genes in SD conditions but not in LD
sd.specific.rhythmic.genes <- setdiff(complete.sd.rhythmic.genes,complete.ld.rhythmic.genes)
length(sd.specific.rhythmic.genes)

#Rhythmic genes in LD conditions but not in SD
ld.specific.rhytmic.genes <-setdiff(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(ld.specific.rhytmic.genes)

#Rhythmic genes in both LD and SD
complete.ld.sd.rhythmic.genes<-intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(complete.ld.sd.rhythmic.genes)


#Save lists of rhythmic genes
ld.rhytmic.table<-matrix(data = complete.ld.rhythmic.genes,ncol = 1)
write.table(ld.rhytmic.table,file="ld.rhytmic.table.combat.covar.tsv")
sd.rhytmic.table<-matrix(data = complete.sd.rhythmic.genes,ncol = 1)
write.table(sd.rhytmic.table,file="sd.rhytmic.table.combat.covar.tsv")

#### Plots of rhythmic genes

#Venn diagram of LD and SD rhythmic genes
library(VennDiagram)
ld.rhythmic.combat<- read.table("ld.rhytmic.table.combat.covar.tsv", header = T)$V1
sd.rhythmic.combat <- read.table("sd.rhytmic.table.combat.covar.tsv", header = T)$V1
ld.rhythmic.2 <- read.table("ld.rhytmic.table.2.tsv", header = T)$V1
sd.rhythmic.2 <- read.table("sd.rhytmic.table.2.tsv", header = T)$V1

png("venn_2_vs_combat_covar_ld.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(ld.rhythmic.combat), area2 = length(ld.rhythmic.2),
                   cross.area = length(intersect(ld.rhythmic.combat,ld.rhythmic.2)),
                   lwd = 3,category = c("LD_combat","LD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

png("venn_2_vs_combat_covar_sd.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(sd.rhythmic.combat), area2 = length(sd.rhythmic.2),
                   cross.area = length(intersect(sd.rhythmic.combat,sd.rhythmic.2)),
                   lwd = 3,category = c("SD_combat","SD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()


# Combat separated
library(sva)
sva.gene.expression.sd <- genes.expression[,c(sd.zt.i)]
sva.gene.expression.ld <- genes.expression[,c(ld.zt.i)]

pheno.ld <- data.frame(sample=1:ncol(sva.gene.expression.ld),
                       photoperiod=rep(c("LD"), each=18),
                       zt=rep(c("ZT00","ZT04","ZT08","ZT12","ZT16","ZT20"), each=3),
                       batch=c(1,1,2,1,1,1,rep(c(1,1,2),4))
)

pheno.sd <- data.frame(sample=1:ncol(sva.gene.expression.sd),
                    photoperiod=rep(c("SD"), each=18),
                    zt=rep(c("ZT00","ZT04","ZT08","ZT12","ZT16","ZT20"), each=3),
                    batch=c(rep(c(1,1,2),6))
)



rownames(pheno.ld) <- colnames(sva.gene.expression.ld)
head(sva.gene.expression.ld)
batch.ld <- pheno.ld$batch
edata.ld <- sva.gene.expression.ld
modcombat.ld = model.matrix(~1, data=pheno.ld)
combat_edata.ld = ComBat(dat=edata.ld, batch=batch.ld, mod=modcombat.ld, par.prior=TRUE, prior.plots=FALSE, ref.batch = 1)
write.table(combat_edata.ld,file = "combat_gene_expression_ld.tsv",sep="\t")

rownames(pheno.sd) <- colnames(sva.gene.expression.sd)
head(sva.gene.expression.sd)
batch.sd <- pheno.sd$batch
edata.sd <- sva.gene.expression.sd
modcombat.sd = model.matrix(~1, data=pheno.sd)
combat_edata.sd = ComBat(dat=edata.sd, batch=batch.sd, mod=modcombat.sd, par.prior=TRUE, prior.plots=FALSE, ref.batch = 1)
write.table(combat_edata.sd,file = "combat_gene_expression_sd.tsv",sep="\t")


# RAIN
library(rain)

# Read data
genes.expression.ld <- read.table(file ="combat_gene_expression_ld.tsv",sep="\t",header=T,as.is=T)
genes.expression.sd <- read.table(file ="combat_gene_expression_sd.tsv",sep="\t",header=T,as.is=T)

sum(rownames(genes.expression.ld) == rownames(genes.expression.sd))

genes.expression <- cbind(genes.expression.ld, genes.expression.sd)
head(genes.expression)
genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)

#LD
ld.gene.expression <- genes.expression[,ld.zt.i]

results.ld <- rain(t(ld.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.ld$pVal < 0.05)/number.genes
rhythmic.genes.ld <- rownames(subset(results.ld, pVal < 0.05))
length(rhythmic.genes.ld)

results.ld.12 <- rain(t(ld.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.ld.12$pVal < 0.05)/number.genes
rhythmic.genes.ld.12 <- rownames(subset(results.ld.12, pVal < 0.05))
length(rhythmic.genes.ld.12)

complete.ld.rhythmic.genes <- unique(c(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(complete.ld.rhythmic.genes)/number.genes

non.rhythmic.ld.genes <- setdiff(mp.genes,complete.ld.rhythmic.genes)
length(non.rhythmic.ld.genes)
length(non.rhythmic.ld.genes)/number.genes
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))/number.genes

#SD
sd.gene.expression <- genes.expression[,sd.zt.i]
results.sd <- rain(t(sd.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.sd$pVal < 0.05)/number.genes
rhythmic.genes.sd <- rownames(subset(results.sd, pVal < 0.05))
length(rhythmic.genes.sd)

#results.sd[order(results.sd$pVal),]

results.sd.12 <- rain(t(sd.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.sd.12$pVal < 0.05)/number.genes
rhythmic.genes.sd.12 <- rownames(subset(results.sd.12, pVal < 0.05))
length(rhythmic.genes.sd.12)

complete.sd.rhythmic.genes <- unique(c(rhythmic.genes.sd,rhythmic.genes.sd.12))
length(complete.sd.rhythmic.genes)
length(complete.sd.rhythmic.genes)/number.genes

non.rhythmic.sd.genes <- setdiff(mp.genes,complete.sd.rhythmic.genes)
length(non.rhythmic.sd.genes)
length(non.rhythmic.sd.genes)/number.genes

#Rhythmic genes in SD conditions but not in LD
sd.specific.rhythmic.genes <- setdiff(complete.sd.rhythmic.genes,complete.ld.rhythmic.genes)
length(sd.specific.rhythmic.genes)

#Rhythmic genes in LD conditions but not in SD
ld.specific.rhytmic.genes <-setdiff(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(ld.specific.rhytmic.genes)

#Rhythmic genes in both LD and SD
complete.ld.sd.rhythmic.genes<-intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(complete.ld.sd.rhythmic.genes)


#Save lists of rhythmic genes
ld.rhytmic.table<-matrix(data = complete.ld.rhythmic.genes,ncol = 1)
write.table(ld.rhytmic.table,file="ld.rhytmic.table.combat.sep.tsv")
sd.rhytmic.table<-matrix(data = complete.sd.rhythmic.genes,ncol = 1)
write.table(sd.rhytmic.table,file="sd.rhytmic.table.combat.sep.tsv")

#### Plots of rhythmic genes

#Venn diagram of LD and SD rhythmic genes
library(VennDiagram)
ld.rhythmic.combat<- read.table("ld.rhytmic.table.combat.sep.tsv", header = T)$V1
sd.rhythmic.combat <- read.table("sd.rhytmic.table.combat.sep.tsv", header = T)$V1
ld.rhythmic.2 <- read.table("ld.rhytmic.table.2.tsv", header = T)$V1
sd.rhythmic.2 <- read.table("sd.rhytmic.table.2.tsv", header = T)$V1

png("venn_2_vs_combat_sep_ld.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(ld.rhythmic.combat), area2 = length(ld.rhythmic.2),
                   cross.area = length(intersect(ld.rhythmic.combat,ld.rhythmic.2)),
                   lwd = 3,category = c("LD_combat","LD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

png("venn_2_vs_combat_sep_sd.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(sd.rhythmic.combat), area2 = length(sd.rhythmic.2),
                   cross.area = length(intersect(sd.rhythmic.combat,sd.rhythmic.2)),
                   lwd = 3,category = c("SD_combat","SD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# Create table
genes.expression.ld <- read.table(file ="combat_gene_expression_ld.tsv",sep="\t",header=T,as.is=T)
genes.expression.sd <- read.table(file ="combat_gene_expression_sd.tsv",sep="\t",header=T,as.is=T)

sum(rownames(genes.expression.ld) == rownames(genes.expression.sd))

genes.expression <- cbind(genes.expression.ld, genes.expression.sd)
head(genes.expression)

write.table(genes.expression,file = "combat_gene_expression_ld_sd_sep.tsv",sep="\t")


genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)

# Pruebas

# Read data
genes.expression <- read.table(file ="combat_gene_expression_ld_sd_sep.tsv",sep="\t",header=T,as.is=T)
dim(genes.expression)
head(genes.expression)
genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)


# Gene expression matrix for samples under light-dark cycles
ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ld.zt.i <- sapply(X = ld.zt,FUN = function(x){ paste(x,1:3,sep="_")})
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt.i <- sapply(X = sd.zt,FUN = function(x){ paste(x,1:3,sep="_")})

ld.sd.gene.expression <- genes.expression[,c(ld.zt.i,sd.zt.i)]
head(ld.sd.gene.expression)
dim(ld.sd.gene.expression)

ld.gene.expression <- genes.expression[,c(ld.zt.i)]
sd.gene.expression <- genes.expression[,c(sd.zt.i)]

plot.ld.sd(gene.id = "Mp3g09460",gene.name = "MpTOC1", ld.gene.expression,sd.gene.expression)

#Examples taken from a paper
paper.ids <- c("Mp3g09460","Mp2g16560","Mp1g10870","Mp1g16510","Mp1g13750","Mp3g14360","Mp2g14310","Mp4g21220","Mp3g14360")
paper.names <- c("MpTOC1","MpPRR","MpELF3","MpELF4","MpGI","MpZTL","MpRVE","MpLUX","MpFKF")

for (i in 1:length(paper.ids))
{
  png(paste0(paper.names[i],".png"), width = 2000, height = 1500, res=300)
  plot.ld.sd(gene.id = paper.ids[i],gene.name = paper.names[i], ld.gene.expression,sd.gene.expression)
  dev.off()
}

png("Mp3g08490.png", width = 2000, height = 1500, res=300)
plot.ld.sd(gene.id = "Mp3g08490",gene.name = " ", ld.gene.expression,sd.gene.expression)
dev.off()

plot.ld.sd(gene.id = "Mp6g11470",gene.name = " ", ld.gene.expression,sd.gene.expression)

# Before pretreatment
genes.expression <- read.table(file ="marchantia_gene_expression.tsv",sep="\t",header=T,as.is=T)
dim(genes.expression)
head(genes.expression)
genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)


# Gene expression matrix for samples under light-dark cycles
ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ld.zt.i <- sapply(X = ld.zt,FUN = function(x){ paste(x,1:3,sep="_")})
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt.i <- sapply(X = sd.zt,FUN = function(x){ paste(x,1:3,sep="_")})

ld.sd.gene.expression <- genes.expression[,c(ld.zt.i,sd.zt.i)]
head(ld.sd.gene.expression)
dim(ld.sd.gene.expression)

ld.gene.expression <- genes.expression[,c(ld.zt.i)]
sd.gene.expression <- genes.expression[,c(sd.zt.i)]

plot.ld.sd(gene.id = "Mp3g09460",gene.name = "MpTOC1", ld.gene.expression,sd.gene.expression)

#Examples taken from a paper
paper.ids <- c("Mp3g09460","Mp2g16560","Mp1g10870","Mp1g16510","Mp1g13750","Mp3g14360","Mp2g14310","Mp4g21220","Mp3g14360")
paper.names <- c("MpTOC1","MpPRR","MpELF3","MpELF4","MpGI","MpZTL","MpRVE","MpLUX","MpFKF")

for (i in 1:length(paper.ids))
{
  png(paste0(paper.names[i],".png"), width = 2000, height = 1500, res=300)
  plot.ld.sd(gene.id = paper.ids[i],gene.name = paper.names[i], ld.gene.expression,sd.gene.expression)
  dev.off()
}

plot.ld.sd(gene.id = "Mp3g08490",gene.name = " ", ld.gene.expression,sd.gene.expression)
plot.ld.sd(gene.id = "Mp6g11470",gene.name = " ", ld.gene.expression,sd.gene.expression)

# Normalyzed
genes.expression <- read.table(file ="separated_groups_norm_expression.tsv",sep="\t",header=T,as.is=T)
rownames(genes.expression) <- genes.expression$genes
genes.expression <- genes.expression[,-1]
dim(genes.expression)
head(genes.expression)

genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)

# Add genes not represented in the normalization
gene.expression.comp <- read.table(file ="marchantia_gene_expression.tsv",sep="\t",header=T,as.is=T)
head(gene.expression.comp)

gene.expression.comp.ld.sd <- gene.expression.comp[,c(ld.zt.i, sd.zt.i)]

added.expression <- gene.expression.comp.ld.sd[which(!(rownames(gene.expression.comp.ld.sd) %in% rownames(genes.expression))),]

genes.expression <- rbind(genes.expression, added.expression)
dim(genes.expression)
genes.expression <- as.matrix(genes.expression)

# Gene expression matrix for samples under light-dark cycles
ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ld.zt.i <- sapply(X = ld.zt,FUN = function(x){ paste(x,1:3,sep="_")})
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt.i <- sapply(X = sd.zt,FUN = function(x){ paste(x,1:3,sep="_")})

ld.sd.gene.expression <- genes.expression[,c(ld.zt.i,sd.zt.i)]
head(ld.sd.gene.expression)
dim(ld.sd.gene.expression)

ld.gene.expression <- genes.expression[,c(ld.zt.i)]
sd.gene.expression <- genes.expression[,c(sd.zt.i)]

plot.ld.sd(gene.id = "Mp3g09460",gene.name = "MpTOC1", ld.gene.expression,sd.gene.expression)

#Examples taken from a paper
paper.ids <- c("Mp3g09460","Mp2g16560","Mp1g10870","Mp1g16510","Mp1g13750","Mp3g14360","Mp2g14310","Mp4g21220","Mp3g14360")
paper.names <- c("MpTOC1","MpPRR","MpELF3","MpELF4","MpGI","MpZTL","MpRVE","MpLUX","MpFKF")

for (i in 1:length(paper.ids))
{
  png(paste0(paper.names[i],".png"), width = 2000, height = 1500, res=300)
  plot.ld.sd(gene.id = paper.ids[i],gene.name = paper.names[i], ld.gene.expression,sd.gene.expression)
  dev.off()
}

plot.ld.sd(gene.id = "Mp3g08490",gene.name = " ", ld.gene.expression,sd.gene.expression)
plot.ld.sd(gene.id = "Mp6g11470",gene.name = " ", ld.gene.expression,sd.gene.expression)


# ADJUST COMPLETE DATASET
library(sva)

genes.expression <- read.table(file ="marchantia_gene_expression.tsv",sep="\t",header=T,as.is=T)
dim(genes.expression)
head(genes.expression)
genes.expression<-as.matrix(genes.expression)
head(genes.expression[,1:36])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)


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


sva.gene.expression.ld <- genes.expression[,c(ld.zt.i, dd.ld.zt.i, ll.ld.zt.i)]

pheno.ld <- data.frame(sample=1:ncol(sva.gene.expression.ld),
                       photoperiod=c(rep("LD",18), rep("DD_LD",12), rep("LL_LD",12)),
                       zt=c(rep(c("ZT00","ZT04","ZT08","ZT12","ZT16","ZT20"), each=3),
                            rep(rep(c("ZT00","ZT04","ZT08","ZT12","ZT16","ZT20"), each=2),2)),
                       batch=c(1,1,2,1,1,1,rep(c(1,1,2),4), rep(2,12), rep(c(1,2),6))
)

rownames(pheno.ld) <- colnames(sva.gene.expression.ld)
head(sva.gene.expression.ld)
batch.ld <- pheno.ld$batch
edata.ld <- sva.gene.expression.ld
modcombat.ld = model.matrix(~1, data=pheno.ld)
combat_edata.ld = ComBat(dat=edata.ld, batch=batch.ld, mod=modcombat.ld, par.prior=TRUE, prior.plots=FALSE, ref.batch = 1)
write.table(combat_edata.ld,file = "combat_gene_expression_ld_comp.tsv",sep="\t")

modcombat.ld2 = model.matrix(~as.factor(photoperiod), data=pheno.ld)
combat_edata.ld = ComBat(dat=edata.ld, batch=batch.ld, mod=modcombat.ld2, par.prior=F, prior.plots=T, ref.batch = 1)
write.table(combat_edata.ld,file = "combat_gene_expression_ld_comp_covar_np.tsv",sep="\t")


sva.gene.expression.sd <- genes.expression[,c(sd.zt.i, dd.sd.zt.i, ll.sd.zt.i)]


pheno.sd <- data.frame(sample=1:ncol(sva.gene.expression.sd),
                       photoperiod=c(rep("SD",18), rep("DD_SD",12), rep("LL_SD",12)),
                       zt=c(rep(c("ZT00","ZT04","ZT08","ZT12","ZT16","ZT20"), each=3),
                            rep(rep(c("ZT00","ZT04","ZT08","ZT12","ZT16","ZT20"), each=2),2)),
                       batch=c(rep(c(1,1,2),6), rep(2,12), rep(c(1,2),6))
)



rownames(pheno.sd) <- colnames(sva.gene.expression.sd)
head(sva.gene.expression.sd)
batch.sd <- pheno.sd$batch
edata.sd <- sva.gene.expression.sd
modcombat.sd = model.matrix(~1, data=pheno.sd)
combat_edata.sd = ComBat(dat=edata.sd, batch=batch.sd, mod=modcombat.sd, par.prior=TRUE, prior.plots=F, ref.batch = 1)
write.table(combat_edata.sd,file = "combat_gene_expression_sd_comp.tsv",sep="\t")

modcombat.sd2 = model.matrix(~as.factor(photoperiod), data=pheno.sd)
combat_edata.sd = ComBat(dat=edata.sd, batch=batch.sd, mod=modcombat.sd2, par.prior=F, prior.plots=F, ref.batch = 1)
write.table(combat_edata.sd,file = "combat_gene_expression_sd_comp_covar_np.tsv",sep="\t")

# RAIN
library(rain)

# Read data
genes.expression.ld <- read.table(file ="combat_gene_expression_ld_comp_covar_np.tsv",sep="\t",header=T,as.is=T)
genes.expression.sd <- read.table(file ="combat_gene_expression_sd_comp_covar_np.tsv",sep="\t",header=T,as.is=T)

sum(rownames(genes.expression.ld) == rownames(genes.expression.sd))

genes.expression <- cbind(genes.expression.ld, genes.expression.sd)
head(genes.expression)
genes.expression<-as.matrix(genes.expression)
head(genes.expression[,c(1:18, 43:60)])
mp.genes <- rownames(genes.expression)
number.genes <- length(mp.genes)

#LD
ld.gene.expression <- genes.expression.ld[,c(ld.zt.i)]

results.ld <- rain(t(ld.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.ld$pVal < 0.05)/number.genes
rhythmic.genes.ld <- rownames(subset(results.ld, pVal < 0.05))
length(rhythmic.genes.ld)

results.ld.12 <- rain(t(ld.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.ld.12$pVal < 0.05)/number.genes
rhythmic.genes.ld.12 <- rownames(subset(results.ld.12, pVal < 0.05))
length(rhythmic.genes.ld.12)

complete.ld.rhythmic.genes <- unique(c(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(complete.ld.rhythmic.genes)/number.genes

non.rhythmic.ld.genes <- setdiff(mp.genes,complete.ld.rhythmic.genes)
length(non.rhythmic.ld.genes)
length(non.rhythmic.ld.genes)/number.genes
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))
length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12))/number.genes

#SD
sd.gene.expression <- genes.expression.sd[,sd.zt.i]
results.sd <- rain(t(sd.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.sd$pVal < 0.05)/number.genes
rhythmic.genes.sd <- rownames(subset(results.sd, pVal < 0.05))
length(rhythmic.genes.sd)

#results.sd[order(results.sd$pVal),]

results.sd.12 <- rain(t(sd.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.sd.12$pVal < 0.05)/number.genes
rhythmic.genes.sd.12 <- rownames(subset(results.sd.12, pVal < 0.05))
length(rhythmic.genes.sd.12)

complete.sd.rhythmic.genes <- unique(c(rhythmic.genes.sd,rhythmic.genes.sd.12))
length(complete.sd.rhythmic.genes)
length(complete.sd.rhythmic.genes)/number.genes

non.rhythmic.sd.genes <- setdiff(mp.genes,complete.sd.rhythmic.genes)
length(non.rhythmic.sd.genes)
length(non.rhythmic.sd.genes)/number.genes

#Rhythmic genes in SD conditions but not in LD
sd.specific.rhythmic.genes <- setdiff(complete.sd.rhythmic.genes,complete.ld.rhythmic.genes)
length(sd.specific.rhythmic.genes)

#Rhythmic genes in LD conditions but not in SD
ld.specific.rhytmic.genes <-setdiff(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(ld.specific.rhytmic.genes)

#Rhythmic genes in both LD and SD
complete.ld.sd.rhythmic.genes<-intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(complete.ld.sd.rhythmic.genes)


#Save lists of rhythmic genes
ld.rhytmic.table<-matrix(data = complete.ld.rhythmic.genes,ncol = 1)
write.table(ld.rhytmic.table,file="ld.rhytmic.table.combat.sep.comp.covar.np.tsv")
sd.rhytmic.table<-matrix(data = complete.sd.rhythmic.genes,ncol = 1)
write.table(sd.rhytmic.table,file="sd.rhytmic.table.combat.sep.comp.covar.np.tsv")

#### Plots of rhythmic genes

#Venn diagram of LD and SD rhythmic genes
library(VennDiagram)
ld.rhythmic.combat<- read.table("ld.rhytmic.table.combat.sep.comp.covar.np.tsv", header = T)$V1
sd.rhythmic.combat <- read.table("sd.rhytmic.table.combat.sep.comp.covar.np.tsv", header = T)$V1
ld.rhythmic.2 <- read.table("ld.rhytmic.table.2.tsv", header = T)$V1
sd.rhythmic.2 <- read.table("sd.rhytmic.table.2.tsv", header = T)$V1

png("venn_2_vs_combat_sep_ld.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(ld.rhythmic.combat), area2 = length(ld.rhythmic.2),
                   cross.area = length(intersect(ld.rhythmic.combat,ld.rhythmic.2)),
                   lwd = 3,category = c("LD_combat","LD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

png("venn_2_vs_combat_sep_sd.png", width = 1500, height = 1500, res=300)
grid.newpage()
draw.pairwise.venn(area1 = length(sd.rhythmic.combat), area2 = length(sd.rhythmic.2),
                   cross.area = length(intersect(sd.rhythmic.combat,sd.rhythmic.2)),
                   lwd = 3,category = c("SD_combat","SD_2"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# Pruebas al completo y al separado

plot.ld.sd(gene.id = "Mp3g09460",gene.name = "MpTOC1", as.matrix(ld.gene.expression),as.matrix(sd.gene.expression))

#Examples taken from a paper
paper.ids <- c("Mp3g09460","Mp2g16560","Mp1g10870","Mp1g16510","Mp1g13750","Mp3g14360","Mp2g14310","Mp4g21220","Mp3g14360")
paper.names <- c("MpTOC1","MpPRR","MpELF3","MpELF4","MpGI","MpZTL","MpRVE","MpLUX","MpFKF")

for (i in 1:length(paper.ids))
{
  png(paste0(paper.names[i],".png"), width = 2000, height = 1500, res=300)
  plot.ld.sd(gene.id = paper.ids[i],gene.name = paper.names[i], ld.gene.expression,sd.gene.expression)
  dev.off()
}

plot.ld.sd(gene.id = "Mp3g08490",gene.name = " ", as.matrix(ld.gene.expression),as.matrix(sd.gene.expression))
plot.ld.sd(gene.id = "Mp6g11470",gene.name = " ", as.matrix(ld.gene.expression),as.matrix(sd.gene.expression))

# Free-running
genes.expression.ld <- read.table(file ="combat_gene_expression_ld_comp_covar.tsv",sep="\t",header=T,as.is=T)
genes.expression.sd <- read.table(file ="combat_gene_expression_sd_comp_covar.tsv",sep="\t",header=T,as.is=T)

genes.expression <- cbind(genes.expression.ld, genes.expression.sd)

plot.ld.sd(gene.id = "Mp7g03540",gene.name = " ", as.matrix(genes.expression.ld),as.matrix(genes.expression.sd))
plot.sd.ll(gene.id = "Mp7g03540", gene.name = " ", genes.expression = as.matrix(genes.expression))
plot.sd.dd(gene.id = "Mp7g03540", gene.name = " ", genes.expression = as.matrix(genes.expression))

"Mp7g03540"  "Mp5g03810" "Mp3g09450"

