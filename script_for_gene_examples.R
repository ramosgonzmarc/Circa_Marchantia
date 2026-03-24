
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


# Load gene information
complete.ld.rhythmic.genes <- read.table("tables/ld.rhythmic.table.tsv", header = T)$V1
complete.sd.rhythmic.genes <- read.table("tables/sd.rhythmic.table.tsv", header = T)$V1

complete.sd.dd.rhythmic.genes <- read.table("tables/sd.dd.rhythmic.table.tsv", sep="\t")$x
complete.ld.dd.rhythmic.genes <- read.table("tables/ld.dd.rhythmic.table.tsv", sep="\t")$x
complete.sd.ll.rhythmic.genes <- read.table("tables/sd.ll.rhythmic.table.tsv", sep="\t")$x
complete.ld.ll.rhythmic.genes <- read.table("tables/ld.ll.rhythmic.table.tsv", sep="\t")$x

circacompare.sd.dd <- read.table("tables/circacompare_sd_dd.tsv", sep="\t")
circacompare.sd.ll <- read.table("tables/circacompare_sd_ll.tsv", sep="\t")
circacompare.ld.dd <- read.table("tables/circacompare_ld_dd.tsv", sep="\t")
circacompare.ld.ll <- read.table("tables/circacompare_ld_ll.tsv", sep="\t")

circa_sd_ll <- read.table("tables/estimates_sd_ll.tsv", sep="\t")
circa_sd_dd <- read.table("tables/estimates_sd_dd.tsv", sep="\t")
circa_ld_ll <- read.table("tables/estimates_ld_ll.tsv", sep="\t")
circa_ld_dd <- read.table("tables/estimates_ld_dd.tsv", sep="\t")

ld_ll_error <- read.table("tables/estimates_ld_ll_with_errors.tsv")
sd_ll_error <- read.table("tables/estimates_sd_ll_with_errors.tsv")
ld_dd_error <- read.table("tables/estimates_ld_dd_with_errors.tsv")
sd_dd_error <- read.table("tables/estimates_sd_dd_with_errors.tsv")


circacompare.ld.sd <- read.table("tables/circacompare_intersection.tsv", sep="\t")
ld_sd_error <- read.table("tables/circacompare_ld_sd_complete_with_errors.tsv")

# Example of genes retaining rhythmicity in LL and not DD
source("functions_circa_marchantia.R")
circacompare.ld.ll.high <- subset(circacompare.ld.ll, rhythmic_p < 0.005)
circacompare.sd.ll.high <- subset(circacompare.sd.ll, rhythmic_p < 0.005)
ordered_ld_ll <- rownames(circacompare.ld.ll.high[order(circacompare.ld.ll.high$rhythmic_p),])
ordered_sd_ll <- rownames(circacompare.sd.ll.high[order(circacompare.sd.ll.high$rhythmic_p),])
filtered_ld_sd_ll <- intersect(ordered_sd_ll, ordered_ld_ll)

examples_ll_no_dd <- setdiff(intersect(complete.ld.ll.rhythmic.genes, complete.sd.ll.rhythmic.genes),
        union(complete.sd.dd.rhythmic.genes, complete.ld.dd.rhythmic.genes))

possible_examples <- filtered_ld_sd_ll[which(filtered_ld_sd_ll %in% examples_ll_no_dd )]

phase_delayed_ll <- intersect(rownames(circa_sd_ll[circa_sd_ll$Phase.difference.estimate < 0,]), 
          rownames(circa_ld_ll[circa_ld_ll$Phase.difference.estimate > 0,]))

phase_delayed_ll_and_no_dd <- intersect(phase_delayed_ll, examples_ll_no_dd)
ordered_ld_ll[which(ordered_ld_ll %in% phase_delayed_ll_and_no_dd)]

plot.ld.ll("Mp8g12930", "", genes.expression)
plot.ld.dd("Mp8g12930", "", genes.expression)
plot.sd.ll("Mp8g12930", "", genes.expression)
plot.sd.dd("Mp8g12930", "", genes.expression)


par(mfrow=c(2,2))
for (x in possible_examples)
{
  plot.ld.ll(x, "", genes.expression)
  plot.ld.dd(x, "", genes.expression)
  plot.sd.ll(x, "", genes.expression)
  plot.sd.dd(x, "", genes.expression)
}

# Mp2g06190, Mp3g09890 patrón interesante de dobles picos alternos, 
# Mp8g12930 el mejor hasta ahora


# Examples for bona fide genes
bona_fide <- intersect(intersect(intersect(complete.sd.dd.rhythmic.genes, complete.ld.dd.rhythmic.genes),
                                 complete.sd.ll.rhythmic.genes), complete.ld.ll.rhythmic.genes)




par(mfrow=c(2,2))
for (x in bona_fide)
{
  plot.ld.ll(x, "", genes.expression)
  plot.ld.dd(x, "", genes.expression)
  plot.sd.ll(x, "", genes.expression)
  plot.sd.dd(x, "", genes.expression)
}



# Mp7g17400
plot.ld.ll("Mp7g17400", "", genes.expression)
plot.ld.dd("Mp7g17400", "", genes.expression)
plot.sd.ll("Mp7g17400", "", genes.expression)
plot.sd.dd("Mp7g17400", "", genes.expression)




# Examples for phase delays under LD LL
circacompare.ld.ll.high <- subset(ld_ll_error, Presence.of.rhythmicity..p.value..for.ll < 0.005)
circacompare.ld.high <- subset(ld_ll_error, Presence.of.rhythmicity..p.value..for.ld < 0.00001)
ordered_ld_ll <- rownames(circacompare.ld.ll.high[order(circacompare.ld.ll.high$Presence.of.rhythmicity..p.value..for.ll),])
ordered_ld <- rownames(circacompare.ld.high[order(circacompare.ld.high$Presence.of.rhythmicity..p.value..for.ld),])
filtered_ld_ll <- intersect(ordered_ld, ordered_ld_ll)

phase_delayed_ll <- intersect(filtered_ld_ll, 
                              rownames(circa_ld_ll[circa_ld_ll$Phase.difference.estimate > 2,]))

library(ggplot2)
plot.ld.ll("Mp6g07690", "", genes.expression)
plot.ld.ll.wave.error("Mp1g09490", "", genes.expression, period=1, 
                      wave_parameters.df = ld_ll_error, data.points = T, 
                      color.1 = "royalblue3", color.2="skyblue")

ggsave("images/examples/ld_ll_delays.png")

# Examples for phase anticipation under LD DD
circacompare.ld.dd.high <- subset(ld_dd_error, Presence.of.rhythmicity..p.value..for.dd < 0.05)
circacompare.ld.high <- subset(ld_dd_error, Presence.of.rhythmicity..p.value..for.ld < 0.001)
ordered_ld_dd <- rownames(circacompare.ld.dd.high[order(circacompare.ld.dd.high$Presence.of.rhythmicity..p.value..for.dd),])
ordered_ld <- rownames(circacompare.ld.high[order(circacompare.ld.high$Presence.of.rhythmicity..p.value..for.ld),])
filtered_ld_dd <- intersect(ordered_ld, ordered_ld_dd)

phase_anticipated_dd <- intersect(filtered_ld_dd, 
                              rownames(circa_ld_dd[circa_ld_dd$Phase.difference.estimate > 3,]))
intersect(phase_anticipated_dd, 
          rownames(circa_ld_dd[circa_ld_dd$Mesor.difference.estimate > -10 & circa_ld_dd$Mesor.difference.estimate < 10 &
                               circa_ld_dd$Amplitude.difference.estimate > 10,]))


plot.ld.dd.wave.error("Mp5g21020", "", genes.expression, period=1, 
                      wave_parameters.df = ld_dd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="darkblue")

ggsave("images/examples/ld_dd_anticipation.png")

# Examples for phase delays under SD LL
circacompare.sd.ll.high <- subset(sd_ll_error, Presence.of.rhythmicity..p.value..for.ll < 0.005)
circacompare.sd.high <- subset(sd_ll_error, Presence.of.rhythmicity..p.value..for.sd < 0.00001)
ordered_sd_ll <- rownames(circacompare.sd.ll.high[order(circacompare.sd.ll.high$Presence.of.rhythmicity..p.value..for.ll),])
ordered_sd <- rownames(circacompare.sd.high[order(circacompare.sd.high$Presence.of.rhythmicity..p.value..for.sd),])
filtered_sd_ll <- intersect(ordered_sd, ordered_sd_ll)

phase_delayed_ll <- intersect(filtered_sd_ll, 
                              rownames(circa_sd_ll[circa_sd_ll$Phase.difference.estimate < -3,]))

intersect(phase_delayed_ll, 
          rownames(circa_sd_ll[circa_sd_ll$Mesor.difference.estimate > -10 & circa_sd_ll$Mesor.difference.estimate < 10 &
                                 circa_sd_ll$Amplitude.difference.estimate > 5,]))


library(ggplot2)

plot.sd.ll.wave.error("Mp5g04500", "", genes.expression, period=1, 
                      wave_parameters.df = sd_ll_error, data.points = T, 
                      color.1 = "tomato3", color.2="salmon")
ggsave("images/examples/sd_ll_delays.png")

# Examples for phase delays under SD dd
circacompare.sd.dd.high <- subset(sd_dd_error, Presence.of.rhythmicity..p.value..for.dd < 0.01)
circacompare.sd.high <- subset(sd_dd_error, Presence.of.rhythmicity..p.value..for.sd < 0.005)
ordered_sd_dd <- rownames(circacompare.sd.dd.high[order(circacompare.sd.dd.high$Presence.of.rhythmicity..p.value..for.dd),])
ordered_sd <- rownames(circacompare.sd.high[order(circacompare.sd.high$Presence.of.rhythmicity..p.value..for.sd),])
filtered_sd_dd <- intersect(ordered_sd, ordered_sd_dd)

phase_delayed_dd <- intersect(filtered_sd_dd, 
                              rownames(circa_sd_dd[circa_sd_dd$Phase.difference.estimate < -2,]))

library(ggplot2)

plot.sd.dd.wave.error("Mp1g19330", "", genes.expression, period=1, 
                      wave_parameters.df = sd_dd_error, data.points = T, 
                      color.1 = "tomato3", color.2="red4")
ggsave("images/examples/sd_dd_delays.png")

# Example for LD/SD cycle: anticipated phase and increased mesor in SD
circacompare.ld.high <- subset(ld_sd_error, Presence.of.rhythmicity..p.value..for.ld < 0.00005)
circacompare.sd.high <- subset(ld_sd_error, Presence.of.rhythmicity..p.value..for.sd < 0.00005)
ordered_ld <- rownames(circacompare.ld.high[order(circacompare.ld.high$Presence.of.rhythmicity..p.value..for.ld),])
ordered_sd <- rownames(circacompare.sd.high[order(circacompare.sd.high$Presence.of.rhythmicity..p.value..for.sd),])
filtered_ld_sd <- intersect(ordered_sd, ordered_ld)

phase_delayed_ld <- intersect(intersect(intersect(filtered_ld_sd, 
                              rownames(ld_sd_error[ld_sd_error$Phase.difference.estimate < -4,])),
                              rownames(ld_sd_error[ld_sd_error$Mesor.difference.estimate > 10,])),
                              rownames(ld_sd_error[ld_sd_error$Amplitude.difference.estimate < 5 & ld_sd_error$Amplitude.difference.estimate > -5,]))

library(ggplot2)
plot.ld.sd.wave.error("Mp1g20380", "", genes.expression, period=1, 
                      wave_parameters.df = ld_sd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="tomato3")

ggsave("images/examples/ld_sd_anticipation.png")

# Bona fide morning LD peak absent in SD
bona_fide <- intersect(intersect(intersect(complete.sd.dd.rhythmic.genes, complete.ld.dd.rhythmic.genes),
                                 complete.sd.ll.rhythmic.genes), complete.ld.ll.rhythmic.genes)

# ld_12 <- read.table("tables/circacompare_ld_12.tsv", sep="\t")
# 
# ld_12_genes <- rownames(subset(ld_12, rhythmic_p < 0.05))
# 
# bona_fide_12_ld <- intersect(bona_fide, ld_12_genes)
# 
# plot.ld.sd("Mp8g17710", "", genes.expression)

bona_fide_morning_ld_evening_sd <- intersect(bona_fide, 
                                             rownames(ld_sd_error[ld_sd_error$ld.peak.time.hours < 12 &
                                                                    ld_sd_error$ld.peak.time.hours > 2 &
                                                                    ld_sd_error$sd.peak.time.hours > 12,]))
x <- "Mp6g07940"
plot.ld.sd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_sd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="tomato3")
ggsave("images/examples/bona_fide_ld_sd.png")

plot.ld.ll.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_ll_error, data.points = T, 
                      color.1 = "royalblue3", color.2="skyblue")
ggsave("images/examples/bona_fide_ld_ll.png")

plot.ld.dd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_dd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="darkblue")
ggsave("images/examples/bona_fide_ld_dd.png")

plot.sd.ll.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = sd_ll_error, data.points = T, 
                      color.1 = "tomato3", color.2="salmon")
ggsave("images/examples/bona_fide_sd_ll.png")

plot.sd.dd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = sd_dd_error, data.points = T, 
                      color.1 = "tomato3", color.2="red4")
ggsave("images/examples/bona_fide_sd_dd.png")


# Examples for histograms
ld.light.response <- read.table("tables/LD_light_activated_dark_repressed.tsv", sep = "\t", header = T)$x
ld.dark.response <- read.table("tables/LD_light_repressed_dark_activated.tsv", sep = "\t", header = T)$x
sd.light.response <- read.table("tables/SD_light_activated_dark_repressed.tsv", sep = "\t", header = T)$x
sd.dark.response <- read.table("tables/SD_light_repressed_dark_activated.tsv", sep = "\t", header = T)$x

# Example for first peak in SD, peak in LD, dark response
intersect(ld.dark.response, rownames(ld_sd_error[ld_sd_error$ld.peak.time.hours > 14 &
                                                 ld_sd_error$sd.peak.time.hours < 10 ,]))

plot.ld.sd.wave.error("Mp1g02250", "", genes.expression, period=1, 
                      wave_parameters.df = ld_sd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="tomato3")
ggsave("images/examples/peak_1_sd_peak_ld_dark_resp.png")


# Example for second peak in SD, peak in LD, dark response
intersect(ld.dark.response, rownames(ld_sd_error[ld_sd_error$ld.peak.time.hours > 14 &
                                                   ld_sd_error$sd.peak.time.hours > 14 ,]))


plot.ld.sd.wave.error("Mp1g01630", "", genes.expression, period=1, 
                      wave_parameters.df = ld_sd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="tomato3")
ggsave("images/examples/peak_2_sd_peak_ld_dark_resp.png")

# Example for first peak in SD, peak in LD, light response in SD
intersect(sd.light.response, rownames(ld_sd_error[ld_sd_error$ld.peak.time.hours > 15 &
                                                   ld_sd_error$sd.peak.time.hours < 9 ,]))

plot.ld.sd.wave.error("Mp2g11080", "", genes.expression, period=1, 
                      wave_parameters.df = ld_sd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="tomato3")
ggsave("images/examples/peak_1_sd_peak_ld_light_resp_sd.png")

# Example for LL circadian in SD and LD with morning phase in LD and night phase in SD

ll_circadian_both <- intersect(complete.sd.ll.rhythmic.genes, complete.ld.ll.rhythmic.genes)



ll_circadian_morning_ld_evening_sd <- intersect(complete.sd.ll.rhythmic.genes, 
                                             rownames(ld_sd_error[ld_sd_error$ld.peak.time.hours < 4 &
                                                                    ld_sd_error$ld.peak.time.hours > 1 &
                                                                    ld_sd_error$sd.peak.time.hours > 12 &
                                                                    ld_sd_error$sd.peak.time.hours < 18,]))

plot.ld.sd.wave.error("Mp2g10760", "", genes.expression, period=1, 
                      wave_parameters.df = ld_sd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="tomato3")
ggsave("images/examples/night_sd_peak_morning_ld_peak_ll_circadian_ld_sd.png")

plot.ld.ll.wave.error("Mp2g10760", "", genes.expression, period=1, 
                      wave_parameters.df = ld_ll_error, data.points = T, 
                      color.1 = "royalblue3", color.2="skyblue")
ggsave("images/examples/night_sd_peak_morning_ld_peak_ll_circadian_ld_ll.png")

plot.sd.ll.wave.error("Mp2g10760", "", genes.expression, period=1, 
                      wave_parameters.df = sd_ll_error, data.points = T, 
                      color.1 = "tomato3", color.2="salmon")
ggsave("images/examples/night_sd_peak_morning_ld_peak_ll_circadian_sd_ll.png")


# HY5, TRIHELIX11, TRIHELIX31, DET1
x <- "Mp1g16800"
plot.ld.sd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_sd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="tomato3")
ggsave("images/examples/hy5_ld_sd.png")

plot.ld.ll.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_ll_error, data.points = T, 
                      color.1 = "royalblue3", color.2="skyblue")
ggsave("images/examples/hy5_ld_ll.png")

plot.ld.dd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_dd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="darkblue")
ggsave("images/examples/hy5_ld_dd.png")

plot.sd.ll.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = sd_ll_error, data.points = T, 
                      color.1 = "tomato3", color.2="salmon")
ggsave("images/examples/hy5_sd_ll.png")

plot.sd.dd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = sd_dd_error, data.points = T, 
                      color.1 = "tomato3", color.2="red4")
ggsave("images/examples/hy5_sd_dd.png")

x <- "Mp2g08360" 
plot.ld.sd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_sd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="tomato3")
ggsave("images/examples/trihelix11_ld_sd.png")

plot.ld.ll.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_ll_error, data.points = T, 
                      color.1 = "royalblue3", color.2="skyblue")
ggsave("images/examples/trihelix11_ld_ll.png")

plot.ld.dd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_dd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="darkblue")
ggsave("images/examples/trihelix11_ld_dd.png")

plot.sd.ll.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = sd_ll_error, data.points = T, 
                      color.1 = "tomato3", color.2="salmon")
ggsave("images/examples/trihelix11_sd_ll.png")

plot.sd.dd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = sd_dd_error, data.points = T, 
                      color.1 = "tomato3", color.2="red4")
ggsave("images/examples/trihelix11_sd_dd.png")

x <- "Mp4g09730" 
plot.ld.sd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_sd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="tomato3")
ggsave("images/examples/trihelix31_ld_sd.png")

plot.ld.ll.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_ll_error, data.points = T, 
                      color.1 = "royalblue3", color.2="skyblue")
ggsave("images/examples/trihelix31_ld_ll.png")

plot.ld.dd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_dd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="darkblue")
ggsave("images/examples/trihelix31_ld_dd.png")

plot.sd.ll.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = sd_ll_error, data.points = T, 
                      color.1 = "tomato3", color.2="salmon")
ggsave("images/examples/trihelix31_sd_ll.png")

plot.sd.dd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = sd_dd_error, data.points = T, 
                      color.1 = "tomato3", color.2="red4")
ggsave("images/examples/trihelix31_sd_dd.png")


x <- "Mp8g00970"
plot.ld.sd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_sd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="tomato3")
ggsave("supplementary_data_for_paper/s5/det1_ld_sd.png")

plot.ld.ll.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_ll_error, data.points = T, 
                      color.1 = "royalblue3", color.2="skyblue")
ggsave("supplementary_data_for_paper/s5/det1_ld_ll.png")

plot.ld.dd.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = ld_dd_error, data.points = T, 
                      color.1 = "royalblue3", color.2="darkblue")
ggsave("supplementary_data_for_paper/s5/det1_ld_dd.png")

plot.sd.ll.wave.error(x, "", genes.expression, period=1, 
                      wave_parameters.df = sd_ll_error, data.points = T, 
                      color.1 = "tomato3", color.2="salmon")
ggsave("supplementary_data_for_paper/s5/det1_sd_ll.png")




# Two peaks in SD, only one in LD
sd_12_errors <- read.table("tables/circacompare_sd_12_with_errors.tsv", sep = "\t")
ld_24_errors <- read.table("tables/circacompare_ld_24_with_errors.tsv", sep = "\t")

sd_12_df <- read.table("tables/circacompare_sd_12.tsv", sep = "\t")
sd_12 <- rownames(subset(sd_12_df, rhythmic_p < 0.0005))

ld_12_df <- read.table("tables/circacompare_ld_12.tsv", sep = "\t")
ld_no_12 <- rownames(subset(ld_12_df, rhythmic_p > 0.4))

intersect(sd_12, ld_no_12)

plot.ld.sd("Mp1g03010", "", genes.expression)

library(ggplot2)
plot.ld.sd.wave.error.12("Mp1g03010", "", genes.expression, period1=1, period2=2, 
                      wave_parameters.df1 = ld_24_errors, wave_parameters.df2 = sd_12_errors, 
                      data.points = T, 
                      color.1 = "royalblue3", color.2="tomato3")
ggsave("images/examples/sd_two_peaks.png")





