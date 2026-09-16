.libPaths("/home/marcos/R/x86_64-pc-linux-gnu-library/4.1")

library(GenomicFeatures)
txdb_mp <- makeTxDbFromGFF('annotation.gff')
trans_mp <- transcripts(txdb_mp)
as.data.frame(trans_mp)

library(ChIPseeker)
#promoters <- getPromoters(TxDb = txdb_mp, upstream = 2500, downstream = 0)
transcriptCoordsByGene.GRangesList <- transcriptsBy(txdb_mp, by = "gene") 

library(Rsamtools)
library(Biostrings)

# Extract lenghts of 5' UTR
#complete_len <- transcriptLengths(txdb_mp, with.cds_len=FALSE,
#                  with.utr5_len=T, with.utr3_len=FALSE)
#mean(complete_len$utr5_len)

# With UTRs
# promoter.seqs <- getPromoterSeq(transcriptCoordsByGene.GRangesList,
#                                 FaFile("genome.fa", index="genome.fa.fai"),
#                                 upstream=1800, downstream=535)

# Without UTRs
promoter.seqs <- getPromoterSeq(transcriptCoordsByGene.GRangesList,
                                FaFile("genome.fa", index="genome.fa.fai"),
                                upstream=1800, downstream=0)


promoter.unlist <- unlist(promoter.seqs)
promoter.table <- data.frame(width=width(promoter.unlist), seq=as.character(promoter.unlist), names=names(promoter.unlist))
promoter.table$names <- sapply(strsplit(promoter.table$names, ".chr"), function(x) x[[1]])

promoter.table.clean <- promoter.table[!duplicated(promoter.table$names),]

library(seqinr)
write.fasta(as.list(promoter.table.clean$seq),promoter.table.clean$names, "marchantia_promoters_1.8kb_without_utrs.fa", as.string = T)

# Create fastas
library(seqinr)
mp_promoters <- read.fasta("marchantia_promoters_1.8kb_without_utrs.fa", seqtype = "DNA", forceDNAtolower = F)

ld_complete <- read.table("../tables/ld.rhythmic.table.tsv", header = T)$V1
promoters_ld <- mp_promoters[ld_complete]
write.fasta(getSequence(promoters_ld), getName(promoters_ld), "promoters_1.8kb_no_utr/promoters_ld_complete.fa")

sd_complete <- read.table("../tables/sd.rhythmic.table.tsv", header = T)$V1
promoters_sd <- mp_promoters[sd_complete]
write.fasta(getSequence(promoters_sd), getName(promoters_sd), "promoters_1.8kb_no_utr/promoters_sd_complete.fa")

cycle_union <- unique(c(ld_complete, sd_complete))
promoters_cycle <- mp_promoters[cycle_union]
write.fasta(getSequence(promoters_cycle), getName(promoters_cycle), "promoters_1.8kb_no_utr/promoters_cycle_union.fa")

# Free running
# LL
ld_ll <- read.table("../tables/ld.ll.rhythmic.table.tsv", sep = "\t")$x
sd_ll <- read.table("../tables/sd.ll.rhythmic.table.tsv", sep = "\t")$x
ll_inter <- intersect(ld_ll, sd_ll)
ll_union <- unique(c(ld_ll, sd_ll))

ld_ll_fasta <- mp_promoters[ld_ll]
write.fasta(getSequence(ld_ll_fasta), getName(ld_ll_fasta), "promoters_1.8kb_no_utr/promoters_ld_ll.fa")

ll_inter_fasta <- mp_promoters[ll_inter]
write.fasta(getSequence(ll_inter_fasta), getName(ll_inter_fasta), "promoters_1.8kb_no_utr/promoters_ll_inter.fa")

ll_union_fasta <- mp_promoters[ll_union]
write.fasta(getSequence(ll_union_fasta), getName(ll_union_fasta), "promoters_1.8kb_no_utr/promoters_ll_union.fa")

# DD
ld_dd <- read.table("../tables/ld.dd.rhythmic.table.tsv", sep = "\t")$x
sd_dd <- read.table("../tables/sd.dd.rhythmic.table.tsv", sep = "\t")$x
dd_inter <- intersect(ld_dd, sd_dd)
dd_union <- unique(c(ld_dd, sd_dd))

ld_dd_fasta <- mp_promoters[ld_dd]
write.fasta(getSequence(ld_dd_fasta), getName(ld_dd_fasta), "promoters_1.8kb_no_utr/promoters_ld_dd.fa")

dd_inter_fasta <- mp_promoters[dd_inter]
write.fasta(getSequence(dd_inter_fasta), getName(dd_inter_fasta), "promoters_1.8kb_no_utr/promoters_dd_inter.fa")

dd_union_fasta <- mp_promoters[dd_union]
write.fasta(getSequence(dd_union_fasta), getName(dd_union_fasta), "promoters_1.8kb_no_utr/promoters_dd_union.fa")

free_running_union <- unique(c(dd_union, ll_union))
free_running_union_fasta <- mp_promoters[free_running_union]
write.fasta(getSequence(free_running_union_fasta), getName(free_running_union_fasta),
            "promoters_1.8kb_no_utr/promoters_free_running_union.fa")

# LL 24h
ld_ll_24 <- read.table("../tables/circacompare_ld_ll.tsv", sep = "\t")
head(ld_ll_24)
ld_ll_24 <- rownames(subset(ld_ll_24, rhythmic_p < 0.05))

sd_ll_24 <- read.table("../tables/circacompare_sd_ll.tsv", sep = "\t")
head(sd_ll_24)
sd_ll_24 <- rownames(subset(sd_ll_24, rhythmic_p < 0.05))

ll_24 <- intersect(ld_ll_24, sd_ll_24)
promoters_ll_24 <- mp_promoters[ll_24]
write.fasta(getSequence(promoters_ll_24), getName(promoters_ll_24), "promoters_ll_24.fa")

# DD 24h
ld_dd_24 <- read.table("../tables/circacompare_ld_dd.tsv", sep = "\t")
head(ld_dd_24)
ld_dd_24 <- rownames(subset(ld_dd_24, rhythmic_p < 0.05))

sd_dd_24 <- read.table("../tables/circacompare_sd_dd.tsv", sep = "\t")
head(sd_dd_24)
sd_dd_24 <- rownames(subset(sd_dd_24, rhythmic_p < 0.05))

dd_24 <- intersect(ld_dd_24, sd_dd_24)
promoters_dd_24 <- mp_promoters[dd_24]
write.fasta(getSequence(promoters_dd_24), getName(promoters_dd_24), "promoters_dd_24.fa")

# LD/SD ZTs
# LD, first cluster by time point, only for those with 24h period
ld_complete_circa <- read.table("../tables/circacompare_ld_24.tsv", header = T)
head(ld_complete_circa)
ld_rhythmic <- subset(ld_complete_circa, rhythmic_p < 0.05)
head(ld_rhythmic)

ld_zt0 <- rownames(ld_rhythmic)[ld_rhythmic$peak_time_hours < 2 | ld_rhythmic$peak_time_hours >= 22]
promoters_ld_zt0 <- mp_promoters[ld_zt0]
write.fasta(getSequence(promoters_ld_zt0), getName(promoters_ld_zt0), "promoters_ld_zt0.fa")

ld_zt4 <- rownames(ld_rhythmic)[ld_rhythmic$peak_time_hours < 6 & ld_rhythmic$peak_time_hours >= 2]
promoters_ld_zt4 <- mp_promoters[ld_zt4]
write.fasta(getSequence(promoters_ld_zt4), getName(promoters_ld_zt4), "promoters_ld_zt4.fa")

ld_zt8 <- rownames(ld_rhythmic)[ld_rhythmic$peak_time_hours < 10 & ld_rhythmic$peak_time_hours >= 6]
promoters_ld_zt8 <- mp_promoters[ld_zt8]
write.fasta(getSequence(promoters_ld_zt8), getName(promoters_ld_zt8), "promoters_ld_zt8.fa")

ld_zt12 <- rownames(ld_rhythmic)[ld_rhythmic$peak_time_hours < 14 & ld_rhythmic$peak_time_hours >= 10]
promoters_ld_zt12 <- mp_promoters[ld_zt12]
write.fasta(getSequence(promoters_ld_zt12), getName(promoters_ld_zt12), "promoters_ld_zt12.fa")

ld_zt16 <- rownames(ld_rhythmic)[ld_rhythmic$peak_time_hours < 18 & ld_rhythmic$peak_time_hours >= 14]
promoters_ld_zt16 <- mp_promoters[ld_zt16]
write.fasta(getSequence(promoters_ld_zt16), getName(promoters_ld_zt16), "promoters_ld_zt16.fa")

ld_zt20 <- rownames(ld_rhythmic)[ld_rhythmic$peak_time_hours < 22 & ld_rhythmic$peak_time_hours >= 18]
promoters_ld_zt20 <- mp_promoters[ld_zt20]
write.fasta(getSequence(promoters_ld_zt20), getName(promoters_ld_zt20), "promoters_ld_zt20.fa")

# SD, first cluster by time point, only for those with 24h period
sd_complete_circa <- read.table("../tables/circacompare_sd_24.tsv", header = T)
head(sd_complete_circa)
sd_rhythmic <- subset(sd_complete_circa, rhythmic_p < 0.05)
head(sd_rhythmic)

sd_zt0 <- rownames(sd_rhythmic)[sd_rhythmic$peak_time_hours < 2 | sd_rhythmic$peak_time_hours >= 22]
promoters_sd_zt0 <- mp_promoters[sd_zt0]
write.fasta(getSequence(promoters_sd_zt0), getName(promoters_sd_zt0), "promoters_sd_zt0.fa")

sd_zt4 <- rownames(sd_rhythmic)[sd_rhythmic$peak_time_hours < 6 & sd_rhythmic$peak_time_hours >= 2]
promoters_sd_zt4 <- mp_promoters[sd_zt4]
write.fasta(getSequence(promoters_sd_zt4), getName(promoters_sd_zt4), "promoters_sd_zt4.fa")

sd_zt8 <- rownames(sd_rhythmic)[sd_rhythmic$peak_time_hours < 10 & sd_rhythmic$peak_time_hours >= 6]
promoters_sd_zt8 <- mp_promoters[sd_zt8]
write.fasta(getSequence(promoters_sd_zt8), getName(promoters_sd_zt8), "promoters_sd_zt8.fa")

sd_zt12 <- rownames(sd_rhythmic)[sd_rhythmic$peak_time_hours < 14 & sd_rhythmic$peak_time_hours >= 10]
promoters_sd_zt12 <- mp_promoters[sd_zt12]
write.fasta(getSequence(promoters_sd_zt12), getName(promoters_sd_zt12), "promoters_sd_zt12.fa")

sd_zt16 <- rownames(sd_rhythmic)[sd_rhythmic$peak_time_hours < 18 & sd_rhythmic$peak_time_hours >= 14]
promoters_sd_zt16 <- mp_promoters[sd_zt16]
write.fasta(getSequence(promoters_sd_zt16), getName(promoters_sd_zt16), "promoters_sd_zt16.fa")

sd_zt20 <- rownames(sd_rhythmic)[sd_rhythmic$peak_time_hours < 22 & sd_rhythmic$peak_time_hours >= 18]
promoters_sd_zt20 <- mp_promoters[sd_zt20]
write.fasta(getSequence(promoters_sd_zt20), getName(promoters_sd_zt20), "promoters_sd_zt20.fa")


# Period 12 vs complete
ld_period_12 <- read.table("../tables/circacompare_ld_12.tsv", sep = "\t")
head(ld_period_12)
ld_period_12 <- rownames(subset(ld_period_12, rhythmic_p < 0.05))
promoters_ld_period_12 <- mp_promoters[ld_period_12]
write.fasta(getSequence(promoters_ld_period_12), getName(promoters_ld_period_12), "promoters_ld_period_12.fa")

sd_period_12 <- read.table("../tables/circacompare_sd_12.tsv", sep = "\t")
head(sd_period_12)
sd_period_12 <- rownames(subset(sd_period_12, rhythmic_p < 0.05))
promoters_sd_period_12 <- mp_promoters[sd_period_12]
write.fasta(getSequence(promoters_sd_period_12), getName(promoters_sd_period_12), "promoters_sd_period_12.fa")

################### Analysis of promoters in LD and SD using GeneSCARAB #########################
files_ld <- sapply(strsplit(list.files("results_circa_genescarab/ld_circa_genescarab/"), ".txt"), function(x) x[[1]])
list_ld <- list()

for (x in files_ld)
{
  new_table <- read.table(paste0("results_circa_genescarab/ld_circa_genescarab/", x, ".txt"), header = T, sep="\t")
  new_entry <- unique(subset(new_table, MotifScore > 9)$FASTA.ID)
  list_ld[[x]] <- new_entry
}

names(which(sapply(list_ld, function(x) "Mp8g00970" %in% x)))

files_sd <- sapply(strsplit(list.files("results_circa_genescarab/sd_circa_genescarab/"), ".txt"), function(x) x[[1]])
list_sd <- list()

for (x in files_sd)
{
  new_table <- read.table(paste0("results_circa_genescarab/sd_circa_genescarab/", x, ".txt"), header = T, sep="\t")
  new_entry <- unique(subset(new_table, MotifScore > 9)$FASTA.ID)
  list_sd[[x]] <- new_entry
}

# GeneSCARAB for differential analysis
source("../../../../PSEA_R/GeneSCARAB/PSEA_functions_package.R")

# Load phases
circa_table_genescarab <- read.table("../../circa_marchantia_paper/tables/circacompare_intersection.tsv", header = T, sep = "\t")
total_phases_table_ld <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["ld.peak.time.hours"]]))
total_phases_table_sd <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["sd.peak.time.hours"]]))

# Generate phase lists for LD and SD
phases.list.ld <- gene.list.to.phases(list_ld, total_phases_table_ld)
phases.list.sd <- gene.list.to.phases(list_sd, total_phases_table_sd)

# Remove entries with no associated rhytmic genes
phases.list.ld.clean <- phases.list.ld[which(sapply(phases.list.ld, nrow) != 0)]
phases.list.sd.clean <- phases.list.sd[which(sapply(phases.list.sd, nrow) != 0)]

# Test preferential concentration of motifs during the diel cycle
library(circular)

go.circa.table.ld <- complete_circular_table(phases.list.ld.clean,hr.on.large.sets = F, hr.on.large.sets.th=400,
                                             iter.hr = 999, force.hr=F, force.hr.th=0.05, 
                                             rao.on.large.sets = F, rao.on.large.sets.th=400,
                                             iter.rao = 999, force.rao=F, force.rao.th=0.05)
write.table(go.circa.table.ld,file="motif_table_ld.tsv", sep="\t")

go.circa.table.sd <- complete_circular_table(phases.list.sd.clean,hr.on.large.sets = F, hr.on.large.sets.th=400,
                                             iter.hr = 999, force.hr=F, force.hr.th=0.05, 
                                             rao.on.large.sets = F, rao.on.large.sets.th=400,
                                             iter.rao = 999, force.rao=F, force.rao.th=0.05)
write.table(go.circa.table.sd,file="motif_table_sd.tsv", sep="\t")

# Read and explore results
go.circa.table.ld <- read.table("motif_table_ld.tsv", sep="\t", header = T)
go.circa.table.sd <- read.table("motif_table_sd.tsv", sep="\t", header = T)

sum(go.circa.table.ld$rayleigh_p_value_adj < 0.05 | go.circa.table.ld$kuiper_p_value_adj < 0.05) # 427
sum(go.circa.table.sd$rayleigh_p_value_adj < 0.05 | go.circa.table.sd$kuiper_p_value_adj < 0.05) # 419
sum(go.circa.table.ld$hr_p_value_adj < 0.05, na.rm = T) # 0
sum(go.circa.table.sd$hr_p_value_adj < 0.05, na.rm = T) # 0

unimod_ld <- rownames(go.circa.table.ld)[which(go.circa.table.ld$rayleigh_p_value_adj < 0.05 | go.circa.table.ld$kuiper_p_value_adj < 0.05)]
unimod_sd <- rownames(go.circa.table.sd)[which(go.circa.table.sd$rayleigh_p_value_adj < 0.05 | go.circa.table.sd$kuiper_p_value_adj < 0.05)]

#multimod_ld <- rownames(go.circa.table.ld)[which(go.circa.table.ld$hr_p_value_adj < 0.05)]
#multimod_sd <- rownames(go.circa.table.sd)[which(go.circa.table.sd$hr_p_value_adj < 0.05)]

library(VennDiagram)
png("ld_sd_unimod_motif.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(unimod_ld), area2 = length(unimod_sd),
                   cross.area = length(intersect(sapply(strsplit(unimod_ld, "_"), function(x) x[[1]]),
                                                 sapply(strsplit(unimod_sd, "_"), function(x) x[[1]]))),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# png("ld_sd_multimod.png", res=300, width = 1500, height = 1500)
# grid.newpage()
# draw.pairwise.venn(area1 = length(multimod_ld), area2 = length(multimod_sd),
#                    cross.area = length(intersect(sapply(strsplit(multimod_ld, "_"), function(x) x[[1]]),
#                                                  sapply(strsplit(multimod_sd, "_"), function(x) x[[1]]))),
#                    lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
#                    alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
# dev.off()

# For the TFBS groups not following an uniform distribution, test them against the general distribution of phases
# per photoperiod
non_uniform_gos_ld <- unique(c(unimod_ld, multimod_ld))
phases.list.ld.non.uni <- phases.list.ld.clean[non_uniform_gos_ld]
gen.dist.ld <- test_against_gen_dist(phase.list = phases.list.ld.non.uni, total.phase.table = total_phases_table_ld)
write.table(gen.dist.ld,file="motif_general_dist_ld.tsv", sep="\t")

non_uniform_gos_sd <- unique(c(unimod_sd, multimod_sd))
phases.list.sd.non.uni <- phases.list.sd.clean[non_uniform_gos_sd]
gen.dist.sd <- test_against_gen_dist(phase.list = phases.list.sd.non.uni, total.phase.table = total_phases_table_sd)
write.table(gen.dist.sd,file="motif_general_dist_sd.tsv", sep="\t")

gen.dist.ld <- read.table("motif_general_dist_ld.tsv", sep="\t", header = T)
gen.dist.sd <- read.table("motif_general_dist_sd.tsv", sep="\t", header = T)

sum(gen.dist.ld$gen_dist_p_value_adj < 0.05) # 27 non-uniform motif terms present a different distribution than that of global
# mRNA accumulation in LD
sum(gen.dist.sd$gen_dist_p_value_adj < 0.05) # 11 non-uniform motif terms present a different distribution than that of global
# mRNA accumulation in SD

go_dist_gen_ld <- rownames(gen.dist.ld)[which(gen.dist.ld$gen_dist_p_value_adj < 0.05)]
go_dist_gen_sd <- rownames(gen.dist.sd)[which(gen.dist.sd$gen_dist_p_value_adj < 0.05)]

phases.list.sd.clean["bZIP44_SD"]
phases.list.sd.clean["ebox_SD"]
phases.list.sd.clean["ABI5_SD"]
phases.list.sd.clean["abf1_SD"]

png("ld_sd_general_dist_motif.png", res=300, width = 1500, height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(go_dist_gen_ld), area2 = length(go_dist_gen_sd),
                   cross.area = length(intersect(sapply(strsplit(go_dist_gen_ld,"_"), function(x) x[[1]]),
                                                 sapply(strsplit(go_dist_gen_sd,"_"), function(x) x[[1]]))),
                   lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c("blue","red"),
                   alpha = 0.3,cex = 2,cat.cex = 2,cat.dist = c(0.035,0.035), cat.pos = c(210,150))
dev.off()

# Test shared TFBS groups not following an uniform distribution for differences by photoperiod
shared_non_uniform <- intersect(sapply(strsplit(non_uniform_gos_ld, "_"), function(x) x[[1]]),
                                sapply(strsplit(non_uniform_gos_sd, "_"), function(x) x[[1]]))
phases.list.shared.ld.non.uni <- phases.list.ld.non.uni[paste0(shared_non_uniform, "_LD")]
phases.list.shared.sd.non.uni <- phases.list.sd.non.uni[paste0(shared_non_uniform, "_SD")]

go_diff_ld_sd <- mapply(function(x, y) summary(manova(cbind(sin(c(x$phase*pi/12, y$phase*pi/12)),cos(c(x$phase*pi/12, y$phase*pi/12))) 
                                                      ~ c(rep("ld", nrow(x)), rep("sd", nrow(y)))))$stats[1,6],
                        phases.list.shared.ld.non.uni, phases.list.shared.sd.non.uni)

# Remove NULL results (due to a single phase in the set of one of them)
go_diff_ld_sd_clean <- go_diff_ld_sd[names(which(sapply(go_diff_ld_sd, length) != 0))]
go_diff_ld_sd_clean_bh <- p.adjust(unlist(go_diff_ld_sd_clean), method = "BH")

go_diff_ld_sd_table <- data.frame(motif_diff_p_value=unlist(go_diff_ld_sd_clean), 
                                  motif_diff_p_value_adj=go_diff_ld_sd_clean_bh)
rownames(go_diff_ld_sd_table) <- sapply(strsplit(rownames(go_diff_ld_sd_table),"_"), function(x) x[[1]])
go_diff_ld_sd_table["RVE1",]
go_diff_ld_sd_table["At1g14580",]

write.table(go_diff_ld_sd_table,file="go_diff_ld_sd_table_motif.tsv", sep="\t")

go_diff_ld_sd_table <- read.table("go_diff_ld_sd_table_motif.tsv", sep="\t", header = T)

sum(go_diff_ld_sd_table$motif_diff_p_value_adj < 0.05) # 391 change their distribution between LD and SD
change_gos <- rownames(go_diff_ld_sd_table)[which(go_diff_ld_sd_table$motif_diff_p_value_adj < 0.05)]

# Difference in medians for changing GOs
diff_gos_medians <- go.circa.table.ld[change_gos,"median"] - go.circa.table.sd[change_gos,"median"]
names(diff_gos_medians) <- change_gos
sd_anticipated <- change_gos[which((0 < diff_gos_medians & diff_gos_medians < 12) | (diff_gos_medians < -12))]
# 376 groups are anticipated in SD

library(circlize)
library(ComplexHeatmap)
library(MetBrewer)

circ_ex_ld <- go.circa.table.ld[go_dist_gen_ld,]
circular_boxplot(circ_ex_ld[1:14,], "prueba_motifs_ld_diff_1.png", color.palette = "Tam",
                 tr.height = 0.06)
circular_boxplot(circ_ex_ld[15:27,], "prueba_motifs_ld_diff_2.png", color.palette = "Tam",
                 tr.height = 0.06)

circ_ex_sd <- go.circa.table.sd[go_dist_gen_sd,]
circular_boxplot(circ_ex_sd, "prueba_motifs_sd_diff.png", color.palette = "Tam",
                 tr.height = 0.075)
# bZIP plots
circ_ex_ld <- go.circa.table.ld[non_uniform_gos_ld[grep("bZIP", non_uniform_gos_ld)],]
circular_boxplot(circ_ex_ld, "motif_bzip_ld.png", color.palette = "Tam",
                 tr.height = 0.07)

circ_ex_sd <- go.circa.table.sd[non_uniform_gos_sd[grep("bZIP", non_uniform_gos_sd)],]
circular_boxplot(circ_ex_sd, "motif_bzip_sd.png", color.palette = "Tam",
                 tr.height = 0.07)
