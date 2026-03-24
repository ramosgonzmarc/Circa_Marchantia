#BiocManager::install("GenomicFeatures")
#BiocManager::install("txdbmaker")
library(GenomicFeatures)
library(txdbmaker)
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

# Response to light
ld_light <- read.table("../tables/LD_light_activated_dark_repressed.tsv", sep="\t", header = T)$x
sd_light <- read.table("../tables/SD_light_activated_dark_repressed.tsv", sep="\t", header = T)$x

inter_light <- intersect(ld_light, sd_light)
promoters_inter_light <- mp_promoters[inter_light]
write.fasta(getSequence(promoters_inter_light), getName(promoters_inter_light), 
            "promoters_inter_light_response.fa")

promoters_ld_light <- mp_promoters[ld_light]
write.fasta(getSequence(promoters_ld_light), getName(promoters_ld_light), 
            "promoters_ld_light_response.fa")

promoters_sd_light <- mp_promoters[sd_light]
write.fasta(getSequence(promoters_sd_light), getName(promoters_sd_light), 
            "promoters_sd_light_response.fa")

# Response to dark
ld_dark <- read.table("../tables/LD_light_repressed_dark_activated.tsv", sep="\t", header = T)$x
sd_dark <- read.table("../tables/SD_light_repressed_dark_activated.tsv", sep="\t", header = T)$x

inter_dark <- intersect(ld_dark, sd_dark)
promoters_inter_dark <- mp_promoters[inter_dark]
write.fasta(getSequence(promoters_inter_dark), getName(promoters_inter_dark), 
            "promoters_inter_dark_response.fa")

promoters_ld_dark <- mp_promoters[ld_dark]
write.fasta(getSequence(promoters_ld_dark), getName(promoters_ld_dark), 
            "promoters_ld_dark_response.fa")

promoters_sd_dark <- mp_promoters[sd_dark]
write.fasta(getSequence(promoters_sd_dark), getName(promoters_sd_dark), 
            "promoters_sd_dark_response.fa")



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
names(which(sapply(list_ld, function(x) "Mp1g16800" %in% x)))



files_sd <- sapply(strsplit(list.files("results_circa_genescarab/sd_circa_genescarab/"), ".txt"), function(x) x[[1]])
list_sd <- list()

for (x in files_sd)
{
  new_table <- read.table(paste0("results_circa_genescarab/sd_circa_genescarab/", x, ".txt"), header = T, sep="\t")
  new_entry <- unique(subset(new_table, MotifScore > 9)$FASTA.ID)
  list_sd[[x]] <- new_entry
}

# GeneSCARAB for differential analysis
source("../../../../GeneSCARAB/PSEA_functions_package.R")

# Load phases
circa_table_genescarab <- read.table("../../circa_marchantia_paper/tables/circacompare_intersection.tsv", header = T, sep = "\t")
total_phases_table_ld <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["ld.peak.time.hours"]]))
total_phases_table_sd <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["sd.peak.time.hours"]]))

# Generate phase lists for LD and SD
phases.list.ld <- gene.list.to.phases(list_ld, total_phases_table_ld)
phases.list.sd <- gene.list.to.phases(list_sd, total_phases_table_sd)

# Remove entries with no associated rhythmic genes in both photoperiods
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
go.circa.table.ld <- read.table("motif_tables/motif_table_ld_9.tsv", sep="\t", header = T)
go.circa.table.sd <- read.table("motif_tables/motif_table_sd_9.tsv", sep="\t", header = T)

sum(go.circa.table.ld$rayleigh_p_value_adj < 0.05 | go.circa.table.ld$kuiper_p_value_adj < 0.05) # 427
sum(go.circa.table.sd$rayleigh_p_value_adj < 0.05 | go.circa.table.sd$kuiper_p_value_adj < 0.05) # 419
sum(go.circa.table.ld$hr_p_value_adj < 0.05, na.rm = T) # 0
sum(go.circa.table.sd$hr_p_value_adj < 0.05, na.rm = T) # 2
subset(go.circa.table.sd, hr_p_value_adj < 0.05)

unimod_ld <- rownames(go.circa.table.ld)[which(go.circa.table.ld$rayleigh_p_value_adj < 0.05 | go.circa.table.ld$kuiper_p_value_adj < 0.05)]
unimod_sd <- rownames(go.circa.table.sd)[which(go.circa.table.sd$rayleigh_p_value_adj < 0.05 | go.circa.table.sd$kuiper_p_value_adj < 0.05)]

write.table(sapply(strsplit(unimod_ld, "_LD"), function(x) x[1]), 
            "motif_tables/motifs_unimodal_ld.tsv", sep = "\t", row.names = F)
write.table(sapply(strsplit(unimod_sd, "_SD"), function(x) x[1]), 
            "motif_tables/motifs_unimodal_sd.tsv", sep = "\t", row.names = F)


multimod_ld <- rownames(go.circa.table.ld)[which(go.circa.table.ld$hr_p_value_adj < 0.05)]
multimod_sd <- rownames(go.circa.table.sd)[which(go.circa.table.sd$hr_p_value_adj < 0.05)]
go.circa.table.sd["E2FA_SD",]


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

gen.dist.ld <- read.table("motif_tables/motif_general_dist_ld.tsv", sep="\t", header = T)
gen.dist.sd <- read.table("motif_tables/motif_general_dist_sd.tsv", sep="\t", header = T)

sum(gen.dist.ld$gen_dist_p_value_adj < 0.05) # 27 non-uniform motif terms present a different distribution than that of global
# mRNA accumulation in LD
sum(gen.dist.sd$gen_dist_p_value_adj < 0.05) # 11 non-uniform motif terms present a different distribution than that of global
# mRNA accumulation in SD

go_dist_gen_ld <- rownames(gen.dist.ld)[which(gen.dist.ld$gen_dist_p_value_adj < 0.05)]
go_dist_gen_sd <- rownames(gen.dist.sd)[which(gen.dist.sd$gen_dist_p_value_adj < 0.05)]

write.table(sapply(strsplit(go_dist_gen_ld, "_LD"), function(x) x[1]), 
            "motif_tables/motif_diff_general_ld.tsv", sep = "\t", row.names = F)
write.table(sapply(strsplit(go_dist_gen_sd, "_SD"), function(x) x[1]), 
            "motif_tables/motif_diff_general_sd.tsv", sep = "\t", row.names = F)


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

go_diff_ld_sd_table <- read.table("motif_tables/go_diff_ld_sd_table_motif.tsv", sep="\t", header = T)

sum(go_diff_ld_sd_table$motif_diff_p_value_adj < 0.05) # 391 change their distribution between LD and SD
change_gos <- rownames(go_diff_ld_sd_table)[which(go_diff_ld_sd_table$motif_diff_p_value_adj < 0.05)]

# Difference in medians for changing GOs
diff_gos_medians <- go.circa.table.ld[change_gos,"median"] - go.circa.table.sd[change_gos,"median"]
names(diff_gos_medians) <- change_gos
sd_anticipated <- change_gos[which((0 < diff_gos_medians & diff_gos_medians < 12) | (diff_gos_medians < -12))]
# 376 groups are anticipated in SD
setdiff(change_gos, sd_anticipated)


# Plots for non-uniform non-experimentally distributed genes by TFBS

# In both photoperiods
library(circlize)
library(ComplexHeatmap)
library(MetBrewer)

dist_gen_shared <- intersect(unlist(strsplit(go_dist_gen_ld, "_LD")),
                             unlist(strsplit(go_dist_gen_sd, "_SD")))


circular_boxplot(go.circa.table.ld[paste0(dist_gen_shared, "_LD"),],
                 "non_exp_motifs_ld_shared.png", color.palette = "Tam",
                 tr.height = 0.08, width = 1600, height = 1600, res = 180)

circular_boxplot(go.circa.table.sd[paste0(dist_gen_shared, "_SD"),],
                 "non_exp_motifs_sd_shared.png", color.palette = "Tam",
                 tr.height = 0.08, width = 1600, height = 1600, res = 180)


circular_histogram(phases.list.ld[paste0(dist_gen_shared, "_LD")],
                   "non_exp_motifs_ld_shared_hist.png", color.palette = "Tam",
                   tr.height = 0.08, nbins = 48, width = 1600, height = 1600, res = 180)

circular_histogram(phases.list.sd[paste0(dist_gen_shared, "_SD")],
                   "non_exp_motifs_sd_shared_hist.png", color.palette = "Tam",
                   tr.height = 0.08, nbins = 48, width = 1600, height = 1600, res = 180)

# Only in LD
dist_gen_ld_specific <- setdiff(unlist(strsplit(go_dist_gen_ld, "_LD")),
                             unlist(strsplit(go_dist_gen_sd, "_SD")))


circular_boxplot(go.circa.table.ld[paste0(dist_gen_ld_specific, "_LD"),],
                 "non_exp_motifs_ld_specific.png", color.palette = "Tam",
                 tr.height = 0.04)

circular_histogram(phases.list.ld[paste0(dist_gen_ld_specific, "_LD")],
                   "non_exp_motifs_ld_specific_hist.png", color.palette = "Tam",
                   tr.height = 0.04, nbins = 48)


# Only in SD
dist_gen_sd_specific <- setdiff(unlist(strsplit(go_dist_gen_sd, "_SD")),
                                unlist(strsplit(go_dist_gen_ld, "_LD")))


circular_boxplot(go.circa.table.sd[paste0(dist_gen_sd_specific, "_SD"),],
                 "non_exp_motifs_sd_specific.png", color.palette = "Tam",
                 tr.height = 0.10)

circular_histogram(phases.list.sd[paste0(dist_gen_sd_specific, "_SD")],
                   "non_exp_motifs_sd_specific_hist.png", color.palette = "Tam",
                   tr.height = 0.15, nbins = 48)


# Create plots for genes with few associated TFBS

# Create a table showing the distribution of promoters by gene
names(list_ld) <- sapply(strsplit(names(list_ld), "_LD"), function(x) x[1])
names(list_sd) <- sapply(strsplit(names(list_sd), "_SD"), function(x) x[1])

complete_list <- c(list_ld, list_sd)

unique_genes <- unique(unlist(complete_list))

motifs_by_genes <- sapply(complete_list, function(vec) unique_genes %in% vec)
rownames(motifs_by_genes) <- unique_genes

# Combine columns with OR
if (any(duplicated(names(complete_list)))) {
  unique_names <- unique(names(complete_list))
  motifs_by_genes <- sapply(unique_names, function(n) {
    rows <- motifs_by_genes[, names(complete_list) == n, drop = FALSE]
    apply(rows, 1, any)
  })
}

motifs_by_genes_df <- as.data.frame(motifs_by_genes)

motifs_by_genes_df[1:6,1:6]
dim(motifs_by_genes_df)

write.table(motifs_by_genes_df,file="motif_tables/motifs_by_genes_9.tsv", sep="\t")

# Read table and subset
motifs_by_genes_df <- read.table("motif_tables/motifs_by_genes_9.tsv", sep = "\t", header = T)

# Subset by bZIP48 and bZIP28
bzip48_df <- subset(motifs_by_genes_df, bZIP48)

# Subset by low number of motifs in promoter
few_motifs_bzip48 <- bzip48_df[apply(bzip48_df, MARGIN = 1, FUN = sum) < 19, ]
dim(few_motifs_bzip48)

# Recompute tables for plots
circa_table_genescarab <- read.table("../../circa_marchantia_paper/tables/circacompare_intersection.tsv", header = T, sep = "\t")
total_phases_table_ld <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["ld.peak.time.hours"]]))
total_phases_table_sd <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["sd.peak.time.hours"]]))

list_ld_red <- lapply(list_ld, function(x) x[x %in% rownames(few_motifs_bzip48)])
list_sd_red <- lapply(list_sd, function(x) x[x %in% rownames(few_motifs_bzip48)])

phases.list.ld <- gene.list.to.phases(list_ld_red, total_phases_table_ld)
phases.list.sd <- gene.list.to.phases(list_sd_red, total_phases_table_sd)

# Remove entries with no associated rhythmic genes in both photoperiods
phases.list.ld.clean <- phases.list.ld[which(sapply(phases.list.ld, nrow) != 0)]
phases.list.sd.clean <- phases.list.sd[which(sapply(phases.list.sd, nrow) != 0)]

# Recompute shared plot
# go.circa.table.ld <- complete_circular_table(phases.list.ld.clean,hr.on.large.sets = F, hr.on.large.sets.th=400,
#                                              iter.hr = 999, force.hr=F, force.hr.th=0.05, 
#                                              rao.on.large.sets = F, rao.on.large.sets.th=400,
#                                              iter.rao = 999, force.rao=F, force.rao.th=0.05)
# 
# go.circa.table.sd <- complete_circular_table(phases.list.sd.clean,hr.on.large.sets = F, hr.on.large.sets.th=400,
#                                              iter.hr = 999, force.hr=F, force.hr.th=0.05, 
#                                              rao.on.large.sets = F, rao.on.large.sets.th=400,
#                                              iter.rao = 999, force.rao=F, force.rao.th=0.05)


circular_histogram(phases.list.ld.clean[dist_gen_shared[!is.na(names(phases.list.ld.clean[dist_gen_shared]))]],
                   "non_exp_motifs_ld_shared_hist_red_3.png", color.palette = "Tam",
                   tr.height = 0.08, nbins = 24)

circular_dotplot(phases.list.ld.clean[dist_gen_shared[!is.na(names(phases.list.ld.clean[dist_gen_shared]))]],
                   "non_exp_motifs_ld_shared_dot_red_3.png", color.palette = "Tam",
                   tr.height = 0.08)

circular_histogram(phases.list.sd.clean[dist_gen_shared[!is.na(names(phases.list.sd.clean[dist_gen_shared]))]],
                   "non_exp_motifs_sd_shared_hist_red_3.png", color.palette = "Tam",
                   tr.height = 0.08, nbins = 24)

circular_dotplot(phases.list.sd.clean[dist_gen_shared[!is.na(names(phases.list.sd.clean[dist_gen_shared]))]],
                   "non_exp_motifs_sd_shared_dot_red_3.png", color.palette = "Tam",
                   tr.height = 0.08)

################# VIOLIN PLOTS OF TFS EXPRESSION AND PROMOTERS ############

genes.expression <- read.table(file ="../gene_expression_marchantia.tsv",sep="\t",header=T,as.is=T)
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

# Define general variables
ld.gene.expression <- genes.expression[,ld.zt.i]
sd.gene.expression <- genes.expression[,sd.zt.i]
ld.sd.gene.expression <- genes.expression[,c(ld.zt.i,sd.zt.i)]

ld.dd.gene.expression <- genes.expression[,dd.ld.zt.i]
sd.dd.gene.expression <- genes.expression[,dd.sd.zt.i]
ld.ll.gene.expression <- genes.expression[,ll.ld.zt.i]
sd.ll.gene.expression <- genes.expression[,ll.sd.zt.i]

# Prepare scaled variables
scaled.sd.gene.expression <- t(scale(t(sd.gene.expression)))
scaled.ld.gene.expression <- t(scale(t(ld.gene.expression)))
scaled.ld.gene.expression <- scaled.ld.gene.expression[!is.nan(scaled.ld.gene.expression[,1]),]
scaled.sd.gene.expression <- scaled.sd.gene.expression[!is.nan(scaled.sd.gene.expression[,1]),]

# Define variables for goi
hy5.gene <- c("Mp1g16800")
det1.gene <- c("Mp8g00970")
prr.gene <- c("Mp2g16560")
rve.gene <- c("Mp2g14310")
toc1.gene <- c("Mp3g09460")
pif.gene <- c("Mp3g17350")



# Define subset og genes exhibiting motifs of interest
motifs_by_genes_df <- read.table("motif_tables/motifs_by_genes_9.tsv", sep = "\t", header = T)

bzip44_df <- subset(motifs_by_genes_df, bZIP48)
bzip44_df <- bzip44_df[apply(bzip44_df, MARGIN = 1, FUN = sum) < 17, ]
colnames(bzip44_df)[rev(order(colSums(bzip44_df, )))]
bzip44.genes <- intersect(intersect(rownames(bzip44_df), rownames(scaled.ld.gene.expression)),
                          rownames(scaled.sd.gene.expression))


cca1_df <- subset(motifs_by_genes_df, cca)
colnames(cca1_df)[rev(order(colSums(cca1_df, )))]
cca1_df <- cca1_df[apply(cca1_df, MARGIN = 1, FUN = sum) < 12, ]
cca1.genes <- intersect(intersect(rownames(cca1_df), rownames(scaled.ld.gene.expression)),
                        rownames(scaled.sd.gene.expression))



# Si filtramos sólo por genes rítmicos en ambos fotoperiodos
complete.ld.rhythmic.genes <- read.table("../tables/ld.rhythmic.table.tsv", header = T)$V1
complete.sd.rhythmic.genes <- read.table("../tables/sd.rhythmic.table.tsv", header = T)$V1
rhythm <- unique(c(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes))
bzip44.genes <- bzip44.genes[bzip44.genes %in% rhythm]
ebox.genes <- ebox.genes[ebox.genes %in% rhythm]
cca1.genes <- cca1.genes[cca1.genes %in% rhythm]

# Si filtro sólo por los que siguen ciclando en alguna de las condiciones continuas
complete.sd.dd.rhythmic.genes <- read.table("../tables/sd.dd.rhythmic.table.tsv", sep="\t")$x
complete.ld.dd.rhythmic.genes <- read.table("../tables/ld.dd.rhythmic.table.tsv", sep="\t")$x
complete.sd.ll.rhythmic.genes <- read.table("../tables/sd.ll.rhythmic.table.tsv", sep="\t")$x
complete.ld.ll.rhythmic.genes <- read.table("../tables/ld.ll.rhythmic.table.tsv", sep="\t")$x
union.free <- unique(c(complete.sd.dd.rhythmic.genes, complete.sd.ll.rhythmic.genes,
                       complete.ld.dd.rhythmic.genes,complete.ld.ll.rhythmic.genes))

# bzip44.genes <- bzip44.genes[bzip44.genes %in% union.free]
# ebox.genes <- ebox.genes[ebox.genes %in% union.free]
# cca1.genes <- cca1.genes[cca1.genes %in% union.free]

# Plot LD violin
# Scale and prepare expression for plotting
{
scaled.ld.gene.expression.bzip44 <- scaled.ld.gene.expression[bzip44.genes,]
scaled.ld.gene.expression.bzip44.mean <- colMeans(scaled.ld.gene.expression.bzip44)
scaled.ld.gene.expression.bzip44.mean.one.day <- tapply(scaled.ld.gene.expression.bzip44.mean, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.ld.gene.expression.bzip44.mean.one.day <- (scaled.ld.gene.expression.bzip44.mean.one.day - 
                                                   mean(scaled.ld.gene.expression.bzip44.mean.one.day))/sd(scaled.ld.gene.expression.bzip44.mean.one.day)

scaled.ld.gene.expression.cca1 <- scaled.ld.gene.expression[cca1.genes,]
scaled.ld.gene.expression.cca1.mean <- colMeans(scaled.ld.gene.expression.cca1)
scaled.ld.gene.expression.cca1.mean.one.day <- tapply(scaled.ld.gene.expression.cca1.mean, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.ld.gene.expression.cca1.mean.one.day <- (scaled.ld.gene.expression.cca1.mean.one.day - 
                                                  mean(scaled.ld.gene.expression.cca1.mean.one.day))/sd(scaled.ld.gene.expression.cca1.mean.one.day)

scaled.ld.gene.expression.hy5 <- scaled.ld.gene.expression[hy5.gene,]
scaled.ld.gene.expression.hy5.mean.one.day <- tapply(scaled.ld.gene.expression.hy5, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.ld.gene.expression.hy5.mean.one.day <- (scaled.ld.gene.expression.hy5.mean.one.day - 
                                                    mean(scaled.ld.gene.expression.hy5.mean.one.day))/sd(scaled.ld.gene.expression.hy5.mean.one.day)

scaled.ld.gene.expression.det1 <- scaled.ld.gene.expression[det1.gene,]
scaled.ld.gene.expression.det1.mean.one.day <- tapply(scaled.ld.gene.expression.det1, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.ld.gene.expression.det1.mean.one.day <- (scaled.ld.gene.expression.det1.mean.one.day - 
                                                 mean(scaled.ld.gene.expression.det1.mean.one.day))/sd(scaled.ld.gene.expression.det1.mean.one.day)

scaled.ld.gene.expression.rve <- scaled.ld.gene.expression[rve.gene,]
scaled.ld.gene.expression.rve.mean.one.day <- tapply(scaled.ld.gene.expression.rve, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.ld.gene.expression.rve.mean.one.day <- (scaled.ld.gene.expression.rve.mean.one.day - 
                                                 mean(scaled.ld.gene.expression.rve.mean.one.day))/sd(scaled.ld.gene.expression.rve.mean.one.day)

scaled.ld.gene.expression.prr <- scaled.ld.gene.expression[prr.gene,]
scaled.ld.gene.expression.prr.mean.one.day <- tapply(scaled.ld.gene.expression.prr, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.ld.gene.expression.prr.mean.one.day <- (scaled.ld.gene.expression.prr.mean.one.day - 
                                                 mean(scaled.ld.gene.expression.prr.mean.one.day))/sd(scaled.ld.gene.expression.prr.mean.one.day)

scaled.ld.gene.expression.toc1 <- scaled.ld.gene.expression[toc1.gene,]
scaled.ld.gene.expression.toc1.mean.one.day <- tapply(scaled.ld.gene.expression.toc1, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.ld.gene.expression.toc1.mean.one.day <- (scaled.ld.gene.expression.toc1.mean.one.day - 
                                                 mean(scaled.ld.gene.expression.toc1.mean.one.day))/sd(scaled.ld.gene.expression.toc1.mean.one.day)

scaled.ld.gene.expression.pif <- scaled.ld.gene.expression[pif.gene,]
scaled.ld.gene.expression.pif.mean.one.day <- tapply(scaled.ld.gene.expression.pif, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.ld.gene.expression.pif.mean.one.day <- (scaled.ld.gene.expression.pif.mean.one.day - 
                                                  mean(scaled.ld.gene.expression.pif.mean.one.day))/sd(scaled.ld.gene.expression.pif.mean.one.day)

#Define complete set of variables for violin plot
genes.scaled.profile.hy5 <- as.numeric(scaled.ld.gene.expression.hy5.mean.one.day)
genes.scaled.profile.det1 <- as.numeric(scaled.ld.gene.expression.det1.mean.one.day)
genes.scaled.profile.rve <- as.numeric(scaled.ld.gene.expression.rve.mean.one.day)
genes.scaled.profile.prr <- as.numeric(scaled.ld.gene.expression.prr.mean.one.day)
genes.scaled.profile.toc1 <- as.numeric(scaled.ld.gene.expression.toc1.mean.one.day)
genes.scaled.profile.pif <- as.numeric(scaled.ld.gene.expression.pif.mean.one.day)
genes.scaled.profile.bzip44 <- as.numeric(scaled.ld.gene.expression.bzip44.mean.one.day)
genes.scaled.profile.cca1 <- as.numeric(scaled.ld.gene.expression.cca1.mean.one.day)

genes.hy5.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.hy5,genes.scaled.profile.hy5[1]),n = 10000)
genes.det1.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.det1,genes.scaled.profile.det1[1]),n = 10000)
genes.rve.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.rve,genes.scaled.profile.rve[1]),n = 10000)
genes.prr.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.prr,genes.scaled.profile.prr[1]),n = 10000)
genes.toc1.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.toc1,genes.scaled.profile.toc1[1]),n = 10000)
genes.pif.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.pif,genes.scaled.profile.pif[1]),n = 10000)
genes.bzip44.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.bzip44,genes.scaled.profile.bzip44[1]),n = 10000)
genes.cca1.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.cca1,genes.scaled.profile.cca1[1]),n = 10000)


if(min(genes.hy5.splines$y) < 0)
{
  genes.hy5.splines$y <- genes.hy5.splines$y - min(genes.hy5.splines$y)
}

if(min(genes.det1.splines$y) < 0)
{
  genes.det1.splines$y <- genes.det1.splines$y - min(genes.det1.splines$y)
}

if(min(genes.rve.splines$y) < 0)
{
  genes.rve.splines$y <- genes.rve.splines$y - min(genes.rve.splines$y)
}

if(min(genes.prr.splines$y) < 0)
{
  genes.prr.splines$y <- genes.prr.splines$y - min(genes.prr.splines$y)
}

if(min(genes.toc1.splines$y) < 0)
{
  genes.toc1.splines$y <- genes.toc1.splines$y - min(genes.toc1.splines$y)
}

if(min(genes.pif.splines$y) < 0)
{
  genes.pif.splines$y <- genes.pif.splines$y - min(genes.pif.splines$y)
}

if(min(genes.bzip44.splines$y) < 0)
{
  genes.bzip44.splines$y <- genes.bzip44.splines$y - min(genes.bzip44.splines$y)
}

if(min(genes.cca1.splines$y) < 0)
{
  genes.cca1.splines$y <- genes.cca1.splines$y - min(genes.cca1.splines$y)
}
}

library(MetBrewer)
library(RColorBrewer)
colors.ld <- c(brewer.pal(9, "Blues")[2:3], met.brewer("OKeeffe1",11)[7], brewer.pal(9, "Blues")[4],met.brewer("OKeeffe1",11)[8:11])

Cairo::CairoPNG("violin_ld_tfs_with_pif.png", height = 1200, width = 1200, res = 200)

{
par(lwd=2)
i <- 0
plot(x = genes.cca1.splines$x,y = genes.cca1.splines$y,type="l",ylim=c(-8,50),axes=F,ylab="",xlab="",lwd=2, xlim = c(0,24))
lines(x = genes.cca1.splines$x, y = -genes.cca1.splines$y,lwd=2)
polygon(x=c(genes.cca1.splines$x,rev(genes.cca1.splines$x)),y = c(genes.cca1.splines$y,rev(-genes.cca1.splines$y)),col=colors.ld[8])

i <- i+7
lines(x = genes.bzip44.splines$x, y = i+genes.bzip44.splines$y,lwd=2)
lines(x = genes.bzip44.splines$x, y = i-genes.bzip44.splines$y,lwd=2)
polygon(x=c(genes.bzip44.splines$x,rev(genes.bzip44.splines$x)),y = c(i + genes.bzip44.splines$y,rev(i - genes.bzip44.splines$y)),col=colors.ld[7])

i <- i+7
lines(x = genes.pif.splines$x, y = i+genes.pif.splines$y,lwd=2)
lines(x = genes.pif.splines$x, y = i-genes.pif.splines$y,lwd=2)
polygon(x=c(genes.pif.splines$x,rev(genes.pif.splines$x)),y = c(i + genes.pif.splines$y,rev(i - genes.pif.splines$y)),col=colors.ld[6])

i <- i+7
lines(x = genes.hy5.splines$x, y = i+genes.hy5.splines$y,lwd=2)
lines(x = genes.hy5.splines$x, y = i-genes.hy5.splines$y,lwd=2)
polygon(x=c(genes.hy5.splines$x,rev(genes.hy5.splines$x)),y = c(i + genes.hy5.splines$y,rev(i - genes.hy5.splines$y)),col=colors.ld[5])

i <- i+7
lines(x = genes.det1.splines$x, y = i+genes.det1.splines$y,lwd=2)
lines(x = genes.det1.splines$x, y = i-genes.det1.splines$y,lwd=2)
polygon(x=c(genes.det1.splines$x,rev(genes.det1.splines$x)),y = c(i + genes.det1.splines$y,rev(i - genes.det1.splines$y)),col=colors.ld[4])

i <- i+7
lines(x = genes.toc1.splines$x, y = i+genes.toc1.splines$y,lwd=2)
lines(x = genes.toc1.splines$x, y = i-genes.toc1.splines$y,lwd=2)
polygon(x=c(genes.toc1.splines$x,rev(genes.toc1.splines$x)),y = c(i + genes.toc1.splines$y,rev(i - genes.toc1.splines$y)),col=colors.ld[3])

i <- i+7
lines(x = genes.rve.splines$x, y = i+genes.rve.splines$y,lwd=2)
lines(x = genes.rve.splines$x, y = i-genes.rve.splines$y,lwd=2)
polygon(x=c(genes.rve.splines$x,rev(genes.rve.splines$x)),y = c(i + genes.rve.splines$y,rev(i - genes.rve.splines$y)),col=colors.ld[2])

i <- i+7
lines(x = genes.prr.splines$x, y = i+genes.prr.splines$y,lwd=2)
lines(x = genes.prr.splines$x, y = i-genes.prr.splines$y,lwd=2)
polygon(x=c(genes.prr.splines$x,rev(genes.prr.splines$x)),y = c(i + genes.prr.splines$y,rev(i - genes.prr.splines$y)),col=colors.ld[1])

polygon(x = c(0,0,16,16),y=c(-8,-5,-5,-8),border="blue",lwd=2)
polygon(x = c(16,16,24,24),y=c(-8,-5,-5,-8),lwd=2,col="blue",border="blue") #510 420
}

dev.off()



# Plot SD violin
# Scale and prepare expression for plotting 
{
scaled.sd.gene.expression.bzip44 <- scaled.sd.gene.expression[bzip44.genes,]
scaled.sd.gene.expression.bzip44.mean <- colMeans(scaled.sd.gene.expression.bzip44)
scaled.sd.gene.expression.bzip44.mean.one.day <- tapply(scaled.sd.gene.expression.bzip44.mean, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.sd.gene.expression.bzip44.mean.one.day <- (scaled.sd.gene.expression.bzip44.mean.one.day - 
                                                    mean(scaled.sd.gene.expression.bzip44.mean.one.day))/sd(scaled.sd.gene.expression.bzip44.mean.one.day)

scaled.sd.gene.expression.cca1 <- scaled.sd.gene.expression[cca1.genes,]
scaled.sd.gene.expression.cca1.mean <- colMeans(scaled.sd.gene.expression.cca1)
scaled.sd.gene.expression.cca1.mean.one.day <- tapply(scaled.sd.gene.expression.cca1.mean, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.sd.gene.expression.cca1.mean.one.day <- (scaled.sd.gene.expression.cca1.mean.one.day - 
                                                  mean(scaled.sd.gene.expression.cca1.mean.one.day))/sd(scaled.sd.gene.expression.cca1.mean.one.day)

scaled.sd.gene.expression.hy5 <- scaled.sd.gene.expression[hy5.gene,]
scaled.sd.gene.expression.hy5.mean.one.day <- tapply(scaled.sd.gene.expression.hy5, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.sd.gene.expression.hy5.mean.one.day <- (scaled.sd.gene.expression.hy5.mean.one.day - 
                                                 mean(scaled.sd.gene.expression.hy5.mean.one.day))/sd(scaled.sd.gene.expression.hy5.mean.one.day)

scaled.sd.gene.expression.det1 <- scaled.sd.gene.expression[det1.gene,]
scaled.sd.gene.expression.det1.mean.one.day <- tapply(scaled.sd.gene.expression.det1, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.sd.gene.expression.det1.mean.one.day <- (scaled.sd.gene.expression.det1.mean.one.day - 
                                                  mean(scaled.sd.gene.expression.det1.mean.one.day))/sd(scaled.sd.gene.expression.det1.mean.one.day)

scaled.sd.gene.expression.rve <- scaled.sd.gene.expression[rve.gene,]
scaled.sd.gene.expression.rve.mean.one.day <- tapply(scaled.sd.gene.expression.rve, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.sd.gene.expression.rve.mean.one.day <- (scaled.sd.gene.expression.rve.mean.one.day - 
                                                 mean(scaled.sd.gene.expression.rve.mean.one.day))/sd(scaled.sd.gene.expression.rve.mean.one.day)

scaled.sd.gene.expression.prr <- scaled.sd.gene.expression[prr.gene,]
scaled.sd.gene.expression.prr.mean.one.day <- tapply(scaled.sd.gene.expression.prr, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.sd.gene.expression.prr.mean.one.day <- (scaled.sd.gene.expression.prr.mean.one.day - 
                                                 mean(scaled.sd.gene.expression.prr.mean.one.day))/sd(scaled.sd.gene.expression.prr.mean.one.day)

scaled.sd.gene.expression.toc1 <- scaled.sd.gene.expression[toc1.gene,]
scaled.sd.gene.expression.toc1.mean.one.day <- tapply(scaled.sd.gene.expression.toc1, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.sd.gene.expression.toc1.mean.one.day <- (scaled.sd.gene.expression.toc1.mean.one.day - 
                                                  mean(scaled.sd.gene.expression.toc1.mean.one.day))/sd(scaled.sd.gene.expression.toc1.mean.one.day)

scaled.sd.gene.expression.pif <- scaled.sd.gene.expression[pif.gene,]
scaled.sd.gene.expression.pif.mean.one.day <- tapply(scaled.sd.gene.expression.pif, factor(rep(seq(0,20,4), each = 3)), mean, na.rm = TRUE)
scaled.sd.gene.expression.pif.mean.one.day <- (scaled.sd.gene.expression.pif.mean.one.day - 
                                                  mean(scaled.sd.gene.expression.pif.mean.one.day))/sd(scaled.sd.gene.expression.pif.mean.one.day)


#Define complete set of variables for violin plot
genes.scaled.profile.hy5 <- as.numeric(scaled.sd.gene.expression.hy5.mean.one.day)
genes.scaled.profile.det1 <- as.numeric(scaled.sd.gene.expression.det1.mean.one.day)
genes.scaled.profile.rve <- as.numeric(scaled.sd.gene.expression.rve.mean.one.day)
genes.scaled.profile.prr <- as.numeric(scaled.sd.gene.expression.prr.mean.one.day)
genes.scaled.profile.toc1 <- as.numeric(scaled.sd.gene.expression.toc1.mean.one.day)
genes.scaled.profile.pif <- as.numeric(scaled.sd.gene.expression.pif.mean.one.day)
genes.scaled.profile.bzip44 <- as.numeric(scaled.sd.gene.expression.bzip44.mean.one.day)
genes.scaled.profile.cca1 <- as.numeric(scaled.sd.gene.expression.cca1.mean.one.day)

genes.hy5.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.hy5,genes.scaled.profile.hy5[1]),n = 10000)
genes.det1.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.det1,genes.scaled.profile.det1[1]),n = 10000)
genes.rve.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.rve,genes.scaled.profile.rve[1]),n = 10000)
genes.prr.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.prr,genes.scaled.profile.prr[1]),n = 10000)
genes.toc1.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.toc1,genes.scaled.profile.toc1[1]),n = 10000)
genes.pif.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.pif,genes.scaled.profile.pif[1]),n = 10000)
genes.bzip44.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.bzip44,genes.scaled.profile.bzip44[1]),n = 10000)
genes.cca1.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile.cca1,genes.scaled.profile.cca1[1]),n = 10000)


if(min(genes.hy5.splines$y) < 0)
{
  genes.hy5.splines$y <- genes.hy5.splines$y - min(genes.hy5.splines$y)
}

if(min(genes.det1.splines$y) < 0)
{
  genes.det1.splines$y <- genes.det1.splines$y - min(genes.det1.splines$y)
}

if(min(genes.rve.splines$y) < 0)
{
  genes.rve.splines$y <- genes.rve.splines$y - min(genes.rve.splines$y)
}

if(min(genes.prr.splines$y) < 0)
{
  genes.prr.splines$y <- genes.prr.splines$y - min(genes.prr.splines$y)
}

if(min(genes.toc1.splines$y) < 0)
{
  genes.toc1.splines$y <- genes.toc1.splines$y - min(genes.toc1.splines$y)
}

if(min(genes.pif.splines$y) < 0)
{
  genes.pif.splines$y <- genes.pif.splines$y - min(genes.pif.splines$y)
}

if(min(genes.bzip44.splines$y) < 0)
{
  genes.bzip44.splines$y <- genes.bzip44.splines$y - min(genes.bzip44.splines$y)
}

if(min(genes.cca1.splines$y) < 0)
{
  genes.cca1.splines$y <- genes.cca1.splines$y - min(genes.cca1.splines$y)
}
}

library(MetBrewer)
library(RColorBrewer)
colors.sd <- c(brewer.pal(9, "Reds")[c(2,3,4,5,6,7,8,9)])

Cairo::CairoPNG("violin_sd_tfs_with_pif.png", height = 1200, width = 1200, res = 200)

{
par(lwd=2)
i <- 0
plot(x = genes.cca1.splines$x,y = genes.cca1.splines$y,type="l",ylim=c(-8,50),axes=F,ylab="",xlab="",lwd=2, xlim = c(0,24))
lines(x = genes.cca1.splines$x, y = -genes.cca1.splines$y,lwd=2)
polygon(x=c(genes.cca1.splines$x,rev(genes.cca1.splines$x)),y = c(genes.cca1.splines$y,rev(-genes.cca1.splines$y)),col=colors.sd[8])

i <- i+7
lines(x = genes.bzip44.splines$x, y = i+genes.bzip44.splines$y,lwd=2)
lines(x = genes.bzip44.splines$x, y = i-genes.bzip44.splines$y,lwd=2)
polygon(x=c(genes.bzip44.splines$x,rev(genes.bzip44.splines$x)),y = c(i + genes.bzip44.splines$y,rev(i - genes.bzip44.splines$y)),col=colors.sd[7])

i <- i+7
lines(x = genes.pif.splines$x, y = i+genes.pif.splines$y,lwd=2)
lines(x = genes.pif.splines$x, y = i-genes.pif.splines$y,lwd=2)
polygon(x=c(genes.pif.splines$x,rev(genes.pif.splines$x)),y = c(i + genes.pif.splines$y,rev(i - genes.pif.splines$y)),col=colors.sd[6])

i <- i+7
lines(x = genes.hy5.splines$x, y = i+genes.hy5.splines$y,lwd=2)
lines(x = genes.hy5.splines$x, y = i-genes.hy5.splines$y,lwd=2)
polygon(x=c(genes.hy5.splines$x,rev(genes.hy5.splines$x)),y = c(i + genes.hy5.splines$y,rev(i - genes.hy5.splines$y)),col=colors.sd[5])

i <- i+7
lines(x = genes.det1.splines$x, y = i+genes.det1.splines$y,lwd=2)
lines(x = genes.det1.splines$x, y = i-genes.det1.splines$y,lwd=2)
polygon(x=c(genes.det1.splines$x,rev(genes.det1.splines$x)),y = c(i + genes.det1.splines$y,rev(i - genes.det1.splines$y)),col=colors.sd[4])

i <- i+7
lines(x = genes.toc1.splines$x, y = i+genes.toc1.splines$y,lwd=2)
lines(x = genes.toc1.splines$x, y = i-genes.toc1.splines$y,lwd=2)
polygon(x=c(genes.toc1.splines$x,rev(genes.toc1.splines$x)),y = c(i + genes.toc1.splines$y,rev(i - genes.toc1.splines$y)),col=colors.sd[3])

i <- i+7
lines(x = genes.rve.splines$x, y = i+genes.rve.splines$y,lwd=2)
lines(x = genes.rve.splines$x, y = i-genes.rve.splines$y,lwd=2)
polygon(x=c(genes.rve.splines$x,rev(genes.rve.splines$x)),y = c(i + genes.rve.splines$y,rev(i - genes.rve.splines$y)),col=colors.sd[2])

i <- i+7
lines(x = genes.prr.splines$x, y = i+genes.prr.splines$y,lwd=2)
lines(x = genes.prr.splines$x, y = i-genes.prr.splines$y,lwd=2)
polygon(x=c(genes.prr.splines$x,rev(genes.prr.splines$x)),y = c(i + genes.prr.splines$y,rev(i - genes.prr.splines$y)),col=colors.sd[1])

polygon(x = c(0,0,8,8),y=c(-8,-5,-5,-8),border="red",lwd=2)
polygon(x = c(8,8,24,24),y=c(-8,-5,-5,-8),lwd=2,col="red",border="red") #510 420
}

dev.off()

# Enriquecimiento en genes de respuesta a luz
ld_light <- read.table("../tables/LD_light_activated_dark_repressed.tsv", sep="\t", header = T)$x
sd_light <- read.table("../tables/SD_light_activated_dark_repressed.tsv", sep="\t", header = T)$x

light_df <- motifs_by_genes_df[unique(c(ld_light, sd_light)),]
light_motifs <- colnames(light_df)[rev(order(colSums(light_df, )))]

ld_dark <- read.table("../tables/LD_light_repressed_dark_activated.tsv", sep="\t", header = T)$x
sd_dark <- read.table("../tables/SD_light_repressed_dark_activated.tsv", sep="\t", header = T)$x

dark_df <- motifs_by_genes_df[unique(c(ld_dark, sd_dark)),]
dark_motifs <- colnames(dark_df)[rev(order(colSums(dark_df, )))]


setdiff(light_motifs[1:40], dark_motifs[1:40])
setdiff(dark_motifs[1:40], light_motifs[1:40])




############### Exploration of genes associated to each TFBS at each time point ###############
files_ld <- sapply(strsplit(list.files("results_circa_genescarab/ld_circa_genescarab/"), ".txt"), function(x) x[[1]])
list_ld <- list()

for (x in files_ld)
{
  new_table <- read.table(paste0("results_circa_genescarab/ld_circa_genescarab/", x, ".txt"), header = T, sep="\t")
  new_entry <- unique(subset(new_table, MotifScore > 7)$FASTA.ID)
  list_ld[[x]] <- new_entry
}

files_sd <- sapply(strsplit(list.files("results_circa_genescarab/sd_circa_genescarab/"), ".txt"), function(x) x[[1]])
list_sd <- list()

for (x in files_sd)
{
  new_table <- read.table(paste0("results_circa_genescarab/sd_circa_genescarab/", x, ".txt"), header = T, sep="\t")
  new_entry <- unique(subset(new_table, MotifScore > 7)$FASTA.ID)
  list_sd[[x]] <- new_entry
}

names(which(sapply(list_ld, function(x) "Mp8g00970" %in% x)))
names(which(sapply(list_ld, function(x) "Mp1g16800" %in% x)))

ld_complete_circa <- read.table("../tables/circacompare_ld_24.tsv", header = T)
head(ld_complete_circa)
ld_rhythmic <- subset(ld_complete_circa, rhythmic_p < 0.05)
head(ld_rhythmic)

sd_complete_circa <- read.table("../tables/circacompare_sd_24.tsv", header = T)
head(sd_complete_circa)
sd_rhythmic <- subset(sd_complete_circa, rhythmic_p < 0.05)
head(sd_rhythmic)


# Get examples of genes clusters showing a specific motif in their promoter
ld_zt0 <- rownames(ld_rhythmic)[ld_rhythmic$peak_time_hours < 2 | ld_rhythmic$peak_time_hours >= 22]
ld_zt4 <- rownames(ld_rhythmic)[ld_rhythmic$peak_time_hours < 6 & ld_rhythmic$peak_time_hours >= 2]
ld_zt8 <- rownames(ld_rhythmic)[ld_rhythmic$peak_time_hours < 10 & ld_rhythmic$peak_time_hours >= 6]
ld_zt12 <- rownames(ld_rhythmic)[ld_rhythmic$peak_time_hours < 14 & ld_rhythmic$peak_time_hours >= 10]
ld_zt16 <- rownames(ld_rhythmic)[ld_rhythmic$peak_time_hours < 18 & ld_rhythmic$peak_time_hours >= 14]
ld_zt20 <- rownames(ld_rhythmic)[ld_rhythmic$peak_time_hours < 22 & ld_rhythmic$peak_time_hours >= 18]

sd_zt0 <- rownames(sd_rhythmic)[sd_rhythmic$peak_time_hours < 2 | sd_rhythmic$peak_time_hours >= 22]
sd_zt4 <- rownames(sd_rhythmic)[sd_rhythmic$peak_time_hours < 6 & sd_rhythmic$peak_time_hours >= 2]
sd_zt8 <- rownames(sd_rhythmic)[sd_rhythmic$peak_time_hours < 10 & sd_rhythmic$peak_time_hours >= 6]
sd_zt12 <- rownames(sd_rhythmic)[sd_rhythmic$peak_time_hours < 14 & sd_rhythmic$peak_time_hours >= 10]
sd_zt16 <- rownames(sd_rhythmic)[sd_rhythmic$peak_time_hours < 18 & sd_rhythmic$peak_time_hours >= 14]
sd_zt20 <- rownames(sd_rhythmic)[sd_rhythmic$peak_time_hours < 22 & sd_rhythmic$peak_time_hours >= 18]

ld_zt0[ld_zt0 %in% list_ld$bZIP16_LD]

# If we intersect that with a GO of interest
library(org.Mpolymorpha71.eg.db)
columns(org.Mpolymorpha71.eg.db)

photo_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0015979",
                      keytype = "GO", columns = c("GID", "GO"))

ribosome_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0005840",
                      keytype = "GO", columns = c("GID", "GO"))

protein_catabolism_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0030163",
                                   keytype = "GO", columns = c("GID", "GO"))
  
splicing_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0000398",
                                   keytype = "GO", columns = c("GID", "GO"))

chloroplast_org_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0009658",
                                   keytype = "GO", columns = c("GID", "GO"))

transcription_reg_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0006355",
                                keytype = "GO", columns = c("GID", "GO"))

translation_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0006412",
                                  keytype = "GO", columns = c("GID", "GO"))

carbon_fix_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0015977",
                           keytype = "GO", columns = c("GID", "GO"))

light_harv_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0009765",
                           keytype = "GO", columns = c("GID", "GO"))

intra_transport_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0046907",
                                keytype = "GO", columns = c("GID", "GO"))

rna_degr_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0006401",
                         keytype = "GO", columns = c("GID", "GO"))

prot_proc_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0016485",
                          keytype = "GO", columns = c("GID", "GO"))

proteasome_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0010498",
                          keytype = "GO", columns = c("GID", "GO"))

signal_genes <- select(org.Mpolymorpha71.eg.db, keys = "GO:0007165",
                       keytype = "GO", columns = c("GID", "GO"))

# LD
intersect(ribosome_genes$GID, ld_zt0[ld_zt0 %in% list_ld$SPL3_LD]) # SBP Ribosomal S17 Mp5g05080
intersect(chloroplast_org_genes$GID, ld_zt0[ld_zt0 %in% list_ld$bZIP16_LD]) # bZIP CHUP1 Mp7g03540

intersect(ribosome_genes$GID, ld_zt4[ld_zt4 %in% list_ld$AT2G01818_LD]) # PLATZ Mp4g05120 Ribosomal L6
intersect(ribosome_genes$GID, ld_zt4[ld_zt4 %in% list_ld$bZIP16_LD]) # bZIP Mp1g05730 Ribosomal protein L9

intersect(light_harv_genes$GID, ld_zt8[ld_zt8 %in% list_ld$At2g45680_LD]) # TCP Mp7g06780 LHCB2
intersect(light_harv_genes$GID, ld_zt8[ld_zt8 %in% list_ld$bZIP16_LD]) # bzip LHCA4 Mp3g07840
intersect(photo_genes$GID, ld_zt8[ld_zt8 %in% list_ld$bZIP16_LD]) # bzip Mp7g18150 psaL
intersect(translation_genes$GID, ld_zt8[ld_zt8 %in% list_ld$At1g64620_LD]) #  C2C2dof Mp1g01690 Transamidase

intersect(transcription_reg_genes$GID, ld_zt12[ld_zt12 %in% list_ld$NTM2_LD]) # NAC MpVg00340 MpMYCY
intersect(splicing_genes$GID, ld_zt12[ld_zt12 %in% list_ld$At5g08330_LD]) # TCP MpCDC5 Mp1g10310
intersect(proteasome_genes$GID, ld_zt12[ld_zt12 %in% list_ld$TGA5_LD]) # bZIP 20S proteasome subunit alpha 7 Mp8g09510

intersect(transcription_reg_genes$GID, ld_zt16[ld_zt16 %in% list_ld$WRKY22_LD]) # WRKY Mp2g07720 MpBZIP7
intersect(signal_genes$GID, ld_zt16[ld_zt16 %in% list_ld$AGL95_LD]) # ND MpNBS-LRR12 Mp4g21990

intersect(translation_genes$GID, ld_zt20[ld_zt20 %in% list_ld$WRKY15_LD]) # WRKY Ribosomal L18 Mp3g09540
intersect(signal_genes$GID, ld_zt20[ld_zt20 %in% list_ld$MYB93_LD]) # MYB MpPAS4 Mp4g11010
ld_zt20[ld_zt20 %in% list_ld$AT2G15740_LD] # C2H2 Mp2g17860 ALDOLASE

# SD
sd_zt0[sd_zt0 %in% list_sd$At2g45680_SD] # TCP Mp1g03160 F-1,6-BP
intersect(photo_genes$GID, sd_zt0[sd_zt0 %in% list_sd$ABI5_SD]) # bZIP PsbP Mp2g17390

intersect(photo_genes$GID, sd_zt4[sd_zt4 %in% list_sd$At1g72010_SD]) # TCP PsbR Mp3g05840
intersect(light_harv_genes$GID, sd_zt4[sd_zt4 %in% list_sd$ABI5_SD]) # ABI5 LHCB2 Mp7g06790

intersect(transcription_reg_genes$GID, sd_zt8[sd_zt8 %in% list_sd$AT3G25990_SD]) # Trihelix MpVg00340 MpMYCY
intersect(prot_proc_genes$GID, sd_zt8[sd_zt8 %in% list_sd$ANAC071_SD]) # NAC peroxisomal leader peptide-processing protease Mp6g19230
intersect(proteasome_genes$GID, sd_zt8[sd_zt8 %in% list_sd$At1g64620_SD]) # C2C2dof 20S proteasome subunit beta 4 Mp4g04320


intersect(transcription_reg_genes$GID, sd_zt12[sd_zt12 %in% list_sd$WRKY46_SD]) # MpLRL Mp2g20695

intersect(chloroplast_org_genes$GID, sd_zt16[sd_zt16 %in% list_sd$TRP2_SD]) #  MYBrel-AAACCC Thioredoxin Mp1g20510
intersect(chloroplast_org_genes$GID, sd_zt16[sd_zt16 %in% list_sd$MYB93_SD]) # MYB Ribokinase Mp7g17740
intersect(ribosome_genes$GID, sd_zt16[sd_zt16 %in% list_sd$AT1G44830_SD]) # AP2EREBP Ribosomal L26 Mp2g12870


intersect(photo_genes$GID, sd_zt20[sd_zt20 %in% list_sd$At5g65130_SD]) # AP2EREBP Mp5g23000 Psp29

library(enrichplot)
library(clusterProfiler)
prueba <- sd_zt0[sd_zt0 %in% list_sd$At2g45680_SD]
ego.low.exp.genes <- enrichGO(gene = prueba,
                              OrgDb         = org.Mpolymorpha71.eg.db,
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.10,
                              qvalueCutoff  = 0.10,
                              readable      = F,
                              keyType = "GID")

subset(as.data.frame(ego.low.exp.genes), ID == "GO:0005985")


# Free running
ld_ll <- read.table("../tables/ld.ll.rhythmic.table.tsv", sep = "\t")$x
sd_ll <- read.table("../tables/sd.ll.rhythmic.table.tsv", sep = "\t")$x
ll_inter <- intersect(ld_ll, sd_ll)
ll_union <- unique(c(ld_ll, sd_ll))

ld_dd <- read.table("../tables/ld.dd.rhythmic.table.tsv", sep = "\t")$x
sd_dd <- read.table("../tables/sd.dd.rhythmic.table.tsv", sep = "\t")$x
dd_inter <- intersect(ld_dd, sd_dd)
dd_union <- unique(c(ld_dd, sd_dd))

intersect(photo_genes$GID, ll_inter[ll_inter %in% list_ld$At1g72010_LD]) # TCP Mp5g04200 PsaD
intersect(chloroplast_org_genes$GID, ll_inter[ll_inter %in% list_ld$ABI5_LD]) # bZIP CHUP1 Mp7g03540
intersect(photo_genes$GID, ll_inter[ll_inter %in% list_ld$ABI5_LD]) # bZIP PsaL Mp7g18150
intersect(light_harv_genes$GID, ll_inter[ll_inter %in% list_ld$LCL1_LD]) # MYB LHCB5 Mp4g10900

intersect(transcription_reg_genes$GID, dd_inter[dd_inter %in% list_ld$dof24_LD]) # C2C2dof Mp1g16800 HY5
intersect(signal_genes$GID, dd_inter[dd_inter %in% list_ld$AT3G51470_LD]) # SBP MpPAS3 Mp2g15090
dd_inter[dd_inter %in% list_ld$GATA6_LD] # C2C2gata Mp5g12555 intra transport


# Response to light/dark
# Response to light. In summary: ribosomes and translation
ld_light <- read.table("../tables/LD_light_activated_dark_repressed.tsv", sep="\t", header = T)$x
sd_light <- read.table("../tables/SD_light_activated_dark_repressed.tsv", sep="\t", header = T)$x
inter_light <- intersect(ld_light, sd_light)

ld_dark <- read.table("../tables/LD_light_repressed_dark_activated.tsv", sep="\t", header = T)$x
sd_dark <- read.table("../tables/SD_light_repressed_dark_activated.tsv", sep="\t", header = T)$x
inter_dark <- intersect(ld_dark, sd_dark)

intersect(photo_genes$GID, ld_light[ld_light %in% list_ld$MYB51_LD]) # MYB psaG/psaK Mp4g06880

intersect(ribosome_genes$GID, sd_light[sd_light %in% list_sd$DEAR3_SD]) # AP2EREBP CCAGC Ribosomal L5 Mp1g19970
intersect(ribosome_genes$GID, sd_light[sd_light %in% list_sd$RAP21_SD]) # AP2EREBP GTCGG Ribosomal LP2 Mp1g18620

intersect(ribosome_genes$GID, inter_light[inter_light %in% list_ld$IDD2_LD]) # C2H2 MpPRPL Mp1g18940


# Response to dark: splicing and transcrption regulation
intersect(splicing_genes$GID, ld_dark[ld_dark %in% list_ld$AT1G28160_LD]) # AP2EREBP CWC26 Mp7g05140
intersect(transcription_reg_genes$GID, ld_dark[ld_dark %in% list_ld$MYB61_LD]) # MYB Mp5g18910 MpBHLH40


intersect(protein_catabolism_genes$GID, sd_dark[sd_dark %in% list_sd$SGR5_SD]) # C2H2 ClpS Mp7g08950
intersect(transcription_reg_genes$GID, sd_dark[sd_dark %in% list_sd$MYB51_SD]) # MYB MpIAA Mp6g05000

intersect(transcription_reg_genes$GID, inter_dark[inter_dark %in% list_ld$SGR5_LD]) # C2H2 MpBHLH21 Mp3g11900
intersect(transcription_reg_genes$GID, inter_dark[inter_dark %in% list_ld$HDG1_LD]) # Homeobox Mp1g24950 MpWRKY11
intersect(transcription_reg_genes$GID, inter_dark[inter_dark %in% list_ld$MYB118_LD]) # MYB Mp5g18820 MpBHLH39


# SD Period 12h
sd_period_12 <- read.table("../tables/circacompare_sd_12.tsv", sep = "\t")
head(sd_period_12)
sd_period_12 <- rownames(subset(sd_period_12, rhythmic_p < 0.05))

intersect(photo_genes$GID, sd_period_12[sd_period_12 %in% list_sd$GBF5_SD]) #bZIP  Mp6g21460 PsbX
