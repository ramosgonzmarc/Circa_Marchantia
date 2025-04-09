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
