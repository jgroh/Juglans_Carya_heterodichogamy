library(data.table)
library(HardyWeinberg)
library(seqinr)

# This script shows representative commands used to analyse phased G-locus variants from J. hindsii
# similar commands were used for analyses in other species

# After variant calling and phasing across G locus, converted to IMPUTE format using vcftools

# ----- SNP positions ----- 
leg <- fread("~/workspace/heterodichogamy/hindsii/Putah_Gloc_phased.impute.legend")

# ID column in this table is both redundant and unecessary
leg[, ID := NULL]

# ----- haps ----- Each column is a biallelic SNP position. Columns are haplotypes.
haps0 <- fread("~/workspace/heterodichogamy/hindsii/Putah_Gloc_phased.impute.hap")
colnames(haps0)
n <- ncol(haps0)
setnames(haps0, paste0("V", 1:n), paste0(rep(seq(1:(n/2)), each = 2), c("_1", "_2")))

d0 <- melt(cbind(leg, haps0), id.vars = c("pos", "allele0", "allele1"), value.name = 'allele')
d0[, c("indiv.id", "hap.id") := tstrsplit(variable, split = "_", fixed = T)]
d0[, variable := NULL]

# ----- individual IDs -----
indiv <- fread("~/workspace/heterodichogamy/hindsii/Putah_Gloc_phased.impute.hap.indv", header = F, col.names = c("ID"))
indiv[, indiv.id := as.character(seq(1, .N))]

d1 <- merge(d0, indiv, by = 'indiv.id', all = T)
d1[, indiv.id := NULL]


# ---- merge with genotype (assigned by phenotype or coverage) -----
# Relevant data provided in Table S2
pheno <- fread("~/workspace/heterodichogamy/data/phenotypes.txt", col.names = c("ID", "phenotype", "species"))

haps <- merge( pheno[species == 'hindsii' & grepl('JHIN_PC', ID), .(ID, phenotype)], d1, by = 'ID', all = T)
haps[phenotype == 'protogynous', genotype := 'Gg'] # coverage indicates all protogynous individuals are het
haps[phenotype == 'protandrous', genotype := 'gg']
haps[ID == 'JHIN_PC_102', genotype := 'Gg'] # no phenotype, assigned on basis of coverage
haps[ID == 'JHIN_PC_052', genotype := 'gg']


# ----- visualize distribution of nonreference alleles -----
# (expect to be bimodal)
a <- haps[hap.id == 1, .(nonref_alleles_hap1 = sum(allele)), by = .(ID,genotype)]
b <- haps[hap.id == 2, .(nonref_alleles_hap2 = sum(allele)), by = .(ID,genotype)]
nonref_alleles <- merge(a,b)
nonref_allhaps <- melt(nonref_alleles, id.vars = c("ID", 'genotype'), value.name = "Non ref alleles")

ggplot(nonref_allhaps, aes(x = `Non ref alleles`)) + geom_histogram() + 
  theme_classic() + 
  theme(aspect.ratio = 1)



# ----- Assign haplotype ID based on fake SNP that codes for SV -----
haps[, hapRank := seq(1, .N), by = .(pos)]

# This was manually added to the vcf at a previous step to represent each individual's SV genotype
G.IDs <- haps[pos == 31370000 & allele == 0, hapRank]
g.IDs <- haps[pos == 31370000 & allele == 1, hapRank]
haps[hapRank %in% G.IDs, hap := 'G']
haps[hapRank %in% g.IDs, hap := 'g']


# ----- visualize haplotypes -----
haps[, SNPrank := seq(1, .N), by = .(ID, hap.id)]


ggplot(haps, 
       aes(x = SNPrank, y = hapRank)) + 
  facet_wrap(~hap) +
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none') 


# JHIN_PC_099 was identified as potential hybrid, and has highly divergent g haplotype, excluded from sample
# 
# ggplot(haps[hapRank!= 57], 
#        aes(x = SNPrank, y = hapRank)) + 
#   facet_wrap(~hap) +
#   geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
#   scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
#   theme_classic() + 
#   labs(x = '', y = '', fill = 'Allele') + 
#   theme(axis.text = element_blank(),
#         legend.position = 'none') + 
#   geom_vline(aes(xintercept = 144))


# ====== Sequence Analysis =====


syn_vs_nonsyn <- function(CDS_dt, rnk, nuc, codon_position){
  # Function to classify polymorphisms as synonymous or nonsynonymous
  
  CDS <- copy(CDS_dt)
  if(codon_position == 1){
    refaa <- translate(CDS[seq_rank %in% rnk:(rnk+2), refseq])
    qryaa <- translate(c(nuc, CDS[seq_rank %in% (rnk+1):(rnk+2), refseq]))
  }
  if(codon_position == 2){
    refaa <- translate(CDS[seq_rank %in% (rnk-1):(rnk+1), refseq])
    qryaa <- translate( c( CDS[seq_rank %in% (rnk-1), refseq], nuc, CDS[seq_rank %in% (rnk+1), refseq]) )
  }
  if(codon_position == 3){
    refaa <- translate(CDS[seq_rank %in% (rnk-2):(rnk), refseq])
    qryaa <- translate( c( CDS[seq_rank %in% (rnk-2):(rnk-1), refseq], nuc ) )
  }
  if(refaa == qryaa){
    return('S')
  } else if(refaa != qryaa){
    return('N')
  }
}

# ----- coding sequence coordinates for TPPD1 -----

# coordinates manually subset from annotation files
coords01 <- fread('~/workspace/heterodichogamy/hindsii/Jcali_primary_TPPD-1.gff', sep = '\t')
coords01 <- coords01[, c(4,5,7)]
setnames(coords01, c("start", "end", 'strand'))
setkey(coords01, "start")

seq01 <- read.fasta('~/workspace/heterodichogamy/01_G_locus_structure/TPPD/CDS/Jcali_G_corrected.fasta')
seq01 <- as.character(seq01[[1]])
#length(seq01) 

# combine positions and CDS sequence
if(coords01[1, strand == '-']){
  cds01 <- data.table(coords01[, seq(start,end), by = start][, .(pos = rev(V1))], refseq=seq01)
} else{
  cds01 <- data.table(coords01[, seq(start,end), by = start][, .(pos = V1)], refseq=seq01)
}


# assign codon position
cds01[, seq_rank := seq(1, .N)]
cds01[seq_rank %% 3 == 1, codon_position := 1]
cds01[seq_rank %% 3 == 2, codon_position := 2]
cds01[seq_rank %% 3 == 0, codon_position := 3]

# get allele counts by haplotype at each site
cnts01 <- haps[pos >= coords01[, min(start)] & pos <= coords01[, max(end)], .(count = sum(allele)), by = .(pos, hap, allele0, allele1)]
cnts01 <- dcast(cnts01, pos + allele0 + allele1 ~ hap, value.var = 'count')

# *conditional reverse complement*
if(coords01[1, strand == '-']){
  cnts01[, c('allele0', 'allele1') := lapply(.SD, comp), .SDcols = c('allele0', 'allele1')]
}

final01 <- merge(cnts01, cds01, by = 'pos')
setkey(final01, seq_rank)

# write out positions of fixed SNPs for other analyses
#fwrite(final01[g -G >= 68, .(pos)], file = '~/workspace/heterodichogamy/hindsii/fixed_TPPD1-SNPs.txt', quote = F, col.names = F, row.names = F)

# classify
vals <- vector()
for(i in 1:nrow(final01)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = cds01, rnk = final01$seq_rank[i], nuc = final01$allele1[i], codon_position = final01$codon_position[i])
}

final01[, type := vals]

# classify polymorphisms
# how many of each morph/ genotype/ hap?
haps[SNPrank == 1,
     table(hap)]

haps[SNPrank == 1 & hap.id == 1,
     table(genotype)]



# ----- classify polymorphisms -----
final01[(G >= 22 & g <= 1) | (G <= 1 & g >= 68), fixed_poly := 'fixed']
final01[G > 1 & G < 22 & g > 1 & g < 68, fixed_poly := 'polymorphic_both']
final01[(G > 1 & G < 22) & (g <= 1 | g >= 68), fixed_poly := 'polymorphic_G']
final01[(G <= 1 | G >= 22) & (g > 1 & g < 68), fixed_poly := 'polymorphic_g']


final01[, fixed_poly_MK := fixed_poly]
final01[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_g']

MK_tbl <- table(final01$fixed_poly_MK, final01$type)
MK_tbl

# ========= NDR1 ==========
coords02 <- fread('~/workspace/heterodichogamy/hindsii/Jcali_primary_NDR1.gff', sep = '\t')
coords02 <- coords02[, c(4,5,7)]
setnames(coords02, c("start", "end", 'strand'))
setkey(coords02, "start")

cnts02 <- haps[pos >= coords02[, min(start)] & pos <= coords02[, max(end)], .(count = sum(allele)), by = .(pos, hap, allele0, allele1)]

cnts02 <- dcast(cnts02, pos + allele0 + allele1 ~ hap, value.var = 'count')

cnts02
# no fixed differences in NDR1/HIN1-like coding sequence

# Is this polymoprhism nonsynonymous or synonymous
coords02 <- fread('~/workspace/heterodichogamy/hindsii/Jcali_primary_NDR1.gff', sep = '\t')
coords02 <- coords02[, c(4,5,7)]
setnames(coords02, c("start", "end", 'strand'))
setkey(coords02, "start")

seq02 <- read.fasta('~/workspace/heterodichogamy/01_G_locus_structure/NDR1/Jcali_primary_NDR1_CDS.fasta')
seq02 <- as.character(seq02[[1]])
#length(seq01) 

# combine positions and CDS sequence
if(coords02[1, strand == '-']){
  cds02 <- data.table(coords02[, seq(start,end), by = start][, .(pos = rev(V1))], refseq=seq02)
} else{
  cds02 <- data.table(coords02[, seq(start,end), by = start][, .(pos = V1)], refseq=seq02)
}


# assign codon position
cds02[, seq_rank := seq(1, .N)]
cds02[seq_rank %% 3 == 1, codon_position := 1]
cds02[seq_rank %% 3 == 2, codon_position := 2]
cds02[seq_rank %% 3 == 0, codon_position := 3]


# get allele counts by haplotype at each site
cnts02 <- haps[pos >= coords02[, min(start)] & pos <= coords02[, max(end)], .(count = sum(allele)), by = .(pos, hap, allele0, allele1)]
cnts02 <- dcast(cnts02, pos + allele0 + allele1 ~ hap, value.var = 'count')


# *conditional reverse complement* (not needed in this case)
if(coords02[1, strand == '-']){
  cnts02[, c('allele0', 'allele1') := lapply(.SD, comp), .SDcols = c('allele0', 'allele1')]
}

final02 <- merge(cnts02, cds02, by = 'pos')
setkey(final02, seq_rank)

final02

# classify
vals <- vector()
for(i in 1:nrow(final02)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = cds02, rnk = final02$seq_rank[i], nuc = final02$allele1[i], codon_position = final02$codon_position[i])
}

final02[, type := vals]


