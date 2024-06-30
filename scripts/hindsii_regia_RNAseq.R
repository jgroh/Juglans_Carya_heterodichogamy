library(data.table)
library(tximport)
library(forcats)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(tidyverse)
library(apeglm)


hindsii_pg <- paste0("JHIN_", c("001", "013", "023"))
hindsii_pa <- paste0("JHIN_", c("002", "006", "014"))

regia_pg <- paste0("JREG_", c("SHARK", "PLACN", "CHICO", "ERLST"))
regia_pa <- paste0("JREG_", c("CHAND", "PI159", "KLPNG"))

meta <- data.table(ID = c(hindsii_pg, hindsii_pa, regia_pg, regia_pa), 
                   species = c(rep("hindsii", 6), rep("regia", 7)),
                   type = c(rep("PG", 3), rep("PA", 3), rep("PG", 4), rep("PA", 3)))

meta[, type := factor(type, levels = c("PA", "PG"))]
meta.regia <- meta[species == 'regia']
meta.hindsii <- meta[species == 'hindsii']

# # ----- hindsii raw read depth from STAR across G-locus (not standardized) -----
# hindsii_dp <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/gene_expression/hindsii_regia/coverage/hindsii/", full.names = T),
#        function(x){
#         z <- fread(x, col.names = c("chr", "pos", "dp"))
#         z[, ID := gsub(".txt.gz", "", basename(x))]
#        }
# ))
# 
# hindsii_dp[ID %in% hindsii_pg, type := 'PG']
# hindsii_dp[ID %in% hindsii_pa, type := 'PA']
# 
# ggplot(hindsii_dp, aes(x = pos, y = dp)) +
#   facet_wrap(~ID) +
#   geom_area()
# 
# ggplot(hindsii_dp[pos > 30710000 & pos < 30714000], aes(x = pos, y = dp)) +
#   facet_wrap(~ID) +
#   geom_area()  +
#   geom_vline(aes(xintercept = 30710259), color = 'red') +
#   geom_vline(aes(xintercept = 30713073), color = 'red')


# ----- regia raw read depth across G-locus (not standardized) -----
# regia_dp <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/gene_expression/hindsii_regia/coverage/regia/", full.names = T),
#                                function(x){
#                                  z <- fread(x, col.names = c("chr", "pos", "dp"))
#                                  z[, ID := gsub(".txt.gz", "", basename(x))]
#                                }
# ))
# 
# regia_dp[ID %in% regia_pg, type := 'PG']
# regia_dp[ID %in% regia_pa, type := 'PA']
# 
# ggplot(regia_dp[pos > 31884268], aes(x = pos, y = dp)) +
#   facet_wrap(~ID) +
#   geom_area()  


# ------ import transcript quantification files ------
#Create a named vector pointing to the quantification files

setwd("~/workspace/heterodichogamy/gene_expression/hindsii_regia/")
quants_files <- list.files("quants/", recursive = T, pattern = ".*quant\\.sf$", full.names = T)
all_names <- gsub("_quant/quant.sf", "", gsub("quants//", "", quants_files))
names(quants_files) <- all_names

regia_files <- quants_files[grepl("JREG", quants_files)]
hindsii_files <- quants_files[grepl("JHIN", quants_files)]

regia_names <- gsub("_quant/quant.sf", "", gsub("quants//", "", regia_files))
names(regia_files) <- regia_names
hindsii_names <- gsub("_quant/quant.sf", "", gsub("quants//", "", hindsii_files))
names(hindsii_files) <- hindsii_names

tx2gene <- fread("~/workspace/heterodichogamy/gene_expression/hindsii_regia/txt2gene.txt", header = F, col.names = c("tx", "gene"))

# The tximport package has a single function for importing transcript-level estimates. 
# The type argument is used to specify what software was used for estimation. 
# A simple list with matrices, "abundance", "counts", and "length", is returned, 
# where the transcript level information is summarized to the gene-level.

txi.all <- tximport(quants_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = 'scaledTPM')
head(txi.all$counts)

txi.regia <- tximport(regia_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = 'scaledTPM')
txi.hindsii <- tximport(hindsii_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = 'scaledTPM')

# *** critical to make rowns of metadata match columns of count matrix
# From DESeqDataSetFromTximport help page:
# Rows of colData correspond to columns of countData
all_order_indices <- match(names(as.data.table(txi.all$counts)), meta$ID)
meta <- meta[all_order_indices]

regia_order_indices <- match(names(as.data.table(txi.regia$counts)), meta.regia$ID)
meta.regia <- meta.regia[regia_order_indices]

hindsii_order_indices <- match(names(as.data.table(txi.hindsii$counts)), meta.hindsii$ID)
meta.hindsii <- meta.hindsii[hindsii_order_indices]


# ------ DESeq2 ------
ddsTxi.all <- DESeqDataSetFromTximport(txi.all,
                                         colData = meta,
                                         design = ~ species + type)


ddsTxi.regia <- DESeqDataSetFromTximport(txi.regia,
                                   colData = meta.regia,
                                   design = ~ type)

ddsTxi.hindsii <- DESeqDataSetFromTximport(txi.hindsii,
                                         colData = meta.hindsii,
                                         design = ~ type)


# ----- all -----
dds.all <- DESeq(ddsTxi.all)
res.all <- results(dds.all, name="type_PG_vs_PA")

topn <- 10
all_top <- rownames(head(res.all[order(res.all$pvalue),], topn))
all_top

all_TPPD1_cnts <- plotCounts(dds.all, gene='gene-LOC108984907', intgroup="type", returnData = T) 
smpls <- rownames(all_TPPD1_cnts)

setDT(all_TPPD1_cnts)
all_TPPD1_cnts[, mean(count), by = type]

all_TPPD1_cnts[, ID := smpls]


# inspect results for TPPD-1
res.all[which(rownames(res.all) == 'gene-LOC108984907'),]

# ---- regia -----
dds.regia <- DESeq(ddsTxi.regia)

res.regia <- results(dds.regia, name="type_PG_vs_PA")
res.regia.shrunk <- lfcShrink(dds.regia, coef="type_PG_vs_PA", res=res.regia)

# inspect results for TPPD-1
res.regia[which(rownames(res.regia) == 'gene-LOC108984907'),]
res.regia.shrunk[which(rownames(res.regia.shrunk) == 'gene-LOC108984907'),]

# regia plot gene counts
regia_TPPD1_cnts <- plotCounts(dds.regia, gene='gene-LOC108984907', intgroup="type", returnData = T) 

# ----- hindsii -----
dds.hindsii <- DESeq(ddsTxi.hindsii)

res.hindsii <- results(dds.hindsii, name="type_PG_vs_PA")
res.hindsii.shrunk <- lfcShrink(dds.hindsii, coef="type_PG_vs_PA", res=res.regia)

# results for TPPD-1
res.hindsii[which(rownames(res.hindsii) == 'gene-LOC108984907'),]
res.hindsii.shrunk[which(rownames(res.hindsii.shrunk) == 'gene-LOC108984907'),]

# plot gene counts
hindsii_TPPD1_cnts <- plotCounts(dds.hindsii, gene='gene-LOC108984907', intgroup="type", returnData = T) 


# ----- compare X most significant DEGs between species ------

topn <- 100
regia_top <- rownames(head(res.regia[order(res.regia$pvalue),], topn))
hindsii_top <- rownames(head(res.hindsii[order(res.hindsii$pvalue),], topn))

intersect(regia_top, hindsii_top)


# ----- interesting candidate genes----- 
# 1. gene-LOC109002076
# annotation: SWEET4-like 

# paper on homolog 
# https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-020-2266-0

plotCounts(dds.regia, gene='gene-LOC109002076', intgroup="type")
plotCounts(dds.hindsii, gene='gene-LOC109002076', intgroup="type")

#2. gene-LOC108983049
# uncharacterized
plotCounts(dds.regia, gene='gene-LOC108983049', intgroup="type")
plotCounts(dds.hindsii, gene='gene-LOC108983049', intgroup="type")


#3. gene-LOC108983651
# annotation: agamous-like MADS-box protein TM6
plotCounts(dds.regia, gene='gene-LOC108983651', intgroup="type")
plotCounts(dds.hindsii, gene='gene-LOC108983651', intgroup="type")



# ----- regia allele depth -----
adr <- fread("~/workspace/heterodichogamy/gene_expression/hindsii_regia/regia_allele_depths_TPPD-1.txt", header = F, sep = " ")
regia_fixed <- fread("~/workspace/heterodichogamy/gene_expression/Jreg_Dang_etal_2016/Gloc_fixed_SNPs_Chandler.txt", col.names = 'pos')

# reformat operations
setnames(adr, paste0("V", 1:3), c("pos", "ref", "alt"),)
adr <- melt(adr, id.vars = c("pos", "ref", "alt"))

st <- 'alignment_files/JregCha/'
adr[, value := gsub("filtered.bam=", "", gsub(st, "", value))]

adr[, c("species", "ID", "AD") := tstrsplit(value, split = "_")]
adr[, value := NULL][, variable := NULL][]
adr[, c("ref_cnt", "alt_cnt") := tstrsplit(AD, split = ",")]

adr[, ref_cnt := as.numeric(ref_cnt)]
adr[, alt_cnt := as.numeric(alt_cnt)]

adr[, AD := NULL]
adr[, species := NULL]
adr[, pos := gsub("NC_049911.1:", "", pos)]

adr <- adr[pos %in% regia_fixed$pos]

adr

adr[, prop.G := alt_cnt/(ref_cnt+alt_cnt), by = .(ID, pos)]
prop.G.regia <- adr[, .(G = weighted.mean(prop.G, w = alt_cnt + ref_cnt, na.rm = T)), by = ID]

prop.G.regia <- prop.G.regia[, .(ID = paste0("JREG_", ID), G, g = 1-G, species = 'regia')]



# ----- hindsii allele depth -----
adh <- fread("~/workspace/heterodichogamy/gene_expression/hindsii_regia/hindsii_allele_depths_TPPD-1_Jcali_alt.txt", header = F, sep = " ")
hindsii_fixed <- fread("~/workspace/heterodichogamy/gene_expression/hindsii_regia/calls/hindsii/Jcali_alt_TPPD1_fixed_diffs.txt", select = 2, col.names = 'pos')

# reformat operations
setnames(adh, paste0("V", 1:3), c("pos", "ref", "alt"),)
adh <- data.table::melt(adh, id.vars = c("pos", "ref", "alt"))

sth <- 'alignment_files/Jcali_alt/'
adh[, value := gsub("filtered.bam=", "", gsub(sth, "", value))]

adh[, c("species", "ID", "AD") := tstrsplit(value, split = "_")]
adh[, value := NULL][, variable := NULL][]

adh[, c("ref_cnt", "alt_cnt", "BLANK") := tstrsplit(AD, split = ",")]
adh[, "BLANK" := NULL]

adh[, ref_cnt := as.numeric(ref_cnt)]
adh[, alt_cnt := as.numeric(alt_cnt)]

adh[, AD := NULL]
adh[, species := NULL]
adh[, pos := gsub("JAKSXL010000006.1:", "", pos)]

adh <- adh[pos %in% hindsii_fixed$pos]

# PG samples
adh[ID == '001']
adh[ID == '013']
adh[ID == '023']

# PA samples
adh[ID == '006'] 
adh[ID == '014']
adh[ID == '002']


adh[, prop.G := alt_cnt/(ref_cnt+alt_cnt), by = .(ID, pos)]
prop.G.hindsii <- adh[, .(G = weighted.mean(prop.G, w = ref_cnt + alt_cnt, na.rm = T)), by = ID]

# software automatically gives 0.5 for log transform, even if no transcription detected, so manually adjust
prop.G.hindsii[ID == '006', G := 0] 
prop.G.hindsii[ID == '014', G := 0]
prop.G.hindsii
prop.G.hindsii[, g := 1-G]

prop.G.hindsii <- prop.G.hindsii[, .(ID = paste0("JHIN_", ID), G, g, species = 'hindsii')]
prop.G.hindsii

# ---------- Main Plot -------
all_TPPD1_cnts

# close inspection of data indicated count value for this sample is error rather than real signal
all_TPPD1_cnts[ID == "JHIN_006", count := 0]

all_TPPD1_cnts[grepl("JHIN", ID), species := 'hindsii']
all_TPPD1_cnts[grepl("JREG", ID), species := 'regia']

prop.G <- rbind(prop.G.hindsii, prop.G.regia)
prop.G
all_TPPD1_cnts

plot_data <- merge(prop.G, all_TPPD1_cnts)

plot_data[, G := count*G]
plot_data[, g := count*g]

plot_data <- melt(plot_data[, .(ID, type, species, G, g)], id.vars = c("ID", "type", "species"), value.name = 'count', variable.name = 'hap')
plot_data

ggplot(plot_data, aes(x = reorder(ID, count), y = count, fill = hap)) + 
  geom_bar(stat = 'identity') + 
  facet_wrap(~species, scales = 'free_x') + 
  scale_fill_manual(values = c("#413a6e", "#cfb582"))  + 
  theme_classic() + 
  labs(x = "", y = "Normalized count", fill = '') + 
  theme(aspect.ratio = 1,
        strip.text = element_text(face = 'italic'),
        axis.text.x = element_blank(),
        
        strip.text.x = element_blank(),
        axis.ticks.length.x = unit(.0, 'cm'))

# ------ NDR1 -------
all_NDR1_cnts <- plotCounts(dds.all, gene='gene-LOC108984910', intgroup="type", returnData = T) 
smpls <- rownames(all_TPPD1_cnts)

all_NDR1_cnts[, mean(count), by = type]
setDT(all_TPPD1_cnts)
all_TPPD1_cnts[, ID := smpls]


# ------ vizualize TPM  -----
cnts0 <- data.table(txi.all$counts)
cnts0[, gene := rownames(txi.all$counts)]

cnts <- melt(cnts0, id.vars = 'gene', variable.name = 'Run', value.name = 'TPM')

# merge with metadata
cnts[Run %in% c(regia_pg, hindsii_pg), type := 'PG']
cnts[Run %in% c(regia_pa, hindsii_pa), type := 'PA']

# subset to TPPD-1
TPPD1_cnts <- cnts[gene == 'gene-LOC108984907']


TPPD1_cnts[, c("species", "ID") := tstrsplit(Run, "_")]
TPPD1_cnts[species == 'JHIN', species := 'J. hindsii']
TPPD1_cnts[species == 'JREG', species := 'J. regia']

TPPD1_cnts$ID <- fct_reorder(TPPD1_cnts$ID, TPPD1_cnts$TPM)

ggplot(TPPD1_cnts, aes(x = ID, y = TPM, fill = type)) + 
  geom_bar(stat = 'identity') + 
  facet_wrap(~species, scales = 'free_x') + 
  scale_fill_manual(values = c("#cfb582","#413a6e"))  + 
  theme_classic() + 
  labs(x = "") + 
  theme(aspect.ratio = 1,
        strip.text = element_text(face = 'italic'),
        axis.text.x = element_blank()) 





