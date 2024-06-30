library(data.table)

# This script calculates Dxy across a specified reference chromosome from anchorwave alignments
# whole genome alignments were done with anchorwave proali -R1 -Q1. 

# ===== Read Inputs =====
args <- commandArgs(trailingOnly = TRUE)
# Arguments:
  # 1. MAF file for chromosome 11 (now whole genome)
  # 2. ref gnom id
  # 3. qry gnom id
  # 4. ref chrom id (as written in maf)
  # 5. window size

# read alignment. Use sep value that's not in file to force fread to behave like readLines
mafLines <- fread(args[1], sep = "?", skip = 1,
                  header = F)
ref <- args[2]
qry <- args[3]
chrom <- args[4]
window_size <- as.numeric(args[5])


# ===== format maf alignment into data.table =====

# subset to alignment blocks
mafLines <- mafLines[grep("^s", V1)]

# set column names
# https://genome.ucsc.edu/FAQ/FAQformat.html#format5
mafDT <- mafLines[, tstrsplit(V1, "\\s+")]
setnames(mafDT, c("s", "chr", "start0", "length", "strand", "srcLength", "seq"))

# subset to focal chromosome. Find which rows have the chromosome ID, and also take the next row which has the aligned chromosome
# assumes the maf format
#mafDT[1:10, 1:5]
mafDT[seq(2, nrow(mafDT), by = 2), chr := paste0(chr, '_qry')]

indices <- which(mafDT$chr == chrom)

# in case chromosome names are the same between alignments, subset odd numbered ones as these come first
# in the maf alignment blocks. (not needed anymore, renamed query chroms in step above)
# indices <- indices[indices %% 2 == 1]

# then the following line is the query chromosome in the alignmetn
indices <- c(indices, indices + 1)
indices <- indices[order(indices)]
mafDT <- mafDT[indices]
    
# add grouping variable for alignment id
mafDT[, alignmentID := rep(1:(.N / 2), each = 2)]
# mafDT[, c(1:6,8)]

# subset to focal columns
mafDT <- mafDT[, .(alignmentID, chr, start0 = as.numeric(start0), length = as.numeric(length), seq)]
# mafDT[, c(1:4)]

# function that takes alignment block and returns data table with aligned seqs in long format. First row in alignment block is ref, 2nd is qry.
tr <- function(DT, by){
  outDT <- data.table(alignmentID = by,
                      refChr = DT[1, chr], 
                      qryChr = DT[2, chr], 
                      refStart = DT[1, start0], 
                      qryStart = DT[2, start0],
                      refLength = DT[1, length], 
                      qryLength = DT[2, length],
                      refSeq = unlist(strsplit(DT[1, seq], "")),
                      qrySeq = unlist(strsplit(DT[2, seq], ""))
                      )
  
  # For calculating Dxy against a reference, 
  # sites in the alignment that are absent in the reference are not informative, so we disregard these. 
  outDT <- outDT[refSeq != "-"]
  
  # set coordinates of reference. 
  outDT[, refPos := seq_len(.N)]
  outDT[, refPos := refPos + refStart]
  outDT
}


seqDT <- mafDT[, tr(DT = .SD, by = .BY), by = alignmentID]


# ====== Tally Substitutions ======

# assign logical to substitutions
seqDT[, sub := ifelse(refSeq != qrySeq, T, F)]

# gaps in the query sequence are treated as NA 
# (we account for the total number of non-gapped bps in the calc of Dxy)
seqDT[qrySeq == "-", sub := NA]


# ===== Assign each bp to its containing window =====
round_up <- function(x) {
  result <- ceiling(x / window_size) * window_size
  return(result)
}

round_down <- function(x) {
  result <- floor(x / window_size) * window_size
  return(result)
}


seqDT[, window := cut(refPos, 
                      breaks = seq(from = round_down(refStart[1]), to = round_up(refStart[1] + refLength[1]), by = window_size) , 
                      include.lowest=T, 
                      labels = seq(from = round_down(refStart[1]), 
                                   to = round_up(refStart[1] + refLength[1]) - window_size, by = window_size) + 0.5*window_size), 
      by = alignmentID]


# ===== compute Dxy =====
Dxy <- seqDT[, .(Dxy_num = sum(sub, na.rm=T), 
                 Dxy_denom = sum(!is.na(sub))), by = .(alignmentID, refChr, qryChr, window)]
Dxy[, Dxy := Dxy_num/Dxy_denom]
Dxy[Dxy_num == 0 & Dxy_denom == 0, Dxy := NA]

Dxy[, refGnom := ref]
Dxy[, qryGnom := qry]

# ===== Output =====
fwrite(Dxy, file = "", quote = F, sep = "\t", row.names = F, na = "NA")

