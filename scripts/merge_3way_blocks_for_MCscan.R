library(data.table)


# ----- 5 way ----- 
AB <- fread('~/workspace/heterodichogamy/HC_locus_structure/MCscan/Lakota_v1.Csin.i1.blocks', 
            header=F, sep = '\t',
            col.names = c("A", "B"))

AC <- fread('~/workspace/heterodichogamy/HC_locus_structure/MCscan/Lakota_v1.Pawnee.i1.blocks', 
            header=F, sep = '\t',
            col.names = c("A", "C"))

AD <- fread("~/workspace/heterodichogamy/HC_locus_structure/MCscan/Lakota_v1.Ccat.i1.blocks",
            header = F, sep = '\t', 
            col.names = c("A", "D"))

AE <- fread("~/workspace/heterodichogamy/HC_locus_structure/MCscan/Lakota_v1.Jregia.i1.blocks",
            header = F, sep = '\t', 
            col.names = c("A", "E"))


ABC.blocks <- merge(AB, AC, all= T, by = 'A')
ABCD.blocks <- merge(ABC.blocks, AD, all=T, by = 'A')
ABCDE.blocks <- merge(ABCD.blocks, AE, all=T, by = 'A')


ABCDE.blocks[is.na(A), A := '.']
ABCDE.blocks[is.na(B), B := '.']
ABCDE.blocks[is.na(C), C := '.']
ABCDE.blocks[is.na(D), D := '.']
ABCDE.blocks[is.na(E), E := '.']


fwrite(ABCDE.blocks[, .(A, B, C, D, E)], 
       file = '~/workspace/heterodichogamy/HC_locus_structure/MCscan/5way.blocks.manual.edit',
       col.names = F,
       row.names = F, 
       quote = F,
       sep = '\t')





# ----- 3 way -----

AB <- fread('~/workspace/heterodichogamy/HC_locus_structure/MCscan/Lakota_v1.Pawnee.i1.blocks', 
             header=F, sep = '\t',
             col.names = c("A", "B"))

AC <- fread('~/workspace/heterodichogamy/HC_locus_structure/MCscan/Lakota_v1.Jregia.i1.blocks', 
            header=F, sep = '\t',
            col.names = c("A", "C"))

ABC.blocks <- merge(AB, AC, all= T, by = 'A')

ABC.blocks[is.na(A), A := '.']
ABC.blocks[is.na(B), B := '.']
ABC.blocks[is.na(C), C := '.']

fwrite(ABC.blocks[, .(A, B, C)], 
       file = '~/workspace/heterodichogamy/HC_locus_structure/MCscan/3way.blocks.manual.edit',
       col.names = F,
       row.names = F, 
       quote = F,
       sep = '\t')



