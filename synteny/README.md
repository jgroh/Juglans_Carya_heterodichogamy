# Commands to perform synteny analysis in MCscan, adapted from https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)

Do the following for A->B and B->C:

1)	Generate bed for each genome
	python -m jcvi.formats.gff bed --type=mRNA --key={} --primary_only {gff} -o {genome}.bed
	For regia, Csin, Ccat, use --key=Parent
	For Lakota_v1 and Pawnee, use --key=ID

2)	Generate CDS files
	python -m jcvi.formats.gff load {genome}.gff {genome}.fasta -o {genome}.cds --id_attribute=ID
	Again use id_attribute=Parent for regia, ID for pecan genomes

3)	Synteny search 
	python -m jcvi.compara.catalog ortholog {A} {B} --no_strip_names --cscore=0.99

4)	“Compute layout for gene-level matchings”
	python -m jcvi.compara.synteny mcscan A.bed A.B.lifted.anchors --iter=1 -o A.B.i1.blocks

5)	manual merge of blocks file into 3-way or 5-way blocks file
	use R script merge_3way_blocks_for_MCscan.R

6)	manually edit Layout file

7)	merge bed files
	python -m jcvi.formats.bed merge A.bed B.bed C.bed -o ABC.bed5

8)	Create plot
	python -m jcvi.graphics.synteny {blocks} ABC.bed {layout}

