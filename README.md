# TSSanalysis
# Extension of 5'UTR region of each gene from CAGE data and public gene model
# Input sorted SAM file describing only 1st read of paired-end reads prepared by samtools
# 
# ref: Tokizawa M, Kusunoki K, Koyama H, Kurotani A, Sakurai T, Suzuki Y, Kurata T, Yamamoto YY. "Identification of Arabidopsis genic and non-genic promoters by pair-end sequencing of TSS tags". Plant J, 2017.
#
# argument 1: input SAM file name
# argument 2: input GFF file name
# argument 3: input feature name such as "gene", "protein", "mRNA", "pseudogenic_transcript", "transposable_element_gene"
# argument 4: Length [bp]
# argument 5: loop count of extension of 5'end of each gene
# 
# Output file:
# 1st-9th columns: same to gff format ("Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute")
# 10th column: Position of extended 5'end UTR (original positions are shown in case of no extension)
