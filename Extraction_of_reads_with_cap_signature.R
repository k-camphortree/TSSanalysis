# Extraction of paired-end reads with cap-signiture (Only reads with G at the 5' end of the first (R1) reads) for identification of TSS-tags
# Input directory including TSS-seq paired-end fastq files

# ref: Tokizawa M, Kusunoki K, Koyama H, Kurotani A, Sakurai T, Suzuki Y, Kurata T, Yamamoto YY. "Identification of Arabidopsis genic and non-genic promoters by pair-end sequencing of TSS tags". Plant J, 2017.

# Output files:
# fastq files would be divided into two files: 
#  WithCap_(fastq_file_name)
#  WithoutCap_(fastq_file_name)
# 

#Argument
in_dir = commandArgs(trailingOnly=TRUE)[1]
out_dir = commandArgs(trailingOnly=TRUE)[2]

#Extraction of paired-end reads with cap-signiture
if(substr(in_dir, nchar(in_dir), nchar(in_dir))!="/"){
  in_dir <- paste(in_dir, "/", sep="")
}
if(substr(in_dir, nchar(out_dir), nchar(out_dir))!="/"){
  out_dir <- paste(out_dir, "/", sep="")
}
library(ShortRead)
fastqR1<-list.files(in_dir, pattern="R1")
fastqR2<-list.files(in_dir, pattern="R2")
for (j in 1:length(fastqR1)){
  fastaR1 <- readFastq(paste(in_dir,fastqR1[j],sep=""))
  objR1<-as.matrix(subseq(sread(fastaR1),start=1,end=1))[1:length(fastaR1)]=="G"
  fastaR1withCap <- fastaR1[objR1]
  fastaR1withoutCap <- fastaR1[!objR1]
  writeFastq(fastaR1withCap, file=paste(out_dir, "onlyG_", fastqR1[j], sep=""))
  writeFastq(fastaR1withoutCap, file=paste(out_dir, "exceptG_",fastqR1[j],sep=""))  
  fastaR2 <- readFastq(paste(in_dir, fastqR2[j], sep=""))
  objR2withCap <- is.element(gsub("2:N:0","",id(fastaR2)),gsub("1:N:0","",id(fastaR1withCap)))
  objR2withoutCap <- is.element(gsub("2:N:0","",id(fastaR2)),gsub("1:N:0","",id(fastaR1withoutCap)))
  fastaR2withCap <- fastaR2[obj2RwithCap]
  fastaR2withoutCap <- fastaR2[objR2withoutCap]
  writeFastq(fastaR2withCap, file=paste(out_dir, "onlyG_",fastq2[j],sep=""))
  writeFastq(fastaR2withoutCap, file=paste(out_dir, "exceptG_",fastq2[j],sep=""))
}
