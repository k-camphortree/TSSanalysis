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
# 10th column: Position of extended 5'end UTR (Result "Inf" means no extension)
# 

#Argument
in_f = commandArgs(trailingOnly=TRUE)[1]
gffFile = commandArgs(trailingOnly=TRUE)[2]
feature = commandArgs(trailingOnly=TRUE)[3]
ReadLength = as.numeric(commandArgs(trailingOnly=TRUE)[4])
loop = as.numeric(commandArgs(trailingOnly=TRUE)[5])

#Calculation of range covered by each reads from SAM file
SamR1 <- read.table(in_f, header=F, sep="\t", skip=max(grep("@",readLines(in_f))), fill=TRUE, quote="", comment.char="")
SamR1Range <- cbind(SamR1, rep(0, nrow(SamR1)), rep(0, nrow(SamR1)))
colnames(SamR1Range)[c(3, ncol(SamR1)+1, ncol(SamR1)+2)] <- c("Chr", "TSS", "Anchor")
SamR1Range[SamR1[,9]>0, ncol(SamR1)+1] <- SamR1[SamR1[,9]>0, 4]+1
SamR1Range[SamR1[,9]>0, ncol(SamR1)+2] <- SamR1[SamR1[,9]>0, 4]+SamR1[SamR1[,9]>0,9]-1
SamR1Range[SamR1[,9]<0, ncol(SamR1)+1] <- SamR1[SamR1[,9]<0, 8]-SamR1[SamR1[,9]<0,9]-2
SamR1Range[SamR1[,9]<0, ncol(SamR1)+2] <- SamR1[SamR1[,9]<0, 8]

#BUndle overlapped TSS paired tags
BundledTagRange <- matrix(c(), ncol=5, nrow=1)
BundledTagRange <- matrix(,ncol=5)
colnames(BundledTagRange) <- c("Chr", "TSS", "Anchor", "Count", "Direction")
for (j in names(table(SamR1Range$Chr))){
  SamR1RangeChr<-subset(SamR1Range,Chr==j)
  Reach_plus <- tapply(SamR1RangeChr[SamR1RangeChr[,9]>0,]$Anchor, SamR1RangeChr[SamR1RangeChr[,9]>0,]$TSS, max)
  Count_plus <- tapply(SamR1RangeChr[SamR1RangeChr[,9]>0,]$Anchor, SamR1RangeChr[SamR1RangeChr[,9]>0,]$TSS, function(x)(sum(hist(x,plot=F)$counts)))
  Comb_plus <- data.frame(cbind(j, names(Reach_plus), as.numeric(Reach_plus), as.numeric(Count_plus), "Plus"))
  colnames(Comb_plus)<-c("Chr", "TSS", "Anchor", "Count", "Direction")
  Reach_minus <- tapply(SamR1RangeChr[SamR1RangeChr[,9]<0,]$Anchor,SamR1RangeChr[SamR1RangeChr[,9]<0,]$TSS, min)
  Count_minus <- tapply(SamR1RangeChr[SamR1RangeChr[,9]<0,]$Anchor,SamR1RangeChr[SamR1RangeChr[,9]<0,]$TSS, function(x)(sum(hist(x,plot=F)$counts)))
  Comb_minus <- data.frame(cbind(j, names(Reach_minus), as.numeric(Reach_minus), as.numeric(Count_minus), "Minus"))
  colnames(Comb_minus) <- c("Chr", "TSS", "Anchor", "Count", "Direction")
  Comb_both <- rbind(Comb_plus, Comb_minus)
  BundledTagRange <- rbind(BundledTagRange, Comb_both)
}
BundledTagRange <- BundledTagRange[-1,]
BundledTagRangeLength <- cbind(BundledTagRange, (as.numeric(BundledTagRange$Anchor)-as.numeric(BundledTagRange$TSS)+1))
colnames(BundledTagRangeLength)[6] <- "Length"
BundledTagRangeLength$TSS <- as.numeric(as.character(BundledTagRangeLength$TSS))
BundledTagRangeLength$Anchor <- as.numeric(as.character(BundledTagRangeLength$Anchor))
BundledTagRangeLength$Count <- as.numeric(as.character(BundledTagRangeLength$Count))

#Update Gene_model
gff <- read.table(gffFile, header=F, sep="\t", quote="")
colnames(gff) <- c("Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute")
SelectedGff <- subset(gff, Feature==feature)
UpdatedGff <- data.frame()
for (j in names(table(SamR1Range$Chr))){
  BundledTagRangeLengthChr <- subset(BundledTagRangeLength, Chr==j)
  SelectedGffChr <- subset(SelectedGff, Chr==j)
  BundledTagRangeLengthChrPlus <- subset(BundledTagRangeLengthChr, Direction=="Plus")
  SelectedGffChrPlus <- subset(SelectedGffChr, Strand=="+")
  FarthestAnchorTagPlus <- c()
  for (n in 1:nrow(SelectedGffChrPlus)){
    AnchorTagPlus <- subset(BundledTagRangeLengthChrPlus, Anchor>=(SelectedGffChrPlus[n,4]+ReadLength-1) & Anchor<=SelectedGffChrPlus[n,5])
    FarthestAnchorTagPlus <- c(FarthestAnchorTagPlus,min(AnchorTagPlus$TSS))
  }
  SelectedGffChrPlusAnchor <- cbind(SelectedGffChrPlus,FarthestAnchorTagPlus)
  colnames(SelectedGffChrPlusAnchor)[10] <- "ExtendedTSS"
  BundledTagRangeLengthChrMinus <- subset(BundledTagRangeLengthChr, Direction=="Minus")
  SelectedGffChrMinus <- subset(SelectedGffChr, Strand=="-")
  FarthestAnchorTagMinus <- c()
  for (n in 1:nrow(SelectedGffChrMinus)){
    AnchorTagMinus <- subset(BundledTagRangeLengthChrMinus, Anchor>=SelectedGffChrMinus[n,4] & Anchor<=(SelectedGffChrMinus[n,5]-ReadLength+1))
    FarthestAnchorTagMinus <- c(FarthestAnchorTagMinus,max(AnchorTagMinus$TSS))
  }
  SelectedGffChrMinusAnchor <- cbind(SelectedGffChrMinus, FarthestAnchorTagMinus)
  colnames(SelectedGffChrMinusAnchor)[10] <- "ExtendedTSS"
  UpdatedGff <- rbind(UpdatedGff, SelectedGffChrPlusAnchor, SelectedGffChrMinusAnchor)
}
UpdatedGffeachloop <- UpdatedGff

#Loop extension of 5'end
for (k in 2:loop){
  UpdatedGff <- UpdatedGff[-c(1:nrow(UpdatedGff)),]
  colnames(UpdatedGff) <- c("Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute", "ExtendedTSS")
  for (j in names(table(SamR1Range$Chr))){
    BundledTagRangeLengthChr <- subset(BundledTagRangeLength, Chr==j)
    SelectedGffChr <- subset(UpdatedGffeachloop, Chr==j)
    BundledTagRangeLengthChrPlus <- subset(BundledTagRangeLengthChr, Direction=="Plus")
    SelectedGffChrPlus <- subset(SelectedGffChr, Strand=="+")
    FarthestAnchorTagPlus <- c()
    for (n in 1:nrow(SelectedGffChrPlus)){
      AnchorTagPlus <- subset(BundledTagRangeLengthChrPlus, Anchor>=(SelectedGffChrPlus[n,10]+ReadLength-1) & Anchor<=SelectedGffChrPlus[n,5])
      FarthestAnchorTagPlus <- c(FarthestAnchorTagPlus,min(AnchorTagPlus$ExtendedTSS))
    }
    SelectedGffChrPlusAnchor <- cbind(SelectedGffChrPlus[,1:9],FarthestAnchorTagPlus)
    colnames(SelectedGffChrPlusAnchor)[10] <- "ExtendedTSS"
    BundledTagRangeLengthChrMinus <- subset(BundledTagRangeLengthChr, Direction=="Minus")
    SelectedGffChrMinus <- subset(SelectedGffChr, Strand=="-")
    FarthestAnchorTagMinus <- c()
    for (n in 1:nrow(SelectedGffChrMinus)){
      AnchorTagMinus <- subset(BundledTagRangeLengthChrMinus, Anchor>=SelectedGffChrMinus[n,4] & Anchor<=(SelectedGffChrMinus[n,10]-ReadLength+1))
      FarthestAnchorTagMinus <- c(FarthestAnchorTagMinus,max(AnchorTagMinus$ExtendedTSS))
    }
    SelectedGffChrMinusAnchor <- cbind(SelectedGffChrMinus[,1:9], FarthestAnchorTagMinus)
    colnames(SelectedGffChrMinusAnchor)[10] <- "ExtendedTSS"
    UpdatedGff <- rbind(UpdatedGff, SelectedGffChrPlusAnchor, SelectedGffChrMinusAnchor)
  }
  UpdatedGffeachloop <- UpdatedGff
}
write.table(UpdatedGff, paste("ExtendedTSS_", gffFile, "_", loop,"timesloop.tab", sep=""), quote=F, sep="\t", row.names=F)
