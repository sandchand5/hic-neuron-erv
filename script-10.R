#Step 0. Create syntenic hg38 files fo mm10 subcompartments using UCSC Genome Browser.

#Step 1. Read/format syntenic [syn] file for subcompartment of interest.
setwd()
head(syn)
colnames(syn) <- c("Chr","Coords")

library(tidyr)
syn <- syn %>% separate(Coords, c("Start", "End"), sep="-")
head(syn)
syn[,2] <- as.numeric(syn[,2])
syn[,3] <- as.numeric(syn[,3])

#Step 2. Bin human hg38 genome.
chrom <- read.table("hg38-chrom-auto-sizes.txt",header=FALSE,sep="\t")
head(chrom)
library(GenomicRanges)
chrSizes<- chrom[,2]
names(chrSizes) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
bins   <- tileGenome(chrSizes, tilewidth=2.5e5, cut.last.tile.in.chrom=T)
bins

#Step 3. Classify bins as syntenic to subcompartment of interest or not.
syn.bins <- subsetByOverlaps(bins,syn)

#Step 4. Read in hg38 Repeatmasker [repeats] file.
setwd()
head(repeats)
ERVK <- with(repeats[grep("ERVK",repeats$type),] , GRanges(chr, IRanges(start, end)))
ERV1 <- with(repeats[grep("ERV1",repeats$type),] , GRanges(chr, IRanges(start, end)))
ERVL <- with(repeats[grep("ERVL",repeats$type),] , GRanges(chr, IRanges(start, end)))
LINE <- with(repeats[grep("LINE",repeats$class),] , GRanges(chr, IRanges(start, end)))
SINE <- with(repeats[grep("SINE",repeats$class),] , GRanges(chr, IRanges(start, end)))

#Step 5. Set up Fisher 2x2 testing with 'high' (>99%ile) and 'low' repeat densities. Example below:
dens <- countOverlaps(bins,ERVK)
bins$overlaps <- dens
bins.subset <- bins[bins$overlaps >= quantile(dens,0.99) ]
bins.subset

a<-sum(countOverlaps(bins.subset,syn.bins))
b<-length(bins.subset)
c<-length(syn.bins)
d<-length(bins)

library(exact2x2)
data <- matrix(c(a,(b-a),(c-a),(d-c-b+a)),nrow=2,ncol=2)
fisher.exact(data,alternative="two.sided")  $p.value 