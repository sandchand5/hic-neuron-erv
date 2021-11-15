#Step 1. Set working directory and read in human subcompartments as [subcomp] (Rao, et. al).
##Note: Coordinates were in hg19 and were converted to hg38 using the UCSC liftOver tool.
setwd()
subcomp <- read.table("GSE63525_GM12878_subcompartments.bed",header=F,sep="\t")
subcomp <- subcomp[,1:4]
A1 <- subcomp[which(subcomp[,4] == "A1"),]
A2 <- subcomp[which(subcomp[,4] == "A2"),]
B1 <- subcomp[which(subcomp[,4] == "B1"),]
B2 <- subcomp[which(subcomp[,4] == "B2"),]
B3 <- subcomp[which(subcomp[,4] == "B3"),]
B4 <- subcomp[which(subcomp[,4] == "B4"),]

write.table(A1,"a1-rao-coords.bed",quote=F,row.names=F,col.names=F)
write.table(A2,"a2-rao-coords.bed",quote=F,row.names=F,col.names=F)
write.table(B1,"b1-rao-coords.bed",quote=F,row.names=F,col.names=F)
write.table(B2,"b2-rao-coords.bed",quote=F,row.names=F,col.names=F)
write.table(B3,"b3-rao-coords.bed",quote=F,row.names=F,col.names=F)
write.table(B4,"b4-rao-coords.bed",quote=F,row.names=F,col.names=F)

#####UCSC Genome Browser liftOver tool#####

#Step 2. Set working directory and read in hg38 liftOver subcompartments.
setwd()
A1 <- read.table("a1-rao-coords-liftedover.bed",header=F,sep="\t")
A2 <- read.table("a2-rao-coords-liftedover.bed",header=F,sep="\t")
B1 <- read.table("b1-rao-coords-liftedover.bed",header=F,sep="\t")
B2 <- read.table("b2-rao-coords-liftedover.bed",header=F,sep="\t")
B3 <- read.table("b3-rao-coords-liftedover.bed",header=F,sep="\t")
B4 <- read.table("b4-rao-coords-liftedover.bed",header=F,sep="\t")

A1 <- with(A1,GRanges(A1[,1],IRanges(A1[,2],A1[,3])))
A2 <- with(A2,GRanges(A2[,1],IRanges(A2[,2],A2[,3])))
B1 <- with(B1,GRanges(B1[,1],IRanges(B1[,2],B1[,3])))
B2 <- with(B2,GRanges(B2[,1],IRanges(B2[,2],B2[,3])))
B3 <- with(B3,GRanges(B3[,1],IRanges(B3[,2],B3[,3])))
B4 <- with(B4,GRanges(B4[,1],IRanges(B4[,2],B4[,3])))

#Step 3. Set working directory and read in mm10 subcompartments as [file].
setwd()
head(file)
library(GenomicRanges)
all <- with(file, GRanges(chr, IRanges(start, end)))
clustera1 <- with(subset(file,file[,4]=="a1") , GRanges(chr, IRanges(start, end)))
clustera2 <- with(subset(file,file[,4]=="a2") , GRanges(chr, IRanges(start, end)))
clusterb1 <- with(subset(file,file[,4]=="b1") , GRanges(chr, IRanges(start, end)))
clusterb2 <- with(subset(file,file[,4]=="b2") , GRanges(chr, IRanges(start, end)))
head(all)

#Step 4. Calculate observed percent overlaps. Example below:
sum(width(subsetByOverlaps(A1,clustera1)))

#####Divide by calculated expected overlap to create [matrix] of obs/exp.#####

#Step 5. Create [matrix] plot using corrplot.
library(corrplot)
corrplot(as.matrix(matrix),is.cor=F,cl.lim = c(-25, 25),col=colorRampPalette(c("blue","white","red"))(200))
