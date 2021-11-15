#Step 1. Read in file of subcompartment designations.
setwd()
head(file)
colnames(file)<-c("chr","start","end","id")
library(GenomicRanges)
clustera1 <- with(subset(file,file[,4]=="a1") , GRanges(chr, IRanges(start, end)))
clustera2 <- with(subset(file,file[,4]=="a2") , GRanges(chr, IRanges(start, end)))
clusterb1 <- with(subset(file,file[,4]=="b1") , GRanges(chr, IRanges(start, end)))
clusterb2 <- with(subset(file,file[,4]=="b2") , GRanges(chr, IRanges(start, end)))

#####ENCODE concordance
#Step 2. Read in ENCODE files using GenomicRanges to calculate subcompartment percent coverage.
##Example below:
encode <- read.table("cortex-rnapol2.bed",header=TRUE,sep="\t")
encode[,3] <- encode[,2] + 1
encode.gr <- with(encode, GRanges(encode[,1],IRanges(encode[,2],encode[,3])))
sum(width(subsetByOverlaps(clustera1,encode.gr)))/sum(width(clustera1))
sum(width(subsetByOverlaps(clustera2,encode.gr)))/sum(width(clustera2))
sum(width(subsetByOverlaps(clusterb1,encode.gr)))/sum(width(clusterb1))
sum(width(subsetByOverlaps(clusterb2,encode.gr)))/sum(width(clusterb2))

#####Histone enrichment concordance
#Step 2. Read in histone enrichment files to calculate subcompartment association.
##Example below:
histone <- read.table("",header=TRUE,sep="\t")
up <- histone[which(histone$Event == "Up"),1:3]
down <- histone[which(histone$Event == "Down"),1:3]

#Step 3. Bin the mm10 genome into 250kb blocks.
library(GenomicRanges)
chrSizes<- c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566)
names(chrSizes) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19")
bins   <- tileGenome(chrSizes, tilewidth=2.5e5, cut.last.tile.in.chrom=T)
bins

#Step 4. Perform Fisher's 2x2 testing for histone enrichments per subcompartment (>99th%ile)
library(dplyr)
library(exact2x2)
neuron<-NULL

clustera1.yes <- countOverlaps(clustera1,file)
clustera1.not <- countOverlaps(bins[!bins %in% (clustera1),],file)

clustera2.yes <- countOverlaps(clustera2,file)
clustera2.not <- countOverlaps(bins[!bins %in% (clustera2),],file)

clusterb1.yes <- countOverlaps(clusterb1,file)
clusterb1.not <- countOverlaps(bins[!bins %in% (clusterb1),],file)

clusterb2.yes <- countOverlaps(clusterb2,file)
clusterb2.not <- countOverlaps(bins[!bins %in% (clusterb2),],file)

samples.y <- mget(ls(pattern="yes"))
samples.n <- mget(ls(pattern="not"))

for (i in 1:length(samples.y)) {
 for (j in 1:length(samples.n)) {
    if (i != j) {
      next
    }
    else {
      yy <- sum(as.numeric(unlist(samples.y[i])) > quantile(overlaps, .99))
      yn <- sum(as.numeric(unlist(samples.y[i])) <= quantile(overlaps, .99))
      ny <- sum(as.numeric(unlist(samples.n[j])) > quantile(overlaps, .99))
      nn <- sum(as.numeric(unlist(samples.n[j])) <= quantile(overlaps, .99))
      data <- matrix(c(nn,ny,yn,yy), 2, 2, byrow=TRUE)
      neuron <- cbind(neuron,fisher.exact(data, y = NULL, or = 1, alternative = "greater",
                                                    tsmethod = "minlike", conf.int = TRUE, conf.level = 0.95,
                                                    tol = 0.00001, midp=FALSE)$p.value)
    }
  }
}

head(neuron)