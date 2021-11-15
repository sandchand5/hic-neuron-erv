#Step 1. Select [subcomp] file and Hi-C interaction [file].
setwd()

#Step 2. Create separate [anchor] and [target] data frames.


#Step 3. Read in 1Mb genome bins.
setwd("~/Desktop")
library(GenomicRanges)
chrSizes<- c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566)
names(chrSizes) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19")
bins   <- tileGenome(chrSizes, tilewidth=1e6, cut.last.tile.in.chrom=T)

#Step 4. Classify each Hi-C interaction by subcompartment designations of anchor and target.
seg1 <- NULL
for (n in 1:length(anchor)) {
  counts <- cbind.data.frame(countOverlaps(anchor[n],clustera1),countOverlaps(anchor[n],clustera2),countOverlaps(anchor[n],clusterb1),countOverlaps(anchor[n],clusterb2))
  seg1 <- rbind.data.frame(seg1,counts)
}

seg2 <- NULL
for (n in 1:length(target)) {
  counts <- cbind.data.frame(countOverlaps(target[n],clustera1),countOverlaps(target[n],clustera2),countOverlaps(target[n],clusterb1),countOverlaps(target[n],clusterb2))
  seg2 <- rbind.data.frame(seg2,counts)
}

unique(seg1)
anchorval <- NULL
for (m in 1:nrow(seg1)) {
      if (max(seg1[m,]) == seg1[m,4]) {
        val <- "d"
      }
      else {
        if (max(seg1[m,]) == seg1[m,2]) {
          val <- "c"
        }
        else {
          if (max(seg1[m,]) == seg1[m,1]) {
            val <- "b"
          }
          else {
            if (max(seg1[m,]) == seg1[m,3]) {
              val <- "a"
            }
      }
    }
  }
  anchorval <- rbind(anchorval,val)
}
  anchor$subcomp <- anchorval

targetval <- NULL
for (m in 1:nrow(seg2)) {
      if (max(seg2[m,]) == seg2[m,4]) {
        val <- "d"
      }
      else {
        if (max(seg2[m,]) == seg2[m,2]) {
          val <- "c"
        }
        else {
          if (max(seg2[m,]) == seg2[m,1]) {
            val <- "b"
          }
          else {
            if (max(seg2[m,]) == seg2[m,3]) {
              val <- "a"
            }
          }
        }
      }
  targetval <- rbind(targetval,val)
}
  target$subcomp <- targetval

finalvalue <- cbind.data.frame(anchorval,targetval)
clusters <- c("a","b","c","d")
table(apply(finalvalue,1,function(x) paste(sort(x),collapse='-')))