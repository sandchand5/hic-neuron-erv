#Step 1. Read in aligned PacBio .bed files, filter for SCORE > 60, and convert to GRanges (gr).
##Note. Keep read IDs associated with coordinates.
library(GenomicRanges)
setwd()
read.table()

#Step 2. Calculate autosomal integrations (10kb resolution).
chrSizes<- c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566,171031299,91744698)
names(chrSizes) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")
bins   <- tileGenome(chrSizes, tilewidth=1e4, cut.last.tile.in.chrom=T)
subsetByOverlaps(gr,bins)

#Step 3. Find number of [IAPEzi] and non-IAPEzi [IAP] insertions at known coordinates.
##Note. Read in mm10 Repeatmasker ERVK track.

length(subsetByOverlaps(gr,iapezi))
length(subsetByOverlaps(gr,iap)) - length(subsetByOverlaps(gr,iapezi))

#Step 4. Find de novo insertions and export ids to blast (UNIX).
length(subsetByOverlaps(bins,subsetByOverlaps(gr,iap,invert=TRUE)))
subsetByOverlaps(gr,iap,invert=TRUE)

#####BLAST#####

#Step 5. Read in .bed files of [de novo] reads and [file] of subcompartment [cluster] designations.
#Step 6. Calculate observed and expected de novo reads per subcompartment.
length(subsetByOverlaps(denovo,cluster))
