library(h5r)
library(Gviz)
library(GenomicFeatures)
library(org.Hs.eg.db)

readwig <- function(h, c) {
  return (matrix(readH5Data(getH5Dataset(h, c$name)), c$dim[1], c$dim[2], byrow=TRUE))
}

mygene <- 'DLX3'
EG <- get(mygene, org.Hs.egSYMBOL2EG)
myensg <- get(EG, org.Hs.egENSEMBL)
mychr <- get(EG, org.Hs.egCHR)
mystrand <- '-'
if (!exists("h5contents")) {
  h <- H5File('../python/wigs.h5', 'r')
  h5contents <- listH5Contents(h)
}
data <- readwig(h, h5contents[[myensg]])

if (!exists("trdb")) {
  trdb <- makeTranscriptDbFromBiomart(biomart="ensembl",
                                      dataset="hsapiens_gene_ensembl",
                                      filter=list("ensembl_gene_id"=myensg))
}

T <- c(0, 5, 10, 20, 40, 80, 160, 320, 640, 1280)

options(ucscChromosomeNames=FALSE)
gtrack <- GeneRegionTrack(trdb, name=mygene)

tracks <- vector("list", 12)
tracks[[1]] <- gtrack
tracks[[12]] <- GenomeAxisTrack()
mystart <- min(start(transcripts(trdb)))
myend <- max(end(transcripts(trdb)))
##starts <- seq(mystart, myend, by=200)
starts <- seq(myend, mystart, by=-200)
##ends <- starts[2:(dim(data)[1]+1)]-1
starts <- starts[1:dim(data)[1]]
for (k in seq(10)) {
  mydata <- rep(0, length(starts))
  mydata[1:dim(data)[1]] <- data[,k]
  tracks[[k+1]] <- DataTrack(data=mydata, start=starts-200, width=200,
                             strand=mystrand, chromosome=mychr, genome="hg19",
                             name=paste(T[k], "min"), size=2)
}
plotTracks(tracks[c(1, seq(4, 9))], type='l', reverseStrand=TRUE)

##dtrack <- DataTrack(data=t(data), start=starts[1:dim(data)[1]]-200, width=200, strand=mystrand, chromosome=mychr, genome="hg19")
##plotTracks(dtrack)

pdf('dlx3_premrna.pdf', width=87/25.4, height=100/25.4)
gpar(fontsize=12)
plotTracks(tracks[c(1, seq(4, 9))], type='l', reverseStrand=TRUE)
dev.off()
