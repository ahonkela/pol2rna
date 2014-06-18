library(h5r)
library(Gviz)
library(GenomicFeatures)

readwig <- function(h, c) {
  return (matrix(readH5Data(getH5Dataset(h, c$name)), c$dim[1], c$dim[2], byrow=TRUE))
}

h <- H5File('../python/wigs.h5', 'r')
contents <- listH5Contents(h)

LHX4 <- 'ENSG00000121454'
data <- readwig(h, contents[[LHX4]])

trdb <- makeTranscriptDbFromBiomart(biomart="ensembl",
                                    dataset="hsapiens_gene_ensembl",
                                    filter=list("ensembl_gene_id"=LHX4))

options(ucscChromosomeNames=FALSE)
gtrack <- GeneRegionTrack(trdb, name="LHX4")

mystart <- min(start(transcripts(trdb)))
myend <- max(end(transcripts(trdb)))
starts <- seq(mystart, myend, by=200)
ends <- starts[2:(dim(data)[1]+1)]-1
starts <- starts[1:dim(data)[1]]
dtrack <- DataTrack(data=data[,1], start=starts, end=ends,
                    strand='+', chromosome=1, genome='hg19', name='pre-mRNA')
plotTracks(list(gtrack, dtrack), transcriptAnnotation = "transcript", type='l')
