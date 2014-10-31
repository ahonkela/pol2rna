# Reads exon ranges from UCSC server and stores them in RData file
.libPaths("/share/mi/workspace/jtpelto/synergy/bioconductor")
library(GenomicFeatures)

# Creates a TranscriptDb object from transcript annotations available at the UCSC Genome Browser.
txdb <- makeTranscriptDbFromUCSC(genome="hg19", tablename="ensGene") 

# Stores exon data as a GRangesList object.
exonRanges <- exonsBy(txdb, "gene") 

save(file="exons.RData", exonRanges)
saveFeatures(txdb, "exons.sqlite")






# set partid to integer from 1 to 10

.libPaths("/share/mi/workspace/jtpelto/synergy/bioconductor")
library(GenomicFeatures)
library(R.matlab)

load('exons.RData')
# load('counts_all.RData')

nparts=10;

ngenes = length(exonRanges);

istart=floor(ngenes*(partid-1)/nparts+1)
iend=floor(ngenes*partid/nparts);
genedata = matrix(nrow=(iend-istart+1), ncol=5);

i = istart;
while (i <=  iend) {
  temp1 <- exonRanges[[i]]
  # gene name
  genedata[i-istart+1,1] <- names(exonRanges)[i]
  # chromosome name
  genedata[i-istart+1,2] <- toString(head(seqnames(temp1),1))
  # strand id
  genedata[i-istart+1,3] <- toString(head(strand(temp1),1))
  # start of exons
  genedata[i-istart+1,4] <- min(start(temp1))
  # end of exons
  genedata[i-istart+1,5] <- max(end(temp1))
        if (i %% 10 == 1) {
        print(c("Iteration", i));
}
i <- i + 1
}

colnames(genedata) <- c('ensemblid','chr','strand','start','end')

# genes = cbind(genedata, rawcounts)
# rownames(genes) <- rownames(rawcounts)

filename1<-sprintf("genes_dec2012_part%d.RData",partid)
save(file=filename1, genedata)

filename2<-sprintf("genes_dec2012_part%d.txt",partid)
write.table(genedata, filename2)

filename3<-sprintf("genes_dec2012_part%d.mat",partid)
writeMat(filename3, genedata=genedata)

