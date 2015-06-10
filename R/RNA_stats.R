FONTSIZE=6
DATAPATH <- '~/projects/synergy/rnaseq/'

files <- paste(DATAPATH, seq(2749, 2758), '.genemeans', sep='')
annfile <- paste(DATAPATH, 'Homo_sapiens.GRCh37.68.cdna.new_ref.2749.counts.genes', sep='')

times <- c('0000', '0005', '0010', '0020', '0040',
           '0080', '0160', '0320', '0640', '1280')

files2 <- paste(DATAPATH, 'MCF7_uniq2_t',
                times,
                'min_RNA_v68_2012-03.bam.counts', sep='')

files3 <- paste(DATAPATH, 'MCF7_t',
                times,
                'min_RNA_v68_2012-03.summarycounts', sep='')

times.short <- c('0', '5', '10', '20', '40', '80', '160', '320', '640', '1280')

finalgenes <- unlist(read.table('../matlab/final_genes.txt'))

annotation <- read.table(annfile)
row.names(annotation) <- annotation[,1]
expmat <- c()
for (k in seq_along(files)) {
  t <- read.table(files[k], header=FALSE)
  expmat <- cbind(expmat, t[,'V1'])
}
expmat <- cbind(annotation[,2:3], expmat)
row.names(expmat) <- annotation[,1]

for (k in seq_along(files2)) {
  d <- dim(expmat)[2]+1
  t <- read.table(files2[k], header=FALSE, row.names=2)
  v <- t[grepl('.*_unspliced', row.names(t)),1]
  names(v) <- row.names(t)[grepl('.*_unspliced', row.names(t))]
  expmat[names(v),d] <- v
  expmat[is.na(expmat[,d]),d] <- 0
}

alignstats <- c()
for (k in seq_along(files3)) {
  t <- read.table(files3[k], header=TRUE)
  alignstats <- rbind(alignstats, t)
}
relativestats <- sweep(as.matrix(alignstats), 1, apply(alignstats, 1, sum), '/')
apply(relativestats, 2, mean)

mrnas <- apply(expmat[1:39071, 3:12], 1, mean)
names(mrnas) <- row.names(expmat)[1:39071]

mrnasum <- apply(expmat[1:39071, 3:12], 2, sum)
premrnasum <- apply(expmat[39072:78142, 3:12], 2, sum)

premrnas <- apply(expmat[39072:78142, 3:12], 1, mean)
names(premrnas) <- row.names(expmat)[39072:78142]
premrnagenes <- sapply(strsplit(names(premrnas), '_'), function (x) x[[1]])
names(premrnas) <- premrnagenes

premrnas2 <- apply(expmat[39072:78142, 13:22], 1, mean)
names(premrnas2) <- premrnagenes

mrnacov <- mrnas / expmat[1:39071, 2]
premrnacov <- premrnas / expmat[39072:78142, 2]

pdf('rnaseq_stats.pdf', width=87/25.4, height=100/25.4)
par(ps=FONTSIZE, cex=1)
par(mar=c(1.0, 0.8, 0, 0.8)+1.0)
par(mgp=c(0.6, 0.1, 0))
par(tck=-0.015)
par(mfrow=c(4, 2))
hist(log(premrnas,10), seq(-0.5, 7, by=0.5), main='All genes', xlab=expression('log'[10]*'(pre-mRNA count)'))
hist(log(premrnacov,10), seq(-6, 3.5, by=0.5), main='All genes', xlab=expression('log'[10]*'(pre-mRNA coverage)'))

hist(log(mrnas,10), seq(-0.5, 7, by=0.5), main='All genes', xlab=expression('log'[10]*'(mRNA count)'))
hist(log(mrnacov,10), seq(-6, 3.5, by=0.5), main='All genes', xlab=expression('log'[10]*'(mRNA coverage)'))

hist(log(premrnas[finalgenes],10), seq(-0.5, 7, by=0.5), main='Selected genes', xlab=expression('log'[10]*'(pre-mRNA count)'))
hist(log(premrnacov[finalgenes],10), seq(-6, 3.5, by=0.5), main='Selected genes', xlab=expression('log'[10]*'(pre-mRNA coverage)'))

hist(log(mrnas[finalgenes],10), seq(-0.5, 7, by=0.5), main='Selected genes', xlab=expression('log'[10]*'(mRNA count)'))
hist(log(mrnacov[finalgenes],10), seq(-6, 3.5, by=0.5), main='Selected genes', xlab=expression('log'[10]*'(mRNA coverage)'))
dev.off()

mrnagenes <- annotation[1:39071,1]
multiexon_genes <- which(annotation[1:39071,2]>1)
multiexon_genenames <- mrnagenes[annotation[1:39071,2]>1]
singleexon_genes <- mrnagenes[annotation[1:39071,2]==1]
multiexon_premrnas <- paste(multiexon_genenames, '_unspliced', sep='')

mrnasum2 <- apply(expmat[multiexon_genes, 3:12], 2, sum)
premrnasum2 <- apply(expmat[paste(multiexon_genenames, '_unspliced', sep=''), 3:12], 2, sum)

avelength.mrna <- mean(annotation[1:39071,3])
cat('Average mRNA length', avelength.mrna, '\n')
avelength.premrna <- mean(annotation[39072:78142,3])
cat('Average pre-mRNA length', avelength.premrna, '\n')
aveinvlength.mrna <- mean(1/annotation[multiexon_genes,3])
cat('Average mRNA inv length', aveinvlength.mrna, '\n')
aveinvlength.premrna <- mean(1/annotation[multiexon_premrnas,3])
cat('Average pre-mRNA inv length', aveinvlength.premrna, '\n')
averatio.premrna <- mean(premrnasum / (mrnasum + premrnasum))
alignmentmeans <- apply(relativestats, 2, mean)
alignmentmeans['pre.mRNA'] + alignmentmeans['both']*aveinvlength.premrna/(aveinvlength.premrna+aveinvlength.mrna)

bitseqstats <- cbind(mrnasum / (mrnasum + premrnasum), premrnasum / (mrnasum + premrnasum))
colnames(bitseqstats) <- c('mRNA', 'pre.mRNA')
bitseqstats2 <- cbind(mrnasum2 / (mrnasum2 + premrnasum2), premrnasum2 / (mrnasum2 + premrnasum2))
bitseqstats2 <- cbind(bitseqstats2, relativestats[,'pre.mRNA'] + relativestats[,'both']*aveinvlength.premrna/(aveinvlength.premrna+aveinvlength.mrna))
colnames(bitseqstats2) <- c('mRNA', 'pre.mRNA', 'pre.mRNA.pred')

row.names(bitseqstats2) <- paste(times.short, 'min')
row.names(relativestats) <- paste(times.short, 'min')

library(xtable)
tbl <- xtable(bitseqstats2)
show(tbl)
tbl2 <- xtable(relativestats[,1:3])
digits(tbl2) <- 3
show(tbl2)
