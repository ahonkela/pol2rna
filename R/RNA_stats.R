DATAPATH <- '~/projects/synergy/rnaseq/'

files <- paste(DATAPATH, seq(2749, 2758), '.genemeans', sep='')
annfile <- paste(DATAPATH, 'Homo_sapiens.GRCh37.68.cdna.new_ref.2749.counts.genes', sep='')

files2 <- paste(DATAPATH, 'MCF7_uniq2_t',
                c('0000', '0005', '0010', '0020', '0040',
                  '0080', '0160', '0320', '0640', '1280'),
                'min_RNA_v68_2012-03.bam.counts', sep='')

finalgenes <- unlist(read.table('../matlab/final_genes.txt'))

annotation <- read.table(annfile)
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


mrnas <- apply(expmat[1:39071, 3:12], 1, mean)
names(mrnas) <- row.names(expmat)[1:39071]

premrnas <- apply(expmat[39072:78142, 3:12], 1, mean)
names(premrnas) <- row.names(expmat)[39072:78142]
premrnagenes <- sapply(strsplit(names(premrnas), '_'), function (x) x[[1]])
names(premrnas) <- premrnagenes

premrnas2 <- apply(expmat[39072:78142, 13:22], 1, mean)
names(premrnas2) <- premrnagenes

mrnacov <- mrnas / expmat[1:39071, 2]
premrnacov <- premrnas / expmat[39072:78142, 2]

pdf('rnaseq_stats.pdf', 5, 5)
par(mfrow=c(2, 2))
hist(log(premrnas,10), seq(-0.5, 7, by=0.5), main='All genes', xlab=expression('log'[10]*'(pre-mRNA count)'))
hist(log(premrnacov,10), seq(-6, 3.5, by=0.5), main='All genes', xlab=expression('log'[10]*'(pre-mRNA coverage)'))

hist(log(mrnas,10), seq(-0.5, 7, by=0.5), main='All genes', xlab=expression('log'[10]*'(mRNA count)'))
hist(log(mrnacov,10), seq(-6, 3.5, by=0.5), main='All genes', xlab=expression('log'[10]*'(mRNA coverage)'))

par(mfrow=c(2, 2))
hist(log(premrnas[finalgenes],10), seq(-0.5, 7, by=0.5), main='Selected genes', xlab=expression('log'[10]*'(pre-mRNA count)'))
hist(log(premrnacov[finalgenes],10), seq(-6, 3.5, by=0.5), main='Selected genes', xlab=expression('log'[10]*'(pre-mRNA coverage)'))

hist(log(mrnas[finalgenes],10), seq(-0.5, 7, by=0.5), main='Selected genes', xlab=expression('log'[10]*'(mRNA count)'))
hist(log(mrnacov[finalgenes],10), seq(-6, 3.5, by=0.5), main='Selected genes', xlab=expression('log'[10]*'(mRNA coverage)'))
dev.off()

