DATAPATH <- '~/projects/synergy/rnaseq/'

files <- paste(DATAPATH, seq(2749, 2758), '.genemeans', sep='')
annfile <- paste(DATAPATH, 'Homo_sapiens.GRCh37.68.cdna.new_ref.2749.counts.genes', sep='')

finalgenes <- unlist(read.table('../matlab/final_genes.txt'))

annotation <- read.table(annfile)
expmat <- c()
for (k in seq_along(files)) {
  t <- read.table(files[k], header=FALSE)
  expmat <- cbind(expmat, t[,'V1'])
}
expmat <- cbind(annotation, expmat)
mrnas <- apply(expmat[1:39071, 4:13], 1, mean)
names(mrnas) <- expmat[1:39071, 1]

premrnas <- apply(expmat[39072:78142, 4:13], 1, mean)
names(premrnas) <- expmat[39072:78142, 1]
premrnagenes <- sapply(strsplit(names(premrnas), '_'), function (x) x[[1]])
names(premrnas) <- premrnagenes

mrnacov <- mrnas / expmat[1:39071, 3]
premrnacov <- premrnas / expmat[39072:78142, 3]

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

