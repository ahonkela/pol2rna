delays.pol2 <- read.table('pol2max_and_meddelays_2013-08-30.txt', row.names=1, header=TRUE)
delays.premrna <- read.table('pol2max_and_meddelays_2013-11-05.txt', row.names=1, header=TRUE)
premrna.fits <- read.table('../python/premrna_halfdiff_2014-06-11.txt', row.names=1, header=FALSE)
names(premrna.fits) <- 'premrna_trend'
##delays <- read.table('pol2max_and_delays_2013-03-11.txt', row.names=1, header=TRUE)
delays.orig <- delays.pol2
delays.joint <- merge(delays.pol2, delays.premrna, by=0)
row.names(delays.joint) <- delays.joint[,'Row.names']
F <- as.character(delays.joint[,'gene.x']) == as.character(delays.joint[,'gene.y'])
stopifnot(all(F[!is.na(F)]))
stopifnot(all(delays.joint['corr.x'] == delays.joint['corr.y']))
delays.orig <- delays.joint[!names(delays.joint) %in% c('Row.names', 'gene.y', 'corr.y')]
n <- names(delays.orig)
n[grep('corr.x', n)] <- 'corr'
n[grep('gene.x', n)] <- 'gene'
names(delays.orig) <- n
delays2 <- merge(delays.orig, premrna.fits, by=0)
row.names(delays2) <- delays2[,'Row.names']
delays2 <- delays2[!names(delays2) %in% c('Row.names')]
delays.orig <- delays2

t <- readLines('intron_lengths.txt.lengths2')
l <- strsplit(t, '\t')
names(l) <- sapply(l, function(x) x[1])
introns <- lapply(l, function(x) as.integer(x[-1:-6]))
strands <- sapply(l, function(x) x[2])
termpos <- sapply(l, function(x) as.integer(x[3]))
lengths <- sapply(l, function(x) as.integer(x[4]))
exonlen5 <- sapply(l, function(x) as.integer(x[5]))
exonlen3 <- sapply(l, function(x) as.integer(x[6]))

## t <- readLines('intron_lengths.txt')
## l <- strsplit(t, '\t')
## names(l) <- sapply(l, function(x) x[1])
## introns <- lapply(l, function(x) as.integer(x[-1]))

numIntrons <- unlist(sapply(introns, length))

introns <- introns[numIntrons > 0]
termpos <- termpos[numIntrons > 0]
strands <- strands[numIntrons > 0]
lengths <- lengths[numIntrons > 0]
exonlen5 <- exonlen5[numIntrons > 0]
exonlen3 <- exonlen3[numIntrons > 0]
numIntrons <- numIntrons[numIntrons > 0]

lastIntronsTr <- unlist(sapply(introns, function(x) x[length(x)]))

maxIntronsTr <- unlist(sapply(introns, function(x) max(x)))
maxmaxIntrons <- sapply(split(maxIntronsTr, substr(names(maxIntronsTr), 1, 15)), max)
maxmaxIntrons <- maxmaxIntrons[names(maxmaxIntrons) %in% row.names(delays.orig)]

maxmaxLastIntrons <- sapply(split(lastIntronsTr, substr(names(lastIntronsTr), 1, 15)), max)

maxmaxLastIntrons <- maxmaxLastIntrons[names(maxmaxLastIntrons) %in% row.names(delays.orig)]

maxTrLengths <- sapply(split(lengths, substr(names(lengths), 1, 15)), max)
maxTrLengths <- maxTrLengths[names(maxTrLengths) %in% row.names(delays.orig)]

maxExon3Lengths <- sapply(split(exonlen3, substr(names(exonlen3), 1, 15)), max)
maxExon3Lengths <- maxExon3Lengths[names(maxExon3Lengths) %in% row.names(delays.orig)]

maxExon5Lengths <- sapply(split(exonlen5, substr(names(exonlen5), 1, 15)), max)
maxExon5Lengths <- maxExon5Lengths[names(maxExon5Lengths) %in% row.names(delays.orig)]

longestLast <- (maxmaxLastIntrons == maxmaxIntrons)
lastProportion <- (maxmaxLastIntrons / maxmaxIntrons)
myfact <- rep(0, length(lastProportion))

delays <- merge(delays.orig, cbind(maxmaxLastIntrons, maxmaxIntrons, longestLast, maxTrLengths, lastProportion, maxExon3Lengths, maxExon5Lengths, myfact), by=0)
row.names(delays) <- delays[,1]
delays <- delays[,-1]
I <- (delays[,'tmax.x'] < 160) & (delays[,'tmax.x'] > 1) & (delays[,'begdev10.x'] < 12) & (delays[,'meddelay.x'] < 120) #& (delays[,'corr'] > 0.5)
I2 <- (delays[,'tmax.x'] < 160) & (delays[,'tmax.x'] > 1) & (delays[,'begdev10.x'] < 12) & (delays[,'meddelay.x'] < 80) #& (delays[,'corr'] > 0.5)
I3 <- (delays[,'tmax.x'] < 160) & (delays[,'tmax.x'] > 1) & (delays[,'begdev10.x'] < 12) & (delays[,'meddelay.x'] < 60) #& (delays[,'corr'] > 0.5)
mydelays <- delays[I,]
mydelays[mydelays['corr']<0.5 & mydelays['meddelay.x']>15, 'myfact'] <- 1
mydelays[mydelays['corr']>0.5 & mydelays['meddelay.x']>15, 'myfact'] <- 2
mydelays[mydelays['corr']<0.5 & mydelays['meddelay.x']<15, 'myfact'] <- 3
mydelays[mydelays['corr']>0.5 & mydelays['meddelay.x']<15, 'myfact'] <- 4
longestLast <- mydelays[,'longestLast']
##lastProportion <- mydelays[,'maxmaxLastIntrons']/mydelays[,'maxTrLengths']
lastProportion <- mydelays[,'lastProportion']/mydelays[,'maxmaxIntrons']
##long <- longestLast[mydelays['meddelay.x']>15,]
##short <- longestLast[mydelays['meddelay.x']<15,]
##longcorr <- longestLast[mydelays['corr']>0.5 & mydelays['meddelay.x']>15]
##longuncorr <- longestLast[mydelays['corr']<0.5 & mydelays['meddelay.x']>15]
##shortcorr <- longestLast[mydelays['corr']>0.5 & mydelays['meddelay.x']<15]
##shortuncorr <- longestLast[mydelays['corr']<0.5 & mydelays['meddelay.x']<15]
longcorr <- mydelays[mydelays['corr']>0.5 & mydelays['meddelay.x']>15, 'lastProportion']
longuncorr <- mydelays[mydelays['corr']<0.5 & mydelays['meddelay.x']>15, 'lastProportion']
shortcorr <- mydelays[mydelays['corr']>0.5 & mydelays['meddelay.x']<15, 'lastProportion']
shortuncorr <- mydelays[mydelays['corr']<0.5 & mydelays['meddelay.x']<15, 'lastProportion']

lencos <- c(0.5, 0.75, 0.9, 0.95)
pdf('delay_survival.pdf', width=87/25.4, height=70/25.4)
par(ps=8, cex=1)
par(mar=c(2, 2, 0, 2)+0.4)
par(mgp=c(1.2, 0.4, 0))
NORM <- 10
MAXVAL <- 0.2
for (l in seq_along(lencos)) {
  lenco <- lencos[l]

  shortlast <- mydelays[mydelays['lastProportion']<lenco,'meddelay.x']
  longlast <- mydelays[mydelays['lastProportion']>lenco,'meddelay.x']

  T <- seq(0, 80)
  shortfreq <- rep(0, length(T))
  longfreq <- rep(0, length(T))
  counts <- matrix(0, length(T), 4)
  pvals <- rep(1, length(T))
  for (k in seq_along(T)) {
    shortfreq[k] <- mean(shortlast > T[k])
    longfreq[k] <- mean(longlast > T[k])
    counts[k, 1] <- sum(shortlast < T[k])
    counts[k, 2] <- sum(shortlast >= T[k])
    counts[k, 3] <- sum(longlast < T[k])
    counts[k, 4] <- sum(longlast >= T[k])
    pvals[k] <- fisher.test(matrix(counts[k,],2))$p.value
    ##cat(sum(shortlast > T[k]) + sum(longlast > T[k]), '\n')
  }
  plot(T, shortfreq, type='l', col='blue', ylim=c(0, 0.2), xlim=c(0, 80), ylab=expression("Fraction of genes with" ~ Delta > t), xlab="t (min)")
  lines(T, longfreq, col='red')
  lines(T, rep(-log(0.05)/log(10)/NORM*MAXVAL, length(T)), col='black', lty=2)
  axis(4, seq(0, MAXVAL, len=6), seq(0, NORM, by=2))
  lines(T, -log(pvals)/log(10)/NORM*MAXVAL, col='black')
  mtext(expression(-log[10](p-value)), side=4, line=1.2)
  legend('topright',
         legend=c(sprintf("f<%.2f, N=%d", lenco, length(shortlast)),
           sprintf("f>%.2f, N=%d", lenco, length(longlast)),
           'p-value'),
         col=c('blue', 'red', 'black'), lty=1)
}

lencos <- c(1e4, 3e4, 1e5)
NORM <- 10
MAXVAL <- 0.2
for (l in seq_along(lencos)) {
  lenco <- lencos[l]

  shortlast <- mydelays[mydelays['maxTrLengths']<lenco,'meddelay.x']
  longlast <- mydelays[mydelays['maxTrLengths']>lenco,'meddelay.x']

  T <- seq(0, 80)
  shortfreq <- rep(0, length(T))
  longfreq <- rep(0, length(T))
  counts <- matrix(0, length(T), 4)
  pvals <- rep(1, length(T))
  for (k in seq_along(T)) {
    shortfreq[k] <- mean(shortlast > T[k])
    longfreq[k] <- mean(longlast > T[k])
    counts[k, 1] <- sum(shortlast < T[k])
    counts[k, 2] <- sum(shortlast >= T[k])
    counts[k, 3] <- sum(longlast < T[k])
    counts[k, 4] <- sum(longlast >= T[k])
    pvals[k] <- fisher.test(matrix(counts[k,],2))$p.value
    ##cat(sum(shortlast > T[k]) + sum(longlast > T[k]), '\n')
  }
  plot(T, shortfreq, type='l', col='blue', ylim=c(0, 0.2), ylab=expression("Fraction of genes with" ~ Delta > t), xlab="t (min)")
  lines(T, longfreq, col='red')
  lines(T, rep(-log(0.05)/log(10)/NORM*MAXVAL, length(T)), col='black', lty=2)
  axis(4, seq(0, MAXVAL, len=6), seq(0, NORM, by=2))
  lines(T, -log(pvals)/log(10)/NORM*MAXVAL, col='black')
  mtext(expression(-log[10](p-value)), side=4, line=1.2)
  legend('topright',
         legend=c(sprintf("m<%.0f, N=%d", lenco, length(shortlast)),
           sprintf("m>%.0f, N=%d", lenco, length(longlast)),
           'p-value'),
         col=c('blue', 'red', 'black'), lty=1)
}
dev.off()


lencos <- c(0.5, 0.75, 0.9, 0.95)
pdf('corr_survival.pdf', width=87/25.4, height=70/25.4)
par(ps=8, cex=1)
par(mar=c(2, 2, 0, 2)+0.4)
par(mgp=c(1.2, 0.4, 0))
NORM <- 10
MAXVAL <- 1
for (l in seq_along(lencos)) {
  lenco <- lencos[l]

  shortlast <- mydelays[mydelays['lastProportion']<lenco,'corr']
  longlast <- mydelays[mydelays['lastProportion']>lenco,'corr']

  T <- seq(-1, 1, len=100)
  shortfreq <- rep(0, length(T))
  longfreq <- rep(0, length(T))
  counts <- matrix(0, length(T), 4)
  pvals <- rep(1, length(T))
  for (k in seq_along(T)) {
    shortfreq[k] <- mean(shortlast < T[k])
    longfreq[k] <- mean(longlast < T[k])
    counts[k, 1] <- sum(shortlast < T[k])
    counts[k, 2] <- sum(shortlast >= T[k])
    counts[k, 3] <- sum(longlast < T[k])
    counts[k, 4] <- sum(longlast >= T[k])
    pvals[k] <- fisher.test(matrix(counts[k,],2))$p.value
    ##cat(sum(shortlast > T[k]) + sum(longlast > T[k]), '\n')
  }
  plot(T, shortfreq, type='l', col='blue', ylim=c(-0.005, 1.005), ylab=expression("Fraction of genes with" ~ rho < r), xlab="r")
  lines(T, longfreq, col='red')
  lines(T, rep(-log(0.05)/log(10)/NORM*MAXVAL, length(T)), col='black', lty=2)
  axis(4, seq(0, MAXVAL, len=6), seq(0, NORM, by=2))
  lines(T, -log(pvals)/log(10)/NORM*MAXVAL, col='black')
  mtext(expression(-log[10](p-value)), side=4, line=1.2)
  legend('topleft',
         legend=c(sprintf("f<%.2f, N=%d", lenco, length(shortlast)),
           sprintf("f>%.2f, N=%d", lenco, length(longlast)),
           'p-value'),
         col=c('blue', 'red', 'black'), lty=1)
}

lencos <- c(1e4, 3e4, 1e5)
NORM <- 30
MAXVAL <- 1
for (l in seq_along(lencos)) {
  lenco <- lencos[l]

  shortlast <- mydelays[mydelays['maxTrLengths']<lenco,'corr']
  longlast <- mydelays[mydelays['maxTrLengths']>lenco,'corr']

  T <- seq(-1, 1, len=100)
  shortfreq <- rep(0, length(T))
  longfreq <- rep(0, length(T))
  counts <- matrix(0, length(T), 4)
  pvals <- rep(1, length(T))
  for (k in seq_along(T)) {
    shortfreq[k] <- mean(shortlast < T[k])
    longfreq[k] <- mean(longlast < T[k])
    counts[k, 1] <- sum(shortlast < T[k])
    counts[k, 2] <- sum(shortlast >= T[k])
    counts[k, 3] <- sum(longlast < T[k])
    counts[k, 4] <- sum(longlast >= T[k])
    pvals[k] <- fisher.test(matrix(counts[k,],2))$p.value
    ##cat(sum(shortlast > T[k]) + sum(longlast > T[k]), '\n')
  }
  plot(T, shortfreq, type='l', col='blue', ylim=c(-0.005, 1.005), ylab=expression("Fraction of genes with" ~ rho < r), xlab="r")
  lines(T, longfreq, col='red')
  lines(T, rep(-log(0.05)/log(10)/NORM*MAXVAL, length(T)), col='black', lty=2)
  axis(4, seq(0, MAXVAL, len=7), seq(0, NORM, by=5))
  lines(T, -log(pvals)/log(10)/NORM*MAXVAL, col='black')
  mtext(expression(-log[10](p-value)), side=4, line=1.2)
  legend('topleft',
         legend=c(sprintf("m<%.0f, N=%d", lenco, length(shortlast)),
           sprintf("m>%.0f, N=%d", lenco, length(longlast)),
           'p-value'),
         col=c('blue', 'red', 'black'), lty=1)
}
dev.off()


## cat(c(sum(longcorr), sum(!longcorr), sum(shortcorr), sum(!shortcorr)), '\n')
## cat(c(sum(longuncorr), sum(!longuncorr), sum(shortuncorr), sum(!shortuncorr)), '\n')
## cat(c(sum(longcorr) / length(longcorr), sum(shortcorr)/length(shortcorr)), '\n')
## cat(c(sum(longuncorr) / length(longuncorr), sum(shortuncorr)/length(shortuncorr)), '\n')
## fisher.test(matrix(c(sum(longcorr), sum(!longcorr), sum(shortcorr), sum(!shortcorr)), 2))
## fisher.test(matrix(c(sum(longuncorr), sum(!longuncorr), sum(shortcorr), sum(!shortcorr)), 2))

## fisher.test(matrix(c(sum(longuncorr) + sum(longcorr),
##                      sum(!longuncorr) + sum(!longcorr),
##                      sum(shortcorr) + sum(shortuncorr),
##                      sum(!shortcorr) + sum(!shortuncorr)), 2))

## boxplot(list(longuncorr=longuncorr, longcorr=longcorr, shortuncorr=shortuncorr, shortcorr=shortcorr))

library(vioplot)
##vioplot(longuncorr, longcorr, shortuncorr, shortcorr, names=expression(Delta > 15 ~ "min," ~ rho < 0.5, Delta > 15 ~ "min," ~ rho > 0.5, Delta < 15 ~ "min," ~ rho < 0.5, Delta < 15 ~ "min," ~ rho > 0.5))

pdf('intron_lengths_med.pdf', width=178/25.4, height=70/25.4)
par(mfrow=c(1, 3))

longcorr <- log(mydelays[mydelays['corr']>0.5 & mydelays['meddelay.x']>15, 'maxTrLengths']) / log(10)
longuncorr <- log(mydelays[mydelays['corr']<0.5 & mydelays['meddelay.x']>15, 'maxTrLengths']) / log(10)
shortcorr <- log(mydelays[mydelays['corr']>0.5 & mydelays['meddelay.x']<15, 'maxTrLengths']) / log(10)
shortuncorr <- log(mydelays[mydelays['corr']<0.5 & mydelays['meddelay.x']<15, 'maxTrLengths']) / log(10)

plot.new()
par(ps=8, cex=1)
par(mar=c(3, 3, 0, 0)+0.4)
par(mgp=c(2, 1, 0))
plot.window(xlim=c(0.5, 4.5), ylim=c(2.8, 6.1))
axis(1, at=seq(4), labels=c(paste('+/+\nn=', length(longuncorr), sep=''),
                     paste('+/-\nn=', length(longcorr), sep=''),
                     paste('-/+\nn=', length(shortuncorr), sep=''),
                     paste('-/-\nn=', length(shortcorr), sep='')), cex.axis=0.9)
axis(2, at=log(c(1000, 10000, 100000, 1000000))/log(10),
     labels=c(1, 10, 100, 1000))
vioplot(longuncorr, longcorr, shortuncorr, shortcorr, col='grey', add=TRUE)
title(ylab="Max transcript length (kb)")
title(xlab=expression(Delta > 15 ~ "min /" ~ rho < 0.5))


longcorr <- log(mydelays[mydelays['corr']>0.5 & mydelays['meddelay.x']>15, 'maxmaxLastIntrons']) / log(10)
longuncorr <- log(mydelays[mydelays['corr']<0.5 & mydelays['meddelay.x']>15, 'maxmaxLastIntrons']) / log(10)
shortcorr <- log(mydelays[mydelays['corr']>0.5 & mydelays['meddelay.x']<15, 'maxmaxLastIntrons']) / log(10)
shortuncorr <- log(mydelays[mydelays['corr']<0.5 & mydelays['meddelay.x']<15, 'maxmaxLastIntrons']) / log(10)

plot.new()
par(ps=8, cex=1)
par(mar=c(3, 3, 0, 0)+0.4)
par(mgp=c(2, 1, 0))
plot.window(xlim=c(0.5, 4.5), ylim=c(1.7, 5.7))
axis(1, at=seq(4), labels=c(paste('+/+\nn=', length(longuncorr), sep=''),
                     paste('+/-\nn=', length(longcorr), sep=''),
                     paste('-/+\nn=', length(shortuncorr), sep=''),
                     paste('-/-\nn=', length(shortcorr), sep='')), cex.axis=0.9)
axis(2, at=log(c(100, 1000, 10000, 100000, 500000))/log(10),
       labels=c(0.1, 1, 10, 100, 500))
vioplot(longuncorr, longcorr, shortuncorr, shortcorr, col='grey', add=TRUE)
title(ylab="Max last intron length (kb)")
title(xlab=expression(Delta > 15 ~ "min /" ~ rho < 0.5))



longcorr <- mydelays[mydelays['corr']>0.5 & mydelays['meddelay.x']>15, 'lastProportion']
longuncorr <- mydelays[mydelays['corr']<0.5 & mydelays['meddelay.x']>15, 'lastProportion']
shortcorr <- mydelays[mydelays['corr']>0.5 & mydelays['meddelay.x']<15, 'lastProportion']
shortuncorr <- mydelays[mydelays['corr']<0.5 & mydelays['meddelay.x']<15, 'lastProportion']

plot.new()
par(ps=8, cex=1)
par(mar=c(3, 3, 0, 0)+0.4)
par(mgp=c(2, 1, 0))
plot.window(xlim=c(0.5, 4.5), ylim=c(-0.02, 1.02))
axis(1, at=seq(4), labels=c(paste('+/+\nn=', length(longuncorr), sep=''),
                     paste('+/-\nn=', length(longcorr), sep=''),
                     paste('-/+\nn=', length(shortuncorr), sep=''),
                     paste('-/-\nn=', length(shortcorr), sep='')), cex.axis=0.9)
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
vioplot(longuncorr, longcorr, shortuncorr, shortcorr, col='grey', add=TRUE)
title(ylab="Last intron fraction")
title(xlab=expression(Delta > 15 ~ "min /" ~ rho < 0.5))


dev.off()



library(plotrix)
J <- (delays[,'tmax.x'] < 160) & (delays[,'tmax.x'] > 1) & (delays[,'begdev10.x'] < 12)

delays.clipped <- delays[J,'meddelay.x']
delays.clipped[delays.clipped > 120] <- 121
h <- hist(delays.clipped, breaks=c(seq(0,130,by=10)), plot=FALSE)
h$counts[1] <- h$counts[1]-1500

par(mfrow=c(1, 1))
plot(h, axes=FALSE, main="Posterior median delays", xlab="Delay (min)", ylab="# of genes")
axis(1, at=c(seq(0, 110, by=40), 125), labels=c(seq(0, 110, by=40), ">120"))
axis.break(2,120,style="zigzag")
axis(2, at=c(0, 50, 100, 150, 200), labels=c(0, 50, 100, 1650, 1700))

pdf('delay_histogram_med.pdf', width=87/25.4, height=70/25.4)
par(ps=8, cex=1)
par(mar=c(2, 2, 1, 0)+0.4)
par(mgp=c(1.5, 0.5, 0))
par(mfrow=c(1, 1))
plot(h, axes=FALSE, main="Posterior median delays", xlab="Delay (min)", ylab="# of genes")
axis(1, at=c(seq(0, 110, by=40), 125), labels=c(seq(0, 110, by=40), ">120"))
axis.break(2,120,style="zigzag")
axis(2, at=c(0, 50, 100, 150, 200), labels=c(0, 50, 100, 1650, 1700))
dev.off()


## plot(h2, axes=FALSE, main="", xlab="Delay (min)", ylab="# of genes")
## axis(1, at=c(seq(0, 110, by=40), 125), labels=c(seq(0, 110, by=40), ">120"))
## axis.break(2,120,style="zigzag")
## axis(2, at=c(0, 50, 100, 150, 200), labels=c(0, 50, 100, 1650, 1700))

## pdf('delay_histogram_premrna.pdf', width=87/25.4, height=70/25.4)
## par(ps=8, cex=1)
## par(mar=c(2, 2, 0, 0)+0.4)
## par(mgp=c(1.5, 0.5, 0))
## plot(h2, axes=FALSE, main="", xlab="Delay (min)", ylab="# of genes")
## axis(1, at=c(seq(0, 110, by=40), 125), labels=c(seq(0, 110, by=40), ">120"))
## axis.break(2,120,style="zigzag")
## axis(2, at=c(0, 50, 100, 150, 200), labels=c(0, 50, 100, 1650, 1700))
## dev.off()

## write.table(delays[,c('maxTrLengths', 'lastProportion')], file='gene_structures.txt', sep='\t', quote=FALSE)

t <- 0:60
v <- t
pvals <- t
for (i in seq_along(t)) {
  v[i] <- mean(mydelays[mydelays["meddelay.x"]>t[i],"premrna_trend"])
  if (i > 1)
    pvals[i] <- wilcox.test(mydelays[mydelays['meddelay.x'] < t[i], 'premrna_trend'], mydelays[mydelays['meddelay.x'] > t[i], 'premrna_trend'])$p.value
}

MAXVAL <- 0.03
NORM <- 3
plot(t, v, xlab="Delay lower bound (min)", ylab="Mean pre-mRNA end accumulation index", type='l', col='blue')
lines(t, rep(-log(0.05)/log(10)/NORM*MAXVAL, length(t)), col='black', lty=2)
axis(4, seq(0, MAXVAL, len=4), seq(0, NORM, by=1))
lines(t, -log(pvals)/log(10)/NORM*MAXVAL, col='black')
mtext(expression(-log[10](p-value)), side=4, line=1.2)

pdf('premrna_halfdiff.pdf', width=87/25.4, height=70/25.4)
par(ps=8, cex=1)
par(mar=c(2, 2, 0, 2)+0.4)
par(mgp=c(1.2, 0.4, 0))
##par(mar=c(2, 2, 1, 0)+0.4)
##par(mgp=c(1.5, 0.5, 0))
par(mfrow=c(1, 1))
plot(t, v, xlab="Delay lower bound (min)", ylab="Mean pre-mRNA end accumulation index", type='l')
MAXVAL <- 0.03
NORM <- 3
plot(t, v, xlab="Delay lower bound (min)", ylab="Mean pre-mRNA end accumulation index", type='l', col='blue')
lines(t, rep(-log(0.05)/log(10)/NORM*MAXVAL, length(t)), col='black', lty=2)
axis(4, seq(0, MAXVAL, len=4), seq(0, NORM, by=1))
lines(t, -log(pvals)/log(10)/NORM*MAXVAL, col='black')
mtext(expression(-log[10](p-value)), side=4, line=1.2)
dev.off()
