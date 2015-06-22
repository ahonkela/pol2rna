# FIGURE PARAMETERS
FONTSIZE=6
DEVBOUND <- 0.05
IQRBOUND <- 120

HISTWIDTH=62/25.4
HISTHEIGHT=45/25.4

cleanup_merge <- function(tbl) {
  row.names(tbl) <- tbl[,'Row.names']
  tbl[!names(tbl) %in% c('Row.names')]
}


delays.pol2 <- read.table('../matlab/results/pol2max_and_meddelays_final.txt', row.names=1, header=TRUE)
quantiles <- read.table('../matlab/results/hmc_results_to_browser_final.txt', row.names=1, header=TRUE)
delays.orig <- cleanup_merge(merge(delays.pol2, quantiles, by=0))
stopifnot(all(delays.orig[,'meddelay'] == delays.orig[,'X50.']))

premrna.fits0 <- read.table('../python/premrna_halfdiff_2014-06-18.txt', row.names=1, header=FALSE)
names(premrna.fits0) <- 'premrna_trend'
pol2.fits <- read.table('../python/pol2_halfdiff_2015-06-10.txt', row.names=1, header=FALSE)
names(pol2.fits) <- 'pol2_trend'
pol2.5pct <- read.table('../python/pol2_last5pct_2015-06-10.txt', row.names=1, header=FALSE)
pol2.5pct.b <- read.table('../python/pol2_last5pct_b_2015-06-10.txt', row.names=1, header=FALSE)
pol2.5pct.c <- read.table('../python/pol2_last5pct_c_2015-06-10.txt', row.names=1, header=FALSE)
premrna.fits1 <- cleanup_merge(merge(premrna.fits0, pol2.fits, by=0))
premrna.fits <- cleanup_merge(merge(premrna.fits1, pol2.5pct, by=0))
premrna.fits.b <- cleanup_merge(merge(premrna.fits1, pol2.5pct.b, by=0))
premrna.fits.c <- cleanup_merge(merge(premrna.fits1, pol2.5pct.c, by=0))
##delays <- read.table('pol2max_and_delays_2013-03-11.txt', row.names=1, header=TRUE)
##delays.orig <- delays.pol2
delays2 <- cleanup_merge(merge(delays.orig, premrna.fits, by=0))
delays3 <- cleanup_merge(merge(delays.orig, premrna.fits.b, by=0))
delays4 <- cleanup_merge(merge(delays.orig, premrna.fits.c, by=0))
##delays.orig <- delays2

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

maxTrLengthID <- sapply(split(lengths, substr(names(lengths), 1, 15)), which.max)
##maxTrLengthID <- maxTrLengthID[names(maxTrLengthID) %in% row.names(delays.orig)]

maxTrLastIntrons <- lastIntronsTr[substr(names(maxTrLengthID), 17, 100)]
names(maxTrLastIntrons) <- substr(names(maxTrLastIntrons), 1, 15)
maxTrLastIntrons <- maxTrLastIntrons[names(maxTrLastIntrons) %in% row.names(delays.orig)]

maxExon3Lengths <- sapply(split(exonlen3, substr(names(exonlen3), 1, 15)), max)
maxExon3Lengths <- maxExon3Lengths[names(maxExon3Lengths) %in% row.names(delays.orig)]

maxExon5Lengths <- sapply(split(exonlen5, substr(names(exonlen5), 1, 15)), max)
maxExon5Lengths <- maxExon5Lengths[names(maxExon5Lengths) %in% row.names(delays.orig)]

##lastProportion <- (maxmaxLastIntrons / maxmaxIntrons)
lastProportion <- (maxTrLastIntrons / maxTrLengths)

exskips <- read.csv('skipped_exons.csv', header=FALSE, colClasses=c(NA, 'logical'))
exonskips <- exskips[,2]
names(exonskips) <- exskips[,1]
exonskips <- exonskips[names(maxmaxIntrons)]

delays <- cleanup_merge(merge(delays.orig, cbind(maxmaxLastIntrons, maxmaxIntrons, maxTrLengths, lastProportion, maxExon3Lengths, maxExon5Lengths, exonskips), by=0))
I <- (delays[,'tmax'] < 160) & (delays[,'tmax'] > 1) & (delays[,'begdev10'] < DEVBOUND) & (delays[,'meddelay'] < 120) & ((delays[,'X75.'] - delays[,'X25.']) < IQRBOUND) #& (delays[,'corr'] > 0.5)
mydelays <- delays[I,]


plot_delay_survival <- function(mydelays, key, lenco, delaykey='meddelay', MAXVAL=0.25, NORM=10) {
  shortlast <- mydelays[mydelays[key]<lenco,delaykey]
  longlast <- mydelays[mydelays[key]>lenco,delaykey]

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
  }
  plot(T, shortfreq, type='l', col='blue', ylim=c(0, MAXVAL), ylab=NA, xlab=NA, axes=FALSE)
  lines(T, longfreq, col='red')
  lines(T, rep(-log(0.05)/log(10)/NORM*MAXVAL, length(T)), col='black', lty=2)
  axis(side=1, labels = NA)
  axis(side=1, lwd = 0, line = -.3)
  mtext("t (min)", side=1, line=0.3)
  axis(side=2, labels = NA)
  axis(side=2, lwd = 0, line = 0)
  mtext(expression("Fraction of genes with" ~ Delta > t), side=2, line=0.5)
  axis(4, seq(0, MAXVAL, len=(NORM/2)+1), seq(0, NORM, by=2), mgp=c(-0.6, -0.2, 0))
  lines(T, -log(pvals)/log(10)/NORM*MAXVAL, col='black')
  mtext(expression(-log[10](p-value)), side=4, line=0.3)
  return (c(short=length(shortlast), long=length(longlast)))
}



pdf('delay_survival.pdf', width=87/25.4, height=50/25.4)
par(ps=FONTSIZE, cex=1)
par(mar=c(1.0, 0.8, 0, 0.8)+0.4)
par(mgp=c(0.6, 0.1, 0))
par(mfrow=c(1, 2))
par(tck=-0.03)

I <- (delays[,'tmax'] < 160) & (delays[,'tmax'] > 1) & (delays[,'begdev10'] < DEVBOUND) & (delays[,'meddelay'] < 120) & ((delays[,'X75.'] - delays[,'X25.']) < IQRBOUND) #& (delays[,'corr'] > 0.5)
mydelays <- delays[I,]

lenco <- 1e4
Nlast <- plot_delay_survival(mydelays, 'maxTrLengths', lenco, MAXVAL=0.3, NORM=12)
leg2 <- legend(x=c(30, 83.2), y=c(0.18, 0.31),
       legend=c(sprintf("m<%.0f\n(N=%d)", lenco, Nlast['short']),
         sprintf("m>%.0f\n(N=%d)", lenco, Nlast['long']),
         'p-value'),
       col=c('blue', 'red', 'black'), lty=1,
       x.intersp=0.5, y.intersp=1, seg.len=1)

lenco <- 0.2
Nlast <- plot_delay_survival(mydelays, 'lastProportion', lenco)
leg1 <- legend(x=c(30, 83.2), y=c(0.15, 0.26),
       legend=c(sprintf("f<%.2f\n(N=%d)", lenco, Nlast['short']),
         sprintf("f>%.2f\n(N=%d)", lenco, Nlast['long']),
         'p-value'),
       col=c('blue', 'red', 'black'), lty=1,
       x.intersp=0.5, y.intersp=1, seg.len=1, xjust=0.5, yjust=0.5)
dev.off()


plot_delay_survival_titled <- function(mydelays, key, lenco, delaykey='meddelay', MAXVAL=0.25, NORM=10, title="") {
  shortlast <- mydelays[mydelays[key]<lenco,delaykey]
  longlast <- mydelays[mydelays[key]>lenco,delaykey]

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
  }
  plot(T, shortfreq, type='l', col='blue', ylim=c(0, MAXVAL), ylab=NA, xlab=NA, axes=FALSE, main=title)
  lines(T, longfreq, col='red')
  lines(T, rep(-log(0.05)/log(10)/NORM*MAXVAL, length(T)), col='black', lty=2)
  axis(side=1, labels = NA)
  axis(side=1, lwd = 0, line = -.2)
  mtext("t (min)", side=1, line=0.5)
  axis(side=2, labels = NA)
  axis(side=2, lwd = 0, line = 0.1)
  mtext(expression("Fraction of genes with" ~ Delta > t), side=2, line=0.7)
  axis(4, seq(0, MAXVAL, len=(NORM/2)+1), seq(0, NORM, by=2), mgp=c(-0.6, -0.1, 0))
  lines(T, -log(pvals)/log(10)/NORM*MAXVAL, col='black')
  mtext(expression(-log[10](p-value)), side=4, line=0.7)
  return (c(short=length(shortlast), long=length(longlast)))
}



pdf('delay_survival_lenfiltered.pdf', width=130/25.4, height=50/25.4)
par(ps=FONTSIZE, cex=1)
par(mar=c(1.3, 1.3, 0.7, 1.3)+0.4)
par(mgp=c(0.6, 0.1, 0))
par(mfrow=c(1, 3))
par(tck=-0.03)

LBOUNDS <- c(-1, 10000, 30000)
TITLES <- c("All genes", "m > 10 kb", "m > 30 kb")

for (k in seq_along(LBOUNDS)) {
I <- (delays[,'tmax'] < 160) & (delays[,'tmax'] > 1) & (delays[,'begdev10'] < DEVBOUND) & (delays[,'meddelay'] < 120) & (delays[,'maxTrLengths'] > LBOUNDS[k]) & ((delays[,'X75.'] - delays[,'X25.']) < IQRBOUND) #& (delays[,'corr'] > 0.5)
mydelays <- delays[I,]

lenco <- 0.2
Nlast <- plot_delay_survival_titled(mydelays, 'lastProportion', lenco, title=TITLES[k])
leg1 <- legend(x=c(30, 83.2), y=c(0.15, 0.26),
       legend=c(sprintf("f<%.2f\n(N=%d)", lenco, Nlast['short']),
         sprintf("f>%.2f\n(N=%d)", lenco, Nlast['long']),
         'p-value'),
       col=c('blue', 'red', 'black'), lty=1,
       x.intersp=0.5, y.intersp=1, seg.len=1, xjust=0.5, yjust=0.5)
}
dev.off()


pdf('delay_survival_devbounds.pdf', width=87/25.4, height=150/25.4)
par(ps=FONTSIZE, cex=1)
par(mar=c(1.3, 1.3, 0.7, 1.3)+0.4)
par(mgp=c(0.6, 0.1, 0))
par(mfrow=c(3, 2))
par(tck=-0.03)

for (DBOUND in c(0.05, 0.1, 0.01)) {
I <- (delays[,'tmax'] < 160) & (delays[,'tmax'] > 1) & (delays[,'begdev10'] < DBOUND) & (delays[,'meddelay'] < 120) & ((delays[,'X75.'] - delays[,'X25.']) < IQRBOUND) #& (delays[,'corr'] > 0.5)
mydelays <- delays[I,]

mybound <- sprintf('%.2f', DBOUND)
lenco <- 1e4
Nlast <- plot_delay_survival_titled(mydelays, 'maxTrLengths', lenco, MAXVAL=0.3, NORM=12, title=bquote(D < .(mybound)))
leg2 <- legend(x=c(30, 83.2), y=c(0.18, 0.31),
       legend=c(sprintf("m<%.0f\n(N=%d)", lenco, Nlast['short']),
         sprintf("m>%.0f\n(N=%d)", lenco, Nlast['long']),
         'p-value'),
       col=c('blue', 'red', 'black'), lty=1,
       x.intersp=0.5, y.intersp=1, seg.len=1)

lenco <- 0.2
Nlast <- plot_delay_survival_titled(mydelays, 'lastProportion', lenco, title=bquote(D < .(mybound)))
leg1 <- legend(x=c(30, 83.2), y=c(0.15, 0.26),
       legend=c(sprintf("f<%.2f\n(N=%d)", lenco, Nlast['short']),
         sprintf("f>%.2f\n(N=%d)", lenco, Nlast['long']),
         'p-value'),
       col=c('blue', 'red', 'black'), lty=1,
       x.intersp=0.5, y.intersp=1, seg.len=1, xjust=0.5, yjust=0.5)
}
dev.off()


I <- (delays[,'tmax'] < 160) & (delays[,'tmax'] > 1) & (delays[,'begdev10'] < DEVBOUND) & (delays[,'meddelay'] < 120) & ((delays[,'X75.'] - delays[,'X25.']) < IQRBOUND) #& (delays[,'corr'] > 0.5)
mydelays <- delays[I,]


pdf('delay_survival_more.pdf', width=87/25.4, height=50/25.4)
par(ps=FONTSIZE, cex=1)
par(mar=c(1.0, 0.8, 0, 0.8)+0.4)
par(mgp=c(0.6, 0.1, 0))
par(mfrow=c(1, 2))
par(tck=-0.03)

mvals <- c(30000, 50000)
fvals <- c(0.3, 0.5)

for (i in seq_along(mvals)) {
  lenco <- mvals[i]
  Nlast <- plot_delay_survival(mydelays, 'maxTrLengths', lenco, MAXVAL=0.3, NORM=12)
  leg2 <- legend(x=c(30, 83.2), y=c(0.18, 0.31),
                 legend=c(sprintf("m<%.0f\n(N=%d)", lenco, Nlast['short']),
                   sprintf("m>%.0f\n(N=%d)", lenco, Nlast['long']),
                   'p-value'),
                 col=c('blue', 'red', 'black'), lty=1,
                 x.intersp=0.5, y.intersp=1, seg.len=1)

  lenco <- fvals[i]
  Nlast <- plot_delay_survival(mydelays, 'lastProportion', lenco)
  leg1 <- legend(x=c(30, 83.2), y=c(0.15, 0.26),
                 legend=c(sprintf("f<%.2f\n(N=%d)", lenco, Nlast['short']),
                   sprintf("f>%.2f\n(N=%d)", lenco, Nlast['long']),
                   'p-value'),
                 col=c('blue', 'red', 'black'), lty=1,
                 x.intersp=0.5, y.intersp=1, seg.len=1, xjust=0.5, yjust=0.5)
}
dev.off()



pdf('delay_survival_exonskip.pdf', width=87/25.4, height=50/25.4)
par(ps=FONTSIZE, cex=1)
par(mar=c(1.0, 0.8, 0, 0.8)+0.4)
par(mgp=c(0.6, 0.1, 0))
par(mfrow=c(1, 1))
par(tck=-0.03)

Nlast <- plot_delay_survival(mydelays, 'exonskips', 0.5)
leg2 <- legend(x=c(30, 83.2), y=c(0.15, 0.26),
       legend=c(sprintf("no skips\n(N=%d)", Nlast['short']),
         sprintf("skips\n(N=%d)", Nlast['long']),
         'p-value'),
       col=c('blue', 'red', 'black'), lty=1,
       x.intersp=0.5, y.intersp=1, seg.len=1)

dev.off()

write.table(row.names(mydelays), file='analysed_genes.txt', quote=FALSE, row.names=FALSE, col.names=FALSE)

library(plotrix)
J <- (delays[,'tmax'] < 160) & (delays[,'tmax'] > 1) & (delays[,'begdev10'] < DEVBOUND) & ((delays[,'X75.'] - delays[,'X25.']) < IQRBOUND)

delays.clipped <- delays[J,'meddelay']
delays.clipped[delays.clipped > 120] <- 121
h.med <- hist(delays.clipped, breaks=c(seq(0,130,by=10)), plot=FALSE)
h.med$counts[1] <- h.med$counts[1]-1350

h.p25 <- hist(mydelays[,'X25.'], breaks=c(seq(0,120,by=10)), plot=FALSE)
h.p25$counts[1] <- h.p25$counts[1]-1575

par(mfrow=c(1, 1))
plot(h.med, axes=FALSE, main="", xlab="Delay (min)", ylab="# of genes")
axis(1, at=c(seq(0, 110, by=40), 125), labels=c(seq(0, 110, by=40), ">120"))
axis.break(2,120,style="zigzag")
axis(2, at=c(0, 50, 100, 150, 200), labels=c(0, 50, 100, 1500, 1550))

pdf('delay_histogram_med.pdf', width=HISTWIDTH, height=HISTHEIGHT)
par(ps=FONTSIZE, cex=1)
par(mar=c(1.5, 1.2, 0, 0)+0.4)
par(mgp=c(1, 0.4, 0))
par(mfrow=c(1, 1))
plot(h.med, axes=FALSE, main="", xlab="Delay (min)", ylab="# of genes")
axis(1, at=c(seq(0, 110, by=40), 125), labels=c(seq(0, 110, by=40), ">120"))
axis.break(2,120,style="zigzag")
axis(2, at=c(0, 50, 100, 150, 200), labels=c(0, 50, 100, 1500, 1550))
## plot(h.p25, axes=FALSE, main="", xlab="Delay (min)", ylab="# of genes", ylim=c(0, 100))
## axis(1, at=seq(0, 120, by=40), labels=seq(0, 120, by=40))
## axis.break(2,60,style="zigzag")
## axis(2, at=c(0, 25, 50, 75, 100), labels=c(0, 25, 50, 75+1575, 100+1575))
dev.off()


I2 <- (delays2[,'tmax'] < 160) & (delays2[,'tmax'] > 1) & (delays2[,'begdev10'] < DEVBOUND) & (delays2[,'meddelay'] < 120) & ((delays2[,'X75.'] - delays2[,'X25.']) < IQRBOUND) #& (delays2[,'corr'] > 0.5)
mydelays2 <- delays2[I2,]
mydelays3 <- delays3[I2,]
mydelays4 <- delays4[I2,]


plot_halfdiff <- function(mydelays2, key, response, ylab) {
  par(ps=FONTSIZE, cex=1)
  par(mar=c(1.0, 0.8, 0, 0.8)+0.4)
  par(mgp=c(0.6, 0.1, 0))
  par(tck=-0.015)

  t <- 5:60
  v <- t
  v2 <- v
  v2[1] <- NA
  pvals <- t
  for (i in seq_along(t)) {
    v[i] <- mean(mydelays2[mydelays2[key]>t[i],response])
    v2[i] <- mean(mydelays2[mydelays2[key]<t[i],response])
    pvals[i] <- wilcox.test(mydelays2[mydelays2[key] < t[i], response], mydelays2[mydelays2[key] > t[i], response])$p.value
  }
  ## Pol-II data has a different sign
  if (response == "pol2_trend") {
    v = -v
    v2 = -v2
  }
  
  MAXVAL <- 0.04
  MINVAL <- -0.01
  NORM <- 7
  plot(t, v, xlab=NA, ylab=NA, type='l', col='blue', ylim=c(MINVAL, MAXVAL), axes=FALSE)
  lines(t, v2, type='l', col='red')
  lines(t, rep(-log(0.05)/log(10)/NORM*(MAXVAL-MINVAL)+MINVAL, length(t)), col='black', lty=2)
  lines(t, -log(pvals)/log(10)/NORM*(MAXVAL-MINVAL)+MINVAL, col='black')
  axis(side=1, labels = NA)
  axis(side=1, lwd = 0, line = -.3)
  mtext("t (min)", side=1, line=0.3)
  axis(side=2, labels = NA)
  axis(side=2, lwd = 0, line = 0)
  mtext(ylab, side=2, line=0.5)
  axis(4, seq(MINVAL, MAXVAL, len=(1+NORM)), seq(0, NORM, by=1), mgp=c(-0.6, -0.2, 0))
  mtext(expression(-log[10](p-value)), side=4, line=0.3)
  if (response == "pol2_trend") {
    legend('right',
           legend=c(expression(Delta > t), expression(Delta < t), 'p-value'),
           col=c('blue', 'red', 'black'), lty=1, x.intersp=0.5, y.intersp=0.5,
           seg.len=1, inset=0.01)
  } else {
    legend('topleft',
           legend=c(expression(Delta > t), expression(Delta < t), 'p-value'),
           col=c('blue', 'red', 'black'), lty=1, x.intersp=0.5, y.intersp=0.5,
           seg.len=1, inset=0.01)
  }
}


par(mfrow=c(1, 1))
plot_halfdiff(mydelays2, "meddelay", "premrna_trend", "Mean pre-mRNA end accumulation index")

pdf('premrna_halfdiff.pdf', width=43/25.4, height=50/25.4)
par(mfrow=c(1, 1))
plot_halfdiff(mydelays2, "meddelay", "premrna_trend", "Mean pre-mRNA end accumulation index")
dev.off()


pdf('pol2_halfdiff.pdf', width=87/25.4, height=70/25.4)
par(mfrow=c(1, 1))
plot_halfdiff(mydelays2, "meddelay", "pol2_trend", "Mean Pol-II end accumulation index")
dev.off()


pdf('delay_comparison.pdf', width=50/25.4, height=50/25.4)
par(mfrow=c(1, 1))
par(ps=FONTSIZE, cex=1)
par(mar=c(1.0, 0.8, 0, 0.8)+0.4)
par(mgp=c(0.6, 0.1, 0))
par(tck=-0.015)
smoothScatter(log(mydelays[,'meddelay'])/log(10),
              log(mydelays[,'maxTrLengths']/2.1/1000)/log(10),
              axes=FALSE,
              xlab=expression("RNA processing delay" ~ Delta ~ "(min)"),
              ylab="Transcriptional delay (min)")
axis(1, at=c(0, log(3)/log(10), 1, log(30)/log(10), 2),
     labels=c(1, 3, 10, 30, 100), mgp=c(-0.6, -0.2, 0))
axis(2, at=c(0, log(3)/log(10), 1, log(30)/log(10), 2),
     labels=c(1, 3, 10, 30, 100))
lines(c(-1, 3), c(-1, 3))
dev.off()


pdf('pol2_last5pct.pdf', width=150/25.4, height=50/25.4)
par(mfrow=c(1, 3))
par(ps=FONTSIZE, cex=1)
par(mar=c(1.0, 0.8, 0, 0.8)+0.4)
par(mgp=c(0.6, 0.1, 0))
par(tck=-0.015)
for (cutoff in c(10, 20, 30)) {
  plot(apply(mydelays2[mydelays2[,'meddelay'] > cutoff, 14:23], 2, mean), type='l', xlab='time point', ylab='pol-II in the last 5% of the gene', main=sprintf('t=%d min', cutoff), col='blue', ylim=c(0.047, 0.061))
  lines(apply(mydelays2[mydelays2[,'meddelay'] < cutoff, 14:23], 2, mean), col='red')
  legend('topleft',
         legend=c(expression(Delta > t), expression(Delta < t)), # 'p-value'),
         col=c('blue', 'red'), lty=1, x.intersp=0.5, y.intersp=0.5,
         seg.len=1, inset=0.01)
}
dev.off()

pdf('pol2_last5pct_b.pdf', width=150/25.4, height=50/25.4)
par(mfrow=c(1, 3))
par(ps=FONTSIZE, cex=1)
par(mar=c(1.0, 0.8, 0, 0.8)+0.4)
par(mgp=c(0.6, 0.1, 0))
par(tck=-0.015)
for (cutoff in c(10, 20, 30)) {
  pvals <- rep(0, 10)
  for (i in seq(10)) {
    pvals[i] <- wilcox.test(mydelays3[mydelays3[,'meddelay'] < cutoff, 13+i], mydelays3[mydelays3[,'meddelay'] > cutoff, 13+i])$p.value
  }
  plot(apply(mydelays3[mydelays3[,'meddelay'] > cutoff, 14:23], 2, mean), type='l', xlab='time point', ylab='pol-II in the last 5% of the gene', main=sprintf('t=%d min', cutoff), col='blue', ylim=c(0.115, 0.131))
  lines(apply(mydelays3[mydelays3[,'meddelay'] < cutoff, 14:23], 2, mean), col='red')
  cat(pvals)
  legend('topleft',
         legend=c(expression(Delta > t), expression(Delta < t)), # 'p-value'),
         col=c('blue', 'red'), lty=1, x.intersp=0.5, y.intersp=0.5,
         seg.len=1, inset=0.01)
}
dev.off()

pdf('pol2_last5pct_c.pdf', width=150/25.4, height=50/25.4)
par(mfrow=c(1, 3))
par(ps=FONTSIZE, cex=1)
par(mar=c(1.0, 0.8, 0, 0.8)+0.4)
par(mgp=c(0.6, 0.1, 0))
par(tck=-0.015)
for (cutoff in c(10, 20, 30)) {
  plot(apply(mydelays4[mydelays4[,'meddelay'] > cutoff, 14:23], 2, mean), type='l', xlab='time point', ylab='pol-II in the last 5% of the gene', main=sprintf('t=%d min', cutoff), col='blue')
  lines(apply(mydelays4[mydelays4[,'meddelay'] < cutoff, 14:23], 2, mean), col='red')
  legend('topleft',
         legend=c(expression(Delta > t), expression(Delta < t)), # 'p-value'),
         col=c('blue', 'red'), lty=1, x.intersp=0.5, y.intersp=0.5,
         seg.len=1, inset=0.01)
}
dev.off()
