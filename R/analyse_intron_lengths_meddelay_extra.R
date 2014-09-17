source('analyse_intron_lengths_meddelay.R')

lencos <- c(0.5, 0.75, 0.9, 0.95)
pdf('corr_survival.pdf', width=87/25.4, height=70/25.4)
par(ps=FONTSIZE, cex=1)
par(mar=c(2, 2, 0, 2)+0.4)
par(mgp=c(1.2, 0.4, 0))
par(mfrow=c(1, 1))
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



library(vioplot)

pdf('intron_lengths_med.pdf', width=178/25.4, height=70/25.4)
par(mfrow=c(1, 3))

longcorr <- log(mydelays[mydelays['corr']>0.5 & mydelays['meddelay.x']>15, 'maxTrLengths']) / log(10)
longuncorr <- log(mydelays[mydelays['corr']<0.5 & mydelays['meddelay.x']>15, 'maxTrLengths']) / log(10)
shortcorr <- log(mydelays[mydelays['corr']>0.5 & mydelays['meddelay.x']<15, 'maxTrLengths']) / log(10)
shortuncorr <- log(mydelays[mydelays['corr']<0.5 & mydelays['meddelay.x']<15, 'maxTrLengths']) / log(10)

plot.new()
par(ps=FONTSIZE, cex=1)
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
par(ps=FONTSIZE, cex=1)
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
par(ps=FONTSIZE, cex=1)
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



pdf('premrna_halfdiff2.pdf', width=87/25.4, height=70/25.4)
par(ps=FONTSIZE, cex=1)
par(mar=c(2, 2, 0, 2)+0.4)
par(mgp=c(1.2, 0.4, 0))
##par(mar=c(2, 2, 1, 0)+0.4)
##par(mgp=c(1.5, 0.5, 0))
par(mfrow=c(1, 1))
for (scale in c(20, 10, 5, 1)) {
  J <- round(mydelays2["meddelay.x"]/scale)
  N <- tabulate(J[,1]+1)
  v <- sapply(split(mydelays2[c("meddelay.x", "premrna_trend")], J), colMeans)
  plot(v[1,], v[2,], xlab="Delay (min)", ylab="Mean pre-mRNA end accumulation index", type='l')
}
dev.off()


pdf('pol2_halfdiff2.pdf', width=87/25.4, height=70/25.4)
par(ps=FONTSIZE, cex=1)
par(mar=c(2, 2, 0, 2)+0.4)
par(mgp=c(1.2, 0.4, 0))
##par(mar=c(2, 2, 1, 0)+0.4)
##par(mgp=c(1.5, 0.5, 0))
par(mfrow=c(1, 1))
for (scale in c(20, 10, 5, 1)) {
  J <- round(mydelays2["meddelay.x"]/scale)
  v <- sapply(split(mydelays2[c("meddelay.x", "pol2_trend")], J), colMeans)
  plot(v[1,], -v[2,], xlab="Delay (min)", ylab="Mean Pol-II end accumulation index", type='l')
}
dev.off()

