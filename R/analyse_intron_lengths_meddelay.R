delays.pol2 <- read.table('pol2max_and_meddelays_2013-08-30.txt', row.names=1, header=TRUE)
delays.premrna <- read.table('pol2max_and_meddelays_2013-11-05.txt', row.names=1, header=TRUE)
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

delays.clipped <- delays[J,'meddelay.y']
delays.clipped[delays.clipped > 120] <- 121
h2 <- hist(delays.clipped, breaks=c(seq(0,130,by=10)), plot=FALSE)
h2$counts[1] <- h2$counts[1]-1500

par(mfrow=c(1, 2))
plot(h, axes=FALSE, main="Pol II - mRNA delays", xlab="Delay (min)", ylab="# of genes")
axis(1, at=c(seq(0, 110, by=40), 125), labels=c(seq(0, 110, by=40), ">120"))
axis.break(2,120,style="zigzag")
axis(2, at=c(0, 50, 100, 150, 200), labels=c(0, 50, 100, 1650, 1700))
plot(h2, axes=FALSE, main="pre-mRNA - mRNA delays", xlab="Delay (min)", ylab="# of genes")
axis(1, at=c(seq(0, 110, by=40), 125), labels=c(seq(0, 110, by=40), ">120"))
axis.break(2,120,style="zigzag")
axis(2, at=c(0, 50, 100, 150, 200), labels=c(0, 50, 100, 1650, 1700))

pdf('delay_histogram_med.pdf', width=178/25.4, height=70/25.4)
par(ps=8, cex=1)
par(mar=c(2, 2, 1, 0)+0.4)
par(mgp=c(1.5, 0.5, 0))
par(mfrow=c(1, 2))
plot(h, axes=FALSE, main="Pol II - mRNA delays", xlab="Delay (min)", ylab="# of genes")
axis(1, at=c(seq(0, 110, by=40), 125), labels=c(seq(0, 110, by=40), ">120"))
axis.break(2,120,style="zigzag")
axis(2, at=c(0, 50, 100, 150, 200), labels=c(0, 50, 100, 1650, 1700))
plot(h2, axes=FALSE, main="pre-mRNA - mRNA delays", xlab="Delay (min)", ylab="# of genes")
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



##library(ggplot2)
##qplot(factor(myfact), lastProportion, data=mydelays, geom = "violin")

## ##glm.longest <- glm(formula=longestLast ~ corr + delay, family=binomial)
## glm.longest <- glm(formula=longestLast ~ corr + delay + corr*delay, family=binomial)
mydelays2 = delays[I2,]
fields <- c('lastProportion', 'maxmaxLastIntrons', 'maxTrLengths', 'maxExon3Lengths', 'maxExon5Lengths', 'maxmaxIntrons')
for (k in fields) {
  mydelays2[,k] <- mydelays2[,k] - mean(mydelays2[,k])
  mydelays2[,k] <- mydelays2[,k] / (2*sd(mydelays2[,k]))
}
##mydelays2[,'meddelay.x'] = mydelays2[,'meddelay.x']/60
##mydelays2[,'maxmaxLastIntrons'] = mydelays2[,'maxmaxLastIntrons'] * sd(mydelays2[,'lastProportion']) / sd(mydelays2[,'maxmaxLastIntrons'])
##mydelays2[,'maxmaxLastIntrons'] = mydelays2[,'maxmaxLastIntrons'] * sd(mydelays2[,'lastProportion']) / sd(mydelays2[,'maxmaxLastIntrons'])
##mydelays2[,'maxmaxLastIntrons'] = mydelays2[,'maxmaxLastIntrons']/1000

## lm.longest <- lm(formula=lastProportion ~ corr + meddelay.x + corr*meddelay.x, data=mydelays2)
## summary(lm.longest)
## anova(lm.longest)

require(arm)
##blm.length <- bayesglm(formula=maxmaxLastIntrons ~ corr + meddelay.x + corr*meddelay.x,
##                       data=mydelays2)
##summary(blm.length)
##anova(blm.length)

## blm.longest <- bayesglm(formula=lastProportion ~ corr + meddelay.x + corr*meddelay.x,
##                         data=mydelays2)
## blm.longest2 <- bayesglm(formula=lastProportion ~ corr + meddelay.x,
##                          data=mydelays2)
## blm.longest3 <- bayesglm(formula=lastProportion ~ corr,
##                          data=mydelays2)
## blm.longest4 <- bayesglm(formula=lastProportion ~ meddelay.x,
##                          data=mydelays2)

blm.delay <- bayesglm(formula=meddelay.x ~ lastProportion + maxmaxLastIntrons + corr,
                      data=mydelays2)
blm.delay2 <- bayesglm(formula=meddelay.x ~ lastProportion + maxmaxLastIntrons,
                      data=mydelays2)
blm.delay3 <- bayesglm(formula=meddelay.x ~ lastProportion,
                      data=mydelays2)
blm.delay4 <- bayesglm(formula=meddelay.x ~ maxmaxLastIntrons,
                      data=mydelays2)
blm.delay5 <- bayesglm(formula=meddelay.x ~ lastProportion + maxmaxLastIntrons +
                       lastProportion*maxmaxLastIntrons,
                       data=mydelays2)
blm.delay5b <- bayesglm(formula=meddelay.x ~ lastProportion + maxTrLengths +
                       lastProportion*maxTrLengths,
                       data=mydelays2)
blm.delay6 <- bayesglm(formula=meddelay.x ~ lastProportion + maxmaxLastIntrons +
                       lastProportion*maxmaxLastIntrons + corr,
                       data=mydelays2)
blm.delay7 <- bayesglm(formula=meddelay.x ~ lastProportion + maxmaxLastIntrons +
                       maxTrLengths,
                       data=mydelays2)

blm.delay8 <- bayesglm(formula=meddelay.x ~ lastProportion + maxmaxLastIntrons +
                       maxTrLengths + maxExon3Lengths + maxExon5Lengths,
                       data=mydelays2)
blm.delay8a <- bayesglm(formula=meddelay.x ~ maxExon3Lengths + maxExon5Lengths,
                        data=mydelays2)


blm.corr <- bayesglm(formula=corr ~ lastProportion + maxmaxLastIntrons + meddelay.x,
                      data=mydelays2)
blm.corr2 <- bayesglm(formula=corr ~ lastProportion + maxmaxLastIntrons,
                      data=mydelays2)
blm.corr3 <- bayesglm(formula=corr ~ lastProportion,
                      data=mydelays2)
blm.corr4 <- bayesglm(formula=corr ~ maxmaxLastIntrons,
                      data=mydelays2)
blm.corr5 <- bayesglm(formula=corr ~ lastProportion + maxmaxLastIntrons +
                      lastProportion*maxmaxLastIntrons,
                      data=mydelays2)
blm.corr6 <- bayesglm(formula=corr ~ lastProportion + maxmaxLastIntrons +
                      lastProportion*maxmaxLastIntrons + meddelay.x,
                      data=mydelays2)


##coefplot(blm.longest2, varnames=expression(Delta, rho), main="")
par(mfrow=c(2, 1))
par(ps=8, cex=1)
par(mar=c(0, 0, 0, 0)+0.4)
coefplot(blm.delay2, varnames=expression("1", beta[f], beta[m]), main="Regression of mean delay")
coefplot(blm.corr2, varnames=expression("1", beta[f], beta[m]), main="Regression of Pol II-pre-mRNA correlation")

pdf('regression_coefs_med.pdf', width=87/25.4, height=80/25.4)
par(mfrow=c(2, 1))
par(ps=8, cex=1)
par(mar=c(0, 0, 0, 0)+0.4)
##par(mgp=c(1.5, 0.5, 0))
coefplot(blm.delay2, varnames=expression("1", beta[f], beta[m]), main="Regression of mean delay")
coefplot(blm.corr2, varnames=expression("1", beta[f], beta[m]), main="Regression of Pol II-pre-mRNA correlation")
dev.off()




## blm.length <- bayesglm(formula=maxmaxLastIntrons ~ corr + meddelay.x + corr*meddelay.x,
##                        data=mydelays2)
## blm.length2 <- bayesglm(formula=maxmaxLastIntrons ~ corr + meddelay.x,
##                         data=mydelays2)
## blm.length3 <- bayesglm(formula=maxmaxLastIntrons ~ corr,
##                         data=mydelays2)
## blm.length4 <- bayesglm(formula=maxmaxLastIntrons ~ meddelay.x,
##                         data=mydelays2)
## summary(blm.length)
## summary(blm.length2)
## summary(blm.length3)
## summary(blm.length4)
