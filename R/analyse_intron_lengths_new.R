delays.pol2 <- read.table('pol2max_and_delays_2013-08-30.txt', row.names=1, header=TRUE)
delays.premrna <- read.table('pol2max_and_delays_2013-11-05.txt', row.names=1, header=TRUE)
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

t <- readLines('intron_lengths.txt.lengths')
l <- strsplit(t, '\t')
names(l) <- sapply(l, function(x) x[1])
introns <- lapply(l, function(x) as.integer(x[-1:-4]))
strands <- sapply(l, function(x) x[2])
termpos <- sapply(l, function(x) as.integer(x[3]))
lengths <- sapply(l, function(x) as.integer(x[4]))

## t <- readLines('intron_lengths.txt')
## l <- strsplit(t, '\t')
## names(l) <- sapply(l, function(x) x[1])
## introns <- lapply(l, function(x) as.integer(x[-1]))

numIntrons <- unlist(sapply(introns, length))

introns <- introns[numIntrons > 0]
termpos <- termpos[numIntrons > 0]
strands <- strands[numIntrons > 0]
lengths <- lengths[numIntrons > 0]
numIntrons <- numIntrons[numIntrons > 0]

lastIntronsTr <- unlist(sapply(introns, function(x) x[length(x)]))

maxIntronsTr <- unlist(sapply(introns, function(x) max(x)))
maxmaxIntrons <- sapply(split(maxIntronsTr, substr(names(maxIntronsTr), 1, 15)), max)
maxmaxIntrons <- maxmaxIntrons[names(maxmaxIntrons) %in% row.names(delays.orig)]

maxmaxLastIntrons <- sapply(split(lastIntronsTr, substr(names(lastIntronsTr), 1, 15)), max)

maxmaxLastIntrons <- maxmaxLastIntrons[names(maxmaxLastIntrons) %in% row.names(delays.orig)]

maxTrLengths <- sapply(split(lengths, substr(names(lengths), 1, 15)), max)
maxTrLengths <- maxTrLengths[names(maxTrLengths) %in% row.names(delays.orig)]

longestLast <- (maxmaxLastIntrons == maxmaxIntrons)
lastProportion <- (maxmaxLastIntrons / maxmaxIntrons)
myfact <- rep(0, length(lastProportion))

delays <- merge(delays.orig, cbind(maxmaxLastIntrons, maxmaxIntrons, longestLast, maxTrLengths, lastProportion, myfact), by=0)
row.names(delays) <- delays[,1]
delays <- delays[,-1]
I <- (delays[,'tmax.x'] < 160) & (delays[,'tmax.x'] > 1) & (delays[,'begdev10.x'] < 12) & (delays[,'delay.x'] < 120) #& (delays[,'corr'] > 0.5)
mydelays <- delays[I,]
mydelays[mydelays['corr']<0.5 & mydelays['delay.x']>15, 'myfact'] <- 1
mydelays[mydelays['corr']>0.5 & mydelays['delay.x']>15, 'myfact'] <- 2
mydelays[mydelays['corr']<0.5 & mydelays['delay.x']<15, 'myfact'] <- 3
mydelays[mydelays['corr']>0.5 & mydelays['delay.x']<15, 'myfact'] <- 4
longestLast <- mydelays[,'longestLast']
##lastProportion <- mydelays[,'maxmaxLastIntrons']/mydelays[,'maxTrLengths']
lastProportion <- mydelays[,'lastProportion']/mydelays[,'maxmaxIntrons']
##long <- longestLast[mydelays['delay.x']>15,]
##short <- longestLast[mydelays['delay.x']<15,]
##longcorr <- longestLast[mydelays['corr']>0.5 & mydelays['delay.x']>15]
##longuncorr <- longestLast[mydelays['corr']<0.5 & mydelays['delay.x']>15]
##shortcorr <- longestLast[mydelays['corr']>0.5 & mydelays['delay.x']<15]
##shortuncorr <- longestLast[mydelays['corr']<0.5 & mydelays['delay.x']<15]
longcorr <- mydelays[mydelays['corr']>0.5 & mydelays['delay.x']>15, 'lastProportion']
longuncorr <- mydelays[mydelays['corr']<0.5 & mydelays['delay.x']>15, 'lastProportion']
shortcorr <- mydelays[mydelays['corr']>0.5 & mydelays['delay.x']<15, 'lastProportion']
shortuncorr <- mydelays[mydelays['corr']<0.5 & mydelays['delay.x']<15, 'lastProportion']



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
par(ps=8, cex=1)
par(mar=c(2, 2, 0, 0)+0.4)
par(mgp=c(1.5, 0.5, 0))
vioplot(longuncorr, longcorr, shortuncorr, shortcorr,
        names=c(paste('+/+, n=', length(longuncorr), sep=''),
          paste('+/-, n=', length(longcorr), sep=''),
          paste('-/+, n=', length(shortuncorr), sep=''),
          paste('-/-, n=', length(shortcorr), sep='')), col='grey')
title(ylab="Last intron fraction")
title(xlab=expression(Delta > 15 ~ "min /" ~ rho < 0.5))

pdf('last_intron_fraction.pdf', width=87/25.4, height=70/25.4)
par(ps=8, cex=1)
par(mar=c(2, 2, 0, 0)+0.4)
par(mgp=c(1.5, 0.5, 0))
vioplot(longuncorr, longcorr, shortuncorr, shortcorr,
        names=c(paste('+/+, n=', length(longuncorr), sep=''),
          paste('+/-, n=', length(longcorr), sep=''),
          paste('-/+, n=', length(shortuncorr), sep=''),
          paste('-/-, n=', length(shortcorr), sep='')), col='grey')
title(ylab="Last intron fraction")
title(xlab=expression(Delta > 15 ~ "min /" ~ rho < 0.5))
dev.off()




library(plotrix)
J <- (delays[,'tmax.x'] < 160) & (delays[,'tmax.x'] > 1) & (delays[,'begdev10.x'] < 12)

delays.clipped <- delays[J,'delay.x']
delays.clipped[delays.clipped > 120] <- 121
h <- hist(delays.clipped, breaks=c(seq(0,130,by=10)), plot=FALSE)
h$counts[1] <- h$counts[1]-1500

plot(h, axes=FALSE, main="", xlab="Delay (min)", ylab="# of genes")
axis(1, at=c(seq(0, 110, by=40), 125), labels=c(seq(0, 110, by=40), ">120"))
axis.break(2,120,style="zigzag")
axis(2, at=c(0, 50, 100, 150, 200), labels=c(0, 50, 100, 1650, 1700))

pdf('delay_histogram.pdf', width=87/25.4, height=70/25.4)
par(ps=8, cex=1)
par(mar=c(2, 2, 0, 0)+0.4)
par(mgp=c(1.5, 0.5, 0))
plot(h, axes=FALSE, main="", xlab="Delay (min)", ylab="# of genes")
axis(1, at=c(seq(0, 110, by=40), 125), labels=c(seq(0, 110, by=40), ">120"))
axis.break(2,120,style="zigzag")
axis(2, at=c(0, 50, 100, 150, 200), labels=c(0, 50, 100, 1650, 1700))
dev.off()


delays.clipped <- delays[J,'delay.y']
delays.clipped[delays.clipped > 120] <- 121
h2 <- hist(delays.clipped, breaks=c(seq(0,130,by=10)), plot=FALSE)
h2$counts[1] <- h2$counts[1]-1500

plot(h2, axes=FALSE, main="", xlab="Delay (min)", ylab="# of genes")
axis(1, at=c(seq(0, 110, by=40), 125), labels=c(seq(0, 110, by=40), ">120"))
axis.break(2,120,style="zigzag")
axis(2, at=c(0, 50, 100, 150, 200), labels=c(0, 50, 100, 1650, 1700))

pdf('delay_histogram_premrna.pdf', width=87/25.4, height=70/25.4)
par(ps=8, cex=1)
par(mar=c(2, 2, 0, 0)+0.4)
par(mgp=c(1.5, 0.5, 0))
plot(h2, axes=FALSE, main="", xlab="Delay (min)", ylab="# of genes")
axis(1, at=c(seq(0, 110, by=40), 125), labels=c(seq(0, 110, by=40), ">120"))
axis.break(2,120,style="zigzag")
axis(2, at=c(0, 50, 100, 150, 200), labels=c(0, 50, 100, 1650, 1700))
dev.off()



##library(ggplot2)
##qplot(factor(myfact), lastProportion, data=mydelays, geom = "violin")

## ##glm.longest <- glm(formula=longestLast ~ corr + delay, family=binomial)
## glm.longest <- glm(formula=longestLast ~ corr + delay + corr*delay, family=binomial)
mydelays2 = mydelays
mydelays2[,'delay.x'] = mydelays2[,'delay.x']/60
mydelays2[,'maxmaxLastIntrons'] = mydelays2[,'maxmaxLastIntrons'] * sd(mydelays2[,'lastProportion']) / sd(mydelays2[,'maxmaxLastIntrons'])
##mydelays2[,'maxmaxLastIntrons'] = mydelays2[,'maxmaxLastIntrons']/1000

## lm.longest <- lm(formula=lastProportion ~ corr + delay.x + corr*delay.x, data=mydelays2)
## summary(lm.longest)
## anova(lm.longest)

require(arm)
##blm.length <- bayesglm(formula=maxmaxLastIntrons ~ corr + delay.x + corr*delay.x,
##                       data=mydelays2)
##summary(blm.length)
##anova(blm.length)

## blm.longest <- bayesglm(formula=lastProportion ~ corr + delay.x + corr*delay.x,
##                         data=mydelays2)
## blm.longest2 <- bayesglm(formula=lastProportion ~ corr + delay.x,
##                          data=mydelays2)
## blm.longest3 <- bayesglm(formula=lastProportion ~ corr,
##                          data=mydelays2)
## blm.longest4 <- bayesglm(formula=lastProportion ~ delay.x,
##                          data=mydelays2)

blm.delay <- bayesglm(formula=delay.x ~ lastProportion + maxmaxLastIntrons + corr,
                      data=mydelays2)
blm.delay2 <- bayesglm(formula=delay.x ~ lastProportion + maxmaxLastIntrons,
                      data=mydelays2)
blm.delay3 <- bayesglm(formula=delay.x ~ lastProportion,
                      data=mydelays2)
blm.delay4 <- bayesglm(formula=delay.x ~ maxmaxLastIntrons,
                      data=mydelays2)
blm.delay5 <- bayesglm(formula=delay.x ~ lastProportion + maxmaxLastIntrons +
                       lastProportion*maxmaxLastIntrons,
                       data=mydelays2)
blm.delay6 <- bayesglm(formula=delay.x ~ lastProportion + maxmaxLastIntrons +
                       lastProportion*maxmaxLastIntrons + corr,
                       data=mydelays2)


##coefplot(blm.longest2, varnames=expression(Delta, rho), main="")
par(ps=8, cex=1)
par(mar=c(0, 0, 0, 0)+0.4)
coefplot(blm.longest2, varnames=expression(1, rho, Delta), xlab=expression("Regression coefficient for ", lambda), main="")

pdf('regression_coefs.pdf', width=87/25.4, height=50/25.4)
par(ps=8, cex=1)
par(mar=c(0, 0, 0, 0)+0.4)
##par(mgp=c(1.5, 0.5, 0))
coefplot(blm.longest2, varnames=expression(1, rho, Delta), xlab=expression("Regression coefficient for ", lambda), main="")
dev.off()


blm.corr <- bayesglm(formula=corr ~ lastProportion + maxmaxLastIntrons + delay.x,
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
                      lastProportion*maxmaxLastIntrons + delay.x,
                      data=mydelays2)



## blm.length <- bayesglm(formula=maxmaxLastIntrons ~ corr + delay.x + corr*delay.x,
##                        data=mydelays2)
## blm.length2 <- bayesglm(formula=maxmaxLastIntrons ~ corr + delay.x,
##                         data=mydelays2)
## blm.length3 <- bayesglm(formula=maxmaxLastIntrons ~ corr,
##                         data=mydelays2)
## blm.length4 <- bayesglm(formula=maxmaxLastIntrons ~ delay.x,
##                         data=mydelays2)
## summary(blm.length)
## summary(blm.length2)
## summary(blm.length3)
## summary(blm.length4)
