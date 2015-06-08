# FIGURE PARAMETERS
FONTSIZE=6
DEVBOUND <- 0.05
##LLBOUND <- -100
RESVARBOUND <- 0.2
POL2VARBOUND <- 0.3

HISTWIDTH=62/25.4
HISTHEIGHT=45/25.4

cleanup_merge <- function(tbl) {
  row.names(tbl) <- tbl[,'Row.names']
  tbl[!names(tbl) %in% c('Row.names')]
}


delays.pol2 <- read.table('../matlab/results/pol2max_and_meddelays_final.txt', row.names=1, header=TRUE)
quantiles <- read.table('../matlab/results/hmc_results_to_browser_final.txt', row.names=1, header=TRUE)
ctds <- read.table('../matlab/results/ctd_delays_2015-05-21_spl1.txt', row.names=1, header=TRUE)
#ctds <- read.table('../matlab/results/ctd_delays_2015-05-15.txt', row.names=1, header=TRUE)
delays.orig0 <- cleanup_merge(merge(delays.pol2, quantiles, by=0))
stopifnot(all(delays.orig0[,'meddelay'] == delays.orig0[,'X50.']))
delays.orig <- cleanup_merge(merge(delays.orig0, ctds, by=0))

premrna.fits0 <- read.table('../python/premrna_halfdiff_2014-06-18.txt', row.names=1, header=FALSE)
names(premrna.fits0) <- 'premrna_trend'
pol2.fits <- read.table('../python/pol2_halfdiff_2014-06-11.txt', row.names=1, header=FALSE)
names(pol2.fits) <- 'pol2_trend'
premrna.fits <- cleanup_merge(merge(premrna.fits0, pol2.fits, by=0))
##delays <- read.table('pol2max_and_delays_2013-03-11.txt', row.names=1, header=TRUE)
##delays.orig <- delays.pol2
delays2 <- cleanup_merge(merge(delays.orig, premrna.fits, by=0))
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
I <- (delays[,'tmax'] < 160) & (delays[,'tmax'] > 1) & (delays[,'begdev10'] < DEVBOUND) & (delays[,'ctd_50.'] < 120) & (delays[,'ctd_resvarfrac'] < RESVARBOUND) & (delays[,'spline_var'] < POL2VARBOUND) #& (delays[,'corr'] > 0.5)
mydelays <- delays[I,]


bounds <- seq(-100, 0, by=10)
corrs <- rep(0, length(bounds))
for (k in seq_along(bounds)) {
  myLLBOUND <- bounds[k]
  I <- (delays[,'tmax'] < 160) & (delays[,'tmax'] > 1) & (delays[,'begdev10'] < DEVBOUND) & (delays[,'ctd_50.'] < 120) & (delays[,'avell'] > myLLBOUND) #& (delays[,'corr'] > 0.5)
  mydelays <- delays[I,]
  corrs[k] <- cor((mydelays[,'X50.']), (mydelays[,'ctd_50.']))
}

I <- (delays[,'tmax'] < 160) & (delays[,'tmax'] > 1) & (delays[,'begdev10'] < DEVBOUND) & (delays[,'ctd_50.'] < 120) & (delays[,'spline_var'] < POL2VARBOUND) #& (delays[,'corr'] > 0.5)
mydelays <- delays[I,]

pdf('spline_delay_scatter.pdf', 87/25.4, 87/25.4)
par(ps=FONTSIZE, cex=1)
par(mar=c(1.5, 1.2, 0, 0)+0.4)
par(mgp=c(1, 0.4, 0))
par(mfrow=c(1, 1))
##colfunc <- colorRampPalette(c("red", "green"))
##bounds <- c(0.2, 0.18, 0.16, 0.14, 0.12, 0.1)
##colors <- colfunc(length(bounds))
foo <- mydelays[mydelays[,'ctd_resvarfrac']<RESVARBOUND,c('X50.', 'ctd_50.', 'ctd_resvarfrac')]
plot(foo[,'X50.'], foo[,'ctd_50.'], type='p', col='black', xlab='GP median delay (min)', ylab='Spline model median delay (min)')
lines(c(0, 120), c(0, 120))
cat(cor(foo[,'X50.'], foo[,'ctd_50.']), '\n')
dev.off()

I <- (delays[,'tmax'] < 160) & (delays[,'tmax'] > 1) & (delays[,'begdev10'] < DEVBOUND) & (delays[,'ctd_50.'] < 120) & (delays[,'ctd_resvarfrac'] < RESVARBOUND) & (delays[,'spline_var'] < POL2VARBOUND) #& (delays[,'corr'] > 0.5)
mydelays <- delays[I,]

