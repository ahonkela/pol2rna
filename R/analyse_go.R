## source('analyse_intron_lengths_meddelay.R')

library(topGO)

ONTOLOGY <- "BP"
ann <- annFUN.org(ONTOLOGY, mapping="org.Hs.eg.db", ID="ensembl")

mygenes <- mydelays[,'meddelay']
names(mygenes) <- row.names(mydelays)

mgenes <- mygenes
##mgenes[mgenes==3] = 1
##mgenes[mgenes==4] = 2

longdelay <- function (x) { return (x > 15); }
shortdelay <- function (x) { return (x < 3); }

GOdata <- new("topGOdata", ontology = ONTOLOGY, allGenes = mgenes, geneSel = shortdelay, annot = annFUN.org, mapping="org.Hs.eg.db", ID="ensembl")
resultFis <- runTest(GOdata, algorithm='classic', statistic='fisher')
tableFis <- GenTable(GOdata, classic=resultFis, topNodes = length(resultFis@score))

resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
##resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
tableKS <- GenTable(GOdata, classic=resultKS, topNodes = length(resultKS@score))
##tableKS.elim <- GenTable(GOdata, elim=resultKS.elim, topNodes = 20)

tableFis['bh'] <- p.adjust(tableFis[,"classic"],method="BH")
tableKS['bh'] <- p.adjust(tableKS[,"classic"],method="BH")

head(tableFis)
head(tableKS)
##show(tableKS.elim)

##allRes <- GenTable(GOdata, classicFisher = resultFis,
##                   classicKS = resultKS, elimKS = resultKS.elim,
##                   orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 20)
