library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggfortify)
library(dplyr)
library(magrittr)
library(tidyr)

data("ciona16_sce")
zeroIndx=rowSums(counts(ciona16_sce))==0
ssData <- ciona16_sce[!isSpike(ciona16_sce) & !zeroIndx,]

phiValues <- assay(ssData,i="phicounts")
corData <- cor(phiValues)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

colSideCols <- brewer.pal(4, "Set1")[as.numeric(colData(ciona16_sce)$Embryo)]
rowSideCols <- brewer.pal(8, "Dark2")[as.numeric(colData(ciona16_sce)$CellType)]
labRow <- colData(ciona16_sce)$CellType

pdf("images/heatmap.pdf")
heatmap.2(corData,Colv = "Rowv",distfun=function(c) as.dist(1 - c),col=hmcol,symm=TRUE,revC=TRUE,RowSideColors=rowSideCols,ColSideColors = colSideCols,trace="none",
          labRow=labRow,labCol=colData(ciona16_sce)$Embryo,srtCol=0,density.info="histogram",key.xlab="Correlation",key.ylab="Count",key.title="Histogram key",xlab="Embryo",ylab="Cell Type")
dev.off()

data <- t(phiValues)
summary(prcomp(data))
gp=autoplot(prcomp(data),label=FALSE,shape = FALSE)
gp2=gp+geom_label(aes(label=colData(ciona16_sce)$CellType,colour=colData(ciona16_sce)$Embryo))+scale_colour_brewer(palette="Set1")
gp2
pdf("images/pca.pdf")
gp2+guides(colour=FALSE)
dev.off()

hiSeqGenesDetected <- data.frame(detected=colSums(counts(ciona16_sce)>0))
hiSeqGenesDetected$embryo <- colData(ciona16_sce)$Embryo
hiSeqGenesDetected$cell <- colData(ciona16_sce)$CellType
hiSeqGenesDetected$sequencer <- "HiSeq"

miSeqGenesDetected <- data.frame(detected=colSums(counts(ciona16MiSeq_sce)>0))
miSeqGenesDetected$embryo <- colData(ciona16MiSeq_sce)$Embryo
miSeqGenesDetected$cell <- colData(ciona16MiSeq_sce)$CellType
miSeqGenesDetected$sequencer <- "MiSeq"

genesDetected <- rbind(hiSeqGenesDetected,miSeqGenesDetected)

pdf("images/genesDetectedBoxPlot.pdf")
p <- ggplot(genesDetected, aes(x=embryo, y=detected,fill = sequencer)) + 
  geom_boxplot()
p
dev.off()

data <- genesDetected %>% group_by(sequencer) 
reorgData <- spread(data,sequencer,detected)
pdf("images/genesDetectedScatter.pdf")
p <- ggplot(reorgData, aes(x=MiSeq, y=HiSeq))
p + geom_point(aes(colour = embryo))
dev.off()
