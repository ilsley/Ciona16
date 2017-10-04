library(SingleCellExperiment)
hiSeqCountData = read.csv("data-raw/hiSeqCountData.csv", row.names = 1, header = TRUE)
miSeqCountData = read.csv("data-raw/miSeqCountData.csv", row.names = 1, header = TRUE)

cellInfo = read.csv("data-raw/phenoData.csv", row.names = 1, header= TRUE,colClasses = "factor")
cellInfo = cellInfo[order(cellInfo$Embryo,cellInfo$CellType),]
sampleNames <- rownames(cellInfo)

geneNames <- read.csv("data-raw/cionaGeneNames.csv",header=FALSE)
colnames(geneNames) <- c("GeneID","GeneName")
geneNames <- geneNames[!is.na(geneNames$GeneID),]
rownames(geneNames)=geneNames$GeneID

ciona16_sce <- SingleCellExperiment(assays=SimpleList(counts=as.matrix(hiSeqCountData[,sampleNames])),colData=cellInfo)
ciona16MiSeq_sce <- SingleCellExperiment(assays=SimpleList(counts=as.matrix(miSeqCountData[,sampleNames])),colData=cellInfo)

addDetails <- function(sce) {
  isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))
  sizeFactors(sce) <- colSums(assay(sce))
  normcounts(sce) <- sweep(counts(sce),MARGIN=2,STATS = sizeFactors(sce),FUN = "/") 
  assay(sce, i="phicounts") <- 2*asin(sqrt(normcounts(sce)))
  rowData=merge(data.frame("GeneID"=rownames(sce)),geneNames,by="GeneID",all.x=TRUE)
  rowData(sce)=rowData
  sce
}

ciona16_sce <- addDetails(ciona16_sce)
ciona16MiSeq_sce <- addDetails(ciona16MiSeq_sce)

save(ciona16_sce,file="data/ciona16_sce.rda")
save(ciona16MiSeq_sce,file="data/ciona16MiSeq_sce.rda")

