data("ciona16_sce")

addPatternContrasts <- function(phenoData,patternFilename) {
  patternText=scan(patternFilename,character())
  patterns=sapply(patternText,function(x) strsplit(x,",",fixed=TRUE))
  patternSets=lapply(patterns,function(patternSet) apply(phenoData, MARGIN=1,function(rowx) ifelse(rowx["CellType"] %in% patternSet,"Set1","Set2")))
  patternLabels=sapply(1:length(patternSets),function(index) paste("pattern",index,sep=""))
  names(patternSets)=patternLabels
  list(patternLabels=patternLabels,phenoData=cbind(phenoData,patternSets))
}

runDESeq <- function(patternFilename,exportFilename) {
  require(DESeq2)
  scedata <- ciona16_sce
  result <-addPatternContrasts(data.frame(lapply(colData(scedata),as.character)),patternFilename)
  patternLabels <- result$patternLabels
  phenoData <- result$phenoData

  results <- lapply(patternLabels,function(patternStr) {
    #dds <- DESeqDataSetFromMatrix(countData=counts(scedata),colData=phenoData,design = as.formula(paste("~ ",patternStr,sep="")))
    dds <- DESeqDataSetFromMatrix(countData=counts(scedata),colData=phenoData,design = as.formula(paste("~ Embryo + ",patternStr,sep="")))
    dds <- DESeq(dds)
    data=results(dds)
    data$pattern=patternStr
    data$ShortKHID=rownames(data)
    data
  })
  write.table(results[[1]],file=exportFilename,append=FALSE,col.names=TRUE,row.names=FALSE,sep=",")
  for(i in 2:length(results)) {write.table(results[[i]],file=exportFilename,append=TRUE,col.names=FALSE,row.names=FALSE,sep=",")}  
}

runDESeq("data-raw/allPatterns.csv","exports/allPatternsDESeqEmbryoPattern.csv")
runDESeq("data-raw/top40Patterns.csv","exports/top40DESeq.csv")
