data(ciona16_sce)
phiValues <-assay(ciona16_sce,i="phicounts")

geneNames <- read.csv("data-raw/cionaGeneNames.csv",header=FALSE)
colnames(geneNames) <- c("GeneID","GeneName")
geneNames <- geneNames[!is.na(geneNames$GeneID),]
rownames(geneNames)=geneNames$GeneID

cellGroups=colData(ciona16_sce)$CellType
lq=function(group) quantile(group,0.25,type=5)
uq=function(group) quantile(group,0.75,type=5)
cellTypesCanonicalOrder=c("a5.3","a5.4","b5.3","b5.4","A5.1","A5.2","B5.1","B5.2")   
patternCode=function(cellGroups) {
 inPattern= cellTypesCanonicalOrder %in% cellGroups
 sum(2^((length(cellTypesCanonicalOrder)-1):0)[inPattern])
}

scorePattern=function(rowx) {
  #patternGroups=cutree(hclust(dist(tapply(rowx,cellGroups,mean)),method="single"),2)
  patternGroups=cutree(hclust(dist(do.call(rbind,tapply(rowx,cellGroups,identity))),method="single"),2)
  indx=patternGroups==1
  group1cells=names(patternGroups[indx])
  group1=rowx[cellGroups %in% group1cells]
  group2=rowx[!cellGroups %in% group1cells]
  if(mean(group1) > mean(group2)) {
    c(patternCode(group1cells),lq(group1)-uq(group2))
  } else {
    c(patternCode(names(patternGroups[!indx])),lq(group2)-uq(group1))
  }
}

result=apply(phiValues,MARGIN=1,scorePattern)
row.names(result) <- c("patternCode","score")
patternScores <- data.frame(t(result))

# patterns for top 40
top40=patternScores[order(patternScores$score,decreasing=TRUE),][1:40,]
newNames=geneNames[row.names(top40),"GeneName"]
indx=!is.na(newNames)
top40$GeneName <- row.names(top40)
top40$GeneName[indx] <- as.character(newNames[indx])
top40Patterns=split(top40,top40$patternCode)
for(patternCode in unique(top40$patternCode)) {
  print(rev(as.integer(intToBits(patternCode))[1:8]))
  print(top40[top40$patternCode==patternCode,])
  
}

