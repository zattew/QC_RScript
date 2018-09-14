library(oligo)
library(limma)
library(factoextra)

load() #eset

df <- cor(exprs(eset))

#prova plot per vedere a che altezza tagliare l'abero per il batch
d <- dist(df, method = "euclidean")
hc1 <- hclust(d, method = "ward.D")
plot(hc1)


pdf() #Clustering in base al batch

for( i in ... ) #colonne in cui ci sono le caratteristiche di interesse (es 7:10)
{
  colnames(df) <- as.character(unlist(pData(eset)[i]))
  rownames(df) <- as.character(unlist(pData(eset)[i]))

  d <- dist(df, method = "euclidean")
  hc1 <- hclust(d, method = "ward.D" )
  
  par(cex=0.6)
  plot(hc1,main=names(pData(eset))[i]
  abline(h= ... ,col="red",lty=2,lwd=2)
  
  text(15, 3, labels = "Cluster1", pos=3,cex=2)
  text(62, 3,"Cluster2",pos=3,cex=2)
  
  eset$batch <- cutree(hc1,h= ...)

  boxplot( as.numeric(unlist(pData(eset)[i])) ~ eset$batch,
          main=names(pData(eset))[i], names=c("cluster1","cluster2"))
          
  legend("topright",legend=paste("p=",round(t.test(as.numeric(unlist(pData(eset)[i])) ~ eset$batch)$p.value,3)),bty="n")
}
dev.off()
