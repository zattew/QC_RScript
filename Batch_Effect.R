library(oligo)
library(limma)
library(factoextra)

load() #eset

df <- cor(exprs(eset))

#prova plot per vedere a che altezza tagliare l'abero per il batch
d <- dist(df, method = "euclidean")
hc1 <- hclust(d, method = "ward.D")
plot(hc1)


pdf( ... ,height=5,width=8) #Clustering in base al batch
par(xpd=F)

for( i in ... ) #colonne in cui ci sono le caratteristiche di interesse (es 7:10)
{
  colnames(df) <- as.character(unlist(pData(eset)[i]))
  rownames(df) <- as.character(unlist(pData(eset)[i]))

  d <- dist(df, method = "euclidean")
  hc1 <- hclust(d, method = "ward.D" )
  
  plot(hc,
     cex=0.6,
     xlab="",
     ylab="Euclidean distance",
     hang=-1,
     main="Distance between samples")
  
  abline(h= ... ,col="red",lty=2,lwd=2)
  text(24,8,"Cluster 1",cex=1.1,font=2)
  text(53.5,8,"Cluster 2",cex=1.1,font=2)

  eset$batch <- cutree(hc1,h= ...)

  table(eset$clustertree,as.character(unlist(pData(eset)[i])))

  ##Test del X^2
  chisq.test(eset$batch,as.character(unlist(eset$description)), correct=FALSE) 
  
  boxplot( as.numeric(unlist(pData(eset)[i])) ~ eset$batch,
          main=names(pData(eset))[i], names=c("cluster1","cluster2"))
          
  legend("topright",legend=paste("p=",round(t.test(as.numeric(unlist(pData(eset)[i])) ~ eset$batch)$p.value,3)),bty="n")
}
dev.off()
