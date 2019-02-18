rm(list=ls());gc()

load() #rawdata
load() #eset

pdf() #Boxplot
par(mfrow=c(2,1),mar=c(4,4.5,4,4),cex.axis=0.8,cex.lab=1.1,xpd=T,las=2) 
boxplot(log2(exprs(rawdata)),
        main="Non-normalized probe level data", 
        ylab="log2(intensity)", 
        cex=0.65,
        range=0,
        xaxt="n",
        col="light grey")
axis(1,at=1:ncol(rawdata),colnames(rawdata),las=2,cex.axis=0.45)

boxplot(exprs(eset),
        main="RMA data", 
        ylab="log2(intensity)", 
        cex=0.65,
        range=0,
        xaxt="n",
        col="grey")
axis(1,at=1:ncol(eset),eset$ID_G.,las=2,cex.axis=0.45)
dev.off()


#---------------------------------------------------------------------------------------------------------------
# HEATMAP
#---------------------------------------------------------------------------------------------------------------

source("C:/Users/milanimatteo/funzioni/heat_cluster.R")
source("C:/Users/milanimatteo/funzioni/heatmap_33.R")

correl1<-cor(log2(exprs(rawdata)))
correl2<-cor(exprs(eset))


pdf() #HEATMAP
par(xpd=TRUE, cex.axis=0.8)

#rawdata
LISTA <- heat.cluster(eset = rawdata,
                      col.cluster = c("Ceppo.topo","gruppi","RIN.class"), #Classi del pDATA
                      col.names.color = c("Ceppo","Gruppo","RIN")) #Nomi delle classi sull'HEATMAP

sample <- as.character(unlist(LISTA$sample.legend))
color <- as.character(unlist(LISTA$color.legend))
vett.color <- as.matrix(LISTA$color.named.vector)


#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="ward.D")}	

par(xpd=TRUE,cex.axis=0.7)
heatmap.3(correl1,
          key = TRUE,
          ColSideColors=vett.color,
          ColSideColorsSize = 2.6,
          main = "Non-normalized probe level data",
          col=redgreen(75),
          symkey=FALSE, 
          density.info="none", 
          trace="none",          
          distfun=mydist,
          hclustfun= myclust,
          cexCol = 0.6,
          Rowv = TRUE,
          labCol=colnames(correl1),
          labRow=FALSE)

legend("topleft",legend=sample, fill=color, border=FALSE,
       bty="n", y= 0.7,x = -0.1, cex=0.6)

#eset
LISTA <- heat.cluster(eset = eset,
                      col.cluster = c(), #Classi del pDATA
                      col.names.color = c()) #Nomi delle classi sull'HEATMAP

sample <- as.character(unlist(LISTA$sample.legend))
color <- as.character(unlist(LISTA$color.legend))
vett.color <- as.matrix(LISTA$color.named.vector)

heatmap.3(correl2,
          key = TRUE,
          ColSideColors=vett.color,
          ColSideColorsSize = 2.6,
          main = "RMA Data",
          col=redgreen(75),
          symkey=FALSE, 
          density.info="none", 
          trace="none",          
          distfun=mydist,
          hclustfun= myclust,
          cexCol = 0.6,
          Rowv = TRUE,
          labCol=colnames(correl2),
          labRow=FALSE)

legend("topleft",legend=sample, fill=color, border=FALSE,
       bty="n", y= 0.7,x = -0.1, cex=0.6)

dev.off()
