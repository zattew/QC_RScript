rm(list=ls());gc()
setwd('C:/Users/milanimatteo/progetti/....')

library(Biobase)
library(gplots)
library(geneplotter)
library(RColorBrewer)
library(limma)
library(gdata)

source('C:/Users/milanimatteo/funzioni/qc.r')

##----------------------------------------------------------------------------------------------------
#Hallmarks
##----------------------------------------------------------------------------------------------------

files<-list.files("directory","gsea_report_for_na",recursive=T,full.names=T)
files<-grep(".xls",files,value=T)

lista <- c("BALBc","NeuT")

conta <- 0
for(j in lista)
{
  fi <- grep(j,files,value=T)
  
  for(i in seq(1,(length(fi)-1),by=2)){
    print(i)
    x<-read.table(fi[i],sep="\t",header=T,as.is=T)
    y<-read.table(fi[i+1],sep="\t",header=T,as.is=T)
    temp<-rbind(x,y)
    temp<-temp[,colnames(temp)%in%c("NAME","NES","FDR.q.val")]
    temp<-temp[order(temp$NAME),]
    rownames(temp)<-temp$NAME
    
    NES <- as.matrix(temp$NES)
    rownames(NES) <- rownames(temp)
    FDR <- temp$FDR.q.val
    
    
    #NeuT
    if ( length(grep("NeuT_12w_Vs_NeuT_06w",fi[i])) != 0 )
    {
      colnames(NES)[1] <- "NeuT_12w_Vs_NeuT_06w"
    }
    
    if ( length(grep("NeuT_24w_Vs_NeuT_06w",fi[i])) != 0 )
    {
      colnames(NES)[1] <- "NeuT_24w_Vs_NeuT_06w"
    }
    
    if ( length(grep("NeuT_24w_Vs_NeuT_12w",fi[i])) != 0 )
    {
      colnames(NES)[1] <- "NeuT_24w_Vs_NeuT_12w"
    }
    
    
    #BALBc
    if ( length(grep("BALBc_12w_Vs_BALBc_06w",fi[i])) != 0 )
    {
      colnames(NES)[1] <- "BALBc_12w_Vs_BALBc_06w"
    }
    
    if ( length(grep("BALBc_24w_Vs_BALBc_06w",fi[i])) != 0 )
    {
      colnames(NES)[1] <- "BALBc_24w_Vs_BALBc_06w"
    }
    
    if ( length(grep("BALBc_24w_Vs_BALBc_12w",fi[i])) != 0 )
    {
      colnames(NES)[1] <- "BALBc_24w_Vs_BALBc_12w"
    }
    
    
    if ( conta == 0)
    {
      z<-NES
      f<-FDR
    }
    else
    {
      z<-cbind(z,NES)
      f<-cbind(f,FDR)
    }
    
    conta <- conta + 1
  }
  
  z <- z[rowSums(f < 0.05) > 1,]
  f <- f[rowSums(f < 0.05) > 1,]
  
  cellnotes <- ifelse(f < 0.05,"","x")
  
  z2<-z
    
    
  if( j == "BALBc")
  {
    z2 <- z2[,c("BALBc_12w_Vs_BALBc_06w",
                "BALBc_24w_Vs_BALBc_12w",
                "BALBc_24w_Vs_BALBc_06w")]
  }
  
  if( j == "NeuT")
  {
    z2 <- z2[,c("NeuT_12w_Vs_NeuT_06w",
                "NeuT_24w_Vs_NeuT_12w",
                "NeuT_24w_Vs_NeuT_06w")]
  }
  
  z2<-as.matrix(z2)
  mycol<- brewer.pal(11,"RdBu")
  mycol <- mycol[11:1]
  
  pdf(paste(j,"pdf",sep=".")) 
  par(cex.axis=0.6) #mar=c(4,3,3,5))
  heatmap.2(z2,
            cellnote = cellnotes,
            notecol = "black",
            density.info="none",
            trace="none",
            col=mycol,
            Colv = F,
            Rowv = F,
            dendrogram="none",
            labRow = gsub("HALLMARK_","",rownames(z2)),
            labCol =colnames(z2),
            cexRow=0.7,
            cexCol = 0.75,
            mar=c(10,18),
            key.xlab = "NES",
            key.title="NES",
            keysize = 1.1,
            lwid=c(0.5,2),
            main=j,
            colsep=0:ncol(z2),
            rowsep=0:nrow(z2),
            sepcolor="black",
            sepwidth=c(0.000001,0.000001)
  )
  
  legend("topright",legend = "X = FDR >= 0.05",
         bty="n",border=F,xpd=T)
  
  dev.off()
  
  conta <- 0
  
}
