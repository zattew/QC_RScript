heat.cluster <- function(eset,col.cluster,col.names.color)
{

  #Define custom dist and hclust functions for use with heatmaps
  mydist=function(c) {dist(c,method="euclidian")}
  myclust=function(c) {hclust(c,method="ward.D")}	

  #Matrice di colori per colonne 
  matrx.color <- rep(NA,ncol(eset))

  #brewer.pal.info
  a <- c(brewer.pal(9,"Set1"),brewer.pal(4,"BuPu"),brewer.pal(8,"Dark2"),brewer.pal(12,"Paired"),
         brewer.pal(9,"Purples"),brewer.pal(9,"BuGn"),brewer.pal(9,"Oranges"))

  prec <- 0
  conta <- 0
  conta1 <- 0
  
  for ( j in 1:length(col.cluster))
  {
  
    temp.color <- matrx.color
    type.temp <- unique(pData(eset)[,colnames(pData(eset))%in%col.cluster[j]])

    for ( i in 1:length(type.temp) )
    {
      pos <- pData(eset)[,colnames(pData(eset))%in%col.cluster[j]]%in%type.temp[i]
      temp.color[pos] <- a[i+prec]
      row <- data.frame(sample=type.temp[i],color=a[i+prec])
      
      if ( conta == 0)
      {
        legend.data <- row
      }
      
      else
      {
        legend.data <- rbind(legend.data,row)
      }
      
      conta <- conta + 1
    }
    
    if ( conta1 == 0)
    {
      col.color <- as.matrix(temp.color)
      
    }
    
    else
    {
      col.color <- cbind(col.color,as.matrix(temp.color))
    }
    
    conta1 <- conta1 + 1
    
    prec <- length(legend.data$color)
  }

  colnames(col.color) <- col.names.color
  oggetto <- list("sample.legend"=legend.data$sample,"color.legend"=legend.data$color,"color.named.vector"=col.color)
  
  return(oggetto)
}
  
