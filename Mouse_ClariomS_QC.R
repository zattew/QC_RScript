#------------------------------------------------------------------------------------
#pos_vs_neg AUC
#------------------------------------------------------------------------------------

#dataset normalizzato
load()

pos_neg<-exprs(eset)[grep("MG", rownames(exprs(eset))),]

neg<-pos_neg[grep("MG-1-neg",rownames(pos_neg)),]
pos<-pos_neg[grep("MG-1-pos",rownames(pos_neg)),]
neg<-gsub("*neg-","",gsub("HTA2-","",rownames(neg)))
pos<-gsub("*pos-","",gsub("HTA2-","",rownames(pos)))

category<-rep(NA,rep(dim(pos_neg)[1]))
category[grep("MG-1-neg",rownames(pos_neg))]<-0 
category[grep("MG-1-pos",rownames(pos_neg))]<-1

auc <- apply(pos_neg,2,function(x){ roc(response=category, predictor = x)$auc  })
ispass<-ifelse(auc>0.7,"","*")
auc.pass <- ispass

#Nome pdf curva AUC
pdf()

par(xpd=TRUE,las=2)
par(xpd=F,las=2)
bp<-barplot2(auc,
             xaxt="n",
             ylim=c(0.3,1),
             xpd=FALSE,
             ylab="AUC",
             main='Pos vs Neg AUC',
             col="light green"
)
axis(side=1,at=bp,labels=paste0(names(auc),ispass),las=2, cex.axis=0.6,tick=F)
abline(h=0.7,lty=2,lwd=2,col=2)
dev.off()

#------------------------------------------------------------------------------------
#Hybridization Control
#------------------------------------------------------------------------------------

grep("AFFX", rownames(exprs(eset)),value=T)
cre<-exprs(eset)[rownames(exprs(eset))=="AFFX-r2-P1-cre-5_at",]
biod<-exprs(eset)[rownames(exprs(eset))=="AFFX-r2-Ec-bioD-5_at",]
bioc<-exprs(eset)[rownames(exprs(eset))=="AFFX-r2-Ec-bioC-5_at",]
biob<-exprs(eset)[rownames(exprs(eset))=="AFFX-r2-Ec-bioB-5_at",]
subset<-t(exprs(eset)[rownames(exprs(eset))%in%c("AFFX-r2-P1-cre-5_at","AFFX-r2-Ec-bioD-5_at","AFFX-r2-Ec-bioC-5_at","AFFX-r2-Ec-bioB-5_at"),])

cols<-brewer.pal(n=4,name="Set1")

ispass<-subset[,4]>subset[,3]&subset[,3]>subset[,2]&subset[,2]>subset[,1]
ispass<-ifelse(ispass=="TRUE","","*")
hyb.pass <- ispass

#Controllo Ibridazione
pdf()

par(las = 2, cex.axis=1, xpd=TRUE, mar=c(5,5,5,5)) 
matplot(subset,
        pch=16,
        col=cols,
        type="o",
        ylab = "Log2(intensity)",
        main="Hybridization Control",
        mar=c(6,6),
        xaxt="n"
)
axis(side=1,at=1:nrow(subset),labels=paste0(rownames(subset),ispass),las=2, cex.axis=0.45)
legend("topright",
       legend=colnames(subset),
       pch=16,
       bty="n",
       cex=0.7,
       inset = c(-0.1,-0.16),
       col=cols
)
dev.off()

#------------------------------------------------------------------------------------
#Labeling Control
#------------------------------------------------------------------------------------

subset<-t(exprs(eset)[rownames(exprs(eset))%in%c("AFFX-r2-Bs-lys-5_st","AFFX-r2-Bs-phe-5_st","AFFX-r2-Bs-thr-5_s_st","AFFX-r2-Bs-dap-5_st"),])

subset<-subset[,c(2,3,4,1)]
cols<-brewer.pal(n=4,name="Set1")
ispass<-subset[,4]>subset[,3]&subset[,3]>subset[,2]&subset[,2]>subset[,1]
ispass<-ifelse(ispass=="TRUE","","*")
lab.pass <- ispass

#Controllo di Labeling
pdf()

par(las = 2, cex.axis=1, xpd=TRUE, mar=c(5,5,5,5)) 
matplot(subset,
        pch=16,
        col=cols,
        type="o",
        ylab = "Log2(intensity)",
        main="Labeling Control",
        mar=c(6,6),
        xaxt="n"
)
axis(side=1,at=1:nrow(subset),labels=paste0(rownames(subset),ispass),las=2, cex.axis=0.45)
legend("topright",
       legend=colnames(subset),
       pch=16,
       bty="n",
       cex=0.7,
       inset = c(-0.1,-0.16),
       col=cols
)
dev.off()
