library(Biobase)
library(gplots)
library(RColorBrewer)
library(geneplotter)
library(HTqPCR)
library(openxlsx)
library(gdata)

files<-list.files(".../original_raw_data",full.names=T)

# rimuovo le linee inutili dai file di output dell'openarray
dir.create(raw_data)

for(i in 1:length(files)){
  prova<-readLines(files[i], n=3094)
  prova<-prova[22:length(prova)]
  write.table(prova,file=gsub("original_","",files[i]),sep="\t",row.names=F,quote=F,col.names=F)
}

# rimuovo le wells vuote e metto il file in ordine di campione e mirna
# prendo solo le colonne utili

colonne<-c("Well","Well.Position","Sample.Name","Target.Name","Crt","Crt.Mean","Crt.SD","Amp.Score","Cq.Conf","CRTAMPLITUDE","HIGHSD","CRTNOISE","ROX.Signal")

files<-list.files("C:/Users/milanimatteo/progetti/Chiodoni/Dati_OpenArray_2018/raw_data",full.names=T)
for(i in 1:length(files)){
  x<-read.table(files[i],header=T,sep="\t",as.is=T,quote="",comment.char="",fill=T)
  x<-x[,colnames(x)%in%colonne]
  x$Sample.Name<-gsub(" ","",x$Sample.Name)
  x<-x[x$Target.Name!="",]
  x<-x[order(x$Sample.Name,x$Target.Name),]
  write.table(x,file=files[i],sep="\t",row.names=F,quote=F)
}

#recupero i nomi dei campioni

samples<-NULL
filename<-NULL
for(i in 1:length(files)){
  x<-read.table(files[i],header=T,sep="\t",as.is=T,quote="",comment.char="",fill=T)
  samples<-c(samples,unique(substr(x$Sample.Name,1,nchar(x$Sample.Name)-2)))
  filename<-c(filename,rep(gsub(".*raw_data\\/","",files[i]),3))
}
samples<-data.frame(SampleID=samples,Sample.Name=samples,File.Name=filename,stringsAsFactors = F)

# controllo ordine delle colonne nei diversi files

colonne<-matrix(0,ncol=length(colonne))
for(i in 1:length(files)){
  x<-read.table(files[i],header=T,sep="\t",as.is=T,quote="",comment.char="",fill=T)
  colonne<-rbind(colonne,colnames(x))
}
colonne<-colonne[-1,]
sum(apply(colonne,2,function(x)isTRUE(all(x == x[1]))))==ncol(colonne)

# controllo ordine dei miRNA nei diversi files

annot<-matrix(rep(0,846),ncol=1)
for(i in 1:length(files)){
  x<-read.table(files[i],header=T,sep="\t",as.is=T,quote="",comment.char="",fill=T)
  annot<-cbind(annot,x[1:846,"Target.Name"],x[847:1692,"Target.Name"],x[1693:2538,"Target.Name"])
}
annot<-annot[,-1]
sum(apply(annot,1,function(x)isTRUE(all(x == x[1]))))==nrow(annot)

# recupero l'informazione della plate per i miRNA

x<-read.table(files[1],header=T,sep="\t",as.is=T,quote="",comment.char="",fill=T)
annot<-data.frame(Target.Name=x$Target.Name,Plate=substr(x$Sample.Name,nchar(x$Sample.Name),nchar(x$Sample.Name)))
annot<-annot[1:846,]

# creo i pdata
info<-read.table("sample_info.txt",sep="\t",header=T,as.is=T)
colnames(samples)[1]<-"GF_ID"
samples<-merge(samples,info,by="GF_ID",all.x=T)
samples<-samples[order(samples$File.Name),]
rownames(samples)<-samples$Sample.Name
identical(gsub(".*\\/","",files),unique(as.character(samples$File.Name)))

# in flag metto l'AmpScore

dataset<-readCtData(files,
                    header=T,
                    n.features = 846,
                    format="plain",
                    column.info=list(flag=8, feature=4, type=4, position=2, Ct=5),
                    n.data=3,
                    samples=as.vector(samples$Sample.Name),
                    na.value=40
)
pData(dataset)<-samples

# reimporto i dati mettendo in flag il CqConf e poi lo sostituisco in "featureCategory"
# del dataset precedente

x<-readCtData(files,
              header=T,
              n.features = 846,
              format="plain",
              column.info=list(flag=9, feature=4, type=4, position=2, Ct=5),
              n.data=3,
              samples=as.vector(samples$Sample.Name),
              na.value=40
)

featureCategory(dataset)<-flag(x)

# attacca l'annotazione della plate per i miRNA

identical(rownames(exprs(dataset)),as.vector(annot$Target.Name))
fData(dataset)$Plate<-annot$Plate

# converto a matrice i data frame con AmpScore e CqConf

for(i in 1:ncol(dataset))
  flag(dataset)[,i]<-as.numeric(flag(dataset)[,i])

for(i in 1:ncol(dataset))
  featureCategory(dataset)[,i]<-as.numeric(featureCategory(dataset)[,i])

# annotazione miRBASE 21
annot<-read.table("C:/Users/milanimatteo/annotazioni_geni/Mouse_OpenArray_miRBASE21/megaplex-pools-array-card-content.txt",
sep="\t",header=T,as.is=T,quote="",comment.char="")

annot<-annot[,c(3,5,6,7)]
single<-names(which(table(annot$Assay.ID)<2))
annot1<-annot[annot$Assay.ID%in%single|annot$Control.Assay!="",]
annot2<-annot[(!annot$Assay.ID%in%single)&annot$Control.Assay=="",]
annot2<-annot2[order(annot2$Assay.ID),]
annot2<-annot2[grep("mmu",annot2$miRBase_ID_21),]
annot2$miRBase_ID_21<-gsub(",rno.*","",annot2$miRBase_ID_21)
annot3<-annot[annot$miRBase_ID_21%in%""&annot$Control.Assay==""&(!annot$Assay.ID%in%single),]
annot.final<-rbind(annot1,annot2,annot3)
annot.final<-unique(annot.final)
fdata<-fData(dataset)
fdata$Assay.ID<-gsub(".*_","",fdata$featureNames)
fdata$rownames<-as.numeric(rownames(fdata))
fdata<-merge(fdata,annot.final,by="Assay.ID",all.x=T)
fdata<-fdata[order(fdata$rownames),]
fdata<-unique(fdata)
fdata<-fdata[!rownames(fdata)%in%c(388,493,511,531),]
rownames(fdata)<-fdata$rownames


identical(rownames(fdata),rownames(fData(dataset)))
fData(dataset)<-fdata


save(dataset,file="raw_data.RData")

fdata$Control.Assay[fdata$Control.Assay == ""] <- NA
fdata$miRBase_ID_21[fdata$miRBase_ID_21 == ""] <- NA

write.table(fdata,file="fdata.txt",sep="\t",col.names = T,row.names = F)


