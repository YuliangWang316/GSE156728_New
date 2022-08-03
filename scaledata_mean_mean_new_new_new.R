library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
library(tidyverse)
a<-list.files("c:/Users/xjmik/Downloads/data.expression/data/expression/CD4/byDataset/")
y<-list.files("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8/byDataset/")
z<-as.data.frame(t(rep(0,38)))
colnames(z)<-c("a.i.",                   "length.patient." ,      
               "mean.d.JMJD1C."  ,       "median.d.JMJD1C." ,     
               "sum.d.JMJD1C."  ,        "mean.e.IFNG." ,         
               "median.e.IFNG."   ,      "sum.e.IFNG."  ,         
               "mean.f.STAT3_mean."  ,   "median.f.STAT3_mean.",  
               "sum.f.STAT3_mean."  ,    "mean.f.STAT3_median.",  
               "median.f.STAT3_median." ,"sum.f.STAT3_median." ,  
               "mean.f.STAT3_sum.",      "median.f.STAT3_sum.",   
               "sum.f.STAT3_sum.",       "mean.f.PI3K_mean." ,    
               "median.f.PI3K_mean."  ,  "sum.f.PI3K_mean."   ,   
               "mean.f.PI3K_median."  ,  "median.f.PI3K_median." ,
               "sum.f.PI3K_median."   ,  "mean.f.PI3K_sum."  ,    
               "median.f.PI3K_sum."   ,  "sum.f.PI3K_sum."   ,    
               "mean.o.GZMB."         ,  "median.o.GZMB."     ,   
               "sum.o.GZMB."          ,  "mean.p.TNF."         ,  
               "median.p.TNF."        ,  "sum.p.TNF."           , 
               "mean.q.IFNG."         ,  "median.q.IFNG."       , 
               "sum.q.IFNG."          ,  "Treg_frequence"       , 
               "CD8_frequence"        ,  "CD8_ratio" )
for (i in c(1:23,25:length(a))) {
  setwd("c:/Users/xjmik/Downloads/data.expression/data/expression/CD4/byDataset/")
  b<-readRDS(a[i])
  setwd("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8/byDataset/")
  x<-readRDS(y[i])
  count_B<-b@assays$data@listData$norm_exprs
  count_B<-as.matrix(count_B)
  count_B<-as.data.frame(count_B)
  count_B_metadata<-as.data.frame(b@colData@listData)
  count_B<-count_B[,count_B_metadata$cellID]
  colnames(count_B)<-count_B_metadata$cellID.uniq
  genename<-rownames(count_B)
  d<-str_sub(genename[1],1,4)
  if(d == "ENSG"){
    e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
    e<-e[!duplicated(e$SYMBOL),]
    count_B<-count_B[e$ENSEMBL,]
    rownames(count_B)<-e$SYMBOL
  }
  f<-sort(genename)
  if(f[1] == "1"){
    g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
    g<-g[!duplicated(g$SYMBOL),]
    count_B<-count_B[g$ENTREZID,]
    rownames(count_B)<-g$SYMBOL
  }
  remove(b,e,g,d,f,genename)
  e<-readRDS("C:/Users/xjmik/Downloads/data.expression/data/expression/CD4/integration/int.CD4.S35.meta.tb.rds")
  rownames(count_B_metadata)<-count_B_metadata$cellID.uniq
  h<-intersect(e$cellID.uniq,count_B_metadata$cellID.uniq)
  e<-as.data.frame(e)
  rownames(e)<-e$cellID.uniq
  e_new<-e[h,]
  count_B_metadata_new<-count_B_metadata[h,]
  metadata<-cbind(e_new,count_B_metadata_new)
  remove(e,e_new,count_B_metadata,count_B_metadata_new,h)        #
  rownames(metadata)<-colnames(count_B)
  count_C<-x@assays$data@listData$norm_exprs
  count_C<-as.matrix(count_C)
  count_C<-as.data.frame(count_C)
  count_C_metadata<-as.data.frame(x@colData@listData)
  count_C<-count_C[,count_C_metadata$cellID]
  colnames(count_C)<-count_C_metadata$cellID.uniq
  genename<-rownames(count_C)
  d<-str_sub(genename[1],1,4)
  if(d == "ENSG"){
    e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
    e<-e[!duplicated(e$SYMBOL),]
    count_C<-count_C[e$ENSEMBL,]
    rownames(count_C)<-e$SYMBOL
  }
  f<-sort(genename)
  if(f[1] == "1"){
    g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
    g<-g[!duplicated(g$SYMBOL),]
    count_C<-count_C[g$ENTREZID,]
    rownames(count_C)<-g$SYMBOL
  }
  remove(x,e,g,d,f,genename)
  e<-readRDS("C:/Users/xjmik/Downloads/data.expression/data/expression/CD8/integration/int.CD8.S35.meta.tb.rds")
  rownames(count_C_metadata)<-count_C_metadata$cellID.uniq
  h<-intersect(e$cellID.uniq,count_C_metadata$cellID.uniq)
  e<-as.data.frame(e)
  rownames(e)<-e$cellID.uniq
  e_new<-e[h,]
  count_C_metadata_new<-count_C_metadata[h,]
  metadata_C<-cbind(e_new,count_C_metadata_new)
  remove(e,e_new,count_C_metadata,count_C_metadata_new,h)        #
  rownames(metadata_C)<-colnames(count_C)
  b<-intersect(rownames(count_B),rownames(count_C))
  count_B_new<-count_B[b,]
  count_C_new<-count_C[b,]
  count<-cbind(count_B_new,count_C_new)
  remove(count_B,count_B_new,count_C,count_C_new,b)
  b<-intersect(colnames(metadata),colnames(metadata_C))
  metadata_new<-metadata[,b]
  metadata_C_new<-metadata_C[,b]
  metadata<-rbind(metadata_new,metadata_C_new)
  remove(metadata_C,metadata_C_new,metadata_new,b)
  pbmc<-CreateSeuratObject(counts = count,meta.data = metadata)
  Idents(pbmc)<-pbmc@meta.data$loc
  remove(count,metadata)
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  Idents(pbmc)<-pbmc@meta.data$meta.cluster
  remove(all.genes)
  b<-read.table("C:/Users/xjmik/Downloads/STAT3.txt",sep = "\t",header = TRUE)
  c<-read.table("C:/Users/xjmik/Downloads/PI3K.txt",sep = "\t",header = TRUE)
  b_new<-b[2:length(rownames(b)),]
  c_new<-c[2:length(rownames(c)),]
  b_new<-intersect(b_new,rownames(pbmc))
  c_new<-intersect(c_new,rownames(pbmc))
  pbmc<-AddModuleScore(pbmc,features = b_new,name = "STAT3")
  pbmc<-AddModuleScore(pbmc,features = c_new,name = "PI3K")
  remove(b,c,b_new,c_new)
  b<-levels(Idents(pbmc))
  c<-b[which(b == "CD4.c18.Treg.RTKN2" | b == "CD4.c19.Treg.S1PR1" | b == "CD4.c20.Treg.TNFRSF9" | b == "CD4.c21.Treg.OAS1" | b == "CD4.c23.Mix.NME1" | b == "CD4.c24.Mix.NME2")]
  Treg<-subset(pbmc,ident=c)
  d<-b[grep("CD8",b)]
  CD8<-subset(pbmc,ident=d)
  e<-b[grep("CD4",b)]
  CD4<-subset(pbmc,ident=e)
  remove(b,c,d,e)
  Idents(Treg)<-Treg@meta.data$loc
  Treg_T<-subset(Treg,ident="T")
  Idents(Treg_T)<-Treg_T@meta.data$patient.uid
  patient_Treg<-unique(Treg_T@meta.data$patient.uid)
  remove(Treg)
  Idents(CD8)<-CD8@meta.data$loc
  CD8_T<-subset(CD8,ident="T")
  Idents(CD8_T)<-CD8_T@meta.data$patient.uid
  patient_CD8_T<-unique(CD8_T@meta.data$patient.uid)
  remove(CD8)
  Idents(CD4)<-CD4@meta.data$loc
  CD4_T<-subset(CD4,ident="T")
  Idents(CD4_T)<-CD4_T@meta.data$patient.uid
  patient_CD4_T<-unique(CD4_T@meta.data$patient.uid)
  remove(CD4)
  Idents(pbmc)<-pbmc@meta.data$loc
  pbmc_T<-subset(pbmc,ident="T")
  Idents(pbmc_T)<-pbmc_T@meta.data$patient.uid
  patient_pbmc_T<-unique(pbmc_T@meta.data$patient.uid)
  remove(pbmc)
  patient_1<-intersect(patient_CD8_T,patient_Treg)
  patient_2<-intersect(patient_CD4_T,patient_pbmc_T)
  patient<-intersect(patient_1,patient_2)
  remove(patient_CD8_T,patient_Treg,patient_1,patient_2,patient_CD4_T,patient_pbmc_T)
  for (j in 1:length(patient)) {
    assign(paste("Treg_",j,"_TIL",sep = ""),subset(Treg_T,idents = patient[j]))
    assign("b",slot(get(paste("Treg_",j,"_TIL",sep = "")),"assays"))
    assign("c",as.data.frame(t(as.data.frame(b$RNA@scale.data))))
    for (r in 1:length(colnames(c))) {
      noise<-rnorm(n=length(x=c[,r]))/1e+05
      c[,r]<-c[,r]+noise
    }
    remove(noise,r)
    assign("d",dplyr::select(c,one_of("JMJD1C")))
    assign("e",dplyr::select(c,one_of("IFNG")))
    assign("f",slot(get(paste("Treg_",j,"_TIL",sep = "")),"meta.data"))
    STAT3<-dplyr::select(f,starts_with("STAT3"))
    for (r in 1:length(colnames(STAT3))) {
      noise<-rnorm(n=length(x=STAT3[,r]))/1e+05
      STAT3[,r]<-STAT3[,r]+noise
    }
    remove(noise,r)
    f$STAT3_mean<-apply(STAT3,MARGIN = 1,mean)
    f$STAT3_median<-apply(STAT3,MARGIN = 1,median)
    f$STAT3_sum<-apply(STAT3,MARGIN = 1,sum)
    PI3K<-dplyr::select(f,starts_with("PI3K"))
    for (r in 1:length(colnames(PI3K))) {
      noise<-rnorm(n=length(x=PI3K[,r]))/1e+05
      PI3K[,r]<-PI3K[,r]+noise
    }
    remove(noise,r)
    f$PI3K_mean<-apply(PI3K,MARGIN = 1,mean)
    f$PI3K_median<-apply(PI3K,MARGIN = 1,median)
    f$PI3K_sum<-apply(PI3K,MARGIN = 1,sum)
    remove(b,c,PI3K,STAT3)
    assign("g",length(colnames(get(paste("Treg_",j,"_TIL",sep = "")))))
    assign(paste("CD4_",j,"_TIL",sep = ""),subset(CD4_T,idents = patient[j]))
    assign("h",length(colnames(get(paste("CD4_",j,"_TIL",sep = "")))))
    assign(paste("pbmc_",j,"_TIL",sep = ""),subset(pbmc_T,idents = patient[j]))
    assign("k",length(colnames(get(paste("pbmc_",j,"_TIL",sep = "")))))
    assign(paste("CD8_",j,"_TIL",sep = ""),subset(CD8_T,idents = patient[j]))
    assign("l",length(colnames(get(paste("CD8_",j,"_TIL",sep = "")))))
    assign("m",slot(get(paste("CD8_",j,"_TIL",sep = "")),"assays"))
    assign("n",as.data.frame(t(as.data.frame(m$RNA@scale.data))))
    for (r in 1:length(colnames(n))) {
      noise<-rnorm(n=length(x=n[,r]))/1e+05
      n[,r]<-n[,r]+noise
    }
    remove(noise,r)
    assign("o",dplyr::select(n,one_of("GZMB")))
    assign("p",dplyr::select(n,one_of("TNF")))
    assign("q",dplyr::select(n,one_of("IFNG")))
    remove(m,n)
    Treg_frequence<-g/h
    CD8_frequence<-l/k
    CD8_ratio<-l/g
    remove(g,h,k,l)
    assign(paste("scaldata_",j,"_TIL",sep = ""),data.frame(a[i],length(patient),
                                                           mean(d$JMJD1C),median(d$JMJD1C),sum(d$JMJD1C),
                                                           mean(e$IFNG),median(e$IFNG),sum(e$IFNG),
                                                           mean(f$STAT3_mean),median(f$STAT3_mean),sum(f$STAT3_mean),
                                                           mean(f$STAT3_median),median(f$STAT3_median),sum(f$STAT3_median),
                                                           mean(f$STAT3_sum),median(f$STAT3_sum),sum(f$STAT3_sum),
                                                           mean(f$PI3K_mean),median(f$PI3K_mean),sum(f$PI3K_mean),
                                                           mean(f$PI3K_median),median(f$PI3K_median),sum(f$PI3K_median),
                                                           mean(f$PI3K_sum),median(f$PI3K_sum),sum(f$PI3K_sum),
                                                           mean(o$GZMB),median(o$GZMB),sum(o$GZMB),
                                                           mean(p$TNF),median(p$TNF),sum(p$TNF),
                                                           mean(q$IFNG),median(q$IFNG),sum(q$IFNG),
                                                           Treg_frequence,CD8_frequence,CD8_ratio))
    remove(list =paste("Treg_",j,"_TIL",sep = ""))
    remove(list =paste("CD4_",j,"_TIL",sep = ""))
    remove(list =paste("CD8_",j,"_TIL",sep = ""))
    remove(list =paste("pbmc_",j,"_TIL",sep = ""))
    remove(d,e,f,o,p,q,CD8_frequence,CD8_ratio,Treg_frequence)
  }
  scaldata<-scaldata_1_TIL
  if (length(patient)>=2){
    for (j in 2:length(patient)) {
      scaldata<-rbind(scaldata,get(paste("scaldata_",j,"_TIL",sep = "")))
      remove(list = paste("scaldata_",j,"_TIL",sep = ""))
    }
  }
  remove(scaldata_1_TIL,j,CD4_T,CD8_T,pbmc_T,Treg_T)
  colnames(scaldata)<-c("a.i.",                   "length.patient." ,      
                        "mean.d.JMJD1C."  ,       "median.d.JMJD1C." ,     
                        "sum.d.JMJD1C."  ,        "mean.e.IFNG." ,         
                        "median.e.IFNG."   ,      "sum.e.IFNG."  ,         
                        "mean.f.STAT3_mean."  ,   "median.f.STAT3_mean.",  
                        "sum.f.STAT3_mean."  ,    "mean.f.STAT3_median.",  
                        "median.f.STAT3_median." ,"sum.f.STAT3_median." ,  
                        "mean.f.STAT3_sum.",      "median.f.STAT3_sum.",   
                        "sum.f.STAT3_sum.",       "mean.f.PI3K_mean." ,    
                        "median.f.PI3K_mean."  ,  "sum.f.PI3K_mean."   ,   
                        "mean.f.PI3K_median."  ,  "median.f.PI3K_median." ,
                        "sum.f.PI3K_median."   ,  "mean.f.PI3K_sum."  ,    
                        "median.f.PI3K_sum."   ,  "sum.f.PI3K_sum."   ,    
                        "mean.o.GZMB."         ,  "median.o.GZMB."     ,   
                        "sum.o.GZMB."          ,  "mean.p.TNF."         ,  
                        "median.p.TNF."        ,  "sum.p.TNF."           , 
                        "mean.q.IFNG."         ,  "median.q.IFNG."       , 
                        "sum.q.IFNG."          ,  "Treg_frequence"       , 
                        "CD8_frequence"        ,  "CD8_ratio" )
  
  z<-rbind(z,scaldata)
  
}
write.table(z,"c:/Users/xjmik/Desktop/scaledata_patientnoise.txt",sep = "\t")
#b<-cor(scaldata$JMJD1C,scaldata$IFNG,method = "spearman")
#c<-cor.test(scaldata$JMJD1C,scaldata$IFNG,method = "spearman",exact = TRUE)
#d<-data.frame(b,c$p.value,a[i],i,length(rownames(scaldata)))
#colnames(d)<-c("R","P","Name","i","Sample_Size")
#z<-rbind(z,d)
#remove(b,c,d,scaldata,Treg_T)
setwd("c:/Users/xjmik/Downloads/data.expression/data/expression/CD4/byDataset/")
b<-readRDS(a[24])
setwd("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8/byDataset/")
x<-readRDS(y[24])
count_B<-b@assays@data@listData$norm_exprs
count_B<-as.matrix(count_B)
count_B<-as.data.frame(count_B)
count_B_metadata<-as.data.frame(b@colData@listData)
count_B<-count_B[,count_B_metadata$cellID]
colnames(count_B)<-count_B_metadata$cellID.uniq
genename<-rownames(count_B)
d<-str_sub(genename[1],1,4)
if(d == "ENSG"){
  e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
  e<-e[!duplicated(e$SYMBOL),]
  count_B<-count_B[e$ENSEMBL,]
  rownames(count_B)<-e$SYMBOL
}
f<-sort(genename)
if(f[1] == "1"){
  g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
  g<-g[!duplicated(g$SYMBOL),]
  count_B<-count_B[g$ENTREZID,]
  rownames(count_B)<-g$SYMBOL
}
remove(b,e,g,d,f,genename)
e<-readRDS("C:/Users/xjmik/Downloads/data.expression/data/expression/CD4/integration/int.CD4.S35.meta.tb.rds")
rownames(count_B_metadata)<-count_B_metadata$cellID.uniq
h<-intersect(e$cellID.uniq,count_B_metadata$cellID.uniq)
e<-as.data.frame(e)
rownames(e)<-e$cellID.uniq
e_new<-e[h,]
count_B_metadata_new<-count_B_metadata[h,]
metadata<-cbind(e_new,count_B_metadata_new)
remove(e,e_new,count_B_metadata,count_B_metadata_new,h)        #
rownames(metadata)<-colnames(count_B)
count_C<-x@assays@data@listData$norm_exprs
count_C<-as.matrix(count_C)
count_C<-as.data.frame(count_C)
count_C_metadata<-as.data.frame(x@colData@listData)
count_C<-count_C[,count_C_metadata$cellID]
colnames(count_C)<-count_C_metadata$cellID.uniq
genename<-rownames(count_C)
d<-str_sub(genename[1],1,4)
if(d == "ENSG"){
  e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
  e<-e[!duplicated(e$SYMBOL),]
  count_C<-count_C[e$ENSEMBL,]
  rownames(count_C)<-e$SYMBOL
}
f<-sort(genename)
if(f[1] == "1"){
  g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
  g<-g[!duplicated(g$SYMBOL),]
  count_C<-count_C[g$ENTREZID,]
  rownames(count_C)<-g$SYMBOL
}
remove(x,e,g,d,f,genename)
e<-readRDS("C:/Users/xjmik/Downloads/data.expression/data/expression/CD8/integration/int.CD8.S35.meta.tb.rds")
rownames(count_C_metadata)<-count_C_metadata$cellID.uniq
h<-intersect(e$cellID.uniq,count_C_metadata$cellID.uniq)
e<-as.data.frame(e)
rownames(e)<-e$cellID.uniq
e_new<-e[h,]
count_C_metadata_new<-count_C_metadata[h,]
metadata_C<-cbind(e_new,count_C_metadata_new)
remove(e,e_new,count_C_metadata,count_C_metadata_new,h)        #
rownames(metadata_C)<-colnames(count_C)
b<-intersect(rownames(count_B),rownames(count_C))
count_B_new<-count_B[b,]
count_C_new<-count_C[b,]
count<-cbind(count_B_new,count_C_new)
remove(count_B,count_B_new,count_C,count_C_new,b)
b<-intersect(colnames(metadata),colnames(metadata_C))
metadata_new<-metadata[,b]
metadata_C_new<-metadata_C[,b]
metadata<-rbind(metadata_new,metadata_C_new)
remove(metadata_C,metadata_C_new,metadata_new,b)
pbmc<-CreateSeuratObject(counts = count,meta.data = metadata)
Idents(pbmc)<-pbmc@meta.data$loc
remove(count,metadata)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
Idents(pbmc)<-pbmc@meta.data$meta.cluster
remove(all.genes)
b<-read.table("C:/Users/xjmik/Downloads/STAT3.txt",sep = "\t",header = TRUE)
c<-read.table("C:/Users/xjmik/Downloads/PI3K.txt",sep = "\t",header = TRUE)
b_new<-b[2:length(rownames(b)),]
c_new<-c[2:length(rownames(c)),]
b_new<-intersect(b_new,rownames(pbmc))
c_new<-intersect(c_new,rownames(pbmc))
pbmc<-AddModuleScore(pbmc,features = b_new,name = "STAT3")
pbmc<-AddModuleScore(pbmc,features = c_new,name = "PI3K")
remove(b,c,b_new,c_new)
b<-levels(Idents(pbmc))
c<-b[which(b == "CD4.c18.Treg.RTKN2" | b == "CD4.c19.Treg.S1PR1" | b == "CD4.c20.Treg.TNFRSF9" | b == "CD4.c21.Treg.OAS1" | b == "CD4.c23.Mix.NME1" | b == "CD4.c24.Mix.NME2")]
Treg<-subset(pbmc,ident=c)
d<-b[grep("CD8",b)]
CD8<-subset(pbmc,ident=d)
e<-b[grep("CD4",b)]
CD4<-subset(pbmc,ident=e)
remove(b,c,d,e)
Idents(Treg)<-Treg@meta.data$loc
Treg_T<-subset(Treg,ident="T")
Idents(Treg_T)<-Treg_T@meta.data$patient.uid
patient_Treg<-unique(Treg_T@meta.data$patient.uid)
remove(Treg)
Idents(CD8)<-CD8@meta.data$loc
CD8_T<-subset(CD8,ident="T")
Idents(CD8_T)<-CD8_T@meta.data$patient.uid
patient_CD8_T<-unique(CD8_T@meta.data$patient.uid)
remove(CD8)
Idents(CD4)<-CD4@meta.data$loc
CD4_T<-subset(CD4,ident="T")
Idents(CD4_T)<-CD4_T@meta.data$patient.uid
patient_CD4_T<-unique(CD4_T@meta.data$patient.uid)
remove(CD4)
Idents(pbmc)<-pbmc@meta.data$loc
pbmc_T<-subset(pbmc,ident="T")
Idents(pbmc_T)<-pbmc_T@meta.data$patient.uid
patient_pbmc_T<-unique(pbmc_T@meta.data$patient.uid)
remove(pbmc)
patient_1<-intersect(patient_CD8_T,patient_Treg)
patient_2<-intersect(patient_CD4_T,patient_pbmc_T)
patient<-intersect(patient_1,patient_2)
remove(patient_CD8_T,patient_Treg,patient_1,patient_2,patient_CD4_T,patient_pbmc_T)
for (j in 1:length(patient)) {
  assign(paste("Treg_",j,"_TIL",sep = ""),subset(Treg_T,idents = patient[j]))
  assign("b",slot(get(paste("Treg_",j,"_TIL",sep = "")),"assays"))
  assign("c",as.data.frame(t(as.data.frame(b$RNA@scale.data))))
  for (r in 1:length(colnames(c))) {
    noise<-rnorm(n=length(x=c[,r]))/1e+05
    c[,r]<-c[,r]+noise
  }
  remove(noise,r)
  assign("d",dplyr::select(c,one_of("JMJD1C")))
  assign("e",dplyr::select(c,one_of("IFNG")))
  assign("f",slot(get(paste("Treg_",j,"_TIL",sep = "")),"meta.data"))
  STAT3<-dplyr::select(f,starts_with("STAT3"))
  for (r in 1:length(colnames(STAT3))) {
    noise<-rnorm(n=length(x=STAT3[,r]))/1e+05
    STAT3[,r]<-STAT3[,r]+noise
  }
  remove(noise,r)
  f$STAT3_mean<-apply(STAT3,MARGIN = 1,mean)
  f$STAT3_median<-apply(STAT3,MARGIN = 1,median)
  f$STAT3_sum<-apply(STAT3,MARGIN = 1,sum)
  PI3K<-dplyr::select(f,starts_with("PI3K"))
  for (r in 1:length(colnames(PI3K))) {
    noise<-rnorm(n=length(x=PI3K[,r]))/1e+05
    PI3K[,r]<-PI3K[,r]+noise
  }
  remove(noise,r)
  f$PI3K_mean<-apply(PI3K,MARGIN = 1,mean)
  f$PI3K_median<-apply(PI3K,MARGIN = 1,median)
  f$PI3K_sum<-apply(PI3K,MARGIN = 1,sum)
  remove(b,c,PI3K,STAT3)
  assign("g",length(colnames(get(paste("Treg_",j,"_TIL",sep = "")))))
  assign(paste("CD4_",j,"_TIL",sep = ""),subset(CD4_T,idents = patient[j]))
  assign("h",length(colnames(get(paste("CD4_",j,"_TIL",sep = "")))))
  assign(paste("pbmc_",j,"_TIL",sep = ""),subset(pbmc_T,idents = patient[j]))
  assign("k",length(colnames(get(paste("pbmc_",j,"_TIL",sep = "")))))
  assign(paste("CD8_",j,"_TIL",sep = ""),subset(CD8_T,idents = patient[j]))
  assign("l",length(colnames(get(paste("CD8_",j,"_TIL",sep = "")))))
  assign("m",slot(get(paste("CD8_",j,"_TIL",sep = "")),"assays"))
  assign("n",as.data.frame(t(as.data.frame(m$RNA@scale.data))))
  for (r in 1:length(colnames(n))) {
    noise<-rnorm(n=length(x=n[,r]))/1e+05
    n[,r]<-n[,r]+noise
  }
  remove(noise,r)
  assign("o",dplyr::select(n,one_of("GZMB")))
  assign("p",dplyr::select(n,one_of("TNF")))
  assign("q",dplyr::select(n,one_of("IFNG")))
  remove(m,n)
  Treg_frequence<-g/h
  CD8_frequence<-l/k
  CD8_ratio<-l/g
  remove(g,h,k,l)
  assign(paste("scaldata_",j,"_TIL",sep = ""),data.frame(a[24],length(patient),
                                                         mean(d$JMJD1C),median(d$JMJD1C),sum(d$JMJD1C),
                                                         mean(e$IFNG),median(e$IFNG),sum(e$IFNG),
                                                         mean(f$STAT3_mean),median(f$STAT3_mean),sum(f$STAT3_mean),
                                                         mean(f$STAT3_median),median(f$STAT3_median),sum(f$STAT3_median),
                                                         mean(f$STAT3_sum),median(f$STAT3_sum),sum(f$STAT3_sum),
                                                         mean(f$PI3K_mean),median(f$PI3K_mean),sum(f$PI3K_mean),
                                                         mean(f$PI3K_median),median(f$PI3K_median),sum(f$PI3K_median),
                                                         mean(f$PI3K_sum),median(f$PI3K_sum),sum(f$PI3K_sum),
                                                         mean(o$GZMB),median(o$GZMB),sum(o$GZMB),
                                                         mean(p$TNF),median(p$TNF),sum(p$TNF),
                                                         mean(q$IFNG),median(q$IFNG),sum(q$IFNG),
                                                         Treg_frequence,CD8_frequence,CD8_ratio))
  remove(list =paste("Treg_",j,"_TIL",sep = ""))
  remove(list =paste("CD4_",j,"_TIL",sep = ""))
  remove(list =paste("CD8_",j,"_TIL",sep = ""))
  remove(list =paste("pbmc_",j,"_TIL",sep = ""))
  remove(d,e,f,o,p,q,CD8_frequence,CD8_ratio,Treg_frequence)
}
scaldata<-scaldata_1_TIL
if (length(patient)>=2){
  for (j in 2:length(patient)) {
    scaldata<-rbind(scaldata,get(paste("scaldata_",j,"_TIL",sep = "")))
    remove(list = paste("scaldata_",j,"_TIL",sep = ""))
  }
}
remove(scaldata_1_TIL,j,CD4_T,CD8_T,pbmc_T,Treg_T)
colnames(scaldata)<-c("a.i.",                   "length.patient." ,      
                      "mean.d.JMJD1C."  ,       "median.d.JMJD1C." ,     
                      "sum.d.JMJD1C."  ,        "mean.e.IFNG." ,         
                      "median.e.IFNG."   ,      "sum.e.IFNG."  ,         
                      "mean.f.STAT3_mean."  ,   "median.f.STAT3_mean.",  
                      "sum.f.STAT3_mean."  ,    "mean.f.STAT3_median.",  
                      "median.f.STAT3_median." ,"sum.f.STAT3_median." ,  
                      "mean.f.STAT3_sum.",      "median.f.STAT3_sum.",   
                      "sum.f.STAT3_sum.",       "mean.f.PI3K_mean." ,    
                      "median.f.PI3K_mean."  ,  "sum.f.PI3K_mean."   ,   
                      "mean.f.PI3K_median."  ,  "median.f.PI3K_median." ,
                      "sum.f.PI3K_median."   ,  "mean.f.PI3K_sum."  ,    
                      "median.f.PI3K_sum."   ,  "sum.f.PI3K_sum."   ,    
                      "mean.o.GZMB."         ,  "median.o.GZMB."     ,   
                      "sum.o.GZMB."          ,  "mean.p.TNF."         ,  
                      "median.p.TNF."        ,  "sum.p.TNF."           , 
                      "mean.q.IFNG."         ,  "median.q.IFNG."       , 
                      "sum.q.IFNG."          ,  "Treg_frequence"       , 
                      "CD8_frequence"        ,  "CD8_ratio" )
write.table(scaldata,"c:/Users/xjmik/Desktop/scaledata2_patientnoise.txt",sep = "\t")