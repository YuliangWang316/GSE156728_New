library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
setwd("c:/Users/xjmik/Downloads/data.expression/data/expression/CD4/byDataset/")
a<-list.files("c:/Users/xjmik/Downloads/data.expression/data/expression/CD4/byDataset/")
count<-data.frame()
metadata<-data.frame()
gene<-character()
c<-readRDS(a[1])
count<-c@assays$data@listData$norm_exprs
count<-as.matrix(count)
count<-as.data.frame(count)
metadata<-as.data.frame(c@colData@listData)
count<-count[,metadata$cellID]
colnames(count)<-metadata$cellID.uniq
gene<-rownames(count)
for (i in c(2:23,25:length(a))) {
  b<-readRDS(a[i])
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
    genename<-e$SYMBOL
  }
  f<-sort(genename)
  if(f[1] == "1"){
    g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
    g<-g[!duplicated(g$SYMBOL),]
    count_B<-count_B[g$ENTREZID,]
    rownames(count_B)<-g$SYMBOL
    genename<-g$SYMBOL
  }
  gene<-intersect(gene,genename)
  count_B<-count_B[gene,]
  count<-count[gene,]
  count<-cbind(count,count_B)
  metadata<-rbind(metadata,count_B_metadata)
  remove(count_B,count_B_metadata,b,d)
}
b<-readRDS(a[24])
count_B<-b@assays@data$norm_exprs
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
  genename<-e$SYMBOL
}
f<-sort(genename)
if(f[1] == "1"){
  g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
  g<-g[!duplicated(g$SYMBOL),]
  count_B<-count_B[g$ENTREZID,]
  rownames(count_B)<-g$SYMBOL
  genename<-g$SYMBOL
}
gene<-intersect(gene,genename)
count_B<-count_B[gene,]
count<-count[gene,]
count<-cbind(count,count_B)
metadata<-rbind(metadata,count_B_metadata)
remove(count_B,count_B_metadata,b,d)
remove(c,e,g,a,f,gene,genename,i)
count_CD4<-count
metadata_CD4<-metadata
remove(count,metadata)

library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
setwd("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8/byDataset/")
a<-list.files("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8/byDataset/")
count<-data.frame()
metadata<-data.frame()
gene<-character()
c<-readRDS(a[1])
count<-c@assays$data@listData$norm_exprs
count<-as.matrix(count)
count<-as.data.frame(count)
metadata<-as.data.frame(c@colData@listData)
count<-count[,metadata$cellID]
colnames(count)<-metadata$cellID.uniq
gene<-rownames(count)
for (i in c(2:23,25:length(a))) {
  b<-readRDS(a[i])
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
    genename<-e$SYMBOL
  }
  f<-sort(genename)
  if(f[1] == "1"){
    g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
    g<-g[!duplicated(g$SYMBOL),]
    count_B<-count_B[g$ENTREZID,]
    rownames(count_B)<-g$SYMBOL
    genename<-g$SYMBOL
  }
  gene<-intersect(gene,genename)
  count_B<-count_B[gene,]
  count<-count[gene,]
  count<-cbind(count,count_B)
  metadata<-rbind(metadata,count_B_metadata)
  remove(count_B,count_B_metadata,b,d)
}
b<-readRDS(a[24])
count_B<-b@assays@data$norm_exprs
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
  genename<-e$SYMBOL
}
f<-sort(genename)
if(f[1] == "1"){
  g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
  g<-g[!duplicated(g$SYMBOL),]
  count_B<-count_B[g$ENTREZID,]
  rownames(count_B)<-g$SYMBOL
  genename<-g$SYMBOL
}
gene<-intersect(gene,genename)
count_B<-count_B[gene,]
count<-count[gene,]
count<-cbind(count,count_B)
metadata<-rbind(metadata,count_B_metadata)
remove(count_B,count_B_metadata,b,d)
remove(c,e,g,a,f,gene,genename,i)
count_CD8<-count
metadata_CD8<-metadata
remove(count,metadata)

metadata<-rbind(metadata_CD4,metadata_CD8)
count<-cbind(count_CD4,count_CD8)
e<-readRDS("C:/Users/xjmik/Downloads/data.expression/data/expression/CD4/integration/int.CD4.S35.meta.tb.rds")
f<-readRDS("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8/integration/int.CD8.S35.meta.tb.rds")
g<-rbind(e,f)

remove(count_CD4,count_CD8,metadata_CD4,metadata_CD8,e,f)
h<-intersect(g$cellID.uniq,metadata$cellID.uniq)
g<-as.data.frame(g)
rownames(g)<-g$cellID.uniq
rownames(metadata)<-metadata$cellID.uniq
g_new<-g[h,]
metadata_new<-metadata[h,]
remove(metadata,g,h)
metadata<-cbind(g_new,metadata_new)
remove(g_new,metadata_new)
metadata<-metadata[colnames(count),]
library(Seurat)
pbmc<-CreateSeuratObject(counts = count,meta.data = metadata)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
remove(metadata,count)
#pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 50)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:30)
Idents(pbmc)<-pbmc@meta.data$meta.cluster
VlnPlot(pbmc,features = c("FOXP3"),pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = c("CD8A"),pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = c("CD8B"),pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = c("CD4"),pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = c("IFNG"),pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = c("JMJD1C"),pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = c("IL2RA"),pt.size = 0,sort = TRUE)
FeaturePlot(pbmc,features = c("JMJD1C","FOXP3"),label = TRUE)
DimPlot(pbmc,label = TRUE)
Treg<-subset(pbmc,ident = c("CD4.c18.Treg.RTKN2","CD4.c19.Treg.S1PR1","CD4.c20.Treg.TNFRSF9","CD4.c21.Treg.OAS1","CD4.c23.Mix.NME1","CD4.c24.Mix.NME2"))
Idents(Treg)<-Treg@meta.data$loc
Treg_T<-subset(Treg,ident = "T")
J1C_IFNG<-FetchData(Treg_T,vars = c("JMJD1C","IFNG"))
for (r in 1:length(colnames(J1C_IFNG))) {
  noise<-rnorm(n=length(x=J1C_IFNG[,r]))/1e+05
  J1C_IFNG[,r]<-J1C_IFNG[,r]+noise
}
remove(noise,r)
library(ggplot2)
library(ggpubr)
ggplot(data = J1C_IFNG,aes(x=JMJD1C,IFNG)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=TRUE,size=1.5,color="red")+stat_cor(data = J1C_IFNG,method = "spearman") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
