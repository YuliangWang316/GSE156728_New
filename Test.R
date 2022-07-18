library(dplyr)
a<-read.table("c:/Users/xjmik/Documents/Rawdata&scaledata3.txt",sep = "\t",header = TRUE,row.names = 1)
a<-a[-1,]
b<-unique(a$a.i.)
c<-b[11:15]
e<-a[which(a$a.i.== c[1]),]
for (k in 2:length(c)) {
  d<-a[which(a$a.i. == c[k],),]
  e<-rbind(e,d)
}     
remove(d,k)
library(ggplot2)
library(ggpubr)
ggplot(data = e,aes(x=sum.d.JMJD1C.,median.e.IFNG.)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=TRUE,size=1.5,color="red")+stat_cor(data = e,method = "spearman") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 

write.table(e,file = "4-9_scale.txt",sep = "\t")

