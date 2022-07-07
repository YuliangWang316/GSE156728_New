library(dplyr)
a<-read.table("c:/Users/xjmik/Desktop/Rawdata3.txt",sep = "\t",header = TRUE,row.names = 1)
a<-a[-1,]
b<-unique(a$a.i.)
q<-data.frame(0,0,0,0,0,0)
colnames(q)<-c("i","j","d","f","o","n")
for (i in 1:34) {
  for (j in 2:35) {
    if(i<j){
      c<-b[i:j]
      e<-a[which(a$a.i.== c[1]),]
      for (k in 2:length(c)) {
        d<-a[which(a$a.i. == c[k],),]
        e<-rbind(e,d)
      }
      remove(d,k)
      for (d in 3:5) {
        for (f in 6:8) {
          n<-cor.test(e[,d],e[,f],method = "spearman",exact = TRUE)
          o<-cor(e[,d],e[,f],method = "spearman")
          p<-data.frame(i,j,d,f,o,n$p.value)
          colnames(p)<-c("i","j","d","f","o","n")
          q<-rbind(q,p)
          remove(n,o,p)
        }
        for (g in 9:17) {
          n<-cor.test(e[,d],e[,g],method = "spearman",exact = TRUE)
          o<-cor(e[,d],e[,g],method = "spearman")
          p<-data.frame(i,j,d,g,o,n$p.value)
          colnames(p)<-c("i","j","d","f","o","n")
          q<-rbind(q,p)
          remove(n,o,p)
        }
        for (h in 18:26) {
          n<-cor.test(e[,d],e[,h],method = "spearman",exact = TRUE)
          o<-cor(e[,d],e[,h],method = "spearman")
          p<-data.frame(i,j,d,h,o,n$p.value)
          colnames(p)<-c("i","j","d","f","o","n")
          q<-rbind(q,p)
          remove(n,o,p)
        }
        for (k in 27:29) {
          n<-cor.test(e[,d],e[,k],method = "spearman",exact = TRUE)
          o<-cor(e[,d],e[,k],method = "spearman")
          p<-data.frame(i,j,d,k,o,n$p.value)
          colnames(p)<-c("i","j","d","f","o","n")
          q<-rbind(q,p)
          remove(n,o,p)
        }
        for (l in 30:32) {
          n<-cor.test(e[,d],e[,l],method = "spearman",exact = TRUE)
          o<-cor(e[,d],e[,l],method = "spearman")
          p<-data.frame(i,j,d,l,o,n$p.value)
          colnames(p)<-c("i","j","d","f","o","n")
          q<-rbind(q,p)
          remove(n,o,p)
        }
        for (m in 33:35) {
          n<-cor.test(e[,d],e[,m],method = "spearman",exact = TRUE)
          o<-cor(e[,d],e[,m],method = "spearman")
          p<-data.frame(i,j,d,m,o,n$p.value)
          colnames(p)<-c("i","j","d","f","o","n")
          q<-rbind(q,p)
          remove(n,o,p)
        }
        n<-cor.test(e[,d],e[,36],method = "spearman",exact = TRUE)
        o<-cor(e[,d],e[,36],method = "spearman")
        p<-data.frame(i,j,d,36,o,n$p.value)
        colnames(p)<-c("i","j","d","f","o","n")
        q<-rbind(q,p)
        remove(n,o,p)
        n<-cor.test(e[,d],e[,37],method = "spearman",exact = TRUE)
        o<-cor(e[,d],e[,37],method = "spearman")
        p<-data.frame(i,j,d,37,o,n$p.value)
        colnames(p)<-c("i","j","d","f","o","n")
        q<-rbind(q,p)
        remove(n,o,p)
        n<-cor.test(e[,d],e[,38],method = "spearman",exact = TRUE)
        o<-cor(e[,d],e[,38],method = "spearman")
        p<-data.frame(i,j,d,38,o,n$p.value)
        colnames(p)<-c("i","j","d","f","o","n")
        q<-rbind(q,p)
        remove(n,o,p)
      }
      remove(d,f,g,h,k,l,m)
      for (f in 6:8) {
        for (g in 9:17) {
          n<-cor.test(e[,f],e[,g],method = "spearman",exact = TRUE)
          o<-cor(e[,f],e[,g],method = "spearman")
          p<-data.frame(i,j,f,g,o,n$p.value)
          colnames(p)<-c("i","j","d","f","o","n")
          q<-rbind(q,p)
          remove(n,o,p)
        }
        for (h in 18:26) {
          n<-cor.test(e[,f],e[,h],method = "spearman",exact = TRUE)
          o<-cor(e[,f],e[,h],method = "spearman")
          p<-data.frame(i,j,f,h,o,n$p.value)
          colnames(p)<-c("i","j","d","f","o","n")
          q<-rbind(q,p)
          remove(n,o,p)
        }
      }
      remove(f,g,h,c,e)
    }
  }
}
write.table(q,"newcombo.txt",sep = "\t")
remove(q,i,j)
q<-data.frame(0,0,0,0)
colnames(q)<-c("d","f","o","n")
for (i in 1:35) {
  e<-a[which(a$a.i.== b[1]),]
  for (d in 3:5) {
    for (f in 6:8) {
      n<-cor.test(e[,d],e[,f],method = "spearman",exact = TRUE)
      o<-cor(e[,d],e[,f],method = "spearman")
      p<-data.frame(d,f,o,n$p.value)
      colnames(p)<-c("d","f","o","n")
      q<-rbind(q,p)
      remove(n,o,p)
    }
    for (g in 9:17) {
      n<-cor.test(e[,d],e[,g],method = "spearman",exact = TRUE)
      o<-cor(e[,d],e[,g],method = "spearman")
      p<-data.frame(d,g,o,n$p.value)
      colnames(p)<-c("d","f","o","n")
      q<-rbind(q,p)
      remove(n,o,p)
    }
    for (h in 18:26) {
      n<-cor.test(e[,d],e[,h],method = "spearman",exact = TRUE)
      o<-cor(e[,d],e[,h],method = "spearman")
      p<-data.frame(d,h,o,n$p.value)
      colnames(p)<-c("d","f","o","n")
      q<-rbind(q,p)
      remove(n,o,p)
    }
    for (k in 27:29) {
      n<-cor.test(e[,d],e[,k],method = "spearman",exact = TRUE)
      o<-cor(e[,d],e[,k],method = "spearman")
      p<-data.frame(d,k,o,n$p.value)
      colnames(p)<-c("d","f","o","n")
      q<-rbind(q,p)
      remove(n,o,p)
    }
    for (l in 30:32) {
      n<-cor.test(e[,d],e[,l],method = "spearman",exact = TRUE)
      o<-cor(e[,d],e[,l],method = "spearman")
      p<-data.frame(d,l,o,n$p.value)
      colnames(p)<-c("d","f","o","n")
      q<-rbind(q,p)
      remove(n,o,p)
    }
    for (m in 33:35) {
      n<-cor.test(e[,d],e[,m],method = "spearman",exact = TRUE)
      o<-cor(e[,d],e[,m],method = "spearman")
      p<-data.frame(d,m,o,n$p.value)
      colnames(p)<-c("d","f","o","n")
      q<-rbind(q,p)
      remove(n,o,p)
    }
    n<-cor.test(e[,d],e[,36],method = "spearman",exact = TRUE)
    o<-cor(e[,d],e[,36],method = "spearman")
    p<-data.frame(d,36,o,n$p.value)
    colnames(p)<-c("d","f","o","n")
    q<-rbind(q,p)
    remove(n,o,p)
    n<-cor.test(e[,d],e[,37],method = "spearman",exact = TRUE)
    o<-cor(e[,d],e[,37],method = "spearman")
    p<-data.frame(d,37,o,n$p.value)
    colnames(p)<-c("d","f","o","n")
    q<-rbind(q,p)
    remove(n,o,p)
    n<-cor.test(e[,d],e[,38],method = "spearman",exact = TRUE)
    o<-cor(e[,d],e[,38],method = "spearman")
    p<-data.frame(d,38,o,n$p.value)
    colnames(p)<-c("d","f","o","n")
    q<-rbind(q,p)
    remove(n,o,p)
  }
  remove(d,f,g,h,k,l,m)
  for (f in 6:8) {
    for (g in 9:17) {
      n<-cor.test(e[,f],e[,g],method = "spearman",exact = TRUE)
      o<-cor(e[,f],e[,g],method = "spearman")
      p<-data.frame(f,g,o,n$p.value)
      colnames(p)<-c("d","f","o","n")
      q<-rbind(q,p)
      remove(n,o,p)
    }
    for (h in 18:26) {
      n<-cor.test(e[,f],e[,h],method = "spearman",exact = TRUE)
      o<-cor(e[,f],e[,h],method = "spearman")
      p<-data.frame(f,h,o,n$p.value)
      colnames(p)<-c("d","f","o","n")
      q<-rbind(q,p)
      remove(n,o,p)
    }
  }
  remove(f,g,h,c,e)
}
write.table(q,"newsingle.txt",sep = "\t")
