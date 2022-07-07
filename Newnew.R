a<-read.table("c:/Users/xjmik/Documents/newcombo.txt",sep = "\t",header = TRUE,row.names = 1)
a<-a[-1,]
a<-a[which(a$n != "NA"),]
a<-a[which(a$n<0.05),]
a$p<-paste0(a$i,"-",a$j)
CD8ratio<-a[which(a$f=="38" & a$o<0),]
CD8<-a[which(a$f=="37"& a$o<0),]
Treg<-a[which(a$f=="36"& a$o>0),]

Treg_CD8ratio<-intersect(Treg$p,CD8ratio$p)
Treg_CD8<-intersect(Treg$p,CD8$p)


Treg_1c_mean<-a[which(a$d == "3"  & a$o <0),]
Treg_1c_mean_IFNG<-Treg_1c_mean[which(Treg_1c_mean$f == "6" | Treg_1c_mean$f == "7" | Treg_1c_mean$f == "8"),]
Treg_1c_mean_STAT3<-Treg_1c_mean[which(Treg_1c_mean$f == "9" | Treg_1c_mean$f == "10" | Treg_1c_mean$f == "11" | Treg_1c_mean$f == "12" | Treg_1c_mean$f == "13" | Treg_1c_mean$f == "14" | Treg_1c_mean$f == "15" | Treg_1c_mean$f == "16" | Treg_1c_mean$f == "17"),]
Treg_1c_mean_PI3K<-Treg_1c_mean[which(Treg_1c_mean$f == "18" | Treg_1c_mean$f == "19" | Treg_1c_mean$f == "20" | Treg_1c_mean$f == "21" | Treg_1c_mean$f == "22" | Treg_1c_mean$f == "23" | Treg_1c_mean$f == "24" | Treg_1c_mean$f == "25" | Treg_1c_mean$f == "26"),]
Treg_1c_mean_CD8GZMB_TNF_IFNG<-Treg_1c_mean[which(Treg_1c_mean$f == "27" | Treg_1c_mean$f == "28" | Treg_1c_mean$f == "29" | Treg_1c_mean$f == "30" | Treg_1c_mean$f == "31" | Treg_1c_mean$f == "32" | Treg_1c_mean$f == "33" | Treg_1c_mean$f == "34" | Treg_1c_mean$f == "35" ),]
PI3K_STAT3<-intersect(Treg_1c_mean_PI3K$p,Treg_1c_mean_STAT3$p)
IFNG_PI3K_STAT3<-intersect(Treg_1c_mean_IFNG$p,PI3K_STAT3)
Treg_CD8_IFNG_PI3K_STAT3<-intersect(Treg_CD8,IFNG_PI3K_STAT3)
Treg_CD8ratio_IFNG_PI3K_STAT3<-intersect(Treg_CD8ratio,IFNG_PI3K_STAT3)
Treg_CD8_IFNG_PI3K_STAT3_CD8GZMB_TNF_IFNG<-intersect(Treg_CD8_IFNG_PI3K_STAT3,Treg_1c_mean_CD8GZMB_TNF_IFNG)
Treg_CD8ratio_IFNG_PI3K_STAT3_CD8GZMB_TNF_IFNG<-intersect(Treg_CD8ratio_IFNG_PI3K_STAT3,Treg_1c_mean_CD8GZMB_TNF_IFNG)

Treg_1c_median<-a[which(a$d == "4"  & a$o <0),]
Treg_1c_median_IFNG<-Treg_1c_median[which(Treg_1c_median$f == "6" | Treg_1c_median$f == "7" | Treg_1c_median$f == "8"),]
Treg_1c_median_STAT3<-Treg_1c_median[which(Treg_1c_median$f == "9" | Treg_1c_median$f == "10" | Treg_1c_median$f == "11" | Treg_1c_median$f == "12" | Treg_1c_median$f == "13" | Treg_1c_median$f == "14" | Treg_1c_median$f == "15" | Treg_1c_median$f == "16" | Treg_1c_median$f == "17"),]
Treg_1c_median_PI3K<-Treg_1c_median[which(Treg_1c_median$f == "18" | Treg_1c_median$f == "19" | Treg_1c_median$f == "20" | Treg_1c_median$f == "21" | Treg_1c_median$f == "22" | Treg_1c_median$f == "23" | Treg_1c_median$f == "24" | Treg_1c_median$f == "25" | Treg_1c_median$f == "26"),]
Treg_1c_median_CD8GZMB_TNF_IFNG<-Treg_1c_median[which(Treg_1c_median$f == "27" | Treg_1c_median$f == "28" | Treg_1c_median$f == "29" | Treg_1c_median$f == "30" | Treg_1c_median$f == "31" | Treg_1c_median$f == "32" | Treg_1c_median$f == "33" | Treg_1c_median$f == "34" | Treg_1c_median$f == "35" ),]
PI3K_STAT3<-intersect(Treg_1c_median_PI3K$p,Treg_1c_median_STAT3$p)
IFNG_PI3K_STAT3<-intersect(Treg_1c_median_IFNG$p,PI3K_STAT3)
Treg_CD8_IFNG_PI3K_STAT3<-intersect(Treg_CD8,IFNG_PI3K_STAT3)
Treg_CD8ratio_IFNG_PI3K_STAT3<-intersect(Treg_CD8ratio,IFNG_PI3K_STAT3)

Treg_1c_sum<-a[which(a$d == "5"  & a$o <0),]
Treg_1c_sum_IFNG<-Treg_1c_sum[which(Treg_1c_sum$f == "6" | Treg_1c_sum$f == "7" | Treg_1c_sum$f == "8"),]
Treg_1c_sum_STAT3<-Treg_1c_sum[which(Treg_1c_sum$f == "9" | Treg_1c_sum$f == "10" | Treg_1c_sum$f == "11" | Treg_1c_sum$f == "12" | Treg_1c_sum$f == "13" | Treg_1c_sum$f == "14" | Treg_1c_sum$f == "15" | Treg_1c_sum$f == "16" | Treg_1c_sum$f == "17"),]
Treg_1c_sum_PI3K<-Treg_1c_sum[which(Treg_1c_sum$f == "18" | Treg_1c_sum$f == "19" | Treg_1c_sum$f == "20" | Treg_1c_sum$f == "21" | Treg_1c_sum$f == "22" | Treg_1c_sum$f == "23" | Treg_1c_sum$f == "24" | Treg_1c_sum$f == "25" | Treg_1c_sum$f == "26"),]
Treg_1c_sum_CD8GZMB_TNF_IFNG<-Treg_1c_sum[which(Treg_1c_sum$f == "27" | Treg_1c_sum$f == "28" | Treg_1c_sum$f == "29" | Treg_1c_sum$f == "30" | Treg_1c_sum$f == "31" | Treg_1c_sum$f == "32" | Treg_1c_sum$f == "33" | Treg_1c_sum$f == "34" | Treg_1c_sum$f == "35" ),]
PI3K_STAT3<-intersect(Treg_1c_sum_PI3K$p,Treg_1c_sum_STAT3$p)
IFNG_PI3K_STAT3<-intersect(Treg_1c_sum_IFNG$p,PI3K_STAT3)
Treg_CD8_IFNG_PI3K_STAT3<-intersect(Treg_CD8,IFNG_PI3K_STAT3)
Treg_CD8ratio_IFNG_PI3K_STAT3<-intersect(Treg_CD8ratio,IFNG_PI3K_STAT3)
Treg_CD8_IFNG_PI3K_STAT3_CD8GZMB_TNF_IFNG<-intersect(Treg_CD8_IFNG_PI3K_STAT3,Treg_1c_mean_CD8GZMB_TNF_IFNG)
Treg_CD8ratio_IFNG_PI3K_STAT3_CD8GZMB_TNF_IFNG<-intersect(Treg_CD8ratio_IFNG_PI3K_STAT3,Treg_1c_mean_CD8GZMB_TNF_IFNG)
