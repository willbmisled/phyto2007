library(dplyr)
library(tidyr)
library(vegan)
library(ecodist)

setwd("L:/Public/Betty/tox_data/Community_20160823/community_20160823")

load("phyto.rda")
data1<-read.csv(file="ltm_tox_20160810.csv")


a<-OTU[,-3]


# ###List based on conversations with Nate Smucker and Bryan Milstead
# removeList<-c("Chlorophyta","Chlorophyta coccoid","Chlorophyta Filament",
#               "Chrysophyta","Cyanophyta","Diatoms","Diatoms (dead)",
#               "Dinoflagellates","Euglenophyte","Heterocontophyta","Stigeoclonium",
#               "Unknown alga","Unknown alga flagellated")
# 
# 
# b<-a[-grep(removeList,a$OTU),]
# b<-a$OTU %in% removeList
# new_a<-a[!b,]
# 
# phytoWide<-new_a %>%spread(OTU,sumAbund,0)
phytoWide<-read.csv(file="phytoWide.csv",header=TRUE)

##A species distance matrix using Bray Curtis
phytoDist<-distance(phytoWide[,6:269],"bray")

##Cluster Analysis
phytoClusS<-hclust(phytoDist,"single")
phytoClusA<-hclust(phytoDist,"average")
phytoClusC<-hclust(phytoDist,"complete")

cor(phytoDist, cophenetic(phytoClusS))
#[1] 0.5210736
cor(phytoDist, cophenetic(phytoClusA))
#[1] 0.7986093
cor(phytoDist, cophenetic(phytoClusC))
#[1] 0.6341144


##Looking at it plotted based on averaging
plot(phytoClusA)
rect.hclust(phytoClusA,20)



##Running NMS fit
#phytoNMS.2<-metaMDS(phytoWide[,-1],k=2,trymax=500)
#phytoNMS.3<-metaMDS(phytoWide[,-1],k=3,trymax=500)
#phytoNMS.4<-metaMDS(phytoWide[,6:269],k=4,trymax=1000)

#phytoNMS.5<-metaMDS(phytoWide[,-1],k=5,trymax=500)
#phytoNMS.6<-metaMDS(phytoWide[,-1],k=6,trymax=500)

#save(phytoNMS.4,file="phytoNMS4.rda")
#load(file="phytoNMS4.rda")

fig<-ordiplot(phytoNMS.4,type="none")
#plot(fig$sites,cex=0.0,axes=FALSE,xlab="")
points(fig,"sites",pch=1,cex=0.5,col=phytoWide[,5])



##Need to do envfitting
#rename NLA_ID to SITE_ID in data1

data1<-rename(data1, SITE_ID = NLA_ID)
x<-inner_join(phytoWide,data1,by="SITE_ID")

phytoNMS.4<-metaMDS(x[,6:269],k=4,trymax=1000)

#save(phytoNMS.4,file="phytoNMS4.rda")
#load(file="phytoNMS4.rda")

detection<-x[,270:274]
ef<-envfit(phytoNMS.4,detection[,3:5],na.rm=TRUE)
ef2<-envfit(phytoNMS.4,x[,270:364],na.rm=TRUE)


fig<-ordiplot(phytoNMS.4,type="none",main="det_sax")
points(fig,"sites",pch=16,cex=0.7,col=detection[,5]+1)
plot(ef)

