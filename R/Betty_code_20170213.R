library(dplyr)
library(tidyr)
library(vegan)
library(ecodist)
library(labdsv)
library(randomForest)

setwd("L:/Public/Betty/tox_data/Community_20170206")


load("phytoNMS4.rda")


load("phyto.rda")
data1<-read.csv(file="ltm_tox_20160810.csv", stringsAsFactors = F)

a<-OTU[,-3]
# ###List based on conversations with Nate Smucker and Bryan Milstead
removeList<-c("Chlorophyta","Chlorophyta coccoid","Chlorophyta Filament",
              "Chrysophyta","Cyanophyta","Diatoms","Diatoms (dead)",
              "Dinoflagellates","Euglenophyte","Heterocontophyta","Stigeoclonium",
              "Unknown alga","Unknown alga flagellated","Cuspidothrix.issatschenkoi",
              "Cuspidothrix issatschenkoi", "Eucocconeis","Psammothidium","Diatoma")

b<-a$OTU %in% removeList
new_a<-a[!b,]

phytoWide<-new_a %>%spread(OTU,sumAbund,0)
#write.csv(file="phytoWide.csv",phytoWide)
#phytoWide<-read.csv(file="phytoWide.csv",header=TRUE, stringsAsFactors = F)



##Need to do envfitting
#rename NLA_ID to SITE_ID in data1

data1<-rename(data1, SITE_ID = NLA_ID)

x<-inner_join(phytoWide,data1,by="SITE_ID")

###Cluster Analysis
nms.dist<-vegdist(phytoNMS.4$points[,1:2],method="euclidean")

nmsclus.com<-hclust(nms.dist,"complete")
nmsclus.sin<-hclust(nms.dist,"single")
nmsclus.ave<-hclust(nms.dist,"average")

cor(nms.dist, cophenetic(nmsclus.com))
#[1] 0.7126077
cor(nms.dist, cophenetic(nmsclus.sin))
#[1] 0.5095576
cor(nms.dist, cophenetic(nmsclus.ave))
#[1] 0.7495803

plot(nmsclus.ave)
rect.hclust(nmsclus.ave,3)

grouping_lake<-cutree(nmsclus.ave,3)

lakes_grouping<-cbind(x,grouping_lake)

write.csv(lakes_grouping,file="grouping.csv")

###Plot of NMS color coded for clusters
png("nms4_cluster.png")
fig<-ordiplot(phytoNMS.4,type="none",main="Final NMS")
points(fig,"sites",pch=16,cex=0.7,col=grouping_lake)
dev.off()

###Env. fit
detection<-lakes_grouping[,264:266]
ef<-envfit(phytoNMS.4,detection,na.rm=TRUE)

png("nms4_ef.png")
fig<-ordiplot(phytoNMS.4,type="none",main="Final NMS")
points(fig,"sites",pch=16,cex=0.7,col=grouping_lake)
plot(ef)
dev.off()


##random forest on NMS clusters

##2015 Model
lakes_grouping[,262]<-factor(lakes_grouping[,262])
#lakes_grouping[,263]<-factor(lakes_grouping[,263])
lakes_grouping[,308]<-factor(lakes_grouping[,308])
predVar<-c(264:268,270,272:287,305:344,346:356)




y.cc<-complete.cases(lakes_grouping[,predVar])
lakes.cc<-lakes_grouping[y.cc,]
a<-factor(lakes.cc[,357],labels=c("a","b","c"))
#x<-lakes.cc[,predVar]

lakesRF<-randomForest(lakes.cc[,predVar],a,ntree=10000,importance=TRUE,proximity=TRUE)


###Indicator Species Analysis
lakes_ind<-indval(lakes_grouping[,2:261],lakes_grouping[,357])
write.csv(lakes_ind$indval,file="lakes_indicate.csv")

