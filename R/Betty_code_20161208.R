library(dplyr)
library(tidyr)
library(vegan)
library(ecodist)
library(labdsv)
library(rgl)
library(vegan3d)

setwd("L:/Public/Betty/tox_data/Community_20160823/community_20160823")

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





phytoNMS.4<-metaMDS(x[,2:261],k=4,trymax=1000)


 
save(phytoNMS.4,file="phytoNMS4.rda")
#load(file="phytoNMS4.rda")


detection<-x[,272:274]
ef<-envfit(phytoNMS.4,detection,na.rm=TRUE)
ef2<-envfit(phytoNMS.4,x[,275:364],na.rm=TRUE)


fig<-ordiplot(phytoNMS.4,type="none",main="det_sax")
points(fig,"sites",pch=16,cex=0.7,col=detection[,1]+1)
plot(ef)


##Doing enviro fit on other data (Landscape)
landscape<-x[,268:311]
efLandscape<-envfit(phytoNMS.4,landscape,na.rm=TRUE)

fig<-ordiplot(phytoNMS.4,type="none",main="det_sax")
points(fig,"sites",pch=16,cex=0.7,col=detection[,1]+1)
plot(efLandscape)


##Doing enviro fit on other data (water quality)
WQ<-x[,312:356]
efWaterQuality<-envfit(phytoNMS.4,WQ,na.rm=TRUE)

fig<-ordiplot(phytoNMS.4,type="none",main="det_sax")
points(fig,"sites",pch=16,cex=0.7,col=detection[,1]+1)
plot(efWaterQuality)

###Indicator analysis
##group based on NMS1 locations

xLoc=phytoNMS.4$point[,1]
for (i in 1:1249){
        if(xLoc[i] >= 0.2){
                xLoc[i]<-1
        }
        else{
                xLoc[i]<-0
        }
        }

indicator<-indval(x[,2:261],xLoc)
microIndSpec<-indval(x[,2:261],x[,265])


