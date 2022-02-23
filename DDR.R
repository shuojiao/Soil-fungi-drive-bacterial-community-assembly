Distance decay relationship (DDR)

library(vegan)
library(geosphere) 
library(ggplot2)
library(RColorBrewer)

otu_ba<-read.csv("ASV_ba_rep.csv",row.names=1)
fa<-read.csv("factor_rep.csv",row.names=1)
se_ar<-subset(fa,cls=="Agriculture")
se_fo<-subset(fa,cls=="Forest")
se_de<-subset(fa,cls=="Desert")
se_gr<-subset(fa,cls=="Grass")
se_we<-subset(fa,cls=="Wetland")

beta_ba<-as.matrix(vegdist(otu_ba))
beta<-as.dist(beta_ba[rownames(se_ar),rownames(se_ar)])
geo.dist<-as.dist(distm(fa[rownames(se_ar),3:4]))
df1<-as.data.frame(cbind(bray=beta,dist=geo.dist))
df1$p<-"Agriculture"
beta<-as.dist(beta_ba[rownames(se_fo),rownames(se_fo)])
geo.dist<-as.dist(distm(fa[rownames(se_fo),3:4]))
df2<-as.data.frame(cbind(bray=beta,dist=geo.dist))
df2$p<-"Forest"
beta<-as.dist(beta_ba[rownames(se_we),rownames(se_we)])
geo.dist<-as.dist(distm(fa[rownames(se_we),3:4]))
df3<-as.data.frame(cbind(bray=beta,dist=geo.dist))
df3$p<-"Wetland"
beta<-as.dist(beta_ba[rownames(se_de),rownames(se_de)])
geo.dist<-as.dist(distm(fa[rownames(se_de),3:4]))
df4<-as.data.frame(cbind(bray=beta,dist=geo.dist))
df4$p<-"Desert"
beta<-as.dist(beta_ba[rownames(se_gr),rownames(se_gr)])
geo.dist<-as.dist(distm(fa[rownames(se_gr),3:4]))
df5<-as.data.frame(cbind(bray=beta,dist=geo.dist))
df5$p<-"Grass"
df<-rbind(df1,df2,df3,df4,df5)
df$p<-factor(df$p,levels=c("Agriculture","Forest","Wetland","Grass","Desert"))

col1<-brewer.pal(8,"Dark2")
pdf("DDR.pdf",width=4.4,height=10)
p<-ggplot(df,aes(dist/10000,100-bray*100,color=p))
p + 
    geom_point(size=1.5,alpha=0.3)+
    stat_smooth(method ='lm',color="black")+
    facet_wrap(~p,nrow=5)+
    scale_color_manual(values=col1)+
    xlab("Geographic distance (10 km)")+ylab("Community similarity (%)")+
    theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()



