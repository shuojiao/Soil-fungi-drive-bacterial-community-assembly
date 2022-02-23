library(ggtree)
library(picante)
library(gplots)
library(reshape2)

phylo_ba<-read.tree("ASVs_ba.tre")
dat<-read.csv("domin_ASV_rep",row.names=1)
lm_rep<-read.csv("lm_rep.csv",row.names=1)

c<-match.phylo.data(phylo_ba,t(dat))
aa<-c$phy
dd<-data.frame(id=aa$tip.label,otuname_ba[aa$tip.label,])
dd$aa<-ifelse(dd$Genus=="",rownames(dd),as.vector(dd$Genus))
bb<-otuname_ba[aa$tip.label,]
bb$cc<-as.vector(bb$Phylum)
groupInfo1 <- split(aa$tip.label, as.vector(bb$cc))
aa<- groupOTU(aa, groupInfo1)
heat<-lm_rep
for(i in 1:ncol(heat)){
heat[,i]<-(heat[,i]-min(heat[,i]))/(max(heat[,i])-min(heat[,i]))
}
heat<-t(heat)
c<-t(heat)
for(i in 1:ncol(c)){
c[,i]<-rank(b[,i],ties.method ="first")}
a1<-colnames(c[,c[14,]==17])
a2<-colnames(c[,c[15,]==17])
a3<-colnames(c[,c[16,]==17])
a4<-colnames(c[,c[17,]==17])
a<-c(a1,a2,a3,a4)
d1<-data.frame(ASV=rownames(heat),value=0.0001)
d1[d1$ASV%in%a,]$value<-1
d1$lab<-dd$aa
p<-ggtree(aa,aes(color=group),layout ="fan",branch.length='none',open.angle=90,lwd=0.5)
p<-p %<+% d1+geom_tippoint(aes(size=value))+
        geom_tiplab2(aes(label=lab),offset=19,align=TRUE,linetype=NA,size=4,linesize=0)
pdf("tree_lm.pdf",width=16,height=16)
gheatmap(p,heat,offset=0.4, width=0.65, font.size=4,hjust=0,colnames_angle=0,colnames = T,colnames_position = "top")+
scale_fill_gradientn(colours=colorpanel(100,low="white",high="darkgreen"))+
scale_size_continuous(range=c(0,1.5))
dev.off()
