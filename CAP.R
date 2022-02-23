library(vegan)
library(ggplot2)

otu_ba<-read.csv("ASV_ba_rep.csv",row.names=1)
fa<-read.csv("factor_rep.csv",row.names=1)
env<-read.csv("env_rep.csv",row.names=1)
se_ar<-subset(fa,cls=="Agriculture")
se_fo<-subset(fa,cls=="Forest")
se_de<-subset(fa,cls=="Desert")
se_gr<-subset(fa,cls=="Grass")
se_we<-subset(fa,cls=="Wetland")

beta_ba<-as.matrix(vegdist(otu_ba))
beta<-as.dist(beta_ba[rownames(se_ar),rownames(se_ar)])
dat<-env[rownames(se_ar),]
mod1<-capscale(beta~.,dat,add=T)
mod0<-capscale(beta~1,dat,add=T)
mod<-ordiR2step(mod0,scope=formula(mod1),perm.max=999)
a<-anova(mod,by="term")
anova(mod)
b<-anova(mod,by="axis")
si<-scores(mod)$sites
sp<-as.data.frame(mod$CCA$biplot/1)
df<-data.frame(si,dat)
pdf("CAP_ba_ar.pdf",width=5,height=4)
p<-ggplot(df,aes(CAP1,CAP2))
p + 
    geom_point(size=3,aes(colour=AD))+
    geom_segment(data=sp,aes(x=0, xend=CAP1, y=0, yend=CAP2),color="gray50",lwd=0.5,
                 arrow=arrow(angle=20,length=unit(0.02,"npc"),type="closed"))+
    geom_text(data=sp,aes(x=CAP1,y=CAP2,label=rownames(sp),
              hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))), 
              color="black",size=5)+
    scale_colour_gradient(limits=c(1,14),low="red",high="blue")+
    labs(x="CAP1",y="CAP2")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()
write.table(a,"sp_ba_ar.xls",sep="\t")
write.table(b,"axis_ba_ar.xls",sep="\t")

beta<-as.dist(beta_ba[rownames(se_fo),rownames(se_fo)])
dat<-env[rownames(se_fo),]
mod1<-capscale(beta~.,dat,add=T)
mod0<-capscale(beta~1,dat,add=T)
mod<-ordiR2step(mod0,scope=formula(mod1),perm.max=999)
a<-anova(mod,by="term")
anova(mod)
b<-anova(mod,by="axis")
si<-scores(mod)$sites
sp<-as.data.frame(mod$CCA$biplot/1)
df<-data.frame(si,dat)
pdf("CAP_ba_fo.pdf",width=5,height=4)
p<-ggplot(df,aes(CAP1,CAP2))
p + 
    geom_point(size=3,aes(colour=AD))+
    geom_segment(data=sp,aes(x=0, xend=CAP1, y=0, yend=CAP2),color="gray50",lwd=0.5,
                 arrow=arrow(angle=20,length=unit(0.02,"npc"),type="closed"))+
    geom_text(data=sp,aes(x=CAP1,y=CAP2,label=rownames(sp),
              hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))), 
              color="black",size=5)+
    scale_colour_gradient(limits=c(1,14),low="red",high="blue")+
    labs(x="CAP1",y="CAP2")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()
write.table(a,"sp_ba_fo.xls",sep="\t")
write.table(b,"axis_ba_fo.xls",sep="\t")

beta<-as.dist(beta_ba[rownames(se_we),rownames(se_we)])
dat<-env[rownames(se_we),]
mod1<-capscale(beta~.,dat,add=T)
mod0<-capscale(beta~1,dat,add=T)
mod<-ordiR2step(mod0,scope=formula(mod1),perm.max=999)
a<-anova(mod,by="term")
anova(mod)
b<-anova(mod,by="axis")
si<-scores(mod)$sites
sp<-as.data.frame(mod$CCA$biplot/1)
df<-data.frame(si,dat)
pdf("CAP_ba_we.pdf",width=5,height=4)
p<-ggplot(df,aes(CAP1,CAP2))
p + 
    geom_point(size=3,aes(colour=AD))+
    geom_segment(data=sp,aes(x=0, xend=CAP1, y=0, yend=CAP2),color="gray50",lwd=0.5,
                 arrow=arrow(angle=20,length=unit(0.02,"npc"),type="closed"))+
    geom_text(data=sp,aes(x=CAP1,y=CAP2,label=rownames(sp),
              hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))), 
              color="black",size=5)+
    scale_colour_gradient(limits=c(1,14),low="red",high="blue")+
    labs(x="CAP1",y="CAP2")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()
write.table(a,"sp_ba_we.xls",sep="\t")
write.table(b,"axis_ba_we.xls",sep="\t")

beta<-as.dist(beta_ba[rownames(se_gr),rownames(se_gr)])
dat<-env[rownames(se_gr),]
mod1<-capscale(beta~.,dat,add=T)
mod0<-capscale(beta~1,dat,add=T)
mod<-ordiR2step(mod0,scope=formula(mod1),perm.max=999)
a<-anova(mod,by="term")
anova(mod)
b<-anova(mod,by="axis")
si<-scores(mod)$sites
sp<-as.data.frame(mod$CCA$biplot)
df<-data.frame(si,dat)
pdf("CAP_ba_gr.pdf",width=5,height=4)
p<-ggplot(df,aes(CAP1,CAP2))
p + 
    geom_point(size=3,aes(colour=AD))+
    geom_segment(data=sp,aes(x=0, xend=CAP1, y=0, yend=CAP2),color="gray50",lwd=0.5,
                 arrow=arrow(angle=20,length=unit(0.02,"npc"),type="closed"))+
    geom_text(data=sp,aes(x=CAP1,y=CAP2,label=rownames(sp),
              hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))), 
              color="black",size=5)+
    scale_colour_gradient(limits=c(1,14),low="red",high="blue")+
    labs(x="CAP1",y="CAP2")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()
write.table(a,"sp_ba_gr.xls",sep="\t")
write.table(b,"axis_ba_gr.xls",sep="\t")

beta<-as.dist(beta_ba[rownames(se_de),rownames(se_de)])
dat<-env[rownames(se_de),]
mod1<-capscale(beta~.,dat,add=T)
mod0<-capscale(beta~1,dat,add=T)
mod<-ordiR2step(mod0,scope=formula(mod1),perm.max=999)
a<-anova(mod,by="term")
anova(mod)
b<-anova(mod,by="axis")
si<-scores(mod)$sites
sp<-as.data.frame(mod$CCA$biplot)
df<-data.frame(si,dat)
pdf("CAP_ba_de.pdf",width=5,height=4)
p<-ggplot(df,aes(CAP1,CAP2))
p + 
    geom_point(size=3,aes(colour=AD))+
    geom_segment(data=sp,aes(x=0, xend=CAP1, y=0, yend=CAP2),color="gray50",lwd=0.5,
                 arrow=arrow(angle=20,length=unit(0.02,"npc"),type="closed"))+
    geom_text(data=sp,aes(x=CAP1,y=CAP2,label=rownames(sp),
              hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))), 
              color="black",size=5)+
    scale_colour_gradient(limits=c(1,14),low="blue",high="red")+
    labs(x="CAP1",y="CAP2")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()
write.table(a,"sp_ba_de.xls",sep="\t")
write.table(b,"axis_ba_de.xls",sep="\t")
