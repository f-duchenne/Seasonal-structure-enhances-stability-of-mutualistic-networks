library(data.table)
library(dplyr)
library(igraph)
library(ggplot2)
library(cowplot)
library(ggplotify)
library(gridExtra)
library(viridis)
library(scales)
library(lme4)
library(nlme)
library(piecewiseSEM)
library(qgraph)
library(igraph)
library(glmmTMB)
library(car)
library(randomForest)
library(directlabels)
library(metR)
library(reshape2)
library(ggcorrplot)
library(e1071) 
library(ggalt)
library(ggExtra)
library(DHARMa)
library(tidyr)
library(ggpubr)
setwd(dir="C:/Users/Duchenne/Documents/EPHI/temporal-dynamics/data and results")
#charge simulations results
dat <- fread("results_simulations_network_level.txt")
dat$model2=factor(dat$model2,c("without seasonal structure","with seasonal structure"))



#### BENEFITS FOR SPECIES (Fig. S3)
library(RColorBrewer)
colo=brewer.pal(11,'PiYG')[c(1,2,4,8,10,11)]
bidon=subset(dat2f,N_eq1>0)
png("figure_S3.png",width=1800,height=1000,res=150)
ggplot(data=bidon,aes(x=func_res,fill=as.factor(efficience),col=as.factor(efficience)))+geom_density(alpha=0.4,trim = TRUE)+facet_wrap(~paste0("c = ",interf),scales="free")+
scale_fill_manual(values=colo)+
scale_color_manual(values=colo)+geom_vline(aes(xintercept=efficience/0.8,color=as.factor(efficience)),linetype="dashed")+xlab("Species benefit from mutualistic interactions")+
labs(fill=expression(alpha),col=expression(alpha))+theme(panel.grid=element_blank())
dev.off();

##### STABILITY (FIG. 2):
seuil=0.98
b=dat %>% dplyr::group_by(interf,efficience,model2) %>% dplyr::summarise(moy=length(pers[pers>seuil])/(250*11))

pl1=ggplot(data = b, aes(x=efficience, y=interf, fill=moy,z=moy))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),
plot.subtitle=element_text(size=13),axis.title=element_text(size=14,color="white"))+
scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),
labels = percent,n.breaks=3)+labs(fill='')+ggtitle("a",subtitle ="Network feasibility")+
xlab("")+
ylab("")+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(model2))


b1=subset(dat,model=="model1") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=length(pers[pers>seuil])/(250*11))
b2=subset(dat,model=="model2") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=length(pers[pers>seuil])/(250*11))
b1$moy[is.na(b1$moy)]=0
b2$moy[is.na(b2$moy)]=0
b1$diffp=(b2$moy-b1$moy)/b1$moy
b1$diffp[is.na(b1$diffp)]=0
mid=abs(min(b1$diffp))/(max(b1$diffp)-min(b1$diffp))
b1$legende=""
fwrite(b1,"changes_in_feasibility.txt")


pl2=ggplot(data = b1, aes(x=efficience, y=interf, fill=diffp,z=diffp))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),
plot.subtitle=element_text(size=13),axis.title=element_text(size=14))+
#scale_fill_gradientn(colors=c("darkred","pink","white","dodgerblue"),values=c(0,mid-0.1,mid,1),labels = percent,n.breaks=3)+
scale_fill_gradient2(low ="darkred",high="dodgerblue4",
labels = percent,n.breaks=4)+
labs(fill='')+ggtitle("b",subtitle ="Changes in feasibility")+
xlab("")+
ylab("")+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(legende))

b=dat %>% dplyr::group_by(interf,efficience,model2) %>% dplyr::summarise(moy=mean(pers))

pl3=ggplot(data = b, aes(x=efficience, y=interf, fill=moy,z=moy))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),
plot.subtitle=element_text(size=13),axis.title=element_text(size=14))+
scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),
labels = percent,n.breaks=3)+labs(fill='')+ggtitle("c",subtitle ="Network persistence")+
xlab("")+
ylab("\n")+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(model2))


b1=subset(dat,model=="model1") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(pers))
b2=subset(dat,model=="model2") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(pers))
b1$moy[is.na(b1$moy)]=0
b2$moy[is.na(b2$moy)]=0
b1$diffp=(b2$moy-b1$moy)/b1$moy
b1$diffp[is.na(b1$diffp)]=0
mid=abs(min(b1$diffp))/(max(b1$diffp)-min(b1$diffp))
b1$legende=""
fwrite(b1,"changes_in_persistence.txt")

pl4=ggplot(data = b1, aes(x=efficience, y=interf, fill=diffp,z=diffp))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),
plot.subtitle=element_text(size=13),axis.title=element_text(size=14))+
#scale_fill_gradientn(colors=c("darkred","pink","white","dodgerblue"),values=c(0,mid-0.1,mid,1),labels = percent,n.breaks=3)+
scale_fill_gradient2(low ="darkred",high="dodgerblue4",
labels = percent,n.breaks=4)+
labs(fill='')+ggtitle("d",subtitle ="Changes in persistence")+
xlab("")+
ylab("")+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(legende))

b=dat %>% dplyr::group_by(interf,efficience,model2) %>% dplyr::summarise(moy=log(mean(stab,na.rm=T)))

pl5=ggplot(data = b, aes(x=efficience, y=interf, fill=moy,z=moy))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=13))+
scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))))+
labs(fill='')+ggtitle("e",subtitle ="Network resilience")+
xlab("")+
ylab(expression("Competition for mutualistic partners  "(c)))+coord_cartesian(expand = FALSE)+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(model2))

b1=subset(dat,model=="model1") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(stab,na.rm=T))
b2=subset(dat,model=="model2") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(stab,na.rm=T))
b1$diffp=(b2$moy-b1$moy)/b1$moy
mid=abs(min(b1$diffp,na.rm=T))/(max(b1$diffp,na.rm=T)-min(b1$diffp,na.rm=T))
b1$legende=""
fwrite(b1,"changes_in_resilience.txt")

pl6=ggplot(data = b1, aes(x=efficience, y=interf, fill=diffp,z=diffp))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=13),axis.title.x=element_blank())+
#scale_fill_gradientn(colors=c("darkred","pink","white","dodgerblue"),values=c(0,mid-0.1,mid,1),labels = percent,n.breaks=3)+
scale_fill_gradient2(low ="darkred",high="dodgerblue4",
labels = percent,n.breaks=4)+labs(fill='')+ggtitle("f",subtitle ="Changes in resilience")+
xlab(expression("Mutualism strength "(alpha)))+
ylab("")+coord_cartesian(expand = FALSE)+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(legende))

b=dat %>% dplyr::group_by(interf,efficience,model2) %>% dplyr::summarise(moy=mean(robust,na.rm=T))

pl7=ggplot(data = b, aes(x=efficience, y=interf, fill=moy,z=moy))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=13),plot.background = element_blank())+
scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),
labels = percent,breaks=c(0.5,0.65,0.8,0.95))+
labs(fill='')+ggtitle("g",subtitle ="Network robustness")+
xlab(expression("Mutualism strength "(alpha)))+
ylab("")+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(model2))


b1=subset(dat,model=="model1") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(robust,na.rm=T))
b2=subset(dat,model=="model2") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(robust,na.rm=T))
b1$diffp=(b2$moy-b1$moy)/b1$moy
mid=abs(min(b1$diffp,na.rm=T))/(max(b1$diffp,na.rm=T)-min(b1$diffp,na.rm=T))
b1$legende=""

pl8=ggplot(data = b1, aes(x=efficience, y=interf, fill=diffp,z=diffp))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=13))+
#scale_fill_gradientn(colors=c("darkred","pink","white","dodgerblue"),values=c(0,mid-0.1,mid,1),labels = percent,n.breaks=3)+
scale_fill_gradient2(low ="darkred",high="dodgerblue4",
labels = percent,n.breaks=4)+labs(fill='')+ggtitle("h",subtitle ="Changes in robustness")+
xlab(expression("Mutualism strength "(alpha)))+
ylab("")+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(legende))

plot_grid(pl1,pl2,pl3,pl4,pl5,pl6,pl7,pl8,ncol=2,rel_widths=c(1.6,1),align="hv")


setwd(dir="C:/Users/Duchenne/Documents/EPHI/temporal-dynamics/data and results")
pdf("figure_2.pdf",height=13,width=9.5)
plot_grid(pl1,pl2,pl3,pl4,pl5,pl6,pl7,pl8,ncol=2,rel_widths=c(1.6,1),align="hv")
dev.off();

dim(subset(dat,!is.na(stab)))

#### STABILITY IN THE FAISIBILITY DOMAIN (Fig. 3) 
seuil=0.98

b=dat %>% dplyr::group_by(interf,efficience,model2) %>% dplyr::summarise(moy=log(mean(stab[pers>seuil],na.rm=T)))

pl1=ggplot(data = b, aes(x=efficience, y=interf, fill=moy,z=moy))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=13))+
scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))))+
labs(fill='')+ggtitle("a",subtitle ="Network resilience")+
xlab("")+
ylab(expression("Competition for mutualistic partners  "(c)))+coord_cartesian(expand = FALSE)+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(model2))

b1=subset(dat,model=="model1") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(stab[pers>seuil],na.rm=T))
b2=subset(dat,model=="model2") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(stab[pers>seuil],na.rm=T))
b1$diffp=(b2$moy-b1$moy)/b1$moy
mid=abs(min(b1$diffp,na.rm=T))/(max(b1$diffp,na.rm=T)-min(b1$diffp,na.rm=T))
b1$legende=""

pl2=ggplot(data = b1, aes(x=efficience, y=interf, fill=diffp,z=diffp))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=13))+
#scale_fill_gradientn(colors=c("darkred","pink","white","dodgerblue"),values=c(0,mid-0.1,mid,1),labels = percent,n.breaks=3)+
scale_fill_gradient2(low ="darkred",high="dodgerblue4",
labels = percent,n.breaks=4)+labs(fill='')+ggtitle("b",subtitle ="Changes in resilience")+
xlab("")+
ylab("")+coord_cartesian(expand = FALSE)+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(legende))

b=dat %>% dplyr::group_by(interf,efficience,model2) %>% dplyr::summarise(moy=mean(robust[pers>seuil],na.rm=T))

pl3=ggplot(data = b, aes(x=efficience, y=interf, fill=moy,z=moy))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=13),plot.background = element_blank())+
scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),
labels = percent)+
labs(fill='')+ggtitle("c",subtitle ="Network robustness")+
xlab(expression("Mutualism strength "(alpha)))+
ylab(expression("Competition for mutualistic partners  "(c)))+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(model2))


b1=subset(dat,model=="model1") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(robust[pers>seuil],na.rm=T))
b2=subset(dat,model=="model2") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(robust[pers>seuil],na.rm=T))
b1$diffp=(b2$moy-b1$moy)/b1$moy
mid=abs(min(b1$diffp,na.rm=T))/(max(b1$diffp,na.rm=T)-min(b1$diffp,na.rm=T))
b1$legende=""

pl4=ggplot(data = b1, aes(x=efficience, y=interf, fill=diffp,z=diffp))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=13))+
#scale_fill_gradientn(colors=c("darkred","pink","white","dodgerblue"),values=c(0,mid-0.1,mid,1),labels = percent,n.breaks=3)+
scale_fill_gradient2(low ="darkred",high="dodgerblue4",
labels = percent,n.breaks=4)+labs(fill='')+ggtitle("d",subtitle ="Changes in robustness")+
xlab(expression("Mutualism strength "(alpha)))+
ylab("")+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(legende))


leg1 <- get_legend(pl1)
leg2 <- get_legend(pl3)
# Convert to a ggplot and print
leg1=as_ggplot(leg1)
leg2=as_ggplot(leg2)
pl1=pl1+theme(legend.position="none")
pl3=pl3+theme(legend.position="none")

pdf("figure_3.pdf",height=8,width=10)
plot_grid(pl1,leg1,pl2,pl3,leg2,pl4,ncol=3,rel_widths=c(2.4,0.4,1.75),align="hv")
dev.off();

png("figure_3.png",width=1500,height=1600,res=140)
plot_grid(pl1,leg1,pl2,pl3,leg2,pl4,pl5,leg3,pl6,ncol=3,rel_widths=c(2.2,0.4,1.7),align="hv")
dev.off();


### RESULTS PER SITE:
##### FIGURE S5:
b1=subset(dat,model=="model1") %>% dplyr::group_by(interf,efficience,site,elevation) %>% dplyr::summarise(moy=length(pers[pers>seuil])/(250*11))
b2=subset(dat,model=="model2") %>% dplyr::group_by(interf,efficience,site,elevation) %>% dplyr::summarise(moy=length(pers[pers>seuil])/(250*11))
b1$diffp=(b2$moy-b1$moy)/b1$moy
b1$diffp[is.na(b1$diffp)]=0
mid=abs(min(b1$diffp))/(max(b1$diffp)-min(b1$diffp))
b1$legende=""

png("figureS5.png",width=1200,height=800,res=120)
ggplot(data = b1, aes(x=efficience, y=interf, fill=diffp,z=diffp))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),
plot.subtitle=element_text(size=13),axis.title=element_text(size=14))+
#scale_fill_gradientn(colors=c("darkred","pink","white","dodgerblue"),values=c(0,mid-0.1,mid,1),labels = percent,n.breaks=3)+
scale_fill_gradient2(low ="darkred",high="dodgerblue4",
labels = percent,n.breaks=4)+
labs(fill='')+ggtitle("",subtitle ="Changes in feasibility")+
xlab(expression("Mutualism strength "(alpha)))+
ylab(expression("Competition for mutualistic partners  "(c)))+coord_fixed(ratio=1,expand = FALSE)+facet_wrap(~site)
dev.off();

##### FIGURE S6:
b1=subset(dat,model=="model1") %>% dplyr::group_by(interf,efficience,site,elevation) %>% dplyr::summarise(moy=mean(pers))
b2=subset(dat,model=="model2") %>% dplyr::group_by(interf,efficience,site,elevation) %>% dplyr::summarise(moy=mean(pers))
b1$diffp=(b2$moy-b1$moy)/b1$moy
b1$diffp[is.na(b1$diffp)]=0
mid=abs(min(b1$diffp))/(max(b1$diffp)-min(b1$diffp))
b1$legende=""

ggplot(data=b1,aes(x=elevation,y=diffp))+geom_point()+stat_smooth(method="lm")+facet_wrap(~interf+efficience,scales="free_y")


png("figureS6.png",width=1200,height=800,res=120)
ggplot(data = b1, aes(x=efficience, y=interf, fill=diffp,z=diffp))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),
plot.subtitle=element_text(size=13),axis.title=element_text(size=14))+
#scale_fill_gradientn(colors=c("darkred","pink","white","dodgerblue"),values=c(0,mid-0.1,mid,1),labels = percent,n.breaks=3)+
scale_fill_gradient2(low ="darkred",high="dodgerblue4",
labels = percent,n.breaks=4)+
labs(fill='')+ggtitle("",subtitle ="Changes in persistence")+
xlab(expression("Mutualism strength "(alpha)))+
ylab(expression("Competition for mutualistic partners  "(c)))+coord_fixed(ratio=1,expand = FALSE)+facet_wrap(~site)
dev.off();

##### FIGURE S8:
b1=subset(dat,model=="model1") %>% dplyr::group_by(interf,efficience,site) %>% dplyr::summarise(moy=mean(stab,na.rm=T))
b2=subset(dat,model=="model2") %>% dplyr::group_by(interf,efficience,site) %>% dplyr::summarise(moy=mean(stab,na.rm=T))
b1$diffp=(b2$moy-b1$moy)/b1$moy
mid=abs(min(b1$diffp))/(max(b1$diffp)-min(b1$diffp))
b1$legende=""

png("figureS8.png",width=1200,height=800,res=120)
ggplot(data = b1, aes(x=efficience, y=interf, fill=diffp,z=diffp))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),
plot.subtitle=element_text(size=13),axis.title=element_text(size=14))+
#scale_fill_gradientn(colors=c("darkred","pink","white","dodgerblue"),values=c(0,mid-0.1,mid,1),labels = percent,n.breaks=3)+
scale_fill_gradient2(low ="darkred",high="dodgerblue4",
labels = percent,n.breaks=4)+
labs(fill='')+ggtitle("",subtitle ="Changes in resilience")+
xlab(expression("Mutualism strength "(alpha)))+
ylab(expression("Competition for mutualistic partners  "(c)))+coord_fixed(ratio=1,expand = FALSE)+facet_wrap(~site)
dev.off();

##### FIGURE S9:
b1=subset(dat,model=="model1") %>% dplyr::group_by(interf,efficience,site) %>% dplyr::summarise(moy=mean(robust,na.rm=T))
b2=subset(dat,model=="model2") %>% dplyr::group_by(interf,efficience,site) %>% dplyr::summarise(moy=mean(robust,na.rm=T))
b1$diffp=(b2$moy-b1$moy)/b1$moy
mid=abs(min(b1$diffp))/(max(b1$diffp)-min(b1$diffp))
b1$legende=""

png("figureS9.png",width=1200,height=800,res=120)
ggplot(data = b1, aes(x=efficience, y=interf, fill=diffp,z=diffp))+geom_raster(interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),
plot.subtitle=element_text(size=13),axis.title=element_text(size=14))+
#scale_fill_gradientn(colors=c("darkred","pink","white","dodgerblue"),values=c(0,mid-0.1,mid,1),labels = percent,n.breaks=3)+
scale_fill_gradient2(low ="darkred",high="dodgerblue4",
labels = percent,n.breaks=4)+
labs(fill='')+ggtitle("",subtitle ="Changes in robustness")+
xlab(expression("Mutualism strength "(alpha)))+
ylab(expression("Competition for mutualistic partners  "(c)))+coord_fixed(ratio=1,expand = FALSE)+facet_wrap(~site)
dev.off();

#RELATION SHIPS BETWEEN VARIABLES:
####FIGURE S10:
dat$pers_round=plyr::round_any(dat$pers,0.01)
dat$stab_round=plyr::round_any(log(dat$stab),0.1)
b=dat %>% group_by(stab_round,pers_round,model2) %>% summarise(rob=mean(robust,na.rm=T),n=length(robust))

pl1=ggplot(data=b,aes(y=stab_round,x=pers_round,size=log(n),alpha=log(n)))+geom_point(aes(col=rob),stat="identity",shape=15)+
scale_color_viridis(na.value="white",option = "plasma")+
theme_bw()+theme(panel.grid=element_blank(),strip.background=element_rect(fill="white"),plot.title=element_text(size=14,face="bold",hjust = 0))+
xlab("Network persistence")+ylab("log(resilience)")+labs(col="Robustness")+facet_wrap(~model2)+scale_alpha_continuous(range=c(1,1),guide=F)+
scale_size_continuous(range=c(0.5,2),guide=F)+ggtitle("a")

pl2=ggplot(data=dat,aes(y=compt,x=connectance_I,col=pers))+geom_point(stat="identity",size=0.01)+scale_color_viridis(na.value="white")+
theme_bw()+theme(panel.grid=element_blank(),strip.background=element_rect(fill="white"),plot.title=element_text(size=14,face="bold",hjust = 0))+
xlab("Connectance")+ylab("Interaction overlap")+labs(col="Persistence")+facet_wrap(~model2)+scale_alpha_continuous(range=c(1,1),guide=F)+
scale_size_continuous(range=c(0.5,2),guide=F)+ggtitle("b")

pl3=ggplot(data=dat,aes(y=P+H,x=ntot,col=log(stab)))+geom_point(stat="identity",size=0.01)+scale_color_viridis(na.value="white",option="mako",direction=-1)+
theme_bw()+theme(panel.grid=element_blank(),strip.background=element_rect(fill="white"),plot.title=element_text(size=14,face="bold",hjust = 0))+xlab("Diversity")+
ylab("Total abundance")+labs(col="log(resilience)")+facet_wrap(~model2)+scale_alpha_continuous(range=c(1,1),guide=F)+
scale_size_continuous(range=c(0.5,2),guide=F)+ggtitle("c")

png("figure_S10.png",width=1600,height=2100,res=200)
grid.arrange(pl1,pl2,pl3,ncol=1)
dev.off();


#### FEASIBILITY DOMAIN (Fig. S4)
dat=dat[order(dat$model,dat$interf,dat$efficience,dat$essai,dat$site),]
bidon=subset(dat,model2=="with seasonal structure")
bidon$feas=bidon$feas+(bidon$feas-subset(dat,model2=="without seasonal structure")$feas)
bidon$feas2="without and with seasonsal structure"
bidon$feas2[bidon$feas==0]="non feasible"
bidon$feas2[bidon$feas==(-1)]="only without seasonsal structure"
bidon$feas2[bidon$feas==(2)]="only with seasonsal structure"
bidon$interf=paste0("c=",bidon$interf)
p=ggplot(data=bidon,aes(fill=as.factor(feas2),color=as.factor(feas2),alpha=as.factor(feas2),y=exp(shape_h),x=exp(shape_p)))+
geom_point()+
facet_grid(rows=vars(interf),cols=vars(efficience),
as.table=F,labeller = label_bquote(cols = alpha==.(efficience)))+
ylab(expression(paste(b," parameter of distribution of hummingbird growth rates",sep="")))+
labs(fill="Feasibility",color="Feasibility",alpha="Feasibility")+xlab(expression(paste(b," parameter of distribution of plant growth rates",sep="")))+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+scale_fill_manual(values=c("lightgrey","dodgerblue1","gold3","black"))+
scale_alpha_manual(values=c(0.1,0.4,0.4,0.6))+scale_color_manual(values=c("lightgrey","dodgerblue1","gold3","black"))+
guides(color = guide_legend(override.aes = list(alpha=0.4) ) )+ggtitle("a")+coord_fixed(ratio=1)

df <- data.frame(x=rep(seq(0,1,length=100),3),b=rep(c(0.3,1,15),each=100))
df$y <- mapply(function(x,y){dbeta(x, 1, y)},x=df$x,y=df$b)
p2=ggplot(data=df,aes(x=-x/2,y=y,color=as.factor(b))) + geom_line()+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+ggtitle("b")+ylab("Density")+xlab("Growth rate value")+labs(color="b")+scale_color_viridis_d()

leg1 <- get_legend(p)
leg2 <- get_legend(p2)
# Convert to a ggplot and print
p=p+theme(legend.position="none")
p2=p2+theme(legend.position="none")
png("figure_S4.png",width=1800,height=2500,res=170)
grid.arrange(p,leg1,p2,leg2,ncol=2,heights=c(1.6,1),widths=c(2.7,1))
dev.off();



#### Relationships among variables (Fig. s11)
p1=ggplot(data=dat,aes(y=stab,x=ntot,col=logit(robust,adjust=0.01)))+geom_jitter(stat="identity",size=0.01,width=0.01,alpha=0.4)+scale_color_viridis(na.value="white",option = "plasma")+
theme_bw()+theme(panel.grid=element_blank(),strip.background=element_rect(fill="white"),plot.title=element_text(size=14,face="bold",hjust = 0))+
xlab("Diversity")+ylab("Resilience")+labs(col="Robustness (logit scale)")+facet_wrap(~model2)+
stat_smooth(col="black",fill="lightgrey")+ scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))+
scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))+
annotation_logticks(sides = 'lb')
png("figure_S11.png",width=1800,height=1000,res=170)
p1
dev.off();