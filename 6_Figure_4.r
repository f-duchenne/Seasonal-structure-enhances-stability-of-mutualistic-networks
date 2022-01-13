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
setwd(dir="C:/Users/Duchenne/Documents/EPHI/temporal-dynamics/data and results")
dat <- fread("results_simulations_network_level.txt")
dat$model2=factor(dat$model2,c("without seasonal structure","with seasonal structure"))

###### PATH ANALYSES:

biche=na.omit(dat[,c("nbsp_p_per","nbsp_h_per","interf","efficience","pers","ai_contrib_hp","ai_contrib_hh","ai_contrib_pp","ai_contrib_ph",
"ai_direct_hp","ai_direct_hh","ai_direct_pp","ai_direct_ph","ai_indirect_hp","ai_indirect_hh","ai_indirect_pp","ai_indirect_ph",
"ai_net_hp","ai_net_hh","ai_net_pp","ai_net_ph","mod_I","shape_h","shape_p",
"P","H","ntot","essai","model2","site","stab","mod_r","NODF_r","NODF_I","model","gaintoth","gaintotp","div","skew","nbsp_p","nbsp_h","connectance_I","connectance_r","compt","robust")])
biche$indw=apply(abs(biche[,c("ai_contrib_hh","ai_contrib_pp")]),1,mean)
biche$inda=apply(abs(biche[,c("ai_contrib_hp","ai_contrib_ph")]),1,mean)
biche$n_p=round(biche$nbsp_p_per*biche$nbsp_p)
biche$n_h=round(biche$nbsp_h_per*biche$nbsp_h)
biche$season=0
biche$season[biche$model=="model2"]=1
biche$biomass=biche$P+biche$H
biche$essai2=paste(biche$site,biche$essai)
biche$stab=log(biche$stab)
biche$ntot=log(biche$ntot)
biche$biomass=biche$biomass
biche$connectance_I=logit(biche$connectance_I,adjust=0.01)
biche$compt=logit(biche$compt,adjust=0.01)
biche$robust=logit(biche$robust,adjust=0.01)
biche[,c("biomass","ntot","NODF_I","NODF_r","mod_r","stab","connectance_I","connectance_r","compt","robust")]=
as.data.frame(scale(biche[,c("biomass","ntot","NODF_I","NODF_r","mod_r","stab","connectance_I","connectance_r","compt","robust")]))
biche2=biche

effetf=NULL
effetdivf=NULL
tab=expand.grid(efficience=c(1,2,3),interf=c(3,0.5))
tab$interf[tab$efficience==2 & tab$interf==0.5]=1.5
tab$interf[tab$efficience==3 & tab$interf==0.5]=2
tab$efficience[tab$efficience==3 & tab$interf==2]=1.5
tab$label=letters[1:nrow(tab)]
pdf("path_anal.pdf",width=13,height=11)
par(mfrow=c(3,2))
for(i in 1:nrow(tab)){
biche=subset(biche2,interf==tab$interf[i] & efficience==tab$efficience[i])
modeldiv=lme(ntot~season,data=biche,random=~1|site/essai,control = lmeControl(opt = "optim"))
modelabond=lme(biomass~season+ntot,data=biche,random=~1|site/essai,control = lmeControl(opt = "optim"))
modelnest=lme(connectance_I~season+ntot,data=biche,random=~1|site/essai,control = lmeControl(opt = "optim"))
modelcomp=lme(compt~season+ntot,data=biche,random=~1|site/essai,control = lmeControl(opt = "optim"))
modelstab=lme(stab~(ntot+connectance_I+biomass+compt),data=biche,random=~1|site/essai,control = lmeControl(opt = "optim"))
biche$stab_resi=residuals(modelstab)
modelstab2=lme(stab_resi~season,data=biche,random=~1|site/essai)
modelrobust=lme(robust~(ntot+connectance_I+biomass+compt+stab),data=biche,random=~1|site/essai,control = lmeControl(opt = "optim"))
biche$rob_resi=residuals(modelrobust)
modelrobust2=lme(rob_resi~season,data=biche,random=~1|site/essai)
print(vif(modelstab))
obj=piecewiseSEM::psem(modeldiv,modelnest,modelcomp,modelabond,modelstab,modelrobust)
objb=summary(obj)
left=0
right=30
center=(right-left)/2
haut=20
bas=0
cex=1.5
echelle_fleche=3
l=objb$coefficients[,c("Predictor","Response","Estimate","P.Value")]
l=l %>% separate(Predictor, c("Season", "Predictor"),sep=":")
vec=l$Season[is.na(l$Predictor)]
l$Season[is.na(l$Predictor)]=NA
l$Predictor[is.na(l$Season)]=vec
l=rbind(l,data.frame(Season=NA,Predictor="season",Response="stab",Estimate=modelstab2$coef$fixed[2],P.Value=Anova(modelstab2)[1,3]),
data.frame(Season=NA,Predictor="season",Response="robust",Estimate=modelrobust2$coef$fixed[2],P.Value=Anova(modelrobust2)[1,3]))
l=l[,c("Predictor","Response","Estimate","P.Value","Season")]
l$Estimate[!is.na(l$Season) & l$Response=="biomass"]=l$Estimate[!is.na(l$Season) & l$Response=="biomass"]+l$Estimate[(l$Predictor %in% l$Predictor[!is.na(l$Season)]) & l$Response=="biomass" & is.na(l$Season)]
l$Estimate[!is.na(l$Season) & l$Response=="connectance_I"]=l$Estimate[!is.na(l$Season) & l$Response=="connectance_I"]+l$Estimate[(l$Predictor %in% l$Predictor[!is.na(l$Season)]) & l$Response=="connectance_I" & is.na(l$Season)]
l$Estimate[!is.na(l$Season) & l$Response=="stab"]=l$Estimate[!is.na(l$Season) & l$Response=="stab"]+l$Estimate[(l$Predictor %in% l$Predictor[!is.na(l$Season)]) & l$Response=="stab" & is.na(l$Season)]
l$Estimate[!is.na(l$Season) & l$Response=="compt"]=l$Estimate[!is.na(l$Season) & l$Response=="compt"]+l$Estimate[(l$Predictor %in% l$Predictor[!is.na(l$Season)]) & l$Response=="compt" & is.na(l$Season)]
l[,3]=as.numeric(l[,3])
l[,4]=as.numeric(l[,4])
l[3:4,]=l[4:3,]
l$colo=ifelse(l$Estimate<0,"firebrick4","dodgerblue4")
l$lty=1
l$lty=ifelse((l$Predictor=="season" & l$Response=="stab"),2,l$lty)
l$lty=ifelse((l$Predictor=="season" & l$Response=="robust"),2,l$lty)
l$position=0.5
l$position=ifelse((l$Predictor=="biomass" & l$Response=="robust"),0.7,l$position)
l$position=ifelse((l$Predictor=="connectance_I" & l$Response=="stab"),0.7,l$position)
l$position=ifelse((l$Predictor=="ntot" & l$Response=="connectance_I"),0.65,l$position)
l$position=ifelse((l$Predictor=="ntot" & l$Response=="biomass"),0.4,l$position)
# l$colo=ifelse((l$Response %in% c("ntot","biomass","connectance_I","stab") & l$Predictor=="season"),"grey",l$colo)
l$curv=NA
l$curv[l$Predictor=="season" & l$Response=="robust"]=4.9
l$curv[l$Predictor=="season" & l$Response=="stab"]=-4.9
l$curv[l$Predictor=="ntot" & l$Response=="connectance_I"]=1.8
g <- graph.data.frame(l, directed=T)
g= g %>% set_edge_attr("color", value =l$colo)
coord=data.frame(label=vertex_attr(g, "name"),
x=c(center,center-4,right,left,center+7,left,right),y=c(haut,bas+11,bas+11,bas+11,bas+11,bas,bas),vsize=23)
coord$vsize[coord$label=="biomass"]=32
g= g %>% set_vertex_attr("name", value =c("Seasonal\nstructure",paste0("Diversity\n(",objb$R2$Marginal[objb$R2$Response=="ntot"],"|",objb$R2$Conditional[objb$R2$Response=="ntot"],")"),
paste0("Connectance\n(",objb$R2$Marginal[objb$R2$Response=="connectance_I"],"|",objb$R2$Conditional[objb$R2$Response=="connectance_I"],")"),
paste0("Total abundance\n(",objb$R2$Marginal[objb$R2$Response=="biomass"],"|",objb$R2$Conditional[objb$R2$Response=="biomass"],")"),
paste0("Inter. overlap\n(",objb$R2$Marginal[objb$R2$Response=="compt"],"|",objb$R2$Conditional[objb$R2$Response=="compt"],")"),
paste0("Resilience\n(",objb$R2$Marginal[objb$R2$Response=="stab"],"|",objb$R2$Conditional[objb$R2$Response=="stab"],")"),
paste0("Robustness\n(",objb$R2$Marginal[objb$R2$Response=="robust"],"|",objb$R2$Conditional[objb$R2$Response=="robust"],")")))

EL=as_edgelist(g)
EL=cbind(EL,l[,3]*2)
asi=abs(l[,3])*echelle_fleche
asi[asi<5]=5
asi[asi>=15]=15

qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
border.color="white",label.cex=cex,label.scale=F,lty=l$lty,
edge.label.cex = cex*1.5,edge.label.position=l$position,vsize2=9,vsize=coord$vsize,title=paste0(letters[i]),curve=l$curv,
title.cex=cex,shape="rectangle",edge.labels=T,fade=F,asize=asi,
mar=c(5,5,5,5),knot.border.color="white",curveShape=-0.4)


effet=data.frame(Diversity=c(l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*l$Estimate[(l$Predictor=="ntot" & l$Response=="stab")],
l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*l$Estimate[(l$Predictor=="ntot" & l$Response=="robust")]),
Diversity.indirect=c(l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*(
l$Estimate[(l$Predictor=="ntot" & l$Response=="biomass")]*l$Estimate[(l$Predictor=="biomass" & l$Response=="stab")]+
l$Estimate[(l$Predictor=="ntot" & l$Response=="compt")]*l$Estimate[(l$Predictor=="compt" & l$Response=="stab")]+
l$Estimate[(l$Predictor=="ntot" & l$Response=="connectance_I")]*l$Estimate[(l$Predictor=="connectance_I" & l$Response=="stab")]),
l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*(
l$Estimate[(l$Predictor=="ntot" & l$Response=="biomass")]*l$Estimate[(l$Predictor=="biomass" & l$Response=="robust")]+
l$Estimate[(l$Predictor=="ntot" & l$Response=="compt")]*l$Estimate[(l$Predictor=="compt" & l$Response=="robust")]+
l$Estimate[(l$Predictor=="ntot" & l$Response=="connectance_I")]*l$Estimate[(l$Predictor=="connectance_I" & l$Response=="robust")])),
Abundance=c(l$Estimate[(l$Predictor=="season" & l$Response=="biomass")]*l$Estimate[(l$Predictor=="biomass" & l$Response=="stab")],
l$Estimate[(l$Predictor=="season" & l$Response=="biomass")]*l$Estimate[(l$Predictor=="biomass" & l$Response=="robust")]),
Inter.overlap=c(l$Estimate[(l$Predictor=="season" & l$Response=="compt")]*l$Estimate[(l$Predictor=="compt" & l$Response=="stab")],
l$Estimate[(l$Predictor=="season" & l$Response=="compt")]*l$Estimate[(l$Predictor=="compt" & l$Response=="robust")]),
Connectance=c(l$Estimate[(l$Predictor=="season" & l$Response=="connectance_I")]*l$Estimate[(l$Predictor=="connectance_I" & l$Response=="stab")],
l$Estimate[(l$Predictor=="season" & l$Response=="connectance_I")]*l$Estimate[(l$Predictor=="connectance_I" & l$Response=="robust")]),
varia=c("on resilience","on robustness"),alpha=tab$efficience[i],interf=tab$interf[i])
effetdiv=data.frame(div.abund=l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*c(
l$Estimate[(l$Predictor=="ntot" & l$Response=="biomass")]*l$Estimate[(l$Predictor=="biomass" & l$Response=="stab")],
l$Estimate[(l$Predictor=="ntot" & l$Response=="biomass")]*l$Estimate[(l$Predictor=="biomass" & l$Response=="robust")]),
div.compt=l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*
c(l$Estimate[(l$Predictor=="ntot" & l$Response=="compt")]*l$Estimate[(l$Predictor=="compt" & l$Response=="stab")],
l$Estimate[(l$Predictor=="ntot" & l$Response=="compt")]*l$Estimate[(l$Predictor=="compt" & l$Response=="robust")]),
div.conn=l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*
c(l$Estimate[(l$Predictor=="ntot" & l$Response=="connectance_I")]*l$Estimate[(l$Predictor=="connectance_I" & l$Response=="stab")],
l$Estimate[(l$Predictor=="ntot" & l$Response=="connectance_I")]*l$Estimate[(l$Predictor=="connectance_I" & l$Response=="robust")]),
varia=c("on resilience","on robustness"),alpha=tab$efficience[i],interf=tab$interf[i])
effetf=rbind(effetf,effet)
effetdivf=rbind(effetdivf,effetdiv)
}
dev.off();

objf=NULL
for(i in 1:nrow(tab)){
obj=cbind(tab[i,1:2]-0.25,tab[i,1:2]+0.25)
names(obj)=c("xmin","ymin","xmax","ymax")
obj$group=letters[i]
objf=rbind(objf,obj)
}


b1=subset(dat,model=="model1") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(pers))
b2=subset(dat,model=="model2") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(pers))
b1$moy[is.na(b1$moy)]=0
b2$moy[is.na(b2$moy)]=0
b1$diffp=(b2$moy-b1$moy)/b1$moy
b1$diffp[is.na(b1$diffp)]=0
mid=abs(min(b1$diffp))/(max(b1$diffp)-min(b1$diffp))
b1$legende=""

b3=merge(b1,b2,by=c("interf","efficience"))
b3=subset(b3,moy.x>0.2 & moy.y>0.2 & moy.x<0.9 & moy.y<0.9,select=c(efficience,interf))
objf2=NULL
for(i in 1:nrow(b3)){
obj=cbind(b3[i,1:2]-0.25,b3[i,1:2]+0.25)
names(obj)=c("xmin","ymin","xmax","ymax")
objf2=rbind(objf2,obj)
}

pl2=ggplot()+geom_raster(data = b1, aes(x=efficience, y=interf, fill=diffp,z=diffp),interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),plot.subtitle=element_blank(),strip.text=element_text(size=12),
legend.title=element_text(size=13),axis.title=element_text(size=16))+
#scale_fill_gradientn(colors=c("darkred","pink","white","dodgerblue"),values=c(0,mid-0.1,mid,1),labels = percent,n.breaks=3)+
scale_fill_gradient2(low ="darkred",high="dodgerblue4",labels = percent,n.breaks=4)+
labs(fill='Changes in\npersistence')+ggtitle("g")+
xlab(expression(alpha))+
ylab(expression(italic(c)))+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(legende))+
geom_rect(data=objf2,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),color="red",fill=NA,size=2)+
geom_rect(data=objf,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,group=group),color="black",fill=NA,size=2)+
geom_text(data=tab,aes(x=efficience, y=interf,label=paste0(label)),col="black")


b=subset(dat,model=="model2") %>% dplyr::group_by(interf,efficience,model2) %>% dplyr::summarise(moy=mean(pers))

pl1=ggplot()+geom_raster(data = b, aes(x=efficience, y=interf, fill=moy,z=moy),interpolate=F)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),
plot.subtitle=element_text(size=13),axis.title=element_text(size=14))+
scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),
labels = percent,n.breaks=4,limits=c(0,1))+labs(fill='Persistence')+ggtitle("a",subtitle ="Network persistence")+
xlab(expression(alpha))+
ylab(expression(italic(c)))+coord_fixed(ratio=1,expand = FALSE)+facet_grid(col=vars(model2))+
geom_rect(data=objf2,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),color="red",fill=NA,size=2)+
geom_rect(data=objf,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,group=group),color="black",fill=NA,size=2)+
geom_text(data=tab,aes(x=efficience, y=interf,label=paste0(label)),col="black")

pdf("fig4g.pdf",width=5,height=4)
pl1
dev.off();




################# FIGURE S12


effetf=NULL
effetdivf=NULL
divf=NULL
b1=subset(dat,model=="model1") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(pers))
b2=subset(dat,model=="model2") %>% dplyr::group_by(interf,efficience) %>% dplyr::summarise(moy=mean(pers))
b1$moy[is.na(b1$moy)]=0
b2$moy[is.na(b2$moy)]=0
b1$diffp=(b2$moy-b1$moy)/b1$moy
b1$diffp[is.na(b1$diffp)]=0
mid=abs(min(b1$diffp))/(max(b1$diffp)-min(b1$diffp))
b1$legende=""
b3=merge(b1,b2,by=c("interf","efficience"))
b3=subset(b3,moy.x>0.2 & moy.y>0.2 & moy.x<0.9 & moy.y<0.9,select=c(efficience,interf))
for(i in 1:nrow(b3)){
biche=subset(biche2,interf==b3$interf[i] & efficience==b3$efficience[i])
modeldiv=lme(ntot~season,data=biche,random=~1|site/essai,control = lmeControl(opt = "optim"))
modelabond=lme(biomass~season+ntot,data=biche,random=~1|site/essai,control = lmeControl(opt = "optim"))
modelnest=lme(connectance_I~season+ntot,data=biche,random=~1|site/essai,control = lmeControl(opt = "optim"))
modelcomp=lme(compt~season+ntot,data=biche,random=~1|site/essai,control = lmeControl(opt = "optim"))
modelstab=lme(stab~(ntot+connectance_I+biomass+compt),data=biche,random=~1|site/essai,control = lmeControl(opt = "optim"))
modelrobust=lme(robust~(ntot+connectance_I+biomass+compt+stab),data=biche,random=~1|site/essai,control = lmeControl(opt = "optim"))
print(vif(modelstab))
obj=piecewiseSEM::psem(modeldiv,modelnest,modelcomp,modelabond,modelstab,modelrobust)
objb=summary(obj)

l=objb$coefficients[,c("Predictor","Response","Estimate","P.Value")]
l=l %>% separate(Predictor, c("Season", "Predictor"),sep=":")
vec=l$Season[is.na(l$Predictor)]
l$Season[is.na(l$Predictor)]=NA
l$Predictor[is.na(l$Season)]=vec
l=l[,c("Predictor","Response","Estimate","P.Value","Season")]
l$Estimate[!is.na(l$Season) & l$Response=="biomass"]=l$Estimate[!is.na(l$Season) & l$Response=="biomass"]+l$Estimate[(l$Predictor %in% l$Predictor[!is.na(l$Season)]) & l$Response=="biomass" & is.na(l$Season)]
l$Estimate[!is.na(l$Season) & l$Response=="connectance_I"]=l$Estimate[!is.na(l$Season) & l$Response=="connectance_I"]+l$Estimate[(l$Predictor %in% l$Predictor[!is.na(l$Season)]) & l$Response=="connectance_I" & is.na(l$Season)]
l$Estimate[!is.na(l$Season) & l$Response=="stab"]=l$Estimate[!is.na(l$Season) & l$Response=="stab"]+l$Estimate[(l$Predictor %in% l$Predictor[!is.na(l$Season)]) & l$Response=="stab" & is.na(l$Season)]
l$Estimate[!is.na(l$Season) & l$Response=="compt"]=l$Estimate[!is.na(l$Season) & l$Response=="compt"]+l$Estimate[(l$Predictor %in% l$Predictor[!is.na(l$Season)]) & l$Response=="compt" & is.na(l$Season)]
l[,3]=as.numeric(l[,3])
l[,4]=as.numeric(l[,4])
l[3:4,]=l[4:3,]


effet=data.frame(Diversity=c(l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*l$Estimate[(l$Predictor=="ntot" & l$Response=="stab")],
l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*l$Estimate[(l$Predictor=="ntot" & l$Response=="robust")]),
Diversity.indirect=c(l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*(
l$Estimate[(l$Predictor=="ntot" & l$Response=="biomass")]*l$Estimate[(l$Predictor=="biomass" & l$Response=="stab")]+
l$Estimate[(l$Predictor=="ntot" & l$Response=="compt")]*l$Estimate[(l$Predictor=="compt" & l$Response=="stab")]+
l$Estimate[(l$Predictor=="ntot" & l$Response=="connectance_I")]*l$Estimate[(l$Predictor=="connectance_I" & l$Response=="stab")]),
l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*(
l$Estimate[(l$Predictor=="ntot" & l$Response=="biomass")]*l$Estimate[(l$Predictor=="biomass" & l$Response=="robust")]+
l$Estimate[(l$Predictor=="ntot" & l$Response=="compt")]*l$Estimate[(l$Predictor=="compt" & l$Response=="robust")]+
l$Estimate[(l$Predictor=="ntot" & l$Response=="connectance_I")]*l$Estimate[(l$Predictor=="connectance_I" & l$Response=="robust")])),
Abundance=c(l$Estimate[(l$Predictor=="season" & l$Response=="biomass")]*l$Estimate[(l$Predictor=="biomass" & l$Response=="stab")],
l$Estimate[(l$Predictor=="season" & l$Response=="biomass")]*l$Estimate[(l$Predictor=="biomass" & l$Response=="robust")]),
Inter.overlap=c(l$Estimate[(l$Predictor=="season" & l$Response=="compt")]*l$Estimate[(l$Predictor=="compt" & l$Response=="stab")],
l$Estimate[(l$Predictor=="season" & l$Response=="compt")]*l$Estimate[(l$Predictor=="compt" & l$Response=="robust")]),
Connectance=c(l$Estimate[(l$Predictor=="season" & l$Response=="connectance_I")]*l$Estimate[(l$Predictor=="connectance_I" & l$Response=="stab")],
l$Estimate[(l$Predictor=="season" & l$Response=="connectance_I")]*l$Estimate[(l$Predictor=="connectance_I" & l$Response=="robust")]),
varia=c("on resilience","on robustness"),alpha=b3$efficience[i],interf=b3$interf[i])
effetr=data.frame(Resilience=c(NA,apply(subset(effet,varia=="on resilience")[,1:5],1,sum)*l$Estimate[(l$Predictor=="stab" & l$Response=="robust")]))
effet=cbind(effet,effetr)


effetdiv=data.frame(div.abund=l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*c(
l$Estimate[(l$Predictor=="ntot" & l$Response=="biomass")]*l$Estimate[(l$Predictor=="biomass" & l$Response=="stab")],
l$Estimate[(l$Predictor=="ntot" & l$Response=="biomass")]*l$Estimate[(l$Predictor=="biomass" & l$Response=="robust")]),
div.compt=l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*
c(l$Estimate[(l$Predictor=="ntot" & l$Response=="compt")]*l$Estimate[(l$Predictor=="compt" & l$Response=="stab")],
l$Estimate[(l$Predictor=="ntot" & l$Response=="compt")]*l$Estimate[(l$Predictor=="compt" & l$Response=="robust")]),
div.conn=l$Estimate[(l$Predictor=="season" & l$Response=="ntot")]*
c(l$Estimate[(l$Predictor=="ntot" & l$Response=="connectance_I")]*l$Estimate[(l$Predictor=="connectance_I" & l$Response=="stab")],
l$Estimate[(l$Predictor=="ntot" & l$Response=="connectance_I")]*l$Estimate[(l$Predictor=="connectance_I" & l$Response=="robust")]),
varia=c("on resilience","on robustness"),alpha=b3$efficience[i],interf=b3$interf[i])

div=data.frame(eff=l$Estimate[(l$Predictor=="ntot" & l$Response=="stab")],alpha=b3$efficience[i],interf=b3$interf[i])

effetf=rbind(effetf,effet)
effetdivf=rbind(effetdivf,effetdiv)
divf=rbind(divf,div)
}

divf$ratio=(divf$interf/divf$alpha)/divf$alpha
ggplot(data=divf,aes(x=ratio,color,y=eff))+geom_point(size=3)+scale_color_viridis()+geom_hline(yintercept=0)

names(effetf)[1:4]=c("Diversity (direct)              ","Diversity (indirect)","Total abundance","Inter. overlap")
b=melt(effetf,id.vars=c("alpha","interf","varia"))
b2=b %>% group_by(alpha,interf,varia) %>% summarise(value=sum(value))
bidon=expand.grid(list(alpha=seq(1,3,0.5),interf=seq(0.5,3,0.5),varia=c("on resilience","on robustness"),variable=unique(b$variable)))
b=merge(b,bidon,by=c("alpha","interf","varia","variable"),all.y=T)
b$interf=paste0("c=",b$interf)
b2$interf=paste0("c=",b2$interf)
pl1=ggplot()+
geom_bar(data=subset(b,varia=="on resilience"),aes(fill=variable,y=value,x=as.factor(1)),stat="identity",position="stack")+facet_grid(rows=vars(interf),cols=vars(alpha),
as.table=F,labeller = label_bquote(cols = alpha==.(alpha)))+
geom_point(data=subset(b2,varia=="on resilience"),aes(y=value,x=as.factor(1)),shape = 21,fill="white",size=2)+scale_fill_manual(values=c("dodgerblue4","dodgerblue2",rgb(233,216,166,maxColorValue=255),
"chartreuse4",rgb(202,103,2,maxColorValue=255),"gold2"))+
ylab("Standardized effect of the seasonal structure")+labs(fill="Mediated by")+xlab("")+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
geom_hline(yintercept=0,linetype="dashed")+ggtitle("a",subtitle="On Resilience")

pl2=ggplot()+
geom_bar(data=subset(b,varia=="on robustness"),aes(fill=variable,y=value,x=as.factor(1)),stat="identity",position="stack")+facet_grid(rows=vars(interf),cols=vars(alpha),
as.table=F,labeller = label_bquote(cols = alpha==.(alpha)))+
geom_point(data=subset(b2,varia=="on robustness"),aes(y=value,x=as.factor(1)),shape = 21,fill="white",size =2)+scale_fill_manual(values=c("dodgerblue4","dodgerblue2",rgb(233,216,166,maxColorValue=255),
"chartreuse4",rgb(202,103,2,maxColorValue=255),"gold2"))+
ylab("Standardized effect of the seasonal structure")+labs(fill="Mediated by")+xlab("")+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
geom_hline(yintercept=0,linetype="dashed")+ggtitle("b",subtitle="On Robustness")

leg1 <- get_legend(pl2)
# Convert to a ggplot and print
pl1=pl1+theme(legend.position="none")
pl2=pl2+theme(legend.position="none")

grid.arrange(pl1,pl2,leg1,ncol=3,widths=c(1,1,0.3))


names(effetdivf)[1:3]=c("Diversity - Total abund.","Diversity - Inter. overlap","Diversity - Connectance")
b=melt(effetf,id.vars=c("alpha","interf","varia"))
b3=rbind(melt(effetdivf,id.vars=c("alpha","interf","varia")),subset(b,variable=="Diversity (direct)              "))
bidon=expand.grid(list(alpha=seq(1,3,0.5),interf=seq(0.5,3,0.5),varia=c("on resilience","on robustness"),variable=unique(b3$variable)))
b3=merge(b3,bidon,by=c("alpha","interf","varia","variable"),all.y=T)
b3$interf=paste0("c=",b3$interf)
pl3=ggplot()+
geom_bar(data=subset(b3,varia=="on resilience"),aes(fill=variable,y=value,x=as.factor(1)),stat="identity",position="stack")+
facet_grid(rows=vars(interf),cols=vars(alpha),
as.table=F,labeller = label_bquote(cols = alpha==.(alpha)))+
scale_fill_manual(values=c("lightblue",rgb(10,147,150,maxColorValue=255),rgb(148,210,189,maxColorValue=255),"dodgerblue4"))+
ylab("Standardized effect of the seasonal structure")+labs(fill="Mediated by")+xlab("")+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
geom_hline(yintercept=0,linetype="dashed")+ggtitle("c",subtitle="On Resilience")+scale_y_continuous(breaks = c(-0.1,0,0.1))

pl4=ggplot()+
geom_bar(data=subset(b3,varia=="on robustness"),aes(fill=variable,y=value,x=as.factor(1)),stat="identity",position="stack")+
facet_grid(rows=vars(interf),cols=vars(alpha),
as.table=F,labeller = label_bquote(cols = alpha==.(alpha)))+
scale_fill_manual(values=c("lightblue",rgb(10,147,150,maxColorValue=255),rgb(148,210,189,maxColorValue=255),"dodgerblue4"))+
ylab("Standardized effect of the seasonal structure")+labs(fill="Mediated by")+xlab("")+theme_bw()+
theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
geom_hline(yintercept=0,linetype="dashed")+ggtitle("d",subtitle="On Robustness")+scale_y_continuous(breaks = c(-0.2,0,0.2,0.4))



leg3 <- get_legend(pl3)
# Convert to a ggplot and print
pl3=pl3+theme(legend.position="none")
pl4=pl4+theme(legend.position="none")

grid.arrange(pl1,pl2,leg1,pl3,pl4,leg3,ncol=3,widths=c(1,1,0.8))


png("figS12.png",width=1800,height=1800,res=180)
grid.arrange(pl1,pl2,leg1,pl3,pl4,leg3,ncol=3,widths=c(1,1,0.8))
dev.off();

pdf("figS12.pdf",width=8,height=8)
grid.arrange(pl1,pl2,leg1,pl3,pl4,leg3,ncol=3,widths=c(1,1,0.8))
dev.off();