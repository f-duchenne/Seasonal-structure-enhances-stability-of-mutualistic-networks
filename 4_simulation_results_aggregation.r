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
intdat_raw <- fread("data_for_modelo.txt")

#charge data
sites=unique(intdat_raw$site)
dat=NULL
dat2f=NULL
for(i in 1:length(sites)){
dat1=fread(paste0("res_net_",sites[i],"_c_efficience_cs_0.txt"))
dat2=fread(paste0("resspecies_",sites[i],"_c_efficience_cs_0.txt"))
dat1$site=sites[i]
dat=rbind(dat,dat1)
dat2$site=sites[i]
dat2f=rbind(dat2f,dat2)
}
#add elevation
elev=intdat_raw %>% group_by(site) %>% summarise(elevation=mean(elev,na.rm=T))
dat=merge(dat,elev,by=c("site"),all.x=T,all.y=F)

#add robustness with one extinctions
rob=NULL
for(i in 1:length(sites)){
rob1=fread(paste0("robustesse_",sites[i],"_c_efficience_cs_0.txt"))
rob=rbind(rob,rob1)
}

robf=subset(rob,!is.na(abond_remov)) %>% group_by(interf,site,efficience,essai,model) %>%
summarise(robust=mean((nbsp_p_per2+nbsp_h_per2)/(nbsp_p_per1+nbsp_h_per1-1)),robust_p=mean(nbsp_p_per2/nbsp_p_per1),
robust_abond=mean((P_eq2+H_eq2)/(P_eq1+H_eq1-abond_remov)))
dat=merge(dat,robf,by=c("essai","model","efficience","interf","site"),all=T)
robf0=subset(rob,is.na(abond_remov)) %>% group_by(interf,site,efficience,essai,model) %>%
summarise(robust0=mean((nbsp_p_per2+nbsp_h_per2)/(nbsp_p_per1+nbsp_h_per1)),robust_p0=mean(nbsp_p_per2/nbsp_p_per1),
robust_abond0=mean((P_eq2+H_eq2)/(P_eq1+H_eq1)))
dat=merge(dat,robf0,by=c("essai","model","efficience","interf","site"),all=T)

setwd(dir="C:/Users/Duchenne/Documents/EPHI/temporal-dynamics/data and results/robustness_2_extinctions")
#add robustness with two extinctions
rob=NULL
for(i in 1:length(sites)){
rob1=fread(paste0("robustesse_",sites[i],"_c_efficience_cs_0.txt"))
rob=rbind(rob,rob1)
}
robf2=subset(rob,!is.na(abond_remov)) %>% group_by(interf,site,efficience,essai,model) %>%
summarise(robust2=mean((nbsp_p_per2+nbsp_h_per2)/(nbsp_p_per1+nbsp_h_per1-2)),robust_p2=mean(nbsp_p_per2/nbsp_p_per1),
robust_abond2=mean((P_eq2+H_eq2)/(P_eq1+H_eq1-abond_remov)))
dat=merge(dat,robf2,by=c("essai","model","efficience","interf","site"),all=T)
setwd(dir="C:/Users/Duchenne/Documents/EPHI/temporal-dynamics/data and results")

#calculate robustness final
dat$robust_moy=apply(dat[,c("robust","robust2")],1,mean,na.rm=T)
dat$robust_p_moy=apply(dat[,c("robust_p","robust_p2")],1,mean,na.rm=T)
dat$robust_abond_moy=apply(dat[,c("robust_abond","robust_abond2")],1,mean,na.rm=T)

#Effective competition metrics:
omegab_phen_li=list()
omegap_phen_li=list()
omegab_sansphen_li=list()
omegap_sansphen_li=list()
I_li=list()
Iphen_li=list()
bidonf=NULL
for(i in 1:length(sites)){
setwd(dir="C:/Users/Duchenne/Documents/EPHI/temporal-dynamics/initial")
load(paste(sites[i],1,".RData",sep="_"))
omegab_phen_li[[i]]=omegab_phen
omegap_phen_li[[i]]=omegap_phen
omegab_sansphen_li[[i]]=omegab_sansphen
omegap_sansphen_li[[i]]=omegap_sansphen
I_li[[i]]=I
Iphen_li[[i]]=I*Phen
bidon=data.frame(compt_ini=c(mean(c(omegab_phen,omegap_phen)),mean(c(omegab_sansphen,omegap_sansphen))),connectance_I_ini=c(mean(I*Phen),mean(I)),model=c("model2","model1"),
site=sites[i])
bidonf=rbind(bidonf,bidon)
}
dat=merge(dat,bidonf,by=c("model","site"),all.x=T,all.y=F)


setwd(dir="C:/Users/Duchenne/Documents/EPHI/temporal-dynamics/data and results")
#### Add some phenological information:
#match pheno plant-bird
pheno1=intdat_raw %>% dplyr::group_by(site,plant,mnumb) %>% dplyr::summarise(nm=length(hummingbird))
pheno1=pheno1 %>% dplyr::group_by(site,plant) %>% dplyr::mutate(nmtot=sum(nm))
pheno1$phen=pheno1$nm/pheno1$nmtot
pheno1$pres=0
pheno1$pres[pheno1$phen>=0]=1
summary(pheno1 %>% dplyr::group_by(site,plant) %>% dplyr::summarise(phenlen=sum(pres)))
pheno2=intdat_raw %>% dplyr::group_by(site,hummingbird,mnumb) %>% dplyr::summarise(nm=length(plant))
pheno2=pheno2 %>% dplyr::group_by(site,hummingbird) %>% dplyr::mutate(nmtot=sum(nm))
pheno2$phen=pheno2$nm/pheno2$nmtot
pheno2$pres=0
pheno2$pres[pheno2$phen>=0]=1
summary(pheno2 %>% dplyr::group_by(site,hummingbird) %>% dplyr::summarise(phenlen=sum(pres)))

########
#METRICS AT NETWORK LEVEL FROM METRICS AT SPECIES LEVEL
bidon=dat2f %>% group_by(type,essai,model,efficience,interf,site) %>% summarise(abond=sum(N_eq1),abond_sd=sd(N_eq1)/mean(N_eq1),rmean=mean(rzero),rvar=var(rzero))
bidon=melt(bidon,id.vars=c("essai","model","efficience","interf","site","type"))
bidon$varia=paste0(bidon$type,gsub("abond","",bidon$variable))
bidon=dcast(bidon,essai+model+efficience+site+interf~varia,value.var=c("value"))
dat=merge(dat,bidon,by=c("essai","model","efficience","interf","site"))
which(duplicated(dat, by=c("essai","model","efficience","interf","site"), fromLast=TRUE))
bidon=dat2f %>% group_by(essai,model,efficience,interf,site) %>% summarise(div=vegan::diversity(N_eq1,index = "shannon"),skew=skewness(N_eq1))
dat=merge(dat,bidon,by=c("essai","model","efficience","interf","site"))
which(duplicated(dat, by=c("essai","model","efficience","interf","site"), fromLast=TRUE))

#TRANSFORM SOME VARIABLES
dat$ntot=round(dat$nbsp_h*dat$nbsp_h_per+dat$nbsp_p*dat$nbsp_p_per)
dat$pers=dat$ntot/(dat$nbsp_p+dat$nbsp_h)
dat$model2="without seasonal structure"
dat$model2[dat$model=="model2"]="with seasonal structure"
dat$stab=-1*(dat$valprop)
dat$feas=0
dat$feas[dat$pers>0.98]=1


fwrite(dat,"results_simulations_network_level.txt")