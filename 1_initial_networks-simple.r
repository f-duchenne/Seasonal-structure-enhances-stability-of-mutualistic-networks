library(data.table)
library(dplyr)
library(doParallel)
library(foreach)
library(parallel)
library(ggplot2)
library(viridis)
setwd(dir="C:/Users/Duchenne/Documents/EPHI/temporal-dynamics/data and results")
intdat_raw <- fread("data_for_modelo.txt") #data accessible on Datadryad

sites=unique(intdat_raw$site)
for(i in 1:length(sites)){
sitechoisi=sites[i]

func=function(x) 
curve(func,0,1)

#annual interaction matrix
bidon=subset(intdat_raw,site==sitechoisi)
b=dcast(bidon,plant~hummingbird)
b2=as.matrix(b[,-1])
b2[b2>0]=1
I=b2
rownames(I)=b$plant
nbsp_h=ncol(I)
nbsp_p=nrow(I)
nbsp_h_dep=nbsp_h
nbsp_p_dep=nbsp_p
Ibin=I


#add all possible combination for each month:
b=dcast(bidon,mnumb+plant~hummingbird)
mer=expand.grid(unique(bidon$mnumb),unique(bidon$plant))
names(mer)=c("mnumb","plant")
b=merge(b,mer,by=c("mnumb","plant"),all=T)
b[is.na(b)]=0

#match pheno plant-bird
pheno1=bidon %>% dplyr::group_by(plant,mnumb) %>% dplyr::summarise(nm=sum(value),nsa=length(unique(y)))
pheno1=pheno1 %>% dplyr::group_by(plant) %>% dplyr::mutate(nmtot=sum(nm/nsa))
pheno1$phen=(pheno1$nm/pheno1$nsa)/pheno1$nmtot
pheno2=bidon %>% dplyr::group_by(hummingbird,mnumb) %>% dplyr::summarise(nm=sum(value),nsa=length(unique(y)))
pheno2=pheno2 %>% dplyr::group_by(hummingbird) %>% dplyr::mutate(nmtot=max(nm/nsa))
pheno2$phen=(pheno2$nm/pheno2$nsa)/pheno2$nmtot

#matrix of plant-hummingbird phenological match
Phen=matrix(0,nbsp_p,nbsp_h)
listep=unique(pheno1$plant)
listeh=unique(pheno2$hummingbird)
for(ii in 1:nbsp_p){
for(iii in 1:nbsp_h){
bid1=subset(pheno1,plant==listep[ii])
bid2=subset(pheno2,hummingbird==listeh[iii])
bid=merge(bid1,bid2,by="mnumb",all=T)
bid[is.na(bid)]=0
Phen[ii,iii]=sum(apply(bid[,c("phen.x","phen.y")],1,prod))
}}


#matrix of hummingbird-hummingbird phenological match
Mbird=matrix(0,nbsp_h,nbsp_h)
for(ii in 1:nbsp_h){
for(iii in ii:nbsp_h){
bid1=subset(pheno2,hummingbird==listeh[ii])
bid2=subset(pheno2,hummingbird==listeh[iii])
bid=merge(bid1,bid2,by="mnumb",all=T)
bid[is.na(bid)]=0
Mbird[ii,iii]=mean(apply(bid[,c("phen.x","phen.y")],1,prod))
}}
Mbird=Mbird+t(Mbird)
diag(Mbird)=1

#matrix of plant-plant phenological match
Mplant=matrix(0,nbsp_p,nbsp_p)
for(ii in 1:nbsp_p){
for(iii in ii:nbsp_p){
bid1=subset(pheno1,plant==listep[ii])
bid2=subset(pheno1,plant==listep[iii])
bid=merge(bid1,bid2,by="mnumb",all=T)
bid[is.na(bid)]=0
Mplant[ii,iii]=sum(apply(bid[,c("phen.x","phen.y")],1,prod))
}}
Mplant=Mplant+t(Mplant)
diag(Mplant)=1

#Calculate phenological overlap among hummingbirds for each given plant
Mbird_p=array(0, dim=c(nbsp_p,nbsp_h,nbsp_h))
for(j in 1:nbsp_p){
for(ii in 1:nbsp_h){
for(iii in ii:nbsp_h){
bid0=subset(pheno1,plant==listep[j])
names(bid0)[6]="phen_p"
bid1=subset(pheno2,hummingbird==listeh[ii])
bid2=subset(pheno2,hummingbird==listeh[iii])
bid=merge(bid1,bid2,by="mnumb",all=T)
bid=merge(bid0,bid,by="mnumb",all=T)
bid[is.na(bid)]=0
Mbird_p[j,ii,iii]=sum(apply(bid[,c("phen_p","phen.x","phen.y")],1,prod))
Mbird_p[j,iii,ii]=Mbird_p[j,ii,iii]
}}
}

#Calculate phenological overlap among plants for each given hummingbird
Mplant_p=array(0, dim=c(nbsp_h,nbsp_p,nbsp_p))
for(j in 1:nbsp_h){
for(ii in 1:nbsp_p){
for(iii in ii:nbsp_p){
bid0=subset(pheno2,hummingbird==listeh[j])
names(bid0)[6]="phen_h"
bid1=subset(pheno1,plant==listep[ii])
bid2=subset(pheno1,plant==listep[iii])
bid=merge(bid1,bid2,by="mnumb",all=T)
bid=merge(bid0,bid,by="mnumb",all=T)
bid[is.na(bid)]=0
Mplant_p[j,ii,iii]=sum(apply(bid[,c("phen_h","phen.x","phen.y")],1,prod))
Mplant_p[j,iii,ii]=Mplant_p[j,ii,iii]
}}
}


#Matrix of effective competition among hummingbirds when all abundances equal 1, without seasonal structure
omegab_sansphen=matrix(0,nbsp_h,nbsp_h)
for(i in 1:nbsp_h){
I2=I[,i]*I
I2=matrix(I2,ncol=nbsp_h,nrow=nbsp_p)
vec=(apply(I2,2,sum))/sum(I[,i])
omegab_sansphen[i,]=vec
}

#Matrix of effective competition among hummingbirds when all abundances equal 1, with seasonal structure
omegab_phen=matrix(0,nbsp_h,nbsp_h)
Iphen=I*Phen
for(i in 1:nbsp_h){
I2=I[,i]*I*Mbird_p[,i,]
I2=matrix(I2,ncol=nbsp_h,nrow=nbsp_p)
vec=(apply(I2,2,sum))/sum(Iphen[,i])
omegab_phen[i,]=vec
}

#Matrix of effective competition among plants when all abundances equal 1, without seasonal structure
omegap_sansphen=matrix(0,nbsp_p,nbsp_p)
for(i in 1:nbsp_p){
I2=t(I)[,i]*t(I)
I2=matrix(I2,ncol=nbsp_p,nrow=nbsp_h)
vec=(apply(I2,2,sum))/sum(I[i,])
omegap_sansphen[i,]=vec
}

#Matrix of effective competition among plants when all abundances equal 1, with seasonal structure
omegap_phen=matrix(0,nbsp_p,nbsp_p)
Iphen=I*Phen
for(i in 1:nbsp_p){
I2=t(I)[,i]*t(I)*Mplant_p[,i,]
I2=matrix(I2,ncol=nbsp_p,nrow=nbsp_h)
vec=apply(I2,2,sum)/sum(Iphen[i,])
omegap_phen[i,]=vec
}

#draw parameter values
resultat=foreach(jj=1:250)%dopar%{
#other parameters:
Nini=rep(1,nbsp_p+nbsp_h)
initial=Nini

efficience=rep(1.5,nbsp_p+nbsp_h)
handling=rep(0.8,nbsp_p+nbsp_h)

#growth rates:
shape_h=runif(1,log(0.3),log(15))
shape_p=runif(1,log(0.3),log(15))
#r=c(runif(nbsp_h,-1,0),runif(nbsp_p,0,0.5))
r=c(-1*rbeta(nbsp_h,1,exp(shape_h))/2,-1*rbeta(nbsp_p,1,exp(shape_p))/2)

#save initial state
setwd(dir="C:/Users/Duchenne/Documents/EPHI/temporal-dynamics/initial")
save(I,Phen,r,efficience,handling,nbsp_h,nbsp_p,nbsp_p_dep,nbsp_h_dep,Nini,shape_h,shape_p,Mbird,Mbird_p,Mplant,Mplant_p,omegab_sansphen,omegab_phen,omegap_sansphen,omegap_phen,file=paste(sitechoisi,jj,".RData",sep="_"),version = 2)

}
}



