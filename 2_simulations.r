#setwd(dir="C:/Users/Duchenne/Documents/EPHI/temporal-dynamics/initial")
.libPaths(c("/home/duchenne/R/x86_64-pc-linux-gnu-library/3.6/",.libPaths()))
library(odeintr)
library(data.table)
library(dplyr)
library(doParallel)
library(foreach)
library(parallel)
library(doMPI)
library(igraph)
cl<-startMPIcluster()
registerDoMPI(cl)
#cl<-makeCluster(96)
#registerDoParallel(cl)

conv=1e-14 # convergence criterion
seuil=1e-5 #extinction threshold
fini=250 #number of replicates
comp_vec=seq(0,3,0.5) #competition strength vector
alpha_vec=seq(0,3,0.5) #mutualism strength vector
cs=0 #interspecific competition for space

intdat_raw <- fread("data_for_modelo.txt") #data accessible from the gitlab of Rafi

param="_c_efficience"

#loop over sites
for(site in unique(intdat_raw$site)){
sitechoisi=site

comb <- function(x, ...) {mapply(rbind,x,...,SIMPLIFY=FALSE)}

#prallelized loop over replicates
resultat=foreach(jj=1:fini,.combine=comb)%dopar%{
library(deSolve)
library(rootSolve)
library(bipartite)
set.seed(jj)

final=NULL
resmodel1=NULL
resmodel2=NULL
for(efficience2 in alpha_vec){ #loop over mutualism strength values
### MODEL 1: without seasonal structure
#fonctions
derivs <-function(t, y,parms){
dy=rep(0,1)
for(i in 1:length(Nini)){
if(i<=nbsp_h){
I2=I[,i]*I*y[(1+nbsp_h):(nbsp_h+nbsp_p)]
I2=matrix(I2,ncol=nbsp_h,nrow=nbsp_p)
vec=((apply(I2,2,sum))/sum(y[(1+nbsp_h):(nbsp_h+nbsp_p)]*I[,i]))
vec[which(is.na(vec))]=0
vec[which(is.infinite(vec))]=0
eq1=(r[i]+
efficience[i]*sum(I[,i]*y[(1+nbsp_h):(nbsp_h+nbsp_p)])/(1+handling[i]*sum(I[,i]*y[(1+nbsp_h):(nbsp_h+nbsp_p)])+
interfp*sum(vec*y[1:nbsp_h]))-
sum(CSh[,i]*y[1:nbsp_h]))*y[i]

myenv=new.env(i)
environment(eq1)=myenv
myenv$i=i
myenv$Nini=Nini
myenv$I=I
myenv$I2=I2
myenv$vec=vec
myenv$r=r
myenv$nbsp_h=nbsp_h
myenv$nbsp_p=nbsp_p
myenv$cs=cs
myenv$interfp=interfp
}else{
I2=t(I)[,(i-nbsp_h)]*t(I)*y[1:nbsp_h]
I2=matrix(I2,ncol=nbsp_p,nrow=nbsp_h)
vec=(apply(I2,2,sum))/sum(y[1:nbsp_h]*I[(i-nbsp_h),])
vec[which(is.na(vec))]=0
vec[which(is.infinite(vec))]=0
eq1=(r[i]+
efficience[i]*sum(I[(i-nbsp_h),]*y[1:nbsp_h])/(1+handling[i]*sum(I[(i-nbsp_h),]*y[1:nbsp_h])+
interff*sum(vec*y[(1+nbsp_h):(nbsp_h+nbsp_p)]))-
sum(CSp[,(i-nbsp_h)]*y[(1+nbsp_h):(nbsp_h+nbsp_p)]))*
y[i]

myenv=new.env(i)
environment(eq1)=myenv
myenv$i=i
myenv$I=I
myenv$Nini=Nini
myenv$I2=I2
myenv$vec=vec
myenv$r=r
myenv$nbsp_h=nbsp_h
myenv$nbsp_p=nbsp_p
myenv$cs=cs
myenv$interff=interff
}
if(y[i]>=seuil){dy[i]=eq1}else{dy[i]=0}}
return(list(dy))}

#run it for each value of competition:
for(j in comp_vec){
interff=j
interfp=j
load(paste(sitechoisi,jj,".RData",sep="_"))
efficience=rep(efficience2,nbsp_h+nbsp_p)
initial=Nini
r[r>(-1e-2) & r<0]=-1e-2
r[r<(1e-2) & r>0]=1e-2
CSp=matrix(cs,ncol=nbsp_p,nrow=nbsp_p)
CSh=matrix(cs,ncol=nbsp_h,nrow=nbsp_h)
diag(CSp)=1
diag(CSh)=1


#Solve system until convergence
a=1
b=0
while(a>conv){
if(a==1){out <- ode(y=Nini, times=seq(0,20,0.05), func = derivs,parms = NULL,method="lsoda")}else{
out <- ode(y=Nini, times=seq(0,20,1), func = derivs,parms = NULL,method="lsoda")}
colnames(out)[1]="Time"
Nini=t(out)[2:(1+nbsp_h+nbsp_p),nrow(out)]
Nini[which(Nini<seuil)]=0
if(is.na(mean(Nini))){a=-2}else{
if(b==0){dyna=out}else{
out[,"Time"]=out[,"Time"]+b
dyna=rbind(dyna,out[-1,])}
if(nrow(dyna)>10){a=max(apply(t(dyna[(nrow(dyna)-10):nrow(dyna),2:(1+nbsp_h+nbsp_p)]),1,var))}
if(max(Nini,na.rm=T)>1e5){a=-1}
b=b+20
}
}


#Store basic data at species level
Nini[is.na(Nini)]=0
popf_p=Nini[(nbsp_h+1):(nbsp_h+nbsp_p)]
popf_h=Nini[1:nbsp_h]
nbsp_p_per=length(popf_p[popf_p>seuil])
nbsp_h_per=length(popf_h[popf_h>seuil])
final1=data.frame(species=c(colnames(I),rownames(I)),type=c(rep("H",nbsp_h),rep("P",nbsp_p)),Nini=initial,rzero=r,cs=cs,efficience=efficience,handling=handling,interactions=c(apply(I,2,sum),apply(I,1,sum)),
pheno_inter=NA,pheno_intra=NA,comp=c(apply(omegab_sansphen,1,sum),apply(omegap_sansphen,1,sum)),interf=c(rep(interfp,nbsp_h),rep(interff,nbsp_p)),essai=jj,N_eq1=c(popf_h,popf_p),model="model1",
interactions_eq=NA,pheno_inter_eq=NA,pheno_intra_eq=NA,func_res=NA)

#Estimate saturation of functional response for each species:
for(i in 1:length(Nini)){
if(i<=nbsp_h){
I2=I[,i]*I*Nini[(1+nbsp_h):(nbsp_h+nbsp_p)]
I2=matrix(I2,ncol=nbsp_h,nrow=nbsp_p)
vec=apply(I2,2,sum)/sum(Nini[(1+nbsp_h):(nbsp_h+nbsp_p)]*I[,i])
vec[which(is.na(vec))]=0
vec[which(is.infinite(vec))]=0
final1$func_res[i]=efficience[i]*sum(I[,i]*Nini[(1+nbsp_h):(nbsp_h+nbsp_p)])/(1+handling[i]*sum(I[,i]*Nini[(1+nbsp_h):(nbsp_h+nbsp_p)])+
interfp*sum(vec*Nini[1:nbsp_h]))
}else{
I2=t(I)[,(i-nbsp_h)]*t(I)*Nini[1:nbsp_h]
I2=matrix(I2,ncol=nbsp_p,nrow=nbsp_h)
vec=apply(I2,2,sum)/sum(Nini[1:nbsp_h]*I[(i-nbsp_h),])
vec[which(is.na(vec))]=0
vec[which(is.infinite(vec))]=0
final1$func_res[i]=efficience[i]*sum(I[(i-nbsp_h),]*Nini[1:nbsp_h])/(1+handling[i]*sum(I[(i-nbsp_h),]*Nini[1:nbsp_h])+
interff*sum(vec*Nini[(1+nbsp_h):(nbsp_h+nbsp_p)]))
}}

#If network is not empty then store information at network level
if(nbsp_h_per>0 & nbsp_p_per>0){
### AGREGER:
I=I[which(popf_p>=seuil),which(popf_h>=seuil)]
CSp=CSp[which(popf_p>=seuil),which(popf_p>=seuil)]
CSh=CSh[which(popf_h>=seuil),which(popf_h>=seuil)]
Mbird=Mbird[which(popf_h>=seuil),which(popf_h>=seuil)]
Mplant=Mplant[which(popf_p>=seuil),which(popf_p>=seuil)]
r=r[which(Nini>=seuil)]
handling=handling[Nini>seuil]
efficience=efficience[Nini>seuil]
popf_h=popf_h[which(popf_h>=seuil)]
popf_p=popf_p[which(popf_p>=seuil)]
mati=sqrt(popf_p%*%t(popf_h))
resf=mati*I
nbsp_p=length(popf_p)
nbsp_h=length(popf_h)
Nini=c(popf_h,popf_p)

I=matrix(I,ncol=nbsp_h_per,nrow=nbsp_p_per)
Mplant=matrix(Mplant,ncol=nbsp_p_per,nrow=nbsp_p_per)
Mbird=matrix(Mbird,ncol=nbsp_h_per,nrow=nbsp_h_per)
CSh=matrix(CSh,ncol=nbsp_h_per,nrow=nbsp_h_per)
CSp=matrix(CSp,ncol=nbsp_p_per,nrow=nbsp_p_per)

final1$interactions_eq[final1$N_eq1>=seuil]=c(if(nbsp_h>1){apply(I,2,sum)}else{sum(I)},if(nbsp_p>1){apply(I,1,sum)}else{sum(I)})
final1$pheno_inter_eq[final1$N_eq1>=seuil]=NA
final1$pheno_intra_eq[final1$N_eq1>=seuil]=NA

omegab=matrix(0,nbsp_h,nbsp_h)
for(i in 1:nbsp_h){
I2=I[,i]*I
I2=matrix(I2,ncol=nbsp_h,nrow=nbsp_p)
vec=(apply(I2,2,sum))/sum(I[,i])
omegab[i,]=vec
}

omegap=matrix(0,nbsp_p,nbsp_p)
for(i in 1:nbsp_p){
I2=t(I)[,i]*t(I)
I2=matrix(I2,ncol=nbsp_p,nrow=nbsp_h)
vec=(apply(I2,2,sum))/sum(I[i,])
omegap[i,]=vec
}
omegab=omegab[which(popf_h>=seuil),which(popf_h>=seuil)]
omegap=omegap[which(popf_p>=seuil),which(popf_p>=seuil)]
omegab=matrix(omegab,ncol=nbsp_h_per,nrow=nbsp_h_per)
omegap=matrix(omegap,ncol=nbsp_p_per,nrow=nbsp_p_per)

mod=tryCatch({computeModules(I, method="Beckett")@likelihood},error = function(e) {NA})
mod2=tryCatch({computeModules(resf, method="Beckett")@likelihood},error = function(e) {NA})
jacob=rootSolve::jacobian.full(Nini,derivs,pert = 1e-6)
eig=eigen(jacob,only.values=T)$values

inv=-1*solve(jacob)
jacob2=jacob
inv2=inv
for(zz in 1:nrow(jacob)){
for(zzz in 1:ncol(jacob)){
inv2[zz,zzz]=inv[zz,zzz]/(inv[zz,zz]*inv[zzz,zzz]-inv[zz,zzz]*inv[zzz,zz])
}}
diag(jacob2)=NA
diag(inv2)=NA

ai_direct_hh=mean(jacob2[(1:nbsp_h),(1:nbsp_h)],na.rm=T)
ai_direct_pp=mean(jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)],na.rm=T)
ai_direct_ph=mean(jacob2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)],na.rm=T)
ai_direct_hp=mean(jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)],na.rm=T)
ai_net_hp=mean(inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)],na.rm=T)
ai_net_hh=mean(inv2[(1:nbsp_h),(1:nbsp_h)],na.rm=T)
ai_net_pp=mean(inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)],na.rm=T)
ai_net_ph=mean(inv2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)],na.rm=T)
ai_indirect_ph=mean(inv2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)]-jacob2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)],na.rm=T)
ai_indirect_hp=mean(inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)]-jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)],na.rm=T)
ai_indirect_hh=mean(inv2[(1:nbsp_h),(1:nbsp_h)]-jacob2[(1:nbsp_h),(1:nbsp_h)],na.rm=T)
ai_indirect_pp=mean(inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)]-jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)],na.rm=T)
ai_contrib_ph=mean((inv2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)]-jacob2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)])/(abs(inv2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)]-jacob2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)])+abs(jacob2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)])),na.rm=T)
ai_contrib_hp=mean((inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)]-jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)])/(abs(inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)]-jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)])+abs(jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)])),na.rm=T)
ai_contrib_hh=mean((inv2[(1:nbsp_h),(1:nbsp_h)]-jacob2[(1:nbsp_h),(1:nbsp_h)])/(abs(inv2[(1:nbsp_h),(1:nbsp_h)]-jacob2[(1:nbsp_h),(1:nbsp_h)])+abs(jacob2[(1:nbsp_h),(1:nbsp_h)])),na.rm=T)
ai_contrib_pp=mean((inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)]-jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)])/(abs(inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)]-jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)])+abs(jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)])),na.rm=T)

if(nbsp_h>1 & nbsp_p>1){
gainh=mean(apply(jacob2,1,mean,na.rm=T)[1:nbsp_h])
gainp=mean(apply(jacob2,1,mean,na.rm=T)[(1+nbsp_h):(nbsp_h+nbsp_p)])
gaintoth=mean(apply(inv2,1,mean,na.rm=T)[1:nbsp_h])
gaintotp=mean(apply(inv2,1,mean,na.rm=T)[(1+nbsp_h):(nbsp_h+nbsp_p)])
}else{
gainh=NA
gainp=NA
gaintoth=NA
gaintotp=NA
}

resmodel1=rbind(resmodel1,data.frame(nbsp_p=nbsp_p_dep,nbsp_h=nbsp_h_dep,nbsp_h_per=nbsp_h_per/nbsp_h_dep,nbsp_p_per=nbsp_p_per/nbsp_p_dep,
interf=interff,model="model1",essai=jj,variance=a,valprop=max(Re(eig)),connectance_I=mean(I),
connectance_r=if(nbsp_h>1 & nbsp_p>1){networklevel(resf,index="weighted connectance")[[1]]}else{NA},
NODF_I=if(nbsp_h>1 & nbsp_p>1){networklevel(I,index="NODF")[[1]]}else{NA},
mod_I=mod,cs=cs,shape_h=shape_h,shape_p=shape_p,NODF_r=if(nbsp_h>1 & nbsp_p>1){networklevel(resf,index="weighted NODF")[[1]]}else{NA},mod_r=mod2,gainp=gainp,gainh=gainh,gaintotp=gaintotp,gaintoth=gaintoth,
ai_direct_hh=ai_direct_hh,ai_direct_pp=ai_direct_pp,ai_direct_ph=ai_direct_ph,
ai_direct_hp=ai_direct_hp,ai_net_hp=ai_net_hp,ai_net_hh=ai_net_hh,ai_net_pp=ai_net_pp,ai_net_ph=ai_net_ph,
ai_indirect_hh=ai_indirect_hh,ai_indirect_pp=ai_indirect_pp,ai_indirect_ph=ai_indirect_ph,ai_indirect_hp=ai_indirect_hp,
ai_contrib_hp=ai_contrib_hp,ai_contrib_hh=ai_contrib_hh,ai_contrib_pp=ai_contrib_pp,ai_contrib_ph=ai_contrib_ph,efficience=efficience2,compb=mean(omegab),compp=mean(omegap),compt=mean(c(omegab,omegap))))
}else{
resmodel1=rbind(resmodel1,data.frame(nbsp_p=nbsp_p_dep,nbsp_h=nbsp_h_dep,nbsp_h_per=nbsp_h_per/nbsp_h_dep,nbsp_p_per=nbsp_p_per/nbsp_p_dep,interf=interff,model="model1",essai=jj,variance=a,
valprop=NA,connectance_I=NA,connectance_r=NA,NODF_I=NA,mod_I=NA,cs=cs,shape_h=shape_h,shape_p=shape_p,NODF_r=NA,mod_r=NA,gainp=NA,gainh=NA,gaintotp=NA,gaintoth=NA,
ai_direct_hh=NA,ai_direct_pp=NA,ai_direct_ph=NA,
ai_direct_hp=NA,ai_net_hp=NA,ai_net_hh=NA,ai_net_pp=NA,ai_net_ph=NA,
ai_indirect_hh=NA,ai_indirect_pp=NA,ai_indirect_ph=NA,ai_indirect_hp=NA,
ai_contrib_hp=NA,ai_contrib_hh=NA,ai_contrib_pp=NA,ai_contrib_ph=NA,efficience=efficience2,compb=NA,compp=NA,compt=NA))
}
final=rbind(final,final1)
}


### MODEL 2: with seasonal structure
#fonctions
derivs <-function(t, y,parms){
dy=rep(0,1)
for(i in 1:length(Nini)){
if(i<=nbsp_h){
I2=Ibin[,i]*Ibin*Mbird_p[,i,]*y[(1+nbsp_h):(nbsp_h+nbsp_p)]
I2=matrix(I2,ncol=nbsp_h,nrow=nbsp_p)
vec=apply(I2,2,sum)/sum(y[(1+nbsp_h):(nbsp_h+nbsp_p)]*I[,i])
vec[which(is.na(vec))]=0
vec[which(is.infinite(vec))]=0
eq1=(r[i]+
efficience[i]*sum(I[,i]*y[(1+nbsp_h):(nbsp_h+nbsp_p)])/(1+handling[i]*sum(I[,i]*y[(1+nbsp_h):(nbsp_h+nbsp_p)])+
interfp*sum(vec*y[1:nbsp_h]))-
sum(CSh[,i]*y[1:nbsp_h]))*y[i]

myenv=new.env(i)
environment(eq1)=myenv
myenv$i=i
myenv$Nini=Nini
myenv$I=I
myenv$I2=I2
myenv$vec=vec
myenv$r=r
myenv$nbsp_h=nbsp_h
myenv$nbsp_p=nbsp_p
myenv$cs=cs
myenv$interfp=interfp
}else{
I2=t(Ibin)[,(i-nbsp_h)]*t(Ibin)*Mplant_p[,(i-nbsp_h),]*y[1:nbsp_h]
I2=matrix(I2,ncol=nbsp_p,nrow=nbsp_h)
vec=apply(I2,2,sum)/sum(y[1:nbsp_h]*I[(i-nbsp_h),])
vec[which(is.na(vec))]=0
vec[which(is.infinite(vec))]=0
eq1=(r[i]+
efficience[i]*sum(I[(i-nbsp_h),]*y[1:nbsp_h])/(1+handling[i]*sum(I[(i-nbsp_h),]*y[1:nbsp_h])+
interff*sum(vec*y[(1+nbsp_h):(nbsp_h+nbsp_p)]))-
sum(CSp[,(i-nbsp_h)]*y[(1+nbsp_h):(nbsp_h+nbsp_p)]))*
y[i]

myenv=new.env(i)
environment(eq1)=myenv
myenv$i=i
myenv$I=I
myenv$I2=I2
myenv$Nini=Nini
myenv$vec=vec
myenv$r=r
myenv$nbsp_h=nbsp_h
myenv$nbsp_p=nbsp_p
myenv$cs=cs
myenv$interff=interff
}
if(y[i]>=seuil){dy[i]=eq1}else{dy[i]=0}}
return(list(dy))}

#run it for each value of competition:
for(j in comp_vec){
interff=j
interfp=j
load(paste(sitechoisi,jj,".RData",sep="_"))
efficience=rep(efficience2,nbsp_h+nbsp_p)
initial=Nini
r[r>(-1e-2) & r<0]=-1e-2
r[r<(1e-2) & r>0]=1e-2
Ibin=I
I=I*Phen
CSp=matrix(cs,ncol=nbsp_p,nrow=nbsp_p)
CSh=matrix(cs,ncol=nbsp_h,nrow=nbsp_h)
diag(CSp)=1
diag(CSh)=1

#Solve system until convergence
a=1
b=0
while(a>conv){
if(a==1){out <- ode(y=Nini, times=seq(0,20,0.5),func = derivs,parms = NULL,method="lsoda")}else{
out <- ode(y=Nini, times=seq(0,20,1), func = derivs,parms = NULL,method="lsoda")}
colnames(out)[1]="Time"
# if(a==1){out=integrate_sys(derivs2,Nini,20,0.05)}else{out=integrate_sys(derivs2,Nini,20,1)}
Nini=t(out)[2:(1+nbsp_h+nbsp_p),nrow(out)]
Nini[which(Nini<seuil)]=0
if(is.na(mean(Nini))){a=-2}else{
if(b==0){dyna=out}else{
out[,"Time"]=out[,"Time"]+b
dyna=rbind(dyna,out[-1,])}
if(nrow(dyna)>10){a=max(apply(t(dyna[(nrow(dyna)-10):nrow(dyna),2:(1+nbsp_h+nbsp_p)]),1,var))}
if(max(Nini,na.rm=T)>1e5){a=-1}
b=b+20
}
}

#Store basic data at species level
Nini[is.na(Nini)]=0
popf_p=Nini[(nbsp_h+1):(nbsp_h+nbsp_p)]
popf_h=Nini[1:nbsp_h]
nbsp_p_per=length(popf_p[popf_p>seuil])
nbsp_h_per=length(popf_h[popf_h>seuil])
final2=data.frame(species=c(colnames(I),rownames(I)),type=c(rep("H",nbsp_h),rep("P",nbsp_p)),Nini=initial,rzero=r,cs=cs,efficience=efficience,handling=handling,interactions=c(apply(I,2,sum),apply(I,1,sum)),
pheno_inter=c(apply(Phen,2,sum),apply(Phen,1,sum)),pheno_intra=c(apply(Mbird,2,sum),apply(Mplant,2,sum)),comp=c(apply(omegab_phen,1,sum),apply(omegap_phen,1,sum)),interf=c(rep(interfp,nbsp_h),rep(interff,nbsp_p)),essai=jj,N_eq1=c(popf_h,popf_p),model="model2",interactions_eq=NA,pheno_inter_eq=NA,pheno_intra_eq=NA,func_res=NA)

#Estimate saturation of functional response for each species:
for(i in 1:length(Nini)){
if(i<=nbsp_h){
I2=Ibin[,i]*Ibin*Mbird_p[,i,]*Nini[(1+nbsp_h):(nbsp_h+nbsp_p)]
I2=matrix(I2,ncol=nbsp_h,nrow=nbsp_p)
vec=apply(I2,2,sum)/sum(Nini[(1+nbsp_h):(nbsp_h+nbsp_p)]*I[,i])
vec[which(is.na(vec))]=0
vec[which(is.infinite(vec))]=0
final2$func_res[i]=efficience[i]*sum(I[,i]*Nini[(1+nbsp_h):(nbsp_h+nbsp_p)])/(1+handling[i]*sum(I[,i]*Nini[(1+nbsp_h):(nbsp_h+nbsp_p)])+
interfp*sum(vec*Nini[1:nbsp_h]))
}else{
I2=t(Ibin)[,(i-nbsp_h)]*t(Ibin)*Mplant_p[,(i-nbsp_h),]*Nini[1:nbsp_h]
I2=matrix(I2,ncol=nbsp_p,nrow=nbsp_h)
vec=apply(I2,2,sum)/sum(Nini[1:nbsp_h]*I[(i-nbsp_h),])
vec[which(is.na(vec))]=0
vec[which(is.infinite(vec))]=0
final2$func_res[i]=efficience[i]*sum(I[(i-nbsp_h),]*Nini[1:nbsp_h])/(1+handling[i]*sum(I[(i-nbsp_h),]*Nini[1:nbsp_h])+
interff*sum(vec*Nini[(1+nbsp_h):(nbsp_h+nbsp_p)]))
}}

#If network is not empty then store information at network level
if(nbsp_h_per>0 & nbsp_p_per>0){
### AGREGER:
I=I[which(popf_p>=seuil),which(popf_h>=seuil)]
Ibin=Ibin[which(popf_p>=seuil),which(popf_h>=seuil)]
Phen=Phen[which(popf_p>=seuil),which(popf_h>=seuil)]
CSp=CSp[which(popf_p>=seuil),which(popf_p>=seuil)]
CSh=CSh[which(popf_h>=seuil),which(popf_h>=seuil)]
Mbird=Mbird[which(popf_h>=seuil),which(popf_h>=seuil)]
Mbird_p=Mbird_p[which(popf_p>=seuil),which(popf_h>=seuil),which(popf_h>=seuil)]
Mplant_p=Mplant_p[which(popf_h>=seuil),which(popf_p>=seuil),which(popf_p>=seuil)]
Mplant=Mplant[which(popf_p>=seuil),which(popf_p>=seuil)]
r=r[which(Nini>=seuil)]
handling=handling[Nini>seuil]
efficience=efficience[Nini>seuil]
popf_h=popf_h[which(popf_h>=seuil)]
popf_p=popf_p[which(popf_p>=seuil)]
mati=sqrt(popf_p%*%t(popf_h))
resf=mati*I
nbsp_p=length(popf_p)
nbsp_h=length(popf_h)
Nini=c(popf_h,popf_p)

I=matrix(I,ncol=nbsp_h_per,nrow=nbsp_p_per)
Ibin=matrix(Ibin,ncol=nbsp_h_per,nrow=nbsp_p_per)
Phen=matrix(Phen,ncol=nbsp_h_per,nrow=nbsp_p_per)
Mplant=matrix(Mplant,ncol=nbsp_p_per,nrow=nbsp_p_per)
Mbird=matrix(Mbird,ncol=nbsp_h_per,nrow=nbsp_h_per)
Mbird_p=array(Mbird_p,dim=c(nbsp_p,nbsp_h,nbsp_h))
Mplant_p=array(Mplant_p,dim=c(nbsp_h,nbsp_p,nbsp_p))
CSh=matrix(CSh,ncol=nbsp_h_per,nrow=nbsp_h_per)
CSp=matrix(CSp,ncol=nbsp_p_per,nrow=nbsp_p_per)

final2$interactions_eq[final2$N_eq1>=seuil]=c(if(nbsp_h>1){apply(I,2,sum)}else{sum(I)},if(nbsp_p>1){apply(I,1,sum)}else{sum(I)})
final2$pheno_inter_eq[final2$N_eq1>=seuil]=c(if(nbsp_h>1){apply(Phen,2,sum)}else{sum(Phen)},if(nbsp_p>1){apply(Phen,1,sum)}else{sum(Phen)})
final2$pheno_intra_eq[final2$N_eq1>=seuil]=c(if(nbsp_h>1){apply(Mbird,2,sum)}else{sum(Mbird)},if(nbsp_p>1){apply(Mplant,2,sum)}else{sum(Mplant)})

omegab=matrix(0,nbsp_h,nbsp_h)
for(i in 1:nbsp_h){
I2=Ibin[,i]*Ibin*Mbird_p[,i,]
I2=matrix(I2,ncol=nbsp_h,nrow=nbsp_p)
vec=(apply(I2,2,sum))/sum(I[,i])
omegab[i,]=vec
}

omegap=matrix(0,nbsp_p,nbsp_p)
for(i in 1:nbsp_p){
I2=t(Ibin)[,i]*t(Ibin)*Mplant_p[,i,]
I2=matrix(I2,ncol=nbsp_p,nrow=nbsp_h)
vec=apply(I2,2,sum)/sum(I[i,])
omegap[i,]=vec
}
omegab=omegab[which(popf_h>=seuil),which(popf_h>=seuil)]
omegap=omegap[which(popf_p>=seuil),which(popf_p>=seuil)]
omegab=matrix(omegab,ncol=nbsp_h_per,nrow=nbsp_h_per)
omegap=matrix(omegap,ncol=nbsp_p_per,nrow=nbsp_p_per)

mod=tryCatch({computeModules(I, method="Beckett")@likelihood},error = function(e) {NA})
mod2=tryCatch({computeModules(resf, method="Beckett")@likelihood},error = function(e) {NA})
jacob=rootSolve::jacobian.full(Nini,derivs,pert = 1e-6)
eig=eigen(jacob,only.values=T)$values

inv=-1*solve(jacob)
jacob2=jacob
inv2=inv
for(zz in 1:nrow(jacob)){
for(zzz in 1:ncol(jacob)){
inv2[zz,zzz]=inv[zz,zzz]/(inv[zz,zz]*inv[zzz,zzz]-inv[zz,zzz]*inv[zzz,zz])
}}
diag(jacob2)=NA
diag(inv2)=NA

ai_direct_hh=mean(jacob2[(1:nbsp_h),(1:nbsp_h)],na.rm=T)
ai_direct_pp=mean(jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)],na.rm=T)
ai_direct_ph=mean(jacob2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)],na.rm=T)
ai_direct_hp=mean(jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)],na.rm=T)
ai_net_hp=mean(inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)],na.rm=T)
ai_net_hh=mean(inv2[(1:nbsp_h),(1:nbsp_h)],na.rm=T)
ai_net_pp=mean(inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)],na.rm=T)
ai_net_ph=mean(inv2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)],na.rm=T)
ai_indirect_ph=mean(inv2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)]-jacob2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)],na.rm=T)
ai_indirect_hp=mean(inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)]-jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)],na.rm=T)
ai_indirect_hh=mean(inv2[(1:nbsp_h),(1:nbsp_h)]-jacob2[(1:nbsp_h),(1:nbsp_h)],na.rm=T)
ai_indirect_pp=mean(inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)]-jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)],na.rm=T)
ai_contrib_ph=mean((inv2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)]-jacob2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)])/(abs(inv2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)]-jacob2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)])+abs(jacob2[(1:nbsp_h),(nbsp_h+1):(nbsp_h+nbsp_p)])),na.rm=T)
ai_contrib_hp=mean((inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)]-jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)])/(abs(inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)]-jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)])+abs(jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(1:nbsp_h)])),na.rm=T)
ai_contrib_hh=mean((inv2[(1:nbsp_h),(1:nbsp_h)]-jacob2[(1:nbsp_h),(1:nbsp_h)])/(abs(inv2[(1:nbsp_h),(1:nbsp_h)]-jacob2[(1:nbsp_h),(1:nbsp_h)])+abs(jacob2[(1:nbsp_h),(1:nbsp_h)])),na.rm=T)
ai_contrib_pp=mean((inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)]-jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)])/(abs(inv2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)]-jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)])+abs(jacob2[(nbsp_h+1):(nbsp_h+nbsp_p),(nbsp_h+1):(nbsp_h+nbsp_p)])),na.rm=T)

if(nbsp_h>1 & nbsp_p>1){
gainh=mean(apply(jacob2,1,mean,na.rm=T)[1:nbsp_h])
gainp=mean(apply(jacob2,1,mean,na.rm=T)[(1+nbsp_h):(nbsp_h+nbsp_p)])
gaintoth=mean(apply(inv2,1,mean,na.rm=T)[1:nbsp_h])
gaintotp=mean(apply(inv2,1,mean,na.rm=T)[(1+nbsp_h):(nbsp_h+nbsp_p)])
}else{
gainh=NA
gainp=NA
gaintoth=NA
gaintotp=NA
}

resmodel2=rbind(resmodel2,data.frame(nbsp_p=nbsp_p_dep,nbsp_h=nbsp_h_dep,nbsp_h_per=nbsp_h_per/nbsp_h_dep,nbsp_p_per=nbsp_p_per/nbsp_p_dep,
interf=interff,model="model2",essai=jj,variance=a,valprop=max(Re(eig)),connectance_I=mean(I),
connectance_r=if(nbsp_h>1 & nbsp_p>1){networklevel(resf,index="weighted connectance")[[1]]}else{NA},
NODF_I=if(nbsp_h>1 & nbsp_p>1){networklevel(I,index="NODF")[[1]]}else{NA},
mod_I=mod,cs=cs,shape_h=shape_h,shape_p=shape_p,NODF_r=if(nbsp_h>1 & nbsp_p>1){networklevel(resf,index="weighted NODF")[[1]]}else{NA},mod_r=mod2,gainp=gainp,gainh=gainh,gaintotp=gaintotp,gaintoth=gaintoth,
ai_direct_hh=ai_direct_hh,ai_direct_pp=ai_direct_pp,ai_direct_ph=ai_direct_ph,
ai_direct_hp=ai_direct_hp,ai_net_hp=ai_net_hp,ai_net_hh=ai_net_hh,ai_net_pp=ai_net_pp,ai_net_ph=ai_net_ph,
ai_indirect_hh=ai_indirect_hh,ai_indirect_pp=ai_indirect_pp,ai_indirect_ph=ai_indirect_ph,ai_indirect_hp=ai_indirect_hp,
ai_contrib_hp=ai_contrib_hp,ai_contrib_hh=ai_contrib_hh,ai_contrib_pp=ai_contrib_pp,ai_contrib_ph=ai_contrib_ph,efficience=efficience2,compb=mean(omegab),compp=mean(omegap),compt=mean(c(omegab,omegap))))
}else{
resmodel2=rbind(resmodel2,data.frame(nbsp_p=nbsp_p_dep,nbsp_h=nbsp_h_dep,nbsp_h_per=nbsp_h_per/nbsp_h_dep,nbsp_p_per=nbsp_p_per/nbsp_p_dep,interf=interff,model="model2",essai=jj,variance=a,
valprop=NA,connectance_I=NA,connectance_r=NA,NODF_I=NA,mod_I=NA,cs=cs,shape_h=shape_h,shape_p=shape_p,NODF_r=NA,mod_r=NA,gainp=NA,gainh=NA,gaintotp=NA,gaintoth=NA,
ai_direct_hh=NA,ai_direct_pp=NA,ai_direct_ph=NA,
ai_direct_hp=NA,ai_net_hp=NA,ai_net_hh=NA,ai_net_pp=NA,ai_net_ph=NA,
ai_indirect_hh=NA,ai_indirect_pp=NA,ai_indirect_ph=NA,ai_indirect_hp=NA,
ai_contrib_hp=NA,ai_contrib_hh=NA,ai_contrib_pp=NA,ai_contrib_ph=NA,efficience=efficience2,compb=NA,compp=NA,compt=NA))
}
final=rbind(final,final2)
}
}

resmodel=rbind(resmodel1,resmodel2)
return(list(resmodel,final))
}
setwd(dir="/home/duchenne/EPHI")
fwrite(resultat[[1]],paste0("res_net_",sitechoisi,param,"_cs_",cs,".txt"),row.names=F,sep="\t")
fwrite(resultat[[2]],paste0("resspecies_",sitechoisi,param,"_cs_",cs,".txt"),row.names=F,sep="\t")
}



closeCluster(cl)
mpi.quit()
