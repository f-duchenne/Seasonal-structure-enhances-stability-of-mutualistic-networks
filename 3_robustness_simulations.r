#setwd(dir="C:/Users/Duchenne/Documents/EPHI/temporal-dynamics/initial")
library(odeintr)
library(data.table)
library(dplyr)
library(doParallel)
library(foreach)
library(parallel)
library(doMPI)
cl<-startMPIcluster()
registerDoMPI(cl)
#cl<-makeCluster(96)
#registerDoParallel(cl)

conv=1e-14
seuil=1e-5
fini=250
cs_vec=seq(0,0.03,0.01)
comp_vec=seq(0,3,0.5)
alpha_vec=seq(0,3,0.5)
cs=0

intdat_raw <- fread("data_for_modelo.txt") #data accessible from the Datadryad

param="_c_efficience"

for(site in unique(intdat_raw$site)){ #loop over sites
sitechoisi=site
finali=fread(paste0("resspecies_",sitechoisi,param,"_cs_",cs,".txt"),sep="\t")

resultat=foreach(jj=1:fini,.combine=rbind)%dopar%{ #prallelized loop over replicates
library(deSolve)
library(rootSolve)

set.seed(jj)

finalr=NULL

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


for(j in comp_vec){ #run it for each value of competition:
interff=j
interfp=j
final=subset(finali,interf==j & efficience==efficience2 & model=="model1" & essai==jj)
load(paste(sitechoisi,jj,".RData",sep="_"))
efficience=rep(efficience2,nbsp_h+nbsp_p)
popf_h=final$N_eq1[final$type=="H"]
popf_p=final$N_eq1[final$type=="P"]
initial=final$N_eq1
Nini=initial
r[r>(-1e-2) & r<0]=-1e-2
r[r<(1e-2) & r>0]=1e-2
CSp=matrix(cs,ncol=nbsp_p,nrow=nbsp_p)
CSh=matrix(cs,ncol=nbsp_h,nrow=nbsp_h)
diag(CSp)=1
diag(CSh)=1
nbsp_p_per=length(final$species[final$N_eq1>=1e-5 & final$type=="P"])
nbsp_h_per=length(final$species[final$N_eq1>=1e-5 & final$type=="H"])

if(nbsp_h_per>0 & nbsp_p_per>0){
I=I[which(popf_p>=seuil),which(popf_h>=seuil)]
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
P_eq1=sum(popf_p)
H_eq1=sum(popf_h)
nbsp_p=length(popf_p)
nbsp_h=length(popf_h)
Nini=c(popf_h,popf_p)
initial=Nini

I=matrix(I,ncol=nbsp_h_per,nrow=nbsp_p_per)
Phen=matrix(Phen,ncol=nbsp_h_per,nrow=nbsp_p_per)
Mplant=matrix(Mplant,ncol=nbsp_p_per,nrow=nbsp_p_per)
Mbird=matrix(Mbird,ncol=nbsp_h_per,nrow=nbsp_h_per)
Mbird_p=array(Mbird_p,dim=c(nbsp_p,nbsp_h,nbsp_h))
Mplant_p=array(Mplant_p,dim=c(nbsp_h,nbsp_p,nbsp_p))
CSh=matrix(CSh,ncol=nbsp_h_per,nrow=nbsp_h_per)
CSp=matrix(CSp,ncol=nbsp_p_per,nrow=nbsp_p_per)

for(z in 0:nbsp_h){ #simulate extinction of one hummingbird species
Nini=initial
if(z>0){Nini[z]=0}

#Solve system until convergence
a=1
b=0
while(a>conv){
if(a==1){out <- ode(y=Nini, times=seq(0,20,0.05), func = derivs,parms = NULL,method="lsoda")}else{
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
# matplot(log(dyna[,2:(1+nbsp_h)]),type="l",ylab="abond")
# matplot(log(dyna[,(2+nbsp_h):(1+nbsp_h+nbsp_p)]),type="l",ylab="abond")

Nini[is.na(Nini)]=0
popf_p=Nini[(nbsp_h+1):(nbsp_h+nbsp_p)]
popf_h=Nini[1:nbsp_h]
nbsp_p_per2=length(popf_p[popf_p>seuil])
nbsp_h_per2=length(popf_h[popf_h>seuil])

final1=data.frame(species_remov=ifelse(z>0,final$species[final$N_eq1>=seuil & final$type=="H"][z],NA),nbsp_p_per1=nbsp_p_per,nbsp_h_per1=nbsp_h_per,
nbsp_p_per2=nbsp_p_per2,nbsp_h_per2=nbsp_h_per2,abond_remov=ifelse(z>0,initial[z],NA),P_eq1=P_eq1,H_eq1=H_eq1,P_eq2=sum(popf_p),H_eq2=sum(popf_h),
cs=cs,essai=jj,model="model1",efficience=efficience2,interf=interfp,site=sitechoisi,type_remove=ifelse(z==0,NA,ifelse(z>nbsp_h_per,"P","H")))
finalr=rbind(finalr,final1)
}
}

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


for(j in comp_vec){ #run it for each value of competition
interff=j
interfp=j
final=subset(finali,interf==j & efficience==efficience2 & model=="model2" & essai==jj)
load(paste(sitechoisi,jj,".RData",sep="_"))
efficience=rep(efficience2,nbsp_h+nbsp_p)
popf_h=final$N_eq1[final$type=="H"]
popf_p=final$N_eq1[final$type=="P"]
initial=final$N_eq1
Nini=initial
r[r>(-1e-2) & r<0]=-1e-2
r[r<(1e-2) & r>0]=1e-2
CSp=matrix(cs,ncol=nbsp_p,nrow=nbsp_p)
CSh=matrix(cs,ncol=nbsp_h,nrow=nbsp_h)
diag(CSp)=1
diag(CSh)=1
nbsp_p_per=length(final$species[final$N_eq1>=1e-5 & final$type=="P"])
nbsp_h_per=length(final$species[final$N_eq1>=1e-5 & final$type=="H"])

if(nbsp_h_per>0 & nbsp_p_per>0){
I=I[which(popf_p>=seuil),which(popf_h>=seuil)]
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
P_eq1=sum(popf_p)
H_eq1=sum(popf_h)
nbsp_p=length(popf_p)
nbsp_h=length(popf_h)
Nini=c(popf_h,popf_p)
initial=Nini

I=matrix(I,ncol=nbsp_h_per,nrow=nbsp_p_per)
Phen=matrix(Phen,ncol=nbsp_h_per,nrow=nbsp_p_per)
Mplant=matrix(Mplant,ncol=nbsp_p_per,nrow=nbsp_p_per)
Mbird=matrix(Mbird,ncol=nbsp_h_per,nrow=nbsp_h_per)
Mbird_p=array(Mbird_p,dim=c(nbsp_p,nbsp_h,nbsp_h))
Mplant_p=array(Mplant_p,dim=c(nbsp_h,nbsp_p,nbsp_p))
CSh=matrix(CSh,ncol=nbsp_h_per,nrow=nbsp_h_per)
CSp=matrix(CSp,ncol=nbsp_p_per,nrow=nbsp_p_per)
Ibin=I
I=I*Phen


for(z in 0:nbsp_h){ #simulate extinction of one hummingbird species
Nini=initial
if(z>0){Nini[z]=0}

#Solve system until convergence
a=1
b=0
while(a>conv){
if(a==1){out <- ode(y=Nini, times=seq(0,20,0.05), func = derivs,parms = NULL,method="lsoda")}else{
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
# matplot(log(dyna[,2:(1+nbsp_h)]),type="l",ylab="abond")
# matplot(log(dyna[,(2+nbsp_h):(1+nbsp_h+nbsp_p)]),type="l",ylab="abond")

Nini[is.na(Nini)]=0
popf_p=Nini[(nbsp_h+1):(nbsp_h+nbsp_p)]
popf_h=Nini[1:nbsp_h]
nbsp_p_per2=length(popf_p[popf_p>seuil])
nbsp_h_per2=length(popf_h[popf_h>seuil])

final2=data.frame(species_remov=ifelse(z>0,final$species[final$N_eq1>=seuil & final$type=="H"][z],NA),nbsp_p_per1=nbsp_p_per,nbsp_h_per1=nbsp_h_per,
nbsp_p_per2=nbsp_p_per2,nbsp_h_per2=nbsp_h_per2,abond_remov=ifelse(z>0,initial[z],NA),P_eq1=P_eq1,H_eq1=H_eq1,P_eq2=sum(popf_p),H_eq2=sum(popf_h),
cs=cs,essai=jj,model="model2",efficience=efficience2,interf=interfp,site=sitechoisi,type_remove=ifelse(z==0,NA,ifelse(z>nbsp_h_per,"P","H")))
finalr=rbind(finalr,final2)
}
}
}
}
return(finalr)
}
setwd(dir="/home/duchenne/EPHI")
fwrite(resultat,paste0("robustesse_",sitechoisi,param,"_cs_",cs,".txt"))
}


closeCluster(cl)
mpi.quit()
