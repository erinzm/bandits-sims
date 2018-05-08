### Function that gives a 3 star score RV

threestar=function(p1,p2)
{
 fos=runif(1,0,1)
 return=0.5
 if(fos<p1){return=0}
 if(fos>(p1+p2)){return=1}
 return
}

### KL UCB functions

leftright=function(muhat,lower,upper,threshold)
{
 if((muhat*(1-muhat))!=0)
 {
  shit=(lower+upper)/2
  KL=(muhat*log(muhat/shit))+((1-muhat)*log((1-muhat)/(1-shit)))
  if(KL>=threshold){return=c(lower,shit,shit)}
  if(KL<threshold){return=c(shit,upper,shit)}
  return
 }
 if(muhat==0)
 {
  shit=(lower+upper)/2
  KL=((1-muhat)*log((1-muhat)/(1-shit)))
  if(KL>=threshold){return=c(lower,shit,shit)}
  if(KL<threshold){return=c(shit,upper,shit)}
  return
 }
 if(muhat==1)
 {
  return=c(1,1,1)
 }
 return
}

computeUCB=function(muhat,threshold,accuracy=10^(-6))
{
 lower=muhat
 upper=1
 return=1
 while((upper-lower)>accuracy)
 {
  new=leftright(muhat,lower,upper,threshold)
  lower=new[1]
  upper=new[2]
  return=new[3]
 }
 return
}

### KL UCB functions

D=function(x,y)
{
 return=0
 if(x!=y){return=x*log(x/y)+(1-x)*log((1-x)/(1-y))}
 if(x==0){return=log(1/(1-y))}
 if(x==1){return=log(1/y)}
 return
}

NEWleftright=function(N,muhat,lower,upper,threshold)
{
 if(muhat!=1)
 {
  m=(lower+upper)/2
  x=(N*muhat+m)/(N+1)
  KL=D(x,m)
  if(KL>threshold){return=c(lower,m,m)}
  if(KL<threshold){return=c(m,upper,m)}
  if(KL==threshold){return=c(m,m,m)}
  return
 }
 if(muhat==1)
 {
  return=c(1,1,1)
 }
 return
}

computeNEWUCB=function(N,muhat,threshold,accuracy=10^(-10))
{
 lower=muhat
 upper=1
 return=1
 while((upper-lower)>accuracy)
 {
  new=NEWleftright(N,muhat,lower,upper,threshold)
  lower=new[1]
  upper=new[2]
  return=new[3]
 }
 return
}

### Import data

setwd('c:/Uw - Madison/Notes, Write-ups/KL-UCB/Paper/Simulations')

data=read.csv("NYr_03-09-2017-KL.csv")
data=read.csv("NYr_03-09-2017-lil.csv")
data=read.csv("NYr_old.csv")
data=read.csv("NYr_05-11-2017-KL.csv")

data=read.csv("NYr_cont558-lil.csv")

### Set up arms.
### WATCH OUT! Data before contest 530 is stored differently!
### After 531, it's "probs[i,j]=(data[i,j+6]/data[i,6])"
### Before 530, it's "probs[i,j]=(data[i,j+1]/(data[i,2]+data[i,3]+data[i,4]))"
### Also, data not necessarily sorted before contest 526

n=length(data[,1]) # Number of arms
probs=mat.or.vec(n,2)
for(i in 1:n)
{
# probs[i,1]=(data[i,5]/(data[i,6])) # For 512-lil
 probs[i,1]=(data[i,7]/(data[i,6])) # For 558-lil
# probs[i,2]=(data[i,4]/(data[i,6])) ### We merge somewhat funny & funny in the simulations!
}

mu=rep(0,n)
for(i in 1:n){mu[i]=1-probs[i,1]-0.5*probs[i,2]}
plot(mu)
max(mu)

parametric=rep(0,n)
for(i in 1:n){parametric[i]=max(mu)*(1-((i-1)/n)^(1/3))}
points(parametric,col="blue")

### Sort ams (if necessary)

sortedprobs=mat.or.vec(n,2)
muhelper=mu
for(i in 1:n)
{
 top=which.is.max(muhelper)
 sortedprobs[i,]=probs[top,]
 muhelper[top]=-1
}
probs=sortedprobs
for(i in 1:n){mu[i]=1-probs[i,1]-0.5*probs[i,2]}
plot(mu)
max(mu)

parametric=rep(0,n)
for(i in 1:n){parametric[i]=max(mu)*(1-((i-1)/n)^(1/3))}
points(parametric,col="blue")

### Merge somewhat funny and funny (if necessary)

probs[,2]=rep(0,n)
for(i in 1:n){mu[i]=1-probs[i,1]-0.5*probs[i,2]}
plot(mu)
max(mu)

### Parametric arm setup

n=1000
probs=mat.or.vec(n,2)
alpha=1/2
mumax=1
for(i in 1:n){probs[i,1]=(1-mumax)+mumax*(((i-1)/n)^alpha)}

mu=rep(0,n)
for(i in 1:n){mu[i]=1-probs[i,1]-0.5*probs[i,2]}
plot(mu)
max(mu)


### The sub-Gaussian bounds and the functions to compute the thresholds for the KL-bound

kaufmann=function(t,d){sqrt(0.5*(log(1/d)+3*log(log(1/d))+1.5*log(log(exp(1)*t)))/t)} # Kaufmann's bound

ourSGbound=function(t,d,k,N){sqrt((((N+1)/N)^2/2)*log(k*log2(2*t)/d)/t)} # Our sub-Gaussian bound

KLthreshold=function(t,d,k){log(k*log2(2*t)/d)/t} # The KL-threshold

kappa=function(N,d)
{
 l=log2(N)
 upper=1000

 return=0
 for(j in 1:N)
 {
  s=0
  for(k in l:upper)
  {
   s=s+(k+1)^(-(N+1)/(N))
  }
  s=s+N*(upper+1)^(-1/N)
  return=return+s
  if(l!=0){return=return+(log2(2*j)^(-(N+1)/N))}
 }
 return=return*(d^(1/N))
 return=return^(N/(N+1))
 return
}

### Choose parameters for our confidence bounds

N=8
delta=0.01 # the error tolerance

k=kappa(N,delta) # the kappa constant

### Compute threshold sequences beforehand to save time (this is OK, if we don't have more pulls on any arm then 'maxpulls')

maxpulls=400000
kaufmannthr=rep(0,maxpulls)
for(i in 1:maxpulls){kaufmannthr[i]=kaufmann(i,delta)}
ourSGthr=rep(0,maxpulls)
for(i in 1:maxpulls){ourSGthr[i]=ourSGbound(i,delta,k,N)}
KLthr=rep(0,maxpulls)
for(i in 1:maxpulls){KLthr[i]=KLthreshold(i,delta,k)}

### Experiment parameters

numberofreps=250
accuracy=10^(-6)

m=5
horizon=50000
howoften=horizon/20 # how often do we evaluate whether or not the best arm is in the top m

kaufmannprob=rep(0,(horizon/howoften))
KLprob=rep(0,(horizon/howoften))
ourSGprob=rep(0,(horizon/howoften))
isitgood=0

kaufmannlowpulls=0
KLlowpulls=0
ourSGlowpulls=0

for(j in 1:numberofreps)
{

### Kaufmann lil-UCB

 pulls=rep(0,n) # number of pulls of each arm
 empmeans=rep(0,n) # means of arms
 ucb=rep(0,n) # upper conf bounds
 tiebreak=runif(n,0,1) # used to decide how we break ties in the ranking

 for(i in 1:n)
 {
  pulls[i]=pulls[i]+1
  X=threestar(probs[i,1],probs[i,2])
  empmeans[i]=(((pulls[i]-1)/pulls[i])*empmeans[i])+(X/pulls[i])
  ucb[i]=empmeans[i]+kaufmannthr[1]
 }

 totalpulls=n

 while(totalpulls<horizon)
 {
  where=which.is.max(ucb)
  pulls[where]=pulls[where]+1
  totalpulls=totalpulls+1

  X=threestar(probs[where,1],probs[where,2])
  empmeans[where]=(((pulls[where]-1)/pulls[where])*empmeans[where])+(X/pulls[where])
  ucb[where]=empmeans[where]+kaufmannthr[pulls[where]]

  if((totalpulls%%howoften)==0)
  {
   rankone=1
   for(a in 2:n)
   {
    if(empmeans[1]<empmeans[a]){rankone=rankone+1}
    if(empmeans[1]==empmeans[a]){if(tiebreak[1]<tiebreak[a]){rankone=rankone+1}}
    if(rankone>m){break}
   }
   if(rankone<(m+1)){isitgood=1}else{isitgood=0}
#   if(which.is.max(empmeans)==1){isitgood=1}else{isitgood=0}
   kaufmannprob[(totalpulls/howoften)]=((j-1)/j)*kaufmannprob[(totalpulls/howoften)]+(isitgood/j)
  }

  if(totalpulls==horizon){kaufmannlowpulls=((j-1)/j)*kaufmannlowpulls+(sum(pulls[ceiling(n/2):n])/j)}

 }

### lil-UCB with our SG-bound

 pulls=rep(0,n) # number of pulls of each arm
 empmeans=rep(0,n) # means of arms
 ucb=rep(0,n) # upper conf bounds
 tiebreak=runif(n,0,1) # used to decide how we break ties in the ranking

 for(i in 1:n)
 {
  pulls[i]=pulls[i]+1
  X=threestar(probs[i,1],probs[i,2])
  empmeans[i]=(((pulls[i]-1)/pulls[i])*empmeans[i])+(X/pulls[i])
  ucb[i]=empmeans[i]+ourSGthr[1]
 }

 totalpulls=n

 while(totalpulls<horizon)
 {
  where=which.is.max(ucb)
  pulls[where]=pulls[where]+1
  totalpulls=totalpulls+1

  X=threestar(probs[where,1],probs[where,2])
  empmeans[where]=(((pulls[where]-1)/pulls[where])*empmeans[where])+(X/pulls[where])
  ucb[where]=empmeans[where]+ourSGthr[pulls[where]]

  if((totalpulls%%howoften)==0)
  {
   rankone=1
   for(a in 2:n)
   {
    if(empmeans[1]<empmeans[a]){rankone=rankone+1}
    if(empmeans[1]==empmeans[a]){if(tiebreak[1]<tiebreak[a]){rankone=rankone+1}}
    if(rankone>m){break}
   }
   if(rankone<(m+1)){isitgood=1}else{isitgood=0}
#   if(which.is.max(empmeans)==1){isitgood=1}else{isitgood=0}
   ourSGprob[(totalpulls/howoften)]=((j-1)/j)*ourSGprob[(totalpulls/howoften)]+(isitgood/j)
  }

  if(totalpulls==horizon){ourSGlowpulls=((j-1)/j)*ourSGlowpulls+(sum(pulls[ceiling(n/2):n])/j)}

 }

### KL-UCB

 pulls=rep(0,n) # number of pulls of each arm
 empmeans=rep(0,n) # means of arms
 ucb=rep(0,n) # upper conf bounds
 tiebreak=runif(n,0,1) # used to decide how we break ties in the ranking

 for(i in 1:n)
 {
  pulls[i]=pulls[i]+1
  X=threestar(probs[i,1],probs[i,2])
  empmeans[i]=(((pulls[i]-1)/pulls[i])*empmeans[i])+(X/pulls[i])
  ucb[i]=computeNEWUCB(N,empmeans[i],KLthr[1],accuracy)
 }

 totalpulls=n

 while(totalpulls<horizon)
 {
  where=which.is.max(ucb)
  pulls[where]=pulls[where]+1
  totalpulls=totalpulls+1

  X=threestar(probs[where,1],probs[where,2])
  empmeans[where]=(((pulls[where]-1)/pulls[where])*empmeans[where])+(X/pulls[where])
  ucb[where]=computeNEWUCB(N,empmeans[where],KLthr[pulls[where]],accuracy)

  if((totalpulls%%howoften)==0)
  {
   rankone=1
   for(a in 2:n)
   {
    if(empmeans[1]<empmeans[a]){rankone=rankone+1}
    if(empmeans[1]==empmeans[a]){if(tiebreak[1]<tiebreak[a]){rankone=rankone+1}}
    if(rankone>m){break}
   }
   if(rankone<(m+1)){isitgood=1}else{isitgood=0}
#   if(which.is.max(empmeans)==1){isitgood=1}else{isitgood=0}
   KLprob[(totalpulls/howoften)]=((j-1)/j)*KLprob[(totalpulls/howoften)]+(isitgood/j)
  }

  if(totalpulls==horizon){KLlowpulls=((j-1)/j)*KLlowpulls+(sum(pulls[ceiling(n/2):n])/j)}

 }

}

### Plot

xax=seq(0,horizon,howoften)/10000
plot(xax,c(0,KLprob),type="l",col="blue",ylim=c(0,1),main="P(best arm in top 5)",sub="Parametric, alpha  1/3",xlab="Number of pulls (10 thousands)",ylab="(Empirical) probability - 250 trials",xaxt="n")
lines(xax,c(0,kaufmannprob))
lines(xax,c(0,ourSGprob),col="red")
legend("bottomright",c("Kaufmann lil-UCB","KL-UCB","Our SG lil-UCB"),fill=c('black','blue','red'))
# axis(1,at=seq(0,horizon,howoften*2)/howoften)
axis(1,at=xax)
abline(h=0.9,lty=3)

kaufmannlowpulls
ourSGlowpulls
KLlowpulls

### Export results

results=matrix(c(xax[2:length(xax)],-1,kaufmannprob,kaufmannlowpulls,ourSGprob,ourSGlowpulls,KLprob,KLlowpulls),ncol=4,nrow=(length(KLprob)+1))
write.table(results,"forpaper_250reps_alpha1_small_new.csv",row.names=F,col.names=c("pulls","Kaufmann","our SG","KL"),sep=",")










###
### Perform Thompson sampling
###


threestar=function(p1,p2)
{
 fos=runif(1,0,1)
 return=0.5
 if(fos<p1){return=0}
 if(fos>(p1+p2)){return=1}
 return
}


numberofreps=100
accuracy=10^(-6)

m=5
horizon=200000
howoften=horizon/20 # how often do we evaluate whether or not the best arm is in the top m

Thompsonprob=rep(0,(horizon/howoften))
isitgood=0

Thompsonlowpulls=0

for(j in 1:numberofreps)
{

### Thompson sampling

 theta=rep(0,n)
 a=1
 b=1
 alpha=1

 pulls=rep(0,n) # number of pulls of each arm
 empmeans=rep(0,n) # means of arms
 ucb=rep(0,n) # upper conf bounds
 tiebreak=runif(n,0,1) # used to decide how we break ties in the ranking

 for(i in 1:n)
 {
  pulls[i]=pulls[i]+1
  X=threestar(probs[i,1],probs[i,2])
  empmeans[i]=(((pulls[i]-1)/pulls[i])*empmeans[i])+(X/pulls[i])
 }

 totalpulls=n

 while(totalpulls<horizon)
 {
  if((totalpulls%%10)==0){for(i in 1:n){theta[i]=rbeta(1,((pulls[i]*empmeans[i])+a)/alpha,((pulls[i]*(1-empmeans[i]))+b)/alpha)}}

  where=which.is.max(theta)
  pulls[where]=pulls[where]+1
  totalpulls=totalpulls+1

  X=threestar(probs[where,1],probs[where,2])
  empmeans[where]=(((pulls[where]-1)/pulls[where])*empmeans[where])+(X/pulls[where])
  theta[where]=rbeta(1,((pulls[where]*empmeans[where])+a)/alpha,((pulls[where]*(1-empmeans[where]))+b)/alpha)

  if((totalpulls%%howoften)==0)
  {
   rankone=1
   for(a in 2:n)
   {
    if(empmeans[1]<empmeans[a]){rankone=rankone+1}
    if(empmeans[1]==empmeans[a]){if(tiebreak[1]<tiebreak[a]){rankone=rankone+1}}
    if(rankone>m){break}
   }
   if(rankone<(m+1)){isitgood=1}else{isitgood=0}
#   if(which.is.max(empmeans)==1){isitgood=1}else{isitgood=0}
   Thompsonprob[(totalpulls/howoften)]=((j-1)/j)*Thompsonprob[(totalpulls/howoften)]+(isitgood/j)
  }

  if(totalpulls==horizon){Thompsonlowpulls=((j-1)/j)*Thompsonlowpulls+(sum(pulls[ceiling(n/2):n])/j)}

 }

}


### Plot

xax=seq(0,horizon,howoften)/10000
plot(xax,c(0,Thompsonprob),type="l",col="blue",ylim=c(0,1),main="P(best arm in top 5)",sub="Parametric, alpha  1",xlab="Number of pulls (10 housands)",ylab="(Empirical) probability - 250 trials",xaxt="n")
axis(1,at=xax)
abline(h=0.9,lty=3)



### Export results

xax=seq(0,horizon,howoften)/10000

results=matrix(c(xax[2:length(xax)],-1,Thompsonprob,Thompsonlowpulls),ncol=2,nrow=(length(Thompsonprob)+1))
write.table(results,"forpaper_100reps_558lil_Thompson.csv",row.names=F,col.names=c("pulls","Thompson"),sep=",")














































### Get old results and merge them together

results1=read.csv("top5_N8_delta10-2_100reps_cont512-lil-merged.csv")

results2=read.csv("top5_N8_25reps_cont528_2.csv")

length(results1[,1])-length(results[,1])

aggresults=matrix(rep(0,length(results1[,1]*3)),ncol=3,nrow=length(results1[,1]))

for(i in 1:length(results1[,1]))
{
 for(j in 1:3)
 {
  aggresults[i,j]=(1/2)*results1[i,j+1]+(1/2)*results[i,j+1]
 }
}

kaufmannprob=aggresults[,1]
ourSGprob=aggresults[,2]
KLprob=aggresults[,3]

horizon=400000
howoften=20000 # how often do we evaluate whether or not the best arm is in the top m













n=5000 # Number of arms
probs=mat.or.vec(n,2)
for(i in 1:n)
{
 probs[i,1]=i/n
 probs[i,2]=0
}

mu=rep(0,n)
for(i in 1:n){mu[i]=1-probs[i,1]-0.5*probs[i,2]}
plot(mu)
max(mu)




###
### What do we gain by using KL versus Hoeffding? For one arm
###

### Computing the UCB using binary bisection. The core function is leftright.

leftright=function(muhat,lower,upper,threshold)
{
 if((muhat*(1-muhat))!=0)
 {
  shit=(lower+upper)/2
  KL=(muhat*log(muhat/shit))+((1-muhat)*log((1-muhat)/(1-shit)))
  if(KL>=threshold){return=c(lower,shit,(shit+lower)/2)}
  if(KL<threshold){return=c(shit,upper,(shit+upper)/2)}
  return
 }
 if(muhat==0)
 {
  shit=(lower+upper)/2
  KL=((1-muhat)*log((1-muhat)/(1-shit)))
  if(KL>=threshold){return=c(lower,shit,(shit+lower)/2)}
  if(KL<threshold){return=c(shit,upper,(shit+upper)/2)}
  return
 }
 if(muhat==1)
 {
  return=c(1,1,1)
 }
 return
}

computeUCB=function(muhat,threshold,accuracy=10^(-6))
{
 lower=muhat
 upper=1
 return=1
 while((upper-lower)>accuracy)
 {
  new=leftright(muhat,lower,upper,threshold)
  lower=new[1]
  upper=new[2]
  return=new[3]
 }
 return
}

### Import data

data=read.csv("NYr_03-09-2017-KL.csv")
data=read.csv("NYr_03-09-2017-lil.csv")
data=read.csv("NYr_old.csv")




data=read.csv("NYr_cont558-KL.csv")

### Set up arms

n=length(data[,1]) # Number of arms
probs=mat.or.vec(n,2)
for(i in 1:n)
{
 for(j in 1:2)
 {
  probs[i,j]=(data[i,j+6]/data[i,6])
 }
}

mu=rep(0,n)
for(i in 1:n){mu[i]=1-probs[i,1]-0.5*probs[i,2]}

### The sub-Gaussian bounds and the functions to compute the thresholds for the KL-bound

kaufmann=function(t,d){sqrt(0.5*(log(1/d)+3*log(log(1/d))+1.5*log(log(exp(1)*t/2)))/t)} # Kaufmann's bound

ourSGbound=function(t,d,k,N){sqrt((((N+1)/N)^2/2)*log(k*log2(2*t)/d)/t)} # Our sub-Gaussian bound

KLthreshold=function(t,d,k){log(k*log2(2*t)/d)/t} # The KL-threshold

kappa=function(N,d)
{
 l=log2(N)
 upper=1000

 return=0
 for(j in 1:N)
 {
  s=0
  for(k in l:upper)
  {
   s=s+(k+1)^(-(N+1)/(N))
  }
  s=s+N*(upper+1)^(-1/N)
  return=return+s
  if(l!=0){return=return+(log2(2*j)^(-(N+1)/N))}
 }
 return=return*(d^(1/N))
 return=return^(N/(N+1))
 return
}

### Compute threshold sequences beforehand to save time (this is OK, if we don't have more pulls on any arm then 'maxpulls')

delta=0.001 # the error tolerance

N=8
k=kappa(N,delta) # the kappa constant

maxpulls=50
SGthr=rep(0,maxpulls)
for(i in 1:maxpulls){SGthr[i]=sqrt(log((pi^2/6)*(i^2)/delta)/(2*i))}
KLthr=rep(0,maxpulls)
for(i in 1:maxpulls){KLthr[i]=log((pi^2/6)*(i^2)/delta)/i}
Kaufthr=rep(0,maxpulls)
for(i in 1:maxpulls){Kaufthr[i]=kaufmann(i,delta)}
itlogSGthr=rep(0,maxpulls)
for(i in 1:maxpulls){itlogSGthr[i]=ourSGbound(i,delta,k,N)}
itlogKLthr=rep(0,maxpulls)
for(i in 1:maxpulls){itlogKLthr[i]=KLthreshold(i,delta,k)}


### Choose an arm

index=1000
prob=probs[index,]

### The actual experiment. Note that I need the following functions: computeNEWUCB, NEWleftright, treestar

numberofreps=50
SGUCB=rep(0,maxpulls)
SGtomu=rep(0,maxpulls)
KLUCB=rep(0,maxpulls)
KaufUCB=rep(0,maxpulls)
itlogSGUCB=rep(0,maxpulls)
itlogKLUCB=rep(0,maxpulls)

for(j in 1:numberofreps)
{
 empmeans=rep(0,maxpulls)
 empmeans[1]=threestar(prob[1],prob[2])

 SGUCB[1]=((j-1)/j)*SGUCB[1]+(empmeans[1]+SGthr[1])/j
 KLupper=computeUCB(empmeans[1],KLthr[1])
 KLUCB[1]=((j-1)/j)*KLUCB[1]+KLupper/j
 KaufUCB[1]=((j-1)/j)*KaufUCB[1]+(empmeans[1]+Kaufthr[1])/j
 itlogSGUCB[1]=((j-1)/j)*itlogSGUCB[1]+(empmeans[1]+itlogSGthr[1])/j
 itlogKLupper=computeNEWUCB(N,empmeans[1],itlogKLthr[1])
 itlogKLUCB[1]=((j-1)/j)*itlogKLUCB[1]+itlogKLupper/j

 for(t in 2:maxpulls)
 {
  X=threestar(prob[1],prob[2])
  empmeans[t]=((t-1)/t)*empmeans[t-1]+X/t

  SGUCB[t]=((j-1)/j)*SGUCB[t]+(empmeans[t]+SGthr[t])/j
  KLupper=computeUCB(empmeans[t],KLthr[t])
  KLUCB[t]=((j-1)/j)*KLUCB[t]+KLupper/j
  KaufUCB[t]=((j-1)/j)*KaufUCB[t]+(empmeans[t]+Kaufthr[t])/j
  itlogSGUCB[t]=((j-1)/j)*itlogSGUCB[t]+(empmeans[t]+itlogSGthr[t])/j
  itlogKLupper=computeNEWUCB(N,empmeans[t],itlogKLthr[t])
  itlogKLUCB[t]=((j-1)/j)*itlogKLUCB[t]+itlogKLupper/j

 }
}

### Plot

plot(empmeans,type="l",ylim=c(0,1))
lines(SGUCB,col="red",lty=2)
lines(KLUCB,col="blue",lty=2)
lines(KaufUCB,col="green")
lines(itlogSGUCB,col="red")
lines(itlogKLUCB,col="blue")
abline(h=mu[index],lty=2)
abline(h=mu[1],lty=3)
legend("topright",c("Kaufmann","KL (solid: log log; dash: union bound)","Our SG (solid: log log; dash: union bound)"),fill=c('green','blue','red'))


## Plot ratios to mu

ratio=rep(0,maxpulls)
for(i in 1:maxpulls){ratio[i]=(KaufUCB[i]-mu)/(itlogKLUCB[i]-mu)}
plot(ratio,type="l")
abline(h=1,lty=2)
















###
### Dueling bandits?
###

data=read.csv("508-round2-dueling-responses.csv")

noofqueries=length(data[,1])
queries=mat.or.vec(noofqueries,2) # This is the list of queries int he format: the first item beat the second item
for(i in 1:noofqueries)
{
 if(data[i,4]==data[i,13]){queries[i,1]=data[i,4];queries[i,2]=data[i,10]}else{queries[i,1]=data[i,10];queries[i,2]=data[i,4]}
}

n=max(queries) # Number of arms
probsmat=mat.or.vec(n,n)
for(i in 1:noofqueries)
{
 probsmat[queries[i,1],queries[i,2]]=probsmat[queries[i,1],queries[i,2]]+1
}
for(i in 1:(n-1))
{
 for(j in (i+1):n)
 {
  if(probsmat[j,i]!=0 | probsmat[i,j]!=0){probsmat[i,j]=probsmat[i,j]/(probsmat[i,j]+probsmat[j,i])}
 }
}

probs=mat.or.vec(n,2)
probs[1,1]=sum(probsmat[1,2:n])/n
for(i in 2:(n-1))
{
 probs[i,1]=(sum(probsmat[i,(i+1):n])+(i-1)-sum(probsmat[1:(i-1),i]))/n
}
probs[n,1]=(n-sum(probsmat[1:(n-1),n]))/n

### Reformulate to be the same as the previous format

for(i in 1:n)
{
 probs[i,1]=1-probs[i,1]
}

### Sort the probs vector afterwards, and use the usual bandit script! For that, you need the mean:

mu=rep(0,n)
for(i in 1:n){mu[i]=1-probs[i,1]-0.5*probs[i,2]}








###
### Misc stuff
###

N=16
k=kappa(N) # the kappa constant
delta=0.01 # the error tolerance

### Compute threshold sequences beforehand to save time (this is OK, if we don't have more pulls on any arm then 'maxpulls')

maxpulls=5000
kaufmannthr=rep(0,maxpulls)
for(i in 1:maxpulls){kaufmannthr[i]=kaufmann(i,delta)}
ourSGthr=rep(0,maxpulls)
for(i in 1:maxpulls){ourSGthr[i]=ourSGbound(i,delta,k,N)}
KLthr=rep(0,maxpulls)
for(i in 1:maxpulls){KLthr[i]=KLthreshold(i,delta,k)}

tmax=500
muhat=0.1
t=seq(1,tmax,1)
ratio=rep(0,tmax)
for(i in 1:tmax)
{
 # muhat+kaufmannthr[t]
 ratio[i]=(muhat+ourSGthr[t[i]])/computeNEWUCB(N,muhat,KLthr[t[i]],accuracy)
}
plot(ratio,type="l")
abline(h=1,lty=2)



means=rep(0,n)
for(i in 1:n){means[i]=1-probs[i,1]-0.5*probs[i,2]}
plot(means)


D=function(x,y)
{
 return=0
 if(x!=y){return=x*log(x/y)+(1-x)*log((1-x)/(1-y))}
 return
}

grain=0.001
x=seq(grain,1-grain,grain)
xl=length(x)
hoef=rep(0,xl)
kl=rep(0,xl)
kauf=rep(0,xl)

N=8
muhat=0.1
index=muhat/grain
for(i in 1:xl)
{
 hoef[i]=2*(N/(N+1))^2*(x[i]-muhat)^2
 kl[i]=D((N*muhat+x[i])/(N+1),x[i])
 kauf[i]=(4/3)*(x[i]-muhat)^2
}

szar=kl/hoef
plot(x[index:900],szar[index:900],type="l")
szar=kauf/hoef
lines(x[index:900],szar[index:900],col="red")

for(i in 1:xl)
{
 hoef[i]=2*(x[i]-muhat)^2
 kl[i]=D(muhat,x[i])
}
szar=kl/hoef
plot(x[index:900],szar[index:900],col="blue")













mu=rep(0,n)
for(i in 1:n){mu[i]=1-probs[i,1]-0.5*probs[i,2]}
plot(mu)

squareddecay=rep(0,n)
for(i in 1:n){squareddecay[i]=0.5*(1-(i/n)^(1/7))}
lines(squareddecay,col="blue")




###
### Put error bars and arm profile on top of plots
###

xax=seq(0,horizon,howoften)/10000
plot(xax,c(0,KLprob),type="l",col="blue",ylim=c(0,1),main="P(best arm in top 5), N=8, delta=0.01",sub="Contest 517-lil, somewhat-funny & funny merged",xlab="Number of pulls (10 thousands)",ylab="(Empirical) probability - 250 trials",xaxt="n")
lines(xax,c(0,kaufmannprob))
lines(xax,c(0,ourSGprob),col="red")
legend("bottomright",c("Kaufmann lil-UCB","KL-UCB","Our SG lil-UCB"),fill=c('black','blue','red'))
# axis(1,at=seq(0,horizon,howoften*2)/howoften)
axis(1,at=xax)
abline(h=0.9,lty=3)

### Error bars

kaufmannse=rep(0,length(kaufmannprob))
for(i in 1:length(kaufmannse)){kaufmannse[i]=2*sqrt(kaufmannprob[i]*(1-kaufmannprob[i])/numberofreps)}
ourSGse=rep(0,length(ourSGprob))
for(i in 1:length(ourSGse)){ourSGse[i]=2*sqrt(ourSGprob[i]*(1-ourSGprob[i])/numberofreps)}
KLse=rep(0,length(KLprob))
for(i in 1:length(KLse)){KLse[i]=2*sqrt(KLprob[i]*(1-KLprob[i])/numberofreps)}

delta=max(xax)/200
segments(xax[2:length(xax)],(kaufmannprob-kaufmannse),xax[2:length(xax)],(kaufmannprob+kaufmannse))
segments(xax[2:length(xax)]-delta,(kaufmannprob-kaufmannse),xax[2:length(xax)]+delta,(kaufmannprob-kaufmannse))
segments(xax[2:length(xax)]-delta,(kaufmannprob+kaufmannse),xax[2:length(xax)]+delta,(kaufmannprob+kaufmannse))
segments(xax[2:length(xax)],(ourSGprob-ourSGse),xax[2:length(xax)],(ourSGprob+ourSGse),col="red")
segments(xax[2:length(xax)]-delta,(ourSGprob-ourSGse),xax[2:length(xax)]+delta,(ourSGprob-ourSGse),col="red")
segments(xax[2:length(xax)]-delta,(ourSGprob+ourSGse),xax[2:length(xax)]+delta,(ourSGprob+ourSGse),col="red")
segments(xax[2:length(xax)],(KLprob-KLse),xax[2:length(xax)],(KLprob+KLse),col="blue")
segments(xax[2:length(xax)]-delta,(KLprob-KLse),xax[2:length(xax)]+delta,(KLprob-KLse),col="blue")
segments(xax[2:length(xax)]-delta,(KLprob+KLse),xax[2:length(xax)]+delta,(KLprob+KLse),col="blue")

### Adding the mean profile:

height=0.25
buffer=0.05
xax=seq(0,horizon,howoften)/10000
plot(xax,c(0,KLprob),type="l",col="blue",ylim=c(0,1+height+buffer),main="P(best arm in top 5), N=8, delta=0.01",sub="Contest 517-lil, somewhat-funny & funny merged",xlab="Number of pulls (10 thousands)",ylab="(Empirical) probability - 250 trials",xaxt="n")
lines(xax,c(0,kaufmannprob))
lines(xax,c(0,ourSGprob),col="red")
legend("bottomright",c("Kaufmann lil-UCB","KL-UCB","Our SG lil-UCB"),fill=c('black','blue','red'))
# axis(1,at=seq(0,horizon,howoften*2)/howoften)
axis(1,at=xax)
abline(h=0.9,lty=3)

abline(h=1+buffer)
bases=seq(0,max(xax),(max(xax))/(n-1))
meanprofile=(height/max(mu))*mu
bottom=rep(1,n)+buffer
segments(bases,bottom,bases,(bottom+meanprofile))
# points(bases,(bottom+meanprofile),pch=20)




###
### Exporting to pdf
###

### Import the results

setwd('')

forplot=read.csv("forpaper_250reps_cont558-lil.csv")
# forplot=read.csv("forpaper_250reps_cont512-lil.csv")
# forplot=read.csv("forpaper_250reps_alpha12_small.csv")
# forplot=read.csv("forpaper_250reps_alpha1_small_new.csv")

horizon=max(forplot[,1])*10000
howoften=horizon/20

kaufmannprob=forplot[1:20,2]
ourSGprob=forplot[1:20,3]
KLprob=forplot[1:20,4]

forplot1=read.csv("forpaper_50reps_558lil_Thompson_reshaped.csv")
forplot2=read.csv("forpaper_50reps_558lil_Thompson_reshaped.csv")

Tp1=forplot1[1:20,2]
Tp2=forplot2[1:20,2]
Thompsonprob=(Tp1+Tp2)/2


### Import the mean for the mean-profile

data=read.csv("NYr_cont558-lil.csv")
# data=read.csv("NYr_cont512-lil.csv")

### Set up arms.
### WATCH OUT! Data before contest 530 is stored differently!
### After 531, it's "probs[i,j]=(data[i,j+6]/data[i,6])"
### Before 530, it's "probs[i,j]=(data[i,j+1]/(data[i,2]+data[i,3]+data[i,4]))"
### Also, data not necessarily sorted before contest 526

n=length(data[,1]) # Number of arms
probs=mat.or.vec(n,2)
for(i in 1:n)
{
 probs[i,1]=(data[i,7]/(data[i,6]))
# probs[i,2]=0
}

mu=rep(0,n)
for(i in 1:n){mu[i]=1-probs[i,1]-0.5*probs[i,2]}

### Sort ams (if necessary)

sortedprobs=mat.or.vec(n,2)
muhelper=mu
for(i in 1:n)
{
 top=which.is.max(muhelper)
 sortedprobs[i,]=probs[top,]
 muhelper[top]=-1
}
probs=sortedprobs
for(i in 1:n){mu[i]=1-probs[i,1]-0.5*probs[i,2]}
plot(mu)
max(mu)

### Parametric arm setup (for the mean profile)

n=1000
probs=mat.or.vec(n,2)
alpha=1/2
mumax=1
for(i in 1:n){probs[i,1]=(1-mumax)+mumax*((i/n)^alpha)}

mu=rep(0,n)
for(i in 1:n){mu[i]=1-probs[i,1]-0.5*probs[i,2]}
plot(mu)
max(mu)

### Plot!

pdf("558_new_Thompson.pdf",width=6,height=6)
# pdf("512_new_Thompson.pdf",width=6,height=6)
# pdf("alpha12.pdf",width=6,height=6)
# pdf("alpha13.pdf",width=6,height=6)
# pdf("alpha1_small_new.pdf",width=6,height=6)
# pdf("alpha12_small_new.pdf",width=6,height=6)



height=0.25
buffer=0.05
xax=seq(0,horizon,howoften)/10000
plot(xax,c(0,KLprob),type="l",col="blue",ylim=c(0,1+height+buffer),main="P(best arm in top 5), Contest 512",xlab="Number of samples (10 thousands)",ylab="(Empirical) probability - 250 trials",xaxt="n",cex.lab=1.5,cex.main=1.5,cex.axis=1.2)
lines(xax,c(0,kaufmannprob))
lines(xax,c(0,ourSGprob),col="red")
lines(xax,c(0,Thompsonprob),col='purple')
legend("bottomright",c("Kaufmann lil-UCB","KL-UCB","SG lil-UCB",'Thompson'),fill=c('black','blue','red','purple'))
# axis(1,at=seq(0,horizon,howoften*2)/howoften)
axis(1,at=xax)
abline(h=0.9,lty=3)

abline(h=1+buffer)
bases=seq(0,max(xax),(max(xax))/(n-1))
meanprofile=(height*max(mu))*mu
bottom=rep(1,n)+buffer
segments(bases,bottom,bases,(bottom+meanprofile))
points(bases,(bottom+meanprofile),pch=20)
points(bases,(bottom),pch=20)

dev.off()















###
### Why did we see bigger gains before?
###

### The threshold functions (we used union bounds before, so we'll use that here as well)

delta=0.01 # the error tolerance

maxpulls=400000
unionKLthr=rep(0,maxpulls)
for(i in 1:maxpulls){unionKLthr[i]=log((pi^2/6)*i^2/delta)/i}
unionSGthr=rep(0,maxpulls)
for(i in 1:maxpulls){unionSGthr[i]=sqrt(0.5*log((pi^2/6)*i^2/delta)/i)}

### Import data. The contests where both methods were used is 559 (maybe also 558, but I'm not sure what that data is)

# data=read.csv("NYr_cont559-lil.csv")
data=read.csv("NYr_cont559-KL.csv")

n=length(data[,1]) # Number of arms
probs=mat.or.vec(n,2)
for(i in 1:n)
{
 probs[i,1]=data[i,7]/(data[i,6])
 probs[i,2]=data[i,8]/(data[i,6])
}

mu=rep(0,n)
for(i in 1:n){mu[i]=1-probs[i,1]-0.5*probs[i,2]}
plot(mu)
max(mu)

### Experiment parameters

numberofreps=100
accuracy=10^(-6)

m=5
horizon=200000
howoften=horizon/20 # how often do we evaluate whether or not the best arm is in the top m

KLprob=rep(0,(horizon/howoften))
SGprob=rep(0,(horizon/howoften))
isitgood=0

SGavgpulls=rep(0,n)
KLavgpulls=rep(0,n)

for(j in 1:numberofreps)
{

### SG-UCB

 pulls=rep(0,n) # number of pulls of each arm
 empmeans=rep(0,n) # means of arms
 ucb=rep(0,n) # upper conf bounds
 tiebreak=runif(n,0,1) # used to decide how we break ties in the ranking

 for(i in 1:n)
 {
  pulls[i]=pulls[i]+1
  X=threestar(probs[i,1],probs[i,2])
  empmeans[i]=(((pulls[i]-1)/pulls[i])*empmeans[i])+(X/pulls[i])
  ucb[i]=empmeans[i]+unionSGthr[1]
 }

 totalpulls=n

 while(totalpulls<horizon)
 {
  where=which.is.max(ucb)
  pulls[where]=pulls[where]+1
  totalpulls=totalpulls+1

  X=threestar(probs[where,1],probs[where,2])
  empmeans[where]=(((pulls[where]-1)/pulls[where])*empmeans[where])+(X/pulls[where])
  ucb[where]=empmeans[where]+unionSGthr[pulls[where]]

  if((totalpulls%%howoften)==0)
  {
   rankone=1
   for(a in 2:n)
   {
    if(empmeans[1]<empmeans[a]){rankone=rankone+1}
    if(empmeans[1]==empmeans[a]){if(tiebreak[1]<tiebreak[a]){rankone=rankone+1}}
    if(rankone>m){break}
   }
   if(rankone<(m+1)){isitgood=1}else{isitgood=0}
#   if(which.is.max(empmeans)==1){isitgood=1}else{isitgood=0}
   SGprob[(totalpulls/howoften)]=((j-1)/j)*SGprob[(totalpulls/howoften)]+(isitgood/j)
  }

  if(totalpulls==horizon)
  {
   for(i in 1:n){SGavgpulls[i]=((j-1)/j)*SGavgpulls[i]+(pulls[i]/j)}
  }

 }

### KL-UCB

 pulls=rep(0,n) # number of pulls of each arm
 empmeans=rep(0,n) # means of arms
 ucb=rep(0,n) # upper conf bounds
 tiebreak=runif(n,0,1) # used to decide how we break ties in the ranking

 for(i in 1:n)
 {
  pulls[i]=pulls[i]+1
  X=threestar(probs[i,1],probs[i,2])
  empmeans[i]=(((pulls[i]-1)/pulls[i])*empmeans[i])+(X/pulls[i])
  ucb[i]=computeUCB(empmeans[i],unionKLthr[1],accuracy)
 }

 totalpulls=n

 while(totalpulls<horizon)
 {
  where=which.is.max(ucb)
  pulls[where]=pulls[where]+1
  totalpulls=totalpulls+1

  X=threestar(probs[where,1],probs[where,2])
  empmeans[where]=(((pulls[where]-1)/pulls[where])*empmeans[where])+(X/pulls[where])
  ucb[where]=computeUCB(empmeans[where],unionKLthr[pulls[where]],accuracy)

  if((totalpulls%%howoften)==0)
  {
   rankone=1
   for(a in 2:n)
   {
    if(empmeans[1]<empmeans[a]){rankone=rankone+1}
    if(empmeans[1]==empmeans[a]){if(tiebreak[1]<tiebreak[a]){rankone=rankone+1}}
    if(rankone>m){break}
   }
   if(rankone<(m+1)){isitgood=1}else{isitgood=0}
#   if(which.is.max(empmeans)==1){isitgood=1}else{isitgood=0}
   KLprob[(totalpulls/howoften)]=((j-1)/j)*KLprob[(totalpulls/howoften)]+(isitgood/j)
  }

  if(totalpulls==horizon)
  {
   for(i in 1:n){KLavgpulls[i]=((j-1)/j)*KLavgpulls[i]+(pulls[i]/j)}
  }

 }

}

### Plots

### P(best arm in top 5)

xax=seq(0,horizon,howoften)/10000
plot(xax,c(0,KLprob),type="l",ylim=c(0,1),main="P(best arm in top 5)",sub="Contest 559",xlab="Number of pulls (10 thousands)",ylab="(Empirical) probability - 100 trials",xaxt="n")
lines(xax,c(0,SGprob),col="blue")
legend("bottomright",c("KL-UCB","SG-UCB"),fill=c('black','blue'))
# axis(1,at=seq(0,horizon,howoften*2)/howoften)
axis(1,at=xax)
abline(h=0.9,lty=3)


### Ratio of pulls per arm

ratio=KLavgpulls/SGavgpulls
plot(ratio,main="KL_pulls(arm i) / SG_pulls(arm i)",ylim=c(0,max(ratio)))
abline(h=1,lty=3)

### Add empirical pull ratio

KLszar=read.csv("NYr_cont559-KL.csv")
lilszar=read.csv("NYr_cont559-lil.csv")

# Align shit! If the experiment was run on 559-KL:

KLpulls=KLszar[,6]/sum(KLszar[,6])
lilpulls=rep(0,n)
lilsum=sum(lilszar[,6])
for(i in 1:n)
{
 j=1
 while(KLszar[i,3]!=lilszar[j,3]){j=j+1}
 lilpulls[i]=lilszar[j,6]/lilsum
}

# Align shit! If the experiment was run on 559-lil:

lilpulls=lilszar[,6]/sum(lilszar[,6])
KLpulls=rep(0,n)
KLsum=sum(KLszar[,6])
for(i in 1:n)
{
 j=1
 while(lilszar[i,3]!=KLszar[j,3]){j=j+1}
 KLpulls[i]=KLszar[j,6]/KLsum
}

# The plot itself

ratio=KLavgpulls/SGavgpulls
oldratio=KLpulls/lilpulls
plot(ratio,main="KL_pulls(arm i) / SG_pulls(arm i)",ylim=c(0,max(c(ratio,oldratio))))
lines(oldratio,col="blue")


### Export results

### P(best in top 5)

results=matrix(c(xax[2:length(xax)],SGprob,KLprob),ncol=3,nrow=length(KLprob))
write.table(results,"top5_100reps_unionb_cont559-KL.csv",row.names=F,col.names=c("pulls","SG","KL"),sep=",")

### Pulls

results=matrix(c(SGavgpulls,KLavgpulls),ncol=2,nrow=n)
write.table(results,"pulls_100reps_unionb_cont559-KL.csv",row.names=F,col.names=c("SG","KL"),sep=",")









