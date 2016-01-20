library(untb)
library(plyr)
library(ggplot2)
library(MCMCpack)
library(zoo)
source("coalescent.r")

#likelihood function for psi or theta
Tha<-function(Tha,n,S){
  if (Tha==0) {-10000000000->jllh}
  else if (Tha<0) {-100000000000->jllh}
	else if (Tha>100000) {-100000000000->jllh}
  else {
    
    theta.likelihood(Tha,J=n,S=S)->jllh1
 
  
jllh1}}

#gets max likelihood for psi or theta
opt<-function(n,S){optimize(f=Tha,interval=c(0,100000),maximum=TRUE,n=n,S=S)$maximum}


#function that calculates the expected pairwise differences from a vector genotype frequency spectrum
ExPWD<-function(GFS){

sum(GFS)->n
length(GFS)->S
(n*(n-1))/2->tpw
apw<-NULL
for (i in 1:S){
	GFS[i]->ng
	(ng*(ng-1))/2->spw
	c(apw,spw)->apw
	}
	sum(apw)/tpw->pi
	pi
}


#simulates samples from populations with the same psi as observed data to make the null distribution/confidence interval around Dpsi; input is a single Genotype frequency spectrum
Dpsi<-function(GFS,r){
	sum(GFS)->n
length(GFS)->S
opt(n,S)->th
Ndist<-NULL
for (i in 1:r){
s1 <- simulate.coalescent(sample=n, theta=th,plot.tree=FALSE)
as.matrix(s1$data$seq)->a
t(a)->x
x.unique <- unique(x, MARGIN  = 2)

freq <- apply(x.unique, MARGIN = 2, 
              function(b) sum(apply(x, MARGIN = 2, function(a) all(a == b)))
)
paste(rep("S",length(freq)),1:length(freq),sep="")->Names
rep(Names,freq)->obs
as.count(obs)->cnt
ExPWD(cnt)->sim.pi
c(sim.pi,Ndist)->Ndist
}
ExPWD(GFS)->Opi
mean(Ndist)->Exdist
(Opi-Ndist)/Exdist->Dpsi.out
data.frame(Dpsi.out)->Dpsi.out
	Dpsi.out
}


#produces the posterior probability distribution of Psi given a GFS; input is a single Genotype frequency spectrum
Psi<-function(GFS){
sum(GFS)->n
length(GFS)->S
opt(n,S)->init
post.samp.sexe <- MCMCmetrop1R(Tha, theta.init=c(init),n=n,S=S,thin=100, mcmc=100000, burnin=500, V=NULL,tune=c(1), verbose=500, logfun=TRUE)
SEXe<-(post.samp.sexe)
SEXe
}


#simulates samples from populations with the same psi as observed data to make the null distribution/confidence interval around Dpsi; input is a list of genotype frequency spectrums with names of list elements being used as groups
Dpsi.lst<-function(GFS.lst,r){
	Dw.all<-NULL
for (i in 1:length(GFS.lst)){
	GFS.lst[[i]]->GFS
names(GFS.lst)[i]->GFS.ID
	Dpsi(GFS,r)->Dw
	cbind((Dw),rep(GFS.ID,length(Dw)))->Dw.tab
rbind(Dw.all,Dw.tab)->Dw.all
}
data.frame(Dw.all)->Dw.all
names(Dw.all)<-c("Dw","GFS.ID")
as.numeric(levels(Dw.all$Dw))[Dw.all$Dw]->Dw.all$Dw
Dw.all
	
}

#produces the posterior probability distribution of Psi given a GFS; input is a list of genotype frequency spectrums with names of list elements being used as groups
Psi.lst<-function(GFS.lst){
	W.all<-NULL
for (i in 1:length(GFS.lst)){
	GFS.lst[[i]]->GFS
names(GFS.lst)[i]->GFS.ID
	Psi(GFS)->W
	cbind((W),rep(GFS.ID,dim(W)[1]))->W.tab
rbind(W.all,W.tab)->W.all
}
data.frame(W.all)->W.all
names(W.all)<-c("W","GFS.ID")
as.numeric(levels(W.all$W))[W.all$W]->W.all$W
W.all
}
#function that takes in Se, Ne and Pi used to solve for Se in the optimize function
Pibd<-function(Se,Ne,pi){
		((1-Se)^2)/((Ne*Se*(2-Se))+(1-Se)^2)->pibd
		-abs(pi-pibd)->out
		out
		}
# this solves the Pibd equation for the value of Se		
NSex<-function(pi,Ne){
1-((Ne*pi)/(sqrt(Ne*pi*(pi*(Ne-1)+1))))->ESe
ESe
}
# gets Psi for the GFS then simulates a bunch of pi values for that value of Psi, then uses those values to make a null distribution for Se given a value for Ne, r is the number of simulations you want to run		
ESex<-function(GFS,Ne,r){
	sum(GFS)->n
length(GFS)->S
opt(n,S)->th
Ndist<-NULL
for (i in 1:r){
s1 <- simulate.coalescent(sample=n, theta=th,plot.tree=FALSE)
as.matrix(s1$data$seq)->a
t(a)->x
x.unique <- unique(x, MARGIN  = 2)

freq <- apply(x.unique, MARGIN = 2, 
              function(b) sum(apply(x, MARGIN = 2, function(a) all(a == b)))
)
paste(rep("S",length(freq)),1:length(freq),sep="")->Names
rep(Names,freq)->obs
as.count(obs)->cnt
ExPWD(cnt)->sim.pi
c(sim.pi,Ndist)->Ndist
}


sapply(Ndist,NSex,Ne=Ne)->Sex.rate
data.frame(Sex.rate)->Sex.rate
}
		
#does the same as above but takes a list of GFS and a list of Ne of the same length as input 	
ESex.lst<-function(GFS.lst,Ne.lst,r){
	ESe.all<-NULL
for (i in 1:length(GFS.lst)){
	GFS.lst[[i]]->GFS
names(GFS.lst)[i]->GFS.ID
	ESex(GFS,Ne.lst[i],r=r)->ESe
	cbind((ESe),rep(GFS.ID,dim(ESe)[1]))->ESe.tab
rbind(ESe.all,ESe.tab)->ESe.all
}
data.frame(ESe.all)->ESe.all
names(ESe.all)<-c("ESe","GFS.ID")
ESe.all
}
