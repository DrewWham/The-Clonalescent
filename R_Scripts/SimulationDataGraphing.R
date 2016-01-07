
x <- scan("infile.test.txt", what="", sep="\n")
y <- strsplit(x, "[[:space:]]+")
names(y) <- sapply(y, `[[`, 1)
y <- lapply(y, `[`, -1)
lapply(y,as.numeric)->y

unlist(lapply(y,sum))->n
unlist(lapply(y,length))->S
as.numeric(sapply(strsplit(names(y),"[.]"),'[',1))->Ne
as.numeric(sapply(strsplit(names(y),"[.]"),'[',2))->Ng

ptm <- proc.time()

optW<-function(x){optimize(f=theta.likelihood,interval=c(0,100000),maximum=TRUE,J=sum(x),S=length(x))$maximum}
proc.time() - ptm
unlist(lapply(y,optW))->Psi

unlist(lapply(y,ExPWD))->Fo

cbind(Ne,Ng,n,S,Psi,Fo)->Data.tab
row.names(Data.tab)<-seq(1:dim(Data.tab)[1])
data.frame(Data.tab)->Data.tab
optWexp<-function(x){optimize(f=theta.likelihood,interval=c(0,10000),maximum=TRUE,J=x[1],S=x[2])$maximum}
apply(Data.tab[,c(1,2)],1,optWexp)->Data.tab$PsiExp


Data.tab$PsiExp-Data.tab$Psi->Data.tab$Psierror
abs(Data.tab$Psierror/Data.tab$PsiExp)->Data.tab$PsiRel

Data.tab$S/Data.tab$n->Data.tab$GtoN
 Data.tab$Ng/Data.tab$Ne->Data.tab$GtoNexp
 Data.tab$GtoNexp-Data.tab$GtoN->Data.tab$GtoNerror
abs(Data.tab$GtoNerror/Data.tab$GtoNexp)->Data.tab$GtoNrel

melt(acast(Data.tab,Ne~Ng,mean,value.var="Fo"),value.name="Fe")->Fe.tab
names(Fe.tab)[1:2]<-c("Ne","Ng")
merge(Data.tab,Fe.tab)->Data.tab
(Data.tab$Fo-Data.tab$Fe)/Data.tab$Fe->Data.tab$Dpsi

ENg.ap<-function(x){
	x[1]->W
	x[2]->Ne
		Ek<-NULL
		round(Ne)->Ne
		0:Ne->v
		unlist(lapply(v,FUN=function(x,W){W/(W+x)},W=W))->Ek
		sum(Ek)->Ng
Ng
	}


apply(Data.tab[,c(5,1)],1,ENg.ap)->Data.tab$ENg
