GFS.ID<-"Test2"
Psi(GFS)->out
Ne<-59
Dw<- -10
dw<-0
ifelse(Dw< -2,1,dw)->dw
ifelse(Dw>-2 & Dw< -0.5,2,dw)->dw
ifelse(Dw>-0.5 & Dw<0.5,3,dw)->dw
ifelse(Dw>0.5 & Dw<2,4,dw)->dw
ifelse(Dw> 2,5,dw)->dw
brewer.pal(5,"RdBu")[dw]->cl
unlist(lapply(out,ENg,Ne=Ne))->ndist

round(ndist)->ndist
round(cumsum(rep(1,length(ndist)))/length(ndist),2)*100->color
#cumsum(ndist[order(ndist)])/sum(ndist)->color
ndist[order(ndist)]->val
cbind(val,color)->tab
data.frame(tab)->tab
names(tab)<-c("ID","Value")
tab[!duplicated(tab$ID),]->tab
tab <- as.data.table(tab)
setkey(tab,ID,Value)
tab[CJ(seq(1:Ne)),roll=TRUE]->tab
data.frame(tab)->tab
tab$Value[is.na(tab$Value)]<-0
tab->cdf.tab
as.factor(round(as.numeric(cdf.tab$Value)))->cdf.tab$CI
c(rep(1:45,each=2),seq(50,100,by=5))->graph.pro
colorRampPalette(c(cl,cl,"gray30"))(100)->colors
colors[graph.pro[as.numeric(levels(cdf.tab$CI))+1]]->col

 
 
 NvG.plot<-ggplot(cdf.tab,aes(x=ID,fill=CI))+geom_bar(binwidth=Ne)+xlim(c(0,Ne))+ coord_polar(theta = "x")+scale_fill_manual(values=col)+geom_text(x=0,y=0,label=paste("Ng=",round(mean(ndist))))+geom_text(x=0,y=Ne,label=paste("Ne=",Ne))
 ggsave(paste(GFS.ID,".pdf"))