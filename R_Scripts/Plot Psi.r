library(ggplot)
Psi<-read.csv("Psi.mcmc2",header=FALSE,sep="\t")
names(Psi)<-c("Psi","ENg","Site")
Psi.p<-ggplot( Psi, aes(x=Site,y=Psi,fill=Site))+ geom_violin(position="identity",alpha=0.7)+coord_flip()
ggsave("Psi.pdf",Psi.p)

ENg.p<-ggplot( Psi, aes(x=Site,y=ENg,fill=Site))+ geom_violin(position="identity",alpha=0.7)+ylim(c(0,100))+coord_flip()
ggsave("ENg.pdf",ENg.p,)
