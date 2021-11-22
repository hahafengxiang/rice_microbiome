#1. load dependencies
library(ape)
library(ggplot2)
library(dplyr)
library(ggsci)
library(scales)
############################################################################
#2. input data
#IBS dissimilarity matrix of the SNPs was calculated in PLINK software
rice_pca.dist<- read.table("rice/plink.mdist", quote="\"", comment.char="")
#
rice_names <- read.delim("rice/plink.mdist.id", header=FALSE)
column<-rice_names$V1
row<-rice_names$V2
dimnames(rice_pca.dist)=list(row,column)
rice_pca.dist<-as.matrix(rice_pca.dist)
rice.dist<-dist(rice_pca.dist)

##############################################################################################
#principle component analysis (PCA) for the genetic variation of rice varieties was done
# in PLINK software

#3. load PCA results to R

rice.pca <- read.table("rice/riceQC.eigenvec", quote="\"", comment.char="")
pca.eig <- read.table("rice/riceQC.eigenval", quote="\"", comment.char="")
pca.eig<-mutate(pca.eig,V2=round(V1*100/sum(V1),2))
rice.pca1<-mutate(rice.pca,Group=env.df$population)

xlab=paste("PC1(",pca.eig[1,2],"%)",sep="") 
ylab=paste("PC2(",pca.eig[2,2],"%)",sep="")
zlab=paste("PC3(",pca.eig[3,2],"%)",sep="")


pca<-ggplot(rice.pca1, aes(x =PC1,
                        y = PC2,color=factor(Group)))+
  geom_point(size=2,alpha=0.5)+
  scale_color_manual(values =pal2,name="Population")+
  theme_bw()+  xlab(xlab)+
  ylab(ylab)

pca+
  theme(legend.position  = "right",panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(colour = "black"))+ 
  scale_y_continuous(expand = c(0.05,0.05))+ 
  coord_fixed(ratio=8/16)

######################################################
#  4. Visualization of rice Population structure
#input the ancestry results from ADMIXTURE software
ta1<-read.table("rice_admixture.6.Q")
head(ta1)

bar.plot<-ggplot(ta1,aes(x=Cultivar,y=value*100,fill=factor(Ancestry)))+
  geom_bar(stat = "identity",width = 1)  +    
  scale_fill_manual(values = c("#1B8E77","#7570B3","#E6AB02"))+
  theme_minimal()+
  xlab("Cultivar")+
  ylab("Ancestry")

bar.plot+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background = element_blank(),panel.border = element_blank(),legend.position = "none")+
  theme(axis.line.x = element_blank(),axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank())+
  scale_y_continuous(expand = c(0,0)) + coord_fixed(ratio=1/4)

