  #load dependencies
library(phyloseq)
library(DESeq2)
library(reshape2)
library(dplyr)
library(ape)
library(vegan)
library(agricolae)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(scales)

pal2<-c("#BEBEBE","#DAA520","#0072B2","#3C6950","#B83F2F","#FFF148","#2E3092")
show_col(pal2)

mypal =pal_d3("category20")(20)
show_col(mypal)

###########################################################
#load data

#taxonomy 
taxa.inf <- read.table("otu.sintax", sep = "\t")
#otu table
bac.otu <- read.csv("otutab.csv",sep = "\t",
                       row.names = 1,
                       header = TRUE)
#meta data
env.df<-read.csv("env.csv")
# import phylogenetic tree
bac.tree <- read.tree("clusters.tree")

#############################################################################################################################


# make phyloseq object for all compartments

bac.phylo <- phyloseq(otu_table(bac.otu,
                                taxa_are_rows = TRUE),
                      tax_table(taxa),
                      sample_data(env.df),
                      phy_tree(bac.tree))


bac.deseq <- phyloseq_to_deseq2(physeq = bac.phylo,
                                design = ~position)
bac.deseq.wald <- DESeq(bac.deseq,
                        fitType = "parametric",
                        test = "Wald")

bac.norm <- phyloseq(otu_table(counts(bac.deseq.wald, 
                                      normalized = TRUE),
                               taxa_are_rows = TRUE),
                     tax_table(taxa1),
                     sample_data(env.df),
                     phy_tree(bac.tree))


bac.otu.norm<-data.frame(otu_table(bac.norm))

#############################################################################################################################
# make phyloseq object for each compartment

library(phyloseq)
rhizo.phylo <- phyloseq(otu_table(rhizo.otu,
                                 taxa_are_rows = TRUE),
                       tax_table(taxa),
                       sample_data(rhizo.env),
                       phy_tree(bac.tree))

library(DESeq2)

rhizo.deseq <- phyloseq_to_deseq2(physeq = rhizo.phylo,
                                  design = ~cultivar)
rhizo.deseq.wald <- DESeq(rhizo.deseq,
                          fitType = "parametric",
                          test = "Wald")

rhizo.norm <- phyloseq(otu_table(counts(rhizo.deseq.wald, 
                                        normalized = TRUE),
                                 taxa_are_rows = TRUE),
                       tax_table(taxa),
                       sample_data(rhizo.env),
                       phy_tree(bac.tree))

rhizo.otu.norm<-data.frame(otu_table(rhizo.norm))

#output: rhizo.norm, bulk.norm,root.norm

##################################################
#microbial dissimilarity
rhizo.dist<-distance(rhizo.norm, method="wunifrac", type="samples")
root.dist<-distance(root.norm, method="wunifrac", type="samples")

#############################################################################################################################

#microbiome compositions

source("make_taxa_bar.R")


rhizo.phylum<- make_bar(phyloseq = rhizo.norm, 
                       keepnum = 7,
                       level = 2,
                       sample.class=1:315
)

rhizo.bar<-data.frame(env.df,rhizo.phylum)
rhizo.bar<-data.frame(aggregate(.~position,data =rhizo.bar, mean))
rhizo.bar2 <- data.frame(melt(rhizo.bar))


phylum.plot<-ggplot(bar.df2,aes(x=group, y=value*100,fill=variable))+
  geom_bar(stat = "identity",width = 1)+
  xlab("Cultivar")+
  ylab("Relative abundance (%)")

phylum.plot+scale_fill_manual(values = c(mypal,"light grey"),name="")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank(),panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(hjust = 0.5,vjust = 1,color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_y_continuous(expand = c(0,0))+
  ylab("Relative abundance (%)")
##############################################################

#Cophylogenetic trees

a<-aggregate(.~population,data=rhizo.bar,mean)
rownames(a)<-a[,1]
a<-a[,-1:-6]
a.tree<-bionj(rhizo.dist)

b<-aggregate(.~population,data=root.bar,mean)
rownames(b)<-b[,1]
b<-b[,-1:-6]
b.tree<-bionj(root.dist)

library(phytools)
tr.assoc<-data.frame(Rhizo=rownames(a),Root=rownames(b))
rice.cophylo<-cophylo(a.tree,b.tree,assoc=tr.assoc)

par(mar=c(5.1,4.1,4.1,2.1))
plot(rice.cophylo,ftype="reg",link.type="straight",link.lwd=2,
     link.lty="dashed",link.col=make.transparent(mypal[8],
                                                 0.5))
######################################################################################################################
#phylum ANOVA between populations
library(agricolae)
for (i in 7:13) {
  s<-aov(root.bar[,i]~population,data=root.bar)
  result_1 <- HSD.test(s, "population", group = T)
  print(colnames(root.bar)[i])
  print(summary(s))
  print(result_1)
}

library(reshape2)
a<-aggregate(.~population,data=rhizo.bar,mean)
a1<-melt(a)

#######################################################
#bubble plot
pops<-c("ADMIX","AROMATIC","AUS","IND","TEJ","TRJ")
bubble.plot<-ggplot(a1,aes(y=population,
                           x=variable,size=value*100))+
  geom_point(shape=16,color=mypal[13])+
  theme_bw()+
  ylab("Rice Population")+
  xlab("Phylum")+
  guides(color=F)+
  scale_size_continuous(range = c(5,15),
                        limits = c(0,100),
                        name = "Relative Abundance (%)")+
  scale_y_discrete(guide = guide_axis(position = "left"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",angle = 45, hjust = 1,vjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks.length = unit(0.05, "cm"),
        legend.position = "none")
#############################################################################################################################

#PCoA

bac.pcoa <- ordinate(bac.norm,method = c("PCoA"), distance = "wunifrac")
ord.df <- data.frame(bac.pcoa$vectors[,1:3],sample_names(bac.phylo),sample_data(bac.phylo))
result.bray <-bac.pcoa$values[,"Relative_eig"]
PCo1 = as.numeric(sprintf("%.4f",result.bray[1]))*100
PCo2 = as.numeric(sprintf("%.4f",result.bray[2]))*100
PCo3 = as.numeric(sprintf("%.4f",result.bray[3]))*100
xlab=paste("PCOA1(",PCo1,"%)",sep="") 
ylab=paste("PCOA2(",PCo2,"%)",sep="")
zlab=paste("PCOA3(",PCo3,"%)",sep="")


pca<-ggplot(ord.df, aes(x = Axis.1,
                        y = Axis.2,
                        color =position))+
  geom_point(size = 2,alpha=0.6)+
  theme_bw()+  xlab(xlab)+
  ylab(ylab)+  theme(legend.position = c("top"))

pca+
  scale_y_continuous(expand = c(0,0.005))+scale_color_manual(values = pal1[c(7,2,4)])+
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  theme(legend.position = "top")
#############################################################################################################################

#ANOSIM

unlist(distanceMethodList)
total.dist<-distance(bac.norm, method="wunifrac", type="samples")
total.anosim<-with(env.df, anosim(total.dist, factor(position)))
summary(total.anosim)

#anosim plot
anosim.df<-data.frame(rank=total.anosim$dis.rank,position=total.anosim$class.vec)

p1<-ggplot(anosim.df, aes(x=position, y=rank)) +
  geom_boxplot(notch = F,alpha=0.8,fill=mypal[11])+
  xlab("Populations")+ylab("Weighted-UniFrac Rank")
p1+theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank(),panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour = "black",size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8))+ 
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0.1,0))


#####################################################################################################################
#permanova
library(vegan)
rice.adonis<-adonis(bac.phylum~population*position,data = env.df,
                     permutations = 9999,method = "bray",by=NULL)

otuput <- data.frame(rice.adonis$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')

############################################################################################################################
#mantel test
mantel(rice.dist,root.dist,method = 'spearman', permutations = 9999, na.rm = TRUE)
mantel(rice.dist,rhizo.dist,method = 'spearman', permutations = 9999, na.rm = TRUE)
mantel(root.dist,rhizo.dist,method = 'spearman', permutations = 9999, na.rm = TRUE)

#plot
micro.m = melt(as.matrix(rhizo.dist))
micro.m = micro.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)
names(micro.m)<-c("V1","V2","value1")

phylo.m<-melt(as.matrix(rice.dist))
phylo.m = phylo.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)

scatter.m<-data.frame(micro.m,phylo.m)
scatter.m<-scatter.m[,-4:-5]
scatter.m<-scatter.m[,-1]

scatter.m<-data.frame(aggregate(.~V2,data = scatter.m,mean),population=rice99.env$Sub.population.1)
names(scatter.m)<-c("rice","microbial","genetic","population")

p1 = ggplot(scatter.m, aes(x = as.numeric(microbial), y = as.numeric(genetic),color=population)) +  
  geom_point(size=2,alpha=0.5)+
  scale_color_manual(values = pal1)+
  #geom_smooth(method=lm, se=T,color=mypal[3])+
  xlab("Microbial wUnifrac Dissimilarity")+ylab("Rice Genetic Dissimilarity")

p1+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank(),panel.border = element_blank(),
        axis.line = element_line(color = "black"),axis.text = element_text(size=8,color = "black"),
        axis.title = element_text(size = 8),
        axis.ticks.length = unit(0.05, "cm"),
        legend.position = "right")
##################################################################################################################

#alpha diversity

#shannon
shannon.div <- plot_richness(physeq = bac.norm,
                             x = "position",
                             measures =  c("Shannon"))
div.df <- shannon.div$data
#one-way anova
shannon.aov<-aov(value~position,data = div.df)
summary(shannon.aov)
TukeyHSD(shannon.aov)

result_1 <- HSD.test(shannon.aov, "position", group = T)
print(result_1)
# plot

p<-ggplot(div.df, aes(x=factor(position,levels = c("BulkSoil","Rhizosphere","Root")), y=value)) +
  geom_boxplot(aes(colour=position))+
  geom_jitter(size=1.5,shape=16, position=position_jitter(0.1),aes(colour=position),alpha=0.5)+
  xlab("Compartment")+ylab("Shannon Index")

p+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank(),panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour = "black",size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8))+ 
  scale_color_manual(values = pal1[c(7,2,4)])+
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0.1,0.01))

ggsave("shannon_position.pdf",height = 4,width = 4,useDingbats=F)






#############################################################################################################################
#canonical correlation
require(ggplot2)
require(GGally)
require(CCA)
#pca correlation

p1<-rice99.pca[,3:12]

#canonical correlation
cc1 <- cc(p1,rhizo.pco1_data)
# display the canonical correlations
cc1$cor
# raw canonical coefficients
cc1[3:4]
# compute canonical loadings
cc2 <- comput(p1,rhizo.pco1_data, cc1)
# display canonical loadings
cc2[3:6]
library(CCP)
# tests of canonical dimensions
rho <- cc1$cor
## Define number of observations, number of variables in first set, and number of variables in the second set.
n <- dim(p1)[1]
p <- length(p1)
q <- length(rhizo.pco1_data)
## Calculate p-values using the F-approximations of different test statistics:
p.asym(rho, n, p, q, tstat = "Wilks")
p.asym(rho, n, p, q, tstat = "Hotelling")
p.asym(rho, n, p, q, tstat = "Pillai")
p.asym(rho, n, p, q, tstat = "Roy")
#Canonical correlation analysis (CCA) measures the degree of linear relationship between two sets of variables. 
#The number of correlation coefficients calculated in CCA is equal to the number of variables in the smaller set:
#m = min(p,q). The coefficients are arranged in descending order of magnitude: rho[1] > rho[2] > rho[3] > ... > rho[m].
#Except for tstat = "Roy", the function p.asym calculates m p-values for each test statistic: the first p-value is
#calculated including all canonical correlation coefficients, the second p-value is calculated by excluding rho[1], 
#the third p-value is calculated by excluding rho[1] and rho[2] etc., therewith allowing assessment of the statistical
#significance of each individual correlation coefficient. On principle, Roy's Largest root takes only rho[1] into 
#account, hence one p-value is calculated only.

