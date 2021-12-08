
#multi-trait Manhattan plot for root microbiome beta diversity
library(dplyr)
library(data.table)

a = list.files(pattern = ".assoc.txt")                                     
root.multi<-fread(a[1],header=T) 
root.multi<-a.multi[,c(2,1,3,12)]
ss<-basename(a[1])
ss<-gsub('_GE_GWAS_root.assoc.txt','',ss)
names(root.multi)[4]<-ss

for (i in 2:length(a)){
        new.data = fread(a[i],header=T)
        new.data<-new.data[,12]
        ss<-basename(a[i])
        ss<-gsub('_GE_GWAS_root.assoc.txt','',ss)
        names(new.data)<-ss
        a.multi = cbind(root.multi,new.data)
}

#multi trait ractangle
MVP.Report(root.multi, plot.type="m", multracks=T, 
           threshold=1.07e-8,threshold.lty=2, 
           threshold.lwd=1, threshold.col="grey", 
           amplify=TRUE,bin.size=1e6,
           chr.den.col=NULL, 
           col=c("#BEBEBE","#DAA520","#0072B2","#009E73","#CC79A7","#00A4DE","#2E3092"),
           file.type="tiff",memo="",dpi=600)

#multi trait qq plot
MVP.Report(root.multi, plot.type="q",
           col=c("#BEBEBE","#DAA520","#0072B2","#009E73","#CC79A7","#00A4DE","#2E3092"),
           conf.int=F,box=FALSE,multracks= T,file.type="tiff",memo="",dpi=600)

#Single trait-Manhattan
root.p3<-fread("root_Axis.3.txt.assoc.txt",header = T)

tiff(file="root_pco3_Manh.tiff",compression="lzw",units="in",res=600,height=6,width=10)

MVP.Report(root.p3, plot.type="m", col=c("#DAA520","#009e73"), LOG10=TRUE, ylim=c(0,20),
           threshold=c(1.07e-8), threshold.lty=2, threshold.lwd=1,
           threshold.col="grey50", amplify=F,chr.den.col=NULL,
           signal.cex=c(1,1),signal.pch=c(19,19),
           file.type="TIFF",memo="",dpi=600)
dev.off()

####################
#Single trait qq plot
tiff(file="root_pco3_QQ.tiff",compression="lzw",units="in",res=600,height=5,width=5)

MVP.Report(root.p3,plot.type="q",conf.int=F,box=TRUE,file.type="TIFF",memo="",dpi=600)
dev.off()


#####################################
#Local positions of significant SNP associated with Root PCO3
root.p3<-fread("root_Axis.3.txt.assoc.txt",header = T)
peak1<-filter(root.p3,chr==11&ps>=21844005 &  ps<=21944005)
peak2<-filter(root.p3,chr==11&BP>=9738437 & BP<=9838437 )
peak3<-filter(root.p3,chr==3 & BP>=689512 &  BP<=789512)
#plot

pdf(file="peak1.pdf", width=5,height=5,useDingbats = F)
ggplot(peak1) +
        geom_point(aes(x=BP/1000,y=-log10(P)),size=2,color=pal1[4])+
        geom_abline(intercept = -log10(1.07e-8),slope = 0,lty=2)+
        xlab("Chromosome 11 (kb)")+
        theme_bw()+theme(legend.position = c("top"))+
        theme(panel.background = element_blank(),panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "black"),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"))+
        ylim(0,10)+
        scale_x_continuous(expand = c(0,0),position = "bottom")

dev.off()
