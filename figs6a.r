
#
library(datawizard)

d<-tdata[,-1]
d1<-data.frame(d)
d2<-normalize(d1)

dm<-dist(d2, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)
library(vegan)
library(ape)
library(ggplot2)
library(grid)
library(dplyr)
library(multcomp)
library(patchwork)
library(ggm)
library(ggpubr)

library(ggsci)

############bray_curtis_distance_matrix
data<- read.table("anno.tsv", header = TRUE, sep = "\t" )
data$group<-data$subject
data$sample<-data$col
rownames(data)<-data$col

# dis <- read.table("bray_curtis_distance_matrix.tsv",sep = "\t",header = T)

dis <- as.matrix(dm)
colnames(dis)<-data$col
rownames(dis)<-data$col
cid<-intersect(rownames(dis),rownames(data))

# dis<-dis[cid,]
# dis<-dis[,paste('X',cid,sep='')]

data<-data[cid,]
# rownames(dis)<-data$subid
# colnames(dis)<-data$subid

pcoa<- pcoa(dis, correction = "none", rn = NULL)
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
plotdata <- data.frame(rownames(dis),PC1,PC2,data$group)
colnames(plotdata) <-c("sample","PC1","PC2","group")

groups<-data.frame(data$sample,data$group)
colnames(groups) <-c("V1","V2")
groups <- as.list(groups)


pich=c(21:24)
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
Palette <- c("#000000", "#000000", "#000000", "#000000")
pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)

otu.adonis=adonis(dis~V2,data = groups,distance = "bray")

scale_fill_nejm(alpha=0.3)+
plotdata$group<-factor(plotdata$group,levels=c('WT-AOM/DSS','KO-AOM/DSS'))
p<-ggplot(plotdata, aes(PC1, PC2)) +
  geom_point(aes(colour=group,fill=group),size=4,alpha=0.8)+
  geom_text(aes(x = 0.05,y = 0.35,label = paste("PERMANOVA: p-value = ",otu.adonis$aov.tab$`Pr(>F)`[1],sep = "")),size = 5,hjust = 0)+
  stat_ellipse(aes(x = PC1, y = PC2, fill = group), geom = "polygon", alpha = 0.3, level = 0.90)+ 
  scale_color_nejm(alpha=0.3)+
  labs(title="PCoA - The composition of gut microbiome") + 
  xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  theme(panel.grid =element_blank(),panel.background=element_rect(fill="transparent",colour="black"))+
	theme(legend.title=element_blank(),legend.background=element_blank())+
	theme(plot.title = element_text(hjust = 0.5),plot.subtitle=element_text(hjust = 0.5))


ggsave(p, filename='bray_curtis_distance_matrix.pdf', width=30, height=30, units=c("cm"))


