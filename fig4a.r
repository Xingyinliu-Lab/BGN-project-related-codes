

library(forcats)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyverse)
library(hrbrthemes)
library(psych)
library(pheatmap)
library(psych)
library(igraph)
library(ggm)
library(ggpubr)

tdata<- read.table("faithpd.csv", header = TRUE, sep = ",")
tdata$group<-factor(tdata$group,levels=c('WT-AOM/DSS','KO-AOM/DSS'))

p2 <- ggplot(data = tdata,aes(x =group,y = value,fill =group,color=group)) + 
	geom_violin(alpha=1, outlier.size=0, fill="transparent") +
	scale_color_nejm()+
	geom_point(size = 3, alpha = 0.8,position = position_jitter(0.1))+
	labs(title="",x="", y="Faith pd index")+
	#theme(plot.title=element_text(hjust=0.5),axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) +
	theme(panel.grid =element_blank(),panel.background=element_rect(fill="transparent",colour="black"))+
	theme(legend.position="none")+stat_compare_means(method = "wilcox.test",paired=FALSE)


ggsave(p2, device=cairo_pdf, filename='faithpd_index.pdf', width=25, height=15, units=c("cm"))

wilcox.test(value~group,data=tdata,paired=FALSE,alternative ="two.sided")

        Wilcoxon rank sum test

data:  value by group
W = 57, p-value = 0.7343
alternative hypothesis: true location shift is not equal to 0
