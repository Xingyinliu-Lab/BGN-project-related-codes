library(ggplot2)
library(ggrepel)

resdata =na.omit(resdata)
resdata$significant <- as.factor(resdata$padj<0.05 & abs(resdata$log2FoldChange)>1)
resdata$threshold = factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) >= 1, ifelse(resdata$log2FoldChange>=1,'Up','Down'),'NS'),levels=c('Up','Down','NS'))

p<- ggplot(data=resdata,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
geom_point()+
# ylim(0,8)+
scale_color_manual(values= c("#DC143C","#00008B","#808080"))+
geom_hline(yintercept=-log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
# # 添加竖线padj<0.05
geom_vline(xintercept=c(1,-1),lty=4,lwd=0.6,alpha=0.8)+ 
# # 添加横线|FoldChange|>2
theme_bw()+
labs(title="Volcanoplot",xlab="log2 (fold change)",ylab="-log10 (padj)")+ 
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),plot.title = element_text(hjust = 0.5))+
geom_text_repel(data=subset(resdata, padj < 0.001 & abs(log2FoldChange)> 1.5),aes(label=Gene),col="black",alpha=0.8)
p

p2= p + scale_y_sqrt()
p2
# 
pdf(paste0(Con1,"-",Con2,"-volcanoplot.pdf"),width=12,height=8)
print(p2)
dev.off()


