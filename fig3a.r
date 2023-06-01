difdata <- merge(as.data.frame(diff_gene_Group2),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
names(difdata)[1] <- "Gene"

difdata =na.omit(difdata)

Difdata1 <- difdata[which(difdata$padj < 0.001 & abs(difdata$log2FoldChange)> 1.5),]

library(pheatmap)
library(RColorBrewer)
#RColorBrewer包就提供了很多种适合不同场景进行可视化的调色版：
RColorBrewer::display.brewer.all()
### manu
rdylbu <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))
rdbu <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))
navy <- colorRampPalette(c("navy", "white", "firebrick3"))
jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026'))
warmcold <- colorRampPalette(c(rev(cold(21)), warm(20)))

a_mean<- data.frame("KO_LPS"=rowSums(Difdata1[,c("KO_LPS1","KO_LPS2","KO_LPS3")])/3,"WT_LPS"=rowSums(Difdata1[,c("WT_LPS1","WT_LPS2","WT_LPS3")])/3)
row.names(a_mean)=Difdata1$Gene

b_mean=a_mean
b_mean$dev=abs(a_mean$KO_LPS-a_mean$WT_LPS)
c_mean=b_mean[order(b_mean$dev,decreasing=TRUE),] 
a_mean <- a_mean[rownames(c_mean),]
Difdata1=diff_gene_Group2
row.names(Difdata1)=Difdata1$Gene
Difdata1<- Difdata1[rownames(c_mean),]
a_mean <- Difdata1[1:30,]

Plot<-pheatmap(log2(Difdata1[,8:19]+1), labels_row = Difdata1$Gene, cluster_rows = TRUE, cluster_cols = FALSE ,main = "",color= rdbu(100),scale="row")

pdf(paste0(Con1,"-",Con2,"-heatmap-reorder-top30absdev.pdf"),width=5,height=8)
print(Plot)
dev.off()


p<-pheatmap(log2(a_mean[,8:19]+1),cluster_rows=FALSE,show_rownames=TRUE,cluster_cols=FALSE,main ="",color=rdbu(100), cellwidth=20,cellheight=20,fontsize_row=10,fontsize_col=10,,file=paste0(Con1,"-",Con2,"-heatmap-count_matrix_mean_reorder-top30absdev.pdf"),width=10,height=15,fontsize_number=12, angle_col= "270")
p

