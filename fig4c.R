data_Diff<-read.delim(file="relative_abundance_p0.05.csv",sep=",",row.names= 1,header=TRUE)
library(RColorBrewer)
rdbu <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))
library(pheatmap)
annotation_col <-  read.csv(file='anno.txt',sep="\t",row.names= 1,header=TRUE)
pheatmap(data_Diff,annotation_col = annotation_col, cluster_rows = TRUE,cutree_rows = 2, cluster_cols = FALSE,gaps_col=12,main = "",color= rdbu(100),scale="row",cellwidth=20,cellheight=20,fontsize_row=10,fontsize_col=10,,file="diff-heatmap.pdf",width=40,height=15,fontsize_number=12)
dev.off()

