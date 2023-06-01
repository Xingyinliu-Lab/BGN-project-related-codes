G1="KO_C"
G2="WT_C"
G3="KO_LPS"
G4="WT_LPS"

S12=data.frame(read.table(paste0("./xyy-小鼠转录组/2-筛后sheet/",G1,"-",G2,"-padj0.05-lfc1-resdata.csv"), sep=",", header = TRUE, row.names = 1))
S31=data.frame(read.table(paste0("./xyy-小鼠转录组/2-筛后sheet/",G3,"-",G1,"-padj0.05-lfc1-resdata.csv"),sep=",", header = TRUE, row.names = 1))
S34=data.frame(read.table(paste0("./xyy-小鼠转录组/2-筛后sheet/",G3,"-",G4,"-padj0.05-lfc1-resdata.csv"),sep=",", header = TRUE, row.names = 1))
S42=data.frame(read.table(paste0("./xyy-小鼠转录组/2-筛后sheet/",G4,"-",G2,"-padj0.05-lfc1-resdata.csv"),sep=",", header = TRUE, row.names = 1))

S12$group=paste0(G1,"-",G2)
S12$groupA=G1
S12$groupB=G2
S31$group=paste0(G3,"-",G1)
S31$groupA=G3
S31$groupB=G1
S34$group=paste0(G3,"-",G4)
S34$groupA=G3
S34$groupB=G4
S42$group=paste0(G4,"-",G2)
S42$groupA=G4
S42$groupB=G2

All<-rbind(S12,S31,S34,S42)

library(stringr)
All[str_order(All$Gene), ]-> resdata

resdata$significant <- as.factor(resdata$padj<0.001 & abs(resdata$log2FoldChange)>1.5)
resdata$threshold = factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) >= 1, ifelse(resdata$log2FoldChange>=1,'Up','Down'),'NS'),levels=c('Up','Down','NS'))

C1= rownames(colData)[which(colData$condition==G1)]
C2= rownames(colData)[which(colData$condition==G2)]
C3= rownames(colData)[which(colData$condition==G3)]
C4= rownames(colData)[which(colData$condition==G4)]

resdataS= resdata[which(resdata$significant=="TRUE"),]
resdataS$KO_C = rowMeans(resdataS[,which(colnames(resdataS) %in% C1),])
resdataS$WT_C= rowMeans(resdataS[,which(colnames(resdataS) %in% C2),])
resdataS$KO_LPS= rowMeans(resdataS[,which(colnames(resdataS) %in% C3),])
resdataS$WT_LPS= rowMeans(resdataS[,which(colnames(resdataS) %in% C4),])

resdataS[,c(1,3,7,20:28)]->data
write.csv(data, file="All-Diffdata.csv")


A <- data[which(data$group=="KO_C-WT_C"),]$Gene
B <- data[which(data$group=="KO_LPS-KO_C"),]$Gene
C <- data[which(data$group=="KO_LPS-WT_LPS"),]$Gene
D <- data[which(data$group=="WT_LPS-WT_C"),]$Gene



venn.plot <- venn.diagram(
#数据列表
x = list(
"KO_C-WT_C" = A,
"KO_LPS-KO_C" = B,
"KO_LPS-WT_LPS" = C ,
"WT_LPS-WT_C" = D
),
filename = "Venn_4set.tiff",    #保存路径
col = "transparent",      #指定图形的圆周边缘颜色  transparent 透明           
fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),  #填充颜色
alpha = 0.50,                                      #透明度
label.col = c("orange", "white", "darkorchid4", "white",
"white", "white", "white", "white", "darkblue", "white",
"white", "white", "white", "darkgreen", "white"),
cex = 1.5,    #每个区域label名称的大小
fontfamily = "serif",  #字体
fontface = "bold",     #字体格式
cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),  #分类颜色 
cat.cex = 1.5,      #每个分类名称大小
cat.pos = 0,        #
cat.dist = 0.09,    #
cat.fontfamily = "serif",     #分类字体
cat.default.pos = "outer",
rotation.degree = 270,        #旋转角度
margin = 0.2               #在网格单元中给出图周围空白量的编号
);


# 更换网站画：Oliveros, J.C. (2007-2015) Venny. An interactive tool for comparing lists with Venn's diagrams. https://bioinfogp.cnb.csic.es/tools/venny/index.html
