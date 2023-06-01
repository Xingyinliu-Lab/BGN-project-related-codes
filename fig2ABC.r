
#————————————————————————————————————————————————————————————————————————————————————————————
# 这个文件包含下载目标数据，处理至筛选神经元的部分
#————————————————————————————————————————————————————————————————————————————————————————————


# https://cells.ucsc.edu/?ds=mouse-nervous-system+enteric
# article：Molecular Architecture of the Mouse Nervous System


# 查看ucsc数据集
	# 表达矩阵： exprMatrix.tsv.gz
	# 细胞元注释： meta.tsv
	# 降维坐标：tSNE.coords.tsv.gz
	# 数据集描述：desc.json
	# 细胞浏览器配置：dataset.json


R
# ① 读入数据
require(Seurat)
require(data.table)
mat <- fread("exprMatrix.tsv.gz")
meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
genes = mat[,1][[1]]
geness = gsub(".+[|]", "", genes) 
mat = data.frame(mat[,-1] , row.names=genes)
# new <- genes[duplicated(genes)]
library(readr)
library(tidyverse)
which(str_detect(geness,"Bgn"))
[1] 24611
# ENSMUSG00000031375|Bgn
# 数据集中测到的少于200个基因的细胞（min.features = 200）和少于3个细胞覆盖的基因（min.cells = 3）被过滤掉
so <- CreateSeuratObject(counts = mat, project = "MouseNervousSystemEnteric", meta.data=meta, min.cells = 3, min.features = 200)


# ② 质控
# 质控的参数主要有两个： 1.每个细胞测到的unique feature数目（unique feature代表一个细胞检测到的基因的数目，可以根据数据的质量进行调整） 2.每个细胞检测到的线粒体基因的比例，理论上线粒体基因组与核基因组相比，只占很小一部分。所以线粒体基因表达比例过高的细胞会被过滤。
so=readRDS("MouseNervousSystemEnteric-1.rds")
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# nFeature_RNA代表每个细胞测到的基因数目，nCount代表每个细胞测到所有基因的表达量之和，percent.mt代表测到的线粒体基因的比例。去除线粒体基因表达比例过高的细胞，和一些极值细胞。
 
so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
# 200为固定低表达值，其它数值根据图调整，筛除极值并保留绝大多数数据
# 这个数据里没有线粒体基因
# ③ 标准化
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
#鉴定细胞间表达量高变的基因（feature selection）
#这一步的目的是鉴定出细胞与细胞之间表达量相差很大的基因，用于后续鉴定细胞类型，
#我们使用默认参数，即“vst”方法选取2000个高变基因。
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)


so[["RNA"]]@var.features[1:10]
which(str_detect(VariableFeatures(so),"Bgn"))
# Bgn作为高变基因
Bgn="ENSMUSG00000031375-Bgn"


# ④ 细胞分类 
#Scaling the data
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)
#Perform linear dimensional reduction
so <- RunPCA(so, features = VariableFeatures(object = so))
#Examine and visualize PCA results a few different ways
print(so[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(so, dims = 1:2, reduction = "pca")
DimPlot(so, reduction = "pca")
DimHeatmap(so, dims = 1, cells = 500, balanced = TRUE)

# JackStraw和Elbow都可以决定数据的“维度”。
so <- JackStraw(so, num.replicate = 100)
so <- ScoreJackStraw(so, dims = 1:20)
JackStrawPlot(so, dims = 1:20)
ElbowPlot(so)

so <- FindNeighbors(so, dims = 1:20)
so <- FindClusters(so, resolution = 0.5) 


# 根据先验知识进行一个简单分类
# Ret,Chat,Nos1表达的细胞被判定为神经细胞
# 提取神经细胞进行下一步分析
neo = DotPlot(so, features = c(Ret,Chat,Nos1), assay='RNA' )
neo=neo$data
neo = neo[which((neo$pct.exp>50) & neo$avg.exp.scaled>-0.5),]
neo = neo[which(neo$features.plot == Ret),]$id %>% as.character()

so1 =subset(so, idents = neo)

# 把筛选出来的细胞重新分簇
# 这里经过调整，改为分7簇

all.genes <- rownames(so1)
so1 <- ScaleData(so1, features = all.genes)
#Perform linear dimensional reduction
So1 <- RunPCA(so1, features = VariableFeatures(object = so1))
ElbowPlot(so1)

so1 <- FindNeighbors(so1, dims = 1:20)
so1 <- FindClusters(so1, resolution = 0.5) 
so1 <- RunTSNE(so1, dims = 1:20)
head(so1@reductions$tsne@cell.embeddings)

p2 <- DimPlot(so1, reduction = "tsne")
Pt<-FeaturePlot(so1, features =Bgn)
pb<-VlnPlot(so1, features = Bgn)
# BGN的情况
p = DotPlot(so1, features = genes_to_check) + coord_flip()


# ⑤ 提取各个细胞类型的marker gene
#find all markers of cluster 1
cluster1.markers <- FindMarkers(so, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
#利用 DoHeatmap 命令可以可视化marker基因的表达
so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25)
?FindMarkers

library(dplyr)
top10 <- so.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pp<-DoHeatmap(so, features = top10$gene) + NoLegend()



# 神经元存档
saveRDS(so1, file = "neuro.rds")




#————————————————————————————————————————————————————————————————————————————————————————————
# 这个文件包含神经元注释和输出数据
#————————————————————————————————————————————————————————————————————————————————————————————


library(dplyr)
library(readr)
library(tidyverse)
library(Seurat)
library(data.table)

# 读取
# 从这里做分析

so1=readRDS("neuro.rds")

# BGN的图出两张
Bgn="ENSMUSG00000031375-Bgn"
Pd=DotPlot(so1, features = Bgn) + coord_flip()
pdf("Bgn-DotPlot.pdf",width=9,height=4)
print(Pd)
dev.off()

pb<-VlnPlot(so1, features = Bgn)
pdf("Bgn-VlnPlot.pdf",width=9,height=4)
print(pb)
dev.off()
 
 
 
# 按不同细胞簇输出数据
#利用system.time记录运行时间
for(case in c(1:6)){
assign(paste0("C",case),so1[,so1@meta.data$seurat_clusters %in% case])
system.time({fwrite(x=as.data.frame(get(paste0("C",case))[["RNA"]]@counts), row.names=T,file =paste0("C",case,".csv"))})
}

# 然后是根据目标基因的表达手动分类

# 参考文章的注释方法：

# 神经元分为两个主要部分，包括胆碱能亚群或氮能亚群（cholinergic or nitrergic）。这种广泛的分化与其他几个基因相关。例如，胶质细胞系来源的神经营养因子(GDNF)家族受体α1 (Gfra1)和α2 (Gfra2)的表达分别分离了Nos1和Chat表达的神经元。Gfra1和Gfra2是GDNF受体Ret的共同受体，Ret对ENS的形成是必要的。同样，Chat和Nos1表达亚群也差异表达了转录因子Casz1和Etv1。

# 基于Chat和Tac1的共同表达(Brookes et al.，14.54.51a)，我们对5个假定的兴奋性运动神经元(PEMNs)子集进行了注释，并在树状图的一棵子树中对它们进行了聚类(图S2D)。 PEMN亚群表达运动神经元标志物内源性阿片类物质脑啡肽(enkephalin, Penk)，和/或肌间运动神经元标志物calretinin (Calb2) 。

#基于Nos1和Vip的共同表达，我们注释了7个假定的抑制性运动神经元(PIMNs)亚群。总的来说，73%的vip阳性神经元共同表达Nos1。

# 小鼠分泌运动神经元和血管扩张神经元。Glp2r+假定分泌运动/血管扩张剂(PSVNs)的两个亚群 包括 Vip+非胆碱能亚群和Chat+胆碱能亚群(分别为PSVN1和PSVN2)。 PSVN2也表达在支配上皮和小动脉的神经元中表达的Galanin(Gal)、在分泌运动神经元中表达的neuropeptide Y (Npy)。 PSVN2中也有部分神经元表达谷氨酸脱羧酶2 (Gad2)，可能形成 胆碱能/GABA能 神经元。

# 肠内中间神经元(INs)传递感觉信息，协调运动神经元的活动。 目前已知的六种亚型是:(1)通过ACh、5HT和ATP信号的下降INs，(2)通过ATP信号的下降Nos1+Vip+Grp+Chat- INs，(3)通过ATP信号的下降Vip+Chat+Nos1+ INs，(4)通过Chat+Sst+ INs下降，(5)通过Sst响应的下降Penk+ INs，(6)通过ATP信号的上升Chat+Penk+ INs 

# 肠神经胶质（Plp1+ , Sox10+）



# 结果见judge.xlsx
celltype=c("PINs","PIMNs","PSVNs","PSVNs","PEMNs","PEMNs","PINs")

names(celltype) <- levels(so1)
so2<- RenameIdents(so1, celltype)

p2 <- DimPlot(so2, reduction = "tsne")
p2
pdf("annoted-11-9.pdf",width=5,height=4)
print(p2)
dev.off()


# 按不同细胞类型输出数据
PINs =so1[,so1@meta.data$seurat_clusters %in% c(0,6)]
system.time({fwrite(x=as.data.frame(PINs[["RNA"]]@counts), row.names=T,file =paste0("PINs.csv"))})
PIMNs=so1[,so1@meta.data$seurat_clusters %in% 1]
system.time({fwrite(x=as.data.frame(PIMNs[["RNA"]]@counts), row.names=T,file =paste0("PIMNs.csv"))})
PSVNs=so1[,so1@meta.data$seurat_clusters %in% c(2,3)]
system.time({fwrite(x=as.data.frame(PSVNs[["RNA"]]@counts), row.names=T,file =paste0("PSVNs.csv"))})
PEMNs=so1[,so1@meta.data$seurat_clusters %in% c(4,5)]
system.time({fwrite(x=as.data.frame(PEMNs[["RNA"]]@counts), row.names=T,file =paste0("PEMNs.csv"))})




