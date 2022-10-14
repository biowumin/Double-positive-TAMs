################################################### integrate data from different patients using LIGER
rm(list=ls())
library(Seurat)
library(rliger)
library(tidyverse)
library(SeuratWrappers)
library(UCell)
library(RColorBrewer)
packageVersion("rliger")   ## version 1.0.0
setwd("/media/usb3/Couturier_scRNA/GBM_res")
readRDS("GBM_QC.rds") -> GBM_QC
GBM_tumor <- subset(GBM_QC,subset=orig.ident%in%c("BT389","BT390","BT397","BT400","BT402","BT407","BT409"))  ## whole tumor without sorting
obj.list <- SplitObject(GBM_tumor, split.by = "orig.ident")
tumor_list <- lapply(obj.list,function(x){x <- CreateSeuratObject(x@assays$RNA@counts)})
whole_tumor <- merge(tumor_list[[1]],tumor_list[2:length(tumor_list)])
whole_tumor <- NormalizeData(whole_tumor) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% ScaleData(split.by = "orig.ident", do.center = FALSE)
nFactors <- 25
whole_tumor <- RunOptimizeALS(whole_tumor,k=nFactors,split.by="orig.ident")
whole_tumor <- RunQuantileNorm(whole_tumor, split.by="orig.ident")
whole_tumor <- FindNeighbors(whole_tumor, reduction="iNMF", dims=1:nFactors) %>% FindClusters(resolution = 0.4)
whole_tumor <- RunUMAP(whole_tumor, dims=1:nFactors, reduction="iNMF")
saveRDS(whole_tumor,file="whole_tumor2.RDS")

readRDS("whole_tumor2.RDS") -> whole_tumor
setwd("/media/usb3/Couturier_scRNA/GBM_res/myscore")
load("TAM_tumor_score2.RData")
whole_tumor$tumor_score <- score$tumor_norm
whole_tumor$TAM_score <- score$TAM_norm
whole_tumor$immune_suppress_score <- score$immune_suppress_norm
whole_tumor$PS_score <- score$PS_norm
tiff("whole_tumor_umap.tif",res=800,compression="lzw",units="in",width=6,height = 6)
DimPlot(whole_tumor, reduction = "umap",label=T,pt.size = 1,label.size = 6) + NoLegend() + NoAxes()
dev.off()
tiff("patient_umap.tif",res=800,compression="lzw",units="in",width=6,height = 6)
DimPlot(whole_tumor, reduction = "umap",pt.size = 1,group.by = "orig.ident") + NoAxes() + labs(title="")
dev.off()
## cell annotation
whole_tumor$cell_type <- as.character(whole_tumor$seurat_clusters)
whole_tumor$cell_type[whole_tumor$cell_type%in%c("2","3","7","8","9","10","11","13","14","17","18")] <- "tumor"
whole_tumor$cell_type[whole_tumor$cell_type%in%c("12","16")] <- "oligo"
whole_tumor$cell_type[whole_tumor$cell_type%in%c("19")] <- "T_cell"
whole_tumor$cell_type[whole_tumor$cell_type%in%c("20")] <- "endothelial"
whole_tumor$cell_type[whole_tumor$cell_type%in%c("0","1")] <- "microglia"
whole_tumor$cell_type[whole_tumor$cell_type%in%c("4","5","6","15")] <- "BMDM"
mycolors <- c(brewer.pal(12,'Paired'),"lightgray")
tiff("cell_anno.tif",res=800,compression="lzw",units="in",width=6,height = 6)
DimPlot(whole_tumor, reduction = "umap",label=F,pt.size = 1,group.by = "cell_type",cols=mycolors[c(1,3,5,7,8,9)]) + NoLegend() + NoAxes() + labs(title="")
dev.off()
pdf("cell_anno_marker.pdf")
VlnPlot(whole_tumor, features = c("MCAM","ESAM","PTPRC","CD3D"),group.by = "cell_type",ncol=2,cols=mycolors[c(1,3,5,7,8,9)],pt.size = 0)
VlnPlot(whole_tumor, features = c("P2RY12","TMEM119","LYZ","TGFBI"),group.by = "cell_type",ncol=2,cols=mycolors[c(1,3,5,7,8,9)],pt.size = 0)
VlnPlot(whole_tumor, features = c("SOX9","BCAN","MOG","MAG"),group.by = "cell_type",ncol=2,cols=mycolors[c(1,3,5,7,8,9)],pt.size = 0)
dev.off()

cells <- colnames(subset(whole_tumor,subset=cell_type%in%c("tumor","BMDM","microglia")))
tiff("tumor_score.tif",res=800,compression="lzw",units="in",width=6,height = 6)
pdf("tumor_score.pdf")
FeaturePlot(whole_tumor,cells=cells, features = "tumor_score",pt.size=1,order=F) + scale_color_gradient2(low = "lightgray", high = "orangered",midpoint = 0.4) + NoAxes() + labs(title="")
dev.off()
tiff("TAM_score.tif",res=800,compression="lzw",units="in",width=6,height = 6)
pdf("TAM_score.pdf")
FeaturePlot(whole_tumor,cells=cells, features = "TAM_score",pt.size=1,order=F,min.cutoff = 0.35) + scale_color_gradient2(low = "lightgray", high = "orangered",midpoint = 0.5) + NoAxes() + labs(title="")
dev.off()

library(ggrepel)
library(RColorBrewer)
library(ggExtra)
library(ggplot2)
library(ggpubr)
subset(whole_tumor,subset=cell_type%in%c("tumor","BMDM","microglia")) -> zz
score <- data.frame(TAM_norm=zz$TAM_score,tumor_norm=zz$tumor_score,cell_type=zz$cell_type)
score$cell_type <- ifelse(score$cell_type=="tumor","tumor","TAM")
p <- ggplot(score,aes(x=TAM_norm,y=tumor_norm,colour=cell_type)) + 
  geom_point(size=0.8) + scale_color_manual(values=c("blue","red")) + theme_bw() + geom_hline(yintercept = 0.45)+geom_vline(xintercept = 0.6) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))
p1 <- ggMarginal(p, type="density", margins = "both",groupFill=T,groupColour=T)
ggsave("TAM_tumor_score_density.pdf",p1)

macrophage <- subset(whole_tumor,subset = cell_type%in%c("BMDM","microglia"))
macrophage$status <- ifelse(macrophage$tumor_score>0.45,"dou-pos","non-pos")
table(macrophage$status,macrophage$seurat_clusters)
table(macrophage$cell_type,macrophage$status)
tiff("cell_status.tif",res=800,compression="lzw",units="in",width=6,height = 6)
pdf("cell_status.pdf")
DimPlot(macrophage, group.by = "status",pt.size=1,order=F,cols=c("orangered","lightgray")) + NoLegend() + NoAxes() + labs(title="")
dev.off()
tiff("BMDM_MG_marker.tif",res=800,compression="lzw",units="in",width=6,height = 6)
pdf("BMDM_MG_marker.pdf")
FeaturePlot(macrophage,features = c("P2RY12","TMEM119","TGFBI","AHR"),pt.size=1,order=T,min.cutoff=1,cols=c("lightgray","orangered"))
dev.off()
tiff("immunesuppress_marker.tif",res=800,compression="lzw",units="in",width=6,height = 3)
pdf("immunesuppress_marker.pdf")
FeaturePlot(macrophage,ncol=2,features = c("MRC1","MARCO","MRC1","MARCO"),pt.size=1,order=T,min.cutoff=1.5,cols=c("lightgray","orangered"))
dev.off()