rm(list = ls());gc()

if(!dir.exists('C:/Users/10784/Desktop/SCCE/3.1.epi_dimRedu')){
  dir.create('C:/Users/10784/Desktop/SCCE/3.1.epi_dimRedu',recursive = T)}

setwd('C:/Users/10784/Desktop/SCCE/3.1.epi_dimRedu')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

grep('.rds',list.files('../2.1.annoMain'),value = T)
#[1] "SCCE_B.rds"             "SCCE_Basophil.rds"      "SCCE_Endothelial.rds"   "SCCE_Epithelial.rds"   
#[5] "SCCE_Fibroblast.rds"    "SCCE_mainCellTypes.rds" "SCCE_MB.rds"            "SCCE_Myeloid.rds"      
#[9] "SCCE_T.rds"

# load data ====
epi <- readRDS('../2.1.annoMain/SCCE_Epithelial.rds')
table(epi$Sample)
# Naive Treatment 
# 1820      5733 

## normalize ====
epi <- NormalizeData(object = epi, assay = 'RNA',
                     normalization.method = "LogNormalize",
                     scale.factor = 10000,
                     margin = 1, #If performing CLR normalization, normalize across features (1) or cells (2)
                     verbose = TRUE) # 进行标准化，默认参数

## find variable features ====
epi <- FindVariableFeatures(
  object = epi,
  assay = 'RNA',
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = "auto",
  num.bin = 20,
  binning.method = "equal_width",
  nfeatures = 5000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE
) #注意使用的是norm data

## 线性变化 ====
epi <- ScaleData(epi, features = rownames(epi))

## PCA降维 ====
epi <- RunPCA(epi, features = VariableFeatures(object = epi))

## 确定数据集维度
epi <- JackStraw(epi, num.replicate = 100)
epi <- ScoreJackStraw(epi, dims = 1:20)
#ElbowPlot(epi) # 20

## Cluster ====
epi <- FindNeighbors(epi, dims = 1:20)
epi <- FindClusters(epi, resolution = .1) %>% 
  FindClusters( resolution = .2) %>% FindClusters( resolution = .3) %>% 
  FindClusters( resolution = .4) %>% FindClusters( resolution = .5) %>% 
  FindClusters( resolution = .6) %>% FindClusters( resolution = .7) %>% 
  FindClusters( resolution = .8) %>% FindClusters( resolution = .9) %>% 
  FindClusters( resolution = 1) %>% FindClusters( resolution = 1.1) %>% 
  FindClusters( resolution = 1.2) %>% FindClusters( resolution = 1.3) %>% 
  FindClusters( resolution = 1.4)

## plot
library(clustree)
pdf('epi_clustree.pdf',width = 15,height = 11,onefile = F)
clustree(epi)
dev.off()

# 根据markers热图选resolution ====
lapply(seq(0.3,1,by=0.1), function(reso){
  
  gc()
  ## set resolution
  epi <- FindClusters(epi, resolution = reso)
  
  # UMAP & tSNE 降维 ====
  epi <- RunUMAP(epi, dims = 1:20)
  epi <- RunTSNE(epi, dims = 1:20)
  
  ## plot UMAP & t-SNE
  pdf(paste0('reso=',reso,'_UMAP.pdf'),width = 13,height = 4,onefile = F)
  print(
    DimPlot(epi, reduction = "umap", label = TRUE, 
            group.by = c("seurat_clusters", "Sample", "orig.ident"))
  )
  dev.off()
  pdf(paste0('reso=',reso,'_tSNE.pdf'),width = 13,height = 4,onefile = F)
  print(
    DimPlot(epi, reduction = "tsne", label = TRUE, 
            group.by = c("seurat_clusters", "Sample", "orig.ident"))
  )
  dev.off()
  
  # find markers ====
  allmarkers <- FindAllMarkers(
    object = epi, assay = 'integrated',
    only.pos= TRUE,min.pct = .25, logfc.threshold = 0.25) %>% 
    group_by(cluster)
  
  ## plot heatmap
  gc()
  pdf(paste0('reso=',reso,'_markers_heatmap.pdf'),onefile = F)
  print(DoHeatmap(epi, features = allmarkers$gene, cells = NULL,
                  group.by = "seurat_clusters", group.bar = TRUE,
                  group.colors = NULL, disp.min = -2.5, disp.max = NULL,
                  slot = "scale.data", assay = 'integrated', label = TRUE,
                  size = 5.5, hjust = 0, angle = 45, raster = TRUE,
                  draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02,
                  combine = TRUE ))
  dev.off()
})

epi <- FindClusters(epi, resolution = .6)
# Look at cluster IDs of the first 5 cells
head(Idents(epi), 5)
table(epi$seurat_clusters)

# UMAP & tSNE 降维 ====
epi <- RunUMAP(epi, dims = 1:20)
epi <- RunTSNE(epi, dims = 1:20)

pdf('epiCell_UMAP.pdf',width = 8,height = 3.5,onefile = F)
DimPlot(epi, reduction = "umap", label = TRUE, 
        group.by = c("Sample",'seurat_clusters'))
dev.off()
pdf('epiCell_tSNE.pdf',width = 8,height = 3.5,onefile = F)
DimPlot(epi, reduction = "tsne", label = TRUE, 
        group.by = c("Sample",'seurat_clusters'))
dev.off()

## plot cell ratio ====
library(plyr)
library(ggplot2)

plotData <- table(epi$Sample,epi$seurat_clusters) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Sample','Cluster','Freq')
plotData$Sample <- factor(plotData$Sample, levels = c('Naive','Treatment'))

plotData <- ddply(plotData, "Cluster", transform,
                  percent = Freq / sum(Freq) * 100)

pdf(paste0('tumorCell_Sample_Percent.pdf'),
    width = 7,height = 5,onefile = F)
print(ggplot(plotData, aes(x=Cluster, y=percent, fill=Sample)) +
        geom_bar(stat="identity") +
        scale_fill_manual(values=c('#a8e6cf','#fdffab','#ffd3b6','#ffaaa5',
                                   '#11999e','#3f72af','#ffc7c7','#9896f1',
                                   '#ff7e67','#0dceda','#6ef3d6','#a6d0e4',
                                   '#fbafaf','#f2c6b4','#f3e8cb','#99e1e5',
                                   '#ffc93c','#ff9a3c','#ff6f3c','#155263',
                                   '#c86b85','#e6a4b4','#f3d7ca','#bfcfff',
                                   '#a5dee5','#e0f9b5','#ffcfdf','#ffc7c7',
                                   '#355c7d','#6c5b7b','#c06c84','#f67280',
                                   '#edb1f1','#d59bf6','#8ef6e4','#99ddcc',
                                   '#ffe2e2','#bad7df','#cabbe9','#a1eafb',
                                   '#ff5722','#00adb5','#a2d5f2','#07689f'))+
        theme_bw()+labs(x='Cluster',y='Percent',fill='Sample'))
dev.off() 

plotData <- ddply(plotData, "Sample", transform,
                  percent = Freq / sum(Freq) * 100)

pdf(paste0('tumorCell_Cluster_Percent.pdf')
    ,width = 4,height = 5,onefile = F)
print(ggplot(plotData, aes(x=Sample, y=percent, fill=Cluster)) +
        geom_bar(stat="identity") +
        scale_fill_manual(values=c('#a8e6cf','#fdffab','#ffd3b6','#ffaaa5',
                                   '#11999e','#3f72af','#ffc7c7','#9896f1',
                                   '#ff7e67','#0dceda','#6ef3d6','#a6d0e4',
                                   '#fbafaf','#f2c6b4','#f3e8cb','#99e1e5',
                                   '#ffc93c','#ff9a3c','#ff6f3c','#155263',
                                   '#c86b85','#e6a4b4','#f3d7ca','#bfcfff',
                                   '#a5dee5','#e0f9b5','#ffcfdf','#ffc7c7',
                                   '#355c7d','#6c5b7b','#c06c84','#f67280',
                                   '#edb1f1','#d59bf6','#8ef6e4','#99ddcc',
                                   '#ffe2e2','#bad7df','#cabbe9','#a1eafb',
                                   '#ff5722','#00adb5','#a2d5f2','#07689f'))+
        
        theme_bw()+
        labs(x='Sample',y='Percent',fill='Cluster'))
dev.off() 

saveRDS(epi, 'SCCE_Epithelial_reCluster.rds')