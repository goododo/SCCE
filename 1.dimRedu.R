#esophageal small-cell carcinoma
rm(list = ls());gc()

if(!dir.exists('SCCE/1.integrate')){
  dir.create('SCCE/1.integrate',recursive = T)}

setwd('SCCE/1.integrate')
# library ====
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(data.table)

# load data ====
SCCE <- readRDS('../0.rawdata/SCCE_Raw_SeuratObj_seuratIntegration.rds')

DefaultAssay(SCCE) <- "integrated"
# Run the standard workflow for visualization and clustering
SCCE <- ScaleData(SCCE, verbose = T, features = rownames(SCCE))

# PCA 降维 ====
SCCE <- RunPCA(SCCE, npcs = 50, 
               verbose = T,
               features = VariableFeatures(object = SCCE))

## 确定数据集维度
SCCE <- JackStraw(SCCE, num.replicate = 100)
SCCE <- ScoreJackStraw(SCCE, dims = 1:20)
#JackStrawPlot(SCCE, dims = 1:20)
pdf('inteData_PCA_elbow.pdf',width = 6,height = 4,onefile = F)
ElbowPlot(SCCE) # 20
dev.off()

SCCE <- RunUMAP(SCCE, reduction = "pca", dims = 1:20)

# Cluster ====
SCCE <- FindNeighbors(SCCE,#reduction = "harmony",
                      dims = 1:20)

SCCE <- FindClusters(SCCE, resolution = .1) %>% 
  FindClusters( resolution = .2) %>% FindClusters( resolution = .3) %>%
  FindClusters( resolution = .4) %>% FindClusters( resolution = .5) %>% 
  FindClusters( resolution = .6) %>% FindClusters( resolution = .7) %>% 
  FindClusters( resolution = .8) %>% FindClusters( resolution = .9) %>% 
  FindClusters( resolution = 1) %>% FindClusters( resolution = 1.1) %>% 
  FindClusters( resolution = 1.2) %>% FindClusters( resolution = 1.3) %>% 
  FindClusters( resolution = 1.4)

## plot
library(clustree)
pdf('inteData_clustree.pdf',width = 15,height = 11,onefile = F)
clustree(SCCE)
dev.off()

# 根据markers热图选resolution ====
lapply(seq(0.3,1.4,by=0.1), function(reso){
  
  gc()
  ## set resolution
  SCCE <- FindClusters(SCCE, resolution = reso)
  
  # UMAP & tSNE 降维 ====
  SCCE <- RunUMAP(SCCE, dims = 1:20)
  SCCE <- RunTSNE(SCCE, dims = 1:20)
  
  ## plot UMAP & t-SNE
  pdf(paste0('reso=',reso,'_UMAP.pdf'),width = 9.5,height = 4,onefile = F)
  print(
    DimPlot(SCCE, reduction = "umap", label = TRUE, 
            group.by = c("seurat_clusters", "orig.ident"))
  )
  dev.off()
  pdf(paste0('reso=',reso,'_tSNE.pdf'),width = 9.5,height = 4,onefile = F)
  print(
    DimPlot(SCCE, reduction = "tsne", label = TRUE, 
            group.by = c("seurat_clusters", "orig.ident"))
  )
  dev.off()
  
  # find markers ====
  allmarkers <- FindAllMarkers(
    object = SCCE, assay = 'integrated',
    only.pos= TRUE,min.pct = .25, logfc.threshold = 0.25) %>% 
    group_by(cluster)
  
  ## plot heatmap
  gc()
  pdf(paste0('reso=',reso,'_markers_heatmap.pdf'),onefile = F)
  print(DoHeatmap(SCCE, features = allmarkers$gene, cells = NULL,
                  group.by = "seurat_clusters", group.bar = TRUE,
                  group.colors = NULL, disp.min = -2.5, disp.max = NULL,
                  slot = "scale.data", assay = 'integrated', label = TRUE,
                  size = 5.5, hjust = 0, angle = 45, raster = TRUE,
                  draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02,
                  combine = TRUE ))
  dev.off()
})

## set resolution
SCCE <- FindClusters(SCCE, resolution = 1.1)

# UMAP & tSNE 降维 ====
SCCE <- RunUMAP(SCCE, dims = 1:20)
SCCE <- RunTSNE(SCCE, dims = 1:20)

## plot
pdf('inteData_UMAP.pdf',width = 12,height = 5,onefile = F)
DimPlot(SCCE, reduction = "umap", label = TRUE, 
        group.by = c("seurat_clusters", "orig.ident"))
dev.off()
pdf('inteData_tSNE.pdf',width = 12,height = 5,onefile = F)
DimPlot(SCCE, reduction = "tsne", label = TRUE, 
        group.by = c("seurat_clusters", "orig.ident"))
dev.off()

## plot cell ratio ====
library(plyr)
library(ggplot2)

plotData <- table(SCCE$Sample,SCCE$seurat_clusters) %>% data.frame()# %>% 
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

plotData <- ddply(plotData, 'Sample',transform,
                  percent = Freq / sum(Freq) * 100)

pdf(paste0('tumorCell_Cluster_Percent.pdf')
    ,width = 7,height = 5,onefile = F)
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

saveRDS(SCCE,'SCCE_inteData_dimRedu11.rds')
#29620 features across 10932 samples within 2 assays 
