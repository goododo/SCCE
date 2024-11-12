rm(list = ls());gc()

if(!dir.exists('/home/gzy//SCCE/2.2.annoImmune')){
  dir.create('/home/gzy//SCCE/2.2.annoImmune',recursive = T)}

setwd('/home/gzy//SCCE/2.2.annoImmune')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

# load data ====
immune <- readRDS('../2.1.annoMain/SCCE_Leukocyte.rds')
#DefaultAssay(immune) <- "integrated"

## normalize ====
immune <- NormalizeData(object = immune, assay = 'RNA',
                     normalization.method = "LogNormalize",
                     scale.factor = 10000,
                     margin = 1, #If performing CLR normalization, normalize across features (1) or cells (2)
                     verbose = TRUE) # 进行标准化，默认参数

## find variable features ====
immune <- FindVariableFeatures(
  object = immune,
  assay = 'RNA',
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = "auto",
  num.bin = 20,
  binning.method = "equal_width",
  nfeatures = 2000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE
) #注意使用的是norm data

## 线性变化 ====
immune <- ScaleData(immune, features = rownames(immune))

## PCA降维 ====
immune <- RunPCA(immune, features = VariableFeatures(object = immune))

## 确定数据集维度
immune <- JackStraw(immune, num.replicate = 100)
immune <- ScoreJackStraw(immune, dims = 1:20)
ElbowPlot(immune) # 20

## Cluster ====
immune <- FindNeighbors(immune, dims = 1:20)
immune <- FindClusters(immune, resolution = .1) %>% 
  FindClusters( resolution = .2) %>% FindClusters( resolution = .3) %>% 
  FindClusters( resolution = .4) %>% FindClusters( resolution = .5) %>% 
  FindClusters( resolution = .6) %>% FindClusters( resolution = .7) %>% 
  FindClusters( resolution = .8) %>% FindClusters( resolution = .9) %>% 
  FindClusters( resolution = 1) %>% FindClusters( resolution = 1.1) %>% 
  FindClusters( resolution = 1.2) %>% FindClusters( resolution = 1.3) %>% 
  FindClusters( resolution = 1.4)

## plot
library(clustree)
pdf('immune_clustree.pdf',width = 15,height = 11,onefile = F)
clustree(immune)
dev.off()

# 根据markers热图选resolution ====
lapply(seq(0.3,1.5,by=0.1), function(reso){
  
  gc()
  ## set resolution
  immune <- FindClusters(immune, resolution = reso)
  
  # UMAP & tSNE 降维 ====
  immune <- RunUMAP(immune, dims = 1:20)
  immune <- RunTSNE(immune, dims = 1:20)
  
  ## plot UMAP & t-SNE
  pdf(paste0('reso=',reso,'_UMAP.pdf'),width = 10,height = 4,onefile = F)
  print(
    DimPlot(immune, reduction = "umap", label = TRUE, 
            group.by = c("seurat_clusters", "orig.ident"))
  )
  dev.off()
  pdf(paste0('reso=',reso,'_tSNE.pdf'),width = 10,height = 4,onefile = F)
  print(
    DimPlot(immune, reduction = "tsne", label = TRUE, 
            group.by = c("seurat_clusters", "orig.ident"))
  )
  dev.off()
  
  # find markers ====
  allmarkers <- FindAllMarkers(
    object = immune, assay = 'integrated',
    only.pos= TRUE,min.pct = .25, logfc.threshold = 0.25) %>% 
    group_by(cluster)
  
  ## plot heatmap
  gc()
  pdf(paste0('reso=',reso,'_markers_heatmap.pdf'),onefile = F)
  print(DoHeatmap(immune, features = allmarkers$gene, cells = NULL,
                  group.by = "seurat_clusters", group.bar = TRUE,
                  group.colors = NULL, disp.min = -2.5, disp.max = NULL,
                  slot = "scale.data", assay = 'integrated', label = TRUE,
                  size = 5.5, hjust = 0, angle = 45, raster = TRUE,
                  draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02,
                  combine = TRUE ))
  dev.off()
})

immune <- FindClusters(immune, resolution = .5)
# Look at cluster IDs of the first 5 cells
head(Idents(immune), 5)
table(immune$seurat_clusters)

# UMAP & tSNE 降维 ====
immune <- RunUMAP(immune, dims = 1:20)
immune <- RunTSNE(immune, dims = 1:20)

# cell annotation ====
# https://www.cellsignal.com/pathways/immune-cell-markers-human
# https://www.abcam.com/primary-antibodies/immune-cell-markers-poster
# https://lishensuo.github.io/posts/bioinfo/028%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--%E5%9F%BA%E4%BA%8E%E6%96%87%E7%8C%AE%E7%9A%84%E7%BB%86%E8%83%9E%E7%B1%BB%E5%9E%8B%E6%B3%A8%E9%87%8Amarker/
library(Hmisc)

humanMarkers <- list(
  #Lymphoid
    TCells = c('CD3','CD3D','CD3C','CD3E','TRAG','CD3G','CD2','GNLY'),
    Th = c('CD4'),
    Tcyto = c('CD8','CD8A','CD8B'),
    Treg = c('FOXP3','IL2RA','IKZF2'),
    Activated = c('CD69','CD25'),
    Naive = c('CD45RA','LEF1','SELL','TCF7'),
    Cytotoxic = c('IFNG','PRF1','GZMK','GZMA','GZMB','NKG7'),
    Proliferation = c('MKI67','PCNA','TOP2A'),
    Effector = c('CCR7'),
    Exhaustion = c('PD1','TIM3','LAG3','TIGIT','TOX','TOX2','CTLA4','HAVCR2','PDCD1'),
    #NK = c('CD56','CD16','CD94','NKP46'),
    BCells = c('CD19','MS4A1','CD79A','CD79B'),
    M_B = c('IL10'), #myeloid-like B: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9930167/
    
  Myeloid = c('CD11B','LYZ','CD74'),
    DC = c('CD11C','HLADR','CD83','CD317'),
    Macrophages = c('CD68','MHCII','CD163'),
    M1 = c('CD80','CD86','INOS'),
    #M2 = c('CD206'),
    Monocyte = c('CD14'),
    Mast = c('CD32','CD33','CD117','CD203','CD203C','FCERI','FCERIA'),
    Neutrophil = c('CD18','CD44','CD55'),
    Eosinophil = c('CD45','CD125','CD193','F4/80','F480','F4','F80','EMR1','SIGLEC8'),
    Basophil = c('CD22','CD123'),

  Progenitor = c('CD34','CD38'), # cd34 only: Multi-potent progenitor (MPP); cd38 only: Megakaryocyte–erythroid progenitor
    HSC = c('CD49','CD49B','CD49D','CD90','THY1')
)

## plot
library(ggplot2)
library(RColorBrewer)

pdf('immuneCellMarkers.pdf', width = 25,height = 8,onefile = F)
DotPlot(immune, assay = 'RNA',
        features = humanMarkers, group.by = "seurat_clusters",
        cols = c("#f8fcfb", "#2193b0"),
        col.min = 0)+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
dev.off()

## plot
library(ggplot2)
library(RColorBrewer)

### violin plot
names(humanMarkers)
lapply(names(humanMarkers), function(cellType){
  pdf(paste0(cellType,'_vlnPlot.pdf'),width = 15,height = 15,
      onefile = F)
  print(
    VlnPlot(immune, #assay = 'RNA',
            features = humanMarkers[[cellType]], pt.size = 0)
  )
  dev.off()
})

## add common cell types metaData ====
ImmuneType <- ifelse(immune$seurat_clusters %in% c(0,5,9),'Tcyto',
                   ifelse(immune$seurat_clusters %in% c(1,4),'Th',
                          ifelse(immune$seurat_clusters %in% c(2),'DC',
                                 ifelse(immune$seurat_clusters %in% c(3,6),'Neutrophil',
                                        ifelse(immune$seurat_clusters %in% c(7),'Myeloid',
                                 'B'))))) %>%
  factor(levels = c('Tcyto','Th','B','Myeloid','DC','Neutrophil'))

immune <- AddMetaData(immune,ImmuneType,col.name = 'ImmuneType')
table(immune$ImmuneType)
table(immune$ImmuneType,immune$orig.ident)

## delete not good markers (common) ====
humanMarkers <- list(
  #Lymphoid
    TCells = c('CD3D','CD3E','CD3G','CD2','GNLY'),
    Th = c('CD4'),
    Tcyto = c('CD8','CD8A','CD8B','IFNG','PRF1','GZMK','GZMA','GZMB','NKG7'),
    Treg = c('FOXP3','IL2RA','IKZF2'),
    Proliferation = c('MKI67','PCNA','TOP2A'),
    Exhaustion = c('PD1','TIM3','LAG3','TIGIT','TOX','TOX2','CTLA4','HAVCR2','PDCD1'),
    #NK = c('CD56','CD16','CD94','NKP46'),
    BCells = c('CD19','MS4A1','CD79A','CD79B'),

  Myeloid = c('CD11B','LYZ','CD74'),
    DC = c('CD11C','HLADR','CD83','CD317'),
    Neutrophil = c('CD18','CD44','CD55')
)

### dot plot
pdf('immuneCellMarkers_del.pdf',width = 12,height = 3,onefile = F)
DotPlot(immune, assay = 'RNA',
        features = humanMarkers,group.by = "ImmuneType",
        cols = c("#f8fcfb", "#f12711"),
        col.min = 0)+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
dev.off()

## plot
pdf('immuneCellType_UMAP.pdf',width = 20,height = 5.5,onefile = F)
DimPlot(immune, reduction = "umap", label = TRUE, 
        group.by = c("orig.ident", "ImmuneType",'seurat_clusters'))
dev.off()
pdf('immuneCellType_tSNE.pdf',width = 20,height = 5.5,onefile = F)
DimPlot(immune, reduction = "tsne", label = TRUE, 
        group.by = c("orig.ident", "ImmuneType",'seurat_clusters'))
dev.off()

## plot cell type percent of each patient ====
library(plyr)
library(ggplot2)

plotData <- table(immune$orig.ident,immune$ImmuneType) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Sample','Cluster','Freq')
plotData$Sample <- factor(plotData$Sample, levels = c('1N', '1T'))

plotData <- ddply(plotData, "Cluster", transform,
                  percent = Freq / sum(Freq) * 100)

pdf(paste0('immuneCell_Sample_Percent.pdf'),
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
        theme_bw()+labs(x='Cluster',y='Percent',fill='orig.ident'))
dev.off() 

plotData <- ddply(plotData, "Sample", transform,
                  percent = Freq / sum(Freq) * 100)

pdf(paste0('immuneCell_Cluster_Percent.pdf')
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
        labs(x='orig.ident',y='Percent',fill='Cluster'))
dev.off() 

saveRDS(immune, 'SCCE_Immune_reAnno.rds')