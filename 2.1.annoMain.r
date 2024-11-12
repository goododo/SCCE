rm(list = ls());gc()

if(!dir.exists('C:/Users/10784/Desktop/SCCE/2.1.annoMain')){
  dir.create('C:/Users/10784/Desktop/SCCE/2.1.annoMain',recursive = T)}

setwd('C:/Users/10784/Desktop/SCCE/2.1.annoMain')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

# load data ====
scce <- readRDS('../1.integrate/SCCE_inteData_dimRedu11.rds')
#DefaultAssay(scce) <- "integrated"

# cell annotation ====
# https://www.cellsignal.com/pathways/immune-cell-markers-human
# https://www.abcam.com/primary-antibodies/immune-cell-markers-poster
# https://lishensuo.github.io/posts/bioinfo/028%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--%E5%9F%BA%E4%BA%8E%E6%96%87%E7%8C%AE%E7%9A%84%E7%BB%86%E8%83%9E%E7%B1%BB%E5%9E%8B%E6%B3%A8%E9%87%8Amarker/
library(Hmisc)
# class markers
allMarkers <- list(
  EpithelialCells = c('KRT18','EPCAM','KRT19','KRT5','KRT15','CD24','CLDN4','MUC1',
                      'KRT8','KRT7',
                      'CALD1','AGRN','CTNNB1'),
  Leukocytes = c('CD45','PTPRC'),
  TCells = c('CD3','CD3D','CD3C','CD3E','TRAG','CD3G','CD2','GNLY'),
  
  BCells = c('CD19','MS4A1','CD79A','CD79B','B220','H2AB1','H2EB1','H2-AB1','H2-EB1'),
  
  Myeloid = c('CD11B','LYZ','CD74','LYZ1','LYZ2'),
  DC = c('CD11C','HLADR','CD83','CD317',
         'LY6C2','LAG3','HAVCR1','SELL','SIGLECH','IRF8'), #pDC
  Macrophages = c('CD68','MHCI','MHCII','CD163','MERTK','C1QA','MACRO','FOLR2','SIGLEC1','MS4A7',
                  'CD80','CD86','INOS'),
  Monocyte = c('CD14','CD115','CX3CR1','Ly6C','ITGAM','CSF1R','ADGRE1','CCR2','FCGR3','IFIH1','ISG15'),
  Basophil = c('CD22','CD123'),
  
  Progenitor = c('CD34','CD38','CD117'),
  HSC = c('CD135','CD49','CD49B','CD49D','CD90','THY1','SCA1','CD150'),
  
  Fibroblasts = c('COL1A1','COL3A1','DCN','POSTN','COL1A2','FN1','RUNX2','LUM'),
  
  Endothelial = c('VWF','RAMP2','DARC','PECAM1','KDR','EFNB2','ENG','PLVAP','CDH5')
)

## plot
library(ggplot2)
library(RColorBrewer)

pdf('mainCellMarkers.pdf', width = 20,height = 8,onefile = F)
DotPlot(scce, assay = 'RNA',
        features = allMarkers, group.by = "seurat_clusters",
        cols = c("#f8fcfb", "#2193b0"),
        col.min = 0)+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
dev.off()

#rm(list = ls());gc()
#scce <- readRDS('../1.integrate/SCCE_inteData_dimRedu.rds')

## plot
library(ggplot2)
library(RColorBrewer)

### violin plot
names(allMarkers)
lapply(names(allMarkers), function(cellType){
  pdf(paste0(cellType,'_vlnPlot.pdf'),width = 15,height = 15,
      onefile = F)
  print(
    VlnPlot(scce, #assay = 'RNA',
            features = allMarkers[[cellType]], pt.size = 0)
  )
  dev.off()
})

## add common cell types metaData ====
mainType <- ifelse(
  scce$seurat_clusters %in% c(0,2,3,5:13,18),'Epithelial',
  ifelse(scce$seurat_clusters %in% c(1,4,23,24),'T',
         ifelse(scce$seurat_clusters %in% c(14,20,21),'Myeloid',
                ifelse(scce$seurat_clusters %in% c(15),'B',
                       ifelse(scce$seurat_clusters %in% c(22),'MB',
                       ifelse(scce$seurat_clusters %in% c(16),'Endothelial',
                              ifelse(scce$seurat_clusters %in% c(17),'Fibroblast','Basophil'
         )))))))%>%
  factor(levels = c('Epithelial','T','Myeloid','MB','B','Basophil','Endothelial','Fibroblast'))

scce <- AddMetaData(scce,mainType,col.name = 'mainType')
table(scce$mainType)
table(scce$mainType,scce$Sample)

## plot
pdf('commonCellType_UMAP.pdf',width = 18,height = 5,onefile = F)
DimPlot(scce, reduction = "umap", label = TRUE, 
        group.by = c( "mainType",'seurat_clusters','Sample'))
dev.off()
pdf('commonCellType_tSNE.pdf',width = 18,height = 5,onefile = F)
DimPlot(scce, reduction = "tsne", label = TRUE, 
        group.by = c( "mainType",'seurat_clusters','Sample'))
dev.off()

## delete not good markers (common) ====
allMarkers <- list(
  EpithelialCells = c('KRT18','EPCAM','KRT19','KRT5','KRT15','CD24','CLDN4','MUC1',
                      'KRT8','KRT7'),
  Leukocytes = c('CD45','PTPRC'),
  TCells = c('CD3','CD3D','CD3C','CD3E','TRAG','CD3G','CD2','GNLY'),
  
  BCells = c('CD19','MS4A1','CD79A','CD79B','B220','H2AB1','H2EB1','H2-AB1','H2-EB1'),
  
  Myeloid = c('CD11B','LYZ','CD74','LYZ1','LYZ2'),
  Basophil = c('CD22','CD123'),
  
  Fibroblasts = c('COL1A1','COL3A1','DCN','POSTN','COL1A2','FN1','LUM'),
  
  Endothelial = c('VWF','RAMP2','DARC','PECAM1','KDR','EFNB2','ENG','PLVAP','CDH5')
)

### dot plot
pdf('mainCellMarkers_del.pdf',width = 16,height = 3,onefile = F)
DotPlot(scce, assay = 'RNA',
        features = allMarkers,group.by = "mainType",
        cols = c("#f8fcfb", "#f12711"),
        col.min = 0)+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
dev.off()

## plot cell type percent of each patient ====
library(plyr)
library(ggplot2)

# 以Date为切割变量()对每组数据进行transform()
plotData <- table(scce$Sample,scce$mainType) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Patient','Cluster','Freq')
plotData$Patient <- factor(plotData$Patient, levels = c('Naive','Treatment'))

plotData <- ddply(plotData, "Patient", transform,
                  percent = Freq / sum(Freq) * 100)

pdf('mainCell_Percent_Patient.pdf',width = 4,height = 4,onefile = F)
ggplot(plotData, aes(x=Patient, y=percent, fill=Cluster)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#a8e6cf','#fdffab','#ffd3b6','#ffaaa5',
                             '#11999e','#3f72af','#ffc7c7','#9896f1',
                             '#ff7e67','#0dceda','#6ef3d6','#a6d0e4'))+
  
  theme_bw()+
  labs(x='Tumor Type',y='Percent',fill='Cell Type')
dev.off() 

plotData <- ddply(plotData, 'Cluster',transform,
                  percent = Freq / sum(Freq) * 100)

pdf(paste0('mainCell_Cluster_Percent.pdf')
    ,width = 7,height = 5,onefile = F)
print(ggplot(plotData, aes(x=Cluster, y=percent, fill=Patient)) +
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

saveRDS(scce, 'SCCE_mainCellTypes.rds')
# subclusters ====
lapply(as.character(mainType) %>% unique, function(celltype){
  subdata <- subset(scce, mainType == celltype)
  print(celltype)
  saveRDS(subdata, paste0('SCCE_',celltype,'.rds'))
})