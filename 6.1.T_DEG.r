rm(list = ls());gc()

if(!dir.exists('C:/Users/10784/Desktop/SCCE/7.1.T_DEG')){
  dir.create('C:/Users/10784/Desktop/SCCE/7.1.T_DEG',recursive = T)}

setwd('C:/Users/10784/Desktop/SCCE/7.1.T_DEG')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)

grep('.rds',list.files('../2.1.annoMain'),value = T)
# load data ====
tcell <- readRDS('../2.1.annoMain/SCCE_T.rds')
#24324 features across 828 samples

## normalize ====
tcell <- NormalizeData(object = tcell, assay = 'RNA',
                       normalization.method = "LogNormalize",
                       scale.factor = 10000,
                       margin = 1, #If performing CLR normalization, normalize across features (1) or cells (2)
                       verbose = TRUE) # 进行标准化，默认参数

## find variable features ====
tcell <- FindVariableFeatures(
  object = tcell,
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
tcell <- ScaleData(tcell, features = rownames(tcell))

## PCA降维 ====
tcell <- RunPCA(tcell, features = VariableFeatures(object = tcell))

## 确定数据集维度
tcell <- JackStraw(tcell, num.replicate = 100)
tcell <- ScoreJackStraw(tcell, dims = 1:20)
ElbowPlot(tcell) # 20

## Cluster ====
tcell <- FindNeighbors(tcell, dims = 1:20)
tcell <- FindClusters(tcell, resolution = .1) %>% 
  FindClusters( resolution = .2) %>% FindClusters( resolution = .3) %>% 
  FindClusters( resolution = .4) %>% FindClusters( resolution = .5) %>% 
  FindClusters( resolution = .6) %>% FindClusters( resolution = .7) %>% 
  FindClusters( resolution = .8) %>% FindClusters( resolution = .9) %>% 
  FindClusters( resolution = 1) %>% FindClusters( resolution = 1.1) %>% 
  FindClusters( resolution = 1.2) %>% FindClusters( resolution = 1.3) %>% 
  FindClusters( resolution = 1.4)

## plot
library(clustree)
pdf('tcell_clustree.pdf',width = 15,height = 11,onefile = F)
clustree(tcell)
dev.off()

# 根据markers热图选resolution ====
lapply(seq(0.1,1.4,by=0.1), function(reso){
  
  gc()
  ## set resolution
  tcell <- FindClusters(tcell, resolution = reso)
  
  # UMAP & tSNE 降维 ====
  tcell <- RunUMAP(tcell, dims = 1:20)
  tcell <- RunTSNE(tcell, dims = 1:20)
  
  ## plot UMAP & t-SNE
  pdf(paste0('reso=',reso,'_UMAP.pdf'),width = 10,height = 4,onefile = F)
  print(
    DimPlot(tcell, reduction = "umap", label = TRUE, 
            group.by = c("seurat_clusters", "orig.ident"))
  )
  dev.off()
  pdf(paste0('reso=',reso,'_tSNE.pdf'),width = 10,height = 4,onefile = F)
  print(
    DimPlot(tcell, reduction = "tsne", label = TRUE, 
            group.by = c("seurat_clusters", "orig.ident"))
  )
  dev.off()
  
  # find markers ====
  allmarkers <- FindAllMarkers(
    object = tcell, assay = 'integrated',
    only.pos= TRUE,min.pct = .25, logfc.threshold = 0.25) %>% 
    group_by(cluster)
  
  ## plot heatmap
  gc()
  pdf(paste0('reso=',reso,'_markers_heatmap.pdf'),onefile = F)
  print(DoHeatmap(tcell, features = allmarkers$gene, cells = NULL,
                  group.by = "seurat_clusters", group.bar = TRUE,
                  group.colors = NULL, disp.min = -2.5, disp.max = NULL,
                  slot = "scale.data", assay = 'integrated', label = TRUE,
                  size = 5.5, hjust = 0, angle = 45, raster = TRUE,
                  draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02,
                  combine = TRUE ))
  dev.off()
})

tcell <- FindClusters(tcell, resolution = 1.4)
# Look at cluster IDs of the first 5 cells
head(Idents(tcell), 5)
table(tcell$seurat_clusters)

# UMAP & tSNE 降维 ====
tcell <- RunUMAP(tcell, dims = 1:20)
tcell <- RunTSNE(tcell, dims = 1:20)

# cell annotation ====
# https://www.cellsignal.com/pathways/tcell-cell-markers-human
# https://www.abcam.com/primary-antibodies/tcell-cell-markers-poster
# https://lishensuo.github.io/posts/bioinfo/028%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--%E5%9F%BA%E4%BA%8E%E6%96%87%E7%8C%AE%E7%9A%84%E7%BB%86%E8%83%9E%E7%B1%BB%E5%9E%8B%E6%B3%A8%E9%87%8Amarker/
library(Hmisc)
humanMarkers <- list(
  #Lymphoid
  Th = c('CD4','CCR6'),
  Tcyto = c('CD8','CD8A','CD8B','IFNG','PRF1','GZMK','GZMA','GZMB','NKG7'),
  Treg = c('FOXP3','IL2RA','IKZF2'),
  #gammaDelta = c('TCRGC1','TCRGC2','TCRG-C1','TCRG-C2'),
  Activated = c('CD69','CD25'),
  Naive = c('CD45RA','LEF1','SELL','TCF7'),
  Proliferation = c('MKI67','PCNA','TOP2A'),
  #Memory = c('CD45RO','CD62L'),
  Effector = c('CCR7'),
  Exhaustion = c('PD1','TIM3','LAG3','TIGIT','TOX','TOX2','CTLA4','HAVCR2','PDCD1')
)

## plot
library(ggplot2)
library(RColorBrewer)

pdf('tcellMarkers.pdf', width = 25,height = 5,onefile = F)
DotPlot(tcell, assay = 'RNA',
        features = humanMarkers, group.by = "seurat_clusters",
        cols = c("#f8fcfb", "#2193b0"),
        col.min = 0)+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
dev.off()

### violin plot
names(humanMarkers)
lapply(names(humanMarkers), function(cellType){
  pdf(paste0(cellType,'_vlnPlot.pdf'),width = 15,height = 15,onefile = F)
  print(
    VlnPlot(tcell, assay = 'RNA',
            features = humanMarkers[[cellType]], pt.size = 0)
  )
  dev.off()
})

## add common cell types metaData ====
tcellType <- ifelse(tcell$seurat_clusters %in% c(0,1,4,6:8,12),'Tcyto',
                    ifelse(tcell$seurat_clusters %in% c(13),'PCTL', #Proliferative Cytotoxic T Lymphocytes
                           ifelse(tcell$seurat_clusters %in% c(2,3,10,11),'Th',
                                  ifelse(tcell$seurat_clusters %in% c(5),'ExhrTH','Teff')))) %>%
  factor(levels = c('Tcyto','PCTL','Th','ExhrTH','Teff'))

tcell <- AddMetaData(tcell,tcellType,col.name = 'tcellType')
table(tcell$tcellType)
table(tcell$tcellType,tcell$orig.ident)

## delete not good markers (common) ====
humanMarkers <- list(
  #Lymphoid
  Tcyto = c('CD8','CD8A','CD8B','IFNG','PRF1','GZMK','GZMA','GZMB','NKG7'),
  Th = c('CD4','CCR6'),
  Treg = c('FOXP3','IL2RA','IKZF2'),
  Proliferation = c('MKI67','PCNA','TOP2A'),
  #Memory = c('CD45RO','CD62L'),
  Effector = c('CCR7'),
  Exhaustion = c('TIGIT','TOX2','CTLA4')
)

### dot plot
pdf('tcellCellMarkers_del.pdf',width = 12,height = 3,onefile = F)
DotPlot(tcell, assay = 'RNA',
        features = humanMarkers,group.by = "tcellType",
        cols = c("#f8fcfb", "#f12711"),
        col.min = 0)+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
dev.off()

## plot
pdf('tcellCellType_UMAP.pdf',width = 13,height = 3.5,onefile = F)
DimPlot(tcell, reduction = "umap", label = TRUE, 
        group.by = c("orig.ident", "tcellType",'seurat_clusters'))
dev.off()
pdf('tcellCellType_tSNE.pdf',width = 13,height = 3.5,onefile = F)
DimPlot(tcell, reduction = "tsne", label = TRUE, 
        group.by = c("orig.ident", "tcellType",'seurat_clusters'))
dev.off()

## plot cell type percent of each patient ====
library(plyr)
library(ggplot2)

plotData <- table(tcell$Sample,tcell$tcellType) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Sample','Cluster','Freq')
plotData$Sample <- factor(plotData$Sample, levels = c('Naive', 'Treatment'))

plotData <- ddply(plotData, "Cluster", transform,
                  percent = Freq / sum(Freq) * 100)

pdf(paste0('tcellCell_Sample_Percent.pdf'),
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

pdf(paste0('tcellCell_Cluster_Percent.pdf')
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
        labs(x='orig.ident',y='Percent',fill='Cluster'))
dev.off() 

saveRDS(tcell, 'SCCE_tcell_reAnno.rds')

# find markers ====
library(ggplot2)

## DEGs of different cell types
Idents(tcell)="tcellType"
DFgenes <- FindAllMarkers(tcell,only.pos = TRUE)#
DFgenes <- DFgenes[DFgenes$p_val_adj<0.05,]
write.table(DFgenes,file=paste0("Tcelltype_DEG.csv"),
            sep = ",",col.names = T,row.names = T)

## DEGs of Naive VS Treatment
Idents(tcell)="Sample"
DFgenes=FindMarkers(tcell,ident.1 = "Naive",ident.2 = "Treatment",logfc.threshold = 0.01,min.pct = 0.01)
write.table(DFgenes,file=paste0("DEG_NvsT.csv"),
            sep = ",",col.names = T,row.names = T)

## volcano plot ====
volcanoFUN = function(dataset=NULL,title=NULL,sampleoutpath=NULL,sample=NULL,cut_off_logFC=NULL,
                      cut_off_pvalue=NULL,
                      labelUp=NULL,labelDown=NULL,w=4,h=4.6){
  library(ggrepel)
  # 设置pvalue和logFC的阈值
  cut_off_pvalue = cut_off_pvalue# 0.001
  cut_off_logFC = cut_off_logFC #0.5
  # 根据阈值分别为上调基因设置‘up’，下调基因设置‘Down’，无差异设置‘Stable’，保存到change列
  # 这里的change列用来设置火山图点的颜色
  # dataset = DFIN
  dataset = dataset[dataset$p_val_adj!=1,]
  dataset$gene = rownames(dataset)
  dataset$change = ifelse(dataset$p_val_adj < cut_off_pvalue & abs(dataset$avg_log2FC) >= cut_off_logFC, 
                          ifelse(dataset$avg_log2FC> cut_off_logFC ,labelUp,labelDown),
                          'NoSig') %>% factor(levels = c('SCCE_1N','NoSig','SCCE_1T'))
  ## 将-log10(p-value)>100的值均写为100 
  # dataset[dataset$p_val_adj<1E-20,]$p_val_adj = 1E-20
  # 绘制火山图
  pdf(paste0(sample,"Volcanoplot_",title,".pdf"),width = w,height = h)
  p = ggplot(
    #设置数据
    dataset,aes(x = avg_log2FC,
                y = -log10(p_val_adj),
                colour=change)) +
    geom_point(alpha=1, size=0.5) +
    scale_color_manual(values=c("#BC3C28","grey","#0072B5")
    )+ 
    
    # 辅助线
    geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6) +
    geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="grey",lwd=0.6) +
    
    # 坐标轴
    labs(x=expression("Log"["2"]*"(Fold change)"),
         y=expression("-Log"["10"]*"(adjust P vaule)"),
         title = title)+ 
    theme_bw()+
    
    # 图例
    theme_classic(base_size = 15) + 
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.text.align = 0,
          legend.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(size = 13),
          axis.title = element_text(size=15),
          legend.position = "top",
          legend.text = element_text(size = 13)) +
    guides(colour = guide_legend(override.aes = list(size=3),reverse = T))+ #可使得legend中的圆点变大
    
    geom_text_repel(aes(label=label), 
                    segment.color = "black",
                    show.legend = FALSE,
                    box.padding=unit(0.15, "lines"), 
                    point.padding=unit(0.5, "lines"), 
                    #color="white",
                    max.overlaps = 10000000000)
  print(p)
  dev.off()
}

load('../volcanoFUN.RData')
## plot M450 vs LC456 ====
DFgenes = read.table(file=paste0("DEG_NvsT.csv"),sep = ",")

genes.to.label = subset(DFgenes,p_val_adj < 0.01) %>% .$avg_log2FC %>% abs %>% 
  order %>% tail(30) %>% rownames(subset(DFgenes,p_val_adj < 0.01))[.]
#selected sig genes

DFgenes$label = ifelse((rownames(DFgenes) %in% genes.to.label), as.character(rownames(DFgenes)),"")

volcanoFUN(dataset = DFgenes[DFgenes$p_val_adj!=1,],
           title = "Naive vs Treatment",
           #sampleoutpath = sampleoutpath_volcano,
           cut_off_pvalue=0.05,
           cut_off_logFC=1,
           # sample = "subMac",
           labelUp="Naive",
           labelDown = "Treatment",w=8,h=6)

## plot cell types ====
volcanoFUN_multi = function(dataset=NULL,title=NULL,sampleoutpath=NULL,sample=NULL,
                            cut_off_logFC=NULL,cut_off_pvalue=NULL,
                            lab1=NULL,lab2=NULL,lab3=NULL,lab4=NULL,
                            colorVar=c("#BC3C28","grey","#0072B5"),w=4,h=4.6){
  library(ggrepel)
  # 设置pvalue和logFC的阈值
  cut_off_pvalue = cut_off_pvalue# 0.001
  cut_off_logFC = cut_off_logFC #0.5
  # 根据阈值分别为上调基因设置‘up’，下调基因设置‘Down’，无差异设置‘Stable’，保存到change列
  # 这里的change列用来设置火山图点的颜色
  # dataset = DFIN
  dataset = dataset[dataset$p_val_adj!=1,]
  dataset <- subset(DFgenes,p_val_adj < cut_off_pvalue & abs(avg_log2FC)>=cut_off_logFC)
  dataset$gene = rownames(dataset)
  dataset$change = factor(dataset$cluster,levels = c(lab1,lab2,lab3,lab4))
  ## 将-log10(p-value)>100的值均写为100 
  # dataset[dataset$p_val_adj<1E-20,]$p_val_adj = 1E-20
  # 绘制火山图
  pdf(paste0(sample,"Volcanoplot_",title,".pdf"),width = w,height = h)
  p = ggplot(
    #设置数据
    dataset,aes(x = avg_log2FC,
                y = -log10(p_val_adj),
                colour=change)) +
    geom_point(alpha=1, size=0.5) +
    scale_color_manual(values=colorVar
    )+ 
    
    # 坐标轴
    labs(x=expression("Log"["2"]*"(Fold change)"),
         y=expression("-Log"["10"]*"(adjust P vaule)"),
         title = title)+ 
    theme_bw()+
    
    # 图例
    theme_classic(base_size = 15) + 
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.text.align = 0,
          legend.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(size = 13),
          axis.title = element_text(size=15),
          legend.position = "top",
          legend.text = element_text(size = 13)) +
    guides(colour = guide_legend(override.aes = list(size=3),reverse = T))+ #可使得legend中的圆点变大
    
    geom_text_repel(aes(label=label), 
                    segment.color = "black",
                    show.legend = FALSE,
                    box.padding=unit(0.15, "lines"), 
                    point.padding=unit(0.5, "lines"), 
                    #color="white",
                    max.overlaps = 10000000000)
  print(p)
  dev.off()
}

load('../volcanoFUN_multi.RData')
DFgenes = read.table(file=paste0("Tcelltype_DEG.csv"),sep = ",")

#selected sig genes
genes.to.label = subset(DFgenes,p_val_adj < 0.01 & abs(avg_log2FC)>1) %>% .$avg_log2FC %>%
  order %>% subset(DFgenes,p_val_adj < 0.01 & avg_log2FC>1)[.,] %>% group_by(cluster) %>%
  slice_max(n=5,order_by = abs(avg_log2FC))

DFgenes$label = ifelse((rownames(DFgenes) %in% genes.to.label$gene), as.character(rownames(DFgenes)),"")

volcanoFUN_multi(dataset=DFgenes,
                 title='DEGs in T cell types',
                 #sampleoutpath=NULL,sample=NULL,
                 cut_off_logFC=1,
                 cut_off_pvalue=0.01,
                 lab1='Tcyto',
                 lab2='PCTL',
                 lab3='Th',
                 lab4='ExhrTH',
                 lab5='Teff',
                 colorVar=c('#ff7e67',"#0278ae","#e03e36","#40a798","#700961"),
                 w=8,h=6.5)