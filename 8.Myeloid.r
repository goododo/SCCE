rm(list = ls());gc()

if(!dir.exists('C:/Users/10784/Desktop/SCCE/8.1.Myeloid_DEG')){
  dir.create('C:/Users/10784/Desktop/SCCE/8.1.Myeloid_DEG',recursive = T)}

setwd('C:/Users/10784/Desktop/SCCE/8.1.Myeloid_DEG')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)

grep('.rds',list.files('../2.1.annoMain'),value = T)
# load data ====
Myeloid <- readRDS('../2.1.annoMain/SCCE_Myeloid.rds')
#24324 features across 828 samples

## normalize ====
Myeloid <- NormalizeData(object = Myeloid, assay = 'RNA',
                       normalization.method = "LogNormalize",
                       scale.factor = 10000,
                       margin = 1, #If performing CLR normalization, normalize across features (1) or cells (2)
                       verbose = TRUE) # 进行标准化，默认参数

## find variable features ====
Myeloid <- FindVariableFeatures(
  object = Myeloid,
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
Myeloid <- ScaleData(Myeloid, features = rownames(Myeloid))

## PCA降维 ====
Myeloid <- RunPCA(Myeloid, features = VariableFeatures(object = Myeloid))

## 确定数据集维度
Myeloid <- JackStraw(Myeloid, num.replicate = 100)
Myeloid <- ScoreJackStraw(Myeloid, dims = 1:20)
#ElbowPlot(Myeloid) # 20

## Cluster ====
Myeloid <- FindNeighbors(Myeloid, dims = 1:20)
Myeloid <- FindClusters(Myeloid, resolution = .1) %>% 
  FindClusters( resolution = .2) %>% FindClusters( resolution = .3) %>% 
  FindClusters( resolution = .4) %>% FindClusters( resolution = .5) %>% 
  FindClusters( resolution = .6) %>% FindClusters( resolution = .7) %>% 
  FindClusters( resolution = .8) %>% FindClusters( resolution = .9) %>% 
  FindClusters( resolution = 1) %>% FindClusters( resolution = 1.1) %>% 
  FindClusters( resolution = 1.2) %>% FindClusters( resolution = 1.3) %>% 
  FindClusters( resolution = 1.4)

## plot
library(clustree)
pdf('Myeloid_clustree.pdf',width = 15,height = 11,onefile = F)
clustree(Myeloid)
dev.off()

# 根据markers热图选resolution ====
lapply(seq(0.1,0.7,by=0.1), function(reso){
  
  gc()
  ## set resolution
  Myeloid <- FindClusters(Myeloid, resolution = reso)
  
  # UMAP & tSNE 降维 ====
  Myeloid <- RunUMAP(Myeloid, dims = 1:20)
  Myeloid <- RunTSNE(Myeloid, dims = 1:20)
  
  ## plot UMAP & t-SNE
  pdf(paste0('reso=',reso,'_UMAP.pdf'),width = 10,height = 4,onefile = F)
  print(
    DimPlot(Myeloid, reduction = "umap", label = TRUE, 
            group.by = c("seurat_clusters", "orig.ident"))
  )
  dev.off()
  pdf(paste0('reso=',reso,'_tSNE.pdf'),width = 10,height = 4,onefile = F)
  print(
    DimPlot(Myeloid, reduction = "tsne", label = TRUE, 
            group.by = c("seurat_clusters", "orig.ident"))
  )
  dev.off()
  
  # find markers ====
  allmarkers <- FindAllMarkers(
    object = Myeloid, assay = 'integrated',
    only.pos= TRUE,min.pct = .25, logfc.threshold = 0.25) %>% 
    group_by(cluster)
  
  ## plot heatmap
  gc()
  pdf(paste0('reso=',reso,'_markers_heatmap.pdf'),onefile = F)
  print(DoHeatmap(Myeloid, features = allmarkers$gene, cells = NULL,
                  group.by = "seurat_clusters", group.bar = TRUE,
                  group.colors = NULL, disp.min = -2.5, disp.max = NULL,
                  slot = "scale.data", assay = 'integrated', label = TRUE,
                  size = 5.5, hjust = 0, angle = 45, raster = TRUE,
                  draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02,
                  combine = TRUE ))
  dev.off()
})

Myeloid <- FindClusters(Myeloid, resolution = .6)
# Look at cluster IDs of the first 5 cells
head(Idents(Myeloid), 5)
table(Myeloid$seurat_clusters)

# UMAP & tSNE 降维 ====
Myeloid <- RunUMAP(Myeloid, dims = 1:20)
Myeloid <- RunTSNE(Myeloid, dims = 1:20)

# cell annotation ====
# https://www.cellsignal.com/pathways/Myeloid-cell-markers-human
# https://www.abcam.com/primary-antibodies/Myeloid-cell-markers-poster
# https://lishensuo.github.io/posts/bioinfo/028%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--%E5%9F%BA%E4%BA%8E%E6%96%87%E7%8C%AE%E7%9A%84%E7%BB%86%E8%83%9E%E7%B1%BB%E5%9E%8B%E6%B3%A8%E9%87%8Amarker/
library(Hmisc)
humanMarkers <- list(
  Mono = c('CD14', 'CD16', 'CD64', 'CCR2'),
  Macro = c('CD68', 'CD163', 'CD206', 'CD11B'),
  DC = c('CD11C', 'CD123', 'CD83', 'CD86', 'CD40', 'HLADR'),
  #Neutro = c('CD66B', 'CEACAM8', 'CD15'),
  Eosino = c('CD125', 'CD9', 'SIGLEC8')
)

## plot
library(ggplot2)
library(RColorBrewer)

pdf('MyeloidMarkers.pdf', width = 12,height = 4.5,onefile = F)
DotPlot(Myeloid, assay = 'RNA',
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
    VlnPlot(Myeloid, assay = 'RNA',
            features = humanMarkers[[cellType]], pt.size = 0)
  )
  dev.off()
})

# find markers ====
library(ggplot2)

## DEGs of different cell types
#Idents(Myeloid)="MyeloidType"
DFgenes <- FindAllMarkers(Myeloid,only.pos = TRUE)#
DFgenes <- DFgenes[DFgenes$p_val_adj<0.05,]
write.table(DFgenes,file=paste0("Myeloidtype_DEG.csv"),
            sep = ",",col.names = T,row.names = T)

## DEGs of Naive VS Treatment
Idents(Myeloid)="Sample"
DFgenes=FindMarkers(Myeloid,ident.1 = "Naive",ident.2 = "Treatment",logfc.threshold = 0.01,min.pct = 0.01)
write.table(DFgenes,file=paste0("DEG_NvsT.csv"),
            sep = ",",col.names = T,row.names = T)

## volcano plot ====
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
load('../volcanoFUN_multi.RData')
DFgenes = read.table(file=paste0("Myeloidtype_DEG.csv"),sep = ",")

#selected sig genes
genes.to.label = subset(DFgenes,p_val_adj < 0.01 & abs(avg_log2FC)>1) %>% .$avg_log2FC %>%
  order %>% subset(DFgenes,p_val_adj < 0.01 & avg_log2FC>1)[.,] %>% group_by(cluster) %>%
  slice_max(n=5,order_by = abs(avg_log2FC))

DFgenes$label = ifelse((rownames(DFgenes) %in% genes.to.label$gene), as.character(rownames(DFgenes)),"")

volcanoFUN_multi(dataset=DFgenes,
                 title='DEGs in Myeloid cell types',
                 #sampleoutpath=NULL,sample=NULL,
                 cut_off_logFC=1,
                 cut_off_pvalue=0.01,
                 lab1='0',
                 lab2='1',
                 lab3='2',
                 lab4='3',
                 lab5='4',
                 lab6 = '5',
                 lab7 = '6',
                 lab8 = '7',
                 colorVar=c('#ff7e67',"#0278ae","#e03e36","#40a798","#700961",'#a5dee5','#e0f9b5','#ffcfdf'),
                 w=8,h=6.5)

# DEG ====
if(!dir.exists('C:\\Users\\10784\\Desktop\\SCCE\\8.2.Myeloid_func')){
  dir.create('C:\\Users\\10784\\Desktop\\SCCE\\8.2.Myeloid_func',recursive = T)}

setwd('C:\\Users\\10784\\Desktop\\SCCE\\8.2.Myeloid_func')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)

# functional enrichment ====
library(clusterProfiler)
library(org.Hs.eg.db)

# Naive vs Treatment ====
DFgenes = read.table("../8.1.Myeloid_DEG/DEG_NvsT.csv",sep = ",")
DFgenes$gene <- rownames(DFgenes)
ids <- bitr(DFgenes$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
DFgenes <- merge(DFgenes,ids,by.x='gene',by.y='SYMBOL')

DFgenes.sig <- subset(DFgenes,p_val < 0.05)

geneList <- DFgenes$avg_log2FC
names(geneList) <- rownames(DFgenes)
geneList <- sort(geneList,decreasing = T)

## KEGG ====
### enrich KEGG ====
ekegg <- enrichKEGG(gene = DFgenes.sig$ENTREZID,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05,
                    minGSSize = 3,
                    maxGSSize = 500,
                    qvalueCutoff = 0.2)

write.table(as.data.frame(ekegg),'allRes_enrichKEGG_NvsT.txt',sep = '\t',quote = F,row.names = F)

#### plot ====
plotData <- ekegg@result[
  which(ekegg@result$ID %in% c(
    'hsa04060','hsa04640','hsa04110','hsa04512',
    'hsa04514','hsa04151','hsa04659','hsa04115',
    'hsa04510','hsa04010')),]
plotData <- plotData[order(plotData$subcategory,plotData$Description),]
plotData$subcategory <- factor(plotData$subcategory,levels = unique(plotData$subcategory) %>% sort)
plotData$Description <- factor(plotData$Description,levels = plotData$Description)

pdf('enrich KEGG Res.pdf',width = 8,height = 4,onefile = F)
ggplot(plotData,
       aes(x=-log10(pvalue),y=Description, fill=subcategory)) + #x、y轴定义；根据type填充颜色
  
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  #scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666",'blue') ) +  #柱状图填充颜色
  
  #coord_flip() +  #让柱状图变为纵向
  
  xlab("-log10(P value)") +  #x轴标签
  ylab("KEGG pathway") +  #y轴标签
  labs(title = "KEGG Pathway")+  #设置标题
  
  theme_bw()
dev.off()

### gse KEGG ====
gsekegg <- gseKEGG(geneList = geneList,
                   organism = 'hsa',
                   minGSSize = 3,
                   pvalueCutoff = 1,
                   verbose = T)

write.table(as.data.frame(gsekegg),'allRes_gseKEGG_NvsT.txt',sep = '\t',quote = F,row.names = F)

## GO ====
### enrich GO ====
ego <- enrichGO(gene = DFgenes.sig$ENTREZID,
                universe = names(geneList),
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                #ont  = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable  = T
)
write.table(as.data.frame(ego@result),'allRes_enrichGO_NvsT.txt',sep = '\t',quote = F,row.names = F)

## gse GO ====
gsego <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               #ont          = "CC",
               minGSSize    = 3,
               maxGSSize    = 500,
               pvalueCutoff = 1,
               verbose      = T)
write.table(as.data.frame(gsego),'allRes_gseGO_NvsT.txt',sep = '\t',quote = F,row.names = F)

# clusters ====
DFgenes = read.table("../8.1.Myeloid_DEG/Myeloid_cluster_DEG.csv",sep = ",",header = T,row.names = 1)
ids <- bitr(DFgenes$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
DFgenes <- merge(DFgenes,ids,by.x='gene',by.y='SYMBOL')

DFgenes.sig <- subset(DFgenes,p_val < 0.05)

geneList <- DFgenes$avg_log2FC
names(geneList) <- rownames(DFgenes)
geneList <- sort(geneList,decreasing = T)

## KEGG ====
compare_kegg <- compareCluster(ENTREZID~cluster,
                               data = DFgenes.sig,
                               fun = "enrichKEGG",
                               organism = "hsa", pvalueCutoff = 0.05)

write.table(as.data.frame(compare_kegg),'allRes_enrichKEGG_cluster.txt',sep = '\t',quote = F,row.names = F)

#### plot ====
compare_kegg <- read.table('selectRes_enrichKEGG_cluster.txt',sep = '\t',header = T)

lapply(0:7, function(clu){
  table(subset(compare_kegg,Cluster==clu)$Cluster,subset(compare_kegg,Cluster==clu)$subcategory)
})

plotData <- compare_kegg
plotData <- plotData[order(plotData$Cluster,plotData$subcategory,plotData$Description),]
plotData$cluster <- factor(plotData$cluster,levels = unique(plotData$cluster) %>% sort)
plotData$subcategory <- factor(plotData$subcategory,levels = unique(plotData$subcategory) %>% sort)
plotData$Description <- factor(plotData$Description,levels =  unique(plotData$Description))

pdf('enrich KEGG Cluster Res.pdf',width = 8,height = 4.5,onefile = F)
ggplot(plotData,
       aes(x=cluster,y=Description, color=-log10(p.adjust),size=Count)) + #x、y轴定义；根据type填充颜色
  
  geom_point() +  
  scale_color_gradient(low = "#ed85b0", high = "#8d1f17")+
  
  #coord_flip() +  #让柱状图变为纵向
  
  xlab("-log10(adjust P value)") +  #x轴标签
  ylab("KEGG pathway") +  #y轴标签
  labs(title = "KEGG Pathway")+  #设置标题
  
  theme_bw()
dev.off()

## GO ====
compare_go <- compareCluster(ENTREZID~cluster,
                             data = DFgenes.sig,
                             fun = "enrichGO",
                             OrgDb = "org.Hs.eg.db",
                             ont = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.01,
                             qvalueCutoff = 0.05)

write.table(as.data.frame(compare_go),'allRes_enrichGO_cluster.txt',sep = '\t',quote = F,row.names = F)

## add common cell types metaData ====
MyeloidType <- ifelse(Myeloid$seurat_clusters %in% c(0),'Mye-Th Collaborator',
                      ifelse(Myeloid$seurat_clusters %in% c(1),'RepliRepair Myeloid',
                             ifelse(Myeloid$seurat_clusters %in% c(2),'Prolif Myeloid',
                                    ifelse(Myeloid$seurat_clusters %in% c(3),'ImmuneControl Myeloid',
                                           ifelse(Myeloid$seurat_clusters %in% c(4),'Onco Myeloid',
                                                  ifelse(Myeloid$seurat_clusters %in% c(5),'Signaling Myeloid',
                                                         ifelse(Myeloid$seurat_clusters %in% c(6),'ProlifReg Myeloid',
                                                                'CytokineReg Myeloid'))))))) %>%
  factor(levels = c('Mye-Th Collaborator','RepliRepair Myeloid','Prolif Myeloid','ImmuneControl Myeloid',
                    'Onco Myeloid','Signaling Myeloid','ProlifReg Myeloid','CytokineReg Myeloid'))

Myeloid <- AddMetaData(Myeloid,MyeloidType,col.name = 'MyeloidType')
table(Myeloid$MyeloidType)
table(Myeloid$MyeloidType,Myeloid$orig.ident)

## plot
pdf('MyeloidCellType_UMAP.pdf',width = 14,height = 3.7,onefile = F)
DimPlot(Myeloid, reduction = "umap", label = TRUE, 
        group.by = c("Sample", "MyeloidType",'seurat_clusters'))
dev.off()
pdf('MyeloidCellType_tSNE.pdf',width = 14,height = 3.7,onefile = F)
DimPlot(Myeloid, reduction = "tsne", label = TRUE, 
        group.by = c("Sample", "MyeloidType",'seurat_clusters'))
dev.off()

## plot cell type percent of each patient ====
library(plyr)
library(ggplot2)

plotData <- table(Myeloid$Sample,Myeloid$MyeloidType) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Sample','Cluster','Freq')
plotData$Sample <- factor(plotData$Sample, levels = c('Naive', 'Treatment'))

plotData <- ddply(plotData, "Cluster", transform,
                  percent = Freq / sum(Freq) * 100)

pdf(paste0('MyeloidCell_Sample_Percent.pdf'),
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
        theme_bw()+labs(x='Cluster',y='Percent',fill='Sample'))+
  coord_flip()
dev.off() 

plotData <- ddply(plotData, "Sample", transform,
                  percent = Freq / sum(Freq) * 100)

pdf(paste0('MyeloidCell_Cluster_Percent.pdf')
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

saveRDS(Myeloid, 'SCCE_Myeloid_reAnno.rds')
