rm(list = ls());gc()

if(!dir.exists('SCCE\\4.1.tumor_DEG')){
  dir.create('SCCE\\4.1.tumor_DEG',recursive = T)}

setwd('SCCE\\4.1.tumor_DEG')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

# load data ====
epi <- readRDS('../3.1.epi_dimRedu/SCCE_Epithelial_reCluster.rds')

# find markers ====
library(ggplot2)

## DEGs of Treatment vs Naive ====
Idents(epi)="Sample"
DFgenes=FindMarkers(epi,ident.1 = "Treatment",ident.2 = "Naive",
                    logfc.threshold = 0.01,min.pct = 0.01, min.cells.group = 5,
                    slot = "scale.data", assay = 'integrated')
write.table(DFgenes,file=paste0("DEG_TreatmentvsNaive.csv"),
            sep = ",",col.names = T,row.names = T)

## DEGs of different cell clusters ====
Idents(epi)="seurat_clusters"
DFgenes <- FindAllMarkers(epi,only.pos = TRUE, min.cells.group = 5,slot = "scale.data", assay = 'integrated')
DFgenes <- DFgenes[DFgenes$p_val_adj<0.05,]
write.table(DFgenes,file=paste0("DEG_cluster.csv"),
            sep = ",",col.names = T,row.names = T)

# heatmap ====
DFgenes <- read.table("DEG_TreatmentvsNaive.csv",sep = ",")

genes.to.label <- subset(DFgenes,p_val_adj < 0.01) %>% .$avg_log2FC
genes.to.label <-  order(genes.to.label, decreasing = T) %>% 
  .[c(1:50,(length(genes.to.label)-49):length(genes.to.label))] %>% 
  rownames(subset(DFgenes,p_val_adj < 0.01))[.]#selected sig genes

pdf(paste0('heatmap_DEG_TreatmentvsNaive.pdf'),width = 6.5,height = 12.5,onefile = F)
print(DoHeatmap(epi, features = genes.to.label, cells = NULL,
                group.by = "Sample", group.bar = TRUE,
                group.colors = NULL, disp.min = -2.5, disp.max = NULL,
                label = TRUE, slot = "scale.data", assay = 'integrated', 
                size = 5.5, hjust = 0, angle = 45, raster = TRUE,
                draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02,
                combine = TRUE ) + 
        #scale_fill_gradientn(colors = c("#619DB8","#AECDD7","#E3EEEF",'#FAE7D9','#F0B79A','#C85D4D')))
        scale_fill_gradientn(colors = c("#619DB8","#E3EEEF",'#C85D4D')))
dev.off()

## DEGs of different cell clusters ====
### all ====
DFgenes <- read.table("DEG_cluster.csv",sep = ",")

genes.to.label = subset(DFgenes,p_val_adj < 0.01 & abs(avg_log2FC)>1) %>% .$avg_log2FC %>%
  order %>% subset(DFgenes,p_val_adj < 0.01 & abs(avg_log2FC)>1)[.,] %>% group_by(cluster) %>%
  slice_max(n=5,order_by = abs(avg_log2FC))
genes.to.label = subset(DFgenes,p_val_adj < 0.01) %>% .$avg_log2FC %>%
  order %>% subset(DFgenes,p_val_adj < 0.01)[.,] %>% group_by(cluster) %>%
  slice_max(n=5,order_by = abs(avg_log2FC))

pdf(paste0('heatmap_DEG_cluster.pdf'),width = 12,height = 9,onefile = F)
print(DoHeatmap(epi, features = genes.to.label$gene, cells = NULL,
                group.by = "seurat_clusters", group.bar = TRUE,
                group.colors = NULL,
                disp.min = -0, disp.max = 1,
                label = TRUE, slot = "scale.data", assay = 'integrated', 
                size = 5.5, hjust = 0, angle = 45, raster = TRUE,
                draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02,
                combine = TRUE ) + 
        #scale_fill_gradientn(colors = c("#619DB8","#AECDD7","#E3EEEF",'#FAE7D9','#F0B79A','#C85D4D')))
        scale_fill_gradientn(colors = c("#E3EEEF",'#F0B79A','#C85D4D')))
dev.off()

# volcano plot ====
load('../volcanoFUN.RData')
load('../volcanoFUN_multi.RData')

## DEGs of Treatment vs Naive ====
DFgenes <- read.table("DEG_TreatmentvsNaive.csv",sep = ",") %>% subset(p_val_adj!=1)

genes.to.label <- subset(DFgenes,p_val_adj < 0.01) %>% .$avg_log2FC
genes.to.label <-  order(genes.to.label, decreasing = T) %>% 
  .[c(1:10,(length(genes.to.label)-9):length(genes.to.label))] %>% 
  rownames(subset(DFgenes,p_val_adj < 0.01))[.]#selected sig genes

DFgenes$label = ifelse((rownames(DFgenes) %in% genes.to.label), as.character(rownames(DFgenes)),"")

volcanoFUN(dataset = DFgenes[DFgenes$p_val_adj!=1,],
           title = "Treatment vs Naive",
           #sampleoutpath = sampleoutpath_volcano,
           cut_off_pvalue=0.05,
           cut_off_logFC=1,
           #sample = "subMac",
           labelUp="Treatment",
           labelDown = "Naive",w=6,h=5)
