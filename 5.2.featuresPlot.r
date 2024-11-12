rm(list = ls());gc()

if(!dir.exists('C:/Users/10784/Desktop/SCCE/5.2.featurePlot')){
  dir.create('C:/Users/10784/Desktop/SCCE/5.2.featurePlot',recursive = T)}

setwd('C:/Users/10784/Desktop/SCCE/5.2.featurePlot')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)

# load data ====
epi <- readRDS('../5.1.pseudo/Seurat/epi_withPseudotime.RDS')

# t-SNE & UMAP ====
pdf('epiCell_UMAP.pdf',width = 12.6,height = 3.5,onefile = F)
DimPlot(epi, reduction = "umap", label = TRUE, 
        group.by = c("Sample",'pseudoType','TumorType'))
dev.off()
pdf('epiCell_tSNE.pdf',width = 12.6,height = 3.5,onefile = F)
DimPlot(epi, reduction = "tsne", label = TRUE, 
        group.by = c("Sample",'pseudoType','TumorType'))
dev.off()

# cell ratio ====
library(plyr)
library(ggplot2)

plotData <- table(epi$Sample,epi$pseudoType) %>% data.frame()# %>% 
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

# tumor suppressor genes ====
#https://doctor.get.com.tw/m/Journal/detail.aspx?no=402462
tsg <- c(
  'p53','TP53','RB1','p16','NK4A','APC','p14arf','CDKN2A','WT1','NF1','NF2',
  'VHL','BRCA1','BRCA2','MEN1','LKB1','STK1','PTCH','PTEN','MMAC1','DPC4',
  'ECAD','EXT1','EXT2','TSC1','TSC2','MSH2','MLH1','PMS1','PMS2','MSH6',
  'APC','PTCH1','ATM','STK11','SMAD4','WWOX','PML','CDH1','HNF1A','BAP1','RUNX3','DCC'
) %>% unique %>% toupper

tsg2 <- tsg[which(tsg %in% rownames(epi))]

# feature plot ====
pdf('tumorSupGene_featurePlot.pdf', width = 17,height = 30,onefile = F)
FeaturePlot(epi, tsg, cols = c("#c6ffdd", "#fbd786", "#f7797d"),reduction = 'tsne')
dev.off() 

## dark theme (UGLY)
baseplot <- FeaturePlot(epi, tsg, cols = c("#c6ffdd", "#fbd786", "#f7797d"),reduction = 'umap',
                        pt.size = 0.1, combine = FALSE)
for (i in 1:length(x = baseplot)) {
  baseplot[[i]] <- baseplot[[i]] + DarkTheme()
}
# Can also use lapply
baseplot <- lapply(
  X = baseplot,
  FUN = function(p) {
    return(p + DarkTheme())
  }
)

pdf('tumorSupGene_featurePlotdark.pdf', width = 17,height = 34,onefile = F)
CombinePlots(plots = baseplot, ncol = 5)
dev.off()

# violin plot ====
# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(hrbrthemes)
library(viridis)

# load data ====
naive <- readRDS('../0.rawdata/SCCE_1N_QC_SeuratObj.rds')
treat <- readRDS('../0.rawdata/SCCE_1T_QC_SeuratObj.rds')

epi <- readRDS('../5.1.pseudo/Seurat/epi_withPseudotime.RDS')
extraMeta.epi <- epi@meta.data[,c(5,8,22:24)]

naive.epi <- naive[, which(colnames(naive) %in% (strsplit(colnames(epi),'_1|_2') %>% unlist %>% unique))]
naive.epi <- AddMetaData(
  object = naive.epi,
  metadata = extraMeta.epi[which(colnames(naive.epi) %in% (
    strsplit(colnames(epi),'_1|_2') %>% unlist %>% unique)),],
  col.name = NULL
)

treat.epi <- treat[, which(colnames(treat) %in% (strsplit(colnames(epi),'_1|_2') %>% unlist %>% unique))]
treat.epi <- AddMetaData(
  object = treat.epi,
  metadata = extraMeta.epi[which(colnames(treat.epi) %in% (
    strsplit(colnames(epi),'_1|_2') %>% unlist %>% unique)),],
  col.name = NULL
)

naive.epi.count <- naive.epi@assays[["RNA"]]@layers[["counts"]]
dimnames(naive.epi.count) <- list(rownames(naive.epi),colnames(naive.epi))

treat.epi.count <- treat.epi@assays[["RNA"]]@layers[["counts"]]
dimnames(treat.epi.count) <- list(rownames(treat.epi),colnames(treat.epi))

plotData <- cbind(naive.epi.count[intersect(rownames(naive.epi),rownames(treat.epi)) %>% intersect(tsg),],
                  treat.epi.count[intersect(rownames(naive.epi),rownames(treat.epi)) %>% intersect(tsg),]) %>% 
  as.data.frame %>% t %>% as.data.frame

group <- c(naive.epi@meta.data[["pseudoType"]],treat.epi@meta.data[["pseudoType"]])

plotData <- cbind(Group=as.character(group), plotData)
plotData$Group <- factor(plotData$Group, levels = levels(group))

# Plot
plotData %>%
  ggplot( aes(x=text, y=value, fill=text, color=text)) +
  geom_violin(width=2.1, size=0.2) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none"
  ) +
  coord_flip() + # This switch X and Y axis and allows to get the horizontal version
  xlab("") +
  ylab("Assigned Probability (%)")

plotList <- lapply(c(2:ncol(plotData)), function(x){
  #range01 <- function(exp){(exp-min(exp))/(max(exp)-min(exp))}
  data <- plotData[,c(1,x)]
  data[,2] <- scale(data[,2])
  gene <- colnames(data)[2]
  colnames(data)[2] <- 'Gene'
  ggplot(data, aes(x=Group, y=Gene, fill=Group, color=Group)) +
    geom_violin(width=1.1, size=0.2) +
    scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) +
    theme_ipsum() +
    theme(
      legend.position="none"
    ) +
    coord_flip() + # This switch X and Y axis and allows to get the horizontal version
    xlab("") +
    ylab("Expression Level")+
    ggtitle(gene)
})

library(grid)
library(gridExtra)
library(showtext)
showtext_auto()
pdf('tumorSupGene_vlnPlot.pdf', width = 25,height = 25,onefile = F)
grid.arrange(plotList[[1]], plotList[[2]], plotList[[3]], plotList[[4]], plotList[[5]],
             plotList[[6]], plotList[[7]], plotList[[8]], plotList[[9]], plotList[[10]],
             plotList[[11]], plotList[[12]], plotList[[13]], plotList[[14]], plotList[[15]],
             plotList[[16]], plotList[[17]], plotList[[18]], plotList[[19]], plotList[[20]],
             plotList[[21]], plotList[[22]], plotList[[23]], plotList[[24]], plotList[[25]],
             plotList[[26]], plotList[[27]], plotList[[28]], plotList[[29]], 
             nrow = 6)
dev.off()

# dot plot ====
pdf('tumorSupGene_dotPlot.pdf', width = 8,height = 2.5,onefile = F)
DotPlot(epi, tsg2, cols = c("#ffefba",'#ff0080'),
        col.min = -2.5,col.max = .5, scale = F,
        group.by = 'pseudoType',
        cluster.idents = F)
dev.off()

# heatmap ====
pdf('tumorSupGene_heatmap.pdf', width = 8,height = 3.6,onefile = F)
DoHeatmap(epi, features = tsg, cells = NULL,
          group.by = "pseudoType", group.bar = TRUE,
          group.colors = NULL,
          disp.min = -1, disp.max = 1,
          slot = "scale.data", assay = 'integrated', label = TRUE,
          size = 3, hjust = 0, angle = 45, raster = F,
          draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.01,
          combine = TRUE ) + 
  #scale_fill_gradientn(colors = c("#619DB8","#AECDD7","#E3EEEF",'#FAE7D9','#F0B79A','#C85D4D')))
  scale_fill_gradientn(colors = c("#619DB8","#AECDD7","#E3EEEF",'#F0B79A','#C85D4D'))
dev.off()