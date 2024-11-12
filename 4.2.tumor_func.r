rm(list = ls());gc()

if(!dir.exists('SCCE\\4.2.tumor_func')){
  dir.create('SCCE\\4.2.tumor_func',recursive = T)}

setwd('SCCE\\4.2.tumor_func')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)

# functional enrichment ====
library(clusterProfiler)
library(org.Hs.eg.db)

#Treatment vs Naive ====
DFgenes = read.table("../4.1.tumor_DEG/DEG_TreatmentvsNaive.csv",sep = ",")
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

write.table(as.data.frame(ekegg),'allRes_enrichKEGG_TvsN.txt',sep = '\t',quote = F,row.names = F)

#### plot ====
plotData <- ekegg@result[
  which(ekegg@result$ID %in% c(
    'hsa04060','hsa04514','hsa04512','hsa04151',
    'hsa04110','hsa04510','hsa04015','hsa05205',
    'hsa04010','hsa04657','hsa04360','hsa04115',
    'hsa04014','hsa05202','hsa04668')),]
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

write.table(as.data.frame(gsekegg),'allRes_gseKEGG_TvsN.txt',sep = '\t',quote = F,row.names = F)

#### plot ====
plotData <- gsekegg@result[
  which(ekegg@result$ID %in% c(
    'hsa00983','hsa04662','hsa04620','hsa04622',
    'hsa01232','hsa00980','hsa05204','hsa04657')),]
plotData <- plotData[order(plotData$Description),]
plotData$Description <- factor(plotData$Description,levels = plotData$Description)

pdf('gse KEGG Res.pdf',width = 5,height = 3,onefile = F)
ggplot(plotData,
       aes(x=setSize,y=Description, fill=-log10(pvalue))) + #x、y轴定义；根据type填充颜色
  
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  #scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666",'blue') ) +  #柱状图填充颜色
  scale_fill_gradient(low = "#f27121",high = "#8a2387")+
  
  #coord_flip() +  #让柱状图变为纵向
  
  xlab("Set Size") +  #x轴标签
  ylab("KEGG pathway") +  #y轴标签
  labs(title = "KEGG Pathway")+  #设置标题
  
  theme_bw()
dev.off()

## GO ====
### enrich GO ====
ego <- enrichGO(gene = DFgenes.sig$ENTREZID,
                universe = names(geneList),
                keyType = 'ENTREZID',
                OrgDb = 'org.Hs.eg.db',
                #ont  = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable  = T
)
write.table(as.data.frame(ego@result),'allRes_enrichGO_TvsN.txt',sep = '\t',quote = F,row.names = F)

#### plot ====
ego <- read.table('allRes_enrichGO_TvsN.txt',sep = '\t',header = T)

plotData <- ego[
  which(ego$ID %in% c(
    'GO:0050840','GO:0005201','GO:0050786','GO:0001968','GO:0019838',
    'GO:0005178','GO:0005518','GO:0098641','GO:0098631','GO:0023023')),]
plotData <- plotData[order(plotData$Description),]
plotData$Description <- factor(plotData$Description,levels = plotData$Description)

pdf('enrich GO Res.pdf',width = 5,height = 3,onefile = F)
ggplot(plotData,
       aes(x=Count,y=Description, fill=-log10(pvalue))) + #x、y轴定义；根据type填充颜色
  
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  #scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666",'blue') ) +  #柱状图填充颜色
  scale_fill_gradient(low = "#f27121",high = "#8a2387")+
  
  #coord_flip() +  #让柱状图变为纵向
  
  xlab("Set Size") +  #x轴标签
  ylab("") +  #y轴标签
  labs(title = "GO Terms")+  #设置标题
  
  theme_bw()
dev.off()

## gse GO ====
gsego <- gseGO(geneList     = geneList,
               OrgDb        = 'org.Hs.eg.db',
               #ont          = "CC",
               minGSSize    = 3,
               maxGSSize    = 500,
               pvalueCutoff = 1,
               verbose      = T)
write.table(as.data.frame(gsego),'allRes_gseGO_TvsN.txt',sep = '\t',quote = F,row.names = F)

#### plot ====
gsego <- read.table('allRes_gseGO_TvsN.txt',sep = '\t',header = T)

plotData <- gsego[
  which(gsego$ID %in% c(
    'GO:2000278','GO:0071897','GO:0006259','GO:0009123','GO:0009161',
    'GO:0032392','GO:0019373','GO:0006177','GO:0034728','GO:0065004',
    'GO:0051091','GO:2000573')),]
plotData <- plotData[order(plotData$Description),]
plotData$Description <- factor(plotData$Description,levels = plotData$Description)

pdf('gse GO Res.pdf',width = 7.6,height = 3,onefile = F)
ggplot(plotData,
       aes(x=setSize,y=Description, fill=-log10(pvalue))) + #x、y轴定义；根据type填充颜色
  
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  #scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666",'blue') ) +  #柱状图填充颜色
  scale_fill_gradient(low = "#f27121",high = "#8a2387")+
  
  #coord_flip() +  #让柱状图变为纵向
  
  xlab("Set Size") +  #x轴标签
  ylab("") +  #y轴标签
  labs(title = "GO Terms")+  #设置标题
  
  theme_bw()
dev.off()

# clusters ====
DFgenes = read.table("../4.1.tumor_DEG/DEG_cluster.csv",sep = ",")
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

plotData <- compare_kegg
plotData <- plotData[order(plotData$Cluster,plotData$subcategory,plotData$Description),]
plotData$cluster <- factor(plotData$cluster,levels = unique(plotData$cluster) %>% sort)
plotData$subcategory <- factor(plotData$subcategory,levels = unique(plotData$subcategory) %>% sort)
plotData$Description <- factor(plotData$Description,levels =  unique(plotData$Description))

pdf('enrich KEGG Cluster Res.pdf',width = 8,height = 4,onefile = F)
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
