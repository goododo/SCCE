rm(list = ls());gc()

if(!dir.exists('C:/Users/10784/Desktop/SCCE/10.1.keyGenes_featurePlot')){
  dir.create('C:/Users/10784/Desktop/SCCE/10.1.keyGenes_featurePlot',recursive = T)}

setwd('C:/Users/10784/Desktop/SCCE/10.1.keyGenes_featurePlot')
# library packages ====
keyG <- c("MIF", "ACKR3", "CD74", "CXCR4", "CD44", "MDK", "SDC1",
          "SDC4", "PTPRZ1", "NCL", "ITGA6", "ITGB1", "FN1", "ITGA4",
          "ITGB1", "ITGAV", "ITGB6", "PTN", "HLA-E", "HLA-F", "HLA-C",
          "CD8A", "TIGIT", "PVR", "NECTIN2", "DSG1", "DSC3", "DSC2",
          "HBEGF", "EGFR", "AREG", "LIFR", "IL6ST", "OSM", "MPZL1",
          "LGALS9", "IGF1", "IGF1R", "CD99", "PILRA", "ANXA1", "FPR1")

# load data ====
naive <- readRDS('../0.rawdata/SCCE_1N_QC_SeuratObj.rds')
treat <- readRDS('../0.rawdata/SCCE_1T_QC_SeuratObj.rds')

epi <- readRDS('../5.1.pseudo/Seurat/epi_withPseudotime.RDS')
extraMeta.epi <- epi@meta.data[,c(5,8,22:24)]

myeloid <- readRDS('../8.2.Myeloid_func/SCCE_Myeloid_reAnno.rds')
extraMeta.mye <- myeloid@meta.data[,c(5,8,10)]
extraMeta.mye$MyeloidType <- gsub('-','',extraMeta.mye$MyeloidType)

tcell <- readRDS('../6.1.T_DEG/SCCE_tcell_reAnno.rds')
extraMeta.t <- tcell@meta.data[,c(5,8,22)]

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

naive.mye <- naive[, which(colnames(naive) %in% (strsplit(colnames(myeloid),'_1|_2') %>% unlist %>% unique))]
naive.mye <- AddMetaData(
  object = naive.mye,
  metadata = extraMeta.mye[which(colnames(naive.mye) %in% (
    strsplit(colnames(myeloid),'_1|_2') %>% unlist %>% unique)),],
  col.name = NULL
)

treat.mye <- treat[, which(colnames(treat) %in% (strsplit(colnames(myeloid),'_1|_2') %>% unlist %>% unique))]
treat.mye <- AddMetaData(
  object = treat.mye,
  metadata = extraMeta.mye[which(colnames(treat.mye) %in% (
    strsplit(colnames(myeloid),'_1|_2') %>% unlist %>% unique)),],
  col.name = NULL
)

naive.t <- naive[, which(colnames(naive) %in% (strsplit(colnames(tcell),'_1|_2') %>% unlist %>% unique))]
naive.t <- AddMetaData(
  object = naive.t,
  metadata = extraMeta.t[which(colnames(naive.t) %in% (
    strsplit(colnames(tcell),'_1|_2') %>% unlist %>% unique)),],
  col.name = NULL
)

treat.t <- treat[, which(colnames(treat) %in% (strsplit(colnames(tcell),'_1|_2') %>% unlist %>% unique))]
treat.t <- AddMetaData(
  object = treat.t,
  metadata = extraMeta.t[which(colnames(treat.t) %in% (
    strsplit(colnames(tcell),'_1|_2') %>% unlist %>% unique)),],
  col.name = NULL
)

# feature plot ====
pdf('keyGene_featurePlot.pdf', width = 17,height = 40,onefile = F)
FeaturePlot(epi, keyG, cols = c("#c6ffdd", "#fbd786", "#f7797d"),reduction = 'tsne')
dev.off()

pdf('keyGene_featurePlot_t.pdf', width = 17,height = 40,onefile = F)
FeaturePlot(tcell, keyG, cols = c("#c6ffdd", "#fbd786", "#f7797d"),reduction = 'tsne')
dev.off()

pdf('keyGene_featurePlot_mye.pdf', width = 17,height = 40,onefile = F)
FeaturePlot(myeloid, keyG, cols = c("#c6ffdd", "#fbd786", "#f7797d"),reduction = 'tsne')
dev.off()

# violin plot ====
pdf('keyGene_VlnPlot.pdf', width = 17,height = 35,onefile = F)
VlnPlot(epi,keyG, pt.size = 0,group.by = 'pseudoType')
dev.off()

pdf('keyGene_VlnPlot_t.pdf', width = 17,height = 35,onefile = F)
VlnPlot(tcell,keyG, pt.size = 0,group.by = 'tcellType')
dev.off()

pdf('keyGene_VlnPlot_mye.pdf', width = 17,height = 35,onefile = F)
VlnPlot(myeloid,keyG, pt.size = 0,group.by = 'MyeloidType')
dev.off()

# GSEA ====
if(!dir.exists('C:/Users/10784/Desktop/SCCE/10.2.keyGenes_GSEA')){
  dir.create('C:/Users/10784/Desktop/SCCE/10.2.keyGenes_GSEA',recursive = T)}

setwd('C:/Users/10784/Desktop/SCCE/10.2.keyGenes_GSEA')

# 加载包
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(scater)
library(Seurat)
library(enrichplot)

# 将分组信息添加到Seurat对象中（如果还没有）
Idents(epi) <- epi$pseudoType
Idents(naive.epi) <- naive.epi$pseudoType
Idents(treat.epi) <- treat.epi$pseudoType

# 进行多组比较，找到每个分组相对于其他所有分组的差异基因
degs_all <- FindAllMarkers(epi, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
degs_all.naive <- FindAllMarkers(naive.epi, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
degs_all.treat <- FindAllMarkers(treat.epi, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

degs <- degs_all[keyG,] %>% na.omit()
degs.naive <- degs_all.naive[keyG,]
degs.treat <- degs_all.treat[keyG,]

# 进行GSEA分析
# 提取基因列表并按log2 fold change排序
gene_list <- degs$avg_log2FC
names(gene_list) <- rownames(degs)
gene_list <- sort(gene_list, decreasing = TRUE)

## hallmarks (No Res) ====
# 读取GMT文件
gmt <- read.gmt("h.all.v2023.2.Hs.symbols.gmt")
gmt$term <- gsub('_',' ',gmt$term) %>% gsub('HALLMARK ','',.)

# 定义一个函数，将字符串从全大写转换为首字母大写
capitalize_words <- function(string) {
  # 将字符串转换为小写
  string <- tolower(string)
  # 将每个单词的首字母大写
  string <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", string, perl = TRUE)
  return(string)
}

# 应用函数
gmt$term <- capitalize_words(gmt$term)

# 进行GSEA分析
gsea_results.hm <- GSEA(gene_list, TERM2GENE=gmt, pvalueCutoff=1, verbose=FALSE)

# 保存结果
write.table(as.data.frame(gsea_results.hm), file = "gsea_results_hallmark.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

# 使用dotplot函数可视化
pdf('hallmark_dotplot.pdf', width = 5.5,height = 3.4,onefile = F)
dotplot(gsea_results.hm, showCategory=20)
dev.off()

# 使用gseaplot2函数绘制GSEA曲线图
pdf('hallmark_gseaPlot.pdf', width = 6.5,height = 5.2,onefile = F)
gseaplot2(gsea_results.hm,1:3,base_size=10,ES_geom="line")
dev.off()

# upset plot
pdf('hallmark_upsetPlot.pdf', width = 5.2,height = 3.5,onefile = F)
upsetplot(gsea_results.hm)
dev.off()

# 山脊图
pdf('hallmark_ridgePlot.pdf', width = 6,height = 3.5,onefile = F)
ridgeplot(gsea_results.hm,fill = "pvalue",5)+scale_fill_continuous(type = "viridis")
dev.off()

## Canonical pathways ====
# 读取GMT文件
gmt <- read.gmt("c2.cp.v2023.2.Hs.symbols.gmt")
gmt$term <- gsub('_',' ',gmt$term) %>% gsub('SA |KEGG |PID |REACTOME |WP ','',.)

# 定义一个函数，将字符串从全大写转换为首字母大写
capitalize_words <- function(string) {
  # 将字符串转换为小写
  string <- tolower(string)
  # 将每个单词的首字母大写
  string <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", string, perl = TRUE)
  return(string)
}

# 应用函数
gmt$term <- capitalize_words(gmt$term)

# 进行GSEA分析
gsea_results.cp <- GSEA(gene_list, TERM2GENE=gmt, pvalueCutoff=1, verbose=FALSE)

# 保存结果
write.table(as.data.frame(gsea_results.cp), file = "gsea_results_CanonicalPw.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

# 使用dotplot函数可视化
pdf('CP_dotplot.pdf', width = 6,height = 9,onefile = F)
dotplot(gsea_results.cp, showCategory=20, color = "NES")
dev.off()

# 使用gseaplot2函数绘制GSEA曲线图
pdf('CP_gseaPlot.pdf', width = 7,height = 4.5,onefile = F)
gseaplot2(gsea_results.cp,1:2,pvalue_table = F,base_size=10,ES_geom="line")
dev.off()

# upset plot
pdf('CP_upsetPlot.pdf', width = 5,height = 3,onefile = F)
upsetplot(gsea_results.cp)
dev.off()

# 山脊图
pdf('CP_ridgePlot.pdf', width = 7,height = 3.8,onefile = F)
ridgeplot(gsea_results.cp,fill = "NES",10)+scale_fill_continuous(type = "viridis")
dev.off()

## GO BP ====
# 读取GMT文件
gmt <- read.gmt("c5.go.bp.v2023.2.Hs.symbols.gmt")
gmt$term <- gsub('_',' ',gmt$term) %>% gsub('GOBP ','',.)

# 定义一个函数，将字符串从全大写转换为首字母大写
capitalize_words <- function(string) {
  # 将字符串转换为小写
  string <- tolower(string)
  # 将每个单词的首字母大写
  string <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", string, perl = TRUE)
  return(string)
}

# 应用函数
gmt$term <- capitalize_words(gmt$term)

# 进行GSEA分析
gsea_results.bp <- GSEA(gene_list, TERM2GENE=gmt, pvalueCutoff=1, verbose=T)

# 保存结果
write.table(as.data.frame(gsea_results.bp), file = "gsea_results_GOBP.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

gsea_results.bp <- subset(gsea_results.bp, pvalue<=0.05)
# 使用dotplot函数可视化
pdf('GOBP_dotplot.pdf', width = 6,height = 9,onefile = F)
dotplot(gsea_results.bp,color = "pvalue", showCategory=4)
dev.off()

# 使用gseaplot2函数绘制GSEA曲线图
pdf('GOBP_gseaPlot.pdf', width = 7,height = 4.7,onefile = F)
gseaplot2(gsea_results.bp,1:4,pvalue_table = F,base_size=10,ES_geom="line")
dev.off()

# upset plot
pdf('GOBP_upsetPlot.pdf', width = 5,height = 3,onefile = F)
upsetplot(gsea_results.bp, n=4)
dev.off()

# 山脊图
pdf('GOBP_ridgePlot.pdf', width = 8,height = 4,onefile = F)
ridgeplot(gsea_results.bp,fill = "pvalue",4)+scale_fill_continuous(type = "viridis")
dev.off()