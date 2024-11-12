rm(list = ls());gc()

if(!dir.exists('SCCE/5.pseudo')){
  dir.create('SCCE/5.pseudo',recursive = T)}

setwd('SCCE/5.pseudo')
# library packages ====
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(data.table)
#BiocManager::install('monocle')
library(monocle)
library(data.table)

# load data ====
epi <- readRDS('../3.1.epi_dimRedu/SCCE_Epithelial_reCluster.rds')

# add tumor types ====
TumorType <- paste0('Tumor',c(as.numeric(epi$seurat_clusters))) %>%
  factor(levels = paste0('Tumor',1:13))

epi <- AddMetaData(epi,TumorType,col.name = 'TumorType')
table(epi$TumorType)
table(epi$TumorType,epi$orig.ident)

#(1) count表达矩阵 ====
expr_matrix = GetAssayData(epi, layer = "data")
expr_matrix[1:3,1:3]

#(2) cell meta注释信息 ====
p_data <- epi@meta.data 
head(p_data)
pd <- new('AnnotatedDataFrame', data = p_data) 

#(3) gene meta注释信息 ====
f_data <- data.frame(gene_short_name = row.names(epi),
                     row.names = row.names(epi))
head(f_data)
fd <- new('AnnotatedDataFrame', data = f_data)

#构建cds对象 ====
cds_pre <- newCellDataSet(expr_matrix,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily = negbinomial.size())
##预处理
# Add Size_Factor文库因子 ====
cds_pre <- estimateSizeFactors(cds_pre)
cds_pre$Size_Factor %>% head()
#[1] 1.008993 1.297276 1.169330 1.013851 1.034906 1.023569

# 计算基因表达量的离散度 ====
cds_pre <- estimateDispersions(cds_pre)
head(dispersionTable(cds_pre))

cds_pre

## 策略1：marker gene by Seurat ====
Idents(epi) = "seurat_clusters"
gene_FAM = FindAllMarkers(epi)
gene_sle = gene_FAM %>% 
  dplyr::filter(p_val<0.01) %>% 
  pull(gene) %>% unique()

### 标记所选择的基因 ====
cds <- setOrderingFilter(cds_pre, gene_sle)

#降维(关键步骤)
cds <- reduceDimension(cds, method = 'DDRTree')

#排序,得到轨迹分化相关的若干State
cds <- orderCells(cds)

### plot ====
pdf('monocle_seurat_clusters_clusters.pdf',height = 4,width = 16,onefile = F)
plot_cell_trajectory(cds, color_by = "Pseudotime") +
  facet_wrap("~seurat_clusters", nrow = 1)
dev.off()

dirPath <- 'Seurat'

if(!dir.exists(dirPath)){
  dir.create(dirPath,recursive = T)}

pdf(paste0(dirPath,'/monocle_State_Seurat.pdf'),height = 6,width = 6,onefile = F)
plot_cell_trajectory(cds, color_by = "State")
dev.off()

pdf(paste0(dirPath,'/monocle_epiSubtype_Seurat.pdf'),height = 10,width = 12,onefile = F)
plot_cell_trajectory(cds, color_by = "TumorType")
dev.off()

pdf(paste0(dirPath,'/monocle_epiSubtype_pseudotime_Seurat.pdf'),height = 12,width = 16,onefile = F)
plot_cell_trajectory(cds, color_by = "Pseudotime") +
  facet_wrap("~TumorType", nrow = 3)
dev.off()

pdf(paste0(dirPath,'/monocle_pseudotime_Seurat.pdf'),height = 6,width = 6,onefile = F)
plot_cell_trajectory(cds, color_by = "Pseudotime")
dev.off()

pdf(paste0(dirPath,'/monocle_clusters.pdf'),height = 12,width = 16,
    onefile = F)
plot_cell_trajectory(cds, color_by = "TumorType") +
  facet_wrap("~TumorType", nrow = 3)
dev.off()


## 策略2：high dispersion gene by monocle ====
gene_Disp = dispersionTable(cds_pre)
gene_sle = gene_Disp %>% 
  dplyr::filter(mean_expression >= 0.1,
                dispersion_empirical >= dispersion_fit) %>% 
  pull(gene_id) %>% unique()

### 标记所选择的基因 ====
cds <- setOrderingFilter(cds_pre, gene_sle)

#降维(关键步骤)
cds <- reduceDimension(cds, method = 'DDRTree',auto_param_selection = F)

#排序,得到轨迹分化相关的若干State
cds <- orderCells(cds)

### plot ====
pdf('monocle_seurat_clusters_highDispersion.pdf',height = 4,width = 10,onefile = F)
plot_cell_trajectory(cds, color_by = "Pseudotime") +
  facet_wrap("~seurat_clusters", nrow = 1)
dev.off()

dirPath <- 'highDis'

if(!dir.exists(dirPath)){
  dir.create(dirPath,recursive = T)}

pdf(paste0(dirPath,'/monocle_State_Seurat.pdf'),height = 6,width = 6,onefile = F)
plot_cell_trajectory(cds, color_by = "State")
dev.off()

pdf(paste0(dirPath,'/monocle_epiSubtype_Seurat.pdf'),height = 10,width = 12,onefile = F)
plot_cell_trajectory(cds, color_by = "TumorType")
dev.off()

pdf(paste0(dirPath,'/monocle_epiSubtype_pseudotime_Seurat.pdf'),height = 12,width = 16,onefile = F)
plot_cell_trajectory(cds, color_by = "Pseudotime") +
  facet_wrap("~TumorType", nrow = 3)
dev.off()

pdf(paste0(dirPath,'/monocle_pseudotime_Seurat.pdf'),height = 6,width = 6,onefile = F)
plot_cell_trajectory(cds, color_by = "Pseudotime")
dev.off()

pdf(paste0(dirPath,'/monocle_clusters.pdf'),height = 12,width = 16,
    onefile = F)
plot_cell_trajectory(cds, color_by = "TumorType") +
  facet_wrap("~TumorType", nrow = 3)
dev.off()

## 策略3：variable(high dispersion) gene by Seurat ====
gene_sle <- VariableFeatures(epi)

### 标记所选择的基因 ====
cds <- setOrderingFilter(cds_pre, gene_sle)

#降维(关键步骤)
cds <- reduceDimension(cds, method = 'DDRTree')

#排序,得到轨迹分化相关的若干State
cds <- orderCells(cds)

### plot ====
pdf('monocle_seurat_clusters_variable.pdf',height = 4,width = 16,onefile = F)
plot_cell_trajectory(cds, color_by = "Pseudotime") +
  facet_wrap("~seurat_clusters", nrow = 1)
dev.off()

dirPath <- 'variable'

if(!dir.exists(dirPath)){
  dir.create(dirPath,recursive = T)}

pdf(paste0(dirPath,'/monocle_State_Seurat.pdf'),height = 6,width = 6,onefile = F)
plot_cell_trajectory(cds, color_by = "State")
dev.off()

pdf(paste0(dirPath,'/monocle_epiSubtype_Seurat.pdf'),height = 10,width = 12,onefile = F)
plot_cell_trajectory(cds, color_by = "TumorType")
dev.off()

pdf(paste0(dirPath,'/monocle_epiSubtype_pseudotime_Seurat.pdf'),height = 12,width = 16,onefile = F)
plot_cell_trajectory(cds, color_by = "Pseudotime") +
  facet_wrap("~TumorType", nrow = 3)
dev.off()

pdf(paste0(dirPath,'/monocle_pseudotime_Seurat.pdf'),height = 6,width = 6,onefile = F)
plot_cell_trajectory(cds, color_by = "Pseudotime")
dev.off()

pdf(paste0(dirPath,'/monocle_clusters.pdf'),height = 12,width = 16,
    onefile = F)
plot_cell_trajectory(cds, color_by = "TumorType") +
  facet_wrap("~TumorType", nrow = 3)
dev.off()

setwd('Seurat')
# 鉴定轨迹分化相关基因 =====
diff_pseudo <- differentialGeneTest(cds[gene_sle,], cores = 1, 
                                    fullModelFormulaStr = "~sm.ns(Pseudotime)")
head(diff_pseudo)
table(diff_pseudo$qval<0.05)
#FALSE  TRUE 
#2570  2025

## 选取其中最显著的进行可视化 ====
diff_pseudo_gene <- diff_pseudo %>% 
  dplyr::arrange(qval) %>% 
  rownames()

#p <- plot_pseudotime_heatmap(cds[diff_pseudo_gene,], 
#                             num_clusters = 4,  # default 6
#                             return_heatmap=F)

#pdf('pseudotime_heatmap.pdf',width = 3,height = 4,onefile = F)
#p
#dev.off()

## 获得具体每个cluster的组成基因 ====
#pseudotime_clusters <- cutree(p$tree_row, k = 5) %>% 
#  data.frame(gene = names(.), cluster = . )
#head(pseudotime_clusters)
#table(pseudotime_clusters$cluster)
# 1  2  3  4  5  6  7 
# 5 10  4  5  3  1  1 

#save(list = c('diff_pseudo','diff_pseudo_gene','pseudotime_clusters'),
#     file = 'pseudo_genes.RData')

## 分支点基因变化情况 ====
epi <- AddMetaData(epi,cds$State,col.name = 'pseudotimeState')
Idents(epi) <- epi$pseudotimeState

allMarkers <- FindAllMarkers(object = epi, assay = NULL,
                             only.pos= TRUE,min.pct = .25, logfc.threshold = 0.25) %>% 
  group_by(cluster)

table(allMarkers$cluster)
# 1    2    3    4    5    6    7    8    9   10   11 
#869  409  760  236 1734    0    0 1256  440 1521  493 


allMarkers.sig <- subset(allMarkers, avg_log2FC > 1 & p_val_adj < 0.05)
dim(allMarkers.sig) #2900 7
table(allMarkers.sig$cluster)
#1   2   3   4   5   6   7   8   9  10  11 
#362  34 215   1 959   0   0 476 101 615 137 

saveRDS(allMarkers.sig,'allMarkersSig_pseudotimeState.RDS')
saveRDS(allMarkers,'allMarkers_pseudotimeState.RDS')

# add pseudotime type ====
## pseudoType点基因变化情况 ====
pseudoType <- ifelse(epi$seurat_clusters %in% c(2,9,11),'Initiator',
                     ifelse(epi$seurat_clusters %in% c(3:8,10),'Middle',
                            ifelse(epi$seurat_clusters %in% c(0),'Terminus',
                                   ifelse(epi$seurat_clusters %in% c(1),'IniTerm',
                                          'Universal')))) %>% 
  factor(levels = c('Initiator','Middle','Terminus','IniTerm','Universal'))
epi <- AddMetaData(epi,pseudoType,col.name = 'pseudoType')
table(epi$pseudoType, epi$Sample)

Idents(epi) <- epi$pseudoType

allMarkers <- FindAllMarkers(object = epi, assay = NULL,
                             only.pos= TRUE,min.pct = .25, logfc.threshold = 0.25) %>% 
  group_by(cluster)

table(allMarkers$cluster)
#Initiator    Middle  Terminus   IniTerm Universal 
#887      1457       435       588       312 

allMarkers.sig <- subset(allMarkers, avg_log2FC > 1 & p_val_adj < 0.05)
dim(allMarkers.sig) #1601 7
table(allMarkers.sig$cluster)
#Initiator    Middle  Terminus   IniTerm Universal 
#327       910        58       222        84 

saveRDS(allMarkers.sig,'allMarkersSig_pseudoType.RDS')
saveRDS(allMarkers,'allMarkers_pseudoType.RDS')

#(1) count表达矩阵 ====
expr_matrix = GetAssayData(epi, layer = "data")
expr_matrix[1:3,1:3]

#(2) cell meta注释信息 ====
p_data <- epi@meta.data 
head(p_data)
pd <- new('AnnotatedDataFrame', data = p_data) 

#(3) gene meta注释信息 ====
f_data <- data.frame(gene_short_name = row.names(epi),
                     row.names = row.names(epi))
head(f_data)
fd <- new('AnnotatedDataFrame', data = f_data)

#构建cds对象 ====
cds_pre <- newCellDataSet(expr_matrix,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily = negbinomial.size())
##预处理
# Add Size_Factor文库因子 ====
cds_pre <- estimateSizeFactors(cds_pre)
cds_pre$Size_Factor %>% head()
#[1] 1.008993 1.297276 1.169330 1.013851 1.034906 1.023569

# 计算基因表达量的离散度 ====
cds_pre <- estimateDispersions(cds_pre)
head(dispersionTable(cds_pre))

cds_pre

## 策略1：marker gene by Seurat ====
Idents(epi) = "seurat_clusters"
gene_FAM = FindAllMarkers(epi)
gene_sle = gene_FAM %>% 
  dplyr::filter(p_val<0.01) %>% 
  pull(gene) %>% unique()

### 标记所选择的基因 ====
cds <- setOrderingFilter(cds_pre, gene_sle)

#降维(关键步骤)
cds <- reduceDimension(cds, method = 'DDRTree')

#排序,得到轨迹分化相关的若干State
cds <- orderCells(cds)

## plot ====
pdf('monocle_epiSubtype_pseudoType.pdf',height = 8,width = 7.8,onefile = F)
plot_cell_trajectory(cds, color_by = "pseudoType")
dev.off()

pdf('monocle_epiSubtype_pseudotime_pseudoType.pdf',height = 5.5,width = 7,onefile = F)
plot_cell_trajectory(cds, color_by = "Pseudotime") +
  facet_wrap("~pseudoType", nrow = 2)
dev.off()

pdf('monocle_clusters_pseudoType_pseudoType.pdf',height = 6.2,width = 7,
    onefile = F)
plot_cell_trajectory(cds, color_by = "TumorType") +
  facet_wrap("~pseudoType", nrow = 2)
dev.off()

saveRDS(epi,'epi_withPseudotime.RDS')
