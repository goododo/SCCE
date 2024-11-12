#esophageal small-cell carcinoma
rm(list = ls());gc()

setwd('SCCE/0.rawdata')
# library ====
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(data.table)

list.files()
# load SCCE_1N ====
SCCE_1N <- Read10X(data.dir = "SCCE_1N/")
dim(SCCE_1N) #36601  4035

## QC ====
SCCE_1N <- CreateSeuratObject(counts = SCCE_1N, 
                              project = "SCCE_1N",
                              assay = "RNA",
                              names.field = 1,
                              #names.delim = "_",
                              #meta.data = NULL,
                              min.cells = 3, min.features = 200)
dim(SCCE_1N) #22031  3897
head(SCCE_1N,3)

#saveRDS(SCCE_1N,'SCCE_1N_Raw_SeuratObj.rds')

## Mitochondrion
SCCE_1N[["percent.mt"]] <- PercentageFeatureSet(SCCE_1N, pattern = "^MT-")

### plot
pdf('SCCE_1N_QC_plot.pdf',width = 10,height = 5,onefile = F)
VlnPlot(SCCE_1N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        pt.size = 0,ncol = 3)
dev.off()

plot1 <- FeatureScatter(SCCE_1N, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        smooth = F, plot.cor = TRUE, pt.size = 1)
plot2 <- FeatureScatter(SCCE_1N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        smooth = F, plot.cor = TRUE, pt.size = 1)

pdf('SCCE_1N_QC_plot2.pdf',width = 10,height = 5,onefile = F)
plot1 + plot2
dev.off()

## filter
#retain 200< Feature_RNA <10000 cell，filter >20% Mitochondria
SCCE_1N.qc <- subset(SCCE_1N, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 &
                       percent.mt < 20)

## Normalization ====
SCCE_1N.qc <- NormalizeData(object = SCCE_1N.qc, assay = NULL,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000,
                            margin = 1, #If performing CLR normalization, normalize across features (1) or cells (2)
                            verbose = TRUE) # 进行标准化，默认参数
#22031 features across 3703 samples within 1 assay 

## Find features ====
SCCE_1N.qc <- FindVariableFeatures(
  object = SCCE_1N.qc,
  assay = NULL,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = "auto",
  num.bin = 20,
  binning.method = "equal_width",
  nfeatures = 5000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE
) #注意使用的是norm data

saveRDS(SCCE_1N.qc,'SCCE_1N_QC_SeuratObj.rds')

list.files()
# load SCCE_1T ====
SCCE_1T <- Read10X(data.dir = "SCCE_1T/")
dim(SCCE_1T) #36601  7870

## QC ====
SCCE_1T <- CreateSeuratObject(counts = SCCE_1T, 
                              project = "SCCE_1T",
                              assay = "RNA",
                              names.field = 1,
                              #names.delim = "_",
                              #meta.data = NULL,
                              min.cells = 3, min.features = 200)
dim(SCCE_1T) #22686  7863
head(SCCE_1T,3)

#saveRDS(SCCE_1T,'SCCE_1T_Raw_SeuratObj.rds')

## Mitochondrion
SCCE_1T[["percent.mt"]] <- PercentageFeatureSet(SCCE_1T, pattern = "^MT-")

### plot
pdf('SCCE_1T_QC_plot.pdf',width = 10,height = 5,onefile = F)
VlnPlot(SCCE_1T, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        pt.size = 0,ncol = 3)
dev.off()

plot1 <- FeatureScatter(SCCE_1T, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        smooth = F, plot.cor = TRUE, pt.size = 1)
plot2 <- FeatureScatter(SCCE_1T, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        smooth = F, plot.cor = TRUE, pt.size = 1)

pdf('SCCE_1T_QC_plot2.pdf',width = 10,height = 5,onefile = F)
plot1 + plot2
dev.off()

## filter
#retain 200< Feature_RNA <15000 cell，filter >20% Mitochondria
SCCE_1T.qc <- subset(SCCE_1T, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 &
                       percent.mt < 20)

## Normalization ====
SCCE_1T.qc <- NormalizeData(object = SCCE_1T.qc, assay = NULL,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000,
                            margin = 1, #If performing CLR normalization, normalize across features (1) or cells (2)
                            verbose = TRUE) # 进行标准化，默认参数
#22686 features across 7229 samples within 1 assay 

## Find features ====
SCCE_1T.qc <- FindVariableFeatures(
  object = SCCE_1T.qc,
  assay = NULL,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = "auto",
  num.bin = 20,
  binning.method = "equal_width",
  nfeatures = 5000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE
) #注意使用的是norm data

saveRDS(SCCE_1T.qc,'SCCE_1T_QC_SeuratObj.rds')

# find anchors ====
anchors <- FindIntegrationAnchors(object.list = list(
  "SCCE_1N" = SCCE_1N.qc,
  "SCCE_1T" = SCCE_1T.qc
), dims = 1:30)

SCCE <- IntegrateData(anchors, dims = 1:30,
                      normalization.method = "LogNormalize") # "SCT"

##
barcodes <- SCCE@meta.data[["orig.ident"]] %>% strsplit('SCCE_') %>% 
  unlist() %>% 
  matrix(ncol = 2,byrow = T) %>% as.data.frame()

table(barcodes$V2)
#  1N    1T
# 3703  7229

# add barcode annotation
SCCE <- AddMetaData(SCCE,barcodes$V2,col.name = 'Sample')

saveRDS(SCCE,'SCCE_Raw_SeuratObj_seuratIntegration.rds')
