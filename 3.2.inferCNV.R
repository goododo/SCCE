rm(list = ls());gc()

if(!dir.exists('/home/gzy/SCCE/3.2.epi_inferCNV')){
  dir.create('/home/gzy/SCCE/3.2.epi_inferCNV',recursive = T)}

setwd('/home/gzy/SCCE/3.2.epi_inferCNV')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

# load data ====
epi <- readRDS('../3.1.epi_dimRedu/SCCE_Epithelial_reCluster.rds')

# inferCNV ====
#install.packages("rjags")
#BiocManager::install("infercnv")
library(infercnv)
library(rtracklayer)

## gene reference file
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
geneRef <- rtracklayer::import('gencode.v45.annotation.gtf') %>% as.data.frame() %>% 
  subset(type == 'gene' & #seqnames != 'chrX' & 
           seqnames != 'chrY' & seqnames != 'chrM' &
           gene_type %in% c('lncRNA','miRNA','protein_coding'))
table(duplicated(geneRef$gene_name))

geneRef2 <- geneRef[,c('gene_name','seqnames','start','end','width','gene_type')] %>% 
  subset(gene_name %in% rownames(epi))
table(duplicated(geneRef2$gene_name))

### delete duplicated genes
dup_genes <- geneRef2$gene_name[duplicated(geneRef2$gene_name)] %>% unique()

library(data.table)
addgenes <- lapply(dup_genes, function(x){
  dup_rows <- geneRef2[which(geneRef2$gene_name %in% x),] #%>% arrange(desc(width))
  each <- dup_rows[1,]
  return(each)
}) %>% rbindlist()

geneRef3 <- geneRef2[-which(geneRef2$gene_name %in% dup_genes),]
geneRef3 <- rbind(geneRef3, addgenes) %>% .[,1:4] %>% arrange(seqnames,start)
table(duplicated(geneRef3$gene_name))

rownames(geneRef3) <- geneRef3$gene_name

write.table(geneRef3,'geneOrderingFile.txt',quote = F,sep = '\t',row.names = T)

# inferCNV ====
inferCNV_fun <- function(refType){
  # library packages ====
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(data.table)
  
  # inferCNV ====
  #install.packages("rjags")
  #BiocManager::install("infercnv")
  library(infercnv)
  library(rtracklayer)
  
  lapply(unique(epi$Sample), function(patient){
    
    # load data ====
    refData <- readRDS(paste0('../2.1.annoMain/SCCE_',refType,'.rds'))
    refData <- NormalizeData(object = refData, assay = 'RNA',
                             normalization.method = "LogNormalize",
                             scale.factor = 10000,
                             margin = 1, #If performing CLR normalization, normalize across features (1) or cells (2)
                             verbose = TRUE) # 进行标准化，默认参数
    
    ## subset each sample data
    inferData.each <- subset(epi, Sample == patient)
    refData.each <- subset(refData, Sample == patient)
    
    ## sample annotation file
    sampleAnno <- data.frame(row.names = c(colnames(inferData.each),colnames(refData.each)),
                             patient = c(inferData.each$Sample,refData.each$Sample),
                             #group = c(rep('carcinoma',ncol(inferData.each)),
                             #          rep(refType,ncol(refData.each)))
                             group = c(inferData.each$seurat_clusters,
                                       rep(refType, ncol(refData.each)))
    )
    
    sampleAnno <- data.frame(
      row.names = rownames(sampleAnno),
      group = sampleAnno$group)
    
    ## ref data
    geneRef3 <- read.table('geneOrderingFile.txt',sep = '\t',row.names = 1)
    
    ## raw count data
    rawCount1 <- inferData.each@assays$integrated$data[match(geneRef3$gene_name,rownames(inferData.each@assays$integrated)),]
    rawCount1 <- as.matrix(rawCount1)
    rawCount2 <- refData.each@assays$integrated$data[match(geneRef3$gene_name,rownames(refData.each@assays$integrated)),]
    rawCount2 <- as.matrix(rawCount2)
    rawCount <- cbind(rawCount1,rawCount2)
    rm(rawCount1,rawCount2)
    
    ## ref data
    geneRef3 <- read.table('geneOrderingFile.txt',sep = '\t',row.names = 1)
    geneRef3 <- geneRef3[,-1]
    geneRef3$start <- as.numeric(geneRef3$start)
    geneRef3$end <- as.numeric(geneRef3$end)
    
    ## infer CNV
    infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=rawCount, # 可以直接提供矩阵对象
                                         annotations_file=sampleAnno,
                                         #delim="\t",
                                         gene_order_file=geneRef3,
                                         ref_group_names=refType)
    
    ## create output dir
    if(!dir.exists(paste0('./',refType))){
      dir.create(paste0('./',refType),recursive = T)}
    if(!dir.exists(paste0('./',refType))){
      dir.create(paste0('./',refType),recursive = T)}
    
    options("Seurat.object.assay.version" = "v3")
    # perform infercnv operations to reveal cnv signal
    infercnv_obj <- infercnv::run(infercnv_obj,
                                  cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                  out_dir=paste0('./',refType,'/',patient,'_output'),  # 输出文件夹
                                  cluster_by_groups=T,   # 聚类
                                  denoise=T, # 去噪
                                  HMM=T, # 是否基于HMM预测CNV
                                  HMM_type = "i3",
                                  analysis_mode = "samples",#c("subclusters", "samples", "cells")
                                  k_nn = 50,
                                  z_score_filter = 0.8,
                                  ref_subtract_use_mean_bounds = F
                                  )
  })
}

reftype <- gsub('SCCE_|.rds','',grep('.rds',list.files('../2.1.annoMain'),value = T)) %>%
  grep('Epithelial|CellTypes',.,value = T,invert = T)
reftype

load('../inferCNV_fun.RData')
lapply(reftype, inferCNV_fun) # T cell