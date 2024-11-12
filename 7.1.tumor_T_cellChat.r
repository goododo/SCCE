rm(list = ls());gc()

if(!dir.exists('C:/Users/10784/Desktop/SCCE/7.1.tumor_T_cellChat')){
  dir.create('C:/Users/10784/Desktop/SCCE/7.1.tumor_T_cellChat',recursive = T)}

setwd('C:/Users/10784/Desktop/SCCE/7.1.tumor_T_cellChat')
# library packages ====
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(data.table)
library(ggplot2)

# load data ====
naive <- readRDS('../0.rawdata/SCCE_1N_QC_SeuratObj.rds')
treat <- readRDS('../0.rawdata/SCCE_1T_QC_SeuratObj.rds')

epi <- readRDS('../5.1.pseudo/Seurat/epi_withPseudotime.RDS')
extraMeta.epi <- epi@meta.data[,c(5,8,22:24)]

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

# cell chat ====
#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
options(stringsAsFactors = FALSE)

## Import ligand-receptor interaction database ====
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)
# Show the structure of the database
#dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB

cellChat_fun <- function(cell1,cell2,iden1.group,iden2.group,CellChatDB.use,patient){
  # Part I: Data input & processing and initialization of CellChat object ====
  ## prepare inputs ====
  data.input1 <- cell1@assays$RNA@layers$data
  dimnames(data.input1) <- dimnames(cell1)
  data.input2 <- cell2@assays$RNA@layers$data
  dimnames(data.input2) <- dimnames(cell2)
  
  identity1 <- data.frame(group = iden1.group,
                          type = cell1$pseudoType,
                          row.names = names(cell1$TumorType))
  identity2 <- data.frame(group = iden2.group,
                          type = cell2$tcellType,
                          row.names = names(cell2$tcellType))
  
  if (length(which(dimnames(cell2)[[2]]%in% dimnames(cell1)[[2]]))!=0) {
    data.input2 <- data.input2[,-c(which(dimnames(cell2)[[2]]%in% dimnames(cell1)[[2]]))]
    identity2 <- identity2[-c(which(dimnames(cell2)[[2]]%in% dimnames(cell1)[[2]])),]
  }
  
  data.input <- cbind(data.input1,data.input2)
  identity <- rbind(identity1,identity2) # create a dataframe consisting of the cell labels
  table(identity$group,identity$type) # check the cell labels
  
  rm(data.input1,data.input2,identity1,identity2)
  ## create cell chat obj ====
  cellchat <- createCellChat(object = data.input, meta = identity, group.by = "type")
  
  ### Add cell information into *meta* slot of the object  (Optional)
  #cellchat <- addMeta(cellchat, meta = identity, meta.name = 'labels')
  cellchat <- setIdent(cellchat, ident.use = "type") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  
  ## Preprocessing the expression data ====
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  # cellchat <- projectData(cellchat, PPI.human)
  
  # Part II: Inference of cell-cell communication network ====
  ## Compute the communication probability and infer cellular communication network ====
  cellchat <- computeCommunProb(cellchat,type = "truncatedMean",#c("triMean", "truncatedMean", "thresholdedMean", "median")
                                trim = 0.1,
                                LR.use = NULL,
                                raw.use = TRUE,
                                population.size = FALSE,
                                distance.use = TRUE,
                                interaction.length = 200,
                                scale.distance = 0.01,
                                k.min = 10,
                                nboot = 100,
                                seed.use = 1L,
                                Kh = 0.5,
                                n = 1 )
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 3)
  
  ## Infer the cell-cell communication at a signaling pathway level ====
  cellchat <- computeCommunProbPathway(cellchat)
  
  ## Calculate the aggregated cell-cell communication network ====
  cellchat <- aggregateNet(cellchat)
  
  saveRDS(cellchat,paste0(patient,'_cellChat.RDS'))
  
}
cellChat_fun(cell1 = naive.epi,cell2 = naive.t,iden1.group = 'Carcinoma',iden2.group = 'T',
             CellChatDB.use = CellChatDB,patient = 'Naive_pseudo')
cellChat_fun(cell1 = treat.epi,cell2 = treat.t,iden1.group = 'Carcinoma',iden2.group = 'T',
             CellChatDB.use = CellChatDB,patient = 'Treat_pseudo')

# combind all cell chat res ====
grep('_cellChat.RDS',dir(),value = T)
object.list <- list(
  Naive = readRDS(grep('_cellChat.RDS',dir(),value = T)[2]),
  Treatment = readRDS(grep('_cellChat.RDS',dir(),value = T)[4])
)
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)

# Visualization ====
## Compare the total number of interactions and interaction strength ====
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1:2)#,group.levels = c('Naive','Treat')
)
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1:2),#group.levels = c('Naive','Treat'),
                           measure = 'weight'
)
gg1 + gg2

## 看数量和强度变化 ====
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T,
                   label.edge= F, edge.weight.max = weight.max[2],
                   edge.width.max = 12,
                   title.name = paste0("Number of interactions - ",
                                       names(object.list)[i]))
}

## 加权重
#> Merge the following slots: 'data.signaling','net', 'netP','meta','idents','var.features','DB',and 'LR'.
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net"),
                           attribute = c("idents","count"))

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T,
                   label.edge= T, edge.weight.max = weight.max[2],
                   edge.width.max = 12,
                   title.name = paste0("Number of interactions - ",
                                       names(object.list)[i]))
}

## 热图展示
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

## Compare the overall information flow of each signaling pathway ====
gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2),
               color.use = c('Naive'='#dcedc2','Treatment'='#ffaaa6'),
               stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2),
               color.use = c('Naive'='#dcedc2','Treatment'='#ffaaa6'),
               stacked = F, do.stat = TRUE)
pdf('informationFlow.pdf',width = 8.5,height = 9,onefile = F)
gg1 + gg2
dev.off()

# Naive ====
## net up & down ====
## Chord diagram (cell) ====
pos.dataset = 'Naive'
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
                                       pos.dataset = pos.dataset,
                                       features.name = features.name,
                                       only.pos = FALSE,
                                       thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = .05)

#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = pos.dataset,
                              ligand.logFC = 0.2, receptor.logFC = NULL)
par(mfrow = c(1,1), xpd=TRUE)
pdf(paste0("Up-regulated signaling pathways - ",
           names(object.list)[1],'.pdf'),width = 6,height = 6)
netVisual_chord_cell(object.list[[1]],signaling = NULL,
                     net = net.up, 
                     slot.name = 'netP',
                     color.use = c('Initiator'='#ff7e67','Middle'='#fdffab',
                                   'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                   'Universal'='#11999e','Tcyto'='#e0f9b5',
                                   'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                   'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways - ",
                                         names(object.list)[1]), legend.pos.x = 5)
dev.off()

## Middle ====
input <- 'Naive Middle - '
pdf(paste0("Up-regulated signaling pathways from ", input,
           names(object.list)[1],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[1]], slot.name = 'netP',
                     net = net.up,
                     sources.use = 'Middle',
                     #targets.use = 'Carcinoma',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use =  c('Initiator'='#ff7e67','Middle'='#fdffab',
                                    'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                    'Universal'='#11999e','Tcyto'='#e0f9b5',
                                    'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                    'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways from ", input,
                                         names(object.list)[1]), legend.pos.x = 5)
dev.off()

pdf(paste0("Up-regulated signaling pathways target ", input,
           names(object.list)[1],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[1]], slot.name = 'netP',
                     net = net.up,
                     #sources.use = 'Carcinoma',
                     targets.use = 'Middle',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use = c('Initiator'='#ff7e67','Middle'='#fdffab',
                                   'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                   'Universal'='#11999e','Tcyto'='#e0f9b5',
                                   'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                   'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways target ", input,
                                         names(object.list)[1]), legend.pos.x = 5)
dev.off()

pdf('netVisual_bubble_Naive_from Middle.pdf',width = 6,height = 7,onefile = F)
netVisual_bubble(object.list[[1]],
                 sources.use = 'Middle',
                 #targets.use = 'Carcinoma',
                 signaling = c('MIF','MK','MHC-I','APP'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

pdf('netVisual_bubble_Naive_target Middle.pdf',width = 3.6,height = 3.5,onefile = F)
netVisual_bubble(object.list[[1]],
                 #sources.use = 'Carcinoma',
                 targets.use = 'Middle',
                 signaling = c('MK','CDH','CADM','FN1'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

## Universal ====
input <- 'Naive Universal - '
pdf(paste0("Up-regulated signaling pathways from ", input,
           names(object.list)[1],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[1]], slot.name = 'netP',
                     net = net.up,
                     sources.use = 'Universal',
                     #targets.use = 'Carcinoma',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use =  c('Initiator'='#ff7e67','Middle'='#fdffab',
                                    'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                    'Universal'='#11999e','Tcyto'='#e0f9b5',
                                    'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                    'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways from ", input,
                                         names(object.list)[1]), legend.pos.x = 5)
dev.off()

pdf(paste0("Up-regulated signaling pathways target ", input,
           names(object.list)[1],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[1]], slot.name = 'netP',
                     net = net.up,
                     #sources.use = 'Carcinoma',
                     targets.use = 'Universal',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use = c('Initiator'='#ff7e67','Middle'='#fdffab',
                                   'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                   'Universal'='#11999e','Tcyto'='#e0f9b5',
                                   'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                   'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways target ", input,
                                         names(object.list)[1]), legend.pos.x = 5)
dev.off()

pdf('netVisual_bubble_Naive_from Universal.pdf',width = 6,height = 10,onefile = F)
netVisual_bubble(object.list[[1]],
                 sources.use = 'Universal',
                 #targets.use = 'Carcinoma',
                 signaling = c('MHC-I','MK','LAMININ','LCK','ICAM'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

pdf('netVisual_bubble_Naive_target Universal.pdf',width = 4,height = 4.3,onefile = F)
netVisual_bubble(object.list[[1]],
                 #sources.use = 'Carcinoma',
                 targets.use = 'Universal',
                 signaling = c('MK','FN1','TIGIT','CD99'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

## all Naive babble plot ====
pdf('netVisual_bubble_from_Naive.pdf',width = 5.8,height = 5.8,onefile = F)
netVisual_bubble(object.list[[1]],
                 sources.use = c('Middle','Universal'),
                 #targets.use = 'Carcinoma',
                 signaling = c('MIF','MK','MHC-I','APP','MHC-I','MK','LCK','ICAM'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

pdf('netVisual_bubble_target_Naive.pdf',width = 5.8,height = 4.3,onefile = F)
netVisual_bubble(object.list[[1]],
                 #sources.use = c('Middle','Universal'),
                 targets.use = c('Middle','Universal'),
                 signaling = c('MK','CDH','CADM','FN1','MK','FN1','TIGIT','CD99'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

# Treatment ====
## Chord diagram (cell) ====
pos.dataset = 'Treatment'
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
                                       pos.dataset = pos.dataset,
                                       features.name = features.name,
                                       only.pos = FALSE,
                                       thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = .05)

#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = pos.dataset,
                              ligand.logFC = 0.2, receptor.logFC = NULL)
par(mfrow = c(1,1), xpd=TRUE)
pdf(paste0("Up-regulated signaling pathways - ",
           names(object.list)[2],'.pdf'),width = 6,height = 6)
netVisual_chord_cell(object.list[[1]],signaling = NULL,
                     net = net.up, 
                     slot.name = 'netP',
                     color.use = c('Initiator'='#ff7e67','Middle'='#fdffab',
                                   'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                   'Universal'='#11999e','Tcyto'='#e0f9b5',
                                   'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                   'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways - ",
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

## Initiator ====
input <- 'Treatment Initiator - '
pdf(paste0("Up-regulated signaling pathways from ", input,
           names(object.list)[2],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP',
                     net = net.up,
                     sources.use = 'Initiator',
                     #targets.use = 'Carcinoma',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use =  c('Initiator'='#ff7e67','Middle'='#fdffab',
                                    'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                    'Universal'='#11999e','Tcyto'='#e0f9b5',
                                    'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                    'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways from ", input,
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

pdf(paste0("Up-regulated signaling pathways target ", input,
           names(object.list)[2],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP',
                     net = net.up,
                     #sources.use = 'Carcinoma',
                     targets.use = 'Initiator',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use = c('Initiator'='#ff7e67','Middle'='#fdffab',
                                   'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                   'Universal'='#11999e','Tcyto'='#e0f9b5',
                                   'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                   'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways target ", input,
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

pdf('netVisual_bubble_Treatment_from Initiator.pdf',width = 3,height = 2.5,onefile = F)
netVisual_bubble(object.list[[2]],
                 sources.use = 'Initiator',
                 #targets.use = 'Carcinoma',
                 signaling = c('PTN','GRN'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

pdf('netVisual_bubble_Treatment_target Initiator.pdf',width = 6,height = 6,onefile = F)
netVisual_bubble(object.list[[2]],
                 #sources.use = 'Carcinoma',
                 targets.use = 'Initiator',
                 signaling = c('PTN','DESMOSOME','PARs','CD99','CD96'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

## Middle ====
input <- 'Treatment Middle - '
pdf(paste0("Up-regulated signaling pathways from ", input,
           names(object.list)[2],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP',
                     net = net.up,
                     sources.use = 'Middle',
                     #targets.use = 'Carcinoma',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use =  c('Initiator'='#ff7e67','Middle'='#fdffab',
                                    'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                    'Universal'='#11999e','Tcyto'='#e0f9b5',
                                    'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                    'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways from ", input,
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

pdf(paste0("Up-regulated signaling pathways target ", input,
           names(object.list)[2],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP',
                     net = net.up,
                     #sources.use = 'Carcinoma',
                     targets.use = 'Middle',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use = c('Initiator'='#ff7e67','Middle'='#fdffab',
                                   'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                   'Universal'='#11999e','Tcyto'='#e0f9b5',
                                   'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                   'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways target ", input,
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

pdf('netVisual_bubble_Treatment_from Middle.pdf',width = 6,height = 7,onefile = F)
netVisual_bubble(object.list[[2]],
                 sources.use = 'Middle',
                 #targets.use = 'Carcinoma',
                 signaling = c('PTN','DESMOSOME','MHC-I','CD99'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

pdf('netVisual_bubble_Treatment_target Middle.pdf',width = 6,height = 6,onefile = F)
netVisual_bubble(object.list[[2]],
                 #sources.use = 'Carcinoma',
                 targets.use = 'Middle',
                 signaling = c('PTN','DESMOSOME','PARs','CD99','CD96'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

## Terminus ====
input <- 'Treatment Terminus - '
pdf(paste0("Up-regulated signaling pathways from ", input,
           names(object.list)[2],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP',
                     net = net.up,
                     sources.use = 'Terminus',
                     #targets.use = 'Carcinoma',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use =  c('Initiator'='#ff7e67','Middle'='#fdffab',
                                    'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                    'Universal'='#11999e','Tcyto'='#e0f9b5',
                                    'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                    'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways from ", input,
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

pdf(paste0("Up-regulated signaling pathways target ", input,
           names(object.list)[2],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP',
                     net = net.up,
                     #sources.use = 'Carcinoma',
                     targets.use = 'Terminus',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use = c('Initiator'='#ff7e67','Middle'='#fdffab',
                                   'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                   'Universal'='#11999e','Tcyto'='#e0f9b5',
                                   'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                   'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways target ", input,
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

pdf('netVisual_bubble_Treatment_from Terminus.pdf',width = 3.5,height = 2.7,onefile = F)
netVisual_bubble(object.list[[2]],
                 sources.use = 'Terminus',
                 #targets.use = 'Carcinoma',
                 signaling = c('PTN'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

pdf('netVisual_bubble_Treatment_target Terminus.pdf',width = 6,height = 6,onefile = F)
netVisual_bubble(object.list[[2]],
                 #sources.use = 'Carcinoma',
                 targets.use = 'Terminus',
                 signaling = c('PTN','DESMOSOME','PARs','CD99','CD96'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

## IniTerm ====
input <- 'Treatment IniTerm - '
pdf(paste0("Up-regulated signaling pathways from ", input,
           names(object.list)[2],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP',
                     net = net.up,
                     sources.use = 'IniTerm',
                     #targets.use = 'Carcinoma',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use =  c('Initiator'='#ff7e67','Middle'='#fdffab',
                                    'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                    'Universal'='#11999e','Tcyto'='#e0f9b5',
                                    'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                    'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways from ", input,
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

pdf(paste0("Up-regulated signaling pathways target ", input,
           names(object.list)[2],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP',
                     net = net.up,
                     #sources.use = 'Carcinoma',
                     targets.use = 'IniTerm',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use = c('Initiator'='#ff7e67','Middle'='#fdffab',
                                   'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                   'Universal'='#11999e','Tcyto'='#e0f9b5',
                                   'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                   'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways target ", input,
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

pdf('netVisual_bubble_Treatment_from IniTerm.pdf',width = 3.5,height = 2.7,onefile = F)
netVisual_bubble(object.list[[2]],
                 sources.use = 'IniTerm',
                 #targets.use = 'Carcinoma',
                 signaling = c('SLURP'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

pdf('netVisual_bubble_Treatment_target IniTerm.pdf',width = 6,height = 6,onefile = F)
netVisual_bubble(object.list[[2]],
                 #sources.use = 'Carcinoma',
                 targets.use = 'IniTerm',
                 signaling = c('PTN','DESMOSOME','PARs','CD99','CD96'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

## Universal ====
input <- 'Treatment Universal - '
pdf(paste0("Up-regulated signaling pathways from ", input,
           names(object.list)[2],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP',
                     net = net.up,
                     sources.use = 'Universal',
                     #targets.use = 'Carcinoma',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use =  c('Initiator'='#ff7e67','Middle'='#fdffab',
                                    'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                    'Universal'='#11999e','Tcyto'='#e0f9b5',
                                    'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                    'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways from ", input,
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

pdf(paste0("Up-regulated signaling pathways target ", input,
           names(object.list)[2],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP',
                     net = net.up,
                     #sources.use = 'Carcinoma',
                     targets.use = 'Universal',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use = c('Initiator'='#ff7e67','Middle'='#fdffab',
                                   'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                   'Universal'='#11999e','Tcyto'='#e0f9b5',
                                   'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                   'ExhrTH'='#c06c84','Teff'='#9896f1'),
                     title.name = paste0("Up-regulated signaling pathways target ", input,
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

pdf('netVisual_bubble_Treatment_from Universal.pdf',width = 6,height = 7,onefile = F)
netVisual_bubble(object.list[[2]],
                 sources.use = 'Universal',
                 #targets.use = 'Carcinoma',
                 signaling = c('PTN','DESMOSOME','MHC-I','MPZ'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

pdf('netVisual_bubble_Treatment_target Universal.pdf',width = 6,height = 6,onefile = F)
netVisual_bubble(object.list[[2]],
                 #sources.use = 'Carcinoma',
                 targets.use = 'Universal',
                 signaling = c('PTN','DESMOSOME','PARs','CD99','CD96'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

## all Naive babble plot ====
pdf('netVisual_bubble_from_Treat.pdf',width = 10.5,height = 4.4,onefile = F)
netVisual_bubble(object.list[[2]],
                 sources.use = c('Initiator','Middle','Terminus','IniTerm','Universal'),
                 #targets.use = 'Carcinoma',
                 signaling = c('SLURP','PTN','GRN','DESMOSOME','MHC-I','CD99','MPZ'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()

pdf('netVisual_bubble_target_Treat.pdf',width = 12,height = 4,onefile = F)
netVisual_bubble(object.list[[2]],
                 #sources.use = c('Middle','Universal'),
                 targets.use = c('Initiator','Middle','Terminus','IniTerm','Universal'),
                 signaling = c('PTN','DESMOSOME','PARs','CD99','CD96'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('Naive'='#66c6ba','Treatment'='#ff9a00'),
                 grid.on = F,
                 angle.x = 45, remove.isolate = T, return.data = F)
dev.off()
























intersect(object.list[[1]]@netP$pathways,object.list[[2]]@netP$pathways)

union(object.list[[1]]@netP$pathways,object.list[[2]]@netP$pathways)

## specific in BoM (chord)
interestPW <- c('MHC-I','PTN','MK')
lapply(interestPW, function(pw){
  pdf(paste0('from Naive Middle ',pw,'.pdf'),width = 4,height = 4,onefile = F)
  print(
    netVisual_chord_gene(object.list[[1]], slot.name = 'net',
                         net = net.up,
                         signaling = pw,
                         sources.use = 'Middle',
                         #targets.use = 'Carcinoma',
                         #remove.isolate = T,
                         show.legend = F,
                         lab.cex = 0.5,
                         #small.gap=0, big.gap=10,
                         color.use = c('Initiator'='#ff7e67','Middle'='#fdffab',
                                       'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                       'Universal'='#11999e','Tcyto'='#e0f9b5',
                                       'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                       'ExhrTH'='#c06c84','Teff'='#9896f1'),
                         title.name = paste0(names(object.list)[3],' - ',pw),
                         legend.pos.x = 5)
  )
  dev.off()
})

interestPW <- c('SPP1','FN1','GALECTIN','GRN','CEACAM')
lapply(interestPW, function(pw){
  pdf(paste0('targetCarcinoma_',pw,'.pdf'),width = 4,height = 4,onefile = F)
  print(
    netVisual_chord_gene(object.list[[3]], slot.name = 'net',
                         net = net.up,
                         signaling = pw,
                         #sources.use = 'Carcinoma',
                         targets.use = 'Carcinoma',
                         #remove.isolate = T,
                         show.legend = F,
                         lab.cex = 0.5,
                         #small.gap=0, big.gap=10,
                         color.use = c('Initiator'='#ff7e67','Middle'='#fdffab',
                                       'Terminus'='#ffd3b6','IniTerm'='#ffaaa5',
                                       'Universal'='#11999e','Tcyto'='#e0f9b5',
                                       'PCTL'='#ffc7c7','Th'='#a6d0e4',
                                       'ExhrTH'='#c06c84','Teff'='#9896f1'),
                         title.name = paste0(names(object.list)[3],' - ',pw),
                         legend.pos.x = 5)
  )
  dev.off()
})


