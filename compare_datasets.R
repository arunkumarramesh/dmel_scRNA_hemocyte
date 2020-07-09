library(pheatmap)
library(dplyr)
library(gridExtra)
library(gtools)
library(ggplot2)
library(reshape2)
library(grid)
library(data.table)
library(Matrix)
library(mltools)
library(plyr)
library(Seurat)
library(GO.db)
library(org.Dm.eg.db)
library(cowplot)
library(gridExtra)

#Tattikota et al. (2020) data obtained from GEO, accession GSE146596, converting to Seurat object

expdata_GSE146596 <- read.csv(file = "GSE146596_counts.csv")
x <- data.table(expdata_GSE146596)  # convert x to a data.table
sparseM <- sparsify(x[, !"X"])  # sparsify everything except the name column
rownames(sparseM) <- x$X # set the rownames
expdata_GSE146596 <- sparseM
save(expdata_GSE146596, file ="expdata_GSE146596.Robj")
rownames(expdata_GSE146596) <- mapIds(org.Dm.eg.db, rownames(expdata_GSE146596), column="FLYBASE", keytype="SYMBOL", multiVals="first")
expdata_GSE146596 <- expdata_GSE146596[!is.na(rownames(expdata_GSE146596)),]
GSE146596 <- CreateSeuratObject(counts  = expdata_GSE146596, project = "GSE146596")
GSE146596_meta <- read.csv(file = "GSE146596_metadata_cellclust.csv", header = T, stringsAsFactors = F)
GSE146596 <- subset(GSE146596,cells = GSE146596_meta$cell)
GSE146596$seurat_clusters <- GSE146596_meta$seurat_clusters
Idents(GSE146596) <- GSE146596_meta$seurat_clusters
save(GSE146596, file = "GSE146596.Robj")

#Cattenoz et al. (2020) data obtained from EBI Arrayexpress, accession E‐MTAB‐8698, combined NI and WI and convert to Seurat object

data <- Read10X(data.dir = "/scratch/Arun/Projects/singlecell/filtered_feature_bc_matrix_RRCZ22/filtered_feature_bc_matrix/")
rownames(data) <- mapIds(org.Dm.eg.db, rownames(data), column="FLYBASE", keytype="SYMBOL", multiVals="first")
data <- data[!is.na(rownames(data)),]
RRCZ22 <- CreateSeuratObject(counts  = data, project = "RRCZ22")
RRCZ22_meta <- read.table(file = "NI_cell_cluster_ID.txt", header = T, stringsAsFactors = F)
RRCZ22 <- subset(RRCZ22,cells = RRCZ22_meta$Cell_ID)
RRCZ22$seurat_clusters <- RRCZ22_meta$Cluster_ID
Idents(RRCZ22) <- RRCZ22_meta$Cluster_ID

data <- Read10X(data.dir = "/scratch/Arun/Projects/singlecell/filtered_feature_bc_matrix_RRCZ23/filtered_feature_bc_matrix/")
rownames(data) <- mapIds(org.Dm.eg.db, rownames(data), column="FLYBASE", keytype="SYMBOL", multiVals="first")
data <- data[!is.na(rownames(data)),]
RRCZ23 <- CreateSeuratObject(counts  = data, project = "RRCZ23")
RRCZ23_meta <- read.table(file = "WI_cell_cluster_ID.txt", header = T, stringsAsFactors = F)
RRCZ23 <- subset(RRCZ23,cells = RRCZ23_meta$Cell_ID)
RRCZ23$seurat_clusters <- RRCZ23_meta$Cluster_ID
Idents(RRCZ23) <- RRCZ23_meta$Cluster_ID

expdata_RRCZ2 <- cbind(RRCZ22@assays[["RNA"]]@counts,RRCZ23@assays[["RNA"]]@counts)
RRCZ2 <- CreateSeuratObject(counts  = expdata_RRCZ2, project = "RRCZ2")
RRCZ2$seurat_clusters <- c(RRCZ22$seurat_clusters,RRCZ23$seurat_clusters)
Idents(RRCZ2) <- RRCZ2$seurat_clusters
save(RRCZ2, file = "RRCZ2.Robj")

#integrated_hemocytes as part of singlecell_dmel_hemocytes_jul2020.R

load(file = "GSE146596.Robj")
load("integrated_hemocytes.Robj")
load("RRCZ2.Robj")

RRCZ2 <- NormalizeData(RRCZ2)
RRCZ2 <- FindVariableFeatures(RRCZ2)
RRCZ2 <- ScaleData(RRCZ2)
GSE146596 <- NormalizeData(GSE146596)
GSE146596 <- FindVariableFeatures(GSE146596)
GSE146596 <- ScaleData(GSE146596)

##pairwise CCA alignment for all three datasets
predict.anchors <- FindTransferAnchors(reference = RRCZ2, query = integrated_hemocytes, reduction = 'cca', l2.norm = TRUE)
predictions <- TransferData(anchorset = predict.anchors, refdata = RRCZ2$seurat_clusters,  weight.reduction = "cca", l2.norm = T)
integrated_hemocytes_RRCZ2 <- integrated_hemocytes
integrated_hemocytes_RRCZ2 <- AddMetaData(integrated_hemocytes_RRCZ2, metadata = predictions)
save(integrated_hemocytes_RRCZ2, file = "integrated_hemocytes_RRCZ2.Robj")
rm(integrated_hemocytes_RRCZ2)

predict.anchors <- FindTransferAnchors(reference = RRCZ2, query = GSE146596, reduction = 'cca', l2.norm = TRUE)
predictions <- TransferData(anchorset = predict.anchors, refdata = RRCZ2$seurat_clusters,  weight.reduction = "cca", l2.norm = T)
GSE146596_RRCZ2 <- GSE146596
GSE146596_RRCZ2 <- AddMetaData(GSE146596_RRCZ2, metadata = predictions)
save(GSE146596_RRCZ2, file = "GSE146596_RRCZ2.Robj")
rm(GSE146596_RRCZ2)

predict.anchors <- FindTransferAnchors(reference = GSE146596, query = RRCZ2, reduction = 'cca', l2.norm = TRUE)
predictions <- TransferData(anchorset = predict.anchors, refdata = as.character(GSE146596$seurat_clusters),  weight.reduction = "cca", l2.norm = T)
RRCZ2_GSE146596 <- RRCZ2
RRCZ2_GSE146596 <- AddMetaData(RRCZ2_GSE146596, metadata = predictions)
save(RRCZ2_GSE146596, file = "RRCZ2_GSE146596.Robj")
rm(RRCZ2_GSE146596)

predict.anchors <- FindTransferAnchors(reference = GSE146596, query = integrated_hemocytes, reduction = 'cca', l2.norm = TRUE)
predictions <- TransferData(anchorset = predict.anchors, refdata = as.character(GSE146596$seurat_clusters),  weight.reduction = "cca", l2.norm = T)
integrated_hemocytes_GSE146596 <- integrated_hemocytes
integrated_hemocytes_GSE146596 <- AddMetaData(integrated_hemocytes_GSE146596, metadata = predictions)
save(integrated_hemocytes_GSE146596, file = "integrated_hemocytes_GSE146596.Robj")
rm(integrated_hemocytes_GSE146596)

predict.anchors <- FindTransferAnchors(reference = integrated_hemocytes, query = RRCZ2, reduction = 'cca', l2.norm = TRUE)
predictions <- TransferData(anchorset = predict.anchors, refdata = integrated_hemocytes$seurat_clusters,  weight.reduction = "cca", l2.norm = T)
RRCZ2_integrated_hemocytes <- RRCZ2
RRCZ2_integrated_hemocytes <- AddMetaData(RRCZ2_integrated_hemocytes, metadata = predictions)
save(RRCZ2_integrated_hemocytes, file = "RRCZ2_integrated_hemocytes.Robj")
rm(RRCZ2_integrated_hemocytes)

predict.anchors <- FindTransferAnchors(reference = integrated_hemocytes, query = GSE146596, reduction = 'cca', l2.norm = TRUE)
predictions <- TransferData(anchorset = predict.anchors, refdata = integrated_hemocytes$seurat_clusters,  weight.reduction = "cca", l2.norm = T)
GSE146596_integrated_hemocytes <- GSE146596
GSE146596_integrated_hemocytes <- AddMetaData(GSE146596_integrated_hemocytes, metadata = predictions)
save(GSE146596_integrated_hemocytes, file = "GSE146596_integrated_hemocytes.Robj")
rm(GSE146596_integrated_hemocytes)

predict.anchors <- FindTransferAnchors(reference = integrated_hemocytes, query = integrated_hemocytes, reduction = 'cca', l2.norm = TRUE)
predictions <- TransferData(anchorset = predict.anchors, refdata = integrated_hemocytes$seurat_clusters,  weight.reduction = "cca", l2.norm = T)
integrated_hemocytes_integrated_hemocytes <- integrated_hemocytes
integrated_hemocytes_integrated_hemocytes <- AddMetaData(integrated_hemocytes_integrated_hemocytes, metadata = predictions)
save(integrated_hemocytes_integrated_hemocytes, file = "integrated_hemocytes_integrated_hemocytes.Robj")
rm(integrated_hemocytes_integrated_hemocytes)

predict.anchors <- FindTransferAnchors(reference = GSE146596, query = GSE146596, reduction = 'cca', l2.norm = TRUE)
predictions <- TransferData(anchorset = predict.anchors, refdata = as.character(GSE146596$seurat_clusters),  weight.reduction = "cca", l2.norm = T)
GSE146596_GSE146596 <- GSE146596
GSE146596_GSE146596 <- AddMetaData(GSE146596_GSE146596, metadata = predictions)
save(GSE146596_GSE146596, file = "GSE146596_GSE146596.Robj")
rm(GSE146596_GSE146596)

predict.anchors <- FindTransferAnchors(reference = RRCZ2, query = RRCZ2, reduction = 'cca', l2.norm = TRUE)
predictions <- TransferData(anchorset = predict.anchors, refdata = RRCZ2$seurat_clusters,  weight.reduction = "cca", l2.norm = T)
RRCZ2_RRCZ2 <- RRCZ2
RRCZ2_RRCZ2 <- AddMetaData(RRCZ2_RRCZ2, metadata = predictions)
save(RRCZ2_RRCZ2, file = "RRCZ2_RRCZ2.Robj")
rm(RRCZ2_RRCZ2)

##Get pairwise prediction matrix
load(file = "GSE146596_GSE146596.Robj")
clcmp_GSE146596_GSE146596 <- as.matrix(table(GSE146596_GSE146596$predicted.id,GSE146596_GSE146596$seurat_clusters))
write.csv(clcmp_GSE146596_GSE146596, file = "clcmp_GSE146596_GSE146596.csv")
rm(GSE146596_GSE146596)

load(file = "RRCZ2_RRCZ2.Robj")
clcmp_RRCZ2_RRCZ2 <- as.matrix(table(RRCZ2_RRCZ2$predicted.id,RRCZ2_RRCZ2$seurat_clusters))
write.csv(clcmp_RRCZ2_RRCZ2, file = "clcmp_RRCZ2_RRCZ2.csv")
rm(RRCZ2_RRCZ2)

load(file = "integrated_hemocytes_integrated_hemocytes.Robj")
clcmp_integrated_hemocytes_integrated_hemocytes <- as.matrix(table(integrated_hemocytes_integrated_hemocytes$predicted.id,integrated_hemocytes_integrated_hemocytes$seurat_clusters))
write.csv(clcmp_integrated_hemocytes_integrated_hemocytes, file = "clcmp_integrated_hemocytes_integrated_hemocytes.csv")
rm(integrated_hemocytes_integrated_hemocytes)

load(file = "GSE146596_integrated_hemocytes.Robj")
clcmp_GSE146596_integrated_hemocytes <- as.matrix(table(GSE146596_integrated_hemocytes$predicted.id,GSE146596_integrated_hemocytes$seurat_clusters))
write.csv(clcmp_GSE146596_integrated_hemocytes, file = "clcmp_GSE146596_integrated_hemocytes.csv")
rm(GSE146596_integrated_hemocytes)

load(file = "RRCZ2_integrated_hemocytes.Robj")
clcmp_RRCZ2_integrated_hemocytes <- as.matrix(table(RRCZ2_integrated_hemocytes$predicted.id,RRCZ2_integrated_hemocytes$seurat_clusters))
write.csv(clcmp_RRCZ2_integrated_hemocytes, file = "clcmp_RRCZ2_integrated_hemocytes.csv")
rm(RRCZ2_integrated_hemocytes)

load(file = "integrated_hemocytes_GSE146596.Robj")
clcmp_integrated_hemocytes_GSE146596 <- as.matrix(table(integrated_hemocytes_GSE146596$predicted.id,integrated_hemocytes_GSE146596$seurat_clusters))
write.csv(clcmp_integrated_hemocytes_GSE146596, file = "clcmp_integrated_hemocytes_GSE146596.csv")
rm(integrated_hemocytes_GSE146596)

load(file = "RRCZ2_GSE146596.Robj")
clcmp_RRCZ2_GSE146596 <- as.matrix(table(RRCZ2_GSE146596$predicted.id,RRCZ2_GSE146596$seurat_clusters))
write.csv(clcmp_RRCZ2_GSE146596, file = "clcmp_RRCZ2_GSE146596.csv")
rm(RRCZ2_GSE146596)

load(file = "GSE146596_RRCZ2.Robj")
clcmp_GSE146596_RRCZ2 <- as.matrix(table(GSE146596_RRCZ2$predicted.id,GSE146596_RRCZ2$seurat_clusters))
write.csv(clcmp_GSE146596_RRCZ2, file = "clcmp_GSE146596_RRCZ2.csv")
rm(GSE146596_RRCZ2)

load(file = "integrated_hemocytes_RRCZ2.Robj")
clcmp_integrated_hemocytes_RRCZ2 <- as.matrix(table(integrated_hemocytes_RRCZ2$predicted.id,integrated_hemocytes_RRCZ2$seurat_clusters))
write.csv(clcmp_integrated_hemocytes_RRCZ2, file = "clcmp_integrated_hemocytes_RRCZ2.csv")
rm(integrated_hemocytes_RRCZ2)


##Generate pairwise heatmaps

highmatch <- data.frame()

clcmp <- read.csv(file = "clcmp_RRCZ2_RRCZ2.csv", header = T, row.names = 1)
clcmpprop <- ""
for (i in 1:nrow(clcmp)){
  clcmpproppl <- clcmp[i,]/colSums(clcmp)
  clcmpprop <- rbind(clcmpprop,clcmpproppl)
}
clcmpprop <- clcmpprop[-c(1),]
clcmpprop <- apply(clcmpprop,2,as.numeric)
rownames(clcmpprop) <- rownames(clcmp)
colnames(clcmpprop) <- gsub("PL.","PL-",colnames(clcmpprop))
colnames(clcmpprop) <- gsub("LM.","LM-",colnames(clcmpprop))
clcmpprop <- clcmpprop[,mixedsort(colnames(clcmpprop))]
clcmpprop <- clcmpprop[mixedsort(rownames(clcmpprop)),]
heatmap_RRCZ2_RRCZ2 <- pheatmap(t(clcmpprop), breaks = seq(0,1, by =0.001), color = colorRampPalette(c('#062854','#ffffbf','#d73027'))(1000),cluster_rows = F, cluster_cols = F, main = "R: Cattenoz et al. (2020), Q: Cattenoz et al. (2020)")

clcmp <- read.csv(file = "clcmp_integrated_hemocytes_RRCZ2.csv", header = T, row.names = 1)
clcmpprop <- ""
for (i in 1:nrow(clcmp)){
  clcmpproppl <- clcmp[i,]/colSums(clcmp)
  clcmpprop <- rbind(clcmpprop,clcmpproppl)
}
clcmpprop <- clcmpprop[-c(1),]
clcmpprop <- apply(clcmpprop,2,as.numeric)
rownames(clcmpprop) <- rownames(clcmp)
clcmpprop <- clcmpprop[,mixedsort(colnames(clcmpprop))]
clcmpprop <- clcmpprop[mixedsort(rownames(clcmpprop)),]
highmatch <- rbind(highmatch, melt(clcmpprop)[melt(clcmpprop)$value > 0.5,])
heatmap_integrated_hemocytes_RRCZ2 <- pheatmap(t(clcmpprop), breaks = seq(0,1, by =0.001), color = colorRampPalette(c('#062854','#ffffbf','#d73027'))(1000),cluster_rows = F, cluster_cols = F, main = "R: Cattenoz et al. (2020), Q: This paper")

clcmp <- read.csv(file = "clcmp_GSE146596_RRCZ2.csv", header = T, row.names = 1)
clcmpprop <- ""
for (i in 1:nrow(clcmp)){
  clcmpproppl <- clcmp[i,]/colSums(clcmp)
  clcmpprop <- rbind(clcmpprop,clcmpproppl)
}
clcmpprop <- clcmpprop[-c(1),]
clcmpprop <- apply(clcmpprop,2,as.numeric)
rownames(clcmpprop) <- rownames(clcmp)
colnames(clcmpprop) <- c("PM1","PM2","LM1","PM3","PM4","PM5","CC2","PM8","PM9","PM6","PM7","LM2","PM10","CC1","PM11","PM12","non-hemo")
clcmpprop <- clcmpprop[,mixedsort(colnames(clcmpprop))]
clcmpprop <- clcmpprop[mixedsort(rownames(clcmpprop)),]
highmatch <- rbind(highmatch, melt(clcmpprop)[melt(clcmpprop)$value > 0.5,])
heatmap_GSE146596_RRCZ2 <- pheatmap(t(clcmpprop), breaks = seq(0,1, by =0.001), color = colorRampPalette(c('#062854','#ffffbf','#d73027'))(1000),cluster_rows = F, cluster_cols = F, main = "R: Cattenoz et al. (2020), Q: Tattikota et al. (2020)")


clcmp <- read.csv(file = "clcmp_RRCZ2_integrated_hemocytes.csv", header = T, row.names = 1)
clcmpprop <- ""
for (i in 1:nrow(clcmp)){
  clcmpproppl <- clcmp[i,]/colSums(clcmp)
  clcmpprop <- rbind(clcmpprop,clcmpproppl)
}
clcmpprop <- clcmpprop[-c(1),]
clcmpprop <- apply(clcmpprop,2,as.numeric)
rownames(clcmpprop) <- rownames(clcmp)
colnames(clcmpprop) <- gsub("PL.","PL-",colnames(clcmpprop))
colnames(clcmpprop) <- gsub("LM.","LM-",colnames(clcmpprop))
clcmpprop <- clcmpprop[,mixedsort(colnames(clcmpprop))]
clcmpprop <- clcmpprop[mixedsort(rownames(clcmpprop)),]
highmatch <- rbind(highmatch, melt(clcmpprop)[melt(clcmpprop)$value > 0.5,])
heatmap_RRCZ2_integrated_hemocytes <- pheatmap(t(clcmpprop), breaks = seq(0,1, by =0.001), color = colorRampPalette(c('#062854','#ffffbf','#d73027'))(1000),cluster_rows = F, cluster_cols = F, main = "R: This paper, Q: Cattenoz et al. (2020)")

clcmp <- read.csv(file = "clcmp_integrated_hemocytes_integrated_hemocytes.csv", header = T, row.names = 1)
clcmpprop <- ""
for (i in 1:nrow(clcmp)){
  clcmpproppl <- clcmp[i,]/colSums(clcmp)
  clcmpprop <- rbind(clcmpprop,clcmpproppl)
}
clcmpprop <- clcmpprop[-c(1),]
clcmpprop <- apply(clcmpprop,2,as.numeric)
rownames(clcmpprop) <- rownames(clcmp)
clcmpprop <- clcmpprop[,mixedsort(colnames(clcmpprop))]
clcmpprop <- clcmpprop[mixedsort(rownames(clcmpprop)),]
heatmap_integrated_hemocytes_integrated_hemocytes <- pheatmap(t(clcmpprop), breaks = seq(0,1, by =0.001), color = colorRampPalette(c('#062854','#ffffbf','#d73027'))(1000),cluster_rows = F, cluster_cols = F, main = "R: This paper, Q: This paper")

clcmp <- read.csv(file = "clcmp_GSE146596_integrated_hemocytes.csv", header = T, row.names = 1)
clcmpprop <- ""
for (i in 1:nrow(clcmp)){
  clcmpproppl <- clcmp[i,]/colSums(clcmp)
  clcmpprop <- rbind(clcmpprop,clcmpproppl)
}
clcmpprop <- clcmpprop[-c(1),]
clcmpprop <- apply(clcmpprop,2,as.numeric)
rownames(clcmpprop) <- rownames(clcmp)
colnames(clcmpprop) <- c("PM1","PM2","LM1","PM3","PM4","PM5","CC2","PM8","PM9","PM6","PM7","LM2","PM10","CC1","PM11","PM12","non-hemo")
clcmpprop <- clcmpprop[,mixedsort(colnames(clcmpprop))]
clcmpprop <- clcmpprop[mixedsort(rownames(clcmpprop)),]
highmatch <- rbind(highmatch, melt(clcmpprop)[melt(clcmpprop)$value > 0.5,])
heatmap_GSE146596_integrated_hemocytes <- pheatmap(t(clcmpprop), breaks = seq(0,1, by =0.001), color = colorRampPalette(c('#062854','#ffffbf','#d73027'))(1000),cluster_rows = F, cluster_cols = F, main = "R: This paper, Q: Tattikota et al. (2020)")


clcmp <- read.csv(file = "clcmp_RRCZ2_GSE146596.csv", header = T, row.names = 1)
clcmpprop <- ""
for (i in 1:nrow(clcmp)){
  clcmpproppl <- clcmp[i,]/colSums(clcmp)
  clcmpprop <- rbind(clcmpprop,clcmpproppl)
}
clcmpprop <- clcmpprop[-c(1),]
clcmpprop <- apply(clcmpprop,2,as.numeric)
rownames(clcmpprop) <- rownames(clcmp)
rownames(clcmpprop) <- mapvalues(rownames(clcmpprop),c(0,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18),c("PM1","PM2","LM1","PM3","PM4","PM5","CC2","PM8","PM9","PM6","PM7","LM2","PM10","CC1","PM11","PM12","non-hemo"))
colnames(clcmpprop) <- gsub("PL.","PL-",colnames(clcmpprop))
colnames(clcmpprop) <- gsub("LM.","LM-",colnames(clcmpprop))
clcmpprop <- clcmpprop[,mixedsort(colnames(clcmpprop))]
clcmpprop <- clcmpprop[mixedsort(rownames(clcmpprop)),]
highmatch <- rbind(highmatch, melt(clcmpprop)[melt(clcmpprop)$value > 0.5,])
heatmap_RRCZ2_GSE146596 <- pheatmap(t(clcmpprop), breaks = seq(0,1, by =0.001), color = colorRampPalette(c('#062854','#ffffbf','#d73027'))(1000),cluster_rows = F, cluster_cols = F, main = "R: Tattikota et al. (2020), Q: Cattenoz et al. (2020)")

clcmp <- read.csv(file = "clcmp_integrated_hemocytes_GSE146596.csv", header = T, row.names = 1)
clcmpprop <- ""
for (i in 1:nrow(clcmp)){
  clcmpproppl <- clcmp[i,]/colSums(clcmp)
  clcmpprop <- rbind(clcmpprop,clcmpproppl)
}
clcmpprop <- clcmpprop[-c(1),]
clcmpprop <- apply(clcmpprop,2,as.numeric)
rownames(clcmpprop) <- rownames(clcmp)
rownames(clcmpprop) <- mapvalues(rownames(clcmpprop),c(0,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18),c("PM1","PM2","LM1","PM3","PM4","PM5","CC2","PM8","PM9","PM6","PM7","LM2","PM10","CC1","PM11","PM12","non-hemo"))
clcmpprop <- clcmpprop[,mixedsort(colnames(clcmpprop))]
clcmpprop <- clcmpprop[mixedsort(rownames(clcmpprop)),]
highmatch <- rbind(highmatch, melt(clcmpprop)[melt(clcmpprop)$value > 0.5,])
heatmap_integrated_hemocytes_GSE146596 <- pheatmap(t(clcmpprop), breaks = seq(0,1, by =0.001), color = colorRampPalette(c('#062854','#ffffbf','#d73027'))(1000),cluster_rows = F, cluster_cols = F, main = "R: Tattikota et al. (2020), Q: This paper")

clcmp <- read.csv(file = "clcmp_GSE146596_GSE146596.csv", header = T, row.names = 1)
clcmpprop <- ""
for (i in 1:nrow(clcmp)){
  clcmpproppl <- clcmp[i,]/colSums(clcmp)
  clcmpprop <- rbind(clcmpprop,clcmpproppl)
}
clcmpprop <- clcmpprop[-c(1),]
clcmpprop <- apply(clcmpprop,2,as.numeric)
rownames(clcmpprop) <- rownames(clcmp)
rownames(clcmpprop) <- mapvalues(rownames(clcmpprop),c(0,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18),c("PM1","PM2","LM1","PM3","PM4","PM5","CC2","PM8","PM9","PM6","PM7","LM2","PM10","CC1","PM11","PM12","non-hemo"))
clcmpprop <- clcmpprop[mixedsort(rownames(clcmpprop)),]
colnames(clcmpprop) <- c("PM1","PM2","LM1","PM3","PM4","PM5","CC2","PM8","PM9","PM6","PM7","LM2","PM10","CC1","PM11","PM12","non-hemo")
clcmpprop <- clcmpprop[,mixedsort(colnames(clcmpprop))]
clcmpprop <- clcmpprop[mixedsort(rownames(clcmpprop)),]
heatmap_GSE146596_GSE146596 <- pheatmap(t(clcmpprop), breaks = seq(0,1, by =0.001), color = colorRampPalette(c('#062854','#ffffbf','#d73027'))(1000), cluster_rows = F, cluster_cols = F, main = "R: Tattikota et al. (2020), Q: Tattikota et al. (2020)")


##Combined heatmaps 

heatmap_forlegend <- pheatmap(t(clcmpprop), breaks = seq(0,1, by =0.001), color = colorRampPalette(c('#062854','#ffffbf','#d73027'))(1000), cluster_rows = F, cluster_cols = F, main = "R: Tattikota et al. (2020), Q: Tattikota et al. (2020)", cellheight=30, cellwidth = 30)
y.grob <- textGrob("Query dataset (Q)", 
                   gp=gpar(fontface="bold", col="black", fontsize=18), rot=90)
x.grob <- textGrob("Reference dataset (R)", 
                   gp=gpar(fontface="bold", col="black", fontsize=18))
ccaplot <- grid.arrange(heatmap_GSE146596_GSE146596[[4]][,1:4],heatmap_integrated_hemocytes_GSE146596[[4]][,1:4],heatmap_RRCZ2_GSE146596[[4]][,1:4],heatmap_GSE146596_integrated_hemocytes[[4]][,1:4],heatmap_integrated_hemocytes_integrated_hemocytes[[4]][,1:4],heatmap_RRCZ2_integrated_hemocytes[[4]][,1:4],heatmap_GSE146596_RRCZ2[[4]][,1:4],heatmap_integrated_hemocytes_RRCZ2[[4]][,1:4],heatmap_RRCZ2_RRCZ2[[4]][,1:4], ncol = 3, left = y.grob, bottom = x.grob)
ccalegend <- grid.arrange(heatmap_forlegend[[4]][,5])
lay <- rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,"NA"),
             c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2),
             c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2))
grid.arrange(grobs = list(ccaplot,ccalegend), layout_matrix = lay)
plot_grid(grobs = ccaplot,ccalegend, rel_widths = c(20,1))

pdf(file = "cca_all.pdf", height = 12, width = 15)
grid.arrange(grobs = list(ccaplot,ccalegend), layout_matrix = lay)
dev.off()