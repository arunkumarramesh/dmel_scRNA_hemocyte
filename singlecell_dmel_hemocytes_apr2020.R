library(gridExtra)
library(ggplot2)
library(Seurat)
library(pals)
library(pryr)
library(GO.db)
library(org.Dm.eg.db)
library(reshape)
library(cowplot)
library(scales)
library(slingshot)
library(dplyr)
library(rsample)
library(parsnip)
library(viridis)
library(plyr)
library(rgl)

####SECTION 1 - READING IN DATA FROM CELL RANGER ########

##Sample data
mito.genes <- read.table(file = "mito_genes.dat")
mito.genes <- as.vector(mito.genes[,1])

##NSRef-1_Inf_res
data <- Read10X(data.dir = "NSRef-1_Inf_res/outs/filtered_gene_bc_matrices/D.mel/")
sample <- CreateSeuratObject(counts  = data, min.cells = 3, min.features = 200, project = "NSRef-1_Inf")
mito.genes <- mito.genes[mito.genes %in% rownames(sample)]
sample[["percent.mt"]] <- PercentageFeatureSet(sample, features = mito.genes)
NSRef_1_Inf_obj <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )

##NSRef-2_Inf_res
data <- Read10X(data.dir = "NSRef-2_Inf_res/outs/filtered_gene_bc_matrices/D.mel/")
sample <- CreateSeuratObject(counts  = data, min.cells = 3, min.features = 200, project = "NSRef-2_Inf")
mito.genes <- mito.genes[mito.genes %in% rownames(sample)]
sample[["percent.mt"]] <- PercentageFeatureSet(sample, features = mito.genes)
NSRef_2_Inf_obj <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )

##NSRef-3_Inf_res
data <- Read10X(data.dir = "NSRef-3_Inf_res/outs/filtered_gene_bc_matrices/D.mel/")
sample <- CreateSeuratObject(counts  = data, min.cells = 3, min.features = 200, project = "NSRef-3_Inf")
mito.genes <- mito.genes[mito.genes %in% rownames(sample)]
sample[["percent.mt"]] <- PercentageFeatureSet(sample, features = mito.genes)
NSRef_3_Inf_obj <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )

##NSRef-1_Uninf_res
data <- Read10X(data.dir = "NSRef-1_Uninf_res/outs/filtered_gene_bc_matrices/D.mel/")
sample <- CreateSeuratObject(counts  = data, min.cells = 3, min.features = 200, project = "NSRef-1_Uninf")
mito.genes <- mito.genes[mito.genes %in% rownames(sample)]
sample[["percent.mt"]] <- PercentageFeatureSet(sample, features = mito.genes)
NSRef_1_Uninf_obj <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )

##NSRef-2_Uninf_res
data <- Read10X(data.dir = "NSRef-2_Uninf_res/outs/filtered_gene_bc_matrices/D.mel/")
sample <- CreateSeuratObject(counts  = data, min.cells = 3, min.features = 200, project = "NSRef-2_Uninf")
mito.genes <- mito.genes[mito.genes %in% rownames(sample)]
sample[["percent.mt"]] <- PercentageFeatureSet(sample, features = mito.genes)
NSRef_2_Uninf_obj <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )

##NSRef-3_Uninf_res
data <- Read10X(data.dir = "NSRef-3_Uninf_res/outs/filtered_gene_bc_matrices/D.mel/")
sample <- CreateSeuratObject(counts  = data, min.cells = 3, min.features = 200, project = "NSRef-3_Uninf")
mito.genes <- mito.genes[mito.genes %in% rownames(sample)]
sample[["percent.mt"]] <- PercentageFeatureSet(sample, features = mito.genes)
NSRef_3_Uninf_obj <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )


##C1_Inf_res
data <- Read10X(data.dir = "C1_Inf_res/outs/filtered_gene_bc_matrices/D.mel/")
sample <- CreateSeuratObject(counts  = data, min.cells = 3, min.features = 200, project = "C1_Inf")
mito.genes <- mito.genes[mito.genes %in% rownames(sample)]
sample[["percent.mt"]] <- PercentageFeatureSet(sample, features = mito.genes)
C1_Inf_obj <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )


##C2_Inf_res
data <- Read10X(data.dir = "C2_Inf_res/outs/filtered_gene_bc_matrices/D.mel/")
sample <- CreateSeuratObject(counts  = data, min.cells = 3, min.features = 200, project = "C2_Inf")
mito.genes <- mito.genes[mito.genes %in% rownames(sample)]
sample[["percent.mt"]] <- PercentageFeatureSet(sample, features = mito.genes)
C2_Inf_obj <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )

##C3_Inf_res
data <- Read10X(data.dir = "C3_Inf_res/outs/filtered_gene_bc_matrices/D.mel/")
sample <- CreateSeuratObject(counts  = data, min.cells = 3, min.features = 200, project = "C3_Inf")
mito.genes <- mito.genes[mito.genes %in% rownames(sample)]
sample[["percent.mt"]] <- PercentageFeatureSet(sample, features = mito.genes)
C3_Inf_obj <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )


##C1_Uninf_res
data <- Read10X(data.dir = "C1_Uninf_res/outs/filtered_gene_bc_matrices/D.mel/")
sample <- CreateSeuratObject(counts  = data, min.cells = 3, min.features = 200, project = "C1_Uninf")
mito.genes <- mito.genes[mito.genes %in% rownames(sample)]
sample[["percent.mt"]] <- PercentageFeatureSet(sample, features = mito.genes)
C1_Uninf_obj <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )

##C2_Uninf_res
data <- Read10X(data.dir = "C2_Uninf_res/outs/filtered_gene_bc_matrices/D.mel/")
sample <- CreateSeuratObject(counts  = data, min.cells = 3, min.features = 200, project = "C2_Uninf")
mito.genes <- mito.genes[mito.genes %in% rownames(sample)]
sample[["percent.mt"]] <- PercentageFeatureSet(sample, features = mito.genes)
C2_Uninf_obj <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )

##C3_Uninf_res
data <- Read10X(data.dir = "C3_Uninf_res/outs/filtered_gene_bc_matrices/D.mel/")
sample <- CreateSeuratObject(counts  = data, min.cells = 3, min.features = 200, project = "C3_Uninf")
mito.genes <- mito.genes[mito.genes %in% rownames(sample)]
sample[["percent.mt"]] <- PercentageFeatureSet(sample, features = mito.genes)
C3_Uninf_obj <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )

####SECTION 2 - KEEP only He+ or Srp+  CELLS ########

integrate.list <- list(C1_Inf_obj, C1_Uninf_obj,C3_Inf_obj, C3_Uninf_obj,NSRef_1_Inf_obj, NSRef_1_Uninf_obj,NSRef_3_Inf_obj, NSRef_3_Uninf_obj)
names(integrate.list)=c("C1_Inf_obj", "C1_Uninf_obj","C3_Inf_obj", "C3_Uninf_obj","NSRef_1_Inf_obj", "NSRef_1_Uninf_obj","NSRef_3_Inf_obj", "NSRef_3_Uninf_obj")
for (i in 1:length(x = integrate.list)) {
  integrate.list[[i]] <- NormalizeData(object = integrate.list[[i]], verbose = FALSE)
  integrate.list[[i]] <- FindVariableFeatures(object = integrate.list[[i]],  selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
reference.list <- integrate.list[c("C1_Inf_obj", "C1_Uninf_obj","C3_Inf_obj", "C3_Uninf_obj","NSRef_1_Inf_obj", "NSRef_1_Uninf_obj","NSRef_3_Inf_obj", "NSRef_3_Uninf_obj")]
integrate.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50)
integrate.combined <- IntegrateData(anchorset = integrate.anchors, dims = 1:50)
integrate.combined <- subset(integrate.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
DefaultAssay(integrate.combined) <- "integrated"
new <- GetAssayData(object = integrate.combined, assay = "RNA", slot = "data")[c("FBgn0003507"),]>0
new2 <- GetAssayData(object = integrate.combined, assay = "RNA", slot = "data")[c("FBgn0028430"),]>0
new <- new+new2
new <- new[new >0]
length(new)
cells <- gsub("_.","",names(new))
cells <- cells[!(duplicated(cells) | duplicated(cells, fromLast = TRUE)) ]

####SECTION 3 - ROUND 1 CLUSTERING, REMOVING SPERM CELLS ########

C1_Inf_obj <- subset(C1_Inf_obj,cells = cells)
C1_Uninf_obj <- subset(C1_Uninf_obj,cells = cells)
C3_Inf_obj <- subset(C3_Inf_obj,cells = cells)
C3_Uninf_obj <- subset(C3_Uninf_obj,cells = cells)
NSRef_1_Inf_obj <- subset(NSRef_1_Inf_obj,cells = cells)
NSRef_1_Uninf_obj <- subset(NSRef_1_Uninf_obj,cells = cells)
NSRef_3_Inf_obj <- subset(NSRef_3_Inf_obj,cells = cells)
NSRef_3_Uninf_obj <- subset(NSRef_3_Uninf_obj,cells = cells)
integrate.list <- list(C1_Inf_obj, C1_Uninf_obj,C3_Inf_obj, C3_Uninf_obj,NSRef_1_Inf_obj, NSRef_1_Uninf_obj,NSRef_3_Inf_obj, NSRef_3_Uninf_obj)
names(integrate.list)=c("C1_Inf_obj", "C1_Uninf_obj","C3_Inf_obj", "C3_Uninf_obj","NSRef_1_Inf_obj", "NSRef_1_Uninf_obj","NSRef_3_Inf_obj", "NSRef_3_Uninf_obj")
for (i in 1:length(x = integrate.list)) {
  integrate.list[[i]] <- NormalizeData(object = integrate.list[[i]], verbose = FALSE)
  integrate.list[[i]] <- FindVariableFeatures(object = integrate.list[[i]],  selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
reference.list <- integrate.list[c("C1_Inf_obj", "C1_Uninf_obj","C3_Inf_obj", "C3_Uninf_obj","NSRef_1_Inf_obj", "NSRef_1_Uninf_obj","NSRef_3_Inf_obj", "NSRef_3_Uninf_obj")]
integrate.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50)
integrate.combined <- IntegrateData(anchorset = integrate.anchors, dims = 1:50)
integrate.combined <- subset(integrate.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
DefaultAssay(integrate.combined) <- "integrated"

integrate.combined <- ScaleData(integrate.combined, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"), features = rownames(integrate.combined))
integrate.combined <- RunPCA(integrate.combined, npcs = 50, verbose = FALSE)
integrate.combined <- RunUMAP(integrate.combined, reduction = "pca", dims = 1:50)
integrate.combined <- FindNeighbors(integrate.combined, reduction = "pca", dims = 1:50)
s.genes <- read.csv(file = "s.genes.fly.csv")
g2m.genes <- read.csv(file = "g2m.genes.fly.csv")
s.genes <- as.character(s.genes$FlyBaseID)
g2m.genes <- as.character(g2m.genes$FlyBaseID)
s.genes <- s.genes[s.genes %in% rownames(integrate.combined)]
g2m.genes <- g2m.genes[g2m.genes %in% rownames(integrate.combined)]
integrate.combined <- CellCycleScoring(integrate.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
integrate.combined$CC.Difference <- integrate.combined$S.Score - integrate.combined$G2M.Score
integrate.combined$treatment <- integrate.combined$orig.ident
integrate.combined$treatment <- gsub("NSRef-._","",integrate.combined$treatment)
integrate.combined$treatment <- gsub("C._","",integrate.combined$treatment)
integrate.combined$treatment <- gsub("Inf","Infected",integrate.combined$treatment)
integrate.combined$population <- integrate.combined$orig.ident
integrate.combined$population <- gsub("_Uninf","",integrate.combined$population)
integrate.combined$population <- gsub("_Inf","",integrate.combined$population)
integrate.combined$population <- gsub("-","",integrate.combined$population)
integrate.combined$population <- gsub(".$","",integrate.combined$population)
integrate.combined$population <- gsub("C","control",integrate.combined$population)
integrate.combined$poptreat <- paste(integrate.combined$population,integrate.combined$treatment)

integrate.combined_0.2 <- FindClusters(integrate.combined, resolution = 0.2)
sample.markers <- FindAllMarkers(integrate.combined_0.2,  min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sample.markers, file = "sample.markers.list_round1.csv", quote = F)
write.csv(table(Idents(integrate.combined_0.2)),file = "round1_cellcounts.csv")

jpeg("round1_umap.jpeg", width = 750, height = 350)
DimPlot(integrate.combined_0.2, reduction = "umap", label = TRUE)
dev.off()

integrate.combined_0.2 <- subset(integrate.combined_0.2, subset = seurat_clusters != 9)
cells <- gsub("_.","",colnames(integrate.combined_0.2))
cells <- cells[!(duplicated(cells) | duplicated(cells, fromLast = TRUE)) ]
C1_Inf_obj <- subset(C1_Inf_obj,cells = cells)
C1_Uninf_obj <- subset(C1_Uninf_obj,cells = cells)
C3_Inf_obj <- subset(C3_Inf_obj,cells = cells)
C3_Uninf_obj <- subset(C3_Uninf_obj,cells = cells)
NSRef_1_Inf_obj <- subset(NSRef_1_Inf_obj,cells = cells)
NSRef_1_Uninf_obj <- subset(NSRef_1_Uninf_obj,cells = cells)
NSRef_3_Inf_obj <- subset(NSRef_3_Inf_obj,cells = cells)
NSRef_3_Uninf_obj <- subset(NSRef_3_Uninf_obj,cells = cells)

####SECTION 4 - ROUND 2 CLUSTERING, REMOVE MUSCLE AND FATBODY CELLS ########

integrate.list <- list(C1_Inf_obj, C1_Uninf_obj,C3_Inf_obj, C3_Uninf_obj,NSRef_1_Inf_obj, NSRef_1_Uninf_obj,NSRef_3_Inf_obj, NSRef_3_Uninf_obj)
names(integrate.list)=c("C1_Inf_obj", "C1_Uninf_obj","C3_Inf_obj", "C3_Uninf_obj","NSRef_1_Inf_obj", "NSRef_1_Uninf_obj","NSRef_3_Inf_obj", "NSRef_3_Uninf_obj")
for (i in 1:length(x = integrate.list)) {
  integrate.list[[i]] <- NormalizeData(object = integrate.list[[i]], verbose = FALSE)
  integrate.list[[i]] <- FindVariableFeatures(object = integrate.list[[i]],  selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
reference.list <- integrate.list[c("C1_Inf_obj", "C1_Uninf_obj","C3_Inf_obj", "C3_Uninf_obj","NSRef_1_Inf_obj", "NSRef_1_Uninf_obj","NSRef_3_Inf_obj", "NSRef_3_Uninf_obj")]
integrate.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50)
integrate.combined <- IntegrateData(anchorset = integrate.anchors, dims = 1:50)
integrate.combined <- subset(integrate.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
DefaultAssay(integrate.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
integrate.combined <- ScaleData(integrate.combined, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"), features = rownames(integrate.combined))
integrate.combined <- RunPCA(integrate.combined, npcs = 50, verbose = FALSE)
integrate.combined <- RunUMAP(integrate.combined, reduction = "pca", dims = 1:50)
integrate.combined <- FindNeighbors(integrate.combined, reduction = "pca", dims = 1:50)
s.genes <- read.csv(file = "s.genes.fly.csv")
g2m.genes <- read.csv(file = "g2m.genes.fly.csv")
s.genes <- as.character(s.genes$FlyBaseID)
g2m.genes <- as.character(g2m.genes$FlyBaseID)
s.genes <- s.genes[s.genes %in% rownames(integrate.combined)]
g2m.genes <- g2m.genes[g2m.genes %in% rownames(integrate.combined)]
integrate.combined <- CellCycleScoring(integrate.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
integrate.combined$CC.Difference <- integrate.combined$S.Score - integrate.combined$G2M.Score
integrate.combined$treatment <- integrate.combined$orig.ident
integrate.combined$treatment <- gsub("NSRef-._","",integrate.combined$treatment)
integrate.combined$treatment <- gsub("C._","",integrate.combined$treatment)
integrate.combined$treatment <- gsub("Inf","Infected",integrate.combined$treatment)
integrate.combined$population <- integrate.combined$orig.ident
integrate.combined$population <- gsub("_Uninf","",integrate.combined$population)
integrate.combined$population <- gsub("_Inf","",integrate.combined$population)
integrate.combined$population <- gsub("-","",integrate.combined$population)
integrate.combined$population <- gsub(".$","",integrate.combined$population)
integrate.combined$population <- gsub("C","control",integrate.combined$population)
integrate.combined$poptreat <- paste(integrate.combined$population,integrate.combined$treatment)

integrate.combined_0.3 <- FindClusters(integrate.combined, resolution = 0.3)

sample.markers <- FindAllMarkers(integrate.combined_0.3,  min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sample.markers, file = "sample.markers.list_round2.csv", quote = F)
write.csv(table(Idents(integrate.combined_0.3)),file = "round2_cellcounts.csv")

jpeg("round2_umap.jpeg", width = 750, height = 350)
DimPlot(integrate.combined_0.3, reduction = "umap", label = TRUE)
dev.off()

integrate.combined_0.3 <- subset(integrate.combined_0.3, subset = seurat_clusters != 10)
integrate.combined_0.3 <- subset(integrate.combined_0.3, subset = seurat_clusters != 9)
cells <- gsub("_.","",colnames(integrate.combined_0.3))
cells <- cells[!(duplicated(cells) | duplicated(cells, fromLast = TRUE)) ]
C1_Inf_obj <- subset(C1_Inf_obj,cells = cells)
C1_Uninf_obj <- subset(C1_Uninf_obj,cells = cells)
C3_Inf_obj <- subset(C3_Inf_obj,cells = cells)
C3_Uninf_obj <- subset(C3_Uninf_obj,cells = cells)
NSRef_1_Inf_obj <- subset(NSRef_1_Inf_obj,cells = cells)
NSRef_1_Uninf_obj <- subset(NSRef_1_Uninf_obj,cells = cells)
NSRef_3_Inf_obj <- subset(NSRef_3_Inf_obj,cells = cells)
NSRef_3_Uninf_obj <- subset(NSRef_3_Uninf_obj,cells = cells)


####SECTION 5 - ROUND 3 CLUSTERING ########

integrate.list <- list(C1_Inf_obj, C1_Uninf_obj,C3_Inf_obj, C3_Uninf_obj,NSRef_1_Inf_obj, NSRef_1_Uninf_obj,NSRef_3_Inf_obj, NSRef_3_Uninf_obj)
names(integrate.list)=c("C1_Inf_obj", "C1_Uninf_obj","C3_Inf_obj", "C3_Uninf_obj","NSRef_1_Inf_obj", "NSRef_1_Uninf_obj","NSRef_3_Inf_obj", "NSRef_3_Uninf_obj")
for (i in 1:length(x = integrate.list)) {
  integrate.list[[i]] <- NormalizeData(object = integrate.list[[i]], verbose = FALSE)
  integrate.list[[i]] <- FindVariableFeatures(object = integrate.list[[i]],  selection.method = "vst", nfeatures = 2000, verbose = TRUE)
}
reference.list <- integrate.list[c("C1_Inf_obj", "C1_Uninf_obj","C3_Inf_obj", "C3_Uninf_obj","NSRef_1_Inf_obj", "NSRef_1_Uninf_obj","NSRef_3_Inf_obj", "NSRef_3_Uninf_obj")]
integrate.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50, anchor.features = 2000)
integrate.combined <- IntegrateData(anchorset = integrate.anchors, dims = 1:50)
integrate.combined <- subset(integrate.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
DefaultAssay(integrate.combined) <- "integrated"
s.genes <- read.csv(file = "s.genes.fly.csv")
g2m.genes <- read.csv(file = "g2m.genes.fly.csv")
s.genes <- as.character(s.genes$FlyBaseID)
g2m.genes <- as.character(g2m.genes$FlyBaseID)
s.genes <- s.genes[s.genes %in% rownames(integrate.combined)]
g2m.genes <- g2m.genes[g2m.genes %in% rownames(integrate.combined)]
integrate.combined <- CellCycleScoring(integrate.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
integrate.combined$CC.Difference <- integrate.combined$S.Score - integrate.combined$G2M.Score
integrate.combined$treatment <- integrate.combined$orig.ident
integrate.combined$treatment <- gsub("NSRef-._","",integrate.combined$treatment)
integrate.combined$treatment <- gsub("C._","",integrate.combined$treatment)
integrate.combined$treatment <- gsub("Inf","Infected",integrate.combined$treatment)
integrate.combined$population <- integrate.combined$orig.ident
integrate.combined$population <- gsub("_Uninf","",integrate.combined$population)
integrate.combined$population <- gsub("_Inf","",integrate.combined$population)
integrate.combined$population <- gsub("-","",integrate.combined$population)
integrate.combined$population <- gsub(".$","",integrate.combined$population)
integrate.combined$population <- gsub("C","control",integrate.combined$population)
integrate.combined$poptreat <- paste(integrate.combined$population,integrate.combined$treatment)
integrate.combined_cellcycle <- integrate.combined
integrate.combined <- ScaleData(integrate.combined, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"), features = rownames(integrate.combined))
integrate.combined <- RunPCA(integrate.combined, npcs = 50, verbose = FALSE)
integrate.combined <- RunUMAP(integrate.combined, reduction = "pca", dims = 1:50)
integrate.combined <- FindNeighbors(integrate.combined, reduction = "pca", dims = 1:50)

integrate.combined_0 <- FindClusters(integrate.combined, resolution = 0)
integrate.combined_0.01 <- FindClusters(integrate.combined, resolution = 0.01)
integrate.combined_0.02 <- FindClusters(integrate.combined, resolution = 0.02)
integrate.combined_0.03 <- FindClusters(integrate.combined, resolution = 0.03)
integrate.combined_0.04 <- FindClusters(integrate.combined, resolution = 0.04)
integrate.combined_0.05 <- FindClusters(integrate.combined, resolution = 0.05)
integrate.combined_0.06 <- FindClusters(integrate.combined, resolution = 0.06)
integrate.combined_0.07 <- FindClusters(integrate.combined, resolution = 0.07)
integrate.combined_0.08 <- FindClusters(integrate.combined, resolution = 0.08)
integrate.combined_0.09 <- FindClusters(integrate.combined, resolution = 0.09)
integrate.combined_0.1 <- FindClusters(integrate.combined, resolution = 0.1)
integrate.combined_0.2 <- FindClusters(integrate.combined, resolution = 0.2)
integrate.combined_0.3 <- FindClusters(integrate.combined, resolution = 0.3)
integrate.combined_0.4 <- FindClusters(integrate.combined, resolution = 0.4)
integrate.combined_0.5 <- FindClusters(integrate.combined, resolution = 0.5)
integrate.combined_0.6 <- FindClusters(integrate.combined, resolution = 0.6)
integrate.combined_0.7 <- FindClusters(integrate.combined, resolution = 0.7)
integrate.combined_0.8 <- FindClusters(integrate.combined, resolution = 0.8)
integrate.combined_0.9 <- FindClusters(integrate.combined, resolution = 0.9)
integrate.combined_1 <- FindClusters(integrate.combined, resolution = 1)

integrate.combined_0.02$Res.0 <- integrate.combined_0$seurat_clusters
integrate.combined_0.02$Res.0.01 <- integrate.combined_0.01$seurat_clusters
integrate.combined_0.02$Res.0.02 <- integrate.combined_0.02$seurat_clusters
integrate.combined_0.02$Res.0.03 <- integrate.combined_0.03$seurat_clusters
integrate.combined_0.02$Res.0.04 <- integrate.combined_0.04$seurat_clusters
integrate.combined_0.02$Res.0.05 <- integrate.combined_0.05$seurat_clusters
integrate.combined_0.02$Res.0.06 <- integrate.combined_0.06$seurat_clusters
integrate.combined_0.02$Res.0.07 <- integrate.combined_0.07$seurat_clusters
integrate.combined_0.02$Res.0.08 <- integrate.combined_0.08$seurat_clusters
integrate.combined_0.02$Res.0.09 <- integrate.combined_0.09$seurat_clusters
integrate.combined_0.02$Res.0.1 <- integrate.combined_0.1$seurat_clusters
integrate.combined_0.02$Res.0.2 <- integrate.combined_0.2$seurat_clusters
integrate.combined_0.02$Res.0.3 <- integrate.combined_0.3$seurat_clusters
integrate.combined_0.02$Res.0.4 <- integrate.combined_0.4$seurat_clusters
integrate.combined_0.02$Res.0.5 <- integrate.combined_0.5$seurat_clusters
integrate.combined_0.02$Res.0.6 <- integrate.combined_0.6$seurat_clusters
integrate.combined_0.02$Res.0.7 <- integrate.combined_0.7$seurat_clusters
integrate.combined_0.02$Res.0.8 <- integrate.combined_0.8$seurat_clusters
integrate.combined_0.02$Res.0.9 <- integrate.combined_0.9$seurat_clusters
integrate.combined_0.02$Res.1 <- integrate.combined_1$seurat_clusters

integrate.combined_0.3$poptreat <- gsub("control Infected", "No Selection, Infection", integrate.combined_0.3$poptreat)
integrate.combined_0.3$poptreat <- gsub("control Uninf", "No Selection, No Infection", integrate.combined_0.3$poptreat)
integrate.combined_0.3$poptreat <- gsub("NSRef Infected", "Selection, Infection", integrate.combined_0.3$poptreat)
integrate.combined_0.3$poptreat <- gsub("NSRef Uninf", "Selection, No Infection", integrate.combined_0.3$poptreat)
integrate.combined_0.3$poptreat <- factor(integrate.combined_0.3$poptreat, levels = c("No Selection, No Infection", "No Selection, Infection", "Selection, No Infection", "Selection, Infection"))
integrate.combined_0.3$population <- gsub("control", "No Selection", integrate.combined_0.3$population)
integrate.combined_0.3$population <- gsub("NSRef", "Selection", integrate.combined_0.3$population)
integrate.combined_0.3$population <- factor(integrate.combined_0.3$population, levels = c("No Selection","Selection"))
integrate.combined_0.3$treatment <- gsub("Infected", "Infection", integrate.combined_0.3$treatment)
integrate.combined_0.3$treatment <- gsub("Uninf", "No Infection", integrate.combined_0.3$treatment)
integrate.combined_0.3$treatment <- factor(integrate.combined_0.3$treatment, levels = c("No Infection","Infection"))

##cluster marker dotplot
sample.markers_res0.3 <- FindAllMarkers(integrate.combined_0.3,  min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
sample.markers_res0.3_pos <- sample.markers_res0.3[sample.markers_res0.3$avg_logFC >0,]
sample.markers_list_res0.3_pos_top <- sample.markers_res0.3_pos %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
write.csv(sample.markers_res0.3, file = "sample.markers_res0.3.csv")
write.csv(sample.markers_list_res0.3_pos, file = "sample.markers_list_res0.3_pos.csv")
write.csv(sample.markers_list_res0.3_pos_top, file = "sample.markers_list_res0.3_pos_top.csv")

DotPlot(integrate.combined_0.3, features = sample.markers_list_res0.3_pos_top$gene)

##cell cycle and plasmatocyte markers
DotPlot(integrate.combined_0.3, features = c("FBgn0003124","FBgn0003525","FBgn0261385", "FBgn0259896", "FBgn0029167", "FBgn0243514", "FBgn0011828", "FBgn0000299")) +
  scale_x_discrete(labels=c("Col4a1","Pxn","eater","Hml","NimC1", "scra","stg","polo")) + 
  theme(axis.text.x = element_text(face = "italic"))+
  ylab("Cluster")+
  xlab("")

###basic umap marker plots
jpeg("integrate.combined_0.3_umap.jpeg", width = 750, height = 350)
DimPlot(integrate.combined_0.3, reduction = "umap",split.by = "poptreat", label = TRUE)
dev.off()

##lam markers from evans 2014 suggest cluster 1,2,4,5 have variable lam expression
lam_markers <- c("FBgn0261363","FBgn0032422","FBgn0004657","FBgn0034005")
lam_markers_name <- c("ItgaPS4", "atilla", "mys", "PPO3")
lam_markers_plot_unmod <- list()
for (l in 1:length(lam_markers)){
  lam_markers_plot_unmod[[l]] <- VlnPlot(object = integrate.combined_0.3, features = lam_markers[l], pt.size = 0) +
    xlab("") +
    ggtitle(lam_markers_name[l])+
    NoLegend()+
    theme(plot.title = element_text(size=14, face="bold.italic"))
}
do.call("grid.arrange", c(lam_markers_plot_unmod, ncol=2))


##ROUND 3 with CC correction
integrate.combined_cellcycle <- ScaleData(integrate.combined, vars.to.regress = c("S.Score","G2M.Score", "nCount_RNA","nFeature_RNA","percent.mt"), features = rownames(integrate.combined))
integrate.combined_cellcycle <- RunPCA(integrate.combined_cellcycle, npcs = 50, verbose = FALSE)
integrate.combined_cellcycle <- RunUMAP(integrate.combined_cellcycle, reduction = "pca", dims = 1:50)
integrate.combined_cellcycle <- FindNeighbors(integrate.combined_cellcycle, reduction = "pca", dims = 1:50)
integrate.combined_cellcycle_0.3 <- FindClusters(integrate.combined_cellcycle, resolution = 0.3)
integrate.combined_cellcycle_0.3$poptreat <- gsub("control Infected", "No Selection, Infection", integrate.combined_cellcycle_0.3$poptreat)
integrate.combined_cellcycle_0.3$poptreat <- gsub("control Uninf", "No Selection, No Infection", integrate.combined_cellcycle_0.3$poptreat)
integrate.combined_cellcycle_0.3$poptreat <- gsub("NSRef Infected", "Selection, Infection", integrate.combined_cellcycle_0.3$poptreat)
integrate.combined_cellcycle_0.3$poptreat <- gsub("NSRef Uninf", "Selection, No Infection", integrate.combined_cellcycle_0.3$poptreat)
integrate.combined_cellcycle_0.3$poptreat <- factor(integrate.combined_cellcycle_0.3$poptreat, levels = c("No Selection, No Infection", "No Selection, Infection", "Selection, No Infection", "Selection, Infection"))
integrate.combined_cellcycle_0.3$population <- gsub("control", "No Selection", integrate.combined_cellcycle_0.3$population)
integrate.combined_cellcycle_0.3$population <- gsub("NSRef", "Selection", integrate.combined_cellcycle_0.3$population)
integrate.combined_cellcycle_0.3$population <- factor(integrate.combined_cellcycle_0.3$population, levels = c("No Selection","Selection"))
integrate.combined_cellcycle_0.3$treatment <- gsub("Infected", "Infection", integrate.combined_cellcycle_0.3$treatment)
integrate.combined_cellcycle_0.3$treatment <- gsub("Uninf", "No Infection", integrate.combined_cellcycle_0.3$treatment)
integrate.combined_cellcycle_0.3$treatment <- factor(integrate.combined_cellcycle_0.3$treatment, levels = c("No Infection","Infection"))

###basic umap marker plots
jpeg("integrate.combined_cellcycle_0.3_umap.jpeg", width = 750, height = 350)
DimPlot(integrate.combined_cellcycle_0.3, reduction = "umap",split.by = "poptreat", label = TRUE)
dev.off()

##compare to integrate.combined
cellscompare <- ""
for (i in 0:(length(table(Idents(integrate.combined_cellcycle_0.3)))-1)){
  tempcells <- WhichCells(integrate.combined_cellcycle_0.3, idents = i)
  tempcells_data <- subset(integrate.combined_0.3, cells = tempcells)
  celllist <- Idents(tempcells_data)
  celllist <- factor(celllist, levels = 0:(length(table(Idents(integrate.combined_0.3)))-1))
  cellscompare <- rbind(cellscompare,table(celllist))
  
}
cellscompare <- cellscompare[-1,]
rownames(cellscompare) <- 0:(length(table(Idents(integrate.combined_cellcycle_0.3)))-1)
write.csv(cellscompare, file = "row_cc_0.3_col_all_0.3.csv")


####SECTION 6 - DATA INTEGRATION FOR ESTIMATING MEAN FOLD CHANGE ACROSS DETECTED GENES########

integrat.list_15000f <- list(C1_Inf_obj, C1_Uninf_obj,C3_Inf_obj, C3_Uninf_obj,NSRef_1_Inf_obj, NSRef_1_Uninf_obj,NSRef_3_Inf_obj, NSRef_3_Uninf_obj)
names(integrat.list_15000f)=c("C1_Inf_obj", "C1_Uninf_obj","C3_Inf_obj", "C3_Uninf_obj","NSRef_1_Inf_obj", "NSRef_1_Uninf_obj","NSRef_3_Inf_obj", "NSRef_3_Uninf_obj")
for (i in 1:length(x = integrat.list_15000f)) {
  integrat.list_15000f[[i]] <- NormalizeData(object = integrat.list_15000f[[i]], verbose = FALSE)
  integrat.list_15000f[[i]] <- FindVariableFeatures(object = integrat.list_15000f[[i]],  selection.method = "vst", nfeatures = 15000, verbose = TRUE)
}
reference.list_15000f <- integrat.list_15000f[c("C1_Inf_obj", "C1_Uninf_obj","C3_Inf_obj", "C3_Uninf_obj","NSRef_1_Inf_obj", "NSRef_1_Uninf_obj","NSRef_3_Inf_obj", "NSRef_3_Uninf_obj")]
integrate.anchors_15000f <- FindIntegrationAnchors(object.list = reference.list_15000f, dims = 1:50, anchor.features = 15000)
integrate.combined_15000f <- IntegrateData(anchorset = integrate.anchors_15000f, dims = 1:50)
integrate.combined_15000f <- subset(integrate.combined_15000f, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
DefaultAssay(integrate.combined_15000f) <- "integrated"
s.genes <- read.csv(file = "s.genes.fly.csv")
g2m.genes <- read.csv(file = "g2m.genes.fly.csv")
s.genes <- as.character(s.genes$FlyBaseID)
g2m.genes <- as.character(g2m.genes$FlyBaseID)
s.genes <- s.genes[s.genes %in% rownames(integrate.combined_15000f)]
g2m.genes <- g2m.genes[g2m.genes %in% rownames(integrate.combined_15000f)]
integrate.combined_15000f <- CellCycleScoring(integrate.combined_15000f, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
integrate.combined_15000f$CC.Difference <- integrate.combined_15000f$S.Score - integrate.combined_15000f$G2M.Score
integrate.combined_15000f$treatment <- integrate.combined_15000f$orig.ident
integrate.combined_15000f$treatment <- gsub("NSRef-._","",integrate.combined_15000f$treatment)
integrate.combined_15000f$treatment <- gsub("C._","",integrate.combined_15000f$treatment)
integrate.combined_15000f$treatment <- gsub("Inf","Infected",integrate.combined_15000f$treatment)
integrate.combined_15000f$population <- integrate.combined_15000f$orig.ident
integrate.combined_15000f$population <- gsub("_Uninf","",integrate.combined_15000f$population)
integrate.combined_15000f$population <- gsub("_Inf","",integrate.combined_15000f$population)
integrate.combined_15000f$population <- gsub("-","",integrate.combined_15000f$population)
integrate.combined_15000f$population <- gsub(".$","",integrate.combined_15000f$population)
integrate.combined_15000f$population <- gsub("C","control",integrate.combined_15000f$population)
integrate.combined_15000f$poptreat <- paste(integrate.combined_15000f$population,integrate.combined_15000f$treatment)
integrate.combined_15000f_cellcycle <- integrate.combined_15000f
integrate.combined_15000f <- ScaleData(integrate.combined_15000f, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"), features = rownames(integrate.combined_15000f))
integrate.combined_15000f <- RunPCA(integrate.combined_15000f, npcs = 50, verbose = FALSE)
integrate.combined_15000f <- RunUMAP(integrate.combined_15000f, reduction = "pca", dims = 1:50)
integrate.combined_15000f <- FindNeighbors(integrate.combined_15000f, reduction = "pca", dims = 1:50)

integrate.combined_15000f_identschange <- integrate.combined_15000f
integrate.combined_15000f_identschange <- subset(integrate.combined_15000f_identschange, treatment == "Uninf")
Idents(integrate.combined_15000f_identschange) <- integrate.combined_15000f_identschange$population
Nsref_v_control_15000f_markers <- FindMarkers(integrate.combined_15000f_identschange, ident.1 = "control", ident.2 = "NSRef",min.pct = 0, logfc.threshold = 0, test.use = "MAST")
Nsref_v_control_15000f_markers$gene <- rownames(Nsref_v_control_15000f_markers)

integrate.combined_15000f_identschange2 <- integrate.combined_15000f
integrate.combined_15000f_identschange2 <- subset(integrate.combined_15000f_identschange2, population == "control")
Idents(integrate.combined_15000f_identschange2) <- integrate.combined_15000f_identschange2$treatment
control_Inf_v_uninf_15000f_markers <- FindMarkers(integrate.combined_15000f_identschange2, ident.1 = "Uninf", ident.2 = "Infected", min.pct = 0, logfc.threshold = 0, test.use = "MAST")
control_Inf_v_uninf_15000f_markers$gene <- rownames(control_Inf_v_uninf_15000f_markers)

integrate.combined_15000f_identschange3 <- integrate.combined_15000f
integrate.combined_15000f_identschange3 <- subset(integrate.combined_15000f_identschange3, population == "NSRef")
Idents(integrate.combined_15000f_identschange3) <- integrate.combined_15000f_identschange3$treatment
wasp_Inf_v_uninf_15000f_markers <- FindMarkers(integrate.combined_15000f_identschange3, ident.1 = "Uninf", ident.2 = "Infected", min.pct = 0, logfc.threshold = 0, test.use = "MAST")
wasp_Inf_v_uninf_15000f_markers$gene <- rownames(wasp_Inf_v_uninf_15000f_markers)

selvinfection_15000f <- inner_join(control_Inf_v_uninf_15000f_markers,Nsref_v_control_15000f_markers, by =c("gene"="gene"))
infvuninfec_15000f <- inner_join(control_Inf_v_uninf_15000f_markers,wasp_Inf_v_uninf_15000f_markers, by =c("gene"="gene"))

infvuninfec_15000f$avg_logFC.x <- infvuninfec_15000f$avg_logFC.x*-1.44269407046
infvuninfec_15000f$avg_logFC.y <- infvuninfec_15000f$avg_logFC.y*-1.44269407046

selvinfection_15000f$avg_logFC.x <- selvinfection_15000f$avg_logFC.x*-1.44269407046
selvinfection_15000f$avg_logFC.y <- selvinfection_15000f$avg_logFC.y*-1.44269407046


##sel, inf vs no sel, no inf
integrate.combined_15000f_identschange4 <- integrate.combined_15000f
integrate.combined_15000f_identschange4 <- subset(integrate.combined_15000f_identschange4, subset = poptreat == c("NSRef Infected","control Uninf"))
Idents(integrate.combined_15000f_identschange4) <- integrate.combined_15000f_identschange4$treatment
Sel_Inf_v_NoSel_Uninf_15000f_markers <- FindMarkers(integrate.combined_15000f_identschange4, ident.1 = "Uninf", ident.2 = "Infected", min.pct = 0, logfc.threshold = 0, test.use = "MAST")
Sel_Inf_v_NoSel_Uninf_15000f_markers$gene <- rownames(Sel_Inf_v_NoSel_Uninf_15000f_markers)

selinfvnoseluninfec_15000f <- inner_join(control_Inf_v_uninf_15000f_markers,Sel_Inf_v_NoSel_Uninf_15000f_markers, by =c("gene"="gene"))
selinfvnoseluninfec_15000f$avg_logFC.x <- selinfvnoseluninfec_15000f$avg_logFC.x*-1.44269407046
selinfvnoseluninfec_15000f$avg_logFC.y <- selinfvnoseluninfec_15000f$avg_logFC.y*-1.44269407046

pdf("bulkgeneexp_mod2.pdf",width = 12, height = 3.5)
par(mfrow=c(1,3))
par(mar=c(5,6,4,2)+0.8)
par(mgp=c(4,1,0))
heatscatter(selvinfection_15000f$avg_logFC.x,selvinfection_15000f$avg_logFC.y, ylab =expression('                              No Infection\nFold Change after selection for resistance (log2)'),xlab =expression('Fold change after infection (log2)\n                  No Selection'), main = "", bty="L", cex.lab=1.3, cex.axis=1.3, xlim= c(-0.9,1.5), ylim=c(-1.1,1.3))
abline(coef = c(0,1))
points(x = selvinfection_15000f[selvinfection_15000f$gene == "FBgn0261363",]$avg_logFC.x,
       y = selvinfection_15000f[selvinfection_15000f$gene == "FBgn0261363",]$avg_logFC.y,
       pch = 16, col = "purple")
text("atilla", font = 3, x = -0.05+selvinfection_15000f[selvinfection_15000f$gene == "FBgn0261363",]$avg_logFC.x,
     y = 0.15+selvinfection_15000f[selvinfection_15000f$gene == "FBgn0261363",]$avg_logFC.y, cex = 1.3)

points(x = selvinfection_15000f[selvinfection_15000f$gene == "FBgn0032422",]$avg_logFC.x,
       y = selvinfection_15000f[selvinfection_15000f$gene == "FBgn0032422",]$avg_logFC.y,
       pch = 16, col = "magenta")
text("PPO3", font = 3, x = 0.13+selvinfection_15000f[selvinfection_15000f$gene == "FBgn0032422",]$avg_logFC.x,
     y = -0.18+selvinfection_15000f[selvinfection_15000f$gene == "FBgn0032422",]$avg_logFC.y, cex = 1.3)
title("A", adj = 0)

heatscatter(infvuninfec_15000f$avg_logFC.x,infvuninfec_15000f$avg_logFC.y, ylab =expression('        Selected for Resistance\nFold Change after infection (log2)'),xlab =expression('Fold change after infection (log2)\n                  No Selection'), main = "", bty = "L", cex.lab=1.3, cex.axis=1.3, xlim= c(-0.9,1.5), ylim=c(-1.1,1.3))
abline(coef = c(0,1))
points(x = infvuninfec_15000f[infvuninfec_15000f$gene == "FBgn0261363",]$avg_logFC.x,
       y = infvuninfec_15000f[infvuninfec_15000f$gene == "FBgn0261363",]$avg_logFC.y,
       pch = 16, col = "purple")
text("atilla", font = 3, x = +infvuninfec_15000f[infvuninfec_15000f$gene == "FBgn0261363",]$avg_logFC.x,
     y = -0.16+infvuninfec_15000f[infvuninfec_15000f$gene == "FBgn0261363",]$avg_logFC.y, cex = 1.3)
points(x = infvuninfec_15000f[infvuninfec_15000f$gene == "FBgn0032422",]$avg_logFC.x,
       y = infvuninfec_15000f[infvuninfec_15000f$gene == "FBgn0032422",]$avg_logFC.y,
       pch = 16, col = "magenta")
text("PPO3", font = 3, x = infvuninfec_15000f[infvuninfec_15000f$gene == "FBgn0032422",]$avg_logFC.x,
     y = 0.24+infvuninfec_15000f[infvuninfec_15000f$gene == "FBgn0032422",]$avg_logFC.y, cex = 1.3)
title("B", adj = 0)


heatscatter(selinfvnoseluninfec_15000f$avg_logFC.x,selinfvnoseluninfec_15000f$avg_logFC.y, ylab =expression('                  Fold Change after \nselection for resistance and infection (log2)'),xlab =expression('Fold change after infection (log2)\n                  No Selection'), main = "", bty="L", cex.lab=1.3, cex.axis=1.3, xlim= c(-0.9,1.5), ylim=c(-1.1,1.3))
abline(coef = c(0,1))
points(x = selinfvnoseluninfec_15000f[selinfvnoseluninfec_15000f$gene == "FBgn0261363",]$avg_logFC.x,
       y = selinfvnoseluninfec_15000f[selinfvnoseluninfec_15000f$gene == "FBgn0261363",]$avg_logFC.y,
       pch = 16, col = "purple")
text("atilla", font = 3, x = -0+selinfvnoseluninfec_15000f[selinfvnoseluninfec_15000f$gene == "FBgn0261363",]$avg_logFC.x,
     y = -0.2+selinfvnoseluninfec_15000f[selinfvnoseluninfec_15000f$gene == "FBgn0261363",]$avg_logFC.y, cex = 1.3)
points(x = selinfvnoseluninfec_15000f[selinfvnoseluninfec_15000f$gene == "FBgn0032422",]$avg_logFC.x,
       y = selinfvnoseluninfec_15000f[selinfvnoseluninfec_15000f$gene == "FBgn0032422",]$avg_logFC.y,
       pch = 16, col = "magenta")
text("PPO3", font = 3, x = -0.18+selinfvnoseluninfec_15000f[selinfvnoseluninfec_15000f$gene == "FBgn0032422",]$avg_logFC.x,
     y = 0.18+selinfvnoseluninfec_15000f[selinfvnoseluninfec_15000f$gene == "FBgn0032422",]$avg_logFC.y, cex = 1.3)
title("C", adj = 0)

dev.off()

#write table of logFC
logfc_allgenes <- cbind(infvuninfec_15000f$gene,selvinfection_15000f$avg_logFC.x,selvinfection_15000f$avg_logFC.y,infvuninfec_15000f$avg_logFC.y,selinfvnoseluninfec_15000f$avg_logFC.y)
colnames(logfc_allgenes) <- c("Gene","Infection only FC","Selection only FC","Infection FC in selected populations", "Selection and Infection FC")
write.csv(logfc_allgenes, file = "logfc_allgenes.csv")


####SECTION 7 - ROUND 1 LAMELLOCYTE  SUBCLUSTERING, MOVING MISCLASSIFIED CRYSTAL CELLS ########

mod <- subset(integrate.combined_0.3, idents = c(0,1,2,4,5))
cells_lm <- gsub("_.","",colnames(mod))
cells_lm <- cells_lm[!(duplicated(cells_lm) | duplicated(cells_lm, fromLast = TRUE)) ]
C1_Inf_obj_lm <- subset(C1_Inf_obj,cells = cells_lm)
C1_Uninf_obj_lm <- subset(C1_Uninf_obj,cells = cells_lm)
C3_Inf_obj_lm <- subset(C3_Inf_obj,cells = cells_lm)
C3_Uninf_obj_lm <- subset(C3_Uninf_obj,cells = cells_lm)
NSRef_1_Inf_obj_lm <- subset(NSRef_1_Inf_obj,cells = cells_lm)
NSRef_1_Uninf_obj_lm <- subset(NSRef_1_Uninf_obj,cells = cells_lm)
NSRef_3_Inf_obj_lm <- subset(NSRef_3_Inf_obj,cells = cells_lm)
NSRef_3_Uninf_obj_lm <- subset(NSRef_3_Uninf_obj,cells = cells_lm)
integrate.list_lm <- list(C1_Inf_obj_lm, C1_Uninf_obj_lm, C3_Inf_obj_lm, C3_Uninf_obj_lm, NSRef_1_Inf_obj_lm, NSRef_1_Uninf_obj_lm,NSRef_3_Inf_obj_lm, NSRef_3_Uninf_obj_lm)
names(integrate.list_lm)=c("C1_Inf_obj_lm","C1_Uninf_obj_lm", "C3_Inf_obj_lm","C3_Uninf_obj_lm", "NSRef_1_Inf_obj_lm", "NSRef_1_Uninf_obj_lm","NSRef_3_Inf_obj_lm", "NSRef_3_Uninf_obj_lm")
k.filter <- min(200, min(sapply(integrate.list_lm, ncol)))
for (i in 1:length(x = integrate.list_lm)) {
  integrate.list_lm[[i]] <- NormalizeData(object = integrate.list_lm[[i]], verbose = TRUE)
  integrate.list_lm[[i]] <- FindVariableFeatures(object = integrate.list_lm[[i]],  selection.method = "vst", nfeatures = 2000, verbose = TRUE)
}
reference.list_lm <- integrate.list_lm[c("C1_Inf_obj_lm","C1_Uninf_obj_lm","C3_Inf_obj_lm","C3_Uninf_obj_lm","NSRef_1_Inf_obj_lm", "NSRef_1_Uninf_obj_lm","NSRef_3_Inf_obj_lm", "NSRef_3_Uninf_obj_lm")]
integrate.anchors_lm <- FindIntegrationAnchors(object.list = reference.list_lm, dims = 1:50, anchor.features = 2000, k.filter = k.filter)
integrate.combined_lm <- IntegrateData(anchorset = integrate.anchors_lm, dims = 1:50)
integrate.combined_lm <- subset(integrate.combined_lm, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
DefaultAssay(integrate.combined_lm) <- "integrated"
s.genes <- read.csv(file = "s.genes.fly.csv")
g2m.genes <- read.csv(file = "g2m.genes.fly.csv")
s.genes <- as.character(s.genes$FlyBaseID)
g2m.genes <- as.character(g2m.genes$FlyBaseID)
s.genes <- s.genes[s.genes %in% rownames(integrate.combined_lm)]
g2m.genes <- g2m.genes[g2m.genes %in% rownames(integrate.combined_lm)]
integrate.combined_lm <- CellCycleScoring(integrate.combined_lm, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
integrate.combined_lm$CC.Difference <- integrate.combined_lm$S.Score - integrate.combined_lm$G2M.Score
integrate.combined_lm <- ScaleData(integrate.combined_lm, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"), features = rownames(integrate.combined_lm))
integrate.combined_lm <- RunPCA(integrate.combined_lm, npcs = 50, verbose = FALSE)
integrate.combined_lm <- RunUMAP(integrate.combined_lm, reduction = "pca", dims = 1:50)
integrate.combined_lm <- FindNeighbors(integrate.combined_lm, reduction = "pca", dims = 1:50)
integrate.combined_lm$treatment <- integrate.combined_lm$orig.ident
integrate.combined_lm$treatment <- gsub("NSRef-._","",integrate.combined_lm$treatment)
integrate.combined_lm$treatment <- gsub("C._","",integrate.combined_lm$treatment)
integrate.combined_lm$treatment <- gsub("Inf","Infected",integrate.combined_lm$treatment)
integrate.combined_lm$population <- integrate.combined_lm$orig.ident
integrate.combined_lm$population <- gsub("_Uninf","",integrate.combined_lm$population)
integrate.combined_lm$population <- gsub("_Inf","",integrate.combined_lm$population)
integrate.combined_lm$population <- gsub("-","",integrate.combined_lm$population)
integrate.combined_lm$population <- gsub(".$","",integrate.combined_lm$population)
integrate.combined_lm$population <- gsub("C","control",integrate.combined_lm$population)
integrate.combined_lm$poptreat <- paste(integrate.combined_lm$population,integrate.combined_lm$treatment)

integrate.combined_lm_0.3 <- FindClusters(integrate.combined_lm, resolution = 0.3)
DotPlot(integrate.combined_lm_0.3_round1, features = sample.markers_list_res0.3_pos_top[sample.markers_list_res0.3_pos_top$cluster %in% c(3,6,7,8),]$gene)

DotPlot(integrate.combined_lm_0.3_round1, features = c("FBgn0033367", "FBgn0283437")) +
  scale_x_discrete(labels=rev(c("PPO2", "PPO1")))+
  theme(axis.text.x = element_text(face = "italic"))+
  ylab("Cluster")+
  xlab("")

sample.markers_lm_round1_0.3 <- FindAllMarkers(integrate.combined_lm_0.3_round1,  min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
write.csv(sample.markers_lm_round1_0.3, file = "sample.markers_lm_round1_0.3.csv")


##get crystal cells from round 1
integrate.combined_0.3_cc <- subset(integrate.combined_lm_0.3, idents = 5)
integrate.combined_lm_0.3_round1 <- integrate.combined_lm_0.3
cc_misclass <- factor(rep("CC",length(colnames(integrate.combined_0.3_cc))), levels = "CC")
names(cc_misclass) <- colnames(integrate.combined_0.3_cc)

####SECTION 8 - ROUND 2 LAMELLOCYTE  SUBCLUSTERING ########

mod2 <- subset(integrate.combined_lm_0.3, idents = c(0:4))
cells_lm <- gsub("_.","",colnames(mod2))
cells_lm <- cells_lm[!(duplicated(cells_lm) | duplicated(cells_lm, fromLast = TRUE)) ]
C1_Inf_obj_lm <- subset(C1_Inf_obj,cells = cells_lm)
C1_Uninf_obj_lm <- subset(C1_Uninf_obj,cells = cells_lm)
C3_Inf_obj_lm <- subset(C3_Inf_obj,cells = cells_lm)
C3_Uninf_obj_lm <- subset(C3_Uninf_obj,cells = cells_lm)
NSRef_1_Inf_obj_lm <- subset(NSRef_1_Inf_obj,cells = cells_lm)
NSRef_1_Uninf_obj_lm <- subset(NSRef_1_Uninf_obj,cells = cells_lm)
NSRef_3_Inf_obj_lm <- subset(NSRef_3_Inf_obj,cells = cells_lm)
NSRef_3_Uninf_obj_lm <- subset(NSRef_3_Uninf_obj,cells = cells_lm)
integrate.list_lm <- list(C1_Inf_obj_lm, C1_Uninf_obj_lm, C3_Inf_obj_lm, C3_Uninf_obj_lm, NSRef_1_Inf_obj_lm, NSRef_1_Uninf_obj_lm,NSRef_3_Inf_obj_lm, NSRef_3_Uninf_obj_lm)
names(integrate.list_lm)=c("C1_Inf_obj_lm","C1_Uninf_obj_lm", "C3_Inf_obj_lm","C3_Uninf_obj_lm", "NSRef_1_Inf_obj_lm", "NSRef_1_Uninf_obj_lm","NSRef_3_Inf_obj_lm", "NSRef_3_Uninf_obj_lm")
k.filter <- min(200, min(sapply(integrate.list_lm, ncol)))
for (i in 1:length(x = integrate.list_lm)) {
  integrate.list_lm[[i]] <- NormalizeData(object = integrate.list_lm[[i]], verbose = TRUE)
  integrate.list_lm[[i]] <- FindVariableFeatures(object = integrate.list_lm[[i]],  selection.method = "vst", nfeatures = 2000, verbose = TRUE)
}
reference.list_lm <- integrate.list_lm[c("C1_Inf_obj_lm","C1_Uninf_obj_lm","C3_Inf_obj_lm","C3_Uninf_obj_lm","NSRef_1_Inf_obj_lm", "NSRef_1_Uninf_obj_lm","NSRef_3_Inf_obj_lm", "NSRef_3_Uninf_obj_lm")]
integrate.anchors_lm <- FindIntegrationAnchors(object.list = reference.list_lm, dims = 1:50, anchor.features = 2000, k.filter = k.filter)
integrate.combined_lm <- IntegrateData(anchorset = integrate.anchors_lm, dims = 1:50)
integrate.combined_lm <- subset(integrate.combined_lm, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
DefaultAssay(integrate.combined_lm) <- "integrated"
s.genes <- read.csv(file = "s.genes.fly.csv")
g2m.genes <- read.csv(file = "g2m.genes.fly.csv")
s.genes <- as.character(s.genes$FlyBaseID)
g2m.genes <- as.character(g2m.genes$FlyBaseID)
s.genes <- s.genes[s.genes %in% rownames(integrate.combined_lm)]
g2m.genes <- g2m.genes[g2m.genes %in% rownames(integrate.combined_lm)]
integrate.combined_lm <- CellCycleScoring(integrate.combined_lm, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
integrate.combined_lm$CC.Difference <- integrate.combined_lm$S.Score - integrate.combined_lm$G2M.Score
integrate.combined_lm <- ScaleData(integrate.combined_lm, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"), features = rownames(integrate.combined_lm))
integrate.combined_lm <- RunPCA(integrate.combined_lm, npcs = 50, verbose = FALSE)
integrate.combined_lm <- RunUMAP(integrate.combined_lm, reduction = "pca", dims = 1:50)
integrate.combined_lm <- FindNeighbors(integrate.combined_lm, reduction = "pca", dims = 1:50)
integrate.combined_lm$treatment <- integrate.combined_lm$orig.ident
integrate.combined_lm$treatment <- gsub("NSRef-._","",integrate.combined_lm$treatment)
integrate.combined_lm$treatment <- gsub("C._","",integrate.combined_lm$treatment)
integrate.combined_lm$treatment <- gsub("Inf","Infected",integrate.combined_lm$treatment)
integrate.combined_lm$population <- integrate.combined_lm$orig.ident
integrate.combined_lm$population <- gsub("_Uninf","",integrate.combined_lm$population)
integrate.combined_lm$population <- gsub("_Inf","",integrate.combined_lm$population)
integrate.combined_lm$population <- gsub("-","",integrate.combined_lm$population)
integrate.combined_lm$population <- gsub(".$","",integrate.combined_lm$population)
integrate.combined_lm$population <- gsub("C","control",integrate.combined_lm$population)
integrate.combined_lm$poptreat <- paste(integrate.combined_lm$population,integrate.combined_lm$treatment)

integrate.combined_lm_0.2 <- FindClusters(integrate.combined_lm, resolution = 0.2)
integrate.combined_lm_0.2$poptreat <- gsub("control Infected", "No Selection, Infection", integrate.combined_lm_0.2$poptreat)
integrate.combined_lm_0.2$poptreat <- gsub("control Uninf", "No Selection, No Infection", integrate.combined_lm_0.2$poptreat)
integrate.combined_lm_0.2$poptreat <- gsub("NSRef Infected", "Selection, Infection", integrate.combined_lm_0.2$poptreat)
integrate.combined_lm_0.2$poptreat <- gsub("NSRef Uninf", "Selection, No Infection", integrate.combined_lm_0.2$poptreat)
integrate.combined_lm_0.2$poptreat <- factor(integrate.combined_lm_0.2$poptreat, levels = c("No Selection, No Infection", "No Selection, Infection", "Selection, No Infection", "Selection, Infection"))
integrate.combined_lm_0.2$population <- gsub("control", "No Selection", integrate.combined_lm_0.2$population)
integrate.combined_lm_0.2$population <- gsub("NSRef", "Selection", integrate.combined_lm_0.2$population)
integrate.combined_lm_0.2$population <- factor(integrate.combined_lm_0.2$population, levels = c("No Selection","Selection"))
integrate.combined_lm_0.2$treatment <- gsub("Infected", "Infection", integrate.combined_lm_0.2$treatment)
integrate.combined_lm_0.2$treatment <- gsub("Uninf", "No Infection", integrate.combined_lm_0.2$treatment)
integrate.combined_lm_0.2$treatment <- factor(integrate.combined_lm_0.2$treatment, levels = c("No Infection","Infection"))

lmcellsids <- integrate.combined_lm_0.2$seurat_clusters
lmcellsids <- mapvalues(lmcellsids, from = c("0","1","2","3"), to = c("PLASM1","LAM2","LAM3","LAM1"))
lmcellsids <- factor(lmcellsids, levels = c("PLASM1","LAM1","LAM2","LAM3"))
integrate.combined_lm_0.2$seurat_clusters <- lmcellsids
Idents(integrate.combined_lm_0.2) <- lmcellsids

lm_with_cc <- factor(c(as.character(unlist(lmcellsids)), as.character(unlist(cc_misclass))), levels = c("PLASM1","LAM1","LAM2","LAM3","LAM4","CC"))
names(lm_with_cc) <- c(names(lmcellsids),names(cc_misclass))

##check lm subcluster against top markers from all cells
DotPlot(integrate.combined_lm_0.2, features = sample.markers_list_res0.3_pos_top[sample.markers_list_res0.3_pos_top$cluster %in% c(3,6,7,8),]$gene)
DotPlot(integrate.combined_lm_0.2, features = sample.markers_list_res0.3_pos_top[sample.markers_list_res0.3_pos_top$cluster %in% c(0,1,2,4,5),]$gene)
DotPlot(integrate.combined_lm_0.2, features = sample.markers_list_res0.3_pos_top[sample.markers_list_res0.3_pos_top$cluster %in% c(0,1,2,4,5),]$gene)

###lamellocyte markers dotplot
DotPlot(object = integrate.combined_lm_0.2, features = c("FBgn0034005","FBgn0004657","FBgn0032422","FBgn0261363")) +
  scale_x_discrete(labels=rev(c("ItgaPS4", "atilla", "mys", "PPO3")))

##compare lm to integrate.combined
cellscompare <- ""
clustersnames <- names(table(Idents(integrate.combined_lm_0.2)))
for (i in 1:(length(table(Idents(integrate.combined_lm_0.2))))){
  tempcells <- WhichCells(integrate.combined_lm_0.2, idents = clustersnames[i])
  tempcells_data <- subset(integrate.combined_0.3, cells = tempcells)
  celllist <- Idents(tempcells_data)
  celllist <- factor(celllist, levels = 0:(length(table(Idents(integrate.combined_0.3)))-1))
  cellscompare <- rbind(cellscompare,table(celllist))
  
}
cellscompare <- cellscompare[-1,]
rownames(cellscompare) <- clustersnames
cellscompare
write.csv(cellscompare, file = "row_lm_0.2_col_all_0.3.csv")

##cluster phase
phase_table <- table(integrate.combined_lm_0.2$seurat_clusters,integrate.combined_lm_0.2$Phase)
phase_table <- phase_table/rowSums(phase_table)
phase_table <- melt(phase_table)
colnames(phase_table) <- c("Cluster","Phase","Proportion")
ggplot(data=phase_table, aes(x=factor(Cluster,levels = c( "PLASM1","LAM1","LAM2","LAM3")), y=Proportion, fill = Phase)) +
  geom_bar(stat="identity")+
  xlab("Cluster")+
  ylab("Proportion")

####SECTION 10 - TRAJECTORY INFERENCE FOR LAMELLOCYTE DIFFERENTIATION ########

sds_pca <- slingshot(Embeddings(integrate.combined_lm_0.2, "pca"), clusterLabels = integrate.combined_lm_0.2$seurat_clusters, start.clus = "PLASM1", stretch = 0)
integrate.combined_lm_0.2 <- FindVariableFeatures(object = integrate.combined_lm_0.2)
top_hvg <- HVFInfo(integrate.combined_lm_0.2) %>% 
  mutate(., bc = rownames(.)) %>% 
  arrange(desc(variance.standardized)) %>% 
  top_n(2000, variance.standardized) %>% 
  pull(bc)
# Prepare data for random forest
dat_use <- t(GetAssayData(integrate.combined_lm_0.2, slot = "scale.data")[top_hvg,])
dat_use_df <- cbind(slingPseudotime(sds_pca)[,1], dat_use)
colnames(dat_use_df)[1] <- "pseudotime"
dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])
dat_split <- initial_split(dat_use_df)
dat_train <- training(dat_split)
dat_val <- testing(dat_split)
model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 3) %>%
  fit(pseudotime ~ ., data = dat_train)
var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
var_imp_df <- as.data.frame(var_imp)
var_imp_df$symbol <- mapIds(org.Dm.eg.db, rownames(var_imp_df), column="SYMBOL", keytype="FLYBASE", multiVals="first")
write.csv(var_imp_df, file = "var_imp_df_lm_lineage_2000hvg.csv")


##dimplot lm subcluster
DimPlot(object = integrate.combined_lm_0.2, reduction = "umap",split.by = "poptreat", label = TRUE)

lm_dimplot <- DimPlot(object = integrate.combined_lm_0.2, reduction = "umap",label = F, order = rev(c( "PLASM1","LAM1","LAM2","LAM3")), cols = c('#9970ab',"#F28848","#FCBC2A", "#F0F921"))
lm_dimplot_data <- lm_dimplot$data
lm_dimplot_data_mean <- aggregate(lm_dimplot_data[, 1:2], list(lm_dimplot_data$ident), mean)

lm_dimplot_lineage <- lm_dimplot +
  geom_line(data = lm_dimplot_data_mean, aes(UMAP_1,UMAP_2), size = 1)+
  geom_point(data = lm_dimplot_data_mean, aes(UMAP_1,UMAP_2), size = 2)+ 
  theme(legend.position="right")


#top lineage marker plot
featureplotgenelist <- list()
featuregene <- c("FBgn0013733",	"FBgn0003888")
featuregenename <- c("shot",	"betaTub60D")
for (p in 1:length(featuregene)){
  featureplotgenelist[[p]] <- FeaturePlot(object = integrate.combined_lm_0.2, features = featuregene[p], cols = c("purple4","yellow")) +
    ggtitle(featuregenename[p]) +
    theme(plot.title = element_text(size=14, face="bold.italic"))
}
do.call("grid.arrange", c(featureplotgenelist, ncol=2))


##ppo3 and atilla plot, subclustering
integrate.combined_lm_0.2_forplot = SetAssayData(object = integrate.combined_lm_0.2, slot = "data", assay = "integrated", new.data = log2(exp(as.matrix(GetAssayData(object = integrate.combined_lm_0.2, slot = "data", assay = "integrated")))))

atilla_vlnplot <- VlnPlot(object = integrate.combined_lm_0.2_forplot, cols = c('#9970ab',"#F28848","#FCBC2A", "#F0F921"), features = "FBgn0032422", slot = "data",assay = "integrated", pt.size = 0) +NoLegend() +
  ggtitle("atilla") +
  theme(plot.title = element_text(size=14, face="bold.italic")) +
  ylab("Log (2) expression")+
  xlab("")

ppo3_vlnplot <- VlnPlot(object = integrate.combined_lm_0.2_forplot, cols = c('#9970ab',"#F28848","#FCBC2A", "#F0F921"), features = "FBgn0261363", slot = "data",assay = "integrated", pt.size = 0) +NoLegend() +
  ggtitle("PPO3") +
  theme(plot.title = element_text(size=14, face="bold.italic"))+
  ylab("Log (2) expression")+
  xlab("")

lineage_plot <- plot_grid(lm_dimplot_lineage,atilla_vlnplot,ppo3_vlnplot, labels = c('B','C','D'), label_size = 12, ncol = 3, rel_widths = c(1.8,1,1))

##lineage pca
lm_dimplot_pca_legend <- DimPlot(object = integrate.combined_lm_0.2, reduction = "pca",label = T, order = rev(c( "PLASM1","LAM1","LAM2","LAM3")), cols = c('#9970ab',"#F28848","#FCBC2A", "#F0F921")) + NoLegend()
lm_dimplot_pca <- DimPlot(object = integrate.combined_lm_0.2, reduction = "pca",label = F, order = rev(c( "PLASM1","LAM1","LAM2","LAM3")), cols = c('#9970ab',"#F28848","#FCBC2A", "#F0F921")) + NoLegend()

lm_dimplot_pca_data <- lm_dimplot_pca$data
lm_dimplot_pca_data_mean <- aggregate(lm_dimplot_pca_data[, 1:2], list(lm_dimplot_pca_data$ident), mean)

lm_dimplot_pca_lineage <- lm_dimplot_pca +
  geom_line(data = lm_dimplot_pca_data_mean, aes(PC_1,PC_2), size = 1)+
  geom_point(data = lm_dimplot_pca_data_mean, aes(PC_1,PC_2), size = 2)

plot_grid(lm_dimplot_pca_legend,lm_dimplot_pca_lineage,labels = "AUTO", label_size = 12, ncol = 2)

##pathway analysis first with last LAM
sample_markers_lm_0.2_LAM1v2 <- FindMarkers(integrate.combined_lm_0.2, ident.1 = "LAM1", ident.2 = "LAM3",min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
sample_markers_lm_0.2_LAM1v2$gene <- rownames(sample_markers_lm_0.2_LAM1v2)
sample_markers_lm_0.2_LAM1v2$symbol <- mapIds(org.Dm.eg.db, as.character(sample_markers_lm_0.2_LAM1v2$gene), column="SYMBOL", keytype="FLYBASE", multiVals="first")
write.csv(sample_markers_lm_0.2_LAM1v2, file = "markers_LAM1v2.csv")


####SECTION 11 - INTEGRATING ROUND 2 SUBCLUSTERING WITH ROUND 3 CLUSTERING AND GENERATING MAIN TEXT PLOTS ########

identities <- Idents(integrate.combined_0.3)
identities_mod <- ""
identities_mod <-  factor(identities_mod, levels = c(levels(lm_with_cc),3,6:8))
for (n in 1:length(identities)){
  if(names(identities[n]) %in% names(lm_with_cc)){
    identity_mod <- lm_with_cc[names(identities[n])]
    identities_mod[n] <- identity_mod
  } else {
    identities_mod[n] <- identities[n]
  }
}
names(identities_mod) <- names(identities)
identities_mod <- mapvalues(identities_mod, from = c("3","6","7","8"), to = c("PLASM2","CC","MET","AMP"))
identities_mod <- factor(identities_mod, levels = c("PLASM1","PLASM2","CC","MET","AMP","LAM1","LAM2","LAM3"))
integrate.combined_0.3_mod <- integrate.combined_0.3
Idents(integrate.combined_0.3_mod) <- identities_mod
integrate.combined_0.3_mod$seurat_clusters <- identities_mod

##basic umap and marker plots
dimplot_data <- DimPlot(integrate.combined_0.3_mod, reduction = "umap",split.by = "poptreat", label = F, order = rev(c( "PLASM1","PLASM2","MET","AMP","CC","LAM1","LAM2","LAM3")), cols = c('#9970ab','#8c96c6','#99d8c9','#41ae76',"#DC3220","#F28848","#FCBC2A", "#F0F921"), ncol = 2) + NoLegend() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

dimplot_legend <- DimPlot(integrate.combined_0.3_mod, reduction = "umap",label = T, order = rev(c( "PLASM1","PLASM2","MET","AMP","CC","LAM1","LAM2","LAM3")), cols = c('#9970ab','#8c96c6','#99d8c9','#41ae76',"#DC3220","#F28848","#FCBC2A", "#F0F921"), repel = T, label.size = 5) + NoLegend() + NoAxes()
Noinflabel <- grobTree(rectGrob(gp=gpar(fill="NA", lwd = 0)),
                 textGrob("No Infection",  gp=gpar(col="black", fontsize=16, fontface="bold")))
Inflabel <- grobTree(rectGrob(gp=gpar(fill="NA", lwd = 0)),
                       textGrob("Infection", rot = 0, gp=gpar(col="black", fontsize=16, fontface="bold")))
Nosellabel <- grobTree(rectGrob(gp=gpar(fill="NA", lwd = 0)),
                       textGrob("No Selection",  rot = 90, gp=gpar(col="black", fontsize=16, fontface="bold")))
Sellabel <- grobTree(rectGrob(gp=gpar(fill="NA", lwd = 0)),
                     textGrob("Selection", rot = 90, gp=gpar(col="black", fontsize=16, fontface="bold")))
lay <- rbind(c("NA",2,2,2,2,2,2,3,3,3,3,3,3),
             c(4,1,1,1,1,1,1,1,1,1,1,1,1),
             c(4,1,1,1,1,1,1,1,1,1,1,1,1),
             c(4,1,1,1,1,1,1,1,1,1,1,1,1),
             c(4,1,1,1,1,1,1,1,1,1,1,1,1),
             c(4,1,1,1,1,1,1,1,1,1,1,1,1),
             c(4,1,1,1,1,1,1,1,1,1,1,1,1),
             c(5,1,1,1,1,1,1,1,1,1,1,1,1),
             c(5,1,1,1,1,1,1,1,1,1,1,1,1),
             c(5,1,1,1,1,1,1,1,1,1,1,1,1),
             c(5,1,1,1,1,1,1,1,1,1,1,1,1),
             c(5,1,1,1,1,1,1,1,1,1,1,1,1),
             c(5,1,1,1,1,1,1,1,1,1,1,1,1))
dimplot_data2 <- grid.arrange(grobs = list(dimplot_data,Noinflabel,Inflabel,Nosellabel,Sellabel), layout_matrix = lay)

lay2 <- rbind(c("NA","NA","NA","NA","NA",1,1,1,1,1,1,1,1),
             c("NA","NA","NA","NA","NA",1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,1,1,1,1,1,1,1,1),
             c("NA","NA","NA","NA","NA",1,1,1,1,1,1,1,1),
             c("NA","NA","NA","NA","NA",1,1,1,1,1,1,1,1))

svg("integrate.combined_0.3_mod_umap.svg", width = 9, height = 4.3)
grid.arrange(grobs = list(dimplot_data2,dimplot_legend), layout_matrix = lay2)
dev.off()

##compare to integrate.combined
cellscompare <- ""
clustersnames <- names(table(Idents(integrate.combined_0.3_mod)))
for (i in 1:(length(table(Idents(integrate.combined_0.3_mod))))){
  tempcells <- WhichCells(integrate.combined_0.3_mod, idents = clustersnames[i])
  tempcells_data <- subset(integrate.combined_0.3, cells = tempcells)
  celllist <- Idents(tempcells_data)
  celllist <- factor(celllist, levels = 0:(length(table(Idents(integrate.combined_0.3)))-1))
  cellscompare <- rbind(cellscompare,table(celllist))
  
}
cellscompare <- cellscompare[-1,]
rownames(cellscompare) <- clustersnames
write.csv(cellscompare, file = "row_all_mod_col_all_0.3.csv")

##percent mt
VlnPlot(integrate.combined_0.3_mod, features = "percent.mt")

##cluster phase
phase_table <- table(integrate.combined_0.3_mod$seurat_clusters,integrate.combined_0.3_mod$Phase)
phase_table <- phase_table/rowSums(phase_table)
phase_table <- melt(phase_table)
colnames(phase_table) <- c("Cluster","Phase","Proportion")
ggplot(data=phase_table, aes(x=factor(Cluster,levels = c( "PLASM1","PLASM2","MET","AMP","CC","LAM1","LAM2","LAM3")), y=Proportion, fill = Phase)) +
  geom_bar(stat="identity")+
  xlab("Cluster")+
  ylab("Proportion")

##cluster phase by poptreat
groups <- unique(integrate.combined_0.3_mod$poptreat)
phaseplots <- list()
for (g in 1:length(groups)){
  integrate.combined_0.3_mod_sub <- subset(integrate.combined_0.3_mod, poptreat == groups[g])
  phase_table <- table(integrate.combined_0.3_mod_sub$seurat_clusters,integrate.combined_0.3_mod_sub$Phase)
  phase_table <- phase_table/rowSums(phase_table)
  phase_table <- melt(phase_table)
  colnames(phase_table) <- c("Cluster","Phase","Proportion")
  phaseplots[[g]] <- ggplot(data=phase_table, aes(x=factor(Cluster,levels = c( "PLASM1","PLASM2","MET","AMP","CRY","LAM1","LAM2","LAM3")), y=Proportion, fill = Phase)) +
    geom_bar(stat="identity")+
    xlab("Cluster")+
    ylab("Proportion")+
    ggtitle(groups[g])
  
}
do.call("grid.arrange", c(phaseplots, ncol=2))

##cluster proportions 
integrate.combined_cellcount <- as.data.frame(table(paste(integrate.combined_0.3_mod$poptreat,integrate.combined_0.3_mod$seurat_clusters,sep="_")))
integrate.combined_cellcount <- integrate.combined_cellcount %>% separate(Var1, into= c("Treatment","Cluster"),sep="_")
propall <- vector()
for (t in levels(as.factor(integrate.combined_cellcount$Treatment))){
  hold <- integrate.combined_cellcount[integrate.combined_cellcount$Treatment == t,]
  prop <- hold$Freq/sum(hold$Freq)
  propall <- c(propall,prop)
}
integrate.combined_cellcount$Prop <- propall
integrate.combined_cellcount$Treatment <- factor(integrate.combined_cellcount$Treatment, levels = c("No Selection, No Infection", "No Selection, Infection", "Selection, No Infection", "Selection, Infection"))
integrate.combined_cellcount <- separate(integrate.combined_cellcount, Treatment, into = c("Selection", "Infection"), sep=", ", remove = F)
integrate.combined_cellcount$Selection <- as.factor(integrate.combined_cellcount$Selection)
integrate.combined_cellcount$Infection <- factor(integrate.combined_cellcount$Infection, levels = c("No Infection", "Infection"))

clusterplot <- ggplot(data=integrate.combined_cellcount, aes(x=Treatment, y=Prop, fill = Selection, color = Infection, size= Infection)) +
  geom_bar(stat="identity")+
  facet_wrap(~ factor(Cluster,levels = rev(c("LAM3","LAM2","LAM1","AMP","MET","CC","PLASM2","PLASM1"))), ncol=5,scales = "free", as.table = F) +
  xlab("")+
  ylab("Proportion of cells") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_size_manual(name = "", values = c(0.75,1.5)) +
  scale_color_manual(name = "Treatment", values = c("black", "red"))+ 
  scale_fill_manual(name = "Selection regime", values = c("#cccccc", "#5aa5e6"))+ 
  theme(legend.position = c(0.84, 1.06),
        legend.justification = c(0.84, 1.06)) + 
  guides(color = guide_legend(override.aes = list(fill = "white")), size = F) +
  theme(legend.margin = margin(0.4,0,0,0, unit="cm"))

pdf("clusprop.pdf", height = 7.8, width = 10)
plot_grid(clusterplot, lineage_plot, labels = c('A', ''), label_size = 12, ncol = 1, rel_heights = c(1.7, 1))
dev.off()

##for flymine pathway analysis
universe_sc <- read.csv(file = "universe_sc")
universe_sc_2000 <- rownames(integrate.combined_0.3_mod)
universe_sc_lm_2000 <- rownames(integrate.combined_lm_0.2)
write.table(universe_sc_2000, file = "universe_sc_2000.txt", quote = F, row.names = F, col.names = F)
write.table(universe_sc_lm_2000, file = "universe_sc_2000.txt", quote = F, row.names = F, col.names = F)
write.table(universe_sc_lm_2000, file = "universe_sc_lm_2000.txt", quote = F, row.names = F, col.names = F)

##cell cycle and plasmatocyte markers
DotPlot(integrate.combined_0.3_mod, features = c("FBgn0003124","FBgn0003525","FBgn0261385", "FBgn0259896", "FBgn0029167", "FBgn0243514", "FBgn0011828", "FBgn0000299")) +
  scale_x_discrete(labels=c("Col4a1","Pxn","eater","Hml","NimC1", "scra","stg","polo")) + 
  theme(axis.text.x = element_text(face = "italic"))+
  ylab("Cluster")+
  xlab("")

##lamellocyte markers
##lamellocyte markers
lam_markers <- c("FBgn0261363","FBgn0032422","FBgn0004657","FBgn0034005")
lam_markers_name <- c("ItgaPS4", "atilla", "mys", "PPO3")
lam_markers_plot <- list()

for (l in 1:length(lam_markers)){
  lam_markers_plot[[l]] <- VlnPlot(object = integrate.combined_0.3_mod, features = lam_markers[l], pt.size = 0) +
    xlab("") +
    ggtitle(lam_markers_name[l])+
    NoLegend()+
    theme(plot.title = element_text(size=14, face="bold.italic"))
}
do.call("grid.arrange", c(lam_markers_plot, ncol=2))


##integrate mod cluster markers
##plasm2 vs plasm1
sample.markers_res0.3_PLASM2v1 <- FindMarkers(integrate.combined_0.3_mod, ident.1 = "PLASM1", ident.2 = "PLASM2",min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
write.csv(sample.markers_res0.3_PLASM2v1, file = "sample.markers_res0.3_PLASM2v1.csv")
DotPlot(integrate.combined_0.3_mod, features = rownames(head(sample.markers_res0.3_PLASM2v1[sample.markers_res0.3_PLASM2v1$avg_logFC <0,])))

##cc vs plasm1
sample.markers_res0.3_CCvPLASM1 <- FindMarkers(integrate.combined_0.3_mod, ident.1 = "PLASM1", ident.2 = "CC",min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
write.csv(sample.markers_res0.3_CCvPLASM1, file = "sample.markers_res0.3_CCvPLASM1.csv")

##met vs plasm1
sample.markers_res0.3_METvPLASM1 <- FindMarkers(integrate.combined_0.3_mod, ident.1 = "PLASM1", ident.2 = "MET",min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
write.csv(sample.markers_res0.3_METvPLASM1, file = "sample.markers_res0.3_METvPLASM1.csv")

##amp vs plasm1
sample.markers_res0.3_AMPvPLASM1 <- FindMarkers(integrate.combined_0.3_mod, ident.1 = "PLASM1", ident.2 = "AMP",min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
write.csv(sample.markers_res0.3_AMPvPLASM1, file = "sample.markers_res0.3_AMPvPLASM1.csv")

##combine find markers
cluster_markers <- ""
tempint <- sample.markers_res0.3_PLASM2v1[sample.markers_res0.3_PLASM2v1$avg_logFC >0,] 
tempint$cluster <- "PLASM1"
cluster_markers <- rbind(cluster_markers,tempint)
cluster_markers <- cluster_markers[-c(1),]
tempint <- sample.markers_res0.3_PLASM2v1[sample.markers_res0.3_PLASM2v1$avg_logFC <0,] 
tempint$cluster <- "PLASM2"
tempint$avg_logFC <- -tempint$avg_logFC
cluster_markers <- rbind(cluster_markers,tempint)
tempint <- sample.markers_res0.3_CCvPLASM1[sample.markers_res0.3_CCvPLASM1$avg_logFC <0,] 
tempint$cluster <- "CC"
tempint$avg_logFC <- -tempint$avg_logFC
cluster_markers <- rbind(cluster_markers,tempint)
tempint <- sample.markers_res0.3_METvPLASM1[sample.markers_res0.3_METvPLASM1$avg_logFC <0,] 
tempint$cluster <- "MET"
tempint$avg_logFC <- -tempint$avg_logFC
cluster_markers <- rbind(cluster_markers,tempint)
tempint <- sample.markers_res0.3_AMPvPLASM1[sample.markers_res0.3_AMPvPLASM1$avg_logFC <0,] 
tempint$cluster <- "AMP"
tempint$avg_logFC <- -tempint$avg_logFC
cluster_markers <- rbind(cluster_markers,tempint)
cluster_markers$gene <- rownames(cluster_markers)
cluster_markers$gene <- substr(cluster_markers$gene, start = 1, stop = 11)
cluster_markers$symbol <- mapIds(org.Dm.eg.db, as.character(cluster_markers$gene), column="SYMBOL", keytype="FLYBASE", multiVals="first")
cluster_markers_top <- cluster_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
write.csv(cluster_markers, file = "cluster_markers.csv")
DotPlot(integrate.combined_0.3_mod, features = cluster_markers_top$gene, slot = "data")+
  scale_x_discrete(labels=rev(cluster_markers_top$symbol)) + 
  theme(axis.text.x = element_text(face = "italic"))


##save count matrix

integrated_rna_count_matrix <- GetAssayData(object = integrate.combined_0.3_mod@assays$RNA, slot = 'counts')
integrated_rna_count_matrix <- as.data.frame(t(integrated_rna_count_matrix))
integrated_rna_count_matrix$seurat_clusters <- as.character(integrate.combined_0.3_mod$seurat_clusters)
write.csv(integrated_rna_count_matrix, file = "integrated_rna_count_matrix.csv")

integrated_scaled_hvg_matrix <-  GetAssayData(integrate.combined_0.3_mod@assays$integrated, slot = 'scale.data')
integrated_scaled_hvg_matrix <- as.data.frame(t(integrated_scaled_hvg_matrix))
integrated_scaled_hvg_matrix$seurat_clusters <- as.character(integrate.combined_0.3_mod$seurat_clusters)
write.csv(integrated_scaled_hvg_matrix, file = "integrated_scaled_hvg_matrix.csv")

integrated_subcluster_scaled_hvg_matrix <-  GetAssayData(integrate.combined_lm_0.2@assays$integrated, slot = 'scale.data')
integrated_subcluster_scaled_hvg_matrix <- as.data.frame(t(integrated_subcluster_scaled_hvg_matrix))
integrated_subcluster_scaled_hvg_matrix$seurat_clusters <- as.character(integrate.combined_lm_0.2$seurat_clusters)
write.csv(integrated_subcluster_scaled_hvg_matrix, file = "integrated_subcluster_scaled_hvg_matrix.csv")

integrated_hemocytes <- integrate.combined_0.3_mod
integrated_lamellocyte_subcluster <- integrate.combined_lm_0.2
save(integrated_hemocytes, file = "integrated_hemocytes.Robj")
save(integrated_lamellocyte_subcluster, file = "integrated_lamellocyte_subcluster.Robj")
