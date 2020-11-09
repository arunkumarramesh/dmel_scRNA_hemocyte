library(grid)
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
library(tidymodels)
library(ranger)
library(viridis)
library(plyr)
library(rgl)
library(LSD)
library(patchwork)
library(RColorBrewer)
library(limma)

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

integrate.combined_0.3 <- FindClusters(integrate.combined, resolution = 0.3)
integrate.combined_0.3$poptreat <- gsub("control Infected", "No Parasitism, Infection", integrate.combined_0.3$poptreat)
integrate.combined_0.3$poptreat <- gsub("control Uninf", "No Parasitism, No Infection", integrate.combined_0.3$poptreat)
integrate.combined_0.3$poptreat <- gsub("NSRef Infected", "High Parasitism, Infection", integrate.combined_0.3$poptreat)
integrate.combined_0.3$poptreat <- gsub("NSRef Uninf", "High Parasitism, No Infection", integrate.combined_0.3$poptreat)
integrate.combined_0.3$poptreat <- factor(integrate.combined_0.3$poptreat, levels = c("No Parasitism, No Infection", "No Parasitism, Infection", "High Parasitism, No Infection", "High Parasitism, Infection"))
integrate.combined_0.3$population <- gsub("control", "No Parasitism", integrate.combined_0.3$population)
integrate.combined_0.3$population <- gsub("NSRef", "High Parasitism", integrate.combined_0.3$population)
integrate.combined_0.3$population <- factor(integrate.combined_0.3$population, levels = c("No Parasitism","High Parasitism"))
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
lam_markers <- c("FBgn0034005","FBgn0032422","FBgn0004657","FBgn0261363")
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


####SECTION 6 - ROUND 3 CLUSTERING WITH CELL CYCLE CORRECTION ########
integrate.combined_cellcycle <- ScaleData(integrate.combined, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt", "CC.Difference"), features = rownames(integrate.combined))
integrate.combined_cellcycle <- RunPCA(integrate.combined_cellcycle, npcs = 50, verbose = FALSE)
integrate.combined_cellcycle <- RunUMAP(integrate.combined_cellcycle, reduction = "pca", dims = 1:50)
integrate.combined_cellcycle <- FindNeighbors(integrate.combined_cellcycle, reduction = "pca", dims = 1:50)
integrate.combined_cellcycle_0.3 <- FindClusters(integrate.combined_cellcycle, resolution = 0.3)
integrate.combined_cellcycle_0.3$poptreat <- gsub("control Infected", "No Parasitism, Infection", integrate.combined_cellcycle_0.3$poptreat)
integrate.combined_cellcycle_0.3$poptreat <- gsub("control Uninf", "No Parasitism, No Infection", integrate.combined_cellcycle_0.3$poptreat)
integrate.combined_cellcycle_0.3$poptreat <- gsub("NSRef Infected", "High Parasitism, Infection", integrate.combined_cellcycle_0.3$poptreat)
integrate.combined_cellcycle_0.3$poptreat <- gsub("NSRef Uninf", "High Parasitism, No Infection", integrate.combined_cellcycle_0.3$poptreat)
integrate.combined_cellcycle_0.3$poptreat <- factor(integrate.combined_cellcycle_0.3$poptreat, levels = c("No Parasitism, No Infection", "No Parasitism, Infection", "High Parasitism, No Infection", "High Parasitism, Infection"))
integrate.combined_cellcycle_0.3$population <- gsub("control", "No Parasitism", integrate.combined_cellcycle_0.3$population)
integrate.combined_cellcycle_0.3$population <- gsub("NSRef", "High Parasitism", integrate.combined_cellcycle_0.3$population)
integrate.combined_cellcycle_0.3$population <- factor(integrate.combined_cellcycle_0.3$population, levels = c("No Parasitism","High Parasitism"))
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
write.csv(cellscompare, file = "row_cc_0.3_col_all_0.3.csv", quote = F)

####SECTION 7 - DATA INTEGRATION FOR ESTIMATING MEAN FOLD CHANGE ACROSS DETECTED GENES########

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

#write table of logFC
logfc_allgenes <- cbind(infvuninfec_15000f$gene,selvinfection_15000f$avg_logFC.x,selvinfection_15000f$avg_logFC.y,infvuninfec_15000f$avg_logFC.y,selinfvnoseluninfec_15000f$avg_logFC.y)
colnames(logfc_allgenes) <- c("Gene","Infection only FC","High Parasitism only FC","Infection FC in selected populations", "High Parasitism and Infection FC")
write.csv(logfc_allgenes, file = "logfc_allgenes.csv")

##Figure 2 plots

logfc_allgenes <- read.csv(file = "logfc_allgenes2.csv", row.names = 1)

d=data.frame(x1=c(1,1,2,2), x2=c(2,2,3,3), y1=c(1,2,1,2), y2=c(2,3,2,3), c=c('a','a','b','b'), t=c('a','b','a','b'), r=c("2: Constitutive\n (Evolved)","","3: Constitutive +\n Induced","1: Induced"))
rec <- ggplot() + 
  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color = "black", fill = "white", size=0.5, linetype = 3) +
  geom_text(data=d, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=r), size=4) +
  geom_text(data=d, aes(label = "No Infection", x = 1.5, y = 3.2), size = 4) +
  geom_text(data=d, aes(label = "Infection", x = 2.5, y = 3.2), size = 4) +
  geom_text(data=d, aes(label = "No Parasitism", x = 0.8, y = 2.5, angle = 90,), size = 4) +
  geom_text(data=d, aes(label = "High Parasitism", x = 0.8, y = 1.5, angle = 90), size = 4) +
  geom_text(data=d, aes(label = "Selection regime", x = 0.6, y = 2, angle = 90, fontface  = "bold"), size = 4) +
  geom_text(data=d, aes(label = "Treatment", x = 2, y = 3.4, fontface  = "bold"), size = 4) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+ 
  theme(legend.position = "none")

##atilla and ppo3

pdf("bulkgeneexp_mod_5.pdf",width = 12, height = 4.4)

par(mfrow=c(1,3))
par(mar=c(5,6,4,1)+0.8)
#par(mgp=c(4,1,0))
plot.new()
title("A", adj = 0, cex.main=2)
heatscatter(logfc_allgenes$Infection.only.FC,logfc_allgenes$Selection.only.FC, ylab =expression(' 2: Constitutive Expression (log2 FC)\n     (Evolved with High Parasitism)'),xlab ="", main = "", bty="L", cex.lab=1.6, cex.axis=1.5, xlim= c(-0.9,1.5), ylim=c(-1.1,1.3))
abline(coef = c(0,1))
points(x = logfc_allgenes[logfc_allgenes$Gene == "FBgn0032422",]$Infection.only.FC,
       y = logfc_allgenes[logfc_allgenes$Gene == "FBgn0032422",]$Selection.only.FC,
       pch = 16, col = "black")
text("atilla", font = 3, x = 0.18+logfc_allgenes[logfc_allgenes$Gene == "FBgn0032422",]$Infection.only.FC,
     y = logfc_allgenes[logfc_allgenes$Gene == "FBgn0032422",]$Selection.only.FC, cex = 1.3)

points(x = logfc_allgenes[logfc_allgenes$Gene == "FBgn0261363",]$Infection.only.FC,
       y = logfc_allgenes[logfc_allgenes$Gene == "FBgn0261363",]$Selection.only.FC,
       pch = 16, col = "black")
text("PPO3", font = 3, x = -0.04+logfc_allgenes[logfc_allgenes$Gene == "FBgn0261363",]$Infection.only.FC,
     y = 0.12+logfc_allgenes[logfc_allgenes$Gene == "FBgn0261363",]$Selection.only.FC, cex = 1.3)
points(x = logfc_allgenes[logfc_allgenes$Gene == "FBgn0003748",]$Infection.only.FC,
       y = logfc_allgenes[logfc_allgenes$Gene == "FBgn0003748",]$Selection.only.FC,
       pch = 16, col = "black")
text("Treh", font = 3, x = -0.18+logfc_allgenes[logfc_allgenes$Gene == "FBgn0003748",]$Infection.only.FC,
     y = logfc_allgenes[logfc_allgenes$Gene == "FBgn0003748",]$Selection.only.FC, cex = 1.3)

points(x = logfc_allgenes[logfc_allgenes$Gene == "FBgn0050035",]$Infection.only.FC,
       y = logfc_allgenes[logfc_allgenes$Gene == "FBgn0050035",]$Selection.only.FC,
       pch = 16, col = "black")
text("Tret1-1", font = 3, x = -0.07+logfc_allgenes[logfc_allgenes$Gene == "FBgn0050035",]$Infection.only.FC,
     y = 0.08+logfc_allgenes[logfc_allgenes$Gene == "FBgn0050035",]$Selection.only.FC, cex = 1.3)
title("B", adj = 0, cex.main=2)

heatscatter(logfc_allgenes$Infection.only.FC,logfc_allgenes$Selection.and.Infection.FC, ylab =expression('   3: Constitutive + Induced Expression (log2 FC)\n                  (Evolved with High Parasitism)'),xlab ="", main = "", bty="L", cex.lab=1.6, cex.axis=1.5, xlim= c(-0.9,1.58), ylim=c(-1.1,1.3))
abline(coef = c(0,1))
points(x = logfc_allgenes[logfc_allgenes$Gene == "FBgn0032422",]$Infection.only.FC,
       y = logfc_allgenes[logfc_allgenes$Gene == "FBgn0032422",]$Selection.and.Infection.FC,
       pch = 16, col = "black")
text("atilla", font = 3, x = -0.18+logfc_allgenes[logfc_allgenes$Gene == "FBgn0032422",]$Infection.only.FC,
     y = logfc_allgenes[logfc_allgenes$Gene == "FBgn0032422",]$Selection.and.Infection.FC, cex = 1.3)
points(x = logfc_allgenes[logfc_allgenes$Gene == "FBgn0261363",]$Infection.only.FC,
       y = logfc_allgenes[logfc_allgenes$Gene == "FBgn0261363",]$Selection.and.Infection.FC,
       pch = 16, col = "black")
text("PPO3", font = 3, x = -0.1+logfc_allgenes[logfc_allgenes$Gene == "FBgn0261363",]$Infection.only.FC,
     y = 0.1+logfc_allgenes[logfc_allgenes$Gene == "FBgn0261363",]$Selection.and.Infection.FC, cex = 1.3)
points(x = logfc_allgenes[logfc_allgenes$Gene == "FBgn0003748",]$Infection.only.FC,
       y = logfc_allgenes[logfc_allgenes$Gene == "FBgn0003748",]$Selection.and.Infection.FC,
       pch = 16, col = "black")
text("Treh", font = 3, x = 0.13+logfc_allgenes[logfc_allgenes$Gene == "FBgn0003748",]$Infection.only.FC,
     y = -0.07+logfc_allgenes[logfc_allgenes$Gene == "FBgn0003748",]$Selection.and.Infection.FC, cex = 1.3)
points(x = logfc_allgenes[logfc_allgenes$Gene == "FBgn0050035",]$Infection.only.FC,
       y = logfc_allgenes[logfc_allgenes$Gene == "FBgn0050035",]$Selection.and.Infection.FC,
       pch = 16, col = "black")
text("Tret1-1", font = 3, x = 0.23+logfc_allgenes[logfc_allgenes$Gene == "FBgn0050035",]$Infection.only.FC,
     y = logfc_allgenes[logfc_allgenes$Gene == "FBgn0050035",]$Selection.and.Infection.FC, cex = 1.3)
title("C", adj = 0, cex.main=2)
mtext("1: Induced expression in populations evolved without parasitism (log2 fold change)",
      side=1,outer=T,line=-2.2,cex=1.1, at = 0.68)

vp.Left <- viewport(height=grid::unit(0.9, "npc"), width=grid::unit(0.3, "npc"), 
                    just=c("left","top"), 
                    y=0.95, x=0.01)

# plot the ggplot using the print command
print(rec, vp=vp.Left)
dev.off()

####SECTION 8 - CLUSTERING WITH ALL DETECTED ANCHORS########

integrate.combined_15000f_0.3 <- FindClusters(integrate.combined_15000f, resolution = 0.3)
integrate.combined_15000f_0.3$poptreat <- gsub("control Infected", "No Parasitism, Infection", integrate.combined_15000f_0.3$poptreat)
integrate.combined_15000f_0.3$poptreat <- gsub("control Uninf", "No Parasitism, No Infection", integrate.combined_15000f_0.3$poptreat)
integrate.combined_15000f_0.3$poptreat <- gsub("NSRef Infected", "High Parasitism, Infection", integrate.combined_15000f_0.3$poptreat)
integrate.combined_15000f_0.3$poptreat <- gsub("NSRef Uninf", "High Parasitism, No Infection", integrate.combined_15000f_0.3$poptreat)
integrate.combined_15000f_0.3$poptreat <- factor(integrate.combined_15000f_0.3$poptreat, levels = c("No Parasitism, No Infection", "No Parasitism, Infection", "High Parasitism, No Infection", "High Parasitism, Infection"))
integrate.combined_15000f_0.3$population <- gsub("control", "No Parasitism", integrate.combined_15000f_0.3$population)
integrate.combined_15000f_0.3$population <- gsub("NSRef", "High Parasitism", integrate.combined_15000f_0.3$population)
integrate.combined_15000f_0.3$population <- factor(integrate.combined_15000f_0.3$population, levels = c("No Parasitism","High Parasitism"))
integrate.combined_15000f_0.3$treatment <- gsub("Infected", "Infection", integrate.combined_15000f_0.3$treatment)
integrate.combined_15000f_0.3$treatment <- gsub("Uninf", "No Infection", integrate.combined_15000f_0.3$treatment)
integrate.combined_15000f_0.3$treatment <- factor(integrate.combined_15000f_0.3$treatment, levels = c("No Infection","Infection"))

integrate.combined_cellcount <- as.data.frame(table(paste(integrate.combined_15000f_0.3$poptreat,integrate.combined_15000f_0.3$seurat_clusters,sep="_")))
integrate.combined_cellcount <- integrate.combined_cellcount %>% separate(Var1, into= c("Treatment","Cluster"),sep="_")
propall <- vector()
for (t in levels(as.factor(integrate.combined_cellcount$Treatment))){
  hold <- integrate.combined_cellcount[integrate.combined_cellcount$Treatment == t,]
  prop <- hold$Freq/sum(hold$Freq)
  propall <- c(propall,prop)
}
integrate.combined_cellcount$Prop <- propall
integrate.combined_cellcount$Treatment <- factor(integrate.combined_cellcount$Treatment, levels = c("No Parasitism, No Infection", "No Parasitism, Infection", "High Parasitism, No Infection", "High Parasitism, Infection"))
integrate.combined_cellcount <- separate(integrate.combined_cellcount, Treatment, into = c("Parasitism", "Infection"), sep=", ", remove = F)
integrate.combined_cellcount$Parasitism <- as.factor(integrate.combined_cellcount$Parasitism)
integrate.combined_cellcount$Infection <- factor(integrate.combined_cellcount$Infection, levels = c("No Infection", "Infection"))

clusterplot15000 <- ggplot(data=integrate.combined_cellcount, aes(x=Treatment, y=Prop, fill = Parasitism, alpha = Infection)) +
  geom_bar(stat="identity")+
  facet_wrap(~ as.numeric(Cluster), ncol=4,scales = "free")+  
  xlab("")+
  ylab("Proportion of cells") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_alpha_manual(name = "", values = c(0.5,1),  guide = guide_legend(override.aes = list(fill = c("#0072b2", "#d55e00"), alpha = 0.5)))+ 
  scale_fill_manual(name = "", values = c("#0072b2", "#d55e00")) +
  theme(legend.position="none")

d=data.frame(x1=c(2.3,3.2,2.3,3.2), x2=c(2.5,3.4,2.5,3.4), y1=c(2,2,3,3), y2=c(2.8,2.8,3.8,3.8), t=c('a','b','a','b'), r=c('c','c','d','d'),s=c('','','No \nParasitism','High \nParasitism'),u=c('Infection','','No Infection',''))
exptlegend <- ggplot() + 
  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t,alpha=r))+
  scale_alpha_manual(name = "", values = c(1,0.5))+
  scale_fill_manual(name = "", values = c("#0072b2","#d55e00"))+
  theme_void()+
  theme(legend.position="none")+
  geom_text(data=d, aes(x=x1+(x2-x1)/2, y=(y1+(y2-y1)/2)+1.5, label=s), size=4)+
  geom_text(data=d, aes(x=(x1+(x2-x1)/2)-0.8, y=y1+(y2-y1)/2, label=u), size=4)+
  xlim(0.8,4)+
  ylim(0.8,6)+
  ggtitle("Selection Regime")+
  theme(plot.title = element_text(hjust = 0.72, vjust = 0, face="bold"))

gnew15000 <- ggdraw(clusterplot15000) + draw_plot(exptlegend,x = 0.72, y = 0, width = 0.3, height = 0.25)

####SECTION 9 - ROUND 1 LAMELLOCYTE  SUBCLUSTERING, MOVING MISCLASSIFIED CRYSTAL CELLS ########

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

####SECTION 10 - ROUND 2 LAMELLOCYTE  SUBCLUSTERING ########

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
integrate.combined_lm_0.2$poptreat <- gsub("control Infected", "No Parasitism, Infection", integrate.combined_lm_0.2$poptreat)
integrate.combined_lm_0.2$poptreat <- gsub("control Uninf", "No Parasitism, No Infection", integrate.combined_lm_0.2$poptreat)
integrate.combined_lm_0.2$poptreat <- gsub("NSRef Infected", "High Parasitism, Infection", integrate.combined_lm_0.2$poptreat)
integrate.combined_lm_0.2$poptreat <- gsub("NSRef Uninf", "High Parasitism, No Infection", integrate.combined_lm_0.2$poptreat)
integrate.combined_lm_0.2$poptreat <- factor(integrate.combined_lm_0.2$poptreat, levels = c("No Parasitism, No Infection", "No Parasitism, Infection", "High Parasitism, No Infection", "High Parasitism, Infection"))
integrate.combined_lm_0.2$population <- gsub("control", "No Parasitism", integrate.combined_lm_0.2$population)
integrate.combined_lm_0.2$population <- gsub("NSRef", "High Parasitism", integrate.combined_lm_0.2$population)
integrate.combined_lm_0.2$population <- factor(integrate.combined_lm_0.2$population, levels = c("No Parasitism","High Parasitism"))
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
  scale_x_discrete(labels=rev(c("PPO3", "atilla", "mys", "ItgaPS4")))

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
lamphase <- ggplot(data=phase_table, aes(x=factor(Cluster,levels = c( "PLASM1","LAM1","LAM2","LAM3")), y=Proportion, fill = Phase)) +
  geom_bar(stat="identity")+
  xlab("Cluster")+
  ylab("Proportion")

####SECTION 11 - TRAJECTORY INFERENCE FOR LAMELLOCYTE DIFFERENTIATION ########

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

model <- rand_forest(trees = 2000, min_n = 15, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 32) %>%
  fit(pseudotime ~ ., data = dat_use_df)
impt_pval <- importance_pvalues(model$fit, method = "altmann",formula = pseudotime ~ ., data = dat_use_df)
impt_pval_df <- as.data.frame(impt_pval)
impt_pval_df$symbol <- mapIds(org.Dm.eg.db, rownames(impt_pval_df), column="SYMBOL", keytype="FLYBASE", multiVals="first")
impt_pval_df <- read.csv(file = "impt_pval_df.csv", row.names = 1)
impt_pval_df <- impt_pval_df[order(rownames(impt_pval_df)),]
aveexp <- read.csv(file = "aveexp.csv", row.names = 1)
aveexp <- aveexp[order(rownames(aveexp)),]
impt_pval_df <- cbind(impt_pval_df[1:2],impt_pval_df[4],aveexp[1:4])
impt_pval_df <- impt_pval_df[order(impt_pval_df$pvalue,-impt_pval_df$importance),]
write.csv(impt_pval_df, file = "var_imp_df_lm_lineage_2000hvg.csv")

##dimplot lm subcluster
lm_dimplot <- DimPlot(object = integrate.combined_lm_0.2, reduction = "umap",label = F, order = rev(c( "PLASM1","LAM1","LAM2","LAM3")), cols = c('#9970ab',"#F28848","#FCBC2A", "#F0F921"))
lm_dimplot_data <- lm_dimplot$data
lm_dimplot_data_mean <- aggregate(lm_dimplot_data[, 1:2], list(lm_dimplot_data$ident), mean)

lm_dimplot_lineage <- lm_dimplot +
  geom_line(data = lm_dimplot_data_mean, aes(UMAP_1,UMAP_2), size = 1, arrow = arrow())+
  #geom_point(data = lm_dimplot_data_mean, aes(UMAP_1,UMAP_2), size = 2)+
  geom_line(data = lm_dimplot_data_mean[1:3,], aes(UMAP_1,UMAP_2), size = 1, arrow = arrow())+
  geom_line(data = lm_dimplot_data_mean[1:2,], aes(UMAP_1,UMAP_2), size = 1, arrow = arrow())+
  #geom_line(data = lm_dimplot_data_mean[1,], aes(UMAP_1,UMAP_2), size = 1, arrow = arrow())+
  theme(legend.position="right")

#top lineage marker plot
integrate.combined_lm_0.2_forplot = SetAssayData(object = integrate.combined_lm_0.2, slot = "data", assay = "integrated", new.data = log2(exp(as.matrix(GetAssayData(object = integrate.combined_lm_0.2, slot = "data", assay = "integrated")))))
integrate.combined_lm_0.2_ave <- integrate.combined_lm_0.2_forplot
integrate.combined_lm_0.2_ave$seurat_clusters <- gsub("LAM3","LAM2",integrate.combined_lm_0.2_ave$seurat_clusters)
Idents(integrate.combined_lm_0.2_ave) <- integrate.combined_lm_0.2_ave$seurat_clusters
aveexp <- AverageExpression(object = integrate.combined_lm_0.2_ave, features = head(as.character(rownames(impt_pval_df)), n = 100), assays = "integrated", slot = "scale.data")
aveexp <- aveexp$integrated
aveexp <- aveexp[order(-aveexp$LAM2),]
scaleddata_forplot <- integrate.combined_lm_0.2_forplot@assays[["integrated"]]@scale.data
scaleddata_forplot <- scaleddata_forplot[,order(scaleddata_forplot[which(rownames(scaleddata_forplot) %in% rownames(impt_pval_df)[101]),])]

top100_part1<- DoHeatmap2(integrate.combined_lm_0.2_forplot, features = head(as.character(rownames(aveexp)), n = 50), angle = 0, size = 6.5, label = T, group.bar = T, draw.lines = T, group.by = "ident", cells = colnames(scaleddata_forplot),  lines.width = 100, group.bar.height = 0, hjust = 0.5, slot = "scale.data") +
  #ggtitle("Top 100 lineage markers") +
  theme(axis.text.y = element_text(face="italic")) +
  guides(color=FALSE) +
  guides(fill=FALSE) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=14))+
  scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),na.value = "white")
top100_part1$data$Feature<- mapvalues(top100_part1$data$Feature, from = head(as.character(rownames(aveexp)), n = 50), to = mapIds(org.Dm.eg.db, head(as.character(rownames(aveexp)), n = 50), column="SYMBOL", keytype="FLYBASE", multiVals="first"))

top100_part2<- DoHeatmap2(integrate.combined_lm_0.2_forplot, features = tail(as.character(rownames(aveexp)), n = 50), angle = 0, size = 6.5, label = T, group.bar = T, draw.lines = T, group.by = "ident", cells = colnames(scaleddata_forplot),  lines.width = 100, group.bar.height = 0, hjust = 0.5, slot = "scale.data") +
  #ggtitle("Top 50 lineage markers") +
  theme(axis.text.y = element_text(face="italic")) +
  guides(color=FALSE) +
  guides(fill=guide_legend(title="Scaled \nLog (2) expression"))+
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=14))+
  scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),na.value = "white")
top100_part2$data$Feature<- mapvalues(top100_part2$data$Feature, from = tail(as.character(rownames(aveexp)), n = 50), to = mapIds(org.Dm.eg.db, tail(as.character(rownames(aveexp)), n = 50), column="SYMBOL", keytype="FLYBASE", multiVals="first"))
title.grob <- textGrob("Top 100 lineage markers",  gp=gpar(fontface="bold", col="black", fontsize=18))

top100markersplot <- grid.arrange(grobs = list(top100_part1,top100_part2), layout_matrix = t(matrix(c(1,1,1,1,1,2,2,2,2,2,2,2))), top = title.grob)

pdf(file = "top100lineagemarkers.pdf", height = 12, width = 14)
grid.arrange(grobs = list(top100_part1,top100_part2), layout_matrix = t(matrix(c(1,1,1,1,1,2,2,2,2,2,2,2))), top = title.grob)
dev.off()

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

lineage_plot <- plot_grid(lm_dimplot_lineage,atilla_vlnplot,ppo3_vlnplot, labels = c('A','B','C'), label_size = 14, ncol = 3, rel_widths = c(1.8,1,1))

write.csv(atilla_vlnplot$data, file = "atilla_sourcedata.csv", quote = F)
write.csv(ppo3_vlnplot$data, file = "ppo3_sourcedata.csv", quote = F)

##lineage pca
lm_dimplot_pca_legend <- DimPlot(object = integrate.combined_lm_0.2, reduction = "pca",label = T, order = rev(c( "PLASM1","LAM1","LAM2","LAM3")), cols = c('#9970ab',"#F28848","#FCBC2A", "#F0F921")) + NoLegend()
lm_dimplot_pca <- DimPlot(object = integrate.combined_lm_0.2, reduction = "pca",label = F, order = rev(c( "PLASM1","LAM1","LAM2","LAM3")), cols = c('#9970ab',"#F28848","#FCBC2A", "#F0F921"))

lm_dimplot_pca_data <- lm_dimplot_pca$data
lm_dimplot_pca_data_mean <- aggregate(lm_dimplot_pca_data[, 1:2], list(lm_dimplot_pca_data$ident), mean)

lm_dimplot_pca_lineage <- lm_dimplot_pca +
  geom_line(data = lm_dimplot_pca_data_mean, aes(PC_1,PC_2), size = 1)+
  geom_point(data = lm_dimplot_pca_data_mean, aes(PC_1,PC_2), size = 2)

##pathway analysis first with last LAM
sample_markers_lm_0.2_LAM1v3 <- FindMarkers(integrate.combined_lm_0.2, ident.1 = "LAM1", ident.2 = "LAM3",min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
sample_markers_lm_0.2_LAM1v3$gene <- rownames(sample_markers_lm_0.2_LAM1v3)
sample_markers_lm_0.2_LAM1v3$symbol <- mapIds(org.Dm.eg.db, as.character(sample_markers_lm_0.2_LAM1v3$gene), column="SYMBOL", keytype="FLYBASE", multiVals="first")
write.csv(sample_markers_lm_0.2_LAM1v3, file = "markers_LAM1v3.csv")


##GO plots for genes enriched in PLASM1 and LAM3
#plots only for Log(e)>1, FDR corrected and REVIGO reduction
#supramolecular fiber organization and muscle attachment removed due to complete gene overlap with others
integrate.combined_lm_0.2_forplot = SetAssayData(object = integrate.combined_lm_0.2, slot = "scale.data", assay = "integrated", new.data = log2(exp(as.matrix(GetAssayData(object = integrate.combined_lm_0.2, slot = "scale.data", assay = "integrated")))))

LAM3vPLASM1 <- FindMarkers(integrate.combined_lm_0.2, ident.1 = "LAM3", ident.2 = "PLASM1",min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
LAM3vPLASM1$gene <- mapIds(org.Dm.eg.db, as.character(row.names(LAM3vPLASM1)), column="SYMBOL", keytype="FLYBASE", multiVals="first")
LAM3vPLASM1$avg_logFC <- log2(exp(LAM3vPLASM1$avg_logFC))
write.csv(LAM3vPLASM1, file = "LAM3vPLASM1.csv")

universelist <- mapIds(org.Dm.eg.db, as.character(VariableFeatures(integrate.combined_lm_0.2)), column="ENTREZID", keytype="FLYBASE", multiVals="first")
genelist <- mapIds(org.Dm.eg.db, rownames(LAM3vPLASM1[LAM3vPLASM1$avg_logFC > 1.4,]), column="ENTREZID", keytype="FLYBASE", multiVals="first")
go.fisher <-  goana(genelist, universe = universelist, species = "Dm", prior.prob = NULL, covariate=NULL,plot = TRUE)
go.fisher$ADJ<- p.adjust(go.fisher$P.DE, method = "fdr")
go.fisher.sig <- subset(go.fisher, go.fisher$ADJ < 0.05)

for (r in 1:length(rownames(go.fisher.sig))){
  x <- org.Dm.egGO2ALLEGS
  Rkeys(x) <- rownames(go.fisher.sig)[r]
  EG <- mappedLkeys(x)
  genes <- paste(mapIds(org.Dm.eg.db, intersect(EG, genelist), column="FLYBASE", keytype="ENTREZID", multiVals="first"),collapse=" ")
  go.fisher.sig$Genes[r] <- genes
}

write.csv(go.fisher.sig,file = "go.fisher.sig_uplam3.csv")

scaleddata_forplot <- integrate.combined_lm_0.2_forplot@assays[["integrated"]]@data
scaleddata_forplot <- scaleddata_forplot[,order(scaleddata_forplot[which(rownames(scaleddata_forplot) %in% "FBgn0032422"),])]

x <- org.Dm.egGO2ALLEGS
Rkeys(x) <- "GO:0097435"
EG <- mappedLkeys(x)
aveexp <- AverageExpression(object = integrate.combined_lm_0.2_forplot, features = mapIds(org.Dm.eg.db, intersect(EG, genelist), column="FLYBASE", keytype="ENTREZID", multiVals="first"), assays = "integrated", slot = "scale.data")
aveexp <- aveexp$integrated
aveexp <- aveexp[order(-aveexp$LAM3),]
go0097435<- DoHeatmap2(integrate.combined_lm_0.2_forplot, features = rownames(aveexp), angle = 0, size = 5.5, label = T, group.bar = T, draw.lines = T, group.by = "ident", cells = colnames(scaleddata_forplot),  lines.width = 100, group.bar.height = 0, hjust = 0.5, slot = "scale.data") +
  ggtitle("Supramolecular fiber organization") +
  theme(axis.text.y = element_text(face="italic")) +
  guides(color=FALSE) +
  guides(fill=FALSE) +
  #guides(fill=guide_legend(title="Log (2) expression"))+
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=14))+
  #scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),na.value = "white")
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),na.value = "white")
go0097435$data$Feature <- mapvalues(go0097435$data$Feature, from =  rownames(aveexp), to = mapIds(org.Dm.eg.db, rownames(aveexp), column="SYMBOL", keytype="FLYBASE", multiVals="first"))

x <- org.Dm.egGO2ALLEGS
Rkeys(x) <- "GO:0030029"
EG <- mappedLkeys(x)
aveexp <- AverageExpression(object = integrate.combined_lm_0.2_forplot, features = mapIds(org.Dm.eg.db, intersect(EG, genelist), column="FLYBASE", keytype="ENTREZID", multiVals="first"), assays = "integrated", slot = "scale.data")
aveexp <- aveexp$integrated
aveexp <- aveexp[order(-aveexp$LAM3),]
go0030029<- DoHeatmap2(integrate.combined_lm_0.2_forplot, features = rownames(aveexp), angle = 0, size = 5.5, label = T, group.bar = T, draw.lines = T, group.by = "ident", cells = colnames(scaleddata_forplot),  lines.width = 100, group.bar.height = 0, hjust = 0.5, slot = "scale.data") +
  ggtitle("Actin filament-based process") +
  theme(axis.text.y = element_text(face="italic")) +
  guides(color=FALSE) +
  guides(fill=FALSE) +
  guides(fill=guide_legend(title="         Scaled \nLog (2) expression"))+
  theme(legend.position="bottom")+
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=14))+
  #scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),na.value = "white")+
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),na.value = "white")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
 go0030029$data$Feature <- mapvalues(go0030029$data$Feature, from =  rownames(aveexp), to = mapIds(org.Dm.eg.db, rownames(aveexp), column="SYMBOL", keytype="FLYBASE", multiVals="first"))

x <- org.Dm.egGO2ALLEGS
Rkeys(x) <- "GO:0016203"
EG <- mappedLkeys(x)
aveexp <- AverageExpression(object = integrate.combined_lm_0.2_forplot, features = mapIds(org.Dm.eg.db, intersect(EG, genelist), column="FLYBASE", keytype="ENTREZID", multiVals="first"), assays = "integrated", slot = "scale.data")
aveexp <- aveexp$integrated
aveexp <- aveexp[order(-aveexp$LAM3),]
go0016203<- DoHeatmap2(integrate.combined_lm_0.2_forplot, features = rownames(aveexp), angle = 0, size = 5.5, label = T, group.bar = T, draw.lines = T, group.by = "ident", cells = colnames(scaleddata_forplot),  lines.width = 100, group.bar.height = 0, hjust = 0.5, slot = "scale.data") +
  ggtitle("Muscle attachment") +
  theme(axis.text.y = element_text(face="italic")) +
  guides(color=FALSE) +
  guides(fill=FALSE) +
  #guides(fill=guide_legend(title="Log (2) expression"))+
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=14))+
  #scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),na.value = "white")
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),na.value = "white")
go0016203$data$Feature <- mapvalues(go0016203$data$Feature, from =  rownames(aveexp), to = mapIds(org.Dm.eg.db, rownames(aveexp), column="SYMBOL", keytype="FLYBASE", multiVals="first"))

x <- org.Dm.egGO2ALLEGS
Rkeys(x) <- "GO:0022610"
EG <- mappedLkeys(x)
aveexp <- AverageExpression(object = integrate.combined_lm_0.2_forplot, features = mapIds(org.Dm.eg.db, intersect(EG, genelist), column="FLYBASE", keytype="ENTREZID", multiVals="first"), assays = "integrated", slot = "scale.data")
aveexp <- aveexp$integrated
aveexp <- aveexp[order(-aveexp$LAM3),]
go0022610<- DoHeatmap2(integrate.combined_lm_0.2_forplot, features = rownames(aveexp), angle = 0, size = 5.5, label = T, group.bar = T, draw.lines = T, group.by = "ident", cells = colnames(scaleddata_forplot),  lines.width = 100, group.bar.height = 0, hjust = 0.5, slot = "scale.data") +
  ggtitle("Biological adhesion") +
  theme(axis.text.y = element_text(face="italic")) +
  guides(color=FALSE) +
  guides(fill=FALSE) +
  #guides(fill=guide_legend(title="Log (2) expression"))+
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=14))+
  #scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),na.value = "white")
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),na.value = "white")
go0022610$data$Feature <- mapvalues(go0022610$data$Feature, from =  rownames(aveexp), to = mapIds(org.Dm.eg.db, rownames(aveexp), column="SYMBOL", keytype="FLYBASE", multiVals="first"))

x <- org.Dm.egGO2ALLEGS
Rkeys(x) <- "GO:0007160"
EG <- mappedLkeys(x)
aveexp <- AverageExpression(object = integrate.combined_lm_0.2_forplot, features = mapIds(org.Dm.eg.db, intersect(EG, genelist), column="FLYBASE", keytype="ENTREZID", multiVals="first"), assays = "integrated", slot = "scale.data")
aveexp <- aveexp$integrated
aveexp <- aveexp[order(-aveexp$LAM3),]
go0007160<- DoHeatmap2(integrate.combined_lm_0.2_forplot, features = rownames(aveexp), angle = 0, size = 5.5, label = T, group.bar = T, draw.lines = T, group.by = "ident", cells = colnames(scaleddata_forplot),  lines.width = 100, group.bar.height = 0, hjust = 0.5, slot = "scale.data") +
  ggtitle("Cell-matrix adhesion") +
  theme(axis.text.y = element_text(face="italic")) +
  guides(color=FALSE) +
  guides(fill=FALSE) +
  #guides(fill=guide_legend(title="Log (2) expression"))+
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=14))+
  scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),na.value = "white")#  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),na.value = "white")
go0007160$data$Feature <- mapvalues(go0007160$data$Feature, from =  rownames(aveexp), to = mapIds(org.Dm.eg.db, rownames(aveexp), column="SYMBOL", keytype="FLYBASE", multiVals="first"))

genelist <- mapIds(org.Dm.eg.db, rownames(LAM3vPLASM1[LAM3vPLASM1$avg_logFC < (-1.4),]), column="ENTREZID", keytype="FLYBASE", multiVals="first")
go.fisher <-  goana(genelist, universe = universelist, species = "Dm", prior.prob = NULL, covariate=NULL,plot = TRUE)
go.fisher$ADJ<- p.adjust(go.fisher$P.DE, method = "fdr")
go.fisher.sig <- subset(go.fisher, go.fisher$ADJ < 0.05)
for (r in 1:length(rownames(go.fisher.sig))){
  x <- org.Dm.egGO2ALLEGS
  Rkeys(x) <- rownames(go.fisher.sig)[r]
  EG <- mappedLkeys(x)
  genes <- paste(mapIds(org.Dm.eg.db, intersect(EG, genelist), column="FLYBASE", keytype="ENTREZID", multiVals="first"),collapse=" ")
  go.fisher.sig$Genes[r] <- genes
}
write.csv(go.fisher.sig,file = "go.fisher.sig_upplasm1.csv")

x <- org.Dm.egGO2ALLEGS
Rkeys(x) <- "GO:0043062"
EG <- mappedLkeys(x)
aveexp <- AverageExpression(object = integrate.combined_lm_0.2_forplot, features = mapIds(org.Dm.eg.db, intersect(EG, genelist), column="FLYBASE", keytype="ENTREZID", multiVals="first"), assays = "integrated", slot = "scale.data")
aveexp <- aveexp$integrated
aveexp <- aveexp[order(-aveexp$LAM3),]
go0043062<- DoHeatmap2(integrate.combined_lm_0.2_forplot, features = rownames(aveexp), angle = 0, size = 5.5, label = T, group.bar = T, draw.lines = T, group.by = "ident", cells = colnames(scaleddata_forplot),  lines.width = 100, group.bar.height = 0, hjust = 0.5, slot = "scale.data") +
  ggtitle("Extracellular structure organization") +
  theme(axis.text.y = element_text(face="italic")) +
  guides(color=FALSE) +
  guides(fill=FALSE) +
  #guides(fill=guide_legend(title="Log (2) expression"))+
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=14))+
  #scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)),na.value = "white")
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),na.value = "white")
go0043062$data$Feature <- mapvalues(go0043062$data$Feature, from =  rownames(aveexp), to = mapIds(org.Dm.eg.db, rownames(aveexp), column="SYMBOL", keytype="FLYBASE", multiVals="first"))

goright <- plot_grid(go0022610,go0043062, ncol = 1, rel_heights = c(1.55,1))
goleftright <- plot_grid(go0030029,goright, ncol = 2, rel_heights = c(1,1.3))

##Figure 4
pdf("lineage_go_plot.pdf", height = 9.5, width = 11)
plot_grid(lineage_plot, goleftright, labels = c('', 'D'), label_size = 14, ncol = 1, rel_heights = c(1, 2.2))
dev.off()

setEPS()
postscript("lineage_go_plot.eps", height = 9.5, width = 11)
plot_grid(lineage_plot, goleftright, labels = c('', 'D'), label_size = 14, ncol = 1, rel_heights = c(1, 2.2))
dev.off()

####SECTION 12 - TRAJECTORY INFERENCE FOR LAMELLOCYTE DIFFERENTIATION WITH CELL CYCLE CORRECTION########

integrate.combined_lm_0.2_cellcycle <- ScaleData(integrate.combined_lm, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt", "CC.Difference"), features = rownames(integrate.combined_lm))
integrate.combined_lm_0.2_cellcycle <- RunPCA(integrate.combined_lm_0.2_cellcycle, npcs = 50, verbose = FALSE)
integrate.combined_lm_0.2_cellcycle <- RunUMAP(integrate.combined_lm_0.2_cellcycle, reduction = "pca", dims = 1:50)
integrate.combined_lm_0.2_cellcycle <- FindNeighbors(integrate.combined_lm_0.2_cellcycle, reduction = "pca", dims = 1:50)
integrate.combined_lm_0.2_cellcycle <- FindClusters(integrate.combined, resolution = 0.2)
seurat_clusters_cc <- Idents(integrate.combined_lm_0.2_cellcycle)
seurat_clusters_cc <- mapvalues(seurat_clusters_cc, from = c("0","1","2","3"), to = c("PLASM1","LAM2","LAM3","LAM1"))
seurat_clusters_cc <- factor(seurat_clusters_cc, levels = c("PLASM1","LAM1","LAM2","LAM3"))
Idents(integrate.combined_lm_0.2_cellcycle) <- seurat_clusters_cc 
integrate.combined_lm_0.2_cellcycle$seurat_clusters <- seurat_clusters_cc

lm_dimplot_cc <- DimPlot(object = integrate.combined_lm_0.2_cellcycle, reduction = "umap",label = F, order = rev(c( "PLASM1","LAM1","LAM2","LAM3")), cols = c('#9970ab',"#F28848","#FCBC2A", "#F0F921"))
lm_dimplot_data_cc <- lm_dimplot_cc$data
lm_dimplot_data_mean_cc <- aggregate(lm_dimplot_data_cc[, 1:2], list(lm_dimplot_data_cc$ident), mean)

lm_dimplot_lineage_cc <- lm_dimplot_cc +
  geom_line(data = lm_dimplot_data_mean_cc, aes(UMAP_1,UMAP_2), size = 1)+
  geom_point(data = lm_dimplot_data_mean_cc, aes(UMAP_1,UMAP_2), size = 2)+ 
  theme(legend.position="right")

top_row <- plot_grid(lamphase, lm_dimplot_lineage_cc, labels = c('a', 'b'), label_size = 24)
plot_grid(top_row, NULL,top100markersplot, labels = c('','','c'), label_size = 24, ncol = 1, rel_heights = c(1,0.3,4))


cellscompare <- ""
clustersnames <- names(table(Idents(integrate.combined_lm_0.2)))
for (i in 1:(length(table(Idents(integrate.combined_lm_0.2))))){
  tempcells <- WhichCells(integrate.combined_lm_0.2, idents = clustersnames[i])
  tempcells_data <- subset(integrate.combined_lm_0.2_cellcycle, cells = tempcells)
  celllist <- Idents(tempcells_data)
  celllist <- factor(celllist, levels = clustersnames)
  cellscompare <- rbind(cellscompare,table(celllist))
  
}
cellscompare <- cellscompare[-1,]
rownames(cellscompare) <- clustersnames
cellscompare

write.csv(cellscompare, file ="lam subcluster with and without cc compare.csv", quote = F )


####SECTION 13 - INTEGRATING ROUND 2 SUBCLUSTERING WITH ROUND 3 CLUSTERING  ########

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
identities_mod <- factor(identities_mod, levels = c("PLASM1","PLASM2","CC","MET","AMP","LAM1","LAM2","LAM3"))
identities_mod <- mapvalues(identities_mod, from = c("3","6","7","8"), to = c("PLASM2","CC","MET","AMP"))
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
                       textGrob("No Parasitism",  rot = 90, gp=gpar(col="black", fontsize=16, fontface="bold")))
Sellabel <- grobTree(rectGrob(gp=gpar(fill="NA", lwd = 0)),
                     textGrob("High Parasitism", rot = 90, gp=gpar(col="black", fontsize=16, fontface="bold")))
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
integrate.combined_cellcount$Treatment <- factor(integrate.combined_cellcount$Treatment, levels = c("No Parasitism, No Infection", "No Parasitism, Infection", "High Parasitism, No Infection", "High Parasitism, Infection"))
integrate.combined_cellcount <- separate(integrate.combined_cellcount, Treatment, into = c("Parasitism", "Infection"), sep=", ", remove = F)
integrate.combined_cellcount$Parasitism <- factor(integrate.combined_cellcount$Parasitism,levels = c("No Parasitism","High Parasitism"))
integrate.combined_cellcount$Infection <- factor(integrate.combined_cellcount$Infection, levels = c("No Infection", "Infection"))

d=data.frame(x1=c(2.3,3.2,2.3,3.2), x2=c(2.5,3.4,2.5,3.4), y1=c(2,2,3,3), y2=c(2.8,2.8,3.8,3.8), t=c('a','b','a','b'), r=c('c','c','d','d'),s=c('','','No \nParasitism','High \nParasitism'),u=c('Infection','','No Infection',''))
exptlegend <- ggplot() + 
  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t,alpha=r))+
  scale_alpha_manual(name = "", values = c(1,0.5))+
  scale_fill_manual(name = "", values = c("#0072b2","#d55e00"))+
  theme_void()+
  theme(legend.position="none")+
  geom_text(data=d, aes(x=x1+(x2-x1)/2, y=(y1+(y2-y1)/2)+1.5, label=s), size=4)+
  geom_text(data=d, aes(x=(x1+(x2-x1)/2)-0.8, y=y1+(y2-y1)/2, label=u), size=4)+
  xlim(0.8,4)+
  ylim(0.8,6)+
  ggtitle("Selection Regime")+
  theme(plot.title = element_text(hjust = 0.72, vjust = 0, face="bold"))

clusterplot <- ggplot(data=integrate.combined_cellcount, aes(x=Treatment, y=Prop, fill = Parasitism, alpha = Infection)) +
  geom_bar(stat="identity")+
  facet_wrap(~ factor(Cluster,levels = rev(c("LAM3","LAM2","LAM1","AMP","MET","CC","PLASM2","PLASM1"))), ncol=5,scales = "free", as.table = F) +
  xlab("")+
  ylab("Proportion of cells") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_alpha_manual(name = "", values = c(0.5,1),  guide = guide_legend(override.aes = list(fill = c("#0072b2", "#d55e00"), alpha = 0.5)))+ 
  scale_fill_manual(name = "", values = c("#0072b2", "#d55e00")) +
  theme(legend.position="none")

g <- ggplot_gtable(ggplot_build(clusterplot))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c('#9970ab','#8c96c6',"#DC3220",'#99d8c9','#41ae76',"#F28848","#FCBC2A", "#F0F921")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

gnew <- ggdraw(clusterplot) + draw_plot(exptlegend,x = 0.67, y = 0.55, width = 0.3, height = 0.45)

lay2 <- rbind(c("NA","NA","NA","NA","NA",1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,1,1,1,1,1,1,1,1),
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
             c(3,3,3,3,3,3,3,3,3,3,3,3,3),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3))

pdf("integrate.combined_0.3_mod_umap.pdf", width = 10, height = 10)
grid.arrange(arrangeGrob(dimplot_data2, top = grid::textGrob("B", x = 0, hjust = 0,gp=gpar(fontsize=15,fontface=2))),arrangeGrob(dimplot_legend), arrangeGrob(gnew, top = grid::textGrob("C", x = 0, hjust = 0,gp=gpar(fontsize=15,fontface=2))), top = grid::textGrob("A", x = 0, hjust = 0,vjust = 2, gp=gpar(fontsize=15,fontface=2)),layout_matrix = lay2)
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

##cluster phase
phase_table <- table(integrate.combined_0.3_mod$seurat_clusters,integrate.combined_0.3_mod$Phase)
phase_table <- phase_table/rowSums(phase_table)
phase_table <- melt(phase_table)
colnames(phase_table) <- c("Cluster","Phase","Proportion")
allcellsphaseplot <- ggplot(data=phase_table, aes(x=factor(Cluster,levels = c( "PLASM1","PLASM2","MET","AMP","CC","LAM1","LAM2","LAM3")), y=Proportion, fill = Phase)) +
  geom_bar(stat="identity")+
  xlab("Cluster")+
  ylab("Proportion")+
  ggtitle("All cells")

##cluster phase by poptreat
groups <- unique(integrate.combined_0.3_mod$poptreat)
phaseplots <- list()
for (g in 1:length(groups)){
  integrate.combined_0.3_mod_sub <- subset(integrate.combined_0.3_mod, poptreat == groups[g])
  phase_table <- table(integrate.combined_0.3_mod_sub$seurat_clusters,integrate.combined_0.3_mod_sub$Phase)
  phase_table <- phase_table/rowSums(phase_table)
  phase_table <- melt(phase_table)
  colnames(phase_table) <- c("Cluster","Phase","Proportion")
  phaseplots[[g]] <- ggplot(data=phase_table, aes(x=factor(Cluster,levels = c( "PLASM1","PLASM2","CC","MET","AMP","LAM1","LAM2","LAM3")), y=Proportion, fill = Phase)) +
    geom_bar(stat="identity")+
    xlab("Cluster")+
    ylab("Proportion")+
    ggtitle(groups[g])
  
}
phaseplotbypoptreat <- do.call("grid.arrange", c(phaseplots, ncol=2))
plot_grid(allcellsphaseplot,phaseplotbypoptreat, ncol=1)

##compare to integrate.combined cellcycle
cellscompare <- ""
clustersnames <- names(table(Idents(integrate.combined_0.3_mod)))
for (i in 1:(length(table(Idents(integrate.combined_0.3_mod))))){
  tempcells <- WhichCells(integrate.combined_0.3_mod, idents = clustersnames[i])
  tempcells_data <- subset(integrate.combined_cellcycle_0.3, cells = tempcells)
  celllist <- Idents(tempcells_data)
  celllist <- factor(celllist, levels = 0:(length(table(Idents(integrate.combined_cellcycle_0.3)))-1))
  cellscompare <- rbind(cellscompare,table(celllist))
  
}
cellscompare <- cellscompare[-1,]
rownames(cellscompare) <- clustersnames
write.csv(cellscompare, file = "row_all_mod_col_all_cc_0.3.csv")

##compare to integrate.combined 15000 features
cellscompare <- ""
clustersnames <- names(table(Idents(integrate.combined_0.3_mod)))
for (i in 1:(length(table(Idents(integrate.combined_0.3_mod))))){
  tempcells <- WhichCells(integrate.combined_0.3_mod, idents = clustersnames[i])
  tempcells_data <- subset(integrate.combined_15000f_0.3, cells = tempcells)
  celllist <- Idents(tempcells_data)
  celllist <- factor(celllist, levels = 0:(length(table(Idents(integrate.combined_15000f_0.3)))-1))
  cellscompare <- rbind(cellscompare,table(celllist))
  
}
cellscompare <- cellscompare[-1,]
rownames(cellscompare) <- clustersnames
write.csv(cellscompare, file = "row_all_mod_col_7716features_0.3.csv", quote = F)

####SECTION 14 - CLUSTER DEFINITIONS   ########

##for flymine pathway analysis
universe_sc <- read.csv(file = "universe_sc")
universe_sc_2000 <- rownames(integrate.combined_0.3_mod)
universe_sc_lm_2000 <- rownames(integrate.combined_lm_0.2)
write.table(universe_sc_2000, file = "universe_sc_2000.txt", quote = F, row.names = F, col.names = F)
write.table(universe_sc_lm_2000, file = "universe_sc_2000.txt", quote = F, row.names = F, col.names = F)
write.table(universe_sc_lm_2000, file = "universe_sc_lm_2000.txt", quote = F, row.names = F, col.names = F)

##cell cycle and plasmatocyte markers
ccplasmmarkers <- DotPlot(integrate.combined_0.3_mod, features = c("FBgn0003124","FBgn0003525","FBgn0261385", "FBgn0259896", "FBgn0029167", "FBgn0243514", "FBgn0011828", "FBgn0000299")) +
  scale_x_discrete(labels=c("Col4a1","Pxn","eater","Hml","NimC1", "scra","stg","polo")) + 
  theme(axis.text.x = element_text(face = "italic"))+
  ylab("Cluster")+
  xlab("")

##lamellocyte markers
##lamellocyte markers
lam_markers <- c("FBgn0261363","FBgn0032422","FBgn0004657","FBgn0034005")
lam_markers_name <- c("PPO3", "atilla", "mys", "ItgaPS4")
lam_markers_plot <- list()

for (l in 1:length(lam_markers)){
  lam_markers_plot[[l]] <- VlnPlot(object = integrate.combined_0.3_mod, features = lam_markers[l], pt.size = 0) +
    xlab("") +
    ggtitle(lam_markers_name[l])+
    NoLegend()+
    theme(plot.title = element_text(size=14, face="bold.italic"))
}
lammarkersall <- grid.arrange(lam_markers_plot[[1]],lam_markers_plot[[2]],lam_markers_plot[[3]],lam_markers_plot[[4]],ncol=2)

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
alltopmarkers <- DotPlot(integrate.combined_0.3_mod, features = cluster_markers_top$gene)+
  scale_x_discrete(labels=rev(cluster_markers_top$symbol)) + 
  theme(axis.text.x = element_text(face = "italic"))+
  ylab("Cluster")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(integrate.combined_0.3_mod, features = unique(g2m.genes))

plot_grid(ccplasmmarkers,alltopmarkers,lammarkersall,ncol=1, rel_heights = c(1,1,2.5), labels = "auto")

####SECTION 15 - SAVE SEURAT OBJECTS AND COUNT MATRICES   ########

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


####SECTION 16 - SAMPLE LIBRARY METRICS  ########

integrate.combined_0.3_mod_formetrics <- integrate.combined_0.3_mod
integrate.combined_0.3_mod_formetrics$old.ident <- gsub("NSRef","HP", integrate.combined_0.3_mod_formetrics$old.ident)
integrate.combined_0.3_mod_formetrics$old.ident <- gsub("^C","NP-", integrate.combined_0.3_mod_formetrics$old.ident)
integrate.combined_0.3_mod_formetrics$old.ident <- gsub("Inf","I", integrate.combined_0.3_mod_formetrics$old.ident)
integrate.combined_0.3_mod_formetrics$old.ident <- gsub("Uninf","NI", integrate.combined_0.3_mod_formetrics$old.ident)

Idents(integrate.combined_0.3_mod_formetrics) <- integrate.combined_0.3_mod_formetrics$old.ident
metricplot <- VlnPlot(integrate.combined_0.3_mod_formetrics, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2, pt.size = 0)
write.csv(table(integrate.combined_0.3_mod_formetrics$old.ident), file = "cell number by library.csv", quote = F)


numcellslib <- as.data.frame(table(integrate.combined_0.3_mod_formetrics$old.ident))
colnames(numcellslib) <- c("Library", "number.of.cells")
numcellslib <- rbind(numcellslib[5:8,],numcellslib[1:4,])
numcellplot <- ggplot(numcellslib, aes(factor(Library, levels = as.character(unique(numcellslib$Library))),number.of.cells, fill = factor(Library, levels = as.character(unique(numcellslib$Library)))))+
  geom_bar(stat="identity") +
  ggtitle("Number of cells")+ 
  ylab("")+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Identity")

plot_grid(numcellplot,metricplot, nrow=1, rel_widths = c(1,3))

numgenslib <- c(9930, 11067, 10743, 10511, 9380, 9017, 9729, 10384)
names(numgenslib) <- numcellslib$Var1

####SECTION 17 - CLUSTER COMPARISONS BY REPLICATE AND FOLLOWING SUBSAMPLING CELLS ########

##rep 3 only

integrate.list <- list(C3_Inf_obj, C3_Uninf_obj,NSRef_3_Inf_obj, NSRef_3_Uninf_obj)
names(integrate.list)=c("C3_Inf_obj", "C3_Uninf_obj","NSRef_3_Inf_obj", "NSRef_3_Uninf_obj")
for (i in 1:length(x = integrate.list)) {
  integrate.list[[i]] <- NormalizeData(object = integrate.list[[i]], verbose = FALSE)
  integrate.list[[i]] <- FindVariableFeatures(object = integrate.list[[i]],  selection.method = "vst", nfeatures = 2000, verbose = TRUE)
}
reference.list <- integrate.list[c("C3_Inf_obj", "C3_Uninf_obj","NSRef_3_Inf_obj", "NSRef_3_Uninf_obj")]
integrate.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50, anchor.features = 2000)
integrate.combined_rep3 <- IntegrateData(anchorset = integrate.anchors, dims = 1:50)
integrate.combined_rep3 <- subset(integrate.combined_rep3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
DefaultAssay(integrate.combined_rep3) <- "integrated"
s.genes <- read.csv(file = "s.genes.fly.csv")
g2m.genes <- read.csv(file = "g2m.genes.fly.csv")
s.genes <- as.character(s.genes$FlyBaseID)
g2m.genes <- as.character(g2m.genes$FlyBaseID)
s.genes <- s.genes[s.genes %in% rownames(integrate.combined_rep3)]
g2m.genes <- g2m.genes[g2m.genes %in% rownames(integrate.combined_rep3)]
integrate.combined_rep3 <- CellCycleScoring(integrate.combined_rep3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
integrate.combined_rep3$CC.Difference <- integrate.combined_rep3$S.Score - integrate.combined_rep3$G2M.Score
integrate.combined_rep3$treatment <- integrate.combined_rep3$orig.ident
integrate.combined_rep3$treatment <- gsub("NSRef-._","",integrate.combined_rep3$treatment)
integrate.combined_rep3$treatment <- gsub("C._","",integrate.combined_rep3$treatment)
integrate.combined_rep3$treatment <- gsub("Inf","Infected",integrate.combined_rep3$treatment)
integrate.combined_rep3$population <- integrate.combined_rep3$orig.ident
integrate.combined_rep3$population <- gsub("_Uninf","",integrate.combined_rep3$population)
integrate.combined_rep3$population <- gsub("_Inf","",integrate.combined_rep3$population)
integrate.combined_rep3$population <- gsub("-","",integrate.combined_rep3$population)
integrate.combined_rep3$population <- gsub(".$","",integrate.combined_rep3$population)
integrate.combined_rep3$population <- gsub("C","control",integrate.combined_rep3$population)
integrate.combined_rep3$poptreat <- paste(integrate.combined_rep3$population,integrate.combined_rep3$treatment)
integrate.combined_rep3_cellcycle <- integrate.combined_rep3
integrate.combined_rep3 <- ScaleData(integrate.combined_rep3, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"), features = rownames(integrate.combined_rep3))
integrate.combined_rep3 <- RunPCA(integrate.combined_rep3, npcs = 50, verbose = FALSE)
integrate.combined_rep3 <- RunUMAP(integrate.combined_rep3, reduction = "pca", dims = 1:50)
integrate.combined_rep3 <- FindNeighbors(integrate.combined_rep3, reduction = "pca", dims = 1:50)

integrate.combined_rep3_0.3 <- FindClusters(integrate.combined_rep3, resolution = 0.5)
integrate.combined_rep3_0.3$poptreat <- gsub("control Infected", "No Parasitism, Infection", integrate.combined_rep3_0.3$poptreat)
integrate.combined_rep3_0.3$poptreat <- gsub("control Uninf", "No Parasitism, No Infection", integrate.combined_rep3_0.3$poptreat)
integrate.combined_rep3_0.3$poptreat <- gsub("NSRef Infected", "High Parasitism, Infection", integrate.combined_rep3_0.3$poptreat)
integrate.combined_rep3_0.3$poptreat <- gsub("NSRef Uninf", "High Parasitism, No Infection", integrate.combined_rep3_0.3$poptreat)
integrate.combined_rep3_0.3$poptreat <- factor(integrate.combined_rep3_0.3$poptreat, levels = c("No Parasitism, No Infection", "No Parasitism, Infection", "High Parasitism, No Infection", "High Parasitism, Infection"))
integrate.combined_rep3_0.3$population <- gsub("control", "No Parasitism", integrate.combined_rep3_0.3$population)
integrate.combined_rep3_0.3$population <- gsub("NSRef", "High Parasitism", integrate.combined_rep3_0.3$population)
integrate.combined_rep3_0.3$population <- factor(integrate.combined_rep3_0.3$population, levels = c("No Parasitism","High Parasitism"))
integrate.combined_rep3_0.3$treatment <- gsub("Infected", "Infection", integrate.combined_rep3_0.3$treatment)
integrate.combined_rep3_0.3$treatment <- gsub("Uninf", "No Infection", integrate.combined_rep3_0.3$treatment)
integrate.combined_rep3_0.3$treatment <- factor(integrate.combined_rep3_0.3$treatment, levels = c("No Infection","Infection"))

#DimPlot(integrate.combined_rep3_0.3, split.by = "poptreat", label = T)
seurat_clusters_rep3 <- as.data.frame(integrate.combined_rep3_0.3$seurat_clusters)
seurat_clusters_mod <- as.data.frame(integrate.combined_0.3_mod$seurat_clusters)
seurat_clusters_rep3$cell <- rownames(seurat_clusters_rep3)
seurat_clusters_mod$cell <- rownames(seurat_clusters_mod)
newclusters <- left_join(seurat_clusters_rep3, seurat_clusters_mod, by = "cell")
Idents(integrate.combined_rep3_0.3) <- newclusters$`integrate.combined_0.3_mod$seurat_clusters`
integrate.combined_rep3_0.3$seurat_clusters <- Idents(integrate.combined_rep3_0.3)

cellscompare <- ""
clustersnames <- names(table(Idents(integrate.combined_0.3_mod)))
for (i in 1:(length(table(Idents(integrate.combined_0.3_mod))))){
  tempcells <- WhichCells(integrate.combined_0.3_mod, idents = clustersnames[i])
  tempcells_data <- subset(integrate.combined_rep3_0.3, cells = tempcells)
  celllist <- Idents(tempcells_data)
  celllist <- factor(celllist, levels = 0:(length(table(Idents(integrate.combined_rep3_0.3)))-1))
  cellscompare <- rbind(cellscompare,table(celllist))
  
}
cellscompare <- cellscompare[-1,]
rownames(cellscompare) <- clustersnames
cellscompare

integrate.combined_cellcount <- as.data.frame(table(paste(integrate.combined_rep3_0.3$poptreat,integrate.combined_rep3_0.3$seurat_clusters,sep="_")))
integrate.combined_cellcount <- integrate.combined_cellcount %>% separate(Var1, into= c("Treatment","Cluster"),sep="_")
propall <- vector()
for (t in levels(as.factor(integrate.combined_cellcount$Treatment))){
  hold <- integrate.combined_cellcount[integrate.combined_cellcount$Treatment == t,]
  prop <- hold$Freq/sum(hold$Freq)
  propall <- c(propall,prop)
}
integrate.combined_cellcount$Prop <- propall
integrate.combined_cellcount$Treatment <- factor(integrate.combined_cellcount$Treatment, levels = c("No Parasitism, No Infection", "No Parasitism, Infection", "High Parasitism, No Infection", "High Parasitism, Infection"))
integrate.combined_cellcount <- separate(integrate.combined_cellcount, Treatment, into = c("Parasitism", "Infection"), sep=", ", remove = F)
integrate.combined_cellcount$Parasitism <- as.factor(integrate.combined_cellcount$Parasitism)
integrate.combined_cellcount$Infection <- factor(integrate.combined_cellcount$Infection, levels = c("No Infection", "Infection"))

rep3dimplot <- DimPlot(integrate.combined_rep3_0.3,order = rev(c( "PLASM1","PLASM2","MET","AMP","CC","LAM1","LAM2","LAM3")), cols = c('#9970ab','#8c96c6','#99d8c9','#41ae76',"#DC3220","#F28848","#FCBC2A", "#F0F921"), label = T)+
  ggtitle("B) Replicate 3")+
  theme(plot.title = element_text(hjust = 0))

clusterplot_rep3 <- ggplot(data=integrate.combined_cellcount, aes(x=Treatment, y=Prop, fill = Parasitism, alpha = Infection)) +
  geom_bar(stat="identity")+
  facet_wrap(~ factor(Cluster,levels = rev(c("LAM3","LAM2","LAM1","AMP","MET","CC","PLASM2","PLASM1"))), ncol=5,scales = "free", as.table = F) +
  xlab("")+
  ylab("Proportion of cells") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_alpha_manual(name = "", values = c(0.5,1),  guide = guide_legend(override.aes = list(fill = c("#0072b2", "#d55e00"), alpha = 0.5)))+ 
  scale_fill_manual(name = "", values = c("#0072b2", "#d55e00")) +
  theme(legend.position="none")

gnew_rep3 <- ggdraw(clusterplot_rep3) + draw_plot(exptlegend,x = 0.67, y = 0.6, width = 0.3, height = 0.35)

rep3dimplot_clus <- plot_grid(rep3dimplot,gnew_rep3, ncol=2, rel_widths = c(1,2))

##rep1

integrate.list <- list(C1_Inf_obj, C1_Uninf_obj,NSRef_1_Inf_obj, NSRef_1_Uninf_obj)
names(integrate.list)=c("C1_Inf_obj", "C1_Uninf_obj","NSRef_1_Inf_obj", "NSRef_1_Uninf_obj")
for (i in 1:length(x = integrate.list)) {
  integrate.list[[i]] <- NormalizeData(object = integrate.list[[i]], verbose = FALSE)
  integrate.list[[i]] <- FindVariableFeatures(object = integrate.list[[i]],  selection.method = "vst", nfeatures = 2000, verbose = TRUE)
}
reference.list <- integrate.list[c("C1_Inf_obj", "C1_Uninf_obj","NSRef_1_Inf_obj", "NSRef_1_Uninf_obj")]
integrate.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50, anchor.features = 2000)
integrate.combined_rep1 <- IntegrateData(anchorset = integrate.anchors, dims = 1:50)
integrate.combined_rep1 <- subset(integrate.combined_rep1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
DefaultAssay(integrate.combined_rep1) <- "integrated"
s.genes <- read.csv(file = "s.genes.fly.csv")
g2m.genes <- read.csv(file = "g2m.genes.fly.csv")
s.genes <- as.character(s.genes$FlyBaseID)
g2m.genes <- as.character(g2m.genes$FlyBaseID)
s.genes <- s.genes[s.genes %in% rownames(integrate.combined_rep1)]
g2m.genes <- g2m.genes[g2m.genes %in% rownames(integrate.combined_rep1)]
integrate.combined_rep1 <- CellCycleScoring(integrate.combined_rep1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
integrate.combined_rep1$CC.Difference <- integrate.combined_rep1$S.Score - integrate.combined_rep1$G2M.Score
integrate.combined_rep1$treatment <- integrate.combined_rep1$orig.ident
integrate.combined_rep1$treatment <- gsub("NSRef-._","",integrate.combined_rep1$treatment)
integrate.combined_rep1$treatment <- gsub("C._","",integrate.combined_rep1$treatment)
integrate.combined_rep1$treatment <- gsub("Inf","Infected",integrate.combined_rep1$treatment)
integrate.combined_rep1$population <- integrate.combined_rep1$orig.ident
integrate.combined_rep1$population <- gsub("_Uninf","",integrate.combined_rep1$population)
integrate.combined_rep1$population <- gsub("_Inf","",integrate.combined_rep1$population)
integrate.combined_rep1$population <- gsub("-","",integrate.combined_rep1$population)
integrate.combined_rep1$population <- gsub(".$","",integrate.combined_rep1$population)
integrate.combined_rep1$population <- gsub("C","control",integrate.combined_rep1$population)
integrate.combined_rep1$poptreat <- paste(integrate.combined_rep1$population,integrate.combined_rep1$treatment)
integrate.combined_rep1_cellcycle <- integrate.combined_rep1
integrate.combined_rep1 <- ScaleData(integrate.combined_rep1, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"), features = rownames(integrate.combined_rep1))
integrate.combined_rep1 <- RunPCA(integrate.combined_rep1, npcs = 50, verbose = FALSE)
integrate.combined_rep1 <- RunUMAP(integrate.combined_rep1, reduction = "pca", dims = 1:50)
integrate.combined_rep1 <- FindNeighbors(integrate.combined_rep1, reduction = "pca", dims = 1:50)

integrate.combined_rep1_0.3 <- FindClusters(integrate.combined_rep1, resolution = 0.5)
integrate.combined_rep1_0.3$poptreat <- gsub("control Infected", "No Parasitism, Infection", integrate.combined_rep1_0.3$poptreat)
integrate.combined_rep1_0.3$poptreat <- gsub("control Uninf", "No Parasitism, No Infection", integrate.combined_rep1_0.3$poptreat)
integrate.combined_rep1_0.3$poptreat <- gsub("NSRef Infected", "High Parasitism, Infection", integrate.combined_rep1_0.3$poptreat)
integrate.combined_rep1_0.3$poptreat <- gsub("NSRef Uninf", "High Parasitism, No Infection", integrate.combined_rep1_0.3$poptreat)
integrate.combined_rep1_0.3$poptreat <- factor(integrate.combined_rep1_0.3$poptreat, levels = c("No Parasitism, No Infection", "No Parasitism, Infection", "High Parasitism, No Infection", "High Parasitism, Infection"))
integrate.combined_rep1_0.3$population <- gsub("control", "No Parasitism", integrate.combined_rep1_0.3$population)
integrate.combined_rep1_0.3$population <- gsub("NSRef", "High Parasitism", integrate.combined_rep1_0.3$population)
integrate.combined_rep1_0.3$population <- factor(integrate.combined_rep1_0.3$population, levels = c("No Parasitism","High Parasitism"))
integrate.combined_rep1_0.3$treatment <- gsub("Infected", "Infection", integrate.combined_rep1_0.3$treatment)
integrate.combined_rep1_0.3$treatment <- gsub("Uninf", "No Infection", integrate.combined_rep1_0.3$treatment)
integrate.combined_rep1_0.3$treatment <- factor(integrate.combined_rep1_0.3$treatment, levels = c("No Infection","Infection"))

seurat_clusters_rep1 <- as.data.frame(integrate.combined_rep1_0.3$seurat_clusters)
seurat_clusters_mod <- as.data.frame(integrate.combined_0.3_mod$seurat_clusters)
seurat_clusters_rep1$cell <- rownames(seurat_clusters_rep1)
seurat_clusters_mod$cell <- rownames(seurat_clusters_mod)
newclusters <- left_join(seurat_clusters_rep1, seurat_clusters_mod, by = "cell")
Idents(integrate.combined_rep1_0.3) <- newclusters$`integrate.combined_0.3_mod$seurat_clusters`
integrate.combined_rep1_0.3$seurat_clusters <- Idents(integrate.combined_rep1_0.3)

rep1dimplot <- DimPlot(integrate.combined_rep1_0.3, order = rev(c( "PLASM1","PLASM2","MET","AMP","CC","LAM1","LAM2","LAM3")), cols = c('#9970ab','#8c96c6','#99d8c9','#41ae76',"#DC3220","#F28848","#FCBC2A", "#F0F921"), label = T)+
  ggtitle("A) Replicate 1")+
  theme(plot.title = element_text(hjust = 0))

integrate.combined_cellcount <- as.data.frame(table(paste(integrate.combined_rep1_0.3$poptreat,integrate.combined_rep1_0.3$seurat_clusters,sep="_")))
integrate.combined_cellcount <- integrate.combined_cellcount %>% separate(Var1, into= c("Treatment","Cluster"),sep="_")
propall <- vector()
for (t in levels(as.factor(integrate.combined_cellcount$Treatment))){
  hold <- integrate.combined_cellcount[integrate.combined_cellcount$Treatment == t,]
  prop <- hold$Freq/sum(hold$Freq)
  propall <- c(propall,prop)
}
integrate.combined_cellcount$Prop <- propall
integrate.combined_cellcount$Treatment <- factor(integrate.combined_cellcount$Treatment, levels = c("No Parasitism, No Infection", "No Parasitism, Infection", "High Parasitism, No Infection", "High Parasitism, Infection"))
integrate.combined_cellcount <- separate(integrate.combined_cellcount, Treatment, into = c("Parasitism", "Infection"), sep=", ", remove = F)
integrate.combined_cellcount$Parasitism <- as.factor(integrate.combined_cellcount$Parasitism)
integrate.combined_cellcount$Infection <- factor(integrate.combined_cellcount$Infection, levels = c("No Infection", "Infection"))

clusterplot_rep1 <- ggplot(data=integrate.combined_cellcount, aes(x=Treatment, y=Prop, fill = Parasitism, alpha = Infection)) +
  geom_bar(stat="identity")+
  facet_wrap(~ factor(Cluster,levels = rev(c("LAM3","LAM2","LAM1","AMP","MET","CC","PLASM2","PLASM1"))), ncol=5,scales = "free", as.table = F) +
  xlab("")+
  ylab("Proportion of cells") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_alpha_manual(name = "", values = c(0.5,1),  guide = guide_legend(override.aes = list(fill = c("#0072b2", "#d55e00"), alpha = 0.5)))+ 
  scale_fill_manual(name = "", values = c("#0072b2", "#d55e00")) +
  theme(legend.position="none")

gnew_rep1 <- ggdraw(clusterplot_rep1) + draw_plot(exptlegend,x = 0.67, y = 0.6, width = 0.3, height = 0.35)

rep1dimplot_clus <- plot_grid(rep1dimplot,gnew_rep1, ncol=2, rel_widths = c(1,2))

##subsample cell to min library

cl1 <- sample(names(integrate.combined_0.3_mod$old.ident[integrate.combined_0.3_mod$old.ident %in% levels(integrate.combined_0.3_mod$old.ident)[1]]),min(table(integrate.combined_0.3_mod$old.ident)), replace = F)
cl2 <- sample(names(integrate.combined_0.3_mod$old.ident[integrate.combined_0.3_mod$old.ident %in% levels(integrate.combined_0.3_mod$old.ident)[2]]),min(table(integrate.combined_0.3_mod$old.ident)), replace = F)
cl3 <- sample(names(integrate.combined_0.3_mod$old.ident[integrate.combined_0.3_mod$old.ident %in% levels(integrate.combined_0.3_mod$old.ident)[3]]),min(table(integrate.combined_0.3_mod$old.ident)), replace = F)
cl4 <- sample(names(integrate.combined_0.3_mod$old.ident[integrate.combined_0.3_mod$old.ident %in% levels(integrate.combined_0.3_mod$old.ident)[4]]),min(table(integrate.combined_0.3_mod$old.ident)), replace = F)
cl5 <- sample(names(integrate.combined_0.3_mod$old.ident[integrate.combined_0.3_mod$old.ident %in% levels(integrate.combined_0.3_mod$old.ident)[5]]),min(table(integrate.combined_0.3_mod$old.ident)), replace = F)
cl6 <- sample(names(integrate.combined_0.3_mod$old.ident[integrate.combined_0.3_mod$old.ident %in% levels(integrate.combined_0.3_mod$old.ident)[6]]),min(table(integrate.combined_0.3_mod$old.ident)), replace = F)
cl7 <- sample(names(integrate.combined_0.3_mod$old.ident[integrate.combined_0.3_mod$old.ident %in% levels(integrate.combined_0.3_mod$old.ident)[7]]),min(table(integrate.combined_0.3_mod$old.ident)), replace = F)
cl8 <- sample(names(integrate.combined_0.3_mod$old.ident[integrate.combined_0.3_mod$old.ident %in% levels(integrate.combined_0.3_mod$old.ident)[8]]),min(table(integrate.combined_0.3_mod$old.ident)), replace = F)

cellall <- c(cl1,cl2,cl3,cl4,cl5,cl6,cl7,cl8) 

C1_Inf_obj_sub <- subset(C1_Inf_obj,cells = cellall)
C1_Uninf_obj_sub <- subset(C1_Uninf_obj,cells = cellall)
C3_Inf_obj_sub <- subset(C3_Inf_obj,cells = cellall)
C3_Uninf_obj_sub <- subset(C3_Uninf_obj,cells = cellall)
NSRef_1_Inf_obj_sub <- subset(NSRef_1_Inf_obj,cells = cellall)
NSRef_1_Uninf_obj_sub <- subset(NSRef_1_Uninf_obj,cells = cellall)
NSRef_3_Inf_obj_sub <- subset(NSRef_3_Inf_obj,cells = cellall)
NSRef_3_Uninf_obj_sub <- subset(NSRef_3_Uninf_obj,cells = cellall)

integrate.list <- list(C1_Inf_obj_sub, C1_Uninf_obj_sub,C3_Inf_obj_sub, C3_Uninf_obj_sub,NSRef_1_Inf_obj_sub, NSRef_1_Uninf_obj_sub,NSRef_3_Inf_obj_sub, NSRef_3_Uninf_obj_sub)
names(integrate.list)=c("C1_Inf_obj_sub", "C1_Uninf_obj_sub","C3_Inf_obj_sub", "C3_Uninf_obj_sub","NSRef_1_Inf_obj_sub", "NSRef_1_Uninf_obj_sub","NSRef_3_Inf_obj_sub", "NSRef_3_Uninf_obj_sub")
for (i in 1:length(x = integrate.list)) {
  integrate.list[[i]] <- NormalizeData(object = integrate.list[[i]], verbose = FALSE)
  integrate.list[[i]] <- FindVariableFeatures(object = integrate.list[[i]],  selection.method = "vst", nfeatures = 2000, verbose = TRUE)
}
reference.list <- integrate.list[c("C1_Inf_obj_sub", "C1_Uninf_obj_sub","C3_Inf_obj_sub", "C3_Uninf_obj_sub","NSRef_1_Inf_obj_sub", "NSRef_1_Uninf_obj_sub","NSRef_3_Inf_obj_sub", "NSRef_3_Uninf_obj_sub")]
integrate.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50, anchor.features = 2000)
integrate.combined_sub <- IntegrateData(anchorset = integrate.anchors, dims = 1:50)
integrate.combined_sub <- subset(integrate.combined_sub, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
DefaultAssay(integrate.combined_sub) <- "integrated"
s.genes <- read.csv(file = "s.genes.fly.csv")
g2m.genes <- read.csv(file = "g2m.genes.fly.csv")
s.genes <- as.character(s.genes$FlyBaseID)
g2m.genes <- as.character(g2m.genes$FlyBaseID)
s.genes <- s.genes[s.genes %in% rownames(integrate.combined_sub)]
g2m.genes <- g2m.genes[g2m.genes %in% rownames(integrate.combined_sub)]
integrate.combined_sub <- CellCycleScoring(integrate.combined_sub, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
integrate.combined_sub$CC.Difference <- integrate.combined_sub$S.Score - integrate.combined_sub$G2M.Score
integrate.combined_sub$treatment <- integrate.combined_sub$orig.ident
integrate.combined_sub$treatment <- gsub("NSRef-._","",integrate.combined_sub$treatment)
integrate.combined_sub$treatment <- gsub("C._","",integrate.combined_sub$treatment)
integrate.combined_sub$treatment <- gsub("Inf","Infected",integrate.combined_sub$treatment)
integrate.combined_sub$population <- integrate.combined_sub$orig.ident
integrate.combined_sub$population <- gsub("_Uninf","",integrate.combined_sub$population)
integrate.combined_sub$population <- gsub("_Inf","",integrate.combined_sub$population)
integrate.combined_sub$population <- gsub("-","",integrate.combined_sub$population)
integrate.combined_sub$population <- gsub(".$","",integrate.combined_sub$population)
integrate.combined_sub$population <- gsub("C","control",integrate.combined_sub$population)
integrate.combined_sub$poptreat <- paste(integrate.combined_sub$population,integrate.combined_sub$treatment)
integrate.combined_sub_cellcycle <- integrate.combined_sub
integrate.combined_sub <- ScaleData(integrate.combined_sub, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"), features = rownames(integrate.combined_sub))
integrate.combined_sub <- RunPCA(integrate.combined_sub, npcs = 50, verbose = FALSE)
integrate.combined_sub <- RunUMAP(integrate.combined_sub, reduction = "pca", dims = 1:50)
integrate.combined_sub <- FindNeighbors(integrate.combined_sub, reduction = "pca", dims = 1:50)

integrate.combined_sub_0.3 <- FindClusters(integrate.combined_sub, resolution = 0.5)
integrate.combined_sub_0.3$poptreat <- gsub("control Infected", "No Parasitism, Infection", integrate.combined_sub_0.3$poptreat)
integrate.combined_sub_0.3$poptreat <- gsub("control Uninf", "No Parasitism, No Infection", integrate.combined_sub_0.3$poptreat)
integrate.combined_sub_0.3$poptreat <- gsub("NSRef Infected", "High Parasitism, Infection", integrate.combined_sub_0.3$poptreat)
integrate.combined_sub_0.3$poptreat <- gsub("NSRef Uninf", "High Parasitism, No Infection", integrate.combined_sub_0.3$poptreat)
integrate.combined_sub_0.3$poptreat <- factor(integrate.combined_sub_0.3$poptreat, levels = c("No Parasitism, No Infection", "No Parasitism, Infection", "High Parasitism, No Infection", "High Parasitism, Infection"))
integrate.combined_sub_0.3$population <- gsub("control", "No Parasitism", integrate.combined_sub_0.3$population)
integrate.combined_sub_0.3$population <- gsub("NSRef", "High Parasitism", integrate.combined_sub_0.3$population)
integrate.combined_sub_0.3$population <- factor(integrate.combined_sub_0.3$population, levels = c("No Parasitism","High Parasitism"))
integrate.combined_sub_0.3$treatment <- gsub("Infected", "Infection", integrate.combined_sub_0.3$treatment)
integrate.combined_sub_0.3$treatment <- gsub("Uninf", "No Infection", integrate.combined_sub_0.3$treatment)
integrate.combined_sub_0.3$treatment <- factor(integrate.combined_sub_0.3$treatment, levels = c("No Infection","Infection"))

seurat_clusters_sub <- as.data.frame(integrate.combined_sub_0.3$seurat_clusters)
seurat_clusters_mod <- as.data.frame(integrate.combined_0.3_mod$seurat_clusters)
seurat_clusters_sub$cell <- rownames(seurat_clusters_sub)
seurat_clusters_mod$cell <- rownames(seurat_clusters_mod)
newclusters <- left_join(seurat_clusters_sub, seurat_clusters_mod, by = "cell")
Idents(integrate.combined_sub_0.3) <- newclusters$`integrate.combined_0.3_mod$seurat_clusters`
integrate.combined_sub_0.3$seurat_clusters <- Idents(integrate.combined_sub_0.3)

subdimplot <- DimPlot(integrate.combined_sub_0.3, order = rev(c( "PLASM1","PLASM2","MET","AMP","CC","LAM1","LAM2","LAM3")), cols = c('#9970ab','#8c96c6','#99d8c9','#41ae76',"#DC3220","#F28848","#FCBC2A", "#F0F921"), label = T)+
  ggtitle("C) Subsample cells")+
  theme(plot.title = element_text(hjust = 0))

integrate.combined_cellcount <- as.data.frame(table(paste(integrate.combined_sub_0.3$poptreat,integrate.combined_sub_0.3$seurat_clusters,sep="_")))
integrate.combined_cellcount <- integrate.combined_cellcount %>% separate(Var1, into= c("Treatment","Cluster"),sep="_")
propall <- vector()
for (t in levels(as.factor(integrate.combined_cellcount$Treatment))){
  hold <- integrate.combined_cellcount[integrate.combined_cellcount$Treatment == t,]
  prop <- hold$Freq/sum(hold$Freq)
  propall <- c(propall,prop)
}
integrate.combined_cellcount$Prop <- propall
integrate.combined_cellcount$Treatment <- factor(integrate.combined_cellcount$Treatment, levels = c("No Parasitism, No Infection", "No Parasitism, Infection", "High Parasitism, No Infection", "High Parasitism, Infection"))
integrate.combined_cellcount <- separate(integrate.combined_cellcount, Treatment, into = c("Parasitism", "Infection"), sep=", ", remove = F)
integrate.combined_cellcount$Parasitism <- as.factor(integrate.combined_cellcount$Parasitism)
integrate.combined_cellcount$Infection <- factor(integrate.combined_cellcount$Infection, levels = c("No Infection", "Infection"))

clusterplot_sub <- ggplot(data=integrate.combined_cellcount, aes(x=Treatment, y=Prop, fill = Parasitism, alpha = Infection)) +
  geom_bar(stat="identity")+
  facet_wrap(~ factor(Cluster,levels = rev(c("LAM3","LAM2","LAM1","AMP","MET","CC","PLASM2","PLASM1"))), ncol=5,scales = "free", as.table = F) +
  xlab("")+
  ylab("Proportion of cells") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_alpha_manual(name = "", values = c(0.5,1),  guide = guide_legend(override.aes = list(fill = c("#0072b2", "#d55e00"), alpha = 0.5)))+ 
  scale_fill_manual(name = "", values = c("#0072b2", "#d55e00")) +
  theme(legend.position="none")

gnew_sub <- ggdraw(clusterplot_sub) + draw_plot(exptlegend,x = 0.67, y = 0.6, width = 0.3, height = 0.35)

subdimplot_clus <- plot_grid(subdimplot,gnew_sub, ncol=2, rel_widths = c(1,2))

pdf("subsample.pdf", height = 13, width = 15)
plot_grid(rep1dimplot_clus,rep3dimplot_clus,subdimplot_clus, ncol = 1)
dev.off()
