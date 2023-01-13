# conda activate ST
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tidyverse)

sc_ovarian_filtered <- 
  readRDS("outputs/seurat_sc_ovarian_filtered.rds")

sp_ovarian_filtered <- 
  readRDS("outputs/filtered_seurat_sp_ovarian.rds")

# --- Cell Cycle Effect/Regression (scRNA) ---
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sc_ovarian_filtered <- CellCycleScoring(sc_ovarian_filtered, s.features = s.genes, 
                                        g2m.features = g2m.genes, set.ident = TRUE)
# No regression

sc_ovarian_noregress <- SCTransform(sc_ovarian_filtered)
sc_ovarian_noregress <- RunPCA(sc_ovarian_noregress, assay = "SCT", verbose = FALSE)
sc_ovarian_noregress <- FindNeighbors(sc_ovarian_noregress, reduction = "pca", dims = 1:10)
sc_ovarian_noregress <- FindClusters(sc_ovarian_noregress, verbose = FALSE)

sc_ovarian_noregress <- RunUMAP(sc_ovarian_noregress, reduction = "pca", dims = 1:10)
Idents(sc_ovarian_noregress) <- "Phase"
DimPlot(sc_ovarian_noregress, reduction = "umap")



# With regression
sc_ovarian_regress <- SCTransform(sc_ovarian_filtered, 
                                   vars.to.regress = c("S.Score", "G2M.Score"))
sc_ovarian_regress <- RunPCA(sc_ovarian_regress, assay = "SCT", verbose = FALSE)
sc_ovarian_regress <- FindNeighbors(sc_ovarian_regress, reduction = "pca", dims = 1:10)
sc_ovarian_regress <- FindClusters(sc_ovarian_regress, verbose = FALSE)

sc_ovarian_regress <- RunUMAP(sc_ovarian_regress, reduction = "pca", dims = 1:10)
Idents(sc_ovarian_regress) <- "Phase"
DimPlot(sc_ovarian_regress, reduction = "umap")
SpatialDimPlot(sc_ovarian_regress)

# --- Cell Cycle Effect/Regression (ST) ---
sp_ovarian_filtered <- CellCycleScoring(sp_ovarian_filtered, s.features = s.genes, 
                                        g2m.features = g2m.genes, set.ident = TRUE)
# No regression
sp_ovarian_noregress <- SCTransform(sp_ovarian_filtered, assay = "Spatial")
sp_ovarian_noregress <- RunPCA(sp_ovarian_noregress, assay = "SCT", verbose = FALSE)
sp_ovarian_noregress <- FindNeighbors(sp_ovarian_noregress, reduction = "pca", dims = 1:10)
sp_ovarian_noregress <- FindClusters(sp_ovarian_noregress, verbose = FALSE)

sp_ovarian_noregress <- RunUMAP(sp_ovarian_noregress, reduction = "pca", dims = 1:10)
Idents(sp_ovarian_noregress) <- "Phase"
DimPlot(sp_ovarian_noregress, reduction = "umap")
SpatialDimPlot(sp_ovarian_noregress)


# With regression
sp_ovarian_regress <- SCTransform(sp_ovarian_filtered, assay = "Spatial",
                                  vars.to.regress = c("S.Score", "G2M.Score"))
sp_ovarian_regress <- RunPCA(sp_ovarian_regress, assay = "SCT", verbose = FALSE)
sp_ovarian_regress <- FindNeighbors(sp_ovarian_regress, reduction = "pca", dims = 1:10)
sp_ovarian_regress <- FindClusters(sp_ovarian_regress, verbose = FALSE)

sp_ovarian_regress <- RunUMAP(sp_ovarian_regress, reduction = "pca", dims = 1:10)
Idents(sp_ovarian_regress) <- "Phase"
DimPlot(sp_ovarian_regress, reduction = "umap")
SpatialDimPlot(sp_ovarian_regress)

ccplot1 <- DimPlot(sc_ovarian_noregress, reduction = "umap")
ccplot2 <- DimPlot(sc_ovarian_regress, reduction = "umap")
ccplot3 <- DimPlot(sp_ovarian_noregress, reduction = "umap")
ccplot4 <- DimPlot(sp_ovarian_regress, reduction = "umap")

ccpatch <- wrap_plots(ccplot1,
           ccplot2,
           ccplot3,
           ccplot4)

ggsave(filename = "outputs/cell_cycle.jpeg", 
       ccpatch, width = 10,
       height = 6)

