# conda activate ST
rm(list = ls())

library(tidyverse)
library(hdf5r)
library(Seurat) # Seurat_4.1.1
library(beyondcell) # beyondcell_1.3.3
library(CARD) # CARD_1.0.0
library(MuSiC) # MuSiC_0.2.0
library(patchwork) # patchwork_1.1.1
library(ggplot2) # ggplot2_3.3.5
library(dplyr) # dplyr_1.0.9
library(SummarizedExperiment) # SummarizedExperiment_1.24.0
library(SingleCellExperiment) # SingleCellExperiment_1.16.0
library(tidyverse) # tidyverse_1.3.2
library(glue) # glue_1.6.2
library(SingleR) # SingleR_1.8.1
library(celldex) # celldex_1.4.0

# --- Single Cell Data (ovarian GSE158937) ---
setwd("/disk1/beyondcell/ovarian/")
cts <- ReadMtx(mtx = 'inputs/matrix.mtx',
               features = 'inputs/features.tsv',
               cells = 'inputs/barcodes.tsv',
               feature.column = 2)
sc_ovarian <- CreateSeuratObject(counts =  cts)
saveRDS(sc_ovarian, file = 'outputs/raw_seurat_sc_ovarian.rds')

# --- Single Cell QC ---
dim(sc_ovarian@meta.data)
# sc_ovarian before: 7123 single cells
sc_ovarian[["percent.mt"]] <- PercentageFeatureSet(sc_ovarian, pattern = "^MT-")
sc_ovarian[["percent.rb"]] <- PercentageFeatureSet(sc_ovarian, pattern="^RP[LS]")

qc_sc_metrics_before <- VlnPlot(sc_ovarian,
                                features = c("nFeature_RNA", 
                                             "nCount_RNA", 
                                             "percent.mt",
                                             "percent.rb"), ncol = 4)


sc_ovarian_filtered <- subset(sc_ovarian, subset = nFeature_RNA > 200 & 
                                nFeature_RNA < 6000 &
                                nCount_RNA < 60000 &
                                percent.mt < 25 &
                                percent.rb < 40)

qc_sc_metrics_after <- VlnPlot(sc_ovarian_filtered, 
                               features = c("nFeature_RNA", 
                                            "nCount_RNA", 
                                            "percent.mt",
                                            "percent.rb"), ncol = 4)
# sc_ovarian(_filtered) after: 3272 single cells
dim(sc_ovarian_filtered@meta.data)
sc_patch1 <- qc_sc_metrics_before / qc_sc_metrics_after

ggsave(qc_sc_metrics_before / qc_sc_metrics_after,
       filename = "outputs/qc_sc_before_after.jpeg")

# Save the object at this point for cell-cycle check 
saveRDS(sc_ovarian_filtered, 
        file = 'outputs/seurat_sc_ovarian_filtered.rds')

# Load object if necessary
sc_ovarian_filtered <- 
  readRDS("/disk1/beyondcell/ovarian/outputs/seurat_sc_ovarian_filtered.rds")

# Normalization/Clustering/UMAP
sc_ovarian_filtered <- SCTransform(sc_ovarian_filtered)
sc_ovarian_filtered <- RunPCA(sc_ovarian_filtered, assay = "SCT", verbose = FALSE)
DimPlot(sc_ovarian_filtered, reduction = "pca")
ElbowPlot(sc_ovarian_filtered)

# Account for results from elbow plot (10 PCs)
sc_ovarian_filtered <- FindNeighbors(sc_ovarian_filtered, reduction = "pca", dims = 1:10)
sc_ovarian_filtered <- FindClusters(sc_ovarian_filtered, resolution = 0.2, verbose = FALSE)
sc_ovarian_filtered <- RunUMAP(sc_ovarian_filtered, reduction = "pca", dims = 1:10)
Idents(sc_ovarian_filtered) <- "seurat_clusters"
DimPlot(sc_ovarian_filtered, reduction = "umap")

# SingleR Cell Annotation (sc)
sce_filtered <- as.SingleCellExperiment(DietSeurat(sc_ovarian_filtered))
ref <- HumanPrimaryCellAtlasData()

# Remove irrelevant cell types from ref
remove_types <- c("BM & Prog.", "Neurons", "Embryonic_stem_cells", "HSC_-G-CSF",
                  "Osteoblast", "BM", "iPS_cells", "Hepatocytes", 
                  "Neuroepithelial_cell", "Astrocyte", "HSC_CD34+", "MSC")
ref <- ref[, !(ref$label.main %in% remove_types)]


seurat_clusters_factor <- setNames(sc_ovarian_filtered@meta.data$seurat_clusters, 
                                   rownames(sc_ovarian_filtered@meta.data))
pred.grun <- SingleR(test = sce_filtered, 
                     ref = ref, labels = ref$label.main, 
                     de.method = "wilcox",
                     clusters = seurat_clusters_factor)

# Transform SingleR column to add to Seurat metadata
pred.grun <- as.data.frame(pred.grun) %>% 
  select(pruned.labels) %>% 
  mutate(seurat_clusters = rownames(pred.grun))
colnames(pred.grun)[1] <- "immunelabels"
newmetadata <- sc_ovarian_filtered@meta.data %>% select(seurat_clusters) %>%
  mutate(cells = rownames(sc_ovarian_filtered@meta.data)) %>%
  merge(pred.grun)
rownames(newmetadata) <- newmetadata$cells
newmetadata <- newmetadata %>% select(-cells)
sc_ovarian_filtered <- AddMetaData(sc_ovarian_filtered, newmetadata)

#Save object here to use with Copykat.R
saveRDS(sc_ovarian_filtered, file = "outputs/sc_ovarian_singleR.rds")

# Seurat UMAP by exp cluster
Idents(sc_ovarian_filtered) <- "seurat_clusters"
DimPlot(sc_ovarian_filtered, reduction = "umap")

seurat_umap <- DimPlot(sc_ovarian_filtered, reduction = "umap")

# Seurat UMAP by SingleR
Idents(sc_ovarian_filtered) <- "immunelabels"
DimPlot(sc_ovarian_filtered, reduction = "umap")
annotated_umap <- DimPlot(sc_ovarian_filtered, reduction = "umap")

# Seurat UMAP by immune.labels with copykat (Need to run Copykat.R first)
sc_ovarian_annotated_copykat <- 
  readRDS("outputs/sc_seurat_object_ovarian_CopyKAT.rds")

Idents(sc_ovarian_annotated_copykat) <- "immunelabels.copykat"
DimPlot(sc_ovarian_annotated_copykat, reduction = "umap")

annotated_umap2 <- DimPlot(sc_ovarian_annotated_copykat, reduction = "umap")

# Seurat UMAP by Copykat
Idents(sc_ovarian_annotated_copykat) <- "copykat.pred"
DimPlot(sc_ovarian_annotated_copykat, reduction = "umap")

umap_copyk <- DimPlot(sc_ovarian_annotated_copykat, reduction = "umap")

# Seurat UMAP by scaled WFDC2 expression
FeaturePlot(sc_ovarian_annotated_copykat, features = "EPCAM")
sc_WFDC2 <- FeaturePlot(sc_ovarian_annotated_copykat, features = "WFDC2")


sc_umaps <- wrap_plots(seurat_umap, 
                       annotated_umap2,
                       umap_copyk,
                       sc_WFDC2,
                       ncol = 2)

ggsave(filename = "outputs/sc_umaps.jpeg", 
       sc_umaps, width = 10,
       height = 6)

# View most frequent cell types per cluster
most_freq_cells <- sc_ovarian_filtered@meta.data %>% 
  dplyr::group_by(seurat_clusters) %>%
  dplyr::count(immunelabels) %>% 
  dplyr::arrange(seurat_clusters, desc(n)) %>% 
  dplyr::top_n(3)
View(most_freq_cells)

# View # cells per cluster
total_cells <- sc_ovarian_filtered@meta.data %>% 
  dplyr::group_by(seurat_clusters) %>%
  dplyr::count()
View(total_cells)

# Gene Expression exploration
cancer_markers <- c("IFI6", "SLC3A1", "PEG10",
                    "MUC16", "WFDC2", "PROM1",
                    "CD44", "EPCAM", "ALDH1A1")
tcell_markers <- c("CD3E", "CD8A", "FOXP3",
                   "GZMK", "GZMB", "CD4",
                   "FOXP3")
bcell_markers <- c("CD79A", "SDC1")
macro_markers <- c("ITGAM")

tam_markers <- c("CD68", "CD80", "CD86",
                 "CD32", "CD163", "CD24",
                 "CSF1")
markers <- c("WFDC2", "CCL5", "IER3",
             "MALAT1","UBE2C", "IL1B",
             "ESM1", "COL1A1", "IGHG1")
immune_markers <- c("PDGFRA", "RGS5", "PECAM1",
                    "CD14", "CD3E", "CD79A", 
                    "SDC1", "CD79A")
immune_markers2 <- c("CD68", "CD80", "CD86",
                     "CD163", "CD24", "CSF1",
                     "TREM1", "TREM2", "CD8")

FeaturePlot(sc_ovarian_filtered, features = "WFDC2")
FeaturePlot(sc_ovarian_filtered, features = tcell_markers)
FeaturePlot(sc_ovarian_filtered, features = bcell_markers)
FeaturePlot(sc_ovarian_filtered, features = macro_markers)
FeaturePlot(sc_ovarian_filtered, features = tam_markers)
FeaturePlot(sc_ovarian_filtered, features = c("IL1B", "TGFB1"))

# At this point, assuming copykat annotations were applied,
# sc dataset should be fully annotated
saveRDS(sc_ovarian_filtered, file = "outputs/sc_ovarian_annotated.rds")
sc_ovarian_filtered



# --- Spatial Data (Human Ovarian Cancer 10x Genomics) ---
setwd("../ovarian")
spatial_sc <- basename(list.files(pattern = "*bc_matrix.h5$",
                                  recursive = TRUE, 
                                  full.names = TRUE))

image <- Read10X_Image(image.dir = paste0(getwd(), "/inputs/images/spatial"),
                       image.name = "tissue_lowres_image.png")

sp_ovarian <- Load10X_Spatial(data.dir = paste0(getwd(), "/inputs"), 
                              filename = spatial_sc, image = image)

# Generate basic Seurat clusters first
# Visualize metrics
Idents(sp_ovarian) <- 'orig.ident'
sp_plot1 <- VlnPlot(sp_ovarian, features = "nCount_Spatial", pt.size = 0.1)
+ NoLegend()
sp_plot2 <- SpatialFeaturePlot(sp_ovarian, features = "nCount_Spatial") + 
  theme(legend.position = "right")
sp_plot3 <- VlnPlot(sp_ovarian, features = "nFeature_Spatial", pt.size = 0.1) 
+ NoLegend()
sp_plot4 <- SpatialFeaturePlot(sp_ovarian, features = "nFeature_Spatial") + 
  theme(legend.position = "right")
sp_patch1 <- (sp_plot1 | sp_plot2) / (sp_plot3 | sp_plot4)

sp_ovarian[["percent.mt"]] <- PercentageFeatureSet(sp_ovarian, pattern = "^MT-")
sp_ovarian[["percent.rb"]] <- PercentageFeatureSet(sp_ovarian, pattern = "^RP[LS]")

# 10x datasets are usually preprocessed
# Spatial dataset JUST needs to be filtered
# Filter based on qc_metrics
qc_sp_metrics_before <- VlnPlot(sp_ovarian, 
                                features = c("nFeature_Spatial", 
                                             "nCount_Spatial", 
                                             "percent.mt",
                                             "percent.rb"), ncol = 4)

sp_ovarian_filtered <- subset(sp_ovarian, subset = nFeature_Spatial > 200 & 
                                nFeature_Spatial < 11000 &
                                nCount_Spatial < 80000)

qc_sp_metrics_after <- VlnPlot(sp_ovarian_filtered, features = c("nFeature_Spatial", 
                                                                 "nCount_Spatial", 
                                                                 "percent.mt",
                                                                 "percent.rb"), ncol = 4)
sp_patch2 <- qc_sp_metrics_before / qc_sp_metrics_after

ggsave(filename = "outputs/qc_sp_before_after.jpeg", 
       sp_patch2, width = 10,
       height = 6)

# Save the object at this point for cell-cycle check
saveRDS(sp_ovarian_filtered, file = "outputs/filtered_seurat_sp_ovarian.rds")
sp_ovarian_filtered <- 
  readRDS("outputs/filtered_seurat_sp_ovarian.rds")

sp_ovarian_transformed <- SCTransform(sp_ovarian_filtered, assay = 'Spatial')
sp_ovarian_transformed <- RunPCA(sp_ovarian_transformed, assay = "SCT", verbose = FALSE)
ElbowPlot(sp_ovarian_transformed)
sp_ovarian_transformed <- FindNeighbors(sp_ovarian_transformed, reduction = "pca", dims = 1:10)
sp_ovarian_transformed <- FindClusters(sp_ovarian_transformed, verbose = FALSE)

sp_ovarian_transformed <- RunUMAP(sp_ovarian_transformed, reduction = "pca", dims = 1:10)

Idents(sp_ovarian_transformed) <- "seurat_clusters"
DimPlot(sp_ovarian_transformed, reduction = "umap")
SpatialDimPlot(sp_ovarian_transformed)

# UMAPs and Exp Clusters
sp_plot5 <- DimPlot(sp_ovarian_transformed, reduction = "umap", label = TRUE)
sp_plot6 <- SpatialDimPlot(sp_ovarian_transformed, label = TRUE, label.size = 3)
sp_patch3 <- sp_plot5 + sp_plot6
ggsave(filename = "outputs/umap_vs_sp_ovarian.png", 
       sp_patch3, width = 10,
       height = 6)

# Ovarian Cancer Biomarkers
WFDC2 <- SpatialFeaturePlot(sp_ovarian_transformed, features = "WFDC2")

# Other Biomarkers
SpatialFeaturePlot(sp_ovarian_transformed, features = c("EGFR", "HER1"))
SpatialFeaturePlot(sp_ovarian_transformed, features = c("COL1A1", "IGHG1"))
SpatialFeaturePlot(sp_ovarian_transformed, features = tcell_markers)
SpatialFeaturePlot(sp_ovarian_transformed, features = bcell_markers)
SpatialFeaturePlot(sp_ovarian_transformed, features = immune_markers)
SpatialFeaturePlot(sp_ovarian_transformed, features = immune_markers2)

# TME Markers
# Tumor cells: EPCAM+, 
# Immune cells: PTPRC+,
# Stromal cells: EPCAM-, PTPRC- 

sp_patch3 <- SpatialFeaturePlot(sp_ovarian_transformed, 
                                features = c("EPCAM", "PTPRC"))
ggsave(filename = "outputs/sp_tme_markers.jpeg", 
       sp_patch3, width = 10,
       height = 6)

Idents(sp_ovarian_annotated) <- "seurat_clusters"
SpatialDimPlot(sp_ovarian_annotated)
sp_seurat <- SpatialDimPlot(sp_ovarian_annotated)

# Save the transformed object to set up for more annotations 
saveRDS(sp_ovarian_transformed, 
        file = 'outputs/transformed_seurat_sp_ovarian.rds')

# Run Copykat.R on ST dataset to get ploidy status in metadata 
Idents(sp_ovarian_annotated) <- "copykat.pred"
SpatialDimPlot(sp_ovarian_annotated)
sp_copyk <- SpatialDimPlot(sp_ovarian_annotated)

sp_tumor <- wrap_plots(TCs, 
                       sp_seurat,
                       WFDC2, 
                       sp_copyk, 
                       ncol = 2)

ggsave(filename = "outputs/sp_tumor.jpeg", 
       sp_tumor, width = 10,
       height = 6)