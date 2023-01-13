# conda activate scCNV
rm(list = ls()) # R version 4.1.1
library(copykat) # copykat_1.0.5
library(dplyr) # dplyr_1.0.7
library(tibble)
library(forcats) # forcats_0.5.1
library(celldex)
library(SingleR)

# --- Data ---
# Seurat object (scRNA)
sc_ovarian_singler <- readRDS("outputs/sc_ovarian_singleR.rds")

# Raw expression matrix
sc_raw_counts <- as.matrix(sc_ovarian_singler@assays$RNA@counts)

# Run copykat
setwd(outdir)
copyk <- copykat(rawmat = sc_raw_counts, id.type = "Symbol", ngene.chr = 5,
                 win.size = 25, KS.cut = 0.1, sam.name = "sc_ovarian", 
                 distance = "euclidean", norm.cell.names = "", n.cores = 4, 
                 output.seg = FALSE)
# If cell is epithelial and aneuploid, label it as cancer epithelial
condition <- 

copyk@meta.data <- 
          copyk@meta.data %>% 
          mutate(immunelabels.copykat = 
            ifelse(copykat.pred == "Epithelial_cells" &
                   immunelabels == "Aneuploid", "Cancer_epithelial_cells",
                   immunelabels))

saveRDS(copyk, file = "outputs/sc_seurat_object_ovarian_CopyKAT.rds")

# Seurat object (ST)
sp_ovarian_transformed <- readRDS("outputs/transformed_seurat_sp_ovarian.rds")

# Raw expression matrix
raw_counts <- as.matrix(sp_ovarian_transformed@assays$Spatial@counts)

# Run copykat
setwd(outdir)
copyk <- copykat(rawmat = raw_counts, id.type = "Symbol", ngene.chr = 5,
                 win.size = 25, KS.cut = 0.1, sam.name = "sp_ovarian", 
                 distance = "euclidean", norm.cell.names = "", n.cores = 4, 
                 output.seg = FALSE)
# saveRDS(copyk, file = "copykat.rds")

saveRDS(copyk, "sp_seurat_object_ovarian_CopyKAT.rds")

