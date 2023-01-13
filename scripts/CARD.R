# conda activate ST

library(Seurat) # Seurat_4.1.1
library(CARD) # CARD_1.0.0
library(MuSiC) # MuSiC_0.2.0

sc_ovarian_annotated <- 
  readRDS('outputs/sc_ovarian_annotated.rds')

#scRNA with immunelabels.copykat
sc_ovarian_annotated_copykat <- 
  readRDS('outputs/sc_seurat_object_ovarian_CopyKAT.rds')

sp_ovarian_transformed <- 
  readRDS('outputs/transformed_seurat_sp_ovarian.rds')


# Set up the required CARD parameters
sc_count = sc_ovarian_annotated_copykat@assays$RNA@counts
sc_meta = sc_ovarian_annotated_copykat@meta.data[colnames(sc_ovarian_annotated_copykat@assays$RNA@counts),]
spatial_count = sp_ovarian_transformed@assays$Spatial@counts
spatial_location = sp_ovarian_transformed@images$slice1@coordinates[, c('imagerow','imagecol')]

# Check if colData and rowData are aligned
identical(colnames(sc_count), rownames(sc_meta))
identical(colnames(spatial_count), rownames(spatial_location))
colnames(spatial_location) <- c('x', 'y')
View(sc_count)

# Define cell type we want to deconvolute
top_cell_types <- c("Macrophage", "Smooth_muscle_cells", "Endothelial_cells",
                     "B_cell", "Epithelial_cells",
                    "T_cells")

### HAVE ALL CARD INPUTS DEFINED #### 
CARD_obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "immunelabels.copykat",
  ct.select = unique(sc_meta$immunelabels.copykat),
  sample.varname = "orig.ident",
  minCountGene = 100,
  minCountSpot = 5)

CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

p2 <- CARD.visualize.prop(
  proportion = CARD_obj@Proportion_CARD,        
  spatial_location = CARD_obj@spatial_location, 
  ct.visualize = unique(sc_meta$immunelabels.copykat),
  NumCols = 3, colors = c("#f7edd0", "#EA4F88", "#4B2991")) 

p2$layers[[1]]$aes_params$size <- 1
p2 <- p2 + coord_flip() + scale_x_reverse()

p2
write.csv(CARD_obj@Proportion_CARD, 
          "/disk1/beyondcell/ovarian/outputs/prop_CARD_sp_ovarian.csv")

ggsave(filename = "outputs/card_sp_cell_comps.jpeg", p2, 
       width = 10, height = 6)


