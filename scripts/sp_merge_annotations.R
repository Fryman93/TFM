# conda activate ST
rm(list = ls())
library(Seurat) # Seurat_4.1.1
library(tidyverse) # tidyverse_1.3.2
library(dplyr) # dplyr_1.0.9

# --- Merge annotations (Only for ST dataset)---
sp_ovarian_transformed <- 
  readRDS('outputs/transformed_seurat_sp_ovarian.rds')

# CopyKAT
sp_copy <- readRDS("sp_ovarian/seurat_object_ovarian_CopyKAT.rds")
View(sp_copy@meta.data)
ploidy <- sp_copy@meta.data %>% 
  filter(row.names(sp_copy@meta.data) %in% 
           row.names(sp_ovarian_filtered@meta.data)) %>% 
  select(copykat.pred)

sp_ovarian_transformed <- AddMetaData(sp_ovarian_filtered, ploidy)

# CopyKAT Check 
Idents(sp_ovarian_transformed) <- "copykat.pred"
DimPlot(sp_ovarian_transformed, reduction = "umap")
SpatialDimPlot(sp_ovarian_transformed)

# ESTIMATE
estimate_scores <- read.table("ESTIMATE_scores_sp_ovarian.tsv",
                              sep = "\t") 

# Transform estimate_scores to facilitate merge
rownames(estimate_scores) <- estimate_scores[,"V1"]
colnames(estimate_scores) <- gsub("\\.", "-", 
                                  estimate_scores["NAME",])

estimate_scores <- estimate_scores[-1,]
estimate_scores[,1:3] <- NULL
estimate_scores <- as.data.frame(t(estimate_scores))
estimate_scores$NAME <- NULL

# Change data types
estimate_scores$StromalScore <- as.numeric(estimate_scores$StromalScore)
estimate_scores$ImmuneScore <- as.numeric(estimate_scores$ImmuneScore)
estimate_scores$ESTIMATEScore <- as.numeric(estimate_scores$ESTIMATEScore)

# Merge estimate scores
sp_ovarian_transformed@meta.data$StromalScore <- estimate_scores$StromalScore
sp_ovarian_transformed@meta.data$ImmuneScore <- estimate_scores$ImmuneScore
sp_ovarian_transformed@meta.data$ESTIMATEScore <- estimate_scores$ESTIMATEScore

# ESTIMATE Check
stromal <- SpatialFeaturePlot(sp_ovarian_transformed, 
                              features = "StromalScore") + 
  theme(legend.position = "right")

immune <- SpatialFeaturePlot(sp_ovarian_transformed, 
                             features = "ImmuneScore") + 
  theme(legend.position = "right")

ESTIMATE <- SpatialFeaturePlot(sp_ovarian_transformed, 
                               features = "ESTIMATEScore") + 
  theme(legend.position = "right")

estimate_plots <- wrap_plots(stromal, 
                             immune, 
                             ESTIMATE, 
                             ncol = 3)
ggsave(filename = "outputs/estimate_plots.jpeg", 
       estimate_plots, width = 10,
       height = 6)

# CARD
card_prop <- read.csv("outputs/prop_CARD_sp_ovarian.csv")
rownames(card_prop) <- card_prop$X
card_prop <- card_prop %>% subset(select = -c(X))

sp_ovarian_annotated <- AddMetaData(sp_ovarian_transformed, card_prop)
saveRDS(sp_ovarian_annotated, file = "outputs/sp_ovarian_annotated.rds")

