# conda activate ST

library(Seurat)
library(beyondcell)
library(ggplot2)
library(patchwork)

sp_ovarian_annotated <- readRDS("outputs/sp_ovarian_annotated.rds")

# --- Beyondcell (ST only) ---
gs <- GenerateGenesets(SSc)
bc <- bcScore(sp_ovarian_annotated, gs)
 #bc <- bcScore(sp_ovarian_transformed, gs)

bc <- readRDS("outputs/raw_ovarian_bc_object.rds")

# Dimensional reduction
bc <- bcUMAP(bc, pc = 5, res = 0.2)

TCs <- bcClusters(bc,
                  UMAP = "beyondcell", 
                  idents = "bc_clusters_res.0.2", 
                  pt.size = 1.5,
                  spatial = TRUE) #+ 
  ggtitle("Therapeutic Clusters (res = 0.2)")

ggsave(filename = "outputs/sp_TCs", 
       TCs, width = 10,
       height = 6)

# bcranks
bc <- bcRanks(bc, idents = "bc_clusters_res.0.2")
bc_rank_TC0 <- bc4Squares(bc, 
                          idents = "bc_clusters_res.0.2", 
                          lvl = 0,
                          y.cutoff = c(0.1, 0.3, 0.6, 0.9),
                          x.cutoff = c(-70,40))

bc_rank_TC0

#bcSignatures
afatinib_IDs <- FindDrugs(bc, "AFATINIB")$IDs
gefitinib_IDs <- FindDrugs(bc, "GEFITINIB")$IDs

drug_hit_ids <- c("sig_21195", "sig_21407")

drug_hits <- bcSignatures(bc, UMAP = "beyondcell", 
             signatures = list(values = drug_hit_ids), 
             spatial = TRUE, pt.size = 1.5)

wrap_plots(drug_hits)

ggsave(filename = "outputs/sp_ovarian_drughits.jpeg", 
       drug_hits, width = 10,
       height = 6)

# Add nBCS to metadata 
bc <- bcAddMetadata(bc, t(bc@normalized[c("sig_21195", "sig_21407"),]))

saveRDS(bc, file = "outputs/bc_ovarian.rds")

# --- Sub-cluster Analysis of Tumor TC (TC0) ---

TC0 <- bc@meta.data[bc@meta.data$bc_clusters_res.0.2 == 0,]
bc_TC0 <- bcSubset(bc, cells = rownames(TC0)) # Cancer TC
bc_TC0 <- bcUMAP(bc_TC0, pc = 5, res = 0.2)
bc_TC0 <- bcRanks(bc_TC0, idents = "bc_clusters_res.0.2")

# Check, (1559 spots in TC0)
sum(!is.na(bc_TC0@meta.data$bc_clusters_res.0.2))

#Visualize subclusters
bc_TC0_TCs <- bcClusters(bc_TC0, idents = "bc_clusters_res.0.2", spatial = TRUE,
                         label = FALSE, label.size = 4, repel = TRUE) + 
  ggtitle("Therapeutic Sub-Clusters (TC0, 0.2)") #+
#theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
bc_TC0_TCs

ggsave(filename = "outputs/sp_TC0_TCs", 
       bc_TC0_TCs, width = 10,
       height = 6)

bc_TC0 <- bcRanks(bc_TC0, idents = "bc_clusters_res.0.2")

bc_rank_TC0_0 <- bc4Squares(bc_TC0, idents = "bc_clusters_res.0.2", 
                            lvl = 0,
                            y.cutoff = c(0.1, 0.3, 0.6, 0.9),
                            x.cutoff = c(-70,40)) + 
ggtitle("sub-Therapeutic Cluster 0 of TC0 (0.2)") + NoLegend()
bc_rank_TC0_0

bc_rank_TC0_1 <- bc4Squares(bc_TC0, idents = "bc_clusters_res.0.2", 
                            lvl = 1,
                            y.cutoff = c(0.1, 0.3, 0.6, 0.9),
                            x.cutoff = c(-70,40)) + 
  ggtitle("sub-Therapeutic Cluster 1 of TC0 (0.2)") + NoLegend()

bc_rank_TC0_1

wrap_plots(bc_rank_TC0_0, bc_rank_TC0_1)

saveRDS(bc_TC0, file = "outputs/bc_TCO_ovarian.rds")

# Save bc metadata to easily use between different envs
write.csv(bc@meta.data, "outputs/bc_metadata.csv")
write.csv(bc_TC0@meta.data, "outputs/bc_TC0_metadata.csv")





