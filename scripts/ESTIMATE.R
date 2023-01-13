# conda activate ESTIMATE

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, 
                 dependencies=TRUE)
library(estimate)
library(Seurat)
library(beyondcell)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# Normalized ST Seurat object
sp_ovarian_transformed <- 
  readRDS('outputs/transformed_seurat_sp_ovarian.rds')

# Main BC and TC0 BC object (beyondcell.R)
bc <- readRDS("outputs/bc_ovarian.rds")
bc_TC0 <- readRDS("outputs/bc_TCO_ovarian.rds")

# --- ESTIMATE ---

# Create input GCT file
outputGCT(as.data.frame(bc@expr.matrix),
          "ovarian_estimate.gct")

outputGCT(as.data.frame(sp_ovarian_filtered@assays$Spatial@counts),
          "ovarian_estimate.gct")

# filtercommongenes troubleshoot
#https://rdrr.io/rforge/estimate/src/R/filterCommonGenes.R

filterCommonGenes_v2 <- function(input.f,
                                 output.f,
                                 id=c("GeneSymbol", "EntrezID")) {
  
  ## Check arguments
  stopifnot((is.character(input.f) && length(input.f) == 1 && nzchar(input.f)) ||
              (inherits(input.f, "connection") && isOpen(input.f, "r")))
  stopifnot((is.character(output.f) && length(output.f) == 1 && nzchar(output.f)))
  id <- match.arg(id)   
  
  ## Read input data
  input.df <- read.table(input.f,
                         header=TRUE,
                         row.names=1,
                         sep="\t", 
                         quote="",
                         stringsAsFactors=FALSE)
  # SAVE CELL NAMES
  colnames(input.df) <- input.df["NAME",]
  
  merged.df <- merge(common_genes, input.df, by.x=id, by.y="row.names")
  rownames(merged.df) <- merged.df$GeneSymbol
  merged.df <- merged.df[, -1:-ncol(common_genes)]
  print(sprintf("Merged dataset includes %d genes (%d mismatched).",
                nrow(merged.df),
                nrow(common_genes) - nrow(merged.df)))
  outputGCT(merged.df, output.f)
}

filterCommonGenes_v2("ovarian_estimate.gct",
                     "ESTIMATE_output_sp_ovarian.gct",
                     id = "GeneSymbol")

estimateScore("ESTIMATE_output_sp_ovarian.gct",
              "ESTIMATE_scores_sp_ovarian.tsv",
              platform = "illumina")

estimate_scores <- read.table("ESTIMATE_scores_sp_ovarian.tsv",
                              sep = "\t") 

# BC vs ESTIMATE
# Add nBCS to metadata 

bc <- bcAddMetadata(bc, t(bc@normalized[c("sig_21195", "sig_21407"),]))
TC0 <- bc@meta.data[bc@meta.data$bc_clusters_res.0.2 == 0,]
bc_TC0 <- bcSubset(bc, cells = rownames(TC0)) # Cancer TC
bc_TC0 <- bcUMAP(bc_TC0, pc = 5, res = 0.2)

saveRDS(bc_TC0, file = "outputs/bc_TC0_estimate.rds")
bc_TC0 <- readRDS("outputs/bc_TC0_estimate.rds")

write.csv(bc_TC0@meta.data, "outputs/bc_TC0_metadata.csv")

stromal <- SpatialFeaturePlot(sp_ovarian_filtered, 
                              features = "StromalScore") + theme(legend.position = "right")

immune <- SpatialFeaturePlot(sp_ovarian_filtered, 
                             features = "ImmuneScore") + theme(legend.position = "right")

ESTIMATE <- SpatialFeaturePlot(sp_ovarian_filtered, 
                               features = "ESTIMATEScore") + theme(legend.position = "right")

wrap_plots(stromal, immune, ESTIMATE, ncol = 1)

# --- ESTIMATE score vs nBCS---

gef_stromal <- ggplot(bc_TC0@meta.data %>% drop_na(bc_clusters_res.0.2), 
                      aes(x = sig_21195, y = StromalScore, col = bc_clusters_res.0.2)) +
  geom_smooth() + NoLegend() + stat_compare

gef_immune <- ggplot(bc_TC0@meta.data %>% drop_na(bc_clusters_res.0.2), 
                     aes(x = sig_21195, y = ImmuneScore, col = bc_clusters_res.0.2)) +
  geom_smooth() + NoLegend()

gef_estimate <- ggplot(bc_TC0@meta.data %>% drop_na(bc_clusters_res.0.2), 
                       aes(x = sig_21195, y = ESTIMATEScore, col = bc_clusters_res.0.2)) +
  geom_smooth() + NoLegend()

afa_stromal <- ggplot(bc_TC0@meta.data %>% drop_na(bc_clusters_res.0.2), 
                      aes(x = sig_21407, y = StromalScore, col = bc_clusters_res.0.2)) +
  geom_smooth()

afa_immune <- ggplot(bc_TC0@meta.data %>% drop_na(bc_clusters_res.0.2), 
                     aes(x = sig_21407, y = ImmuneScore, col = bc_clusters_res.0.2)) +
  geom_smooth()

afa_estimate <- ggplot(bc_TC0@meta.data %>% drop_na(bc_clusters_res.0.2), 
                       aes(x = sig_21407, y = ESTIMATEScore, col = bc_clusters_res.0.2)) +
  geom_smooth()

wrap_plots(gef_stromal,
           afa_stromal,
           gef_immune,
           afa_immune,
           gef_estimate,
           afa_estimate, ncol = 2)

ggplot(bc_TC0@meta.data %>% drop_na(bc_clusters_res.0.2), 
       aes(x = bc_clusters_res.0.2, y = StromalScore)) +
  geom_boxplot() + NoLegend()

ggplot(bc_TC0@meta.data %>% drop_na(bc_clusters_res.0.2), 
       aes(x = bc_clusters_res.0.2, y = ImmuneScore)) +
  geom_boxplot() + NoLegend()

ggplot(bc_TC0@meta.data %>% drop_na(bc_clusters_res.0.2), 
       aes(x = bc_clusters_res.0.2, y = ESTIMATEScore)) +
  geom_boxplot() + NoLegend()

ggsave(filename = "outputs/sp_estimate_results2.jpeg", 
       wrap_plots(stromal, immune, ESTIMATE, ncol = 1), 
       width = 10,
       height = 6)

write.csv(estimate_scores, file = "outputs/estimate_scores.csv")

Cellcomposition <- function(estimate) {
  ### Takes ESTIMATE output
  # Scale ESTIMATEScore
  estimate <- estimate %>%
    mutate(ESTIMATEScore = scales::rescale(ESTIMATEScore, to = c(0, 1)))
  # Compute signs
  estimate <- estimate %>% 
    mutate(sign = sign(ImmuneScore) + sign(StromalScore),
           sign = case_when(sign == 0 & sign(ImmuneScore) > 0 ~ 1,
                            sign == 0 & sign(StromalScore) > 0 ~ -1,
                            sign != 0 ~ sign),
           greater.than = abs(ImmuneScore) > abs(StromalScore))
  # Compute log2
  estimate <- estimate %>%
    mutate(log2 = 
             log2(case_when(sign == 2 | (sign == 1 & greater.than) | 
                              (sign == -1 & !greater.than) ~ 
                              abs(ImmuneScore/StromalScore),
                            sign == -2 | (sign == 1 & !greater.than) |
                              (sign == -1 & greater.than) ~ 
                              abs(StromalScore/ImmuneScore))))
  
  # Compute composition
  estimate <- estimate %>%
    mutate(Composition = case_when(ESTIMATEScore <= 0.25 ~ "High-Tumor Purity",
                                   ESTIMATEScore > 0.75 ~ "Low-Tumor Purity",
                                   TRUE ~ "Med-Tumor Purity"))
  return(estimate$Composition)
  
}

# Compute cell composition
composition <- Cellcomposition(estimate_scores)

# Rescale scores?
estimate_scores <- apply(estimate_scores, 2, 
                         scales::rescale, to = c(0, 1)) %>%
  as.data.frame %>% mutate(Composition = composition)




