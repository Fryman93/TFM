# conda activate plots

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

setwd("/disk1/beyondcell/ovarian/")

metadata <- read.csv("outputs/bc_TC0_metadata.csv")
metadata$bc_clusters_res.0.2 <- as.factor(metadata$bc_clusters_res.0.2)
rownames(metadata) <- metadata$X
metadata <- metadata %>% subset(select = -c(X))

metadata <- metadata %>% rename(cell = X, subTC = bc_clusters_res.0.2)
mainmetadata <- mainmetadata %>% rename(cell = Row.names, TC = bc_clusters_res.0.2)


mainmetadata <- read.csv("outputs/bc_metadata.csv")
rownames(mainmetadata) <- mainmetadata$X
mainmetadata <- mainmetadata %>% 
  subset(select = -c(X, Macrophage, DC, Fibroblasts, Endothelial_cells,
                     B_cell, Epithelial_cells, T_cells))


card_prop <- read.csv("outputs/prop_CARD_sp_ovarian.csv")
rownames(card_prop) <- card_prop$X
card_prop <- card_prop %>% subset(select = -c(X))

metadata <- metadata %>% 
  subset(select = -c(Macrophage, DC, Fibroblasts, Endothelial_cells,
                     B_cell, Epithelial_cells, T_cells))

View(card_prop[rownames(card_prop) %in% rownames(metadata),])

metadata <- merge(metadata, 
                  card_prop[rownames(card_prop) %in% rownames(metadata),],
                  by = 'row.names', all = TRUE)

mainmetadata <- merge(mainmetadata, 
                  card_prop,
                  by = 'row.names', all = TRUE)

# Comparisons
comp <- list(c("0", "1"))

labels <- list(c("0", "1"), c("0", "2"), c("1", "2"))


# ESTIMATE lineplots of nBCS vs ESTIMATE scores
gef_stromal <- ggplot(metadata %>% drop_na(bc_clusters_res.0.2), 
                      aes(x = sig_21195, 
                          y = StromalScore, 
                          col = bc_clusters_res.0.2)) +
                      geom_smooth() + theme(legend.position="none")

gef_immune <- ggplot(metadata %>% drop_na(bc_clusters_res.0.2), 
                     aes(x = sig_21195, 
                         y = ImmuneScore, 
                         col = bc_clusters_res.0.2)) +
                    geom_smooth() + theme(legend.position="none")

gef_estimate <- ggplot(metadata %>% drop_na(bc_clusters_res.0.2), 
                       aes(x = sig_21195, 
                           y = ESTIMATEScore, 
                           col = bc_clusters_res.0.2)) +
                       geom_smooth() + theme(legend.position="none")

afa_stromal <- ggplot(metadata %>% drop_na(bc_clusters_res.0.2), 
                      aes(x = sig_21407, 
                          y = StromalScore, 
                          col = bc_clusters_res.0.2)) +
                      geom_smooth() #+ theme(legend.position="none")

afa_immune <- ggplot(metadata %>% drop_na(bc_clusters_res.0.2), 
                     aes(x = sig_21407, 
                         y = ImmuneScore, 
                         col = bc_clusters_res.0.2)) +
                     geom_smooth() #+ theme(legend.position="none")

afa_estimate <- ggplot(metadata %>% drop_na(bc_clusters_res.0.2), 
                       aes(x = sig_21407, 
                           y = ESTIMATEScore, 
                           col = bc_clusters_res.0.2)) +
                       geom_smooth() #+ theme(legend.position="none")

wrap_plots(gef_stromal,
           afa_stromal,
           gef_immune,
           afa_immune,
           gef_estimate,
           afa_estimate, ncol = 2)

sub_stromal <- ggplot(metadata %>% drop_na(bc_clusters_res.0.2),
       aes(x = bc_clusters_res.0.2, 
           y = StromalScore)) +
       geom_boxplot() + theme(legend.position="none") +
       stat_compare_means(aes(group = bc_clusters_res.0.2), 
                              label = "p.format",
                              method = "t.test",
                              #method.args = list(var.equal = FALSE,
                              #               alternative = "greater"),
                              label.x = 1.4, label.y = 1000)

sub_immune <- ggplot(metadata %>% drop_na(bc_clusters_res.0.2),
       aes(x = bc_clusters_res.0.2, 
           y = ImmuneScore)) +
  geom_boxplot() + theme(legend.position="none") +
  stat_compare_means(aes(group = bc_clusters_res.0.2), 
                     label = "p.format",
                     method = "t.test",
                     #method.args = list(var.equal = FALSE,
                     #                    alternative = "greater"),
                     label.x = 1.4, label.y = 3000)

sub_estimate <- ggplot(metadata %>% drop_na(bc_clusters_res.0.2),
       aes(x = bc_clusters_res.0.2, 
           y = ESTIMATEScore)) +
  geom_boxplot() + theme(legend.position="none") +
  stat_compare_means(aes(group = bc_clusters_res.0.2), 
                     label = "p.format",
                     method = "t.test",
                     #method.args = list(var.equal = FALSE,
                     #                    alternative = "greater"),
                     label.x = 1.4, label.y = 3500)

sub_estimate_plots <- wrap_plots(sub_stromal, 
                                 sub_immune, 
                                 sub_estimate)

ggsave(filename = "outputs/sub_estimate_plots.jpeg", 
       sub_estimate_plots, width = 10,
       height = 6)

# CARD Boxplots (subTCs)

cellgroup1 <- c("B_cell", "Macrophage", "Smooth_muscle_cells")
cellgroup2 <- c("Cancer_epithelial_cells", "Normal_epithelial_cells", 
                "T_cells", "Endothelial_cells")

cellgroup3 <- c("B_cell", "Macrophage", "Endothelial_cells")
cellgroup4 <- c("Cancer_epithelial_cells", "Normal_epithelial_cells", 
                "T_cells", "Smooth_muscle_cells")

CARD1 <- metadata %>% 
  select(cell, subTC, Normal_epithelial_cells:Smooth_muscle_cells) %>%
  pivot_longer(cols = cellgroup1, names_to = "CARD1", 
               values_to = "Composition")

pcard <- ggboxplot(CARD1, x = "subTC", y = "Composition", fill = "subTC") + 
  stat_compare_means(label = "p.signif", comparisons = comp) +
  theme(legend.position = "bottom") +
  facet_wrap(~CARD1)

pcard

CARD2 <- metadata %>% 
  select(cell, subTC, Normal_epithelial_cells:Smooth_muscle_cells) %>%
  pivot_longer(cols = cellgroup3, names_to = "CARD2", 
               values_to = "Composition")

pcard2 <- ggboxplot(CARD2, x = "subTC", y = "Composition", fill = "subTC") + 
  stat_compare_means(label = "p.signif", comparisons = comp) +
  theme(legend.position = "bottom") +
  facet_wrap(~CARD2)

CARD3 <- metadata %>% 
  select(cell, subTC, Normal_epithelial_cells:Smooth_muscle_cells) %>%
  pivot_longer(cols = cellgroup4, names_to = "CARD3", 
               values_to = "Composition")

pcard3 <- ggboxplot(CARD3, x = "subTC", y = "Composition", fill = "subTC") + 
  stat_compare_means(label = "p.signif", comparisons = comp) +
  theme(legend.position = "bottom") +
  facet_wrap(~CARD3)

ggsave(filename = paste0(outdir, "ESTIMATE_boxplots.pdf"), p1)
ggsave(filename = paste0(outdir, "ESTIMATE_boxplots.png"), p1)

# CARD Boxplots (main TCs)
#boxplot 1
CARDmain1 <- mainmetadata %>% 
  select(cell, TC, B_cell, Macrophage, Smooth_muscle_cells) %>%
  pivot_longer(cols = cellgroup1, names_to = "CARDmain1", 
               values_to = "Composition")

pcardmain1 <- ggboxplot(CARDmain1, x = "TC", y = "Composition", fill = "TC") + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = labels) +
  theme(legend.position = "bottom") +
  facet_wrap(~CARDmain1)


# boxplot 2
CARDmain2 <- mainmetadata %>% 
  select(cell, TC, Normal_epithelial_cells:Smooth_muscle_cells) %>%
  pivot_longer(cols = cellgroup2, names_to = "CARDmain2", 
               values_to = "Composition")

pcardmain2 <- ggboxplot(CARDmain2, x = "TC", y = "Composition", fill = "TC") + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = labels) +
  theme(legend.position = "bottom") +
  facet_wrap(~CARDmain2)

# boxplot 3
CARDmain3 <- mainmetadata %>% 
  select(cell, TC, StromalScore:ESTIMATEScore) %>%
  pivot_longer(cols = StromalScore:ESTIMATEScore, 
               names_to = "CARDmain3", 
               values_to = "Score")

pcardmain3 <- ggboxplot(CARDmain3, x = "TC", y = "Score", fill = "TC") + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = labels) +
  theme(legend.position = "bottom") +
  facet_wrap(~CARDmain3)

ggsave(filename = paste0(outdir, "cellcomp1_boxplots.pdf"), pcardmain1)
ggsave(filename = paste0(outdir, "cellcomp2_boxplots.png"), pcardmain2)



# BCS Boxplots
afa_sub <- ggplot(metadata %>% drop_na(bc_clusters_res.0.2),
              aes(x = bc_clusters_res.0.2, 
                  y = sig_21407,
                  fill = bc_clusters_res.0.2)) +
  geom_boxplot() + theme(legend.position="none") +
  stat_compare_means(aes(group = bc_clusters_res.0.2), 
                     label = "p.format",
                     method = "wilcox.test",
                     #method.args = list(var.equal = FALSE,
                     #                    alternative = "less"),
                     label.x = 1.4)

gef_sub <- ggplot(metadata %>% drop_na(bc_clusters_res.0.2),
                  aes(x = bc_clusters_res.0.2, 
                      y = sig_21195,
                      fill = bc_clusters_res.0.2)) +
  geom_boxplot() + theme(legend.position="none") +
  stat_compare_means(aes(group = bc_clusters_res.0.2), 
                     label = "p.format",
                     method = "wilcox.test",
                     #method.args = list(var.equal = FALSE,
                     #                    alternative = "less"),
                     label.x = 1.4)

wrap_plots(afa_sub, gef_sub)

# --- Data Tests---
df <- read.table("outputs/bc_TC0_metadata.csv", header = TRUE, sep = ",")
dfmain <- read.table("outputs/bc_metadata.csv", header = TRUE, sep = ",")
# --- Code ---
# Create outdir
outdir <- "outputs/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Rename column
df <- df %>% rename(cell = X, subTC = bc_clusters_res.0.2)
dfmain <- dfmain %>% rename(cell = X, TC = bc_clusters_res.0.2)
# Comparisons
comp <- list(c("0", "1"))

# ESTIMATE plots
estimate <- df %>% 
  select(cell, subTC, StromalScore:ESTIMATEScore) %>%
  pivot_longer(cols = StromalScore:ESTIMATEScore, names_to = "ESTIMATE", 
               values_to = "Score")

p1 <- ggboxplot(estimate, x = "subTC", y = "Score", fill = "subTC") + 
  stat_compare_means(label = "p.signif", comparisons = comp) +
  theme(legend.position = "bottom") +
  facet_wrap(~ESTIMATE)
p1
ggsave(filename = paste0(outdir, "ESTIMATE_boxplots.pdf"), p1)
ggsave(filename = paste0(outdir, "ESTIMATE_boxplots.png"), p1)

# Beyondcell plots 
beyondcell <- df %>% 
  select(cell, subTC, sig_21195, sig_21407) %>%
  pivot_longer(cols = sig_21195:sig_21407, names_to = "signature", 
               values_to = "BCS") %>%
  mutate(signature = case_when(signature == "sig_21195" ~ "sig_21195: Gefitinib",
                               signature == "sig_21407" ~ "sig_21407: Afatinib"))

p2 <- ggboxplot(beyondcell, x = "subTC", y = "BCS", fill = "subTC") + 
  stat_compare_means(label = "p.signif", comparisons = comp) +
  ylab("Normalized BCS") +
  theme(legend.position = "bottom") +
  facet_wrap(~signature)
p2
ggsave(filename = paste0(outdir, "Beyondcell_boxplots.pdf"), p2)
ggsave(filename = paste0(outdir, "Beyondcell_boxplots.png"), p2)





