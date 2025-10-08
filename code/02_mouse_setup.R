library(data.table)
library(dplyr)
library(Seurat)
library(stringr)
library(BuenColors)
library(harmony)
source("00_gene_sets.R")

# Import
exp2d5_hdf <- fread("../data/hash/mHSPC_2_d5_hash_ids_Seurat.csv", header = FALSE) %>% filter(!(V2 %in% c("Doublet", "Negative")))
exp2d8_hdf <- fread("../data/hash/mHSPC_2_d8_hash_ids_Seurat.csv", header = FALSE) %>% filter(!(V2 %in% c("Doublet", "Negative")))
exp2d8_hdf$V1 <- gsub("-1", "-2", exp2d8_hdf$V1)

# collate meta
full_meta <- rbind(exp2d5_hdf, exp2d8_hdf)
list_mat_d5 <- Read10X_h5("../data/counts/mHSPC_2_d5_filtered_feature_bc_matrix_corrected.h5")
list_mat_d8 <- Read10X_h5("../data/counts/mHSPC_2_d8_filtered_feature_bc_matrix_corrected.h5")
d8mat <- list_mat_d8[["Gene Expression"]]; colnames(d8mat) <- gsub("-1", "-2", colnames(d8mat))
RNA_mat <- cbind(list_mat_d5[["Gene Expression"]],
                 d8mat
)[,full_meta$V1]

# Do basic seurat handling
so <- CreateSeuratObject(counts = RNA_mat)
so$full_geno <- full_meta$V2
so$mouse_genotype <-  gsub("WT", "aWT", str_split_fixed( full_meta$V2, "_", 3)[,1])
so$gRNA <- gsub("Rosa26", "gWT", str_split_fixed( full_meta$V2, "_", 3)[,2])
so$day <- str_split_fixed( full_meta$V2, "_", 3)[,3]
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
summary(so$percent.mt)
summary(so$nCount_RNA)
summary(so$nFeature_RNA)
so <- subset(so, percent.mt < 5 & nCount_RNA > 2000)

# Do dimensionality reduction
s.genes <- human_to_mouse_gs(cc.genes$s.genes)
g2m.genes <- human_to_mouse_gs(cc.genes$g2m.genes)

so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes)

so <- so %>% RunHarmony(., c("Phase")) 
so <- so %>% FindNeighbors(reduction = "harmony", dims = 1:30)
so <- so %>% RunUMAP(dims = 1:30, reduction  = "harmony") 
so <- so %>% FindClusters(resolution = 0.4)
so$cl_clusters <- case_when(
  so$seurat_clusters %in% c(7,9,12) ~ "HSPC",
  so$seurat_clusters %in% c(0,13) ~ "Stem_like",
  so$seurat_clusters %in% c(1,8) ~ "Neutrophil_like",
  so$seurat_clusters %in% c(2,3,5) ~ "Macrophage_like",
  so$seurat_clusters %in% c(6,10) ~ "Megakaryocyte_like",
  so$seurat_clusters %in% c(4,11) ~ "Baso_Eso_like",
  TRUE ~ "HSPC"
)


###--- here
saveRDS(so, file = "../output/big_seurat_object.rds")

so <- readRDS("../output/big_seurat_object.rds")

###--- here

if(FALSE){
  
  DimPlot(so, group.by = c("full_geno", "seurat_clusters", "gRNA", "day", "Phase"), shuffle = TRUE)
  pumapb <- DimPlot(so, group.by = "cl_clusters", label = FALSE, shuffle = TRUE) + 
    scale_color_manual(values = c("firebrick","lightgrey", "dodgerblue3","purple2", "orange1","#7f7f7f"))+
    theme_void() + ggtitle("") + theme(legend.position = "none")
  
  so$geno2 <- paste0(so$mouse_genotype, "_", so$gRNA)
  df <- data.frame(so@meta.data, so@reductions$umap@cell.embeddings)
  big <- ggplot(df, aes(x = umap_1, y = umap_2, color = full_geno)) +
    geom_point(size = 0.3) +
    scale_color_manual(values = jdb_palette("corona")) +
    facet_grid(day~geno2) + pretty_plot(fontsize = 8)
  cowplot::ggsave2(big, file = "../output/big_out.pdf", width = 7, height = 2.6)
  
  cowplot::ggsave2(pumapb, file = "../output/umap_base.pdf", width = 3.6, height = 3.0)
  cowplot::ggsave2(pumapb, file = "../output/umap_base.png", width = 3.0, height = 3.0, dpi = 400)
  
  mk_plot <- function(gene){
    pu <- FeaturePlot(so, features = c( gene),  
                      pt.size = 0.1, max.cutoff = "q98") + FontSize(main = 0.0001) + 
      theme_void() + theme(legend.position = "none") + ggtitle("") + 
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      scale_color_gradientn(colors = c("lightgrey", jdb_palette("solar_rojos")[c(2:9)]))
    return(pu)
  }
  
  # FeaturePlot(so, features = c("Cd34", "Mecom", "Adgre1", "Elane", "Itga2b",  "Gata2"))
  
  cowplot::ggsave2(
    cowplot::plot_grid(
      mk_plot("Cd34"),
      mk_plot("Mecom"),
      mk_plot("Adgre1"),
      mk_plot("Elane"),
      mk_plot("Itga2b"),
      mk_plot("Gata2"), 
      ncol = 3, scale = 1
    ), file = "../output/markers_umap.png",width = 3*3, height = 1.5*3, dpi = 600)
  
  
  px <- so@meta.data %>% group_by(gRNA, cl_clusters, mouse_genotype,day) %>% # day
    summarize(count = n()) %>%
    group_by(gRNA, mouse_genotype,day) %>% # day
    mutate(pct = count / sum(count)*100) %>%
    ggplot( aes(x = interaction( mouse_genotype, gRNA,day), y=count*100, fill = cl_clusters)) +
    geom_bar(position="fill", stat="identity", color = "black", width = 0.7) +
    scale_fill_manual(values = c("firebrick","lightgrey", "dodgerblue3","purple2", "orange1","#7f7f7f")) +
    pretty_plot(fontsize = 7) + L_border() +
    theme(legend.position = "none") + scale_y_continuous(expand = c(0,0)) +
    labs(x = "", y = "% of cells") #+ facet_wrap(~day)
  px
  cowplot::ggsave2(px, file = "../output/bar_stack_splits.pdf", width = 4, height = 2)
  
  
  # make 
  px <- ggplot(so@meta.data %>% group_by(gRNA, cl_clusters) %>%
                 summarize(count = n()), aes(x = gRNA, y=count, fill = cl_clusters)) + 
    geom_bar(position="fill", stat="identity", color = "black", width = 0.7) +
    scale_fill_manual(values = c("firebrick","lightgrey", "dodgerblue3","purple2", "orange1","#7f7f7f")) +
    pretty_plot(fontsize = 7) + L_border() +
    theme(legend.position = "none") + scale_y_continuous(expand = c(0,0)) +
    labs(x = "", y = "% of cells")
  cowplot::ggsave2(px, file = "../output/bar_stack.pdf", width = 1, height = 1.5)
}
###########

ggplot(so@meta.data, aes(x = mouse_genotype, y = GOBP_REGULATION_OF_MYELOID_LEUKOCYTE_DIFFERENTIATION1, fill = gRNA)) + 
  geom_violin() + geom_boxplot(fill = NA, outlier.shape = NA)

ggplot(so@meta.data, aes(x = mouse_genotype, y = REACTOME_TRANSLATION2, fill = gRNA)) + 
  geom_boxplot()

# Function to compute pvalues
compute_wrs_p <- function( what, day_vec = c("d5", "d8"), genotype_vec = c("aWT", "R3", "R3C8"), df = so@meta.data){
  ss_df <- so@meta.data %>% 
    filter(mouse_genotype %in% genotype_vec & day %in% day_vec)
  
  wt <- wilcox.test(
    ss_df[ss_df$gRNA == "M41T",what],
    ss_df[ss_df$gRNA == "gWT",what]
  )
  wt$p.value
}

compute_wrs_p(what = "GOBP_REGULATION_OF_MYELOID_LEUKOCYTE_DIFFERENTIATION1")

compute_wrs_p(what = "deltaHuman", day_vec = "d5", genotype_vec = "aWT")
compute_wrs_p(what = "deltaHuman", day_vec = "d5", genotype_vec = "R3")
compute_wrs_p(what = "deltaHuman", day_vec = "d5", genotype_vec = "R3C8")
compute_wrs_p(what = "deltaHuman", day_vec = "d8", genotype_vec = "aWT")
compute_wrs_p(what = "deltaHuman", day_vec = "d8", genotype_vec = "R3")
compute_wrs_p(what = "deltaHuman", day_vec = "d8", genotype_vec = "R3C8")


compute_wrs_p(what = "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR4", day_vec = "d5", genotype_vec = "aWT")
compute_wrs_p(what = "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR4", day_vec = "d5", genotype_vec = "R3")
compute_wrs_p(what = "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR4", day_vec = "d5", genotype_vec = "R3C8")
compute_wrs_p(what = "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR4", day_vec = "d8", genotype_vec = "aWT")
compute_wrs_p(what = "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR4", day_vec = "d8", genotype_vec = "R3")
compute_wrs_p(what = "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR4", day_vec = "d8", genotype_vec = "R3C8")


# Add module scores
so <- AddModuleScore(so, l_m_landau, name = l_h_names)
so <- AddModuleScore(so, list(vexas_up_mouse), name = "vexas_up")
so <- AddModuleScore(so, list(vexas_down_mouse), name = "vexas_down")

so$deltaHuman <- so$vexas_up1 - so$vexas_down1

pG <- ggplot(so@meta.data, aes(x = mouse_genotype, y =  deltaHuman , fill = gRNA)) + 
  geom_violin(position = position_dodge(0.9)) + facet_wrap(~day) +
  geom_boxplot(aes(group=interaction(gRNA, mouse_genotype)), width = 0.5, color = "black",
               outlier.shape = NA, position = position_dodge(0.9), fill = NA) +
  scale_fill_manual(values = c("lightblue", "dodgerblue3")) + 
  pretty_plot(fontsize = 8) +  theme(legend.position = "none") + 
  labs(x = "", y = "Ganesan VEXAS score") 

cowplot::ggsave2(pG, 
                 file = "../output/ganesan_supplement_split.pdf", 
                 width = 2.5, height = 1.4)

pM <- ggplot(so@meta.data, aes(x = mouse_genotype, y = GOBP_REGULATION_OF_MYELOID_LEUKOCYTE_DIFFERENTIATION1, fill = gRNA)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(aes(group=interaction(gRNA, mouse_genotype)), width = 0.5, color = "black",
               outlier.shape = NA, position = position_dodge(0.9), fill = NA) +
  scale_fill_manual(values = c("lightblue", "dodgerblue3")) + 
  pretty_plot(fontsize = 8) +   theme(legend.position = "none") + 
  labs(x = "", y = "myeloid differentiation") + facet_wrap(~day)

cowplot::ggsave2(pM, 
                 file = "../output/myeloid_supplement_split.pdf", 
                 width = 2.5, height = 1.4)

pU <- ggplot(so@meta.data, aes(x = mouse_genotype, y = REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR4, fill = gRNA)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(aes(group=interaction(gRNA, mouse_genotype)), width = 0.5, color = "black",
               outlier.shape = NA, position = position_dodge(0.9), fill = NA) +
  scale_fill_manual(values = c("lightblue", "dodgerblue3")) + 
  pretty_plot(fontsize = 8) +   theme(legend.position = "none") + 
  labs(x = "", y = "upr") + facet_wrap(~day)

cowplot::ggsave2(pU, 
                 file = "../output/upr_supplement_split.pdf", 
                 width = 2.5, height = 1.4)

upr1 <- ggplot(so@meta.data, aes(x = mouse_genotype, y = REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR4, fill = gRNA)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(aes(group=interaction(gRNA, mouse_genotype)), width = 0.5, color = "black",
               outlier.shape = NA, position = position_dodge(0.9), fill = NA) +
  scale_fill_manual(values = c("lightblue", "dodgerblue3")) + 
  pretty_plot(fontsize = 8) + L_border() +  theme(legend.position = "none") + 
  labs(x = "", y = "UPR")

my1 <- ggplot(so@meta.data, aes(x = mouse_genotype, y = GOBP_REGULATION_OF_MYELOID_LEUKOCYTE_DIFFERENTIATION1, fill = gRNA)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(aes(group=interaction(gRNA, mouse_genotype)), width = 0.5, color = "black",
               outlier.shape = NA, position = position_dodge(0.9), fill = NA) +
  scale_fill_manual(values = c("lightblue", "dodgerblue3")) +
  pretty_plot(fontsize = 8) + L_border() +  theme(legend.position = "none") + 
  labs(x = "", y = "myeloid")

trans1 <- ggplot(so@meta.data, aes(x = mouse_genotype, y = REACTOME_TRANSLATION2, fill = gRNA)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(aes(group=interaction(gRNA, mouse_genotype)), width = 0.5, color = "black",
               outlier.shape = NA, position = position_dodge(0.9), fill = NA) +
  scale_fill_manual(values = c("lightblue", "dodgerblue3")) + 
  pretty_plot(fontsize = 8) + L_border() +  theme(legend.position = "none") + 
  labs(x = "", y = "Translation")

my1 <- ggplot(so@meta.data, aes(x = mouse_genotype, y = GOBP_REGULATION_OF_MYELOID_LEUKOCYTE_DIFFERENTIATION1, fill = gRNA)) + 
  geom_violin(position = position_dodge(0.9)) + 
  geom_boxplot(aes(group=interaction(gRNA, mouse_genotype)), width = 0.5, color = "black",
               outlier.shape = NA, position = position_dodge(0.9), fill = NA) +
  scale_fill_manual(values = c("lightblue", "dodgerblue3")) +
  pretty_plot(fontsize = 8) + L_border() +  theme(legend.position = "none") + 
  labs(x = "", y = "myeloid")


cowplot::ggsave2(cowplot::plot_grid( trans1, ncol = 1), 
                 file = "../output/violins_2.pdf", 
                 width = 1.2, height = 1.4)


FeaturePlot(so, features = c("Cd14", "Cd68", "Cd80", "Trem2", "Csf1r", "Adgre1"))


px <- so@meta.data %>% group_by(gRNA, cl_clusters, mouse_genotype, day) %>%
  summarize(count = n()) %>%
  group_by(gRNA, mouse_genotype, day) %>%
  mutate(pct = count / sum(count)*100) %>%
  filter(cl_clusters == "Macrophage_like") %>%
  ggplot( aes(x = interaction(gRNA, day), y=pct, color = mouse_genotype, group = mouse_genotype)) +
  geom_point() + geom_line() +
  scale_color_manual(values = c("violet", "purple2", "purple4")) +
  pretty_plot(fontsize = 7) + L_border() +
  labs(x = "", y = "% of macrophages") + theme(legend.position = "none")
cowplot::ggsave2(px,
                 file = "../output/guide.pdf", 
                 width = 1.2, height = 1.2)

# density plot
dfp <- data.frame(so@meta.data)
dfp$umap1 <- so@reductions$umap@cell.embeddings[,1]
dfp$umap2 <- so@reductions$umap@cell.embeddings[,2]
ggplot(dfp, aes(x = umap1, y = umap2, color = gRNA)) + 
  geom_point(size = 0.1) + 
  geom_density_2d()

pGRNA <- DimPlot(so, group.by = c("gRNA"), shuffle = TRUE) + 
  theme_void() + ggtitle("") + theme(legend.position = "none") + 
  scale_color_manual(values = c("lightblue", "dodgerblue3"))
cowplot::ggsave2(pGRNA, file = "../output/umap_grna.png", width = 3.0, height = 3.0, dpi = 400)

pDAY <- DimPlot(so, group.by = c("day"), shuffle = TRUE) + 
  theme_void() + ggtitle("") + theme(legend.position = "none") + 
  scale_color_manual(values = c("pink2", "firebrick"))
cowplot::ggsave2(pDAY, file = "../output/umap_day.png", width = 3.0, height = 3.0, dpi = 400)

pGeno <- DimPlot(so, group.by = c("mouse_genotype"), shuffle = TRUE) + 
  theme_void() + ggtitle("") + theme(legend.position = "none") + 
  scale_color_manual(values = c("violet", "purple2", "purple4"))
cowplot::ggsave2(pGeno, file = "../output/umap_geno.png", width = 3.0, height = 3.0, dpi = 400)

# annotate cells
all_genotypes <- rbind(
  fread("../data/bam/barcoded_genotypes/LSK.d5.txt", header = FALSE),
  fread("../data/bam/barcoded_genotypes/LSK.d8.txt", header = FALSE) %>%
    mutate(V2 = gsub("-1", "-2", V2))
)
wt_cells <- all_genotypes %>% filter(V1 == "T") %>% pull(V2)
mut_cells <- all_genotypes %>% filter(V1 == "C") %>% pull(V2)

so$observed_genotype <- case_when(
  colnames(so) %in% wt_cells ~ "WT", 
  colnames(so) %in% mut_cells ~ "M41T", 
  TRUE ~ "Anone"
)

pobs <- DimPlot(so, group.by = c("observed_genotype"), order = c("WT", "M41T", "Anone")) + 
  theme_void() + ggtitle("")  + theme(legend.position = "none") + 
  scale_color_manual(values = c( "grey", "dodgerblue3", "purple3"))
cowplot::ggsave2(pobs, file = "../output/umap_obs_geno.png", width = 3.0, height = 3.0, dpi = 400)

so$what2 <- paste0(so$mouse_genotype, "_", so$gRNA)
so_d8 <- subset(so, day == "d8")
DotPlot(so_d8, c("Mecom", "Hoxa9", "Hlf", "Mpo", "Elane", "Gfi1", "Prtn3", "Cxcr2", "S100a8", "S100a9", 
              "Cd14", "Cd68", "Csf1r","Itga2b", "Fcer1a", "Gata2", "Kit", "Mpl", "Ly6a"), group.by = "what2")
