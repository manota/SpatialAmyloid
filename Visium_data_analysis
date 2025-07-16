##10x Visium analysis with Seurat--------------------
##Rererence 1 https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html
##Reference 2 https://satijalab.org/seurat/archive/v3.2/spatial_vignette.html
#clear the decks
rm(list=ls())

##Load following packages---------

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the dataset --------------------------------------------

APPKI_8M.data <- '/Volumes/SSD-PUTU3C/Visium loupe file/8M-APPKI-M373/outs'
WT_8M.data <- '/Volumes/SSD-PUTU3C/Visium loupe file/8M-WT-M359/outs'
APPKI_4M.data <- '/Volumes/SSD-PUTU3C/Visium loupe file/4M-APPKI-M409/outs'
WT_4M.data <- '/Volumes/SSD-PUTU3C/Visium loupe file/4M-WT-M408/outs'


# Visualize file list that should contain "filtered_feature_bcmatrix.h5" & "tissue_lowres_image.png"
list.files(APPKI_8M.data)
list.files(WT_8M.data)
list.files(APPKI_4M.data)
list.files(WT_4M.data)

# Load the datasets
APPKI_8M_M373 <- Load10X_Spatial(data.dir = APPKI_8M.data, slice = "appki8m")
WT_8M_M359 <- Load10X_Spatial(data.dir = WT_8M.data, slice = "wt8m")
APPKI_4M_M409 <- Load10X_Spatial(data.dir = APPKI_4M.data, slice = "appki4m")
WT_4M_M408 <- Load10X_Spatial(data.dir = WT_4M.data, slice = "wt4m")
APPKI_8M_M373
WT_8M_M359
APPKI_4M_M409
WT_4M_M408

APPKI_8M_M373@meta.data
WT_8M_M359@meta.data
APPKI_4M_M409@meta.data
WT_4M_M408@meta.data

APPKI_8M_M373@meta.data$strain <- "APPKI_8M"
WT_8M_M359@meta.data$strain <- "WT_8M"
APPKI_4M_M409@meta.data$strain <- "APPKI_4M"
WT_4M_M408@meta.data$strain <- "WT_4M"

APPKI_8M_M373@meta.data$strain_8m <- "APPKI_8M"
WT_8M_M359@meta.data$strain_8m <- "WT_8M"
APPKI_4M_M409@meta.data$strain_4m <- "APPKI_4M"
WT_4M_M408@meta.data$strain_4m <- "WT_4M"

# Merge into one seurat object--------------------------------------------

Brain_8M <- merge(WT_8M_M359, APPKI_8M_M373)
Brain_8M

Brain_8M@meta.data

Brain_4M <- merge(WT_4M_M408, APPKI_4M_M409)
Brain_4M

Brain_4M@meta.data

Brain_4M_8M <- merge(Brain_4M, Brain_8M)
Brain_4M_8M

# Quality control-------------------------------------------------------

Brain_4M_8M <- PercentageFeatureSet(Brain_4M_8M, "^mt-", col.name = "percent_mito")
Brain_4M_8M <- PercentageFeatureSet(Brain_4M_8M, "^Hb.*-", col.name = "percent_hb")
plot1 <- VlnPlot(Brain_4M_8M, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_hb"), pt.size = 0.1, ncol = 2, group.by = "strain") + NoLegend()
plot2 <- SpatialFeaturePlot(Brain_4M_8M, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito","percent_hb"))
plot1
plot2


# Filter low quality spots and contaminated red blood cells (hemoglobin genes)

Brain_4M_8M.filtered <- Brain_4M_8M[, Brain_4M_8M$nFeature_Spatial > 200 & Brain_4M_8M$percent_mito < 35 & Brain_4M_8M$percent_hb < 20]
plot3 <- SpatialFeaturePlot(Brain_4M_8M.filtered, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"))
plot3


# Analysis (SCTransform for Normalization)

Brain_4M_8M.filtered.sct <- SCTransform(Brain_4M_8M.filtered, assay = "Spatial", verbose = TRUE, method = "poisson")
# #240209 Set the variable.features.n = 30000 to obtain the scale.data
# Brain_4M_8M.filtered.sct <- SCTransform(Brain_4M_8M.filtered, assay = "Spatial", verbose = TRUE, method = "poisson",variable.features.n = 30000)

# Dimensionality reduction and clustering
Brain_4M_8M.umap <- RunPCA(Brain_4M_8M.filtered.sct, assay = "SCT", verbose = FALSE)
Brain_4M_8M.umap <- FindNeighbors(Brain_4M_8M.umap, reduction = "pca", dims = 1:30)
Brain_4M_8M.umap <- FindClusters(Brain_4M_8M.umap, verbose = FALSE)
Brain_4M_8M.umap <- RunUMAP(Brain_4M_8M.umap, reduction = "pca", dims = 1:30)

# Umap plot
plot4 <- DimPlot(Brain_4M_8M.umap, reduction = "umap", group.by = "strain", order = c("WT_4M", "APPKI_4M", "WT_8M", "APPKI_8M"))
plot5 <- DimPlot(Brain_4M_8M.umap, reduction = "umap", label = TRUE,label.size = 8, repel = TRUE)
plot6 <- DimPlot(Brain_4M_8M.umap, reduction = "umap", label = FALSE)
plot4 + ggtitle("UMAP Plot") + labs(x = "UMAP 1", y = "UMAP 2") +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5), # タイトルを中央揃え＆文字サイズ
    axis.title.x = element_text(size = 18),           # x軸ラベルの文字サイズ
    axis.title.y = element_text(size = 18),           # y軸ラベルの文字サイズ
    axis.text.x = element_text(size = 18),            # x軸目盛りラベルの文字サイズ
    axis.text.y = element_text(size = 18),             # y軸目盛りラベルの文字サイズ
    legend.text = element_text(size = 14)           # 凡例ラベルの文字サイズを設定
  )
plot5 + ggtitle("UMAP Plot") + labs(x = "UMAP 1", y = "UMAP 2") +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5), # タイトルを中央揃え＆文字サイズ
    axis.title.x = element_text(size = 18),           # x軸ラベルの文字サイズ
    axis.title.y = element_text(size = 18),           # y軸ラベルの文字サイズ
    axis.text.x = element_text(size = 18),            # x軸目盛りラベルの文字サイズ
    axis.text.y = element_text(size = 18),             # y軸目盛りラベルの文字サイズ
    legend.text = element_text(size = 14)           # 凡例ラベルの文字サイズを設定
  )
plot6 + ggtitle("UMAP Plot") + labs(x = "UMAP 1", y = "UMAP 2") +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5), # タイトルを中央揃え＆文字サイズ
    axis.title.x = element_text(size = 18),           # x軸ラベルの文字サイズ
    axis.title.y = element_text(size = 18),           # y軸ラベルの文字サイズ
    axis.text.x = element_text(size = 18),            # x軸目盛りラベルの文字サイズ
    axis.text.y = element_text(size = 18),             # y軸目盛りラベルの文字サイズ
    legend.text = element_text(size = 14)           # 凡例ラベルの文字サイズを設定
  )


plot7 <- DimPlot(Brain_4M_8M.umap, reduction = "umap", split.by = "strain", ncol = 2, label = TRUE, label.size = 6, repel = TRUE)
plot8 <- DimPlot(Brain_4M_8M.umap, reduction = "umap", split.by = "strain", ncol = 2, label = FALSE)
plot7 + ggtitle("UMAP Plot") + labs(x = "UMAP 1", y = "UMAP 2") +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5), # タイトルを中央揃え＆文字サイズ
    axis.title.x = element_text(size = 18),           # x軸ラベルの文字サイズ
    axis.title.y = element_text(size = 18),           # y軸ラベルの文字サイズ
    axis.text.x = element_text(size = 18),            # x軸目盛りラベルの文字サイズ
    axis.text.y = element_text(size = 18),             # y軸目盛りラベルの文字サイズ
    legend.text = element_text(size = 14)           # 凡例ラベルの文字サイズを設定
  )

plot8 + ggtitle("UMAP Plot") + labs(x = "UMAP 1", y = "UMAP 2") +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5), # タイトルを中央揃え＆文字サイズ
    axis.title.x = element_text(size = 18),           # x軸ラベルの文字サイズ
    axis.title.y = element_text(size = 18),           # y軸ラベルの文字サイズ
    axis.text.x = element_text(size = 18),            # x軸目盛りラベルの文字サイズ
    axis.text.y = element_text(size = 18),             # y軸目盛りラベルの文字サイズ
    legend.text = element_text(size = 14)           # 凡例ラベルの文字サイズを設定
  )

# Spatial Dimplot
plot9 <- SpatialDimPlot(Brain_4M_8M.umap, label = TRUE, label.box = FALSE, label.size = 5, ncol = 2)
plot9 &
  theme(
    legend.position = "bottom",                 # 凡例を下部に移動
    legend.text = element_text(size = 14),      # 凡例ラベルの文字サイズ
    legend.title = element_text(size = 12)      # 凡例タイトルの文字サイズ
  )

plot9_2 <- SpatialDimPlot(Brain_4M_8M.umap, label = FALSE, ncol = 2)
plot9_2 &
  theme(
    legend.position = "bottom",                 # 凡例を下部に移動
    legend.text = element_text(size = 14),      # 凡例ラベルの文字サイズ
    legend.title = element_text(size = 12)      # 凡例タイトルの文字サイズ
  )

plot10 <- SpatialDimPlot(Brain_4M_8M.umap, cells.highlight = CellsByIdentities(Brain_4M_8M.umap, idents = c(19)), facet.highlight = TRUE, images = "appki4m")
plot11 <- SpatialDimPlot(Brain_4M_8M.umap, cells.highlight = CellsByIdentities(Brain_4M_8M.umap, idents = c(19)), facet.highlight = TRUE, images = "appki8m")
plot10 + plot11

###色を変える時--------------------------------------------------------------
highlight_colors <- c("white", "black")  #ハイライト色と背景色を定義

plot10 <- SpatialDimPlot(
  Brain_4M_8M.umap, 
  cells.highlight = CellsByIdentities(Brain_4M_8M.umap, idents = c(19)), 
  facet.highlight = TRUE, 
  images = "appki4m", 
  cols.highlight = highlight_colors, alpha = 0.4
)

plot11 <- SpatialDimPlot(
  Brain_4M_8M.umap, 
  cells.highlight = CellsByIdentities(Brain_4M_8M.umap, idents = c(19)), 
  facet.highlight = TRUE, 
  images = "appki8m", 
  cols.highlight = highlight_colors, alpha = 0.4
)

# プロットを並べて表示
plot10 + plot11

#Visualize gene expression of individual genes (For example, Hippocampal marker Hpca, Chroid plexus marker Ttr)
plot12 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Lrp1", "Atox1", "Rpl41", "Nrsn1", "Ehd4"))
plot13 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Tgfbr1", "Mtch1", "Tmem176a", "Rplp1", "Rps19"))
plot14 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Pabpc1", "Rac2", "Asah1", "Hsp90aa1", "Mef2c"))
plot15 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Rps26", "Cx3cr1", "Lrrc4", "Spp1", "Mafb"))
plot16 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Rpl12", "Hspa5", "Syt11", "Gstp1", "Rps15a"))
plot17 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Tmem47", "Itm2b", "2410006H16Rik", "Manf", "Ubc"))
plot18 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Skap2", "St3gal6", "Rpl37a", "Lhfpl2", "Rpl39"))
plot19 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Cmtm3", "Lpcat2", "Rpl30", "Ssr4", "Rps28"))

plot12
plot13
plot14
plot15
plot16
plot17
plot18
plot19

plot20 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Tgfbr1", "Pabpc1", "Asah1", "Rps26", "Cx3cr1"))
plot21 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Spp1", "Siglech", "Hif1a", "B2m", "Rps25"))
plot22 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Slamf9", "Pmepa1", "Slc11a1", "Baiap2l2", "Abcd2"))
plot23 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Ctsb", "Arhgdib", "Axl", "Ifi27l2a", "Il10ra"))
plot24 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Itm2b", "2410006H16Rik", "Ubc", "St3gal6", "Rpl37a"))
plot25 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Cmtm3", "Lpcat2", "Rpl30", "Ssr4", "Rps28"))
plot26 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Csf1", "H2-T23", "Rps21", "Pld4", "Lamp2"))
plot27 <- SpatialFeaturePlot(Brain_4M_8M.umap, features = c("Ank", "Gns", "Lair1", "Lyn", "Lgmn"))

plot20
plot21
plot22
plot23
plot24
plot25
plot26
plot27


# Find all markers among all clusters
All_markers <- FindAllMarkers(Brain_4M_8M.umap)

# Save the gene list
write.table(All_markers,"/Users/okiru/Desktop/230601 APPKI_Visium_Seurat/Spatial_all_markers.txt", sep = "\t")




