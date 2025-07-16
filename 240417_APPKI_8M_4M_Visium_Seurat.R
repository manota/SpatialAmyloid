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
library(stringr)

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

#Visualize gene expression of individual genes (For example, Hippocampal marker Hpca, Chroid plexus marker Ttr)
SpatialFeaturePlot(Brain_4M_8M.filtered.sct, features = c("Cd74", "Itgax", "Tmem119", "Hpca", "Ttr", "Mbp", "Mobp"))
SpatialFeaturePlot(Brain_4M_8M.filtered.sct, features = c("Mbp", "Mog", "Mag", "Cnp", "Plp1", "Mobp", "Olig2", "Olig1"))

# Dimensionality reduction and clustering
Brain_4M_8M.umap <- RunPCA(Brain_4M_8M.filtered.sct, assay = "SCT", verbose = FALSE)
Brain_4M_8M.umap <- FindNeighbors(Brain_4M_8M.umap, reduction = "pca", dims = 1:30)
Brain_4M_8M.umap <- FindClusters(Brain_4M_8M.umap, verbose = FALSE)
Brain_4M_8M.umap <- RunUMAP(Brain_4M_8M.umap, reduction = "pca", dims = 1:30)

#Save cluster numbers
cluster_barcodes <- Brain_4M_8M.umap@active.ident
cluster_strains <-  Brain_4M_8M.umap@meta.data$strain
cluster_barcodes <-as.data.frame(cluster_barcodes)
cluster_strains <-as.data.frame(cluster_strains)
cluster_barcodes_strains <-cbind(cluster_barcodes, cluster_strains)
list_of_cluster_barcodes <- split(cluster_barcodes_strains, cluster_barcodes_strains$cluster_barcodes)

# list_of_data[[1]]

save_parent_path<-'/Users/okiru/Desktop/240124 Visium_Hiroshima University/saved_cluster_barcodes_strain'
sample_sizes<-c()
for(i in 1:length(list_of_cluster_barcodes)){
  barcodes<-rownames(list_of_cluster_barcodes[[i]])
  strains<-list_of_cluster_barcodes[[i]]$cluster_strains
  barcodes_strains<-cbind(as.data.frame(barcodes), as.data.frame(strains))
  print(head(barcodes_strains))
  barcodes_strains$barcodes<-apply(as.data.frame(barcodes_strains$barcodes), 1, function(x) {
    x <- stringr::str_replace(x, "-1_1_", "-tmp-1_")
    x <- stringr::str_replace(x, "-1_2_", "-1_1_")
    x <- stringr::str_replace(x, "-tmp-1_", "-1_2_")
    return(x)
  }) ###1-2をWT、1-1をAPPに対応させるため実行
  print(head(barcodes_strains))
  sample_sizes<-c(sample_sizes,nrow(barcodes_strains))
  cluster_num_string<-as.character(unique(list_of_cluster_barcodes[[i]]$cluster_barcodes)[1])
  #cluster_num_string<-as.character(unique(list_of_cluster_barcodes[[i]]))
  print(cluster_num_string)
  save_name<-file.path(save_parent_path,paste(cluster_num_string,'.csv',sep=''))
  write.csv(barcodes_strains,save_name)
}
print(sum(sample_sizes))

