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
library(jsonlite)
library(imager)


filtering_spots<-function(ob_10X){
  ob_10X <- PercentageFeatureSet(ob_10X, "^mt-", col.name = "percent_mito")
  ob_10X <- PercentageFeatureSet(ob_10X, "^Hb.*-", col.name = "percent_hb")
  ob_10X <- ob_10X [,ob_10X$nFeature_Spatial > 200 & ob_10X$percent_mito < 35 &ob_10X$percent_hb < 20]
  return(ob_10X)
}

SCTtransform_obtain_all_genes<-function(ob_10X){
  ob_10X <- SCTransform(ob_10X, assay = "Spatial", verbose = TRUE, method = "poisson",variable.features.n = 30000)
  gene_nums<-nrow(ob_10X@assays$SCT@data)
  ob_10X<- SCTransform(ob_10X, assay = "Spatial", verbose = TRUE, method = "poisson",variable.features.n = gene_nums)
  return(ob_10X)
}

get_corrected_coord_radius<-function(dir_name,sample_coord){
  json_path<-paste(dir_name,'spatial/scalefactors_json.json',sep='/')
  Json_file <- read_json(json_path)
  # change scale
  scale_factor<- read_json(json_path)
  corrected_y_x<-sample_coord[, 4:5]*scale_factor$tissue_lowres_scalef
  #row is y, col is x
  spot_radius<-scale_factor$spot_diameter_fullres*scale_factor$tissue_lowres_scalef/2
  return(list(corrected_y_x,spot_radius))
}

get_red_intensitiy<-function(dir_name,corrected_y_x){
  image_path<-paste(dir_name,'spatial/tissue_lowres_image.png',sep='/')
  ori_img<-load.image(image_path)
  redlayer_values <- as.vector(ori_img[,,1])
  int_corrected_y_x<-corrected_y_x
  int_corrected_y_x$imagecol <- lapply(corrected_y_x$imagecol, as.integer)
  int_corrected_y_x$imagerow <- lapply(corrected_y_x$imagerow, as.integer)
  
  red_intensities<-matrix(nrow=nrow(int_corrected_y_x),ncol=1)
  j=1
  for (spot_barcode in row.names(int_corrected_y_x)){
    print(j)
    print(spot_barcode)
    tmp_x <- as.numeric(int_corrected_y_x[row.names(int_corrected_y_x)==spot_barcode,]$imagecol)
    tmp_y <- as.numeric(int_corrected_y_x[row.names(int_corrected_y_x)==spot_barcode,]$imagerow)
    tmp_red_intensity <- ori_img[tmp_x, tmp_y, 1]
    red_intensities[j,1]<-tmp_red_intensity
    j=j+1
  }
  red_intensities<-as.data.frame(red_intensities)
  rownames(red_intensities)<-rownames(int_corrected_y_x)
  return(red_intensities)
}


# Load the dataset --------------------------------------------
parent_dir<-'/Volumes/SSD-PUTU3C/Visium loupe file'
save_dir<-'/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data'
dir_names<-list.dirs(parent_dir,recursive = F)
dir_names<-dir_names[-length(dir_names)]  #remove last file 'Gazo file" (-1)
sample_names<-dir_names
sample_names<-basename(dir_names)    #obtain only file name
dir_names<-paste(dir_names,'outs',sep='/')  #connect the dir_names and '/' and 'outs'

pair_dir_names<-list(c(dir_names[1],dir_names[2]),c(dir_names[3],dir_names[4]))
pair_save_names<-c('4month','8month')
 #get scaled_data
i=1
for(pair_dir_name in pair_dir_names){
  pair_dir_name<-as.vector(pair_dir_names[i])
  pair_save_name<-pair_save_names[i]
  app_dir<-pair_dir_name[[1]][1]
  wt_dir<-pair_dir_name[[1]][2]
  app<-Load10X_Spatial(data.dir = app_dir,slice='app')
  wt<-Load10X_Spatial(data.dir = wt_dir,slice='wt')

  app@meta.data$strain <- "app"
  wt@meta.data$strain <- "wt"
  app_wt<-merge(app,wt)
  app_wt<-filtering_spots(app_wt)
  app_wt.sct <- SCTtransform_obtain_all_genes(app_wt)
  scaled_data<-as.data.frame(app_wt.sct@assays$SCT@scale.data)
  
  #forced_changed_barcode_names<-colnames(scaled_data)
  
  app_raw_coord<-app_wt.sct@images$app@coordinates
  wt_raw_coord<-app_wt.sct@images$wt@coordinates
  
  app_coord_radius<-get_corrected_coord_radius(app_dir,app_raw_coord)
  wt_coord_radius<-get_corrected_coord_radius(wt_dir,wt_raw_coord)
  
  #get_corrected_cord_radius is list(corrected_y_x,spot_radius)
  
  app_corrected_y_x<-as.data.frame(app_coord_radius[[1]])
  app_corrected_radius<-as.matrix(app_coord_radius[[2]])
  wt_corrected_y_x<-as.data.frame(wt_coord_radius[[1]])
  wt_corrected_radius<-as.matrix(wt_coord_radius[[2]])
  
  #app_red_intensities<-get_red_intensitiy(app_dir,app_corrected_y_x)
  #wt_red_intensities<-get_red_intensitiy(wt_dir,wt_corrected_y_x)
  #red_intensities<-rbind(app_red_intensities,wt_red_intensities)
  
  image_path<-paste(app_dir,'spatial/tissue_lowres_image.png',sep='/')
  ori_img<-load.image(image_path)
  redlayer_values <- as.vector(ori_img[,,1])
  corrected_y_x_t<-as.data.frame(app_corrected_y_x)
  int_corrected_y_x_t<-corrected_y_x_t
  int_corrected_y_x_t$imagecol <- lapply(corrected_y_x_t$imagecol, as.integer)
  int_corrected_y_x_t$imagerow <- lapply(corrected_y_x_t$imagerow, as.integer)
  
  red_intensities<-matrix(nrow=nrow(int_corrected_y_x_t),ncol=1)
  j=1
  for (spot_barcode in row.names(int_corrected_y_x_t)){
    print(j)
    print(spot_barcode)
    tmp_x <- as.numeric(int_corrected_y_x_t[row.names(int_corrected_y_x_t)==spot_barcode,]$imagecol)
    tmp_y <- as.numeric(int_corrected_y_x_t[row.names(int_corrected_y_x_t)==spot_barcode,]$imagerow)
    tmp_red_intensity <- ori_img[tmp_x, tmp_y, 1]
    red_intensities[j,1]<-tmp_red_intensity
    j=j+1
  }
  red_intensities<-as.data.frame(red_intensities)
  rownames(red_intensities)<-rownames(int_corrected_y_x_t)
  app_red_intensities<-red_intensities
  
  
  image_path<-paste(wt_dir,'spatial/tissue_lowres_image.png',sep='/')
  ori_img<-load.image(image_path)
  redlayer_values <- as.vector(ori_img[,,1])
  corrected_y_x_t<-as.data.frame(wt_corrected_y_x)
  int_corrected_y_x_t<-corrected_y_x_t
  int_corrected_y_x_t$imagecol <- lapply(corrected_y_x_t$imagecol, as.integer)
  int_corrected_y_x_t$imagerow <- lapply(corrected_y_x_t$imagerow, as.integer)
  
  red_intensities<-matrix(nrow=nrow(int_corrected_y_x_t),ncol=1)
  j=1
  for (spot_barcode in row.names(int_corrected_y_x_t)){
    print(j)
    print(spot_barcode)
    tmp_x <- as.numeric(int_corrected_y_x_t[row.names(int_corrected_y_x_t)==spot_barcode,]$imagecol)
    tmp_y <- as.numeric(int_corrected_y_x_t[row.names(int_corrected_y_x_t)==spot_barcode,]$imagerow)
    tmp_red_intensity <- ori_img[tmp_x, tmp_y, 1]
    red_intensities[j,1]<-tmp_red_intensity
    j=j+1
  }
  red_intensities<-as.data.frame(red_intensities)
  rownames(red_intensities)<-rownames(int_corrected_y_x_t)
  wt_red_intensities<-red_intensities
  red_intensities<-rbind(app_red_intensities,wt_red_intensities)
  
  scaled_data
  red_intensities
  app_corrected_y_x
  app_corrected_radius
  wt_corrected_y_x
  wt_corrected_radius
  
  save_base_name<-paste(pair_save_name,'scaledata',sep='_')
  save_base_name<-paste(save_base_name,'csv',sep='.')
  save_file_name<-paste(save_dir,save_base_name,sep='/')
  write.csv(t(scaled_data),save_file_name) ####240423 転置して保存しなおした。
  
  save_base_name<-paste(pair_save_name,'red',sep='_')
  save_base_name<-paste(save_base_name,'csv',sep='.')
  save_file_name<-paste(save_dir,save_base_name,sep='/')
  write.csv(red_intensities,save_file_name)
  
  save_base_name<-paste(pair_save_name,'app_coord',sep='_')
  save_base_name<-paste(save_base_name,'csv',sep='.')
  save_file_name<-paste(save_dir,save_base_name,sep='/')
  write.csv(app_corrected_y_x,save_file_name)
  
  save_base_name<-paste(pair_save_name,'app_radius',sep='_')
  save_base_name<-paste(save_base_name,'csv',sep='.')
  save_file_name<-paste(save_dir,save_base_name,sep='/')
  write.csv(app_corrected_radius,save_file_name)
  
  save_base_name<-paste(pair_save_name,'wt_coord',sep='_')
  save_base_name<-paste(save_base_name,'csv',sep='.')
  save_file_name<-paste(save_dir,save_base_name,sep='/')
  write.csv(wt_corrected_y_x,save_file_name)
  
  save_base_name<-paste(pair_save_name,'wt_radius',sep='_')
  save_base_name<-paste(save_base_name,'csv',sep='.')
  save_file_name<-paste(save_dir,save_base_name,sep='/')
  write.csv(wt_corrected_radius,save_file_name)
  i=i+1
}
