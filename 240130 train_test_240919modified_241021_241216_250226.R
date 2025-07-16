library('glmnet')
library('caret')
library(imager)
library(data.table)
library(openxlsx)
library(ggplot2)
library(dplyr)
#library(stringr)

#Load scale data
# X<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_scaledata.csv',row.names = 1)
# X<-t(X)
# y<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_red.csv',row.names = 1)
# count_1_1 <- sum(grepl("1_1", rownames(X)))
# y[count_1_1+1:nrow(y),1]<-0 #change red intensity of wt into 0


#countXy<-cbind(X,y)

#threshold<-0.48
#y[y<threshold]<-0
#y[y>=threshold]<-1
#weights<-1/table(as.vector(y)$V1)[as.factor(as.vector(y)$V1)]
#print(sum(y))

#train_index<-createDataPartition(Xy$V1,p=0.8,list=FALSE)

#X_train<-X[train_index,]
#X_test<-X[-train_index,]
#y_train<-y[train_index,]
#y_test<-y[-train_index,]

#weights_train<-weights[train_index]


#####240920追加　scaledataのロードに時間かかかるため、fread_fixを使用
fread_fix<-function(X_path){
  X<-fread(X_path,check.names=FALSE,data.table=FALSE)
  rownames(X)<-X[,1]
  X<-X[,-1]
  return(X)
}

#training 8 month
#X_train<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_scaledata.csv',row.names = 1)
#lines<-readLines('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_scaledata.csv',n=10)
#part_data<-read.csv(textConnection((paste(lines,collapse = '\n'))))
X_train<-fread_fix('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_scaledata.csv')  ###faster than read.csv
#X_wt_train<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/8M-WT-M359.csv',row.names = 1)
y_train<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_red.csv',row.names = 1)
#y_wt_train<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/8M-WT-M359_red.csv',row.names = 1)
#y_wt_train[,]<-0

#X_test<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_scaledata.csv',row.names = 1)
X_test<-fread_fix('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_scaledata.csv')  ###faster than read.csv
y_test<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_red.csv',row.names = 1)
# test_coord_y_x<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_app_coord.csv',row.names = 1)
test_coord_y_x<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_wt_coord.csv',row.names = 1)
#X_test<-fread_fix('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_scaledata.csv')  ###faster than read.csv
#y_test<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_red.csv',row.names = 1)
#test_coord_y_x<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_app_coord.csv',row.names = 1)
#test_coord_y_x<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_wt_coord.csv',row.names = 1)

train_coord_y_x<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_app_coord.csv',row.names = 1)
#test_radius<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_app_radius.csv',row.names = 1)
test_radius<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_wt_radius.csv',row.names = 1)
#test_radius<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_app_radius.csv',row.names = 1)
# test_radius<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_wt_radius.csv',row.names = 1)


train_radius<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_app_radius.csv',row.names = 1)
#img <- load.image("/Volumes/SSD-PUTU3C/Visium loupe file/4M-APPKI-M409/outs/spatial/tissue_lowres_image.png")
# img <- load.image("/Volumes/SSD-PUTU3C/Visium loupe file/4M-WT-M408/outs/spatial/tissue_lowres_image.png")
#img <- load.image("/Volumes/SSD-PUTU3C/Visium loupe file/8M-APPKI-M373/outs/spatial/tissue_lowres_image.png")#
#img <- load.image("/Volumes/SSD-PUTU3C/Visium loupe file/8M-WT-M359/outs/spatial/tissue_lowres_image.png")

#cluster19の変動遺伝子リストをアップロード
c19_genes <- read.xlsx("/Users/okiru/Desktop/230601 APPKI_Visium_Seurat/Spatial_all_markers_cluster19 q<0.05.xlsx")
#cluster19のgene listを指定する。
#c19_gene_names<-c19_genes$gene[c19_genes$FC>1.5|c19_genes$FC<1/1.5]
c19_gene_names<-c19_genes$gene
# X_train<-t(X_train)    ###240423 scaledataを転置して保存し直したため不要
# #X_wt_train<-t(X_wt_train)
# X_test<-t(X_test) ###240423 scaledataを転置して保存し直したため不要
#y_train<-rbind(y_train,y_wt_train)
count_1_1 <- sum(grepl("1_1", rownames(X_train)))
y_train[(count_1_1+1):nrow(y_train),1]<-0 #change red intensity of wt into 0 


#common_gene_names<-intersect(intersect(colnames(X_train),colnames(X_wt_train)),colnames(X_test))
common_gene_names<- intersect(colnames(X_train),colnames(X_test))
X_train<-X_train[,common_gene_names]
#X_wt_train<-X_wt_train[,common_gene_names]
X_test<-X_test[,common_gene_names]
#X_train<-rbind(X_train,X_wt_train)

X_train<-X_train[,c19_gene_names]
#X_wt_train<-X_wt_train[,common_gene_names]
X_test<-X_test[,c19_gene_names]


#traing 4 month
#X_test<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/8M-APPKI-M373.csv',row.names = 1)
#y_test<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/8M-APPKI-M373_red.csv',row.names = 1)
#X_train<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/4M-APPKI-M409.csv',row.names = 1)
#y_train<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/4M-APPKI-M409_red.csv',row.names = 1)
#train_coord_y_x<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/4M-APPKI-M409_coord.csv',row.names = 1)
#test_coord_y_x<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/8M-APPKI-M373_coord.csv',row.names = 1)
#train_radius<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/4M-APPKI-M409_radius.csv',row.names = 1)
#test_radius<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/8M-APPKI-M373_radius.csv',row.names = 1)
#img <- load.image("/Volumes/SSD-PUTU3C/Visium loupe file/8M-APPKI-M373/outs/spatial/tissue_lowres_image.png")


#X_train<-t(X_train)
#X_test<-t(X_test)

#common_gene_names<-intersect(colnames(X_train),colnames(X_test))
#X_train<-X_train[,common_gene_names]
#X_test<-X_test[,common_gene_names]

#閾値を指定
threshold<-0.48
y_train[y_train<threshold]<-0
y_train[y_train>=threshold]<-1
y_test[y_test<threshold]<-0
y_test[y_test>=threshold]<-1

#閾値の範囲を指定
#稼動確認用コード（240131太田先生より）
# y1<-as.data.frame(matrix(1,nrow=5,ncol=1))
# y2<-as.data.frame(matrix(0.5,nrow=5,ncol=1))
# y3<-as.data.frame(matrix(0.2,nrow=5,ncol=1))
# y<-rbind(y1,y2,y3)
# threshold_min<-0.4
# threshold_max<-0.8
# 
# print(y)
# y[(y[,1] < threshold_min) | (y[,1] > threshold_max), 1] <- 0　#|の意味は、またはって意味です。threshold_minより小さいか、または、threshold_maxより大きい場合、0を代入
# y[(y[,1] > threshold_min) & (y[,1] < threshold_max), 1] <- 1　#&の意味は、かつって意味です。threshold_minより大きくて、かつ、threshold_maxより小さい場合、0を代入
# print(y)

# #実際のコード（240131太田先生より）
# threshold_min <- 0.48
# threshold_max <- 0.8
# y_train[(y_train<threshold_min)|(y_train>threshold_max)]<-0
# y_test[(y_test<threshold_min)|(y_test>threshold_max)]<-0
# y_train[(y_train>=threshold_min)&(y_train<=threshold_max)]<-1
# y_test[(y_test>=threshold_min)&(y_test>=threshold_max)]<-1


weights_train<-1/table(as.vector(y_train)$V1)[as.factor(as.vector(y_train)$V1)]


#logisticlasso<-glmnet(as.matrix(X_train),y_train,family='binomial',alpha=1)
#logisticlasso<-cv.glmnet(as.matrix(X_train),y_train,family='binomial',alpha=1,weights = weights_train)
# logisticlasso<-cv.glmnet(as.matrix(X_train),as.vector(y_train)$V1,family='binomial',alpha=1,weights = weights_train)
logisticridge<-cv.glmnet(as.matrix(X_train),as.vector(y_train)$V1,family='binomial',alpha=0,weights = weights_train)

# y_test_pred<-predict(logisticlasso,as.matrix(X_test),type='response',s='lambda.min')
y_test_pred<-predict(logisticridge,as.matrix(X_test),type='response',s='lambda.min')
#y_test_pred_bin<-ifelse(y_test_pred>0.5,1,0)
p_threshold<-0.8
y_test_pred_bin<-ifelse(y_test_pred>p_threshold,1,0)

y_wt_test_pred<-y_test_pred[(count_1_1+1):nrow(y_train)]
y_app_test_pred<-y_test_pred[1:count_1_1]
par(mfrow=c(1,2))
hist((y_wt_test_pred),col='red',breaks=20,ylim=c(0,700), xlim=c(0,1))
hist((y_app_test_pred),col='blue',breaks=20,ylim=c(0,700), xlim=c(0,1))
par(mfrow=c(1,1))

TP<-sum(y_test_pred_bin==1&y_test==1)
FP<-sum(y_test_pred_bin==1&y_test==0)
TN<-sum(y_test_pred_bin==0&y_test==0)
FN<-sum(y_test_pred_bin==0&y_test==1)

accuracy<-(TP+TN)/(TP+FP+TN+FN)
precison<-TP/(TP+FP)
recall<-TP/(TP+FN)
F1<-2*precison*recall/(precison+recall)

print(paste('Accuracy=',accuracy))
print(paste('F1=',F1))
print(paste('TP=',TP))
print(paste('FP=',FP))
print(paste('TN=',TN))
print(paste('FN=',FN))

print(paste('total training sample=',nrow(X_train)))
print(paste('training positive sample=',sum(y_train)))
print(paste('total test sample=',nrow(X_test)))
print(paste('test positive sample=',sum(y_test)))

# #範囲を半径19マイクロに広げた結果を表示（250228実施）
y_test_blurred<-read.csv("/Users/okiru/Desktop/240124 Visium_Hiroshima University/250226 Blurred red intensity prediction test/250228_y_blurred_3px_diameter.csv",row.names = 1)

#範囲を半径38マイクロに広げた結果を表示（250226実施）
#y_test_blurred<-read.csv("/Users/okiru/Desktop/240124 Visium_Hiroshima University/250226 Blurred red intensity prediction test/250228_y_blurred_5px_diameter.csv",row.names = 1)


y_test_blurred<-as.data.frame(y_test_blurred[rownames(y_test),])
rownames(y_test_blurred)<-rownames(y_test)
colnames(y_test_blurred)<-"lambda.min"
#y_test_blurred<-as.matrix(y_test_blurred)

#y_test_blurred<-as.matrix(y_test_blurred)
TP<-sum(y_test_pred_bin==1&y_test_blurred==1)
FP<-sum(y_test_pred_bin==1&y_test_blurred==0)
TN<-sum(y_test_pred_bin==0&y_test_blurred==0)
FN<-sum(y_test_pred_bin==0&y_test_blurred==1)

accuracy<-(TP+TN)/(TP+FP+TN+FN)
precison<-TP/(TP+FP)
recall<-TP/(TP+FN)
F1<-2*precison*recall/(precison+recall)

print(paste('Accuracy=',accuracy))
print(paste('F1=',F1))
print(paste('TP=',TP))
print(paste('FP=',FP))
print(paste('TN=',TN))
print(paste('FN=',FN))


# spot_barcodes<-rownames(X_test)
# spot_barcodes<-stringr::str_replace(spot_barcodes,pattern='.1',replace='-1')
# 
# # j=1
# # for(spot_barcode in spot_barcodes){
# #   corespond_test_coord_y_x<-test_coord_y_x[spot_barcode,]
# #   #y_test_pred_j<-y_test_pred[j]
# #   y_test_pred_j<-y_test_pred_bin[j]
# #   img <- img %>%
# #     imager::draw_circle(x = corespond_test_coord_y_x[,2], y = corespond_test_coord_y_x[,1], r = test_radius[1,1], col = "blue",opacity = y_test_pred_j)
# #   j=j+1
# # }
# # plot(img)
# 
# spot_wt_4_coord_y_x<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_wt_coord.csv',row.names = 1)
# spot_app_4_coord_y_x<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_app_coord.csv',row.names = 1)
# spot_wt_8_coord_y_x<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_wt_coord.csv',row.names = 1)
# spot_app_8_coord_y_x<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_app_coord.csv',row.names = 1)
# 
# spot_barcodes_4_month<-rownames(fread_fix('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_scaledata.csv') )
# spot_barcodes_4_month<-stringr::str_replace(spot_barcodes_4_month,pattern='.1',replace='-1')
# spot_barcodes_8_month<-rownames(fread_fix('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_scaledata.csv') )
# spot_barcodes_8_month<-stringr::str_replace(spot_barcodes_8_month,pattern='.1',replace='-1')
# 
# img_wt_4_month <- load.image("/Volumes/SSD-PUTU3C/Visium loupe file/4M-WT-M408/outs/spatial/tissue_lowres_image.png")
# img_app_4_month <- load.image("/Volumes/SSD-PUTU3C/Visium loupe file/4M-APPKI-M409/outs/spatial/tissue_lowres_image.png")
# img_wt_8_month <- load.image("/Volumes/SSD-PUTU3C/Visium loupe file/8M-WT-M359/outs/spatial/tissue_lowres_image.png")
# img_app_8_month <- load.image("/Volumes/SSD-PUTU3C/Visium loupe file/8M-APPKI-M373/outs/spatial/tissue_lowres_image.png")#
# 
# 
# print('Caution: X_train is 8-month data, and X_test is 4-month data')
# #y_pred_4_month<-predict(logisticlasso,as.matrix(X_test),type='response',s='lambda.min')
# y_pred_4_month<-predict(logisticridge,as.matrix(X_test),type='response',s='lambda.min')
# y_pred_4_month_bin<-ifelse(y_pred_4_month>p_threshold,1,0)
# #y_pred_8_month<-predict(logisticlasso,as.matrix(X_train),type='response',s='lambda.min')
# y_pred_8_month<-predict(logisticridge,as.matrix(X_train),type='response',s='lambda.min')
# y_pred_8_month_bin<-ifelse(y_pred_8_month>p_threshold,1,0)
# 
# 
# #グラフィックスデバイスのレイアウトを2×2のグリッドに設定し、4つのプロットを一度に表示する
# par(mfrow=c(2,2))
# #WT 4Mデータの描画
# j=1
# for(spot_barcode in spot_barcodes_4_month){
#   corespond_test_coord_y_x<-spot_wt_4_coord_y_x[spot_barcode,]
#   #y_test_pred_j<-y_test_pred[j]
#   y_test_pred_j<-y_pred_4_month_bin[j]
#   img_wt_4_month <- img_wt_4_month %>%
#     imager::draw_circle(x = corespond_test_coord_y_x[,2], y = corespond_test_coord_y_x[,1], r = test_radius[1,1], col = "gray",opacity = y_test_pred_j)
#   j=j+1
# }
# plot(img_wt_4_month)
# #WT 8Mデータの描画
# j=1
# for(spot_barcode in spot_barcodes_8_month){
#   corespond_test_coord_y_x<-spot_wt_8_coord_y_x[spot_barcode,]
#   #y_test_pred_j<-y_test_pred[j]
#   y_test_pred_j<-y_pred_8_month_bin[j]
#   img_wt_8_month <- img_wt_8_month %>%
#     imager::draw_circle(x = corespond_test_coord_y_x[,2], y = corespond_test_coord_y_x[,1], r = test_radius[1,1], col = "gray",opacity = y_test_pred_j)
#   j=j+1
# }
# plot(img_wt_8_month)
# #APP 4Mデータの描画
# j=1
# for(spot_barcode in spot_barcodes_4_month){
#   corespond_test_coord_y_x<-spot_app_4_coord_y_x[spot_barcode,]
#   #y_test_pred_j<-y_test_pred[j]
#   y_test_pred_j<-y_pred_4_month_bin[j]
#   img_app_4_month <- img_app_4_month %>%
#     imager::draw_circle(x = corespond_test_coord_y_x[,2], y = corespond_test_coord_y_x[,1], r = test_radius[1,1], col = "gray",opacity = y_test_pred_j)
#   j=j+1
# }
# plot(img_app_4_month)
# #APP 8Mデータの描画
# j=1
# for(spot_barcode in spot_barcodes_8_month){
#   corespond_test_coord_y_x<-spot_app_8_coord_y_x[spot_barcode,]
#   #y_test_pred_j<-y_test_pred[j]
#   y_test_pred_j<-y_pred_8_month_bin[j]
#   img_app_8_month <- img_app_8_month %>%
#     imager::draw_circle(x = corespond_test_coord_y_x[,2], y = corespond_test_coord_y_x[,1], r = test_radius[1,1], col = "gray",opacity = y_test_pred_j)
#   j=j+1
# }
# plot(img_app_8_month)
# #レイアウトを元に戻す
# par(mfrow=c(1,1))
# print('Caution: X_train is 8-month data, and X_test is 4-month data')
# 
# 
# # y_test_pred_app<-c()
# # # for(spot_barcode in spot_barcodes){
# # #   y
# # #
# # # }
# # par(mfrow=c(1,2))
# # hist(y_test_pred[1:dim(test_coord_y_x)[1]],col='red',breaks=20,ylim=c(0,700), xlim=c(0,1))
# # hist(y_test_pred[(dim(test_coord_y_x)[1]+1):nrow(y_test_pred)],col='blue',breaks=20,ylim=c(0,700), xlim=c(0,1))
# 
# #Lasso回帰モデルの係数取得
# # lr_coefficients <- coef(logisticlasso, s='lambda.min')
# # #係数のデータフレームへの変換
# # lr_coefficients_mat <- as.matrix(lr_coefficients)
# # df_coeffs <- as.data.frame(lr_coefficients_mat)
# # df_coeffs$genename <- rownames(df_coeffs)
# # #不要な（Intercept)行を除外
# # df_coeffs <- df_coeffs[df_coeffs$genename != "(Intercept)",]
# # #係数の絶対値で降順に並べる
# # df_coeffs <- df_coeffs[order(abs(df_coeffs$s1), decreasing=TRUE),]
# # #閾値を設定して重要な係数を抽出
# # effective_coeffs <- df_coeffs[abs(df_coeffs$s1) >0.0001,]
# # #視覚化（棒グラフ）
# # ggplot(effective_coeffs, aes(x=reorder(genename, s1), y=s1)) + geom_bar(stat = "identity") + theme_minimal() + labs(title = "Effective genes for logistic regression", x="Gene", y="Weights") +
# #   coord_flip()  # 横向きにして遺伝子名が見やすくする
# # 
# # print(effective_coeffs)
# # print(effective_coeffs[order(effective_coeffs$s1),])
# 
# #write.table(effective_coeffs,"/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/effective_coeffs_8M_train_4M_test.txt", sep = "\t")
# #write.table(effective_coeffs,"/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/effective_coeffs_8M_train_4M_WT_test.txt", sep = "\t")
# #write.table(effective_coeffs,"/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/effective_coeffs_8M_train_4M_test_threshold_0_8_0_48.txt", sep = "\t")
# 
# ##c19遺伝子リストに絞って学習した結果
# #write.table(effective_coeffs,"/Users/okiru/Desktop/240124 Visium_Hiroshima University/241216 prediction with c19 genes/Lasso_effective_coeffs_8M_c19genes_train_4M_test.txt", sep = "\t")
# 
# # #DAM遺伝子(FDRp<0.05, FCでは切っていないリスト)のアップロード
# # DAM_genes <- read.xlsx("/Users/okiru/Desktop/DAM marker genes from Karen-Shaul et al 2017 /210907 DAM marker genes FDRp<0.05 .xlsx")
# # #cluster19のgene listを指定する。
# # DAM_gene_names<-DAM_genes$Gene.name
# # head(DAM_gene_names)
# # #遺伝子名がリストに含まれているものを抽出
# # selected_genes_lasso_coeffs <- effective_coeffs[effective_coeffs$genename %in% DAM_gene_names, ]
# # 
# # # 結果の確認
# # head(selected_genes_lasso_coeffs)
# # 
# # write.table(selected_genes_lasso_coeffs,"/Users/okiru/Desktop/240124 Visium_Hiroshima University/241216 prediction with c19 genes/DAM_lasso_effective_coeffs_8M_c19genes_train_4M_test.txt", sep = "\t")
# 
# 
# #Ridge回帰モデルの係数を取得
# ridge_coefficients <- coef(logisticridge, s = "lambda.min")  # lambda.minで最適な正則化パラメータ
# #係数のデータフレームへの変換
# ridge_coefficients_mat <- as.matrix(ridge_coefficients)
# df_ridge_coeffs <- as.data.frame(ridge_coefficients_mat)
# df_ridge_coeffs$genename <- rownames(df_ridge_coeffs)
# 
# # 切片を除外
# df_ridge_coeffs <- df_ridge_coeffs[df_ridge_coeffs$genename != "(Intercept)",]
# 
# # 係数の絶対値で降順にソート
# df_ridge_coeffs <- df_ridge_coeffs[order(abs(df_ridge_coeffs$s1), decreasing = TRUE),]
# 
# # 可視化: 全特徴の係数を表示
# ggplot(df_ridge_coeffs, aes(x = reorder(genename, s1), y = s1)) +
#   geom_bar(stat = "identity") +
#   theme_minimal() +
#   labs(title = "Ridge Regression Coefficients",
#        x = "Gene",
#        y = "Weight Coefficients") +
#   coord_flip()  # 横向きにして遺伝子名が見やすくする
# 
# print(df_ridge_coeffs)
# print(df_ridge_coeffs[order(df_ridge_coeffs$s1),])
# 
# write.table(df_ridge_coeffs,"/Users/okiru/Desktop/240124 Visium_Hiroshima University/241216 prediction with c19 genes/250207_Ridge_effective_coeffs_8M_c19genes_train_4M_test.txt", sep = "\t")
# 
# # 可視化: 上位50位の正と負の重みを表示-----------------------------------
# # 上位50の正の重み
# top_positive <- df_ridge_coeffs %>%
#   arrange(desc(s1)) %>%
#   head(50)
# 
# # 上位50の負の重み
# top_negative <- df_ridge_coeffs %>%
#   arrange(s1) %>%
#   head(50)
# 
# # データを結合
# top_coeffs <- bind_rows(top_positive, top_negative)
# 
# # プロット
# ggplot(top_coeffs, aes(x = reorder(genename, s1), y = s1, fill = s1 > 0)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c("blue", "red"), labels = c("Negative", "Positive")) +
#   theme_minimal() +
#   labs(title = "Top 50 Positive & Negative Ridge Regression Coefficients",
#        x = "Gene",
#        y = "Weight Coefficients",
#        fill = "Coefficient Sign") +
#   coord_flip() +
#   theme(
#     plot.title = element_text(size = 14, face = "bold", color = "darkblue"),
#     axis.title.x = element_text(size = 14, color = "black"),
#     axis.title.y = element_text(size = 14, color = "black"),
#     axis.text.x = element_text(size = 12, color = "black"),
#     axis.text.y = element_text(size = 10, color = "black"),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 10)
#   )
# 
# # 可視化: 上位20位の正と負の重みを表示-----------------------------------
# # 上位20の正の重み
# top_positive <- df_ridge_coeffs %>%
#   arrange(desc(s1)) %>%
#   head(20)
# 
# # 上位20の負の重み
# top_negative <- df_ridge_coeffs %>%
#   arrange(s1) %>%
#   head(20)
# 
# # データを結合
# top_coeffs <- bind_rows(top_positive, top_negative)
# 
# # プロット
# ggplot(top_coeffs, aes(x = reorder(genename, s1), y = s1, fill = s1 > 0)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c("blue", "red"), labels = c("Negative", "Positive")) +
#   theme_minimal() +
#   labs(title = "Top 20 Positive & Negative Ridge Regression Coefficients",
#        x = "Gene",
#        y = "Weight Coefficients",
#        fill = "Coefficient Sign") +
#   coord_flip() +
#   theme(
#     plot.title = element_text(size = 14, face = "bold", color = "darkblue"),
#     axis.title.x = element_text(size = 14, color = "black"),
#     axis.title.y = element_text(size = 14, color = "black"),
#     axis.text.x = element_text(size = 12, color = "black"),
#     axis.text.y = element_text(size = 10, color = "black"),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 10)
#   )
# 
# # #DAM遺伝子(FDRp<0.05, FCでは切っていないリスト)のアップロード
# # DAM_genes <- read.xlsx("/Users/okiru/Desktop/DAM marker genes from Karen-Shaul et al 2017 /210907 DAM marker genes FDRp<0.05 .xlsx")
# 
# #DAM遺伝子(FDRp<0.05, FC1.2 Up&Dnリスト)のアップロード
# DAM_genes <- read.xlsx("/Users/okiru/Documents/論文投稿/APPKI_Visium_paper/250131_Comparison between c19 and DAM genes/250131 DAM marker genes logFC_0.25 UpDn FDRp<0.05_modified.xlsx")
# 
# #cluster19のgene listを指定する。
# DAM_gene_names<-DAM_genes$Gene.name
# head(DAM_gene_names)
# #遺伝子名がリストに含まれているものを抽出
# selected_genes_ridge_coeffs <- df_ridge_coeffs[df_ridge_coeffs$genename %in% DAM_gene_names, ]
# 
# # 結果の確認
# head(selected_genes_ridge_coeffs)
# 
# # write.table(selected_genes_ridge_coeffs,"/Users/okiru/Desktop/240124 Visium_Hiroshima University/241216 prediction with c19 genes/DAM_ridge_effective_coeffs_8M_c19genes_train_4M_test.txt", sep = "\t")
# write.table(selected_genes_ridge_coeffs,"/Users/okiru/Desktop/240124 Visium_Hiroshima University/241216 prediction with c19 genes/250213_DAM_ridge_effective_coeffs_8M_c19genes_train_4M_test.txt", sep = "\t")
# 
# # 可視化: 上位50位の正と負の重みを表示-----------------------------------
# # 上位50の正の重み
# top_positive_DAM <- selected_genes_ridge_coeffs %>%
#   arrange(desc(s1)) %>%
#   head(50)
# 
# # 上位50の負の重み
# top_negative_DAM <- selected_genes_ridge_coeffs %>%
#   arrange(s1) %>%
#   head(50)
# 
# # データを結合
# top_coeffs_DAM <- bind_rows(top_positive_DAM, top_negative_DAM)
# 
# # プロット
# ggplot(top_coeffs_DAM, aes(x = reorder(genename, s1), y = s1, fill = s1 > 0)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c("blue", "red"), labels = c("Negative", "Positive")) +
#   theme_minimal() +
#   labs(title = "Top 50 DAM genes Positive & Negative Ridge Regression Coefficients",
#        x = "DAM Gene",
#        y = "Weight Coefficients",
#        fill = "Coefficient Sign") +
#   coord_flip() +
#   theme(
#     plot.title = element_text(size = 14, face = "bold", color = "darkblue"),
#     axis.title.x = element_text(size = 14, color = "black"),
#     axis.title.y = element_text(size = 14, color = "black"),
#     axis.text.x = element_text(size = 12, color = "black"),
#     axis.text.y = element_text(size = 10, color = "black"),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 10)
#   )
# 
# # 可視化: 上位20位の正と負の重みを表示-----------------------------------
# # 上位20の正の重み
# top_positive_DAM <- selected_genes_ridge_coeffs %>%
#   arrange(desc(s1)) %>%
#   head(20)
# 
# # 上位50の負の重み
# top_negative_DAM <- selected_genes_ridge_coeffs %>%
#   arrange(s1) %>%
#   head(20)
# 
# # データを結合
# top_coeffs_DAM <- bind_rows(top_positive_DAM, top_negative_DAM)
# 
# # プロット
# ggplot(top_coeffs_DAM, aes(x = reorder(genename, s1), y = s1, fill = s1 > 0)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c("blue", "red"), labels = c("Negative", "Positive")) +
#   theme_minimal() +
#   labs(title = "Top 20 DAM genes Positive & Negative Ridge Regression Coefficients",
#        x = "DAM Gene",
#        y = "Weight Coefficients",
#        fill = "Coefficient Sign") +
#   coord_flip() +
#   theme(
#     plot.title = element_text(size = 14, face = "bold", color = "darkblue"),
#     axis.title.x = element_text(size = 14, color = "black"),
#     axis.title.y = element_text(size = 14, color = "black"),
#     axis.text.x = element_text(size = 12, color = "black"),
#     axis.text.y = element_text(size = 10, color = "black"),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 10)
#   )
# 
# #Permutation test用の座標を取得
# # 値が1の行を取得
# indices_4_month <- which(y_pred_4_month_bin[, 1] == 1)
# indices_8_month <- which(y_pred_8_month_bin[, 1] == 1)
# # 行名（スポットバーコード）を取得
# y_pred_4_month_spot_barcodes <- rownames(y_pred_4_month_bin)[indices_4_month]
# y_pred_8_month_spot_barcodes <- rownames(y_pred_8_month_bin)[indices_8_month]
# # y_pred_4_month_spot_barcodesに一致する行名を持つ行を抽出
# y_pred_4_month_app_coord_y_x <- spot_app_4_coord_y_x[rownames(spot_app_4_coord_y_x) %in% y_pred_4_month_spot_barcodes, ]
# y_pred_8_month_app_coord_y_x <- spot_app_8_coord_y_x[rownames(spot_app_8_coord_y_x) %in% y_pred_8_month_spot_barcodes, ]
# #データの保存
# write.csv(y_pred_4_month_app_coord_y_x,"/Users/okiru/Desktop/240124 Visium_Hiroshima University/241216 prediction with c19 genes/y_pred_4_month_app_coord_y_x.csv", row.names = TRUE)
# write.csv(y_pred_8_month_app_coord_y_x,"/Users/okiru/Desktop/240124 Visium_Hiroshima University/241216 prediction with c19 genes/y_pred_8_month_app_coord_y_x.csv", row.names = TRUE)
# 
# 
