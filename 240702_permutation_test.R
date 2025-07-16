app_8_coord<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_app_coord.csv')
app_4_coord<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_app_coord.csv')
red_8<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/8month_red.csv')
red_4<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_red.csv')
c19<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/saved_cluster_barcodes_strain/19.csv', row.names=1)


# サフィックスを追加する関数
add_suffix <- function(x, suff) {
  return(paste0(x, suff))
}

# データフレームの各要素に関数を適用
app_8_coord[,1] <- lapply(app_8_coord, function(col) sapply(col, add_suffix, suff = "_2"))
app_4_coord[,1] <- lapply(app_4_coord, function(col) sapply(col, add_suffix, suff = "_1"))
red_8[,1]<-lapply(red_8, function(col) sapply(col, add_suffix, suff = "_2"))
red_4[,1]<-lapply(red_4, function(col) sapply(col, add_suffix, suff = "_1"))

# 行名を設定
rownames(app_8_coord)<-app_8_coord[,1]
rownames(app_4_coord)<-app_4_coord[,1]
rownames(red_8)<-red_8[,1]
rownames(red_4)<-red_4[,1]

# 最初の列を削除
app_8_coord <- app_8_coord[, -1, drop = FALSE]
app_4_coord <- app_4_coord[, -1, drop = FALSE]
red_8 <- red_8[, -1, drop = FALSE]
red_4 <- red_4[, -1, drop = FALSE]

# 特定の行をフィルタリング
rows_to_leave <- grep("1_1_2", rownames(red_8))
red_8 <- red_8[rows_to_leave, , drop = FALSE]
rows_to_leave <- grep("1_1_1", rownames(red_4))
red_4 <- red_4[rows_to_leave, , drop = FALSE]

# クラスター19の座標を抽出
app8_c19<-app_8_coord[c19$barcodes,]
app8_c19[is.na(app8_c19)]<-0
app4_c19<-app_4_coord[c19$barcodes,]
app4_c19[is.na(app4_c19)]<-0
rownames(app8_c19)<-c19$barcodes
rownames(app4_c19)<-c19$barcodes

# 8ヶ月と4ヶ月のデータを合計
c19_coord<-app8_c19+app4_c19
c19_coord[,3]<-c19[,2]
c19_8_coord<-c19_coord[c19_coord[, 3] == "APPKI_8M", ]
c19_4_coord<-c19_coord[c19_coord[, 3] == "APPKI_4M", ]

# 赤色（Aβ）のフィルタリング
red_threshold<-0.48
row_filter <- apply(red_8, 1, function(row) any(row > red_threshold))
red_8_app <- red_8[row_filter, , drop = FALSE]
row_filter <- apply(red_4, 1, function(row) any(row > red_threshold))
red_4_app <- red_4[row_filter, , drop = FALSE]

# 赤色の座標を抽出
red_8_app_coord<-app_8_coord[rownames(red_8_app),]
red_4_app_coord<-app_4_coord[rownames(red_4_app),]

###using data
c19_8_coord
c19_4_coord
red_8_app_coord
red_4_app_coord
app_8_coord
app_4_coord


# ユークリッド距離を計算する関数
calc_distance<-function(first_coord,second_coord){
  return(sqrt(sum((first_coord-second_coord)^2)))
}

# calculate all patterns (c19とAβとの最短距離の総和を計算）
c19_total_min<-0    ###計算枠のようなものをつくるため0を代入（初期化）
for(i in 1:nrow(c19_8_coord)){
  cell_coord<-c19_8_coord[i,c(1,2)]
  cell_coord<-as.numeric(cell_coord)
  all_distance<-apply(red_8_app_coord,1,function(x) calc_distance(x,cell_coord))
  min_distance<-min(all_distance)
  c19_total_min<-c19_total_min+min_distance
  print(c19_total_min)  #進行状況を表示
}

print(paste("Total minimum distance:", c19_total_min))

#####C19の数と同じスポットをランダムに選んでそれを繰り返した場合の最短距離の総和

sampling_num<-nrow(c19_8_coord)


print(Sys.time())
repeat_times<-10000   ####change yourself 10000回に設定(Takes 17min)
repeat_total_mins <- numeric(repeat_times)
#repeat_total_mins<-matrix(0,repeat_times,1)

# 置換検定の実行
for (rep_i in 1:repeat_times) {
  # ランダムにサンプリング
  sampled_indices <- sample(1:nrow(app_8_coord), sampling_num)
  sampling_coords <- app_8_coord[sampled_indices, ]
  
  # 最小距離の合計を計算
  total_min <- 0
  for (i in 1:nrow(sampling_coords)) {
    cell_coord <- as.numeric(sampling_coords[i, c(1, 2)])
    all_distance <- apply(red_8_app_coord, 1, function(x) calc_distance(x, cell_coord))
    min_distance <- min(all_distance)
    total_min <- total_min + min_distance
  }
  
  # 繰り返しの結果を保存
  repeat_total_mins[rep_i] <- total_min
}

# 現在の時刻の表示
print(Sys.time())

# p値の計算
p_value <- max(1/repeat_times, sum(repeat_total_mins <= c19_total_min) / repeat_times)
print(paste("p value:", p_value))

# for(rep_i in 1:repeat_times){
#   sampled_indices <- sample(1:nrow(app_8_coord), sampling_num)
#   sampling_coords<-app_8_coord[sampled_indices,]
#   total_min<-0
#   for(i in 1:nrow(sampling_coords)){
#     cell_coord<-sampling_coords[i,c(1,2)]
#     cell_coord<-as.numeric(cell_coord)
#     all_distance<-apply(red_8_app_coord,1,function(x) calc_distance(x,cell_coord))
#     min_distance<-min(all_distance)
#     total_min<-total_min+min_distance
#   }
#   repeat_total_mins[rep_i]<-total_min
# }
# print(Sys.time())
# p_value <- sum(repeat_total_mins <= c19_total_min) / repeat_times
# print(paste("p value:", p_value))

####Draw histogram (scale pixels)
library(jsonlite)
scale_factor_app_8m <-read_json('/Volumes/SSD-PUTU3C/Visium loupe file/8M-APPKI-M373/outs/spatial/scalefactors_json.json')
hist(repeat_total_mins/sampling_num/scale_factor_app_8m$tissue_lowres_scalef,xlim=c(150,400),ylim = c(0,2000),main="Permutation Test Distribution", xlab="Averaged minimum distance (pixel)", ylab="Number of Occurrences/10000", col="lightblue", border="black")
abline(v = c19_total_min/sampling_num/scale_factor_app_8m$tissue_lowres_scalef, col = "red", lwd = 2, lty = 2)


####Change scale pixels to microns (Ref. 10xgenomics website: scale factor)
hist(repeat_total_mins / sampling_num / scale_factor_app_8m$tissue_lowres_scalef * microns_pixels,
     xlim = c(100, 350), ylim = c(0, 2500),
     main = "Permutation Test Distribution", 
     xlab = "Averaged minimum distance (micron)", 
     ylab = "Number of Occurrences/10000", 
     col = "lightblue", border = "black",
     cex.main = 1.5,  # タイトルのフォントサイズ
     cex.lab = 1.5,   # 軸ラベルのフォントサイズ
     cex.axis = 1.5   # 目盛りのフォントサイズ
)

abline(v = c19_total_min / sampling_num / scale_factor_app_8m$tissue_lowres_scalef * microns_pixels, 
       col = "red", lwd = 2, lty = 2)







