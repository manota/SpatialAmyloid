app_4_coord<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_app_coord.csv')
red_4<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/learning_data/4month_red.csv')
y_pred_4_month_app_coord_y_x<-read.csv('/Users/okiru/Desktop/240124 Visium_Hiroshima University/241216 prediction with c19 genes/y_pred_4_month_app_coord_y_x.csv')

#特定の行をフィルタリング（Appのみにする）
rows_to_leave <- grep("1_1", red_4[, 1])
red_4 <- red_4[rows_to_leave, , drop = FALSE]
# 1列目のバーコードを行名に設定
rownames(app_4_coord)<-app_4_coord[,1]
rownames(red_4) <- red_4[, 1]
rownames(y_pred_4_month_app_coord_y_x)<-y_pred_4_month_app_coord_y_x[,1]
# 1列目を削除
app_4_coord<-app_4_coord[,-1, drop = FALSE]
red_4 <- red_4[, -1, drop = FALSE]
y_pred_4_month_app_coord_y_x<-y_pred_4_month_app_coord_y_x[,-1, drop = FALSE]

#赤色（Aβ）のフィルタリング
red_threshold<-0.48
row_filter <- apply(red_4, 1, function(row) any(row > red_threshold))
red_4_app <- red_4[row_filter, , drop = FALSE]
#赤色の座標を抽出
red_4_app_coord<-app_4_coord[rownames(red_4_app),]

###Using data
red_4_app_coord
app_4_coord
y_pred_4_month_app_coord_y_x

# ユークリッド距離を計算する関数
calc_distance<-function(first_coord,second_coord){
  return(sqrt(sum((first_coord-second_coord)^2)))
}

# calculate all patterns (PredictionのスポットとAβとの最短距離の総和を計算）
y_pred_4_month_app_min<-0    ###計算枠のようなものをつくるため0を代入（初期化）
for(i in 1:nrow(y_pred_4_month_app_coord_y_x)){
  cell_coord<-y_pred_4_month_app_coord_y_x[i,c(1,2)]
  cell_coord<-as.numeric(cell_coord)
  all_distance<-apply(red_4_app_coord,1,function(x) calc_distance(x,cell_coord))
  min_distance<-min(all_distance)
  y_pred_4_month_app_min<-y_pred_4_month_app_min+min_distance
  print(y_pred_4_month_app_min)  #進行状況を表示
}

print(paste("Total minimum distance:", y_pred_4_month_app_min))

#####Predicitonのスポットの数と同じスポットをランダムに選んでそれを繰り返した場合の最短距離の総和

sampling_num<-nrow(y_pred_4_month_app_coord_y_x)


print(Sys.time())
repeat_times<-10000   ####change yourself 10000回に設定(Takes 17min)
repeat_total_mins <- numeric(repeat_times)

# 置換検定の実行
for (rep_i in 1:repeat_times) {
  # ランダムにサンプリング
  sampled_indices <- sample(1:nrow(app_4_coord), sampling_num)
  sampling_coords <- app_4_coord[sampled_indices, ]
  
  # 最小距離の合計を計算
  total_min <- 0
  for (i in 1:nrow(sampling_coords)) {
    cell_coord <- as.numeric(sampling_coords[i, c(1, 2)])
    all_distance <- apply(red_4_app_coord, 1, function(x) calc_distance(x, cell_coord))
    min_distance <- min(all_distance)
    total_min <- total_min + min_distance
  }
  
  # 繰り返しの結果を保存
  repeat_total_mins[rep_i] <- total_min
}

# 現在の時刻の表示
print(Sys.time())

# p値の計算
p_value <- max(1/repeat_times, sum(repeat_total_mins <= y_pred_4_month_app_min) / repeat_times)
print(paste("p value:", p_value))

####Draw histogram (scale pixels)
library(jsonlite)
scale_factor_app_4m <-read_json('/Volumes/SSD-PUTU3C/Visium loupe file/4M-APPKI-M409/outs/spatial/scalefactors_json.json')
hist(repeat_total_mins/sampling_num/scale_factor_app_4m$tissue_lowres_scalef,xlim=c(150,600),ylim = c(0,2000),main="Permutation Test Distribution", xlab="Averaged minimum distance (pixel)", ylab="Number of Occurrences/10000", col="lightblue", border="black")
abline(v = y_pred_4_month_app_min/sampling_num/scale_factor_app_4m$tissue_lowres_scalef, col = "red", lwd = 2, lty = 2)


####Change scale pixels to microns (Ref. 10xgenomics website: scale factor)
# microns_pixels <- 65/scale_factor_app_4m$spot_diameter_fullres
# hist(repeat_total_mins/sampling_num/scale_factor_app_4m$tissue_lowres_scalef*microns_pixels,xlim=c(100,450), ylim = c(0,3000),main="Permutation Test Distribution", xlab="Averaged minimum distance (micron)", ylab="Number of Occurrences/10000", col="lightblue", border="black")
# abline(v = y_pred_4_month_app_min/sampling_num/scale_factor_app_4m$tissue_lowres_scalef*microns_pixels, col = "red", lwd = 2, lty = 2)

###250130 modified ラベルなどのサイズを変更
microns_pixels <- 65/scale_factor_app_4m$spot_diameter_fullres
hist(repeat_total_mins / sampling_num / scale_factor_app_4m$tissue_lowres_scalef * microns_pixels,
     xlim = c(100, 450), ylim = c(0, 3000),
     main = "Permutation Test Distribution", 
     xlab = "Averaged minimum distance (micron)", 
     ylab = "Number of Occurrences/10000", 
     col = "lightblue", border = "black",
     cex.main = 1.5,  # タイトルのフォントサイズ
     cex.lab = 1.5,   # 軸ラベルのフォントサイズ
     cex.axis = 1.5   # 目盛りのフォントサイズ
)

abline(v = y_pred_4_month_app_min/sampling_num/scale_factor_app_4m$tissue_lowres_scalef*microns_pixels, col = "red", lwd = 2, lty = 2)


