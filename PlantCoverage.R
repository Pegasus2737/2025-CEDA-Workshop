## 資料讀取
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(vegan)
library(infotheo)

# df <- read.csv("Plant/raw.csv", encoding = "UTF-8")
# 篩掉芒草跟次生林的大區資料，不要列入分析
df <- read.csv("Plant/raw-copy.csv", encoding = "UTF-8")

# 1. 各地區樹種分布與DBA總和
area_species <- df %>%
    group_by(大區, 樹種) %>%
    summarise(DBA_sum = sum(DBA), .groups = "drop")
print(head(area_species))

# 2. 各地區樹種DBA長條圖
png("Plant/area_species_dba_barchart.png", width = 1000, height = 800)
ggplot(area_species, aes(x = 樹種, y = DBA_sum, fill = 大區)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "各地區樹種DBA總和", y = "DBA總和")
dev.off()

# 3. 優勢樹種分析（全區DBA總和排序）
dominant_species <- df %>%
    group_by(樹種) %>%
    summarise(DBA_total = sum(DBA), .groups = "drop") %>%
    arrange(desc(DBA_total))
print(head(dominant_species, 10))

# 4. DBA分布統計
summary(df$DBA)
cat("DBA標準差：", sd(df$DBA), "\n")

# 5. DBA分布視覺化
png("Plant/dba_distribution_histogram.png", width = 800, height = 600)
ggplot(df, aes(x = DBA)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    labs(title = "DBA分布直方圖")
dev.off()

# 5.1 Gapped Histogram
set.seed(123)
data <- df$DBA

# 階層式分群
h <- hclust(dist(data))
clus <- 24 # 可依資料調整分群數
cl <- cutree(h, clus)
tcl <- table(cl)
mm <- rbind(tapply(data, cl, min), tapply(data, cl, max))

png("Plant/Plant_gapped_histogram.png", width=2000, height=1500)
plot(NULL, 
     xlim=c(min(data), max(data)),
     ylim=c(0, max(tcl)), 
     ylab="freq.", xlab="DBA", 
     main="Plant DBA Gapped Histogram")
for (i in 1:ncol(mm)){
  rect(xleft=mm[1, i ], 
       ybottom=0,
       xright=mm[2, i ],
       ytop=tcl[i], 
       col=rgb(runif(1), runif(1), runif(1), 0.5), 
       border="blue")
}
dev.off()

# 5.2 ECDF Plot
png("Plant/Plant_ecdf_plot.png", width=8000, height=6000)
Fn <- ecdf(data)
plot(Fn, col=8, main="Plant DBA ECDF Plot")
tdata <- cbind(data, cl)
utdata <- unique(tdata)
mutdata <- tapply(utdata[,1], utdata[,2], range)
for (i in 1:length(mutdata)){
  lines(mutdata[[i]], Fn(mutdata[[i]]),
        col=2, lwd=2)
}
dev.off()

# 6. 熱力圖：地區 x 樹種 DBA
mat <- area_species %>%
    pivot_wider(names_from = 樹種, values_from = DBA_sum, values_fill = 0)
mat2 <- as.data.frame(mat)
rownames(mat2) <- mat2$大區
mat2$大區 <- NULL
png("Plant/area_species_dba_heatmap.png", width = 2000, height = 1500)
pheatmap(as.matrix(mat2),
                 main = "地區 x 樹種 DBA 熱力圖",
                 fontsize_row = 8,
                 fontsize_col = 8)
dev.off()

# 6.1. 把樹種組成群之後再畫熱力圖，並標註每群包含的樹種
set.seed(123)
species_data <- dominant_species$DBA_total
species_hclust <- hclust(dist(species_data))
species_clusters <- cutree(species_hclust, k = 8)
# 加入群組資訊
dominant_species$Cluster <- species_clusters

# 標註每群包含的樹種
cluster_species <- dominant_species %>%
    group_by(Cluster) %>%
    summarise(Species = paste(樹種, collapse = ", "), .groups = "drop")
print(cluster_species)

# 計算每群的DBA總和
clustered_species <- area_species %>%
    left_join(dominant_species %>% select(樹種, Cluster), by = "樹種") %>%
    group_by(大區, Cluster) %>%
    summarise(DBA_cluster_sum = sum(DBA_sum), .groups = "drop")

# 轉換為矩陣格式
cluster_matrix <- clustered_species %>%
    pivot_wider(names_from = Cluster, values_from = DBA_cluster_sum, values_fill = 0)
cluster_matrix2 <- as.data.frame(cluster_matrix)
rownames(cluster_matrix2) <- cluster_matrix2$大區
cluster_matrix2$大區 <- NULL

# 畫熱力圖
png("Plant/area_species_cluster_dba_heatmap.png", width = 1200, height = 1000)
pheatmap(as.matrix(cluster_matrix2),
        main = "地區 x 樹種群組 DBA 熱力圖",
        fontsize_row = 8,
        fontsize_col = 8)
dev.off()

# 7. 群落多樣性指數（Shannon）
area_matrix <- area_species %>%
    pivot_wider(names_from = 樹種, values_from = DBA_sum, values_fill = 0)
area_matrix2 <- as.data.frame(area_matrix)
rownames(area_matrix2) <- area_matrix2$大區
area_matrix2$大區 <- NULL
shannon <- diversity(area_matrix2, index = "shannon")
shannon_df <- data.frame(大區 = rownames(area_matrix2), Shannon = shannon)
print(shannon_df)

# 8. 多樣性指數折線圖
png("Plant/shannon_diversity_lineplot.png", width = 800, height = 600)
ggplot(shannon_df, aes(x = 大區, y = Shannon, group = 1)) +
    geom_line() + geom_point() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "各地區Shannon多樣性指數", y = "Shannon指數")
dev.off()

# 9. 標準化各區域樹種DBA並視覺化
standardized_area_species <- area_species %>%
    group_by(大區) %>%
    mutate(DBA_standardized = DBA_sum / sum(DBA_sum)) %>%
    ungroup()
png("Plant/standardized_area_species_dba_barchart.png", width = 1000, height = 800)
ggplot(standardized_area_species, aes(x = 樹種, y = DBA_standardized, fill = 大區)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "各地區樹種DBA標準化", y = "DBA標準化值")
dev.off()

# 9.1 標準化DBA Gapped Histogram
set.seed(123)
standardized_data <- standardized_area_species$DBA_standardized
h <- hclust(dist(standardized_data))
standardized_data
clus <- 5 # 可依資料調整分群數
cl <- cutree(h, clus)
tcl <- table(cl)
mm <- rbind(tapply(standardized_data, cl, min), tapply(standardized_data, cl, max))
png("Plant/Plant_gapped_histogram_standardized.png", width=2000, height=1500)
plot(NULL, 
     xlim=c(min(standardized_data), max(standardized_data)),
     ylim=c(0, max(tcl)), 
     ylab="freq.", xlab="DBA標準化值", 
     main="Plant 樹種 DBA 標準化 Gapped Histogram")
for (i in 1:ncol(mm)){
  rect(xleft=mm[1, i ], 
       ybottom=0,
       xright=mm[2, i ],
       ytop=tcl[i], 
       col=rgb(runif(1), runif(1), runif(1), 0.5), 
       border="blue")
}
dev.off()

# 9.2 標準化DBA ECDF Plot
png("Plant/Plant_ecdf_plot_standardized.png", width=8000, height=6000)
Fn <- ecdf(standardized_data)
plot(Fn, col=8, main="Plant 樹種 DBA 標準化 ECDF Plot")
tdata <- cbind(standardized_data, cl)
utdata <- unique(tdata)
mutdata <- tapply(utdata[,1], utdata[,2], range)
for (i in 1:length(mutdata)){
  lines(mutdata[[i]], Fn(mutdata[[i]]),
        col=2, lwd=2)
}
dev.off()

# 10. 樹種DBA標準化熱力圖
mat_standardized <- standardized_area_species %>%
    dplyr::select(大區, 樹種, DBA_standardized) %>%
    tidyr::pivot_wider(
        names_from = 樹種,
        values_from = DBA_standardized,
        values_fill = 0
    )
mat_standardized2 <- as.data.frame(mat_standardized)
rownames(mat_standardized2) <- mat_standardized2$大區
mat_standardized2$大區 <- NULL
png("Plant/PlantCoverage_standarized_heatmap.png", width = 2000, height = 1500)
pheatmap(as.matrix(mat_standardized2),
                 main = "地區 x 樹種 DBA 標準化熱力圖",
                 clustering_method = "ward.D2",
                 clustering_distance_cols = "manhattan",
                 clustering_distance_rows = "manhattan",
                 fontsize_row = 8,
                 fontsize_col = 8)
dev.off()

# 10.1. 用曼哈頓距離分樹種群
# 轉置矩陣，讓 row = 樹種, col = 地區
species_matrix <- t(as.matrix(mat_standardized2))

# 計算樹種間的曼哈頓距離
species_dist <- dist(species_matrix, method = "manhattan")

# 階層式分群
species_hclust <- hclust(species_dist, method = "ward.D2")
png("Plant/species_hclust_manhattan.png", width = 1200, height = 800)
plot(species_hclust, main = "樹種分群樹狀圖（曼哈頓距離）")
dev.off()

# 指定群數
species_clusters <- cutree(species_hclust, k = 14)

# 加入群組資訊
species_cluster_df <- data.frame(樹種 = rownames(species_matrix), Cluster = species_clusters)

# 標註每群包含的樹種
cluster_species <- species_cluster_df %>%
    group_by(Cluster) %>%
    summarise(Species = paste(樹種, collapse = ", "), .groups = "drop")
print(cluster_species, n = 15)

# 計算每群的 DBA 標準化總和
standardized_clustered_species <- standardized_area_species %>%
    left_join(species_cluster_df, by = "樹種") %>%
    group_by(大區, Cluster) %>%
    summarise(DBA_cluster_standardized_sum = sum(DBA_standardized), .groups = "drop")

# 轉換為矩陣格式
standardized_cluster_matrix <- standardized_clustered_species %>%
    pivot_wider(names_from = Cluster, values_from = DBA_cluster_standardized_sum, values_fill = 0)
standardized_cluster_matrix2 <- as.data.frame(standardized_cluster_matrix)
rownames(standardized_cluster_matrix2) <- standardized_cluster_matrix2$大區
standardized_cluster_matrix2$大區 <- NULL

# 畫標準化熱力圖
png("Plant/PlantCoverage_standarized_heatmap_cluster_manhattan.png", width = 4000, height =3000)
# 確保 annotation_col 的 rownames 和矩陣的 colnames 順序和數量一致
cluster_names_in_matrix <- colnames(standardized_cluster_matrix2)
annotation_col_df <- cluster_species %>% 
  filter(Cluster %in% cluster_names_in_matrix) %>%
  arrange(match(Cluster, cluster_names_in_matrix))
annotation_col <- data.frame(Species = annotation_col_df$Species, row.names = annotation_col_df$Cluster)


pheatmap(as.matrix(standardized_cluster_matrix2),
         main = "地區 x 樹種群組 DBA 標準化熱力圖（曼哈頓距離）",
         clustering_method = "ward.D2",
         clustering_distance_cols = "manhattan",
         clustering_distance_rows = "manhattan",
         fontsize_row = 8,
         fontsize_col = 8,
         annotation_col = annotation_col)
dev.off()

# 10.2 用曼哈頓距離分地區群
# 計算地區間的曼哈頓距離
area_matrix <- as.matrix(mat_standardized2)
area_dist <- dist(area_matrix, method = "manhattan")

# 階層式分群
area_hclust <- hclust(area_dist, method = "ward.D2")
png("Plant/area_hclust_manhattan.png", width = 1200, height = 800)
plot(area_hclust, main = "地區分群樹狀圖（曼哈頓距離）")
dev.off()

# Alternative: 計算地區間的bray-curtis距離
area_dist <- vegdist(area_matrix, method = "bray")
area_hclust <- hclust(area_dist, method = "ward.D2")
png("Plant/area_hclust_braycurtis.png", width = 1200, height = 800)
plot(area_hclust, main = "地區分群樹狀圖（Bray-Curtis距離）")
dev.off()

# Alternative: 計算地區間的歐氏距離
area_dist <- dist(area_matrix, method = "euclidean")
area_hclust <- hclust(area_dist, method = "ward.D2")
png("Plant/area_hclust_euclidean.png", width = 1200, height = 800)
plot(area_hclust, main = "地區分群樹狀圖（歐氏距離）")
dev.off()

# 指定群數
area_clusters <- cutree(area_hclust, k = 8)

# 加入群組資訊
area_cluster_df <- data.frame(
  大區 = rownames(as.matrix(mat_standardized2)), 
  Cluster = area_clusters
)

# 標註每群包含的地區
cluster_areas <- area_cluster_df %>%
  group_by(Cluster) %>%
  summarise(Areas = paste(大區, collapse = ", "), .groups = "drop")
print(cluster_areas)

# 計算每群的 DBA 標準化總和
std_clustered_areas <- standardized_area_species %>%
  left_join(area_cluster_df, by = "大區") %>%
  group_by(Cluster, 樹種) %>%
  summarise(DBA_cluster_sum = sum(DBA_standardized), .groups = "drop")

# 轉換為矩陣格式
std_area_cluster_matrix <- std_clustered_areas %>%
  pivot_wider(names_from = 樹種, 
              values_from = DBA_cluster_sum, 
              values_fill = 0)
std_area_cluster_matrix2 <- as.data.frame(std_area_cluster_matrix)
rownames(std_area_cluster_matrix2) <- std_area_cluster_matrix2$Cluster
std_area_cluster_matrix2$Cluster <- NULL

# 畫標準化熱力圖
png("Plant/PlantCoverage_standarized_heatmap_area_cluster_manhattan.png",
    width = 4000, height = 3000)
# 確保 annotation_row 的 rownames 和矩陣的 rownames 順序和數量一致
cluster_names_in_area_matrix <- rownames(std_area_cluster_matrix2)
annotation_row_df <- cluster_areas %>%
  filter(Cluster %in% cluster_names_in_area_matrix) %>%
  arrange(match(Cluster, cluster_names_in_area_matrix))
annotation_row <- data.frame(
  Areas = annotation_row_df$Areas, 
  row.names = annotation_row_df$Cluster
)

pheatmap(as.matrix(std_area_cluster_matrix2),
         main = "地區群組 x 樹種 DBA 標準化熱力圖（曼哈頓距離）",
         clustering_method = "ward.D2",
         clustering_distance_cols = "manhattan",
         clustering_distance_rows = "manhattan",
         fontsize_row = 8,
         fontsize_col = 8,
         annotation_row = annotation_row)
dev.off()

# 10.3 樹種跟地區都分群的熱力圖
# 轉換為矩陣格式，行為地區，列為樹種
combined_matrix <- standardized_area_species %>%
    dplyr::select(大區, 樹種, DBA_standardized) %>%
    pivot_wider(names_from = 樹種, values_from = DBA_standardized, values_fill = 0)
combined_matrix2 <- as.data.frame(combined_matrix)
rownames(combined_matrix2) <- combined_matrix2$大區
combined_matrix2$大區 <- NULL

# 計算地區間的曼哈頓距離
combined_area_dist <- dist(combined_matrix2, method = "manhattan")
# 階層式分群
combined_area_hclust <- hclust(combined_area_dist, method = "ward.D2")
png("Plant/combined_area_hclust_manhattan.png", width = 1200, height = 800)
plot(combined_area_hclust, main = "地區分群樹狀圖（合併）")
dev.off()
# 指定群數
combined_area_clusters <- cutree(combined_area_hclust, k = 8)
# 加入群組資訊
combined_area_cluster_df <- data.frame(
    大區 = rownames(combined_matrix2), 
    Cluster = combined_area_clusters
)
# 標註每群包含的地區
combined_cluster_areas <- combined_area_cluster_df %>%
    group_by(Cluster) %>%
    summarise(Areas = paste(大區, collapse = ", "), .groups = "drop")
print(combined_cluster_areas)

# 計算樹種間的曼哈頓距離
combined_species_matrix <- t(as.matrix(combined_matrix2))
combined_species_dist <- dist(combined_species_matrix, method = "manhattan")
# 階層式分群
combined_species_hclust <- hclust(combined_species_dist, method = "ward.D2")
png("Plant/combined_species_hclust_manhattan.png", width = 1200, height = 800)
plot(combined_species_hclust, main = "樹種分群樹狀圖（合併）")
dev.off()
# 指定群數
combined_species_clusters <- cutree(combined_species_hclust, k = 14)
# 加入群組資訊
combined_species_cluster_df <- data.frame(樹種 = rownames(combined_species_matrix), Cluster = combined_species_clusters)
# 標註每群包含的樹種
combined_cluster_species <- combined_species_cluster_df %>%
    group_by(Cluster) %>%
    summarise(Species = paste(樹種, collapse = ", "), .groups = "drop")
print(combined_cluster_species)

# 計算地區群組 x 樹種群組的 DBA 標準化總和
combined_std_clustered <- standardized_area_species %>%
    left_join(combined_area_cluster_df, by = "大區") %>%
    left_join(combined_species_cluster_df, by = "樹種") %>%
    group_by(Cluster.x, Cluster.y) %>%
    summarise(DBA_cluster_sum = sum(DBA_standardized), .groups = "drop")

# 轉換為矩陣格式
combined_std_cluster_matrix <- combined_std_clustered %>%
    pivot_wider(names_from = Cluster.y, 
                values_from = DBA_cluster_sum, 
                values_fill = 0)
combined_std_cluster_matrix2 <- as.data.frame(combined_std_cluster_matrix)
rownames(combined_std_cluster_matrix2) <- combined_std_cluster_matrix2$Cluster.x
combined_std_cluster_matrix2$Cluster.x <- NULL

# 畫標準化熱力圖
png("Plant/PlantCoverage_standarized_heatmap_combined_cluster_manhattan.png",
    width = 4000, height = 3000)
# 準備 annotation
area_cluster_names_in_matrix <- rownames(combined_std_cluster_matrix2)
species_cluster_names_in_matrix <- colnames(combined_std_cluster_matrix2)

annotation_row_df <- combined_cluster_areas %>%
    filter(Cluster %in% area_cluster_names_in_matrix) %>%
    arrange(match(Cluster, area_cluster_names_in_matrix))
annotation_row <- data.frame(
    Areas = annotation_row_df$Areas,
    row.names = annotation_row_df$Cluster
)

annotation_col_df <- combined_cluster_species %>%
    filter(Cluster %in% species_cluster_names_in_matrix) %>%
    arrange(match(Cluster, species_cluster_names_in_matrix))
annotation_col <- data.frame(
    Species = annotation_col_df$Species,
    row.names = annotation_col_df$Cluster
)

pheatmap(as.matrix(combined_std_cluster_matrix2),
         main = "地區群組 x 樹種群組 DBA 標準化熱力圖（曼哈頓距離）",
         clustering_method = "ward.D2",
         clustering_distance_cols = "manhattan",
         clustering_distance_rows = "manhattan",
         fontsize_row = 8,
         fontsize_col = 8,
         annotation_row = annotation_row,
         annotation_col = annotation_col)
dev.off()

# 10.4 把熱力圖裡的數值改成Standarized Gapped Histogram的區間類別
set.seed(123)
standardized_data <- standardized_area_species$DBA_standardized
h <- hclust(dist(standardized_data))
clus <- 10 # 與 Gapped Histogram 相同的分群數
cl <- cutree(h, clus)

# 將分群類別加回資料框
standardized_area_species_with_gaps <- standardized_area_species %>%
    mutate(Gap_Cluster = cl)

# 建立熱力圖矩陣，值為區間類別
gap_matrix <- standardized_area_species_with_gaps %>%
    dplyr::select(大區, 樹種, Gap_Cluster) %>%
    pivot_wider(
        names_from = 樹種,
        values_from = Gap_Cluster,
        values_fill = 0 # 若無該樹種，以0填充
    )

gap_matrix2 <- as.data.frame(gap_matrix)
rownames(gap_matrix2) <- gap_matrix2$大區
gap_matrix2$大區 <- NULL

# 畫出熱力圖
png("Plant/PlantCoverage_standarized_heatmap_gap_clusters.png", width = 2000, height = 1500)
pheatmap(as.matrix(gap_matrix2),
                 main = "地區 x 樹種 Gapped Histogram 區間類別熱力圖",
                 clustering_method = "ward.D2",
                 clustering_distance_cols = "manhattan",
                 clustering_distance_rows = "manhattan",
                 fontsize_row = 8,
                 fontsize_col = 8,
                 display_numbers = TRUE, # 顯示數字
                 number_format = "%.0f", # 顯示整數
                 fontsize_number = 5)
dev.off()

# 11. 各大區的Entropy Histogram
# 定義entropy函數（參考WolfDice的做法）
my.entropy = function(x){p=x/sum(x);sum(-log(p)*p,na.rm=TRUE)}

# 計算每個大區的樹種組成entropy
areas <- unique(df$大區)
area_entropies <- data.frame()

for(area in areas) {
    # 取得該大區的樹種DBA分布
    area_data <- df %>% 
        filter(大區 == area) %>%
        group_by(樹種) %>%
        summarise(DBA_sum = sum(DBA), .groups = "drop")
    
    # 計算entropy
    entropy_value <- my.entropy(area_data$DBA_sum)
    area_entropies <- rbind(area_entropies, data.frame(大區 = area, Entropy = entropy_value))
}

print("各大區的Entropy值：")
print(area_entropies)

# 為每個大區生成null distribution並計算entropy histogram
set.seed(123)
n_simulations <- 10000
all_area_entropy_data <- list()

for(area in areas) {
    # 取得該大區的樹種DBA分布
    area_data <- df %>% 
        filter(大區 == area) %>%
        group_by(樹種) %>%
        summarise(DBA_sum = sum(DBA), .groups = "drop")
    
    # 計算實際的樹種比例
    actual_prob <- area_data$DBA_sum / sum(area_data$DBA_sum)
    n_species <- length(actual_prob)
    total_count <- sum(area_data$DBA_sum)
    
    # 生成null distribution（均勻分布）
    null_prob <- rep(1/n_species, n_species)
    null_samples <- replicate(n_simulations, {
        sample_counts <- rmultinom(1, total_count, null_prob)
        my.entropy(sample_counts)
    })
    
    # 生成alternative distribution（基於實際觀察）
    alt_samples <- replicate(n_simulations, {
        sample_counts <- rmultinom(1, total_count, actual_prob)
        my.entropy(sample_counts)
    })
    
    # 儲存資料
    all_area_entropy_data[[area]] <- list(
        null_entropy = null_samples,
        alt_entropy = alt_samples,
        actual_entropy = my.entropy(area_data$DBA_sum)
    )
}

# 為每個大區繪製entropy histogram
png("Plant/entropy_comparison.png", width = 4000, height = 3000)
par(mfrow = c(5, 4), mar = c(4, 4, 3, 1))

for(area in areas) {
    null_entropy <- all_area_entropy_data[[area]]$null_entropy
    alt_entropy <- all_area_entropy_data[[area]]$alt_entropy
    actual_entropy <- all_area_entropy_data[[area]]$actual_entropy
    
    # 設定圖形範圍
    xlim2 <- range(c(null_entropy, alt_entropy, actual_entropy))
    ax <- pretty(xlim2, n = 11)
    
    # 設定顏色
    col_nul <- adjustcolor("blue", alpha.f = 0.5)
    col_alt <- adjustcolor("red", alpha.f = 0.5)
    
    # 計算histogram
    hgNul <- hist(null_entropy, breaks = ax, plot = FALSE)
    hgAlt <- hist(alt_entropy, breaks = ax, plot = FALSE)
    
    # 繪製histogram
    plot(hgNul, col = col_nul, main = paste("大區:", area), 
         xlab = "Entropy", ylab = "Frequency")
    plot(hgAlt, col = col_alt, add = TRUE)
    
    # 加入實際entropy的垂直線
    abline(v = actual_entropy, col = "black", lwd = 2, lty = 2)
    
    # The legend will be added after the loop.
}

# Create an empty plot for the legend
plot.new()
# Define colors again for the legend, as they were defined inside the loop
col_nul <- adjustcolor("blue", alpha.f = 0.5)
col_alt <- adjustcolor("red", alpha.f = 0.5)
# Add the legend to the new empty plot area
legend("center", 
       legend = c("Null (均勻)", "Alternative (觀察)", "實際值"),
       col = c(col_nul, col_alt, "black"),
       lty = c(1, 1, 2),
       lwd = c(10, 10, 2),
       cex = 1.5, # Increased size for better readability in its own panel
       bty = "n") # No box around the legend

dev.off()