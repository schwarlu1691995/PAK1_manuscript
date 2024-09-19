# Exp08, EGF-stimulated time-course
# data filtered & normalized by EV

load("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/from_Efstathios/full_normalized/TimeCourse_FullProteome_Normalized_Filt_Proc_17052024.rda")


full_timeCourse = data.frame(se_final@assays@data@listData$Norm_Filt)
full_timeCourse$gene = sapply(strsplit(rownames(full_timeCourse), split = "_"), function(x) x[2])



PAK1_targets = read.csv("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/normal_growth_condition/PHOSPHO/Omnipath_KSN_NGC_exp_LTED_vs_WT_filt_PAK1_targets_Extract_TargetProtein_SiteID.tsv", sep = "\t", header = T, stringsAsFactors = F)
TC_full_PAK1_targets = full_timeCourse[full_timeCourse$gene %in% PAK1_targets$Target_Protein,]

pheatmap(TC_full_PAK1_targets[,1:63], scale = "row", cluster_rows = T, cluster_cols = F, breaks = seq(-3,3,6/100))

# calculate means:
full_timeCourse_means = data.frame(matrix(NA, nrow=nrow(full_timeCourse), ncol = 63/3))
colnames(full_timeCourse_means) = paste0(rep(c("LTED_", "TAMR_", "WT_"), each = 7), rep(c(0,2,5,10,20,60,120), 3))
rownames(full_timeCourse_means) = rownames(full_timeCourse)

subtract = rep(1:21, each = 3)

for (i in seq(1,63,3)){
  full_timeCourse_means[,i-(i-subtract[i])] <- apply(full_timeCourse[,i:(i+2)],1,mean,na.rm=T)
}
full_timeCourse_means$gene = sapply(strsplit(rownames(full_timeCourse_means), split = "_"), function(x) x[2])


# calculate sd:
full_timeCourse_sd = data.frame(matrix(NA, nrow=nrow(full_timeCourse), ncol = 63/3))
colnames(full_timeCourse_sd) = paste0(rep(c("LTED_", "TAMR_", "WT_"), each = 7), rep(c(0,2,5,10,20,60,120), 3))
rownames(full_timeCourse_sd) = rownames(full_timeCourse)

subtract = rep(1:21, each = 3)

for (i in seq(1,63,3)){
  full_timeCourse_sd[,i-(i-subtract[i])] <- apply(full_timeCourse[,i:(i+2)],1,sd,na.rm=T)
}
full_timeCourse_sd$gene = sapply(strsplit(rownames(full_timeCourse_sd), split = "_"), function(x) x[2])


# plot PAK1_targets:
TC_full_means_PAK1_targets = full_timeCourse_means[full_timeCourse_means$gene %in% PAK1_targets$Target_Protein,]
TC_full_sd_PAK1_targets = full_timeCourse_sd[full_timeCourse_sd$gene %in% PAK1_targets$Target_Protein,]

TC_full_PAK1_targets_means_sd = merge(TC_full_means_PAK1_targets, TC_full_sd_PAK1_targets, by = "gene", all = T)

pheatmap(TC_full_means_PAK1_targets[,1:21], scale = "row", cluster_rows = T, cluster_cols = F, breaks = seq(-3,3,6/100),color = colorRampPalette(c("#5390c1", "white", "#ea513f"))(100),)

TC_full_PAK1_targets_means_sd$FC_LTED_WT = apply(TC_full_PAK1_targets_means_sd[,c(2:8, 16:22)], 1, function(x) mean(x[1:7]-x[8:14], na.rm=T))
rownames(TC_full_PAK1_targets_means_sd) = TC_full_PAK1_targets_means_sd$gene
TC_full_PAK1_targets_means_sd$FC_TAMR_WT = apply(TC_full_PAK1_targets_means_sd[,c(9:15, 16:22)], 1, function(x) mean(x[1:7]-x[8:14], na.rm=T))
TC_full_PAK1_targets_means_sd$sd_LTED_WT = apply(TC_full_PAK1_targets_means_sd[,c(2:8, 16:22)], 1, function(x) sd(x[1:7]-x[8:14], na.rm=T))
TC_full_PAK1_targets_means_sd$sd_TAMR_WT = apply(TC_full_PAK1_targets_means_sd[,c(9:15, 16:22)], 1, function(x) sd(x[1:7]-x[8:14], na.rm=T))

TC_full_PAK1_targets_means_sd = TC_full_PAK1_targets_means_sd[-9,]

par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
plot(TC_full_PAK1_targets_means_sd$FC_LTED_WT, TC_full_PAK1_targets_means_sd$FC_TAMR_WT,ylim = c(-1,4), xlim = c(-1,4), main = "known PAK1 target proteins\nlog2 FC compared to WT", xlab = "LTED", ylab = "TAMR", pch = 16, cex = 1.5)
abline(h=0, lty = 2)
abline(v=0, lty = 2)
segments(x0 = TC_full_PAK1_targets_means_sd$FC_LTED_WT-TC_full_PAK1_targets_means_sd$sd_LTED_WT/sqrt(3), x1 = TC_full_PAK1_targets_means_sd$FC_LTED_WT+TC_full_PAK1_targets_means_sd$sd_LTED_WT/sqrt(3), y0 = TC_full_PAK1_targets_means_sd$FC_TAMR_WT, col = "#876cad", lwd = 2)
segments(y0 = TC_full_PAK1_targets_means_sd$FC_TAMR_WT-TC_full_PAK1_targets_means_sd$sd_TAMR_WT/sqrt(3), y1 = TC_full_PAK1_targets_means_sd$FC_TAMR_WT+TC_full_PAK1_targets_means_sd$sd_TAMR_WT/sqrt(3), x0 = TC_full_PAK1_targets_means_sd$FC_LTED_WT, col = "#d76d17", lwd = 2)
library(wordcloud)
textplot(TC_full_PAK1_targets_means_sd$FC_LTED_WT, TC_full_PAK1_targets_means_sd$FC_TAMR_WT, rownames(TC_full_PAK1_targets_means_sd), new = F, show.lines = F)


# phospho PAK1 targets
library(PhosR)
load("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/from_Efstathios/phospho_normalized/RE_Analysis_TimeCourse_Phospho_Norm_Filt_Proc_09102023.RData")
phospho.mat = ppe_filt@assays@data$Normalization
phospho.mat = data.frame(phospho.mat)

phospho.mat$gene_site = sapply(strsplit(rownames(phospho.mat), split = ";", fixed = T), function(x) paste0(x[2], "_", x[3]))
table(duplicated(phospho.mat$gene_site))
phospho.mat = phospho.mat[!duplicated(phospho.mat$gene_site),]
rownames(phospho.mat) = phospho.mat$gene_site
phospho.mat = phospho.mat[,-64]


# calculate means:
phospho_timeCourse_means = data.frame(matrix(NA, nrow=nrow(phospho.mat), ncol = 63/3))
colnames(phospho_timeCourse_means) = paste0(rep(c("LTED_", "TAMR_", "WT_"), each = 7), rep(c(0,2,5,10,20,60,120), 3))
rownames(phospho_timeCourse_means) = rownames(phospho.mat)

subtract = rep(1:21, each = 3)

for (i in seq(1,63,3)){
  phospho_timeCourse_means[,i-(i-subtract[i])] <- apply(phospho.mat[,i:(i+2)],1,mean,na.rm=T)
}


# calculate sd:
phospho_timeCourse_sd = data.frame(matrix(NA, nrow=nrow(phospho.mat), ncol = 63/3))
colnames(phospho_timeCourse_sd) = paste0(rep(c("LTED_", "TAMR_", "WT_"), each = 7), rep(c(0,2,5,10,20,60,120), 3))
rownames(phospho_timeCourse_sd) = rownames(phospho.mat)

subtract = rep(1:21, each = 3)

for (i in seq(1,63,3)){
  phospho_timeCourse_sd[,i-(i-subtract[i])] <- apply(phospho.mat[,i:(i+2)],1,sd,na.rm=T)
}

# PAK1 targets:
PAK1_targets$gene_site = paste0(PAK1_targets$Target_Protein, "_", PAK1_targets$Phospho_Site_ID)
phospho_timeCourse_means = phospho_timeCourse_means[rownames(phospho_timeCourse_means) %in% PAK1_targets$gene_site,]
phospho_timeCourse_sd = phospho_timeCourse_sd[rownames(phospho_timeCourse_sd) %in% PAK1_targets$gene_site,]

phospho_timeCourse_means_sd = merge(phospho_timeCourse_means, phospho_timeCourse_sd, by = 0, all = T)
rownames(phospho_timeCourse_means_sd) = phospho_timeCourse_means_sd$Row.names
phospho_timeCourse_means_sd = phospho_timeCourse_means_sd[,-1]

pheatmap(phospho_timeCourse_means_sd[,1:21], scale = "row", cluster_rows = T, cluster_cols = F, breaks = seq(-3,3,6/100),color = colorRampPalette(c("#5390c1", "white", "#ea513f"))(100),)

phospho_timeCourse_means_sd$FC_LTED_WT = apply(phospho_timeCourse_means_sd[,c(1:7, 15:21)], 1, function(x) mean(x[1:7]-x[8:14], na.rm=T))

phospho_timeCourse_means_sd$FC_TAMR_WT = apply(phospho_timeCourse_means_sd[,c(8:14, 15:21)], 1, function(x) mean(x[1:7]-x[8:14], na.rm=T))
phospho_timeCourse_means_sd$sd_LTED_WT = apply(phospho_timeCourse_means_sd[,c(1:7, 15:21)], 1, function(x) sd(x[1:7]-x[8:14], na.rm=T))
phospho_timeCourse_means_sd$sd_TAMR_WT = apply(phospho_timeCourse_means_sd[,c(8:14, 15:21)], 1, function(x) sd(x[1:7]-x[8:14], na.rm=T))



par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
plot(phospho_timeCourse_means_sd$FC_LTED_WT, phospho_timeCourse_means_sd$FC_TAMR_WT,ylim = c(-3,7), xlim = c(-3,7), main = "known PAK1 target phosphosites\nlog2 FC compared to WT", xlab = "LTED", ylab = "TAMR", pch = 16, cex = 1.5)
abline(h=0, lty = 2)
abline(v=0, lty = 2)
segments(x0 = phospho_timeCourse_means_sd$FC_LTED_WT-phospho_timeCourse_means_sd$sd_LTED_WT/sqrt(3), x1 = phospho_timeCourse_means_sd$FC_LTED_WT+phospho_timeCourse_means_sd$sd_LTED_WT/sqrt(3), y0 = phospho_timeCourse_means_sd$FC_TAMR_WT, col = "#876cad", lwd = 2)
segments(y0 = phospho_timeCourse_means_sd$FC_TAMR_WT-phospho_timeCourse_means_sd$sd_TAMR_WT/sqrt(3), y1 = phospho_timeCourse_means_sd$FC_TAMR_WT+phospho_timeCourse_means_sd$sd_TAMR_WT/sqrt(3), x0 = phospho_timeCourse_means_sd$FC_LTED_WT, col = "#d76d17", lwd = 2)
library(wordcloud)
textplot(phospho_timeCourse_means_sd$FC_LTED_WT[!is.na(phospho_timeCourse_means_sd$FC_LTED_WT) & !is.na(phospho_timeCourse_means_sd$FC_TAMR_WT)], phospho_timeCourse_means_sd$FC_TAMR_WT[!is.na(phospho_timeCourse_means_sd$FC_LTED_WT) & !is.na(phospho_timeCourse_means_sd$FC_TAMR_WT)], rownames(phospho_timeCourse_means_sd)[!is.na(phospho_timeCourse_means_sd$FC_LTED_WT) & !is.na(phospho_timeCourse_means_sd$FC_TAMR_WT)], new = F, show.lines = F)

# merge phospho and full:
phospho_timeCourse_means_sd$gene = sapply(strsplit(rownames(phospho_timeCourse_means_sd), "_"), function(x) x[1])

TC_PAK1_targets_full_phospho_merged = merge(phospho_timeCourse_means_sd, TC_full_PAK1_targets_means_sd, by.x = "gene", by.y = 0, all = T)
rownames(TC_PAK1_targets_full_phospho_merged) = rownames(phospho_timeCourse_means_sd)


# scatterplot full vs phospho
# 1. LTED:

par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
plot(TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.y, TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.x,ylim = c(-2,4.5), xlim = c(-2,4.5), main = "known PAK1 target proteins and phosphosites\nLTED compared to WT", xlab = "protein FC", ylab = "phosphosite FC", pch = 16, cex = 1.5)
abline(h=0, lty = 2)
abline(v=0, lty = 2)
segments(x0 = TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.y-(TC_PAK1_targets_full_phospho_merged$sd_LTED_WT.y/sqrt(3)), x1 = TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.y+TC_PAK1_targets_full_phospho_merged$sd_LTED_WT.y/sqrt(3), y0 = TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.x, col = "#876cad", lwd = 1)
segments(y0 = TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.x-TC_PAK1_targets_full_phospho_merged$sd_LTED_WT.x/sqrt(3), y1 = TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.x+TC_PAK1_targets_full_phospho_merged$sd_LTED_WT.x/sqrt(3), x0 = TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.y, col = "#876cad", lwd = 1)
library(wordcloud)
textplot(TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.y[!is.na(TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.y) & !is.na(TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.x)], TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.x[!is.na(TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.x) & !is.na(TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.y)], rownames(TC_PAK1_targets_full_phospho_merged)[!is.na(TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.x) & !is.na(TC_PAK1_targets_full_phospho_merged$FC_LTED_WT.y)], new = F, show.lines = F)

par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
plot(TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.y, TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.x,ylim = c(-2,6), xlim = c(-2,3), main = "known PAK1 target proteins and phosphosites\nTAMR compared to WT", xlab = "protein FC", ylab = "phosphosite FC", pch = 16, cex = 1.5)
abline(h=0, lty = 2)
abline(v=0, lty = 2)
segments(x0 = TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.y-(TC_PAK1_targets_full_phospho_merged$sd_TAMR_WT.y/sqrt(3)), x1 = TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.y+TC_PAK1_targets_full_phospho_merged$sd_TAMR_WT.y/sqrt(3), y0 = TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.x, col = "#d76d17", lwd = 1)
segments(y0 = TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.x-TC_PAK1_targets_full_phospho_merged$sd_TAMR_WT.x/sqrt(3), y1 = TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.x+TC_PAK1_targets_full_phospho_merged$sd_TAMR_WT.x/sqrt(3), x0 = TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.y, col = "#d76d17", lwd = 1)
library(wordcloud)
textplot(TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.y[!is.na(TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.y) & !is.na(TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.x)], TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.x[!is.na(TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.x) & !is.na(TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.y)], rownames(TC_PAK1_targets_full_phospho_merged)[!is.na(TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.x) & !is.na(TC_PAK1_targets_full_phospho_merged$FC_TAMR_WT.y)], new = F, show.lines = F)



# compare expression at time-point 0 (after starvation) with expression in normal growth condition:

TC_full_means_PAK1_targets_TP0 = TC_full_means_PAK1_targets[,seq(1,22,7)]

# Full proteome:
load("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/normal_growth_condition/FULL/RE_Analysis_NormalGrowth_FullProteome_Norm_Filt_Proc_09112023.RData")
normalGrowth_full = se_final@assays@data@listData$Norm_Filt



# calculate means:
full_normalGrowth_means = data.frame(matrix(NA, nrow=nrow(normalGrowth_full), ncol = 3))
colnames(full_normalGrowth_means) = c("LTED", "TAMR","WT")
rownames(full_normalGrowth_means) = rownames(normalGrowth_full)

subtract = rep(1:9, each = 3)

for (i in seq(1,9,3)){
  full_normalGrowth_means[,i-(i-subtract[i])] <- apply(normalGrowth_full[,i:(i+2)],1,mean,na.rm=T)
}
full_normalGrowth_means$gene = sapply(strsplit(rownames(full_normalGrowth_means), split = "_"), function(x) x[2])

# calculate sd:
full_normalGrowth_sd = data.frame(matrix(NA, nrow=nrow(normalGrowth_full), ncol = 3))
colnames(full_normalGrowth_sd) = c("LTED", "TAMR","WT")
rownames(full_normalGrowth_sd) = rownames(normalGrowth_full)

subtract = rep(1:9, each = 3)

for (i in seq(1,9,3)){
  full_normalGrowth_sd[,i-(i-subtract[i])] <- apply(normalGrowth_full[,i:(i+2)],1,sd,na.rm=T)
}
full_normalGrowth_sd$gene = sapply(strsplit(rownames(full_normalGrowth_sd), split = "_"), function(x) x[2])


normalGrowth_full_PAK1_targets_means_sd = merge(full_normalGrowth_means, full_normalGrowth_sd, by = "gene", all = T)


normal_TP0_PAK1Targets_merged = merge(normalGrowth_full_PAK1_targets_means_sd, TC_full_PAK1_targets_means_sd[,c(1,2,9,16,23,30,37)], by = "gene", all = T)
rownames(normal_TP0_PAK1Targets_merged) = normal_TP0_PAK1Targets_merged$gene

normal_TP0_PAK1Targets_merged=normal_TP0_PAK1Targets_merged[,-1]

colnames(normal_TP0_PAK1Targets_merged) =paste(rep(c("LTED_", "TAMR_", "WT_"),4),rep(c("normal","starved"), each = 6), rep(rep(c("mean", "sd"), each = 3),2))

normal_TP0_PAK1Targets_merged = normal_TP0_PAK1Targets_merged[rownames(normal_TP0_PAK1Targets_merged) %in% PAK1_targets$Target_Protein,]


WT_col = "#3fb6a9"
LTED_col = "#d64c9a"
TAMR_col = "#d86d16"


par(mfrow=c(2,2))
par(mar=c(10,5,5,5))

for (i in 1:nrow(normal_TP0_PAK1Targets_merged)){
  pg.bar = barplot(data.matrix(normal_TP0_PAK1Targets_merged[i,c(1,7,2,8,3,9)]), las = 2, main = rownames(normal_TP0_PAK1Targets_merged)[i])
  #segments(x0 = pg.bar, y0 = (normal_TP0_PAK1Targets_merged[i,c(1:3,7:9)]-(normal_TP0_PAK1Targets_merged[i,c(4:6,10:12)], y1 = (normal_TP0_PAK1Targets_merged[i,c(1:3,7:9)]+(normal_TP0_PAK1Targets_merged[i,c(4:6,10:12)])
}


normal_TP0_PAK1Targets_merged$FC_LTED_WT_normal = apply(normal_TP0_PAK1Targets_merged[,c(1,3)],1, function(x) x[1]-x[2])
normal_TP0_PAK1Targets_merged$FC_TAMR_WT_normal = apply(normal_TP0_PAK1Targets_merged[,c(2,3)],1, function(x) x[1]-x[2])
normal_TP0_PAK1Targets_merged$FC_LTED_WT_starved = apply(normal_TP0_PAK1Targets_merged[,c(7,9)],1, function(x) x[1]-x[2])
normal_TP0_PAK1Targets_merged$FC_TAMR_WT_starved = apply(normal_TP0_PAK1Targets_merged[,c(8,9)],1, function(x) x[1]-x[2])

pheatmap(normal_TP0_PAK1Targets_merged[,c(13,15,14,16)], cluster_rows = F, cluster_cols = F, scale = "none", breaks = seq(-2,2,4/100))

for (i in 1:nrow(normal_TP0_PAK1Targets_merged)){
  pg.bar = barplot(data.matrix(normal_TP0_PAK1Targets_merged[i,c(13,15,14,16)]), las = 2, main = rownames(normal_TP0_PAK1Targets_merged)[i])
  #segments(x0 = pg.bar, y0 = (normal_TP0_PAK1Targets_merged[i,c(1:3,7:9)]-(normal_TP0_PAK1Targets_merged[i,c(4:6,10:12)], y1 = (normal_TP0_PAK1Targets_merged[i,c(1:3,7:9)]+(normal_TP0_PAK1Targets_merged[i,c(4:6,10:12)])
}

