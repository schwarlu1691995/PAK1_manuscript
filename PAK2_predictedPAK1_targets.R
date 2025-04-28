

lted_stats.pak2 = LTED_stats[pak2_targets_identified,]

tamr_stats.pak2 = TAMR_stats[pak2_targets_identified,]
pak2_merged = merge(lted_stats.pak2, tamr_stats.pak2, by = 0, all = T)


pak2_merged = na.omit(pak2_merged)
rownames(pak2_merged) = pak2_merged$Row.names

plot(pak2_merged$t.x, pak2_merged$t.y, pch = 16, main = "T-stats of PAK2 target phosphosites", xlim=c(-30,5),ylim = c(-30,5),xlab = "LTED", ylab = "TAMR")
textplot(pak2_merged$t.x, pak2_merged$t.y, rownames(pak2_merged), new = F, show.lines = T)
abline(v=0, lty = 2)
abline(h=0, lty = 2)



lted_stats.pak2 = LTED_stats[grep("PAK2", rownames(LTED_stats)),]

tamr_stats.pak2 = TAMR_stats[grep("PAK2", rownames(TAMR_stats)),]

pak2_merged = merge(lted_stats.pak2, tamr_stats.pak2, by = 0, all = T)


pak2_merged = na.omit(pak2_merged)
rownames(pak2_merged) = pak2_merged$Row.names

plot(pak2_merged$t.x, pak2_merged$t.y, pch = 16, main = "T-stats of PAK2 target phosphosites", xlim=c(-30,5),ylim = c(-30,5),xlab = "LTED", ylab = "TAMR")
textplot(pak2_merged$t.x, pak2_merged$t.y, rownames(pak2_merged), new = F, show.lines = T)
abline(v=0, lty = 2)
abline(h=0, lty = 2)



# predicted PAK1 target sites:

setwd("K:/Ergebnisse/LS_testing/1_PhD/paper/EGF_dynamics/additional_plots/kinase_assignment/")
nature_paper_result = read.csv("nature_paper_result.csv", header = T, stringsAsFactors = F, sep = ",")

nature_pak1_sites = nature_paper_result[nature_paper_result$PAK1_percentile>98 & nature_paper_result$PAK1_rank <6,]
# 288 potential pak1 target sites
nature_pak1_sites$gene_site = paste0(nature_pak1_sites$Gene, "_", nature_pak1_sites$Phosphosite)



lted_stats.pak1 = LTED_stats[rownames(LTED_stats) %in% nature_pak1_sites$gene_site,]

tamr_stats.pak1 = TAMR_stats[rownames(TAMR_stats) %in% nature_pak1_sites$gene_site,]

pak1_merged = merge(lted_stats.pak1, tamr_stats.pak1, by = 0, all = T)
rownames(pak1_merged) = pak1_merged$Row.names
pak1_merged = na.omit(pak1_merged)

barplot(pak1_merged$t.x[order(pak1_merged$t.x, decreasing = T)], las = 2, names.arg = rownames(pak1_merged)[order(pak1_merged$t.x, decreasing = T)], main = "T-stats of potential PAK1 target sites in LTED")
barplot(pak1_merged$t.y[order(pak1_merged$t.y, decreasing = T)], las = 2, names.arg = rownames(pak1_merged)[order(pak1_merged$t.y, decreasing = T)], main = "T-stats of potential PAK1 target sites in TAMR")





plot(pak1_merged$t.x, pak1_merged$t.y, pch = 16, main = "T-stats of predicted PAK1 target phosphosites",xlim = c(-13,5), ylim = c(-13,5), xlab = "LTED", ylab = "TAMR")
textplot(pak1_merged$t.x[abs(pak1_merged$t.x) > 2 | abs(pak1_merged$t.y) > 2], pak1_merged$t.y[abs(pak1_merged$t.x) > 2 | abs(pak1_merged$t.y) > 2], rownames(pak1_merged)[abs(pak1_merged$t.x) > 2 | abs(pak1_merged$t.y) > 2], new = F, show.lines = T)
abline(v=-2, lty = 2, col = "blue")
abline(h=-2, lty = 2, col = "blue")
abline(v=2, lty = 2, col = "red")
abline(h=2, lty = 2, col = "red")








