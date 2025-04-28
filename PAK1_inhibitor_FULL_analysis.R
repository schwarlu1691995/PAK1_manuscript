# Exp61
# MCF7 WT, LTED, TAMR
# PAK1-inhibitor treatment, 1h
# full + phospho MS
# SN17, directDIA, "final_settings", no isoforms


setwd("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp61_MCF7_LTED_TAMR_WT_PAK1-inh_treatment_phosphoMS/")


PAK1_inh_full = read.csv("20240826_143012_20240826_Exp61_PAK1inh_FULL_Report.csv")
table(duplicated(PAK1_inh_full$PG.Genes))
PAK1_inh_full$PG.Genes[duplicated(PAK1_inh_full$PG.Genes)]
PAK1_inh_full = PAK1_inh_full[!duplicated(PAK1_inh_full$PG.Genes),]

rownames(PAK1_inh_full) = PAK1_inh_full$PG.Genes
PAK1_inh_full = PAK1_inh_full[,-c(1:3)]

colnames(PAK1_inh_full) = paste0(rep(c("LTED", "TAMR", "WT"), each = 6),"_", rep(c("DMSO", "inh"), c(3,3)),"_", rep(1:3, 6))

par(mar=c(10,5,5,5))
barplot(apply(PAK1_inh_full, 2, function(x) sum(!is.na(x))), ylim = c(0,8000), las = 2, main = "identified protein groups")

boxplot(log2(PAK1_inh_full), las = 2, main = "log2 protein intensities")

library(pheatmap)
pheatmap(na.omit(PAK1_inh_full), scale = "row")

library(limma)
library(tidyverse)

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

PAK1_inh_full[is.nan.data.frame(PAK1_inh_full)] <- NA

#vsn
PAK1_inh_full_norm = limma::normalizeVSN(PAK1_inh_full)
boxplot((PAK1_inh_full_norm), las = 2, main = "log2 protein intensities")

column_annotation = data.frame(cell_line = rep(c("LTED", "TAMR", "WT"), each = 6), treatment = rep(c("DMSO", "inh", "DMSO", "inh", "DMSO", "inh"), c(3,3,3,3,3,3)), replicate = rep(1:3, 6))
rownames(column_annotation) = colnames(PAK1_inh_full_norm)

pheatmap(na.omit(PAK1_inh_full_norm), scale = "row", annotation_col = column_annotation, show_rownames = F, color = colorRampPalette(c("#5390c1", "white", "#ea513f"))(100), breaks = seq(-3,3,6/100))

# limma stats of inhibitor vs. DMSO:

# unpaired comparison

data_log2 <- PAK1_inh_full_norm
LTED_columns = grep("LTED", colnames(PAK1_inh_full_norm))
TAMR_columns = grep("TAMR", colnames(PAK1_inh_full_norm))
WT_columns = grep("WT", colnames(PAK1_inh_full_norm))

# get design matrix:
condition <- factor(sapply(strsplit(colnames(PAK1_inh_full_norm), split = "_"), function(x) paste0(x[1], "_", x[2])))
table(condition)

design <- model.matrix(~0+ condition)
colnames(design) = gsub("condition","",colnames(design))






# contrast matrix:
condition1 = "WT_inh"
condition2 = "WT_DMSO"

contrast_name <- paste0(condition1, '-', condition2)
contrast_name

# filter for 6/6 values present per condition:
data_filtered <- data_log2[apply(data_log2[,WT_columns],1,function(x) sum(!is.na(x))) == 6, ]
nrow(data_filtered)#6671, 7036, 7255 proteins left

# fit limma linear model:
array_weights <- NULL

fit <- limma::lmFit(data_filtered, design=design, weights=array_weights)

contrast_matrix <- limma::makeContrasts(contrasts=eval(contrast_name),
                                        levels=design)

fit_contrast <- limma::contrasts.fit(fit, contrast_matrix)
# statistics with eBayes method:
stats <- eBayes(fit_contrast, trend=T,robust = T)

# get result in different format:

stats_output <- limma::topTable(stats, coef=1, number=nrow(stats), adjust.method="BH", sort.by="none") %>% 
  drop_na()

write.csv(stats_output, file = "20240826_WT_inh_vs_DMSO_limma_stats.csv", row.names = T)



# plot volcano:
#install.packages(c("wordcloud","tm"),repos="http://cran.r-project.org")
library(wordcloud)


plot(stats_output$logFC[stats_output$logFC > 0], -log10(stats_output$adj.P.Val[stats_output$logFC > 0]), pch = 16, col = ifelse(stats_output$logFC[stats_output$logFC > 0] > 1 & stats_output$adj.P.Val[stats_output$logFC > 0] < 0.05,'#cc006680','#00000030'), main = paste0(contrast_name," \nfull proteome (limma eBayes statistics)"), xlab = "log2 FC", ylab = "adjusted -log10 p-value", xlim = range(stats_output$logFC)*1.2, ylim = range(-log10(stats_output$adj.P.Val))*1.2)
points(stats_output$logFC[stats_output$logFC < 0], -log10(stats_output$adj.P.Val[stats_output$logFC < 0]), pch = 16, col = ifelse(stats_output$logFC[stats_output$logFC < 0] < -1 & stats_output$adj.P.Val[stats_output$logFC < 0] < 0.05,'#0066cc80','#00000030'))
textplot(stats_output$logFC[abs(stats_output$logFC) > 1 & stats_output$adj.P.Val < 0.05], -log10(stats_output$adj.P.Val[abs(stats_output$logFC) > 1 & stats_output$adj.P.Val < 0.05])+0.05, rownames(stats_output)[abs(stats_output$logFC) > 1 & stats_output$adj.P.Val < 0.05], new = F, show.lines = F)
legend("topright", bty = "n", pch = 16, col = "#cc006680", legend = paste0("up in\n", condition1))
legend("topleft", bty = "n", pch = 16, col = "#0066cc80", legend = paste0("up in\n", condition2))




plot(stats_output$logFC[stats_output$logFC > 0], -log10(stats_output$P.Value[stats_output$logFC > 0]), pch = 16, col = ifelse(stats_output$logFC[stats_output$logFC > 0] > 1 & stats_output$P.Value[stats_output$logFC > 0] < 0.05,'#cc006680','#00000030'), main = paste0(contrast_name," \nfull proteome (limma eBayes statistics)"), xlab = "log2 FC", ylab = "-log10 p-value (non adjusted!)", xlim = range(stats_output$logFC)*1.2, ylim = range(-log10(stats_output$P.Value))*1.2)
points(stats_output$logFC[stats_output$logFC < 0], -log10(stats_output$P.Value[stats_output$logFC < 0]), pch = 16, col = ifelse(stats_output$logFC[stats_output$logFC < 0] < -1 & stats_output$P.Value[stats_output$logFC < 0] < 0.05,'#0066cc80','#00000030'))
textplot(stats_output$logFC[abs(stats_output$logFC) > 1 & stats_output$P.Value < 0.05], -log10(stats_output$P.Value[abs(stats_output$logFC) > 1 & stats_output$P.Value < 0.05])+0.05, rownames(stats_output)[abs(stats_output$logFC) > 1 & stats_output$P.Value < 0.05], new = F, show.lines = F)
legend("topright", bty = "n", pch = 16, col = "#cc006680", legend = paste0("up in\n", condition1))
legend("topleft", bty = "n", pch = 16, col = "#0066cc80", legend = paste0("up in\n", condition2))


# no significant full proteome changes!


# Reactome pathways:
library(ReactomePA)

# map gene names to Entrez ID
library('org.Hs.eg.db')
stats_output = read.csv("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp61_MCF7_LTED_TAMR_WT_PAK1-inh_treatment_phosphoMS/20240826_TAMR_inh_vs_DMSO_limma_stats.csv", row.names = 1)

stats_output$entrez_ID = mapIds(org.Hs.eg.db, rownames(stats_output), 'ENTREZID', 'SYMBOL')

require(DOSE)
library(stats)

# GSEA:
gene_list = stats_output$logFC
names(gene_list) = stats_output$entrez_ID

gene_list_sorted = gene_list[order(gene_list, decreasing = T)]
gene_list_sorted = gene_list_sorted[!duplicated(names(gene_list_sorted))]
gene_list_sorted = na.omit(gene_list_sorted)
length(gene_list_sorted)


y <- gsePathway(gene_list_sorted, 
                minGSSize=120, pvalueCutoff=0.2, 
                pAdjustMethod="BH", verbose=FALSE)


#head(summary(x))
dotplot(y, showCategory=20, x = "NES", color = "p.adjust", title = "WT PAK1-inhibitor vs. DMSO (1h), Reactome GSEA")

# no terms enriched! 








