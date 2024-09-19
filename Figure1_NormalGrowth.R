# MCF7 time-course paper
# Plots for figure 1
# A: Volcano of normal growth condition, full proteome TAMR/LTED vs. WT
# Statistics from EV

setwd("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/normal_growth_condition/FULL/")

# Full proteome:
load("RE_Analysis_NormalGrowth_FullProteome_Norm_Filt_Proc_09112023.RData")
full = se_final@assays@data@listData$Norm_Filt


library(pheatmap)
full_cluster = pheatmap(na.omit(full[,c(7:9,1:6)]), scale = "row", cluster_rows = T,
                        cluster_cols = F,
                        colorRampPalette(c("#5390c1", "white", "#ea513f"))(100), show_rownames = F, treeheight_col =0, treeheight_row = 20)

gene_names_with_cluster = sort(cutree(full_cluster$tree_row, k=3))


LTED_cluster = names(gene_names_with_cluster[gene_names_with_cluster==1])
LTED_cluster = sapply(strsplit(LTED_cluster, split = "_"), function(x) x[2])
WT_cluster = names(gene_names_with_cluster[gene_names_with_cluster==3])
WT_cluster = sapply(strsplit(WT_cluster, split = "_"), function(x) x[2])
TAMR_cluster = names(gene_names_with_cluster[gene_names_with_cluster==2])
TAMR_cluster = sapply(strsplit(TAMR_cluster, split = "_"), function(x) x[2])


# statistics:
library(limma)
condition <- factor(sapply(strsplit(colnames(full), split = "_"), function(x) x[2]))
table(condition)

replicate = factor(sapply(strsplit(colnames(full), split = "_"), function(x) x[3]))
table(replicate)



design <- model.matrix(~0+ condition + replicate)
colnames(design) = c(levels(condition), levels(replicate)[-1])
nrow(design)
ncol(full)

# fit limma linear model:

fit <- limma::lmFit(full, design=design, weights=NULL)



# contrast matrix:
condition1 = "TAMR"
condition2 = "WT"

contrast_name <- paste0(condition1, '-', condition2)
contrast_name


contrast_matrix <- limma::makeContrasts(contrasts=eval(contrast_name),
                                        levels=design)

fit_contrast <- limma::contrasts.fit(fit, contrast_matrix)
# statistics with eBayes method:
stats <- eBayes(fit_contrast, trend=T,robust = T)

# get result in different format:
# # change name!!!
library(tidyverse)
stats_output <- limma::topTable(stats, coef=1, number=nrow(stats), adjust.method="BH", sort.by="none") %>% 
  drop_na()

write.csv(stats_output, file = "20240523_LTED_vs_WT_limma_stats.csv", row.names = T)

LTED_stats = read.csv("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/normal_growth_condition/FULL/20240523_LTED_vs_WT_limma_stats.csv", row.names = 1)
TAMR_stats = read.csv("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/normal_growth_condition/FULL/20240523_TAMR_vs_WT_limma_stats.csv", row.names = 1)


rownames(LTED_stats) = sapply(strsplit(rownames(LTED_stats), split = "_"), function(x) x[2])
rownames(TAMR_stats) = sapply(strsplit(rownames(TAMR_stats), split = "_"), function(x) x[2])
sum(LTED_stats$adj.P.Val<0.05) #4591 proteins, 60%
sum(LTED_stats$adj.P.Val<0.05)/nrow(LTED_stats)
sum(LTED_stats$adj.P.Val<0.05 & abs(LTED_stats$logFC) > 1) #559 proteins, 7.4%
sum(LTED_stats$adj.P.Val<0.05 & abs(LTED_stats$logFC) > 1)/nrow(LTED_stats)


sum(TAMR_stats$adj.P.Val<0.05) #3333 proteins, 44%
sum(TAMR_stats$adj.P.Val<0.05)/nrow(TAMR_stats)
sum(TAMR_stats$adj.P.Val<0.05 & abs(TAMR_stats$logFC) > 1) #328 proteins, 4.3%
sum(TAMR_stats$adj.P.Val<0.05 & abs(TAMR_stats$logFC) > 1)/nrow(TAMR_stats)

rownames(full) = sapply(strsplit(rownames(full), "_"), function(x) x[2])
full_sign = full[union(rownames(TAMR_stats)[(TAMR_stats$adj.P.Val<0.05 & abs(TAMR_stats$logFC) > 1)],rownames(LTED_stats)[(LTED_stats$adj.P.Val<0.05 & abs(LTED_stats$logFC) > 1)]),]

pheatmap(na.omit(full_sign[,c(7:9,1:6)]), scale = "row", cluster_rows = T,
                        cluster_cols = F,
                        colorRampPalette(c("#5390c1", "white", "#ea513f"))(100), show_rownames = F, treeheight_col =0, treeheight_row = 20)


# PCA
library(factoextra)
library(ggplot2)
library(ggpubr)
metadata = data.frame(cell_line = rep(c("LTED", "TAMR", "WT"), each = 3))
rownames(metadata) = colnames(full)[1:9]

pca_res = prcomp(t(na.omit(full[,1:9])), scale = T)

plot.new()
f <- fviz_eig(pca_res)
p <- fviz_pca_ind(pca_res,
                  col.ind = metadata$cell_line,
                  repel = TRUE, 
                  geom = c("point", "text"),
                  pointsize = 2,
                  mean.point = FALSE,
                  #addEllipses=TRUE,
                  title = "PCA")

c_var <- fviz_contrib(pca_res, choice = "var", axes = c(1,2),top = 10)
p_cont <- fviz_pca_biplot(pca_res,
                          col.ind = metadata$cell_line,
                          repel = TRUE, 
                          geom.ind = c("point"),
                          pointsize = 4,
                          mean.point = FALSE,
                          col.var = "black",
                          select.var = list(name = NULL, cos2 = 10, contrib = NULL),
                          #addEllipses=TRUE,
                          title = "PCA")

ggarrange(f, p + theme_minimal()+
            scale_shape_manual(values=rep(16, 12)),
          c_var, p_cont + theme_minimal()+
            scale_shape_manual(values=rep(16, 12)),
          
          widths = c(1, 2, 1, 2),
          ncol = 2, nrow = 2)


# # stats from EV, already filtered for p < 0.05 and FC > 0.5
# setwd("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/normal_growth_condition/FULL/from_EV/")
# 
# LTED = read.csv("LTED_vs_WT_NormGrowth_FullProteome_09112023.tsv", sep = "\t", header = T)
# TAMR = read.csv("TAMR_vs_WT_NormGrowth_FullProteome_09112023.tsv", sep = "\t", header = T)
# 
# LTED$Gene = sapply(strsplit(LTED$Protein_ID, split = "_"), function(x) x[2])
# TAMR$Gene = sapply(strsplit(TAMR$Protein_ID, split = "_"), function(x) x[2])
# rownames(LTED) = LTED$Gene
# rownames(TAMR) = TAMR$Gene
# 
# TAMR_merged = merge(TAMR_stats, TAMR, by = 0)
# LTED_merged = merge(LTED_stats, LTED, by = 0)
# # identical FC, T, p-value as Efstathios' statistics
# --> use mine from now on since it's not filtered yet

library(wordcloud)

WT_col = "#6bc4ca"
LTED_col = "#e9457b"
TAMR_col = "#896dae"


LTED_stats$color = rep("#00000050", nrow (LTED_stats))
LTED_stats$color[LTED_stats$adj.P.Val < 0.05 & LTED_stats$logFC < -1] <- paste0(WT_col,"50")
LTED_stats$color[LTED_stats$adj.P.Val < 0.05 & LTED_stats$logFC > 1] <- paste0(LTED_col,"50")

TAMR_stats$color = rep("#00000050", nrow (TAMR_stats))
TAMR_stats$color[TAMR_stats$adj.P.Val < 0.05 & TAMR_stats$logFC < -1] <- paste0(WT_col,"50")
TAMR_stats$color[TAMR_stats$adj.P.Val < 0.05 & TAMR_stats$logFC > 1] <- paste0(TAMR_col,"50")

plot(LTED_stats$logFC, -log10(LTED_stats$adj.P.Val), xlim = range(LTED_stats$logFC)*1.2, 
     ylim = c(0,range(-log10(LTED_stats$adj.P.Val))[2]*1.1), col = LTED_stats$color, pch = 16,
     main = "normal growth condition, full proteome\n LTED vs. WT", xlab = "log2 FC", ylab = "-log10(adjusted p-value)")
textplot(LTED_stats$logFC[LTED_stats$adj.P.Val < 0.05 & abs(LTED_stats$logFC) > 3],
         -log10(LTED_stats$adj.P.Val)[LTED_stats$adj.P.Val < 0.05 & abs(LTED_stats$logFC) > 3],
         rownames(LTED_stats)[LTED_stats$adj.P.Val < 0.05 & abs(LTED_stats$logFC) > 3], new = F, show.lines = T)



plot(TAMR_stats$logFC, -log10(TAMR_stats$adj.P.Val), xlim = range(TAMR_stats$logFC)*1.2, 
     ylim = c(0,range(-log10(TAMR_stats$adj.P.Val))[2]*1.1), col = TAMR_stats$color, pch = 16,
     main = "normal growth condition, full proteome\n TAMR vs. WT", xlab = "log2 FC", ylab = "-log10(adjusted p-value)")
textplot(TAMR_stats$logFC[TAMR_stats$adj.P.Val < 0.05 & abs(TAMR_stats$logFC) > 2.5],
         -log10(TAMR_stats$adj.P.Val)[TAMR_stats$adj.P.Val < 0.05 & abs(TAMR_stats$logFC) > 2.5],
         rownames(TAMR_stats)[TAMR_stats$adj.P.Val < 0.05 & abs(TAMR_stats$logFC) > 2.5], new = F, show.lines = T)


# Reactome GSEA:
#BiocManager::install("ReactomePA")
library(ReactomePA)
require(DOSE)
library(stats)

# overrepresentation analysis of significant proteins:
# # map gene names to Entrez ID
library('org.Hs.eg.db')
entrez_ID_lted = mapIds(org.Hs.eg.db, rownames(LTED_stats), 'ENTREZID', 'SYMBOL')
entrez_ID_lted = as.vector(unname(entrez_ID_lted))
LTED_stats = data.frame(LTED_stats, entrez_ID_lted)

entrez_ID_tamr = mapIds(org.Hs.eg.db, rownames(TAMR_stats), 'ENTREZID', 'SYMBOL')
entrez_ID_tamr = as.vector(unname(entrez_ID_tamr))
TAMR_stats = data.frame(TAMR_stats, entrez_ID_tamr)


LTED_UP = LTED_stats$entrez_ID[LTED_stats$adj.P.Val < 0.05 & LTED_stats$logFC > 0]
TAMR_UP = TAMR_stats$entrez_ID[TAMR_stats$adj.P.Val < 0.05 & TAMR_stats$logFC > 0]
LTED_DOWN = LTED_stats$entrez_ID[LTED_stats$adj.P.Val < 0.05 & LTED_stats$logFC < 0]
TAMR_DOWN = TAMR_stats$entrez_ID[TAMR_stats$adj.P.Val < 0.05 & TAMR_stats$logFC < 0]
LTED_sign = LTED_stats$entrez_ID[LTED_stats$adj.P.Val < 0.05 & abs(LTED_stats$logFC) > 1.5]
TAMR_sign = TAMR_stats$entrez_ID[TAMR_stats$adj.P.Val < 0.05 & abs(TAMR_stats$logFC) > 1.5]

LTED_up_enrich <- enrichPathway(gene=LTED_UP,pvalueCutoff=0.05, readable=T, universe = LTED_stats$entrez_ID, pAdjustMethod = "BH")
TAMR_up_enrich <- enrichPathway(gene=TAMR_UP,pvalueCutoff=0.05, readable=T, universe = TAMR_stats$entrez_ID, pAdjustMethod = "BH")
LTED_down_enrich <- enrichPathway(gene=LTED_DOWN,pvalueCutoff=0.05, readable=T, universe = LTED_stats$entrez_ID, pAdjustMethod = "BH")
TAMR_down_enrich <- enrichPathway(gene=TAMR_DOWN,pvalueCutoff=0.05, readable=T, universe = TAMR_stats$entrez_ID, pAdjustMethod = "BH")
LTED_sign_enrich <- enrichPathway(gene=LTED_sign,pvalueCutoff=0.05, readable=T, universe = LTED_stats$entrez_ID, pAdjustMethod = "BH")
TAMR_sign_enrich <- enrichPathway(gene=TAMR_sign,pvalueCutoff=0.05, readable=T, universe = TAMR_stats$entrez_ID, pAdjustMethod = "BH")



dotplot(LTED_up_enrich,showCategory=10)
dotplot(TAMR_up_enrich,  showCategory=10)
dotplot(LTED_down_enrich, showCategory=10)
dotplot(TAMR_down_enrich, showCategory=10)
dotplot(LTED_sign_enrich, showCategory=10)
dotplot(TAMR_sign_enrich, showCategory=10)

gcSample= list(LTED_UP = LTED_stats$entrez_ID[LTED_stats$adj.P.Val < 0.05 & LTED_stats$logFC > 1], 
               LTED_DOWN =  LTED_stats$entrez_ID[LTED_stats$adj.P.Val < 0.05 & LTED_stats$logFC < 1],
               TAMR_UP = TAMR_stats$entrez_ID[TAMR_stats$adj.P.Val < 0.05 & TAMR_stats$logFC > 1],
               TAMR_DOWN = TAMR_stats$entrez_ID[TAMR_stats$adj.P.Val < 0.05 & TAMR_stats$logFC < 1]
)

ck <- compareCluster(geneCluster = gcSample, fun = enrichPathway)

dotplot(ck, showCategory = 3)


#cholesterol biosynthesis term:


cholesterol = ck@compareClusterResult$geneID[ck@compareClusterResult$Description == "Cholesterol biosynthesis"]
cholesterol = unlist((strsplit(cholesterol, split = "/", fixed = T)))

cholesterol_genes = mapIds(org.Hs.eg.db, cholesterol, 'SYMBOL', 'ENTREZID')

# plot volcano:
LTED_stats$new_color = rep("#00000050", nrow (LTED_stats))
LTED_stats$new_color[rownames(LTED_stats) %in% cholesterol_genes] <- paste0(LTED_col,"50")
TAMR_stats$new_color = rep("#00000050", nrow (TAMR_stats))
TAMR_stats$new_color[rownames(TAMR_stats) %in% cholesterol_genes] <- paste0(TAMR_col,"50")

plot(LTED_stats$logFC, -log10(LTED_stats$adj.P.Val), xlim = range(LTED_stats$logFC)*1.2, 
     ylim = c(0,range(-log10(LTED_stats$adj.P.Val))[2]*1.1), col = LTED_stats$new_color, pch = 16,
     main = "normal growth condition, full proteome\n LTED vs. WT", xlab = "log2 FC", ylab = "-log10(adjusted p-value)")
points(LTED_stats$logFC[LTED_stats$new_color == paste0(LTED_col,"50")], -log10(LTED_stats$adj.P.Val)[LTED_stats$new_color == paste0(LTED_col,"50")], col = LTED_col, pch = 16)

textplot(LTED_stats$logFC[LTED_stats$new_color == paste0(LTED_col,"50")],
         -log10(LTED_stats$adj.P.Val)[LTED_stats$new_color == paste0(LTED_col,"50")],
         rownames(LTED_stats)[LTED_stats$new_color == paste0(LTED_col,"50")], new = F, show.lines = F, col = LTED_col)


plot(TAMR_stats$logFC, -log10(TAMR_stats$adj.P.Val), xlim = range(TAMR_stats$logFC)*1.2, 
     ylim = c(0,range(-log10(TAMR_stats$adj.P.Val))[2]*1.1), col = TAMR_stats$new_color, pch = 16,
     main = "normal growth condition, full proteome\n TAMR vs. WT", xlab = "log2 FC", ylab = "-log10(adjusted p-value)")
points(TAMR_stats$logFC[TAMR_stats$new_color == paste0(TAMR_col,"50")], -log10(TAMR_stats$adj.P.Val)[TAMR_stats$new_color == paste0(TAMR_col,"50")], col = TAMR_col, pch = 16)

textplot(TAMR_stats$logFC[TAMR_stats$new_color == paste0(TAMR_col,"50")],
         -log10(TAMR_stats$adj.P.Val)[TAMR_stats$new_color == paste0(TAMR_col,"50")],
         rownames(TAMR_stats)[TAMR_stats$new_color == paste0(TAMR_col,"50")], new = F, show.lines = F, col = TAMR_col)




# pathway enrichment:
gene_list = LTED_stats$logFC
names(gene_list) = LTED_stats$entrez_ID

gene_list_sorted = gene_list[order(gene_list, decreasing = T)]
gene_list_sorted = gene_list_sorted[!duplicated(names(gene_list_sorted))]


require(ReactomePA)

y <- gsePathway(gene_list_sorted, 
                minGSSize=120, pvalueCutoff=0.05, 
                pAdjustMethod="BH", verbose=FALSE)


dotplot(y, showCategory=17, x = "NES", color = "p.adjust", title = "LTED vs. WT Reactome GSEA")



y_readable <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
y_readable@result$Description

y.heat = heatplot(y_readable,showCategory = "Hemostasis" , foldChange=gene_list_sorted)
y.heat + scale_fill_gradientn(colors = c("#5390c1","white", "#ea513f"),limits = c(-4,4))
                       
                       


# pathway enrichment TAMR:
gene_list = TAMR_stats$logFC
names(gene_list) = TAMR_stats$entrez_ID

gene_list_sorted = gene_list[order(gene_list, decreasing = T)]
gene_list_sorted = gene_list_sorted[!duplicated(names(gene_list_sorted))]


require(ReactomePA)

y <- gsePathway(gene_list_sorted, 
                minGSSize=120, pvalueCutoff=0.05, 
                pAdjustMethod="BH", verbose=FALSE)

length(y@result$ID)
dotplot(y, showCategory=5, x = "NES", color = "p.adjust", title = "TAMR vs. WT Reactome GSEA")



y_readable <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
y_readable@result$Description

y.heat = heatplot(y_readable,showCategory = y_readable@result$Description, foldChange=gene_list_sorted)
y.heat + scale_fill_gradientn(colors = c("#5390c1","white", "#ea513f"),limits = c(-4,4))





# Phospho proteome:
setwd("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/normal_growth_condition/PHOSPHO/")

library(PhosR)
load("Normal_Growth_Phospho_Normalized_Filt_Proc_05122023.RData")
phospho = ppe_filt@assays@data@listData$Normalization

phospho = data.frame(phospho)
phospho$gene_site = sapply(strsplit(rownames(phospho), split = ";"), function(x) paste0(x[2],"_",x[3]))


pheatmap(na.omit(phospho[,c(7:9,1:6)]), scale = "row", cluster_rows = T,
         cluster_cols = F, 
         colorRampPalette(c("#5390c1", "white", "#ea513f"))(100), show_rownames = F, treeheight_col =0, treeheight_row = 20)

phospho_ordered = phospho[order(apply(phospho[,1:9], 1, sd, na.rm=T), decreasing = T),]
pheatmap(na.omit(phospho_ordered[1:2000,c(7:9,1:6)]), scale = "row", cluster_rows = T,
         cluster_cols = F, 
         colorRampPalette(c("#5390c1", "white", "#ea513f"))(100), show_rownames = F, treeheight_col =0, treeheight_row = 20)


# statistics:
library(limma)
condition <- factor(sapply(strsplit(colnames(phospho)[1:9], split = "_"), function(x) x[3]))
table(condition)

replicate = factor(sapply(strsplit(colnames(phospho)[1:9], split = "_"), function(x) x[4]))
table(replicate)



design <- model.matrix(~0+ condition + replicate)
colnames(design) = c(levels(condition), levels(replicate)[-1])
nrow(design)
ncol(phospho[,1:9])

# fit limma linear model:

fit <- limma::lmFit(phospho[,1:9], design=design, weights=NULL)



# contrast matrix:
condition1 = "TAMR"
condition2 = "WT"

contrast_name <- paste0(condition1, '-', condition2)
contrast_name


contrast_matrix <- limma::makeContrasts(contrasts=eval(contrast_name),
                                        levels=design)

fit_contrast <- limma::contrasts.fit(fit, contrast_matrix)
# statistics with eBayes method:
stats <- eBayes(fit_contrast, trend=T,robust = T)

# get result in different format:
# # change name!!!
library(tidyverse)
stats_output <- limma::topTable(stats, coef=1, number=nrow(stats), adjust.method="BH", sort.by="none") %>% 
  drop_na()

write.csv(stats_output, file = "20240523_PHOSPHO_TAMR_vs_WT_limma_stats.csv", row.names = T)

LTED_stats = read.csv("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/normal_growth_condition/PHOSPHO/20240523_PHOSPHO_LTED_vs_WT_limma_stats.csv", row.names = 1)
TAMR_stats = read.csv("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/normal_growth_condition/PHOSPHO/20240523_PHOSPHO_TAMR_vs_WT_limma_stats.csv", row.names = 1)


phospho_stats_merged = merge(LTED_stats, TAMR_stats, by = 0, all = T)
rownames(phospho_stats_merged) = phospho_stats_merged$Row.names
phospho_stats_merged = phospho_stats_merged[,-1]
colnames(phospho_stats_merged) = paste0(rep(c("LTED_", "TAMR_"), each = 6),colnames(phospho_stats_merged))

phospho_stats_merged$gene_site = sapply(strsplit(rownames(phospho_stats_merged), split = ";", fixed = T), function(x) paste0(x[2], "_", x[3]))



#heatmap of PAK1 target sites:
full_merged = merge(LTED_stats, TAMR_stats, by = 0, all = T)
rownames(full_merged) = sapply(strsplit(full_merged$Row.names, "_"), function(x) x[2])
full_merged = full_merged[,-1]



phospho_PAK1_targets = phospho_stats_merged[phospho_stats_merged$gene_site %in% KSN_PAK1_targets$phospho_site,]
rownames(phospho_PAK1_targets) = phospho_PAK1_targets$gene_site


par(mfrow=c(1,2))


p = pheatmap(phospho_PAK1_targets[,c(1,7)], scale = "none", cluster_cols = F,cluster_rows = T, colorRampPalette(c("#5390c1", "white", "#ea513f"))(100), breaks = seq(-2,2,4/100), border_color = NA, treeheight_row = 0)

PAK1_target_genes = unique(sapply(strsplit(phospho_PAK1_targets$gene_site[p$tree_row$order], "_"), function(x) x[1]))

full_PAK1_targets = full_merged[PAK1_target_genes,]

pheatmap(full_PAK1_targets[,c(1,7)], scale = "none", cluster_cols = F,cluster_rows = F, colorRampPalette(c("#5390c1", "white", "#ea513f"))(100), breaks = seq(-2,2,4/100), border_color = NA, treeheight_row = 0)


#log2 FC:


# kinase activities:
setwd("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/normal_growth_condition/PHOSPHO")
LTED_kinases = read.csv("Out_RE_Analysis_NormalGrowth_Condition_Omnipath_Kin_Act_LTED_vs_WT_09112023.tsv", sep = "\t", header = T, stringsAsFactors = F)
TAMR_kinases = read.csv("Out_RE_Analysis_NormalGrowth_Condition_Omnipath_Kin_Act_TAMR_vs_WT_09112023.tsv", sep = "\t", header = T, stringsAsFactors = F)

kinases_merged = merge(LTED_kinases, TAMR_kinases, by = "Omnipath_Kinase", all = T)
kinases_merged = na.omit(kinases_merged)

kinases_merged$col = rep("#000000", nrow(kinases_merged))
kinases_merged$col[abs(kinases_merged$score.x) > 2] <- LTED_col
kinases_merged$col[abs(kinases_merged$score.y) > 2] <- TAMR_col
kinases_merged$col[abs(kinases_merged$score.y) > 2 & abs(kinases_merged$score.x) > 2] <- "blue"

plot(kinases_merged$score.x, kinases_merged$score.y, pch = 16,main="Kinase activities (normal growth condition)",
     xlab = "LTED vs. WT", ylab = "TAMR vs. WT", col = kinases_merged$col, cex = 1.2, xlim = c(-7,7), ylim =c(-4,4))
abline(h=2, lty = 2, col = TAMR_col)
abline(h=-2, lty = 2, col = TAMR_col)
abline(v=2, lty = 2, col = LTED_col)
abline(v=-2, lty = 2, col = LTED_col)

wordcloud::textplot(kinases_merged$score.x[abs(kinases_merged$score.x) > 2 | abs(kinases_merged$score.y) > 2],
                    kinases_merged$score.y[abs(kinases_merged$score.x) > 2 | abs(kinases_merged$score.y) > 2],
                    kinases_merged$Omnipath_Kinase[abs(kinases_merged$score.x) > 2 | abs(kinases_merged$score.y) > 2], new = F, show.lines = T)


# merge kinase activities with FP data:
rownames(kinases_merged) = kinases_merged$Omnipath_Kinase

LTED_kinase_full_merged = merge(kinases_merged, LTED_stats, by = 0, all.x = T)
LTED_kinase_full_merged = na.omit(LTED_kinase_full_merged)
#LTED:
plot(LTED_kinase_full_merged$score.x, LTED_kinase_full_merged$logFC, pch = 16, col = LTED_col, xlab = "LTED vs. WT kinase score", ylab = "LTED vs. WT log2FC")
library(wordcloud)
textplot(LTED_kinase_full_merged$score.x[abs(LTED_kinase_full_merged$score.x) > 2],
                    LTED_kinase_full_merged$logFC[abs(LTED_kinase_full_merged$score.x) > 2],
                    LTED_kinase_full_merged$Omnipath_Kinase[abs(LTED_kinase_full_merged$score.x) > 2], new = F, show.lines = T)
#TAMR:
TAMR_kinase_full_merged = merge(kinases_merged, TAMR_stats, by = 0, all.x = T)
TAMR_kinase_full_merged = na.omit(TAMR_kinase_full_merged)

plot(TAMR_kinase_full_merged$score.y, TAMR_kinase_full_merged$logFC,xlim=range(TAMR_kinase_full_merged$score.x, na.rm=T),ylim=range(LTED_kinase_full_merged$logFC, na.rm=T), pch = 16, cex = 1.5,col = TAMR_col, xlab = "TAMR vs. WT kinase score", ylab = "TAMR vs. WT log2FC")

textplot(TAMR_kinase_full_merged$score.y[abs(TAMR_kinase_full_merged$score.y) > 2 | TAMR_kinase_full_merged$Omnipath_Kinase == "PAK1"],
         TAMR_kinase_full_merged$logFC[abs(TAMR_kinase_full_merged$score.y) > 2 | TAMR_kinase_full_merged$Omnipath_Kinase == "PAK1"],
         TAMR_kinase_full_merged$Omnipath_Kinase[abs(TAMR_kinase_full_merged$score.y) > 2 | TAMR_kinase_full_merged$Omnipath_Kinase == "PAK1"], new = F, show.lines = T)

abline(h=0, lty = 2)
abline(v=0, lty = 2)
#abline(v=-2, lty = 2)

# get Omnipath prior knowledge:
library(OmnipathR)
# import KSN from omnipath
omnipath_ptm <- get_signed_ptms()

# likely erroneous database records: As all databases, the resources constituting OmniPath contain wrong records, some of them more, others less. Here we see all these records are from only one resource called ProtMapper. ProtMapper uses literature mining algorithms (REACH, Sparser), which are known to produce a number of false positives in their output. On this premise:

omnipath_ptm <- omnipath_ptm %>% filter(!grepl('ProtMapper', sources) | n_resources > 1)

# Which means: not from ProtMapper, or from more than one resource (i.e. if it is from ProtMapper, but also another resource confirms, then we keep the record).

omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
KSN <- omnipath_ptm[,c(4,3)]
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
KSN$mor <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)
KSN$likelihood <- 1

#we remove ambiguous modes of regulations
KSN$id <- paste(KSN$substrate_genesymbol,KSN$enzyme_genesymbol,sep ="")
KSN <- KSN[!duplicated(KSN$id),]
KSN <- KSN[,-5]

#rename KSN to fit decoupler format for downstream analysis
names(KSN)[c(1,2)] <- c("phospho_site","kinase")

KSN_PAK1_targets = KSN[KSN$kinase == "PAK1",]

# which kinases have PAK1 as target?
KSN_PAK1 = KSN[grep("PAK1_", KSN$phospho_site),]

upstream_kinases = unique(KSN_PAK1$kinase)
full = data.frame(full)
full$gene = sapply(strsplit(rownames(full), split = "_"), function(x) x[2])

full_upstream = full[rownames(full) %in% upstream_kinases,]

full_upstream$LTED = apply(full_upstream[,1:3], 1, mean, na.rm=T)
full_upstream$TAMR = apply(full_upstream[,4:6], 1, mean, na.rm=T)
full_upstream$WT = apply(full_upstream[,7:9], 1, mean, na.rm=T)



library(pheatmap)
pheatmap(na.omit(full_upstream[,1:9]),cluster_cols = F, scale = "row")

#activities of these kinases:
LTED_upstream_kinases_activities = LTED_kinases[LTED_kinases$Omnipath_Kinase %in% upstream_kinases,]
TAMR_upstream_kinases_activities = TAMR_kinases[TAMR_kinases$Omnipath_Kinase %in% upstream_kinases,]

par(mfrow=c(1,2))
barplot(LTED_upstream_kinases_activities$score[order(LTED_upstream_kinases_activities$score,decreasing  = T)],ylim=c(-5,6), names.arg = LTED_upstream_kinases_activities$Omnipath_Kinase[order(LTED_upstream_kinases_activities$score, decreasing = T)], las = 2, main = "kinase activities upstream of PAK1\n (LTED vs. WT)", col = LTED_col)
barplot(TAMR_upstream_kinases_activities$score[order(TAMR_upstream_kinases_activities$score,decreasing  = T)],ylim=c(-5,6), names.arg = TAMR_upstream_kinases_activities$Omnipath_Kinase[order(TAMR_upstream_kinases_activities$score, decreasing = T)], las = 2, main = "kinase activities upstream of PAK1\n (TAMR vs. WT)", col = TAMR_col)

# merge upstream activities and kinases:
upstream_kinases_merged = merge(full_upstream, LTED_upstream_kinases_activities, by.x = 0, by.y = "Omnipath_Kinase", all = T)
upstream_kinases_merged = merge(upstream_kinases_merged, TAMR_upstream_kinases_activities, by.x = 1, by.y = "Omnipath_Kinase", all = T)

rownames(upstream_kinases_merged) = upstream_kinases_merged$Row.names
upstream_kinases_merged = upstream_kinases_merged[,-1]


plot(upstream_kinases_merged$score.x, upstream_kinases_merged$LTED-upstream_kinases_merged$WT,xlim=c(-6,6),main = "kinases upstream of PAK1",ylim=c(-1,2), col = LTED_col, xlab = "kinase activity",pch=16, ylab = "kinase expression log2FC")
textplot(upstream_kinases_merged$score.x[!is.na(upstream_kinases_merged$score.x) & !is.na(upstream_kinases_merged$LTED-upstream_kinases_merged$WT)],
         c(upstream_kinases_merged$LTED-upstream_kinases_merged$WT)[!is.na(upstream_kinases_merged$score.x) & !is.na(upstream_kinases_merged$LTED-upstream_kinases_merged$WT)],
         rownames(upstream_kinases_merged)[!is.na(upstream_kinases_merged$score.x) & !is.na(upstream_kinases_merged$LTED-upstream_kinases_merged$WT)], new = F)

plot(upstream_kinases_merged$score.y, upstream_kinases_merged$TAMR-upstream_kinases_merged$WT, xlim=c(-6,6),main = "kinases upstream of PAK1",ylim=c(-1,2),col = TAMR_col, xlab = "kinase activity",pch=16, ylab = "kinase expression log2FC")
textplot(upstream_kinases_merged$score.y[!is.na(upstream_kinases_merged$score.y) & !is.na(upstream_kinases_merged$TAMR-upstream_kinases_merged$WT)],
         c(upstream_kinases_merged$TAMR-upstream_kinases_merged$WT)[!is.na(upstream_kinases_merged$score.y) & !is.na(upstream_kinases_merged$TAMR-upstream_kinases_merged$WT)],
         rownames(upstream_kinases_merged)[!is.na(upstream_kinases_merged$score.y) & !is.na(upstream_kinases_merged$TAMR-upstream_kinases_merged$WT)], new = F)


# what does this look like at time point 0?
LTED_kinases_tp0 = read.csv("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/20230912_PHOSPHO/kinase_activities_compared_to_WT/Out_Kin_Act_Total_Per_Timepoint/Out_RE_Analysis_Omnipath_Kin_Act_LTED_0_09112023.tsv", sep = "\t")
TAMR_kinases_tp0 = read.csv("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/20230912_PHOSPHO/kinase_activities_compared_to_WT/Out_Kin_Act_Total_Per_Timepoint/Out_RE_Analysis_Omnipath_Kin_Act_TAMR_0_09112023.tsv", sep = "\t")

LTED_upstream_kinases_activities_tp0 = LTED_kinases_tp0[LTED_kinases_tp0$Omnipath_Kinase %in% upstream_kinases,]
TAMR_upstream_kinases_activities_tp0 = TAMR_kinases_tp0[TAMR_kinases_tp0$Omnipath_Kinase %in% upstream_kinases,]

par(mfrow=c(1,2))
barplot(LTED_upstream_kinases_activities_tp0$score[order(LTED_upstream_kinases_activities_tp0$score,decreasing  = T)],ylim=c(-5,6), names.arg = LTED_upstream_kinases_activities_tp0$Omnipath_Kinase[order(LTED_upstream_kinases_activities_tp0$score, decreasing = T)], las = 2, main = "TP0: kinase activities upstream of PAK1 (LTED vs. WT)")
barplot(TAMR_upstream_kinases_activities_tp0$score[order(TAMR_upstream_kinases_activities_tp0$score,decreasing  = T)],ylim=c(-5,6), names.arg = TAMR_upstream_kinases_activities_tp0$Omnipath_Kinase[order(TAMR_upstream_kinases_activities_tp0$score, decreasing = T)], las = 2, main = "TP0: kinase activities upstream of PAK1 (TAMR vs. WT)")


# what does this look like at time point 10?
LTED_kinases_tp10 = read.csv("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/20230912_PHOSPHO/kinase_activities_compared_to_WT/Out_Kin_Act_Total_Per_Timepoint/Out_RE_Analysis_Omnipath_Kin_Act_LTED_10_09112023.tsv", sep = "\t")
TAMR_kinases_tp10 = read.csv("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/20230912_PHOSPHO/kinase_activities_compared_to_WT/Out_Kin_Act_Total_Per_Timepoint/Out_RE_Analysis_Omnipath_Kin_Act_TAMR_10_09112023.tsv", sep = "\t")

LTED_upstream_kinases_activities_tp10 = LTED_kinases_tp10[LTED_kinases_tp10$Omnipath_Kinase %in% upstream_kinases,]
TAMR_upstream_kinases_activities_tp10 = TAMR_kinases_tp10[TAMR_kinases_tp10$Omnipath_Kinase %in% upstream_kinases,]

par(mfrow=c(1,2))
barplot(LTED_upstream_kinases_activities_tp10$score[order(LTED_upstream_kinases_activities_tp10$score,decreasing  = T)],ylim=c(-5,6), names.arg = LTED_upstream_kinases_activities_tp10$Omnipath_Kinase[order(LTED_upstream_kinases_activities_tp10$score, decreasing = T)], las = 2, main = "TP10: kinase activities upstream of PAK1 (LTED vs. WT)")
barplot(TAMR_upstream_kinases_activities_tp10$score[order(TAMR_upstream_kinases_activities_tp10$score,decreasing  = T)],ylim=c(-5,6), names.arg = TAMR_upstream_kinases_activities_tp10$Omnipath_Kinase[order(TAMR_upstream_kinases_activities_tp10$score, decreasing = T)], las = 2, main = "TP10: kinase activities upstream of PAK1 (TAMR vs. WT)")


# phosphosites on PAK1:
phospho = data.frame(phospho)
phospho$gene_site = sapply(strsplit(rownames(phospho), split = ";", fixed = T), function(x) paste0(x[2], "_", x[3]))

phospho_PAK1 = phospho[grep("PAK1_",phospho$gene_site),]
# only this one site found: S144 (autophospho site)




###################################

PAK1_targets = KSN_PAK1_targets
PAK1_targets$Target_Protein = sapply(strsplit(PAK1_targets$phospho_site, split = "_"), function(x) x[1])

full_PAK1_targets = full[rownames(full) %in% PAK1_targets$Target_Protein,]
full_PAK1_targets = data.frame(full_PAK1_targets)

full_PAK1_targets$LTED = apply(full_PAK1_targets[,1:3],1,mean,na.rm=T)
full_PAK1_targets$TAMR = apply(full_PAK1_targets[,4:6],1,mean,na.rm=T)
full_PAK1_targets$WT = apply(full_PAK1_targets[,7:9],1,mean,na.rm=T)

full_PAK1_targets = full_PAK1_targets[rownames(full_PAK1_targets) %in% phospho_PAK1_targets$gene,]

pheatmap(full_PAK1_targets[,1:9],cluster_rows = F, cluster_cols = F, scale = "row", colorRampPalette(c("#5390c1", "white", "#ea513f"))(100), breaks = seq(-2,2,4/100), border_color = NA, treeheight_row = 0)


# PAK1 targets in phospho data of normal growth condition:
PAK1_targets$gene_site = PAK1_targets$phospho_site

rownames(phospho)[1:10]

rownames(phospho)[!duplicated(phospho$gene_site)] = phospho$gene_site[!duplicated(phospho$gene_site)]


phospho_PAK1_targets = phospho[rownames(phospho) %in% PAK1_targets$gene_site,]
# 30 are found
phospho_PAK1_targets = data.frame(phospho_PAK1_targets)
phospho_PAK1_targets$gene = sapply(strsplit(phospho_PAK1_targets$gene_site, "_"), function(x) x[1])

phospho_PAK1_targets$LTED = apply(phospho_PAK1_targets[,1:3],1,mean,na.rm=T)
phospho_PAK1_targets$TAMR = apply(phospho_PAK1_targets[,4:6],1,mean,na.rm=T)
phospho_PAK1_targets$WT = apply(phospho_PAK1_targets[,7:9],1,mean,na.rm=T)

pheatmap(phospho_PAK1_targets[,10:12],cluster_cols = F,cluster_rows = F, scale = "row")

# calculate means:
normalGrowth_phospho_PAK1_targets_means = data.frame(matrix(NA, nrow=nrow(phospho_PAK1_targets), ncol = 3))
colnames(normalGrowth_phospho_PAK1_targets_means) = c("LTED", "TAMR","WT")
rownames(normalGrowth_phospho_PAK1_targets_means) = rownames(phospho_PAK1_targets)

subtract = rep(1:9, each = 3)

for (i in seq(1,9,3)){
  normalGrowth_phospho_PAK1_targets_means[,i-(i-subtract[i])] <- apply(phospho_PAK1_targets[,i:(i+2)],1,mean,na.rm=T)
}


# merge normalGrowth with timeCourse_TP0:

# Phospho proteome:
setwd("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/from_Efstathios/")

library(PhosR)
load("phospho_normalized/RE_Analysis_TimeCourse_Phospho_Norm_Filt_Proc_09102023.RData")

time_course_phospho = data.frame(ppe_filt@assays@data@listData$Normalization)

time_course_phospho$gene_site = sapply(strsplit(rownames(time_course_phospho), split = ";", fixed = T), function(x) paste0(x[2],"_",x[3]))
time_course_phospho = time_course_phospho[!duplicated(time_course_phospho$gene_site),]

rownames(time_course_phospho) = time_course_phospho$gene_site

time_course_phospho_PAK1_targets = time_course_phospho[rownames(time_course_phospho) %in% PAK1_targets$gene_site,]

# calculate means:
time_course_phospho_PAK1_targets_means = data.frame(matrix(NA, nrow=nrow(time_course_phospho_PAK1_targets), ncol = 63/3))
colnames(time_course_phospho_PAK1_targets_means) = paste0(rep(c("LTED_", "TAMR_", "WT_"), each = 7), rep(c(0,2,5,10,20,60,120), 3))
rownames(time_course_phospho_PAK1_targets_means) = rownames(time_course_phospho_PAK1_targets)

subtract = rep(1:21, each = 3)

for (i in seq(1,63,3)){
  time_course_phospho_PAK1_targets_means[,i-(i-subtract[i])] <- apply(time_course_phospho_PAK1_targets[,i:(i+2)],1,mean,na.rm=T)
}

pheatmap(time_course_phospho_PAK1_targets_means, scale = "row", cluster_rows = T, cluster_cols = F)

# merge normalGrowth with timeCourse_TP0:
phospho_normal_TC0_PAK1Targets_merged = merge(normalGrowth_phospho_PAK1_targets_means, time_course_phospho_PAK1_targets_means[,c(1,8,15)], by = 0, all = T)
rownames(phospho_normal_TC0_PAK1Targets_merged) =phospho_normal_TC0_PAK1Targets_merged$Row.names
phospho_normal_TC0_PAK1Targets_merged = phospho_normal_TC0_PAK1Targets_merged[,-1]

pheatmap(phospho_normal_TC0_PAK1Targets_merged, scale = "none", cluster_rows = T, cluster_cols = F)


phospho_normal_TC0_PAK1Targets_merged$FC_LTED_WT_normal = apply(phospho_normal_TC0_PAK1Targets_merged[,c(1,3)],1, function(x) x[1]-x[2])
phospho_normal_TC0_PAK1Targets_merged$FC_TAMR_WT_normal = apply(phospho_normal_TC0_PAK1Targets_merged[,c(2,3)],1, function(x) x[1]-x[2])
phospho_normal_TC0_PAK1Targets_merged$FC_LTED_WT_starved = apply(phospho_normal_TC0_PAK1Targets_merged[,c(4,6)],1, function(x) x[1]-x[2])
phospho_normal_TC0_PAK1Targets_merged$FC_TAMR_WT_starved = apply(phospho_normal_TC0_PAK1Targets_merged[,c(5,6)],1, function(x) x[1]-x[2])

pheatmap(phospho_normal_TC0_PAK1Targets_merged[,c(7,9,8,10)], cluster_rows = F, cluster_cols = F, scale = "none", breaks = seq(-3,3,6/100))


for (i in 1:nrow(phospho_normal_TC0_PAK1Targets_merged)){
  pg.bar = barplot(data.matrix(phospho_normal_TC0_PAK1Targets_merged[i,c(7,9,8,10)]), las = 2, main = rownames(phospho_normal_TC0_PAK1Targets_merged)[i])
  #segments(x0 = pg.bar, y0 = (normal_TP0_PAK1Targets_merged[i,c(1:3,7:9)]-(normal_TP0_PAK1Targets_merged[i,c(4:6,10:12)], y1 = (normal_TP0_PAK1Targets_merged[i,c(1:3,7:9)]+(normal_TP0_PAK1Targets_merged[i,c(4:6,10:12)])
}



# divide phospho by full protein intensity:


rownames(phospho_PAK1_targets)
rownames(full_PAK1_targets)


phospho_PAK1_targets_norm = data.frame(phospho_PAK1_targets)
phospho_PAK1_targets_norm$gene = sapply(strsplit(rownames((phospho_PAK1_targets_norm)), split = "_"), function(x) x[1])

for (i in 1:nrow(phospho_PAK1_targets_norm)){
  phospho_PAK1_targets_norm[i,1:9] = phospho_PAK1_targets_norm[i,1:9]/full_PAK1_targets[phospho_PAK1_targets_norm$gene[i],1:9]
}

pheatmap(phospho_PAK1_targets_norm[,1:9],cluster_cols = F,cluster_rows = F, scale = "row")

phospho_PAK1_targets_norm$mean_LTED = apply(phospho_PAK1_targets_norm[,1:3], 1, mean, na.rm=T)
phospho_PAK1_targets_norm$mean_TAMR = apply(phospho_PAK1_targets_norm[,4:6], 1, mean, na.rm=T)
phospho_PAK1_targets_norm$mean_WT = apply(phospho_PAK1_targets_norm[,7:9], 1, mean, na.rm=T)


pheatmap(phospho_PAK1_targets_norm[,11:13],cluster_cols = F,cluster_rows = F, scale = "row")

