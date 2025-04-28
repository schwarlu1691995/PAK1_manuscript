# two-way anova to find phosphosites that significantly change over time and are different between the cell lines
# 1. load filtered, normalized phospho data:


norm.mat <- read.csv("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp08_MCF7_EGF_phospho_dynamics/re-analysis_SN17/20230912_PHOSPHO/20240618_filtered_phospho.csv", row.names = 1)
colnames(norm.mat)
norm.mat = norm.mat[,c(43:63, 1:42,64:66)]
head(rownames(norm.mat))
rownames(norm.mat) = norm.mat$gene_site

# keep only sites found in >2/3 of samples:

norm.mat = norm.mat[apply(norm.mat[,1:63], 1, function(x) sum(!is.na(x))) > 41,]

cell_line = sapply(strsplit(colnames(norm.mat)[1:63], split = "_", fixed = T), function(x) x[1])
table(cell_line)
cell_line = factor(cell_line, levels = c("WT", "LTED", "TAMR"))

time = sapply(strsplit(colnames(norm.mat)[1:63], split = "_", fixed = T), function(x) paste0(x[2], "_min"))
table(time)
time = factor(time)



design <- model.matrix(~cell_line*time)
design


phospho.mat = norm.mat[,1:63]
library(limma)
fit <- lmFit(phospho.mat, design)
fit2 <- eBayes(fit)

topTable(fit2, coef=2, adjust="BH") # Cell line
stats_cell_line = limma::topTable(fit2, coef=2, number=nrow(fit2), adjust.method="fdr", sort.by="none")
stats_cell_line_sign = stats_cell_line[stats_cell_line$adj.P.Val < 0.05,]
#1104 sites
stats_time = limma::topTable(fit2, coef=3, number=nrow(fit2), adjust.method="fdr", sort.by="none")
stats_time_sign = stats_time[stats_time$adj.P.Val < 0.05,]
#1593 sites
stats_interaction = limma::topTable(fit2, coef=4, number=nrow(fit2), adjust.method="fdr", sort.by="none")
stats_interaction_sign = stats_interaction[stats_interaction$adj.P.Val < 0.05,]
#279 sites


# get these sites from the original phospho data:
phospho.mat.cl.sign = phospho.mat[rownames(stats_cell_line_sign),]
phospho.mat.time.sign = phospho.mat[rownames(stats_time_sign),]
phospho.mat.interaction.sign = phospho.mat[rownames(stats_interaction_sign),]

load("K:/Ergebnisse/LS_testing/1_PhD/01_experiments/Exp43_NCT_MASTER_CRC_patient_samples/combined/final_output_tables/metadata/plotting_colors.Rdata")

plotting_order = c(1:3, 10:12, 16:18, 4:6, 13:15, 19:21, 7:9, c(1:3, 10:12, 16:18, 4:6, 13:15, 19:21, 7:9)+21, c(1:3, 10:12, 16:18, 4:6, 13:15, 19:21, 7:9)+42)
colnames(phospho.mat)[plotting_order]

pheatmap(na.omit(phospho.mat.cl.sign[,plotting_order]), cluster_rows = T, cluster_cols = F,scale="row", color = plotting_colors$heatmap, breaks = seq(-3,3,6/100))
pheatmap(na.omit(phospho.mat.time.sign[,plotting_order]), cluster_rows = T, cluster_cols = F,scale="row", color = plotting_colors$heatmap, breaks = seq(-3,3,6/100))
pheatmap((phospho.mat.interaction.sign[,plotting_order]), cluster_rows = T, cluster_cols = F,scale="row", color = plotting_colors$heatmap, breaks = seq(-3,3,6/100), fontsize = 6)

phospho.mat.interaction.sign_means = phospho.mat.interaction.sign

for (i in 1:61){
  phospho.mat.interaction.sign_means[,i] <- apply(phospho.mat.interaction.sign[,i:(i+2)],1,mean,na.rm=T)
}

phospho.mat.interaction.sign_means = phospho.mat.interaction.sign_means[,seq(1,63,3)]
pheatmap((phospho.mat.interaction.sign_means[,plotting_order_means]), cluster_rows = T, cluster_cols = F,scale="row", color = plotting_colors$heatmap, breaks = seq(-3,3,6/100), fontsize = 6)


### ORA of proteins:
phospho.mat.interaction.sign_genes = sapply(strsplit(rownames(phospho.mat.interaction.sign),"_"),function(x) x[1])
phospho.mat.interaction.sign_genes = unique(phospho.mat.interaction.sign_genes)
all_genes = unique(norm.mat$gene)

library('org.Hs.eg.db')
phospho.mat.interaction.sign_genes_entrez = mapIds(org.Hs.eg.db, phospho.mat.interaction.sign_genes, 'ENTREZID', 'SYMBOL')
all_genes_entrez = mapIds(org.Hs.eg.db, all_genes, 'ENTREZID', 'SYMBOL')

require(DOSE)
library(stats)


x_up <- enrichPathway(gene=phospho.mat.interaction.sign_genes_entrez,pvalueCutoff=0.05, readable=T, universe = all_genes_entrez, pAdjustMethod = "BH")


dotplot(x_up, showCategory=10)

#RTK signaling:
RTK_proteins= x_up@result$geneID[x_up@result$Description == "Signaling by Receptor Tyrosine Kinases"]
RTK_proteins = unlist(strsplit(RTK_proteins,"/"))

phospho.mat.interaction.sign_RTK = phospho.mat.interaction.sign
phospho.mat.interaction.sign_RTK = phospho.mat.interaction.sign_RTK[sapply(strsplit(rownames(phospho.mat.interaction.sign),"_"),function(x) x[1]) %in% RTK_proteins,]

phospho.mat.interaction.sign_RTK = phospho.mat.interaction.sign_RTK[,plotting_order]
pheatmap((phospho.mat.interaction.sign_RTK), cluster_rows = T, cluster_cols = F,scale="row", color = plotting_colors$heatmap, breaks = seq(-3,3,6/100), fontsize = 6)

plot(rep(c(0,2,5,10,30,60,120), each = 3),phospho.mat.interaction.sign_RTK["FGFR4_S573",1:21], pch = 16, col = WT_color, ylim = c(0,22))
points(rep(c(0,2,5,10,30,60,120), each = 3),phospho.mat.interaction.sign_RTK["FGFR4_S573",22:42], pch = 16, col = LTED_color)
points(rep(c(0,2,5,10,30,60,120), each = 3),phospho.mat.interaction.sign_RTK["FGFR4_S573",43:63], pch = 16, col = TAMR_color)

# filter for presence in Omnipath:
# get prior knowledge:
# remove also KEA and phosphonetwork as sources of the omnipath_ptm

library(OmnipathR)
library(decoupleR)

uniprot_kinases <- OmnipathR::import_omnipath_annotations(resources = "UniProt_keyword") %>%
  dplyr::filter(value == "Kinase" & !grepl("COMPLEX", uniprot)) %>%
  distinct() %>%
  pull(genesymbol) %>%
  unique()
omnipath_ptm <- OmnipathR::get_signed_ptms() %>%
  dplyr::filter(modification %in% c("dephosphorylation","phosphorylation")) %>%
  dplyr::filter(!(stringr::str_detect(sources, "ProtMapper") & n_resources == 1)) %>%
  dplyr::filter(!(stringr::str_detect(sources, "KEA") & n_resources == 1)) %>%
  dplyr::filter(!(stringr::str_detect(sources, "PhosphoNetworks") & n_resources == 1)) %>%
  dplyr::mutate(p_site = paste0(substrate_genesymbol, "_", residue_type, residue_offset),
                mor = ifelse(modification == "phosphorylation", 1, -1)) %>%
  dplyr::transmute(p_site, enzyme_genesymbol, mor) %>%
  dplyr::filter(enzyme_genesymbol %in% uniprot_kinases)

omnipath_ptm$likelihood <- 1

#we remove ambiguous modes of regulations
omnipath_ptm$id <- paste(omnipath_ptm$p_site,omnipath_ptm$enzyme_genesymbol, sep ="")
omnipath_ptm <- omnipath_ptm[!duplicated(omnipath_ptm$id),]
omnipath_ptm <- omnipath_ptm[,-5]


names(omnipath_ptm)[c(1,2)] <- c("target","kinase")

phospho.mat.interaction.sign_omnipath = phospho.mat.interaction.sign[rownames(phospho.mat.interaction.sign) %in% omnipath_ptm$target,]
#60sites
pheatmap((phospho.mat.interaction.sign_omnipath[,plotting_order]), cluster_rows = T, cluster_cols = F,scale="row", color = plotting_colors$heatmap, breaks = seq(-3,3,6/100), fontsize = 6)


sites_annotation = data.frame(site = rownames(phospho.mat.interaction.sign_omnipath), kinase = rep(NA, nrow(phospho.mat.interaction.sign_omnipath)))

for (i in 1:60){
  sites_annotation$kinase[i] = paste(omnipath_ptm$kinase[omnipath_ptm$target == sites_annotation$site[i]], collapse = ";")
}

all_kinases = unlist(strsplit(sites_annotation$kinase, split = ";"))
table(all_kinases)

background = omnipath_ptm$kinase

#hypergeometry test
#how likely is it to find e.g. 19xMAPK1 in this set?

set_size = length(all_kinases)
total_size = length(background)
kinase_size = as.numeric(table(all_kinases)[names(table(all_kinases)) == "MAPK1"])
total_kinase_size = as.numeric(table(background)[names(table(background)) == "MAPK1"])

phyper(q = kinase_size - 1, m = set_size, n = total_size - set_size, k = total_kinase_size, lower.tail = FALSE)


kinase_ora = data.frame(kinase = unique(all_kinases), phyper = rep(NA, length(unique(all_kinases))))
head(kinase_ora)

for (i in 1:nrow(kinase_ora)){
  kinase_size = as.numeric(table(all_kinases)[names(table(all_kinases)) == kinase_ora$kinase[i]])
  total_kinase_size = as.numeric(table(background)[names(table(background)) == kinase_ora$kinase[i]])
  
  kinase_ora$phyper[i] <- phyper(q = kinase_size - 1, m = set_size, n = total_size - set_size, k = total_kinase_size, lower.tail = FALSE)

}

head(kinase_ora)

kinase_ora$phyper.adjust = p.adjust(kinase_ora$phyper, "BH")

kinase_ora$kinase[kinase_ora$phyper.adjust<0.05]
#overrepresented kinases in the significant phosphosites from interaction of cell line + time:  
#"AKT1"    "PDK1"    "RPS6KA1" "MAPK1"   "BRAF"    "MAPK3"   "RPS6KA3" "MAP2K2"  "MAP2K1" 
# "RPS6KB1" "ARAF"


sites_annotation$AKT1 = rep(NA, nrow(sites_annotation))
sites_annotation$AKT1[grep("AKT1", sites_annotation$kinase)] <- "yes"

sites_annotation$PDK1 = rep(NA, nrow(sites_annotation))
sites_annotation$PDK1[grep("PDK1", sites_annotation$kinase)] <- "yes"

sites_annotation$RPS6KA1 = rep(NA, nrow(sites_annotation))
sites_annotation$RPS6KA1[grep("RPS6KA1", sites_annotation$kinase)] <- "yes"

sites_annotation$MAPK1 = rep(NA, nrow(sites_annotation))
sites_annotation$MAPK1[grep("MAPK1", sites_annotation$kinase)] <- "yes"

sites_annotation$BRAF = rep(NA, nrow(sites_annotation))
sites_annotation$BRAF[grep("BRAF", sites_annotation$kinase)] <- "yes"

sites_annotation$MAPK3 = rep(NA, nrow(sites_annotation))
sites_annotation$MAPK3[grep("MAPK3", sites_annotation$kinase)] <- "yes"

sites_annotation$RPS6KA3 = rep(NA, nrow(sites_annotation))
sites_annotation$RPS6KA3[grep("RPS6KA3", sites_annotation$kinase)] <- "yes"

sites_annotation$MAP2K2 = rep(NA, nrow(sites_annotation))
sites_annotation$MAP2K2[grep("MAP2K2", sites_annotation$kinase)] <- "yes"

sites_annotation$MAP2K1 = rep(NA, nrow(sites_annotation))
sites_annotation$MAP2K1[grep("MAP2K1", sites_annotation$kinase)] <- "yes"

sites_annotation$RPS6KB1 = rep(NA, nrow(sites_annotation))
sites_annotation$RPS6KB1[grep("RPS6KB1", sites_annotation$kinase)] <- "yes"

sites_annotation$ARAF = rep(NA, nrow(sites_annotation))
sites_annotation$ARAF[grep("ARAF", sites_annotation$kinase)] <- "yes"

rownames(sites_annotation) = sites_annotation$site

pheatmap((phospho.mat.interaction.sign_omnipath[,plotting_order]), cluster_rows = T, cluster_cols = F,scale="row", color = plotting_colors$heatmap, breaks = seq(-3,3,6/100), fontsize = 6, annotation_row = sites_annotation[,3:13])


phospho.mat.interaction.sign_omnipath_means = phospho.mat.interaction.sign_omnipath

for (i in 1:61){
  phospho.mat.interaction.sign_omnipath_means[,i] <- apply(phospho.mat.interaction.sign_omnipath[,i:(i+2)],1,mean,na.rm=T)
}

phospho.mat.interaction.sign_omnipath_means = phospho.mat.interaction.sign_omnipath_means[,seq(1,63,3)]

plotting_order_means = c(1,4,6,2,5,7,3,c(1,4,6,2,5,7,3)+7, c(1,4,6,2,5,7,3)+14)
pheatmap((phospho.mat.interaction.sign_omnipath_means[,plotting_order_means]), cluster_rows = T, cluster_cols = F,scale="row", color = plotting_colors$heatmap, breaks = seq(-3,3,6/100), fontsize = 6, annotation_row = sites_annotation[,3:13])

#plot target sites of ORA kinases:
ora_kinases = kinase_ora$kinase[kinase_ora$phyper.adjust<0.05]

WT_color = "#3fb6a970"
LTED_color = "#886dae70"
TAMR_color = "#d86d1870"


plot_time_course = function(kinase){
  relevant_df = phospho.mat.interaction.sign_omnipath_means[na.omit(rownames(sites_annotation)[sites_annotation[,kinase] == "yes"]),]
  relevant_df = t(scale(t(relevant_df)))
  
  plot(c(0,2,5,10,30,60,120), relevant_df[1, plotting_order_means[1:7]], type = "l", col = WT_color, lwd = 2, ylim = range(relevant_df, na.rm=T),
       xlab = "time", ylab = "phosphosite intensity", main = paste0("target sites of kinase ", kinase))
  lines(c(0,2,5,10,30,60,120), relevant_df[1, plotting_order_means[8:14]], type = "l", col = LTED_color, lwd = 2)
  lines(c(0,2,5,10,30,60,120), relevant_df[1, plotting_order_means[15:21]], type = "l", col = TAMR_color, lwd = 2)
  
  for (i in 2:nrow(relevant_df)){
    lines(c(0,2,5,10,30,60,120), relevant_df[i, plotting_order_means[1:7]], type = "l", col = WT_color, lwd = 2)
    lines(c(0,2,5,10,30,60,120), relevant_df[i, plotting_order_means[8:14]], type = "l", col = LTED_color, lwd = 2)
    lines(c(0,2,5,10,30,60,120), relevant_df[i, plotting_order_means[15:21]], type = "l", col = TAMR_color, lwd = 2)
  }
  legend("topright", bty="n", lwd = 2, col = c(WT_color, LTED_color, TAMR_color), legend = c("WT", "LTED", "TAMR"))
  }

plot_time_course("MAP2K1")


ora_kinases

phospho.mat.interaction.sign_omnipath_means_scaled = t(scale(t(phospho.mat.interaction.sign_omnipath_means)))


cell_colors <- c(WT = WT_color, LTED = LTED_color, TAMR = TAMR_color)

plot_time_course_means = function(kinase){
  relevant_df = phospho.mat.interaction.sign_omnipath_means[na.omit(rownames(sites_annotation)[sites_annotation[,kinase] == "yes"]),]
  relevant_df = data.frame(t(scale(t(relevant_df))))
  relevant_df_means = apply(relevant_df,2,mean,na.rm=T)
  relevant_df_sd = apply(relevant_df,2,sd,na.rm=T)
  
  ploting_df = data.frame(mean = relevant_df_means[plotting_order_means], sd = relevant_df_sd[plotting_order_means], time = rep(c(0,2,5,10,30,60,120),3), cell = rep(c("WT", "LTED", "TAMR"), each = 7))
  ggplot(data = ploting_df, aes(x = time, group = cell)) + 
    geom_line(aes(y = mean, color = cell), size = 1) + 
    geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = cell), alpha = .2) +
    xlab("Minutes") + 
    theme_bw() +  
    theme(legend.key = element_blank()) + 
    theme(plot.margin=unit(c(3,3,1,1),"cm"))+
    theme(legend.position = c(1.05,.6), legend.direction = "vertical") +
    theme(legend.title = element_blank())+http://127.0.0.1:28365/graphics/5f84af04-2c59-42fc-becf-24b595ef0854.png
    ggtitle(label = paste0("target sites of kinase ", kinase, " (n=", nrow(relevant_df), ")"))+
    ylab(label = "mean scaled phosphosite intensity")+
    scale_color_manual(values = cell_colors)+
    scale_fill_manual(values = cell_colors)+
    annotate(geom="text", x=rep(100, nrow(relevant_df)), y=seq(1,2,1/(nrow(relevant_df)-1)), label=rownames(relevant_df),
             color="black")
  
  }

plot_time_course_means("MAP2K1")

ora_kinases

omnipath_ptm[omnipath_ptm$target %in% rownames(phospho.mat.interaction.sign_omnipath),]
#mor=1 for all


# plot ALL target sites of these ora_kinases:


phospho.mat_means = phospho.mat

for (i in 1:61){
  phospho.mat_means[,i] <- apply(phospho.mat_means[,i:(i+2)],1,mean,na.rm=T)
}

phospho.mat_means = phospho.mat_means[,seq(1,63,3)]
head(phospho.mat_means)

plot_time_course_means_all = function(kinase){
  kinase_sites = omnipath_ptm$target[omnipath_ptm$kinase == kinase]
  kinase_sites = kinase_sites[kinase_sites %in% rownames(phospho.mat_means)]
  relevant_df = phospho.mat_means[kinase_sites,]
  relevant_df = data.frame(t(scale(t(relevant_df))))
  relevant_df_means = apply(relevant_df,2,mean,na.rm=T)
  relevant_df_sd = apply(relevant_df,2,sd,na.rm=T)
  
  ploting_df = data.frame(mean = relevant_df_means[plotting_order_means], sd = relevant_df_sd[plotting_order_means], time = rep(c(0,2,5,10,30,60,120),3), cell = rep(c("WT", "LTED", "TAMR"), each = 7))
  ggplot(data = ploting_df, aes(x = time, group = cell)) + 
    geom_line(aes(y = mean, color = cell), size = 1) + 
    geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = cell), alpha = .2) +
    xlab("Minutes") + 
    theme_bw() +  
    theme(legend.key = element_blank()) + 
    theme(plot.margin=unit(c(3,3,1,1),"cm"))+
    theme(legend.position = c(1.05,.6), legend.direction = "vertical") +
    theme(legend.title = element_blank())+
    ggtitle(label = paste0("ALL target sites of kinase ", kinase, " (n=", nrow(relevant_df), ")"))+
    ylab(label = "mean scaled phosphosite intensity")+
    scale_color_manual(values = cell_colors)+
    scale_fill_manual(values = cell_colors)
  
}

ora_kinases
plot_time_course_means_all("ARAF")


plot_time_course_sites_of_protein = function(protein){
  kinase_sites = rownames(phospho.mat_means)[grep(protein, rownames(phospho.mat_means))]
  
  relevant_df = phospho.mat_means[kinase_sites,]
  relevant_df = data.frame(t(scale(t(relevant_df))))
  relevant_df_means = apply(relevant_df,2,mean,na.rm=T)
  relevant_df_sd = apply(relevant_df,2,sd,na.rm=T)
  
  ploting_df = data.frame(mean = relevant_df_means[plotting_order_means], sd = relevant_df_sd[plotting_order_means], time = rep(c(0,2,5,10,30,60,120),3), cell = rep(c("WT", "LTED", "TAMR"), each = 7))
  ggplot(data = ploting_df, aes(x = time, group = cell)) + 
    geom_line(aes(y = mean, color = cell), size = 1) + 
    #geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = cell), alpha = .2) +
    xlab("Minutes") + 
    theme_bw() +  
    theme(legend.key = element_blank()) + 
    theme(plot.margin=unit(c(3,3,1,1),"cm"))+
    theme(legend.position = c(1.05,.6), legend.direction = "vertical") +
    theme(legend.title = element_blank())+
    ggtitle(label = rownames(relevant_df))+
    ylab(label = "mean scaled phosphosite intensity")+
    scale_color_manual(values = cell_colors)+
    scale_fill_manual(values = cell_colors)
  
}

ora_kinases
plot_time_course_means_all("ARAF")







