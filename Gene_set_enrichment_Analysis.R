# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(readr)
library(DOSE)
library(pathview)
library(VennDiagram)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Set working directory
setwd("~/Downloads/Biocompin")

# ============================
# PART 1: Rubina Analysis
# ============================

rubina_df <- read_csv("Rubina_Genes.csv")
rubina_genes <- rubina_df$Gene.refGene
rubina_entrez <- bitr(rubina_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

write.csv(rubina_entrez, "Rubina_EntrezIDs.csv", row.names = FALSE)

# GO Enrichment
ego_rubina <- enrichGO(gene = rubina_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                       ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2, readable = TRUE)

# KEGG Enrichment
ekegg_rubina <- enrichKEGG(gene = rubina_entrez$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)

# Save plots
png("Rubina_GO_Barplot.png")
barplot(ego_rubina, showCategory = 10, title = "Rubina - GO Biological Process")
dev.off()

png("Rubina_KEGG_Dotplot.png")
dotplot(ekegg_rubina, showCategory = 20, title = "Rubina - KEGG Pathways")
dev.off()

# Filter Stroke-Related Pathways
go_rubina_df <- as.data.frame(ego_rubina)
kegg_rubina_df <- as.data.frame(ekegg_rubina)
stroke_keywords <- "stroke|brain|blood|vascular|platelet|inflamm|coagul|immune|ischemia"
stroke_go_rubina <- go_rubina_df[grep(stroke_keywords, go_rubina_df$Description, ignore.case = TRUE), ]
stroke_kegg_rubina <- kegg_rubina_df[grep(stroke_keywords, kegg_rubina_df$Description, ignore.case = TRUE), ]

write.csv(stroke_go_rubina, "Stroke_GO_Rubina.csv", row.names = FALSE)
write.csv(stroke_kegg_rubina, "Stroke_KEGG_Rubina.csv", row.names = FALSE)

# ============================
# PART 2: Sumitra Analysis
# ============================

sumitra_df <- read_csv("Sumitra_Genes.csv")
sumitra_genes <- sumitra_df$Gene.refGene
sumitra_entrez <- bitr(sumitra_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

write.csv(sumitra_entrez, "Sumitra_EntrezIDs.csv", row.names = FALSE)

# GO Enrichment
ego_sumitra <- enrichGO(gene = sumitra_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                        ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2, readable = TRUE)

# KEGG Enrichment
ekegg_sumitra <- enrichKEGG(gene = sumitra_entrez$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)

# Save plots
png("Sumitra_GO_Barplot.png")
barplot(ego_sumitra, showCategory = 10, title = "Sumitra - GO Biological Process")
dev.off()

png("Sumitra_KEGG_Dotplot.png")
dotplot(ekegg_sumitra, showCategory = 10, title = "Sumitra - KEGG Pathways")
dev.off()

# Filter Stroke-Related Pathways
go_sumitra_df <- as.data.frame(ego_sumitra)
kegg_sumitra_df <- as.data.frame(ekegg_sumitra)
stroke_go_sumitra <- go_sumitra_df[grep(stroke_keywords, go_sumitra_df$Description, ignore.case = TRUE), ]
stroke_kegg_sumitra <- kegg_sumitra_df[grep(stroke_keywords, kegg_sumitra_df$Description, ignore.case = TRUE), ]

write.csv(stroke_go_sumitra, "Stroke_GO_Sumitra.csv", row.names = FALSE)
write.csv(stroke_kegg_sumitra, "Stroke_KEGG_Sumitra.csv", row.names = FALSE)

# ============================
# PART 3: Comparison Analysis
# ============================

# Side-by-side plot for KEGG Stroke Pathways
kegg_df <- data.frame(
  Sample = c(rep("Rubina", nrow(stroke_kegg_rubina)), rep("Sumitra", nrow(stroke_kegg_sumitra))),
  Description = c(stroke_kegg_rubina$Description, stroke_kegg_sumitra$Description),
  GeneRatio = c(stroke_kegg_rubina$GeneRatio, stroke_kegg_sumitra$GeneRatio),
  p.adjust = c(stroke_kegg_rubina$p.adjust, stroke_kegg_sumitra$p.adjust)
)

kegg_df <- kegg_df %>%
  separate(GeneRatio, into = c("num", "denom"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatioValue = num / denom)

png("KEGG_Stroke_Comparison_Dotplot.png")
ggplot(kegg_df, aes(x = GeneRatioValue, y = Description, color = Sample, size = num)) +
  geom_point(alpha = 0.8) +
  labs(title = "Stroke-related KEGG Pathways: Rubina vs Sumitra",
       x = "Gene Ratio", y = "Pathway", size = "Gene Count") +
  theme_minimal() +
  scale_color_manual(values = c("Rubina" = "red", "Sumitra" = "blue"))
dev.off()

# Extract genes from platelet activation (hsa04611)
rubina_hsa04611 <- unlist(strsplit(kegg_rubina_df[kegg_rubina_df$ID == "hsa04611", "geneID"], "/"))
sumitra_hsa04611 <- unlist(strsplit(kegg_sumitra_df[kegg_sumitra_df$ID == "hsa04611", "geneID"], "/"))

# Map to SYMBOL
rubina_symbols <- bitr(rubina_hsa04611, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
sumitra_symbols <- bitr(sumitra_hsa04611, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

rubina_set <- unique(rubina_symbols$SYMBOL)
sumitra_set <- unique(sumitra_symbols$SYMBOL)

# Venn Diagram
png("Platelet_Activation_Venn.png")
venn.plot <- venn.diagram(
  x = list(Rubina = rubina_set, Sumitra = sumitra_set),
  filename = NULL,
  fill = c("tomato", "steelblue"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 1.5,
  main = "Platelet Activation (hsa04611) Gene Overlap"
)
grid.newpage()
grid.draw(venn.plot)
dev.off()

# Export unique/shared gene lists
write.csv(setdiff(rubina_set, sumitra_set), "Rubina_Unique_hsa04611_GeneSymbols.csv", row.names = FALSE)
write.csv(setdiff(sumitra_set, rubina_set), "Sumitra_Unique_hsa04611_GeneSymbols.csv", row.names = FALSE)
write.csv(intersect(rubina_set, sumitra_set), "Shared_hsa04611_GeneSymbols.csv", row.names = FALSE)
