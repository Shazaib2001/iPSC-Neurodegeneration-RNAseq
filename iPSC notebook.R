# ── 1. Install packages ────────────────────────────────────────────────────────

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "limma",          # differential expression engine (voom + linear modelling)
  "edgeR",          # filtering + voom support
  "org.Hs.eg.db",   # human gene annotation (Ensembl / symbol / Entrez mapping)
  "clusterProfiler",# pathway enrichment (GO / KEGG)
  "ReactomePA",     # Reactome pathway analysis
  "AnnotationDbi"   # gene ID mapping utilities
), update = FALSE, ask = FALSE)

install.packages(c(
  "tidyverse",    # data wrangling + ggplot2
  "ggrepel",      # clean gene labels on plots
  "pheatmap",     # heatmaps
  "RColorBrewer", # colour palettes
  "pROC",         # ROC / AUC curves
  "caret",        # ML framework
  "glmnet",       # LASSO / Ridge regression
  "randomForest", # random forest classification
  "e1071"         # SVM models
))


# ── 2. Load libraries ──────────────────────────────────────────────────────────

library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplot2)
library(randomForest)


# ── 3. Load and rename datasets ────────────────────────────────────────────────

als        <- ALS_MN.samples.counts.processed_GSE304550
parkinsons <- Parkinsons_RNA_seq_Norm_Counts_GSE299913

# Fix headers — column names loaded as first row
colnames(als) <- als[1, ]
als           <- als[-1, ]

# Rename Parkinson's columns to short readable labels
colnames(parkinsons) <- c("gene", "Tri1", "Tri2", "Tri3", "Iso1", "Iso2", "Iso3")

dim(als)
dim(parkinsons)


# =============================================================================
# ALS PIPELINE
# Dataset: GSE304550
# Comparison: SCR (scrambled control) vs TDP (TDP-43 knockdown — ALS model)
# Samples: 8 (4 control, 4 ALS) across 2 donor batches (M85, M86)
# =============================================================================


# ── 4. ALS EDA ────────────────────────────────────────────────────────────────

# Convert count columns from character to numeric
als[, 4:19] <- sapply(als[, 4:19], as.numeric)

# Distribution of counts per sample
# All boxes should sit at similar height — flags any failed sequencing runs
boxplot(log2(als[, 4:19] + 1),
        las  = 2,
        main = "ALS - Count Distribution (all 16 samples)",
        ylab = "Log2(count + 1)",
        col  = "steelblue")

# Library size per sample
# Total reads per sample — should be in a comparable range
colSums(als[, 4:19])
barplot(colSums(als[, 4:19]) / 1e6,
        las  = 2,
        main = "ALS - Library Sizes",
        ylab = "Total counts (millions)",
        col  = "steelblue")

# PCA — all 16 samples
# Checks for outliers and batch effects before filtering
als_pca <- prcomp(t(log2(als[, 4:19] + 1)))
plot(als_pca$x[, 1], als_pca$x[, 2],
     main = "ALS PCA — All 16 Samples",
     xlab = "PC1",
     ylab = "PC2",
     pch  = 19,
     col  = "steelblue")
text(als_pca$x[, 1], als_pca$x[, 2],
     labels = colnames(als[, 4:19]),
     pos    = 3,
     cex    = 0.7)

# Filter to SCR and TDP only — remove UPF and TU (side experiments)
als_filtered <- als[, !grepl("UPF|TU", colnames(als))]
dim(als_filtered)

# PCA — filtered 8 samples
# Green = SCR (control), Red = TDP (ALS)
# Checks class separation and confirms batch effect before limma
als_pca_filtered <- prcomp(t(log2(als_filtered[, 4:11] + 1)))
colours_als      <- c("green", "green", "red", "red", "green", "green", "red", "red")

plot(als_pca_filtered$x[, 1], als_pca_filtered$x[, 2],
     main = "ALS PCA — Filtered (SCR vs TDP)",
     xlab = "PC1",
     ylab = "PC2",
     pch  = 19,
     col  = colours_als)
text(als_pca_filtered$x[, 1], als_pca_filtered$x[, 2],
     labels = colnames(als_filtered[, 4:11]),
     pos    = 3,
     cex    = 0.7)
legend("topright",
       legend = c("Control (SCR)", "ALS (TDP)"),
       col    = c("green", "red"),
       pch    = 19)


# ── 5. ALS — Voom + Limma ─────────────────────────────────────────────────────

# Design matrix — accounts for donor batch effect (M85 vs M86)
# limma will test for SCR vs TDP difference while controlling for donor
als_group  <- factor(c("SCR", "SCR", "TDP", "TDP",
                       "SCR", "SCR", "TDP", "TDP"))
als_donor  <- factor(c("M85", "M85", "M85", "M85",
                       "M86", "M86", "M86", "M86"))
als_design <- model.matrix(~ als_donor + als_group)

# Voom — stabilises raw count variance for linear modelling
# Converts counts to logCPM and estimates precision weights per gene
als_v <- voom(counts = als_filtered[, 4:11],
              design = als_design,
              plot   = TRUE)

# Fit linear model, apply empirical Bayes shrinkage, extract results
als_fit     <- lmFit(als_v, als_design)
als_fit     <- eBayes(als_fit)
als_results <- topTable(als_fit, coef = "als_groupTDP", number = Inf)
head(als_results)

# Filter DEGs — significant (adj.P.Val < 0.05) and biologically meaningful (|logFC| > 1)
als_degs <- als_results[als_results$adj.P.Val < 0.05 & abs(als_results$logFC) > 1, ]
nrow(als_degs)
table(als_degs$logFC > 0)

# Attach Ensembl IDs back to DEG results (row numbers used as index)
als_degs$ensembl_id <- als_filtered$ensembl_gene_id[as.numeric(rownames(als_degs))]


# ── 6. ALS — Visualisation ────────────────────────────────────────────────────

# Volcano plot — all genes coloured by significance and direction
# Red = upregulated in ALS, Blue = downregulated, Grey = not significant
plot(als_results$logFC, -log10(als_results$adj.P.Val),
     main = "ALS — Volcano Plot",
     xlab = "log2 Fold Change",
     ylab = "-log10 adjusted p-value",
     pch  = 19,
     cex  = 0.4,
     col  = ifelse(als_results$adj.P.Val < 0.05 & als_results$logFC > 1,  "red",
                   ifelse(als_results$adj.P.Val < 0.05 & als_results$logFC < -1, "blue",
                          "grey")))
abline(h = -log10(0.05), lty = 2)
abline(v = c(-1, 1),     lty = 2)
legend("topright",
       legend = c("Upregulated", "Downregulated", "Not significant"),
       col    = c("red", "blue", "grey"),
       pch    = 19)

# Heatmap — top 50 DEGs
# Confirms samples cluster by group based on gene expression patterns
top50      <- head(als_degs[order(als_degs$adj.P.Val), ], 50)
top50_expr <- als_v$E[rownames(top50), ]

pheatmap(top50_expr,
         scale        = "row",
         main         = "ALS — Top 50 DEGs",
         show_rownames = FALSE,
         fontsize_col  = 8)


# ── 7. ALS — Pathway Enrichment ───────────────────────────────────────────────

# Convert Ensembl IDs to Entrez IDs — required by GO and KEGG databases
als_entrez         <- bitr(als_degs$ensembl_id,
                           fromType = "ENSEMBL",
                           toType   = "ENTREZID",
                           OrgDb    = org.Hs.eg.db)
als_entrez$ENTREZID <- as.character(als_entrez$ENTREZID)

# Note — strict FDR correction (BH) returns zero results due to small sample size
# n=4 per group is insufficient statistical power against 5880 GO terms
# Raw p-values used as exploratory indicator of enriched pathways
# Results should be interpreted cautiously and validated with larger datasets

# GO biological process enrichment
als_go_test <- enrichGO(gene         = als_entrez$ENTREZID,
                        OrgDb        = org.Hs.eg.db,
                        ont          = "BP",
                        pAdjustMethod = "none",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        readable     = TRUE)

# Filter by raw p-value
als_go_results <- as.data.frame(als_go_test)
als_go_sig     <- als_go_results[als_go_results$pvalue < 0.05, ]
nrow(als_go_sig)

# GO dotplot — top 20 enriched pathways
als_go_top20 <- als_go_sig[order(als_go_sig$pvalue), ][1:20, ]

ggplot(als_go_top20, aes(x = Count, y = reorder(Description, Count), colour = pvalue)) +
  geom_point(size = 4) +
  scale_colour_gradient(low = "red", high = "blue") +
  labs(title = "ALS — Top 20 Enriched GO Pathways",
       x     = "Gene Count",
       y     = "Pathway") +
  theme_minimal()

# KEGG pathway enrichment
als_kegg    <- enrichKEGG(gene          = als_entrez$ENTREZID,
                          organism      = "hsa",
                          pAdjustMethod = "none",
                          pvalueCutoff  = 0.05)
als_kegg_df <- as.data.frame(als_kegg)
nrow(als_kegg_df)

# KEGG bar plot
ggplot(als_kegg_df, aes(x = Count, y = reorder(Description, Count), fill = pvalue)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "ALS — KEGG Pathways",
       x     = "Gene Count",
       y     = "Pathway") +
  theme_minimal()


# =============================================================================
# PARKINSON'S PIPELINE
# Dataset: GSE299913
# Comparison: Triplication (SNCA x3 — Parkinson's model) vs Isogenic (corrected control)
# Samples: 6 (3 Triplication, 3 Isogenic) — pre-normalised, no voom needed
# =============================================================================


# ── 8. Parkinson's EDA ────────────────────────────────────────────────────────

# Convert count columns to numeric
parkinsons[, 2:7] <- sapply(parkinsons[, 2:7], as.numeric)

# Distribution of counts per sample
boxplot(log2(parkinsons[, 2:7] + 1),
        las  = 2,
        main = "Parkinson's — Count Distribution",
        ylab = "Log2(count + 1)",
        col  = "steelblue")

# PCA — Red = Triplication (disease), Green = Isogenic (control)
park_pca      <- prcomp(t(log2(parkinsons[, 2:7] + 1)))
colours_park  <- c("red", "red", "red", "green", "green", "green")

plot(park_pca$x[, 1], park_pca$x[, 2],
     main = "Parkinson's PCA",
     xlab = "PC1",
     ylab = "PC2",
     pch  = 19,
     col  = colours_park)
text(park_pca$x[, 1], park_pca$x[, 2],
     labels = colnames(parkinsons[, 2:7]),
     pos    = 3,
     cex    = 0.7)
legend("topright",
       legend = c("Triplication", "Isogenic"),
       col    = c("red", "green"),
       pch    = 19)


# ── 9. Parkinson's — Limma ────────────────────────────────────────────────────

# Design matrix — group only, no batch effect in this dataset
park_group  <- factor(c("Tri", "Tri", "Tri", "Iso", "Iso", "Iso"))
park_design <- model.matrix(~ park_group)

# Log2 transform — data is pre-normalised but not log transformed
# Required to meet limma's assumption of normally distributed values
park_counts             <- parkinsons[, 2:7]
rownames(park_counts)   <- parkinsons$gene
park_counts_log         <- log2(park_counts + 1)

# Fit linear model directly on log transformed values — no voom needed
park_fit     <- lmFit(park_counts_log, park_design)
park_fit     <- eBayes(park_fit)
park_results <- topTable(park_fit, coef = "park_groupTri", number = Inf)
head(park_results)

# Filter DEGs
park_degs <- park_results[park_results$adj.P.Val < 0.05 & abs(park_results$logFC) > 1, ]
nrow(park_degs)
table(park_degs$logFC > 0)


# ── 10. Parkinson's — Visualisation ───────────────────────────────────────────

# Volcano plot
plot(park_results$logFC, -log10(park_results$adj.P.Val),
     main = "Parkinson's — Volcano Plot",
     xlab = "log2 Fold Change",
     ylab = "-log10 adjusted p-value",
     pch  = 19,
     cex  = 0.4,
     col  = ifelse(park_results$adj.P.Val < 0.05 & park_results$logFC > 1,  "red",
                   ifelse(park_results$adj.P.Val < 0.05 & park_results$logFC < -1, "blue",
                          "grey")))
abline(h = -log10(0.05), lty = 2)
abline(v = c(-1, 1),     lty = 2)
legend("topright",
       legend = c("Upregulated", "Downregulated", "Not significant"),
       col    = c("red", "blue", "grey"),
       pch    = 19)

# Heatmap — top 50 DEGs
park_top50      <- head(park_degs[order(park_degs$adj.P.Val), ], 50)
park_top50_expr <- park_counts_log[rownames(park_top50), ]

pheatmap(park_top50_expr,
         scale         = "row",
         main          = "Parkinson's — Top 50 DEGs",
         show_rownames = FALSE,
         fontsize_col  = 8)


# ── 11. Parkinson's — Pathway Enrichment ──────────────────────────────────────

# Parkinson's uses gene symbols — convert directly to Entrez IDs
park_entrez          <- bitr(rownames(park_degs),
                             fromType = "SYMBOL",
                             toType   = "ENTREZID",
                             OrgDb    = org.Hs.eg.db)
park_entrez$ENTREZID <- as.character(park_entrez$ENTREZID)
nrow(park_entrez)

# GO biological process enrichment — strict BH correction passes here
# Stronger signal from SNCA triplication vs TDP-43 knockdown
park_go_test <- enrichGO(gene          = park_entrez$ENTREZID,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         readable      = TRUE)
nrow(park_go_test)

# GO dotplot — top 20 pathways
park_go_results <- as.data.frame(park_go_test)
park_go_top20   <- park_go_results[order(park_go_results$p.adjust), ][1:20, ]

ggplot(park_go_top20, aes(x = Count, y = reorder(Description, Count), colour = p.adjust)) +
  geom_point(size = 4) +
  scale_colour_gradient(low = "red", high = "blue") +
  labs(title = "Parkinson's — Top 20 Enriched GO Pathways",
       x     = "Gene Count",
       y     = "Pathway") +
  theme_minimal()

# KEGG pathway enrichment
park_kegg    <- enrichKEGG(gene          = park_entrez$ENTREZID,
                           organism      = "hsa",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05)
park_kegg_df <- as.data.frame(park_kegg)
nrow(park_kegg_df)

# KEGG bar plot — top 20
park_kegg_top20 <- park_kegg_df[order(park_kegg_df$p.adjust), ][1:20, ]

ggplot(park_kegg_top20, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "Parkinson's — KEGG Pathways",
       x     = "Gene Count",
       y     = "Pathway") +
  theme_minimal()


# =============================================================================
# CROSS-DISEASE COMPARISON
# Identifies shared DEGs and pathways between ALS and Parkinson's
# =============================================================================


# ── 12. DEG overlap ───────────────────────────────────────────────────────────

# Convert ALS Ensembl IDs to gene symbols for comparison
als_symbols <- bitr(als_degs$ensembl_id,
                    fromType = "ENSEMBL",
                    toType   = "SYMBOL",
                    OrgDb    = org.Hs.eg.db)

# Find shared DEGs
shared_genes <- intersect(als_symbols$SYMBOL, rownames(park_degs))
length(shared_genes)
head(shared_genes)

# Find shared KEGG pathways
shared_pathways <- intersect(als_kegg_df$Description, park_kegg_df$Description)
length(shared_pathways)
shared_pathways

# Visualise DEG overlap
venn_data <- data.frame(
  disease = c("ALS only", "Shared", "Parkinson's only"),
  count   = c(length(als_symbols$SYMBOL) - length(shared_genes),
              length(shared_genes),
              length(rownames(park_degs)) - length(shared_genes))
)

ggplot(venn_data, aes(x = disease, y = count, fill = disease)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("ALS only"         = "red",
                               "Shared"           = "purple",
                               "Parkinson's only" = "green")) +
  labs(title = "DEG Overlap — ALS vs Parkinson's",
       x     = "",
       y     = "Number of genes") +
  theme_minimal()


# =============================================================================
# MACHINE LEARNING
# Random forest classifiers trained separately for each disease
# Purpose — validate that DEG signatures are predictive of disease vs control
# =============================================================================


# ── 13. ML — ALS ──────────────────────────────────────────────────────────────

# Extract DEG expression values from voom output
# Transpose so samples are rows and genes are columns — ML format
als_ml_matrix <- t(als_v$E[rownames(als_degs), ])
als_ml_labels  <- factor(c("SCR", "SCR", "TDP", "TDP",
                           "SCR", "SCR", "TDP", "TDP"))
dim(als_ml_matrix)

# Scale — normalises each gene to mean 0, SD 1
# Prevents high expression genes dominating the model
als_ml_scaled <- scale(als_ml_matrix)

# Train/test split — 6 train, 2 test
set.seed(42)
als_train_idx <- sample(1:8, 6)
als_test_idx  <- setdiff(1:8, als_train_idx)

als_train_x <- als_ml_scaled[als_train_idx, ]
als_test_x  <- als_ml_scaled[als_test_idx, ]
als_train_y <- als_ml_labels[als_train_idx]
als_test_y  <- als_ml_labels[als_test_idx]

table(als_train_y)

# Train random forest
set.seed(42)
als_rf <- randomForest(x         = als_train_x,
                       y         = als_train_y,
                       ntree     = 500,
                       importance = TRUE)

# Evaluate
als_pred <- predict(als_rf, als_test_x)
als_pred
als_test_y

# Feature importance — which genes drove classification
als_importance  <- importance(als_rf)
als_top_genes   <- head(als_importance[order(als_importance[, "MeanDecreaseAccuracy"],
                                             decreasing = TRUE), ], 10)

# Convert row numbers to gene symbols
als_top_gene_ids   <- als_filtered$ensembl_gene_id[as.numeric(rownames(als_top_genes))]
als_top_gene_names <- bitr(als_top_gene_ids,
                           fromType = "ENSEMBL",
                           toType   = "SYMBOL",
                           OrgDb    = org.Hs.eg.db)

# Plot feature importance
als_importance_df <- data.frame(
  gene       = als_top_gene_names$SYMBOL,
  importance = als_top_genes[rownames(als_top_genes) %in%
                               match(als_top_gene_names$ENSEMBL,
                                     als_filtered$ensembl_gene_id),
                             "MeanDecreaseAccuracy"]
)

ggplot(als_importance_df, aes(x = importance,
                              y = reorder(gene, importance),
                              fill = importance)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "steelblue", high = "red") +
  labs(title = "ALS — Top 10 Most Important Genes",
       x     = "Mean Decrease Accuracy",
       y     = "Gene") +
  theme_minimal()


# ── 14. ML — Parkinson's ──────────────────────────────────────────────────────

# Extract DEG expression values
park_ml_matrix <- t(park_counts_log[rownames(park_degs), ])
park_ml_labels  <- factor(c("Tri", "Tri", "Tri", "Iso", "Iso", "Iso"))
dim(park_ml_matrix)

# Scale
park_ml_scaled <- scale(park_ml_matrix)

# Balanced train/test split — manually ensures 2 Tri and 2 Iso in training
# Random split gave 3 Iso and 1 Tri which produced unreliable feature importance
park_train_idx <- c(1, 2, 4, 5)  # Tri1, Tri2, Iso1, Iso2
park_test_idx  <- c(3, 6)         # Tri3, Iso3

park_train_x <- park_ml_scaled[park_train_idx, ]
park_test_x  <- park_ml_scaled[park_test_idx, ]
park_train_y <- park_ml_labels[park_train_idx]
park_test_y  <- park_ml_labels[park_test_idx]

table(park_train_y)

# Train random forest
set.seed(42)
park_rf <- randomForest(x         = park_train_x,
                        y         = park_train_y,
                        ntree     = 500,
                        importance = TRUE)

# Evaluate
park_pred <- predict(park_rf, park_test_x)
park_pred
park_test_y

# Feature importance
park_importance <- importance(park_rf)
park_top_genes  <- head(park_importance[order(park_importance[, "MeanDecreaseAccuracy"],
                                              decreasing = TRUE), ], 10)

# Plot feature importance
park_importance_df <- data.frame(
  gene       = rownames(park_top_genes),
  importance = park_top_genes[, "MeanDecreaseAccuracy"]
)

ggplot(park_importance_df, aes(x = importance,
                               y = reorder(gene, importance),
                               fill = importance)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "steelblue", high = "red") +
  labs(title = "Parkinson's — Top 10 Most Important Genes",
       x     = "Mean Decrease Accuracy",
       y     = "Gene") +
  theme_minimal()


# ── 15. Compare ML feature importance across both diseases ────────────────────

# Check if any top predictive genes are shared between ALS and Parkinson's models
als_top_names          <- als_top_gene_names$SYMBOL
park_top_names         <- rownames(park_top_genes)
shared_important_genes <- intersect(als_top_names, park_top_names)

length(shared_important_genes)
shared_important_genes
