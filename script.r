# ================================
# 1. Install & Load Libraries
# ================================
if (!requireNamespace("BiocManager")) install.packages("BiocManager")

packages <- c("edgeR", "limma", "EnhancedVolcano", "pheatmap", "RColorBrewer", "dplyr")
BiocManager::install(packages, ask = FALSE)

lapply(packages, library, character.only = TRUE)

# ================================
# 2. Load Data
# ================================
setwd("~/MS Bioinfo/Semester 2/NGS/EdgeR_data_Project2/Dataset for Project 2")

metadata <- read.delim("metadata.txt")
counts   <- read.delim("Counts.txt", row.names = 1)

# Clean group labels
metadata$Category <- recode(metadata$Category,
                            "Crohn's disease" = "CD",
                            "non inflammatory bowel disease control" = "Control")

metadata$Category <- factor(metadata$Category, levels = c("Control", "CD"))

# ================================
# 3. Create DGE Object
# ================================
y <- DGEList(counts = counts, group = metadata$Category)

# Filter low expression genes
y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]

# Normalize
y <- calcNormFactors(y)

# ================================
# 4. MDS Plot (Points Only)
# ================================

plotMDS(y,
        col = c("blue", "red")[y$samples$group],
        pch = 16,          # solid circles
        cex = 1.5,         # increase point size
        labels = NULL,     # remove labels
        main = "MDS Plot")

legend("topright",
       legend = levels(y$samples$group),
       col = c("blue", "red"),
       pch = 16)

# ================================
# 5. Differential Expression
# ================================
design <- model.matrix(~ Category, data = metadata)

y   <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)

# Results table
res <- topTags(qlf, n = Inf)$table
res <- res[order(res$FDR), ]

write.table(res, "edgeR_DGE_Results.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# Define cutoffs
p_cutoff  <- 1e-20
fc_cutoff <- 5.5

# ================================
# 6. Volcano Plot
# ================================
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'logFC',
                y = 'PValue',
                pCutoff = p_cutoff,
                FCcutoff = fc_cutoff,
                
                title = 'Crohns Disease vs Control',
                subtitle = sprintf("Cutoffs: |log2FC| > %.1f & P-value < %.1e",
                                   fc_cutoff, p_cutoff),
                
                col = c("grey70", "steelblue", "orange", "red2"),
                colAlpha = 0.6,
                pointSize = 1.8,
                labSize = 4,
                legendPosition = "right")

# ================================
# 7. Heatmap of Top 30 DEGs
# ================================

# Get logCPM
logCPM <- cpm(y, log = TRUE)

# Select top 30 DEGs by FDR
top30_genes <- rownames(res)[1:30]

# Subset expression matrix
heatmap_data <- logCPM[top30_genes, ]

# Scale by row (important for visualization)
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# Annotation (group info)
annotation_col <- data.frame(Group = metadata$Category)
rownames(annotation_col) <- colnames(heatmap_data_scaled)

# Colors
ann_colors <- list(Group = c(Control = "blue", CD = "red"))

# Plot heatmap
pheatmap(heatmap_data_scaled,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 6,
         main = "Top 30 Differentially Expressed Genes")
