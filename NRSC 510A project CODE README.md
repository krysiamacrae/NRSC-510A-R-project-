#Set working directory 
setwd("~/Desktop/NRSC 510A/NRSC 510A - project ")

# Install necessary packages
install.packages("Seurat")
install.packages("dplyr")

# Load libraries
library(Seurat)
library(dplyr)

# Define paths to the control and stressed data files
control_matrix <- "/Users/krystynamacrae/Desktop/NRSC 510A/NRSC 510A - project /E15.5 microglia from CD1 embryos /GSM3901919_sample_1_matrix.mtx"
control_features <- "/Users/krystynamacrae/Desktop/NRSC 510A/NRSC 510A - project /E15.5 microglia from CD1 embryos /GSM3901919_sample_1_features.tsv"
control_barcodes <- "/Users/krystynamacrae/Desktop/NRSC 510A/NRSC 510A - project /E15.5 microglia from CD1 embryos /GSM3901919_sample_1_barcodes.tsv"

stressed_matrix <- "/Users/krystynamacrae/Desktop/NRSC 510A/NRSC 510A - project /E15.5 microglia from stressed CD1 embryos/GSM3901920_sample_2_matrix.mtx"
stressed_features <- "/Users/krystynamacrae/Desktop/NRSC 510A/NRSC 510A - project /E15.5 microglia from stressed CD1 embryos/GSM3901920_sample_2_features.tsv"
stressed_barcodes <- "/Users/krystynamacrae/Desktop/NRSC 510A/NRSC 510A - project /E15.5 microglia from stressed CD1 embryos/GSM3901920_sample_2_barcodes.tsv"

# Load the control dataset
control_data <- ReadMtx(mtx = control_matrix, features = control_features, cells = control_barcodes)

# Load the stressed dataset
stressed_data <- ReadMtx(mtx = stressed_matrix, features = stressed_features, cells = stressed_barcodes)

# Create Seurat objects for control and stressed conditions
control_seurat <- CreateSeuratObject(counts = control_data, project = "Control")
stressed_seurat <- CreateSeuratObject(counts = stressed_data, project = "Stressed")

# Add metadata for condition
control_seurat$condition <- "Control"
stressed_seurat$condition <- "Stressed"

combined_seurat <- merge(stressed_seurat, y = control_seurat, 
                         add.cell.ids = c("Stressed", "Control"), 
                         project = "UMAP replicaton")

unique(combined_seurat$condition)

# Calculate the percentage of mitochondrial genes
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^mt-")

# Visualize QC metrics to decide on filtering thresholds
VlnPlot(combined_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter out cells with:
# - Fewer than 200 genes
# - More than 5000 genes (indicative of potential doublets)
# - More than 10% mitochondrial gene expression

combined_seurat <- subset(combined_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Create a scatter plot of nCount_RNA (total UMI count) vs percent.mt (mitochondrial percentage)
plot1 <- FeatureScatter(combined_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")

# Create a scatter plot of nCount_RNA (total UMI count) vs nFeature_RNA (number of detected genes)
plot2 <- FeatureScatter(combined_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Combine and visualize both scatter plots
plot1 + plot2

# Normalize the combined dataset
combined_seurat <- NormalizeData(combined_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Normalize the combined dataset
# combined_seurat <- NormalizeData(combined_seurat)

# Identify highly variable features
combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 2000)

# Plot the variable features
top10 <- head(VariableFeatures(combined_seurat), 10)

# Plot variable features with the mean expression and dispersion
plot1 <- VariableFeaturePlot(combined_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Display the plots
plot2

# Get all gene names (row names of the Seurat object)
all.genes <- rownames(combined_seurat)

# Scale the data for all genes
combined_seurat <- ScaleData(combined_seurat, features = all.genes)

# Perform PCA on the scaled data
combined_seurat <- RunPCA(combined_seurat, features = VariableFeatures(object = combined_seurat))

# Visualizing PCA results in a few ways 

# Check the top genes contributing to the first few principal components
print(combined_seurat[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize the top genes that contribute to the first two principal components
VizDimLoadings(combined_seurat, dims = 1:2, reduction = "pca")

# Plot the first two PCs to visualize cell separation
DimPlot(combined_seurat, reduction = "pca") + NoLegend()

# Visualize the top genes defining the first principal component
DimHeatmap(combined_seurat, dims = 1, cells = 500, balanced = TRUE)

# Visualize the top genes defining multiple principal components (PC1 to PC6)
DimHeatmap(combined_seurat, dims = 1:6, cells = 500, balanced = TRUE)

# Create an elbow plot to visualize the variance explained by each PC
ElbowPlot(combined_seurat)

# Construct the SNN graph using the chosen number of PCs (adjust num_pcs based on the previous step)
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:8)

# Perform clustering using the SNN graph
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

# View the cluster IDs of the first 5 cells 
head(Idents(combined_seurat))

# Visualize the clusters using UMAP 
combined_seurat <- RunUMAP(combined_seurat, dims = 1:8)
DimPlot(combined_seurat, reduction = "umap")

# Visualize the clusters using UMAP 
combined_seurat <- RunUMAP(combined_seurat, dims = 1:8)
DimPlot(combined_seurat, reduction = "umap", group.by = "condition")

#FINDING UPREGULARTED GENES IN DR. ROSIN CLUSTER 3 IN MY MERGED UMAP 

# Replace "GeneA", "GeneB", etc. with actual gene names enriched in cluster 3 from the paper
FeaturePlot(combined_seurat, features = c("Ube2c", "Hist1h2ap", "Birc5", "Hmgb2", "H2afx", "Cenpa", "Mki67", "Pclaf", "Top2a", "Pttg1"), cols = c("lightgrey", "blue"))


#SUBSET CLUSTER 3

# List all unique cluster identities to ensure "3" is one of them
unique(Idents(combined_seurat))

# Subset cluster 3 from the combined Seurat object
cluster_3 <- subset(combined_seurat, idents = 3)

# Run UMAP on the subset if you want a separate visualization
cluster_3 <- RunUMAP(cluster_3, dims = 1:15)
DimPlot(cluster_3, reduction = "umap",  group.by = "condition")

# Run UMAP on the subset if you want a separate visualization
cluster_3 <- RunUMAP(cluster_3, dims = 1:8)
DimPlot(cluster_3, reduction = "umap")

#SUBSET CLUSTER 5

# Subset cluster 5 from the combined Seurat object
cluster_5 <- subset(combined_seurat, idents = 5)

# Run UMAP on the subset if you want a separate visualization
cluster_5 <- RunUMAP(cluster_5, dims = 1:8)
DimPlot(cluster_5, reduction = "umap",  group.by = "condition")

#SUBSET CLUSTER 3 AND 5

# Subset clusters 3 and 5 from the combined Seurat object
clusters_3_and_5 <- subset(combined_seurat, idents = c(3, 5))

# Run UMAP on the subset if desired
clusters_3_and_5 <- RunUMAP(clusters_3_and_5, dims = 1:8)
DimPlot(clusters_3_and_5, reduction = "umap", group.by = "condition")


#WORSE RESOLUTION AND LOOK FOR DR.ROSIN'S CLUSTER 3

# Perform clustering using the SNN graph
combined_seurat <- FindClusters(combined_seurat, resolution = 0.1)

# View the cluster IDs of the first 5 cells 
head(Idents(combined_seurat))

# Visualize the clusters using UMAP 
combined_seurat <- RunUMAP(combined_seurat, dims = 1:8)
DimPlot(combined_seurat, reduction = "umap")

# Visualize the clusters using UMAP + condition 
combined_seurat <- RunUMAP(combined_seurat, dims = 1:8)
DimPlot(combined_seurat, reduction = "umap", group.by = "condition")

# Replace "GeneA", "GeneB", etc. with actual gene names enriched in cluster 3 from the paper
FeaturePlot(combined_seurat, features = c("Ube2c", "Hist1h2ap", "Birc5", "Hmgb2", "H2afx", "Cenpa", "Mki67", "Pclaf", "Top2a", "Pttg1"), cols = c("lightgrey", "blue"))

#SUBSET CLUSTER 0 (after lower resolution)

# List all unique cluster identities to ensure "0" is one of them
unique(Idents(combined_seurat))

# Subset cluster 0 from the combined Seurat object
cluster_0 <- subset(combined_seurat, idents = 0)

# Run UMAP on the subset if you want a separate visualization
cluster_0 <- RunUMAP(cluster_0, dims = 1:15)
DimPlot(cluster_0, reduction = "umap",  group.by = "condition")

# Run UMAP on the subset if you want a separate visualization
cluster_0 <- RunUMAP(cluster_0, dims = 1:15)
DimPlot(cluster_0, reduction = "umap")


#SUBSET CLUSTER 1 (after lower resolution)

# List all unique cluster identities to ensure "1" is one of them
unique(Idents(combined_seurat))

# Subset cluster 1 from the combined Seurat object
cluster_1 <- subset(combined_seurat, idents = 1)

# Run UMAP on the subset if you want a separate visualization
cluster_1 <- RunUMAP(cluster_1, dims = 1:15)
DimPlot(cluster_1, reduction = "umap",  group.by = "condition")

# Run UMAP on the subset if you want a separate visualization
cluster_1 <- RunUMAP(cluster_1, dims = 1:15)
DimPlot(cluster_1, reduction = "umap")


#SUBSET CLUSTER 2 (after lower resolution)

# List all unique cluster identities to ensure "2" is one of them
unique(Idents(combined_seurat))

# Subset cluster 2 from the combined Seurat object
cluster_2 <- subset(combined_seurat, idents = 2)

# Run UMAP on the subset if you want a separate visualization
cluster_2 <- RunUMAP(cluster_2, dims = 1:15)
DimPlot(cluster_2, reduction = "umap",  group.by = "condition")

# Run UMAP on the subset if you want a separate visualization
cluster_2 <- RunUMAP(cluster_2, dims = 1:15)
DimPlot(cluster_2, reduction = "umap")


#SUBSET CLUSTER 3 (after lower resolution)

# List all unique cluster identities to ensure "3" is one of them
unique(Idents(combined_seurat))

# Subset cluster 3 from the combined Seurat object
cluster_3 <- subset(combined_seurat, idents = 3)

# Run UMAP on the subset if you want a separate visualization
cluster_3 <- RunUMAP(cluster_3, dims = 1:15)
DimPlot(cluster_3, reduction = "umap",  group.by = "condition")

# Run UMAP on the subset if you want a separate visualization
cluster_3 <- RunUMAP(cluster_3, dims = 1:15)
DimPlot(cluster_3, reduction = "umap")



#RECLUSTERRING CLUSTER 3

# Identify highly variable features
cluster_3 <- FindVariableFeatures(cluster_3, selection.method = "vst", nfeatures = 2000)

# Plot the variable features
top10 <- head(VariableFeatures(cluster_3), 10)

# Plot variable features with the mean expression and dispersion
plot1 <- VariableFeaturePlot(cluster_3)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Display the plots
plot2

# Get all gene names (row names of the Seurat object)
all.genes <- rownames(cluster_3)

# Scale the data for all genes
cluster_3 <- ScaleData(cluster_3, features = all.genes)

# Perform PCA on the scaled data 
cluster_3 <- RunPCA(cluster_3, features = VariableFeatures(object = cluster_3))

# Visualizing PCA results in a few ways 

# Check the top genes contributing to the first few principal components
print(cluster_3[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize the top genes that contribute to the first two principal components
VizDimLoadings(cluster_3, dims = 1:2, reduction = "pca")

# Plot the first two PCs to visualize cell separation
DimPlot(cluster_3, reduction = "pca") + NoLegend()

# Visualize the top genes defining the first principal component
DimHeatmap(cluster_3, dims = 1, cells = 500, balanced = TRUE)

# Visualize the top genes defining multiple principal components (PC1 to PC6)
DimHeatmap(cluster_3, dims = 1:6, cells = 500, balanced = TRUE)

# Create an elbow plot to visualize the variance explained by each PC
ElbowPlot(cluster_3)

# Construct the SNN graph using the chosen number of PCs (adjust num_pcs based on the previous step)
cluster_3 <- FindNeighbors(cluster_3, dims = 1:7)

# Perform clustering using the SNN graph
cluster_3 <- FindClusters(cluster_3, resolution = 0.5)

# View the cluster IDs of the first 5 cells 
head(Idents(cluster_3))

# Visualize the clusters using UMAP 
cluster_3 <- RunUMAP(cluster_3, dims = 1:7)
DimPlot(cluster_3, reduction = "umap")

# Visualize the clusters using UMAP 
combined_seurat <- RunUMAP(cluster_3, dims = 1:7)
DimPlot(cluster_3, reduction = "umap")
DimPlot(cluster_3, reduction = "umap", group.by = "condition")

all_markers <- FindAllMarkers(cluster_3,
                              only.pos = TRUE,         
                              logfc.threshold = 0.25,  
                              test.use = "wilcox")     

#Find the top 50 marker genes in each cluster of cluster 3
top_genes <- all_markers %>%
  filter(pct.2 < 0.3) %>%
  group_by(cluster) %>%
  slice_max(order_by = pct.1, n = 50)

#Looking for the other marker genes from Dr. Rosin's paper in my UMAP

# Replace "GeneA", "GeneB", etc. with actual gene names enriched in cluster 1 from the paper
FeaturePlot(combined_seurat, features = c("Pf4", "F13a1", "Mrc1", "Lyve1", "Dab2", "Ccl24", "Ms4a7", "Ms4a6c", "Blvrb", "Fcgrt"), cols = c("lightgrey", "blue"))

# Replace "GeneA", "GeneB", etc. with actual gene names enriched in cluster 2 from the paper
FeaturePlot(combined_seurat, features = c("Spp1", "Fabp5", "Ccl4", "Ccl3", "Csf1", "Lpl", "Ctsb", "Fam20c", "Mif", "Fabp3"), cols = c("lightgrey", "blue"))

# Replace "GeneA", "GeneB", etc. with actual gene names enriched in cluster 4 from the paper
FeaturePlot(combined_seurat, features = c("Crybb1", "P2ry12", "Sparc", "Ccr5", "Tgfbr1", "Lpcat2", "Siglech", "C1qa", "Tmem119", "Sall1"), cols = c("lightgrey", "blue"))


#FINDING CLUSTER GENES 

#Most expressed in cluster 0

# Display genes from positions in the top_genes data frame
genes_1_to_10 <- top_genes$gene[1:10]
print(genes_1_to_10)

#Marker genes for cluster 0 (1-10, 0)
VlnPlot(cluster_3, features = genes_1_to_10)

#Most expressed in cluster 1

# Display genes from positions in the top_genes data frame
genes_52_to_65 <- top_genes$gene[52:65]
print(genes_52_to_65)

#Marker genes for cluster 1 (52-65, 1)
VlnPlot(cluster_3, features = genes_52_to_65)

#Most expressed in cluster 2

# Display genes from positions in the top_genes data frame
genes_102_to_120 <- top_genes$gene[102:120]
print(genes_102_to_120)

#Marker genes for cluster 2 (102-120, 2)
VlnPlot(cluster_3, features = genes_102_to_120)


#PATHWAY ENRICHEMENT ANALYSIS (CLUSTER 3)

# Install devtools if not already installed
install.packages("devtools")

# Install enrichR from GitHub
devtools::install_github("wjawaid/enrichR")

library(enrichR)

# List available databases for enrichment analysis
enrichr_databases <- listEnrichrDbs()

#ident 0

#GO

# Select the databases to use (e.g., "GO_Biological_Process_2021", "KEGG_2021_Human")
selected_databases <- c('GO_Biological_Process_2023')

# Run DEenrichRPlot
results_GO_ident0 <- DEenrichRPlot(
  object = cluster_3,               # The Seurat object
  ident.1 = 0,           # Specify the identity group (e.g., control)
  ident.2 = NULL,           # Specify the comparison group (e.g., treated)
  balanced = FALSE,                  # Ensure equal up- and down-regulated genes
  logfc.threshold = 0.25,           # Threshold for log fold change
  assay = "RNA",                    # Assay to use
  max.genes = 50,                   # Max number of genes for enrichment
  test.use = "wilcox",              # Statistical test for DE
  p.val.cutoff = 0.05,         # P-value cutoff for DE genes
  cols = NULL,
  enrich.database = selected_databases,  # The databases to use for enrichment
  num.pathway = 10,             # Number of pathways to display
  return.gene.list = TRUE,
)

# Extract the dataframe from `results`
pos_results_GO_ident0 <- results_GO_ident0$pos  # Assuming "pos" contains the relevant data

# Filter rows where the p-value column (`KEGG_2019_Mouse.P.value`) is less than 0.05
filtered_results_GO_ident0 <- pos_results_GO_ident0[pos_results_GO_ident0$GO_Biological_Process_2023.Adjusted.P.value < 0.05, ]

# View the filtered results
print(filtered_results_GO_ident0)

#DOT PLOT

# Convert adjusted p-values to -log10 for visualization
filtered_results_GO_ident0$log10_pvalue <- -log10(filtered_results_GO_ident0$GO_Biological_Process_2023.Adjusted.P.value)

# Create the dot plot
library(ggplot2)

ggplot(filtered_results_GO_ident0, aes(x = GO_Biological_Process_2023.Odds.Ratio, 
                                       y = reorder(GO_Biological_Process_2023.Term, GO_Biological_Process_2023.Odds.Ratio))) +
  geom_point(aes(size = GO_Biological_Process_2023.Odds.Ratio, color = log10_pvalue)) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(P-value)") +
  theme_minimal() +
  labs(
    title = "GO Biological Process Enrichment Analysis",
    x = "Odds Ratio",
    y = "Pathway",
    size = "Odds Ratio"
  )


#BAR PLOT

# Create the bar plot
ggplot(filtered_results_GO_ident0, aes(x = reorder(GO_Biological_Process_2023.Term, log10_pvalue), y = log10_pvalue)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Flip coordinates to make bars horizontal
  theme_minimal() +
  labs(
    title = "GO Biological Process Enrichment (Top Pathways)",
    x = "Pathway",
    y = "-log10(Adjusted P-value)"
  )



#KEGG
selected_databases <- c('KEGG_2019_Mouse')

# Run DEenrichRPlot
results_kegg_ident0 <- DEenrichRPlot(
  object = cluster_3,               # The Seurat object
  ident.1 = 0,           # Specify the identity group (e.g., control)
  ident.2 = NULL,           # Specify the comparison group (e.g., treated)
  balanced = FALSE,                  # Ensure equal up- and down-regulated genes
  logfc.threshold = 0.25,           # Threshold for log fold change
  assay = "RNA",                    # Assay to use
  max.genes = 50,                   # Max number of genes for enrichment
  test.use = "wilcox",              # Statistical test for DE
  p.val.cutoff = 0.05,         # P-value cutoff for DE genes
  cols = NULL,
  enrich.database = selected_databases,  # The databases to use for enrichment
  num.pathway = 10,             # Number of pathways to display
  return.gene.list = FALSE,
)


# Extract the dataframe from `results`
data_results_kegg_ident0 <- results_kegg_ident0$data  # Assuming "pos" contains the relevant data

# Filter rows where the p-value column (`KEGG_2019_Mouse.P.value`) is less than 0.05
filtered_results_kegg_ident0 <- data_results_kegg_ident0[data_results_kegg_ident0$KEGG_2019_Mouse.Adjusted.P.value < 0.05, ]

# View the filtered results
print(filtered_results_kegg_ident0)

#DOT PLOT

# Convert adjusted p-values to -log10 for color mapping
filtered_results_kegg_ident0$log10_adjusted_pvalue <- -log10(filtered_results_kegg_ident0$KEGG_2019_Mouse.Adjusted.P.value)

# Create the dot plot for KEGG enrichment analysis
library(ggplot2)

ggplot(filtered_results_kegg_ident0, aes(
  x = KEGG_2019_Mouse.Odds.Ratio, 
  y = reorder(KEGG_2019_Mouse.Term, KEGG_2019_Mouse.Odds.Ratio)
)) +
  geom_point(aes(size = KEGG_2019_Mouse.Odds.Ratio, 
                 color = log10_adjusted_pvalue)) +
  scale_color_gradient(low = "blue", high = "red") +  # More significant pathways in red
  theme_minimal() +
  labs(
    title = "KEGG Pathway Enrichment Analysis",
    x = "Odds Ratio",
    y = "Pathway",
    color = "-log10(Adjusted P-value)",
    size = "Odds Ratio"
  )


#BAR PLOT

# Convert adjusted p-values to -log10 for color mapping
filtered_results_kegg_ident0$log10_adjusted_pvalue <- -log10(filtered_results_kegg_ident0$KEGG_2019_Mouse.Adjusted.P.value)

# Create the bar plot for KEGG enrichment analysis
ggplot(filtered_results_kegg_ident0, aes(
  x = reorder(KEGG_2019_Mouse.Term, log10_adjusted_pvalue), 
  y = log10_adjusted_pvalue
)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Flip coordinates to make bars horizontal
  theme_minimal() +
  labs(
    title = "KEGG Pathway Enrichment (Top Pathways)",
    x = "Pathway",
    y = "-log10(Adjusted P-value)"
  )

#ident 1

#GO

selected_databases <- c('GO_Biological_Process_2023')

# Run DEenrichRPlot
results_GO_ident1 <- DEenrichRPlot(
  object = cluster_3,               
  ident.1 = 1,          
  ident.2 = NULL,           
  balanced = FALSE,                  
  logfc.threshold = 0.25,           
  assay = "RNA",                   
  max.genes = 50,                   
  test.use = "wilcox",              
  p.val.cutoff = 0.05,         
  cols = NULL,
  enrich.database = selected_databases,  
  num.pathway = 10,             
  return.gene.list = TRUE,
)

# Extract the dataframe from `results`
pos_results_GO_ident1 <- results_GO_ident1$pos  # Assuming "pos" contains the relevant data

# Filter rows where the p-value column (`KEGG_2019_Mouse.P.value`) is less than 0.05
filtered_results_GO_ident1 <- pos_results_GO_ident1[pos_results_GO_ident1$GO_Biological_Process_2023.Adjusted.P.value < 0.05, ]

# View the filtered results
print(filtered_results_GO_ident1)

#DOT PLOT

# Convert adjusted p-values to -log10 for visualization
filtered_results_GO_ident1$log10_pvalue <- -log10(filtered_results_GO_ident1$GO_Biological_Process_2023.Adjusted.P.value)

# Create the dot plot
library(ggplot2)

ggplot(filtered_results_GO_ident1, aes(x = GO_Biological_Process_2023.Odds.Ratio, 
                                       y = reorder(GO_Biological_Process_2023.Term, GO_Biological_Process_2023.Odds.Ratio))) +
  geom_point(aes(size = GO_Biological_Process_2023.Odds.Ratio, color = log10_pvalue)) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(P-value)") +
  theme_minimal() +
  labs(
    title = "GO Biological Process Enrichment Analysis",
    x = "Odds Ratio",
    y = "Pathway",
    size = "Odds Ratio"
  )


#BAR PLOT

# Create the bar plot
ggplot(filtered_results_GO_ident1, aes(x = reorder(GO_Biological_Process_2023.Term, log10_pvalue), y = log10_pvalue)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Flip coordinates to make bars horizontal
  theme_minimal() +
  labs(
    title = "GO Biological Process Enrichment (Top Pathways)",
    x = "Pathway",
    y = "-log10(Adjusted P-value)"
  )

cluster_3 <- JoinLayers(cluster_3)


#KEGG
selected_databases <- c('KEGG_2019_Mouse')

# Run DEenrichRPlot
results_kegg_ident1 <- DEenrichRPlot(
  object = cluster_3,               
  ident.1 = 1,           
  ident.2 = NULL,           
  balanced = FALSE,                  
  logfc.threshold = 0.25,           
  assay = "RNA",                   
  max.genes = 50,                   
  test.use = "wilcox",              
  p.val.cutoff = 0.05,         
  cols = NULL,
  enrich.database = selected_databases,  
  num.pathway = 10,             
  return.gene.list = FALSE,
)


# Extract the dataframe from `results`
data_results_kegg_ident1 <- results_kegg_ident1$data  # Assuming "pos" contains the relevant data

# Filter rows where the p-value column (`KEGG_2019_Mouse.P.value`) is less than 0.05
filtered_results_kegg_ident1 <- data_results_kegg_ident1[data_results_kegg_ident1$KEGG_2019_Mouse.Adjusted.P.value < 0.05, ]

# View the filtered results
print(filtered_results_kegg_ident1)

#DOT PLOT

# Convert adjusted p-values to -log10 for color mapping
filtered_results_kegg_ident1$log10_adjusted_pvalue <- -log10(filtered_results_kegg_ident1$KEGG_2019_Mouse.Adjusted.P.value)

# Create the dot plot for KEGG enrichment analysis
library(ggplot2)

ggplot(filtered_results_kegg_ident1, aes(
  x = KEGG_2019_Mouse.Odds.Ratio, 
  y = reorder(KEGG_2019_Mouse.Term, KEGG_2019_Mouse.Odds.Ratio)
)) +
  geom_point(aes(size = KEGG_2019_Mouse.Odds.Ratio, 
                 color = log10_adjusted_pvalue)) +
  scale_color_gradient(low = "blue", high = "red") +  # More significant pathways in red
  theme_minimal() +
  labs(
    title = "KEGG Pathway Enrichment Analysis",
    x = "Odds Ratio",
    y = "Pathway",
    color = "-log10(Adjusted P-value)",
    size = "Odds Ratio"
  )


#BAR PLOT

# Convert adjusted p-values to -log10 for color mapping
filtered_results_kegg_ident1$log10_adjusted_pvalue <- -log10(filtered_results_kegg_ident1$KEGG_2019_Mouse.Adjusted.P.value)

# Create the bar plot for KEGG enrichment analysis
ggplot(filtered_results_kegg_ident1, aes(
  x = reorder(KEGG_2019_Mouse.Term, log10_adjusted_pvalue), 
  y = log10_adjusted_pvalue
)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Flip coordinates to make bars horizontal
  theme_minimal() +
  labs(
    title = "KEGG Pathway Enrichment (Top Pathways)",
    x = "Pathway",
    y = "-log10(Adjusted P-value)"
  )

cluster_3 <- JoinLayers(cluster_3)


#ident 2

#GO

selected_databases <- c('GO_Biological_Process_2023')

# Run DEenrichRPlot
results_GO_ident2 <- DEenrichRPlot(
  object = cluster_3,               
  ident.1 = 2,           
  ident.2 = NULL,           
  balanced = FALSE,                 
  logfc.threshold = 0.25,          
  assay = "RNA",                    
  max.genes = 50,                  
  test.use = "wilcox",              
  p.val.cutoff = 0.05,         
  cols = NULL,
  enrich.database = selected_databases,  
  num.pathway = 10,             
  return.gene.list = TRUE,
)

# Extract the dataframe from `results`
pos_results_GO_ident2 <- results_GO_ident2$pos  # Assuming "pos" contains the relevant data

# Filter rows where the p-value column (`KEGG_2019_Mouse.P.value`) is less than 0.05
filtered_results_GO_ident2 <- pos_results_GO_ident2[pos_results_GO_ident2$GO_Biological_Process_2023.Adjusted.P.value < 0.05, ]

# View the filtered results
print(filtered_results_GO_ident2)

#DOT PLOT

# Convert adjusted p-values to -log10 for visualization
filtered_results_GO_ident2$log10_pvalue <- -log10(filtered_results_GO_ident2$GO_Biological_Process_2023.Adjusted.P.value)

# Create the dot plot
library(ggplot2)

ggplot(filtered_results_GO_ident2, aes(x = GO_Biological_Process_2023.Odds.Ratio, 
                                       y = reorder(GO_Biological_Process_2023.Term, GO_Biological_Process_2023.Odds.Ratio))) +
  geom_point(aes(size = GO_Biological_Process_2023.Odds.Ratio, color = log10_pvalue)) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(P-value)") +
  theme_minimal() +
  labs(
    title = "GO Biological Process Enrichment Analysis",
    x = "Odds Ratio",
    y = "Pathway",
    size = "Odds Ratio"
  )


#BAR PLOT

# Create the bar plot
ggplot(filtered_results_GO_ident2, aes(x = reorder(GO_Biological_Process_2023.Term, log10_pvalue), y = log10_pvalue)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Flip coordinates to make bars horizontal
  theme_minimal() +
  labs(
    title = "GO Biological Process Enrichment (Top Pathways)",
    x = "Pathway",
    y = "-log10(Adjusted P-value)"
  )



#KEGG
selected_databases <- c('KEGG_2019_Mouse')

# Run DEenrichRPlot
results_kegg_ident2 <- DEenrichRPlot(
  object = cluster_3,               
  ident.1 = 2,           
  ident.2 = NULL,           
  balanced = FALSE,                  
  logfc.threshold = 0.25,           
  assay = "RNA",                   
  max.genes = 50,                   
  test.use = "wilcox",              
  p.val.cutoff = 0.05,         
  cols = NULL,
  enrich.database = selected_databases,  
  num.pathway = 10,            
  return.gene.list = FALSE,
)


# Extract the dataframe from `results`
data_results_kegg_ident2 <- results_kegg_ident2$data  # Assuming "pos" contains the relevant data

# Filter rows where the p-value column (`KEGG_2019_Mouse.P.value`) is less than 0.05
filtered_results_kegg_ident2 <- data_results_kegg_ident2[data_results_kegg_ident2$KEGG_2019_Mouse.Adjusted.P.value < 0.05, ]

# View the filtered results
print(filtered_results_kegg_ident2)

#DOT PLOT

# Convert adjusted p-values to -log10 for color mapping
filtered_results_kegg_ident2$log10_adjusted_pvalue <- -log10(filtered_results_kegg_ident2$KEGG_2019_Mouse.Adjusted.P.value)

# Create the dot plot for KEGG enrichment analysis
library(ggplot2)

ggplot(filtered_results_kegg_ident2, aes(
  x = KEGG_2019_Mouse.Odds.Ratio, 
  y = reorder(KEGG_2019_Mouse.Term, KEGG_2019_Mouse.Odds.Ratio)
)) +
  geom_point(aes(size = KEGG_2019_Mouse.Odds.Ratio, 
                 color = log10_adjusted_pvalue)) +
  scale_color_gradient(low = "blue", high = "red") +  # More significant pathways in red
  theme_minimal() +
  labs(
    title = "KEGG Pathway Enrichment Analysis",
    x = "Odds Ratio",
    y = "Pathway",
    color = "-log10(Adjusted P-value)",
    size = "Odds Ratio"
  )


#BAR PLOT

# Convert adjusted p-values to -log10 for color mapping
filtered_results_kegg_ident2$log10_adjusted_pvalue <- -log10(filtered_results_kegg_ident2$KEGG_2019_Mouse.Adjusted.P.value)

# Create the bar plot for KEGG enrichment analysis
ggplot(filtered_results_kegg_ident2, aes(
  x = reorder(KEGG_2019_Mouse.Term, log10_adjusted_pvalue), 
  y = log10_adjusted_pvalue
)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Flip coordinates to make bars horizontal
  theme_minimal() +
  labs(
    title = "KEGG Pathway Enrichment (Top Pathways)",
    x = "Pathway",
    y = "-log10(Adjusted P-value)"
  )


# Extract cluster and condition information
cluster_condition_df <- data.frame(
  Cluster = Idents(cluster_3), # Clusters
  Condition = cluster_3$condition # Experimental condition
)

# Generate bar plot
library(ggplot2)
ggplot(cluster_condition_df, aes(x = Cluster, fill = Condition)) +
  geom_bar(position = "fill") +
  labs(title = "Cluster Composition by Condition", y = "Proportion", x = "Cluster") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal()

# Generate pie chart (optional)
cluster_condition_summary <- cluster_condition_df %>%
  group_by(Cluster, Condition) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100)

ggplot(cluster_condition_summary, aes(x = "", y = Percentage, fill = Condition)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  facet_wrap(~ Cluster) +
  theme_void() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Cluster Composition by Condition (Pie Chart)")

library(ggplot2)
library(dplyr)

# Summarize data for pie chart
cluster_condition_summary <- cluster_condition_df %>%
  group_by(Cluster, Condition) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Calculate total number of cells for each cluster
cluster_totals <- cluster_condition_summary %>%
  group_by(Cluster) %>%
  summarise(TotalCells = sum(Count))

# Generate pie chart with cluster labels
ggplot(cluster_condition_summary, aes(x = "", y = Percentage, fill = Condition)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  facet_wrap(~ Cluster) +
  geom_text(data = cluster_totals, aes(x = 1.5, y = 0, label = paste0("Total: ", TotalCells)), 
            inherit.aes = FALSE, color = "black", size = 4, fontface = "bold") +
  theme_void() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Cluster Composition by Condition (Pie Chart)")


