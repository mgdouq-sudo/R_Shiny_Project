# BF530 R Shiny Final Project

**Author:** Mohammad Gharandouq

**Overview:** This Shiny application explores gene expression data from a Huntington's Disease study (GEO, Labadorf et al., 2015), comparing mRNA expression in post-mortem human prefrontal cortex from 20 HD patients and 49 neurologically normal controls. The app uses pre-computed DESeq2 normalized counts and differential expression results provided by the authors. Pathway analysis was performed with fgsea using log2 fold change as the ranking metric against MSigDB Hallmark gene sets, with significance defined at FDR < 0.05.

---

## How to Run

1. Clone or download this repository.
2. Open R or RStudio and set the working directory to the project folder.
3. Install required packages (see Dependencies below).
4. Run the app with:

```r
shiny::runApp()
```

---

## Required Data Files

The following files must be uploaded through the app interface:

- **Samples tab:** A CSV file containing sample metadata (e.g. diagnosis, age of death, PMI, RIN) for all 69 samples.
- **Counts tab:** The DESeq2 normalized counts matrix (~28,000 genes × 69 samples).
- **DE tab:** The DESeq2 differential expression results table.
- **GSEA tab:** The fgsea output results table.

All files are derived from the dataset published by Dr. Labadorf, available on [GEO](https://www.ncbi.nlm.nih.gov/geo/).

---

## Application Tabs

### 1. Sample Information Explorer
Accepts a CSV metadata file. Displays a summary of each column (data type, mean/SD for numeric variables, distinct values for categorical variables), a sortable table of all 69 samples, and violin plots with customizable grouping and variable selection.

### 2. Counts Matrix Explorer
Loads the DESeq2 normalized counts matrix. Allows filtering genes by variance percentile and minimum non-zero samples, with an optional log transformation toggle. Displays a filter summary, diagnostic scatter plots (median count vs. variance and vs. zeros), a clustered heatmap of the top 500 most variable genes, and a PCA plot colored by condition.

### 3. Differential Expression Explorer
Volcano plot visualization of DE results. Supports customizable x/y axes, point colors, and an adjustable p-adjusted significance threshold. Significant genes are highlighted in red. Includes a searchable, sortable table of significant genes.

### 4. GSEA Explorer
Displays fgsea results against MSigDB Hallmark gene sets. Includes a barplot of top pathways by adjusted p-value (filterable by NES direction), a scatter plot of NES vs. −log10 adjusted p-value with labeled significant pathways, and a filterable results table with CSV export.

---

## Dependencies

```r
install.packages(c("shiny", "ggplot2", "DT", "colourpicker", "pheatmap", "dplyr"))
# Bioconductor packages:
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "fgsea"))
```

---

## Notes & Limitations

- The heatmap is capped at the top 500 most variable genes for rendering performance; the full dataset contains over 28,000 genes.
- All `actionButton` inputs are paired with `isolate()` to control when outputs update, avoiding global reactivity conflicts across tabs.
- All input and output IDs are uniquely namespaced to prevent conflicts across the four tab modules.

## References

- Labadorf A, Hoss AG, Lagomarsino V, et al. RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression. PLoS One. 2015;10(12):e0143563. Published 2015 Dec 4. doi:10.1371/journal.pone.0143563
