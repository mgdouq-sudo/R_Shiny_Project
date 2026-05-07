# tab2_counts.R
# Standalone app for Tab 2: Counts Matrix Exploration
# Test this independently before combining into the full app

library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(pheatmap)  # for clustered heatmap - install with install.packages("pheatmap")

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("Counts Matrix Explorer"),
  sidebarLayout(
    sidebarPanel(
      # File upload for normalized counts matrix
      fileInput("counts_file", "Upload Normalized Counts Matrix (CSV or TSV)",
                accept = c(".csv", ".txt", ".tsv")),
      hr(),
      # Slider: filter genes by variance percentile
      # e.g. 50 means keep only genes with variance >= 50th percentile
      sliderInput("var_percentile",
                  "Minimum percentile of variance",
                  min = 0, max = 100, value = 50, step = 5),
      # Slider: filter genes by minimum number of non-zero samples
      sliderInput("nonzero_samples",
                  "Minimum number of non-zero samples",
                  min = 0, max = 69, value = 10, step = 1),
      # Checkbox to log-transform counts for heatmap visualization
      checkboxInput("log_transform", "Log-transform counts for heatmap", value = TRUE),
      submitButton("Apply Filters")
    ),
    mainPanel(
      tabsetPanel(
        # Tab 2a: Filter summary text/table
        tabPanel("Filter Summary", tableOutput("filter_summary")),
        # Tab 2b: Diagnostic scatter plots
        tabPanel("Diagnostic Plots",
                 plotOutput("var_plot"),    # median vs variance
                 plotOutput("zero_plot")),  # median vs zeros
        # Tab 2c: Clustered heatmap of filtered genes
        tabPanel("Heatmap", plotOutput("heatmap", height = "600px")),
        # Tab 2d: PCA scatter plot
        tabPanel("PCA", plotOutput("pca_plot"))
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  #' load_counts
  #' Loads the normalized counts matrix from an uploaded file or default path.
  #' The counts matrix has genes as rows and samples as columns.
  #' Returns a numeric matrix.
  load_counts <- reactive({
    if (is.null(input$counts_file)) {
      # Default file for testing - tab separated
      counts <- read.table(
        "data/GSE64810_mlhd_DESeq2_norm_counts_adjust.txt",
        header = TRUE, sep = "\t", row.names = 1,
        check.names = FALSE
      )
    } else {
      # Detect separator from file extension
      ext <- tools::file_ext(input$counts_file$name)
      sep <- if (ext == "csv") "," else "\t"
      counts <- read.table(
        input$counts_file$datapath,
        header = TRUE, sep = sep, row.names = 1,
        check.names = FALSE
      )
    }
    # Convert to numeric matrix
    as.matrix(counts)
  })
  
  #' filter_counts
  #' Filters genes based on variance percentile and non-zero sample thresholds.
  #' Returns a list with:
  #'   $filtered  - matrix of genes passing both filters
  #'   $pass      - logical vector indicating which genes pass
  #'   $gene_var  - variance of each gene
  #'   $gene_med  - median count of each gene
  #'   $gene_zeros - number of zero counts per gene
  filter_counts <- reactive({
    counts <- load_counts()
    req(counts)
    
    # Calculate per-gene statistics
    gene_var   <- apply(counts, 1, var)
    gene_med   <- apply(counts, 1, median)
    gene_zeros <- apply(counts, 1, function(x) sum(x == 0))
    
    # Calculate variance threshold from percentile slider
    var_threshold <- quantile(gene_var, probs = input$var_percentile / 100)
    
    # Apply both filters: variance AND non-zero samples
    # non-zero samples = total samples - number of zeros
    nonzero_count <- ncol(counts) - gene_zeros
    pass <- gene_var >= var_threshold & nonzero_count >= input$nonzero_samples
    
    list(
      filtered   = counts[pass, ],
      pass       = pass,
      gene_var   = gene_var,
      gene_med   = gene_med,
      gene_zeros = gene_zeros
    )
  })
  
  #' Render filter summary table showing how many genes pass/fail
  output$filter_summary <- renderTable({
    counts  <- load_counts()
    filtered <- filter_counts()
    req(counts, filtered)
    
    n_total  <- nrow(counts)
    n_pass   <- sum(filtered$pass)
    n_fail   <- n_total - n_pass
    
    data.frame(
      Metric  = c(
        "Number of samples",
        "Total number of genes",
        "Genes passing filter",
        "% genes passing filter",
        "Genes not passing filter",
        "% genes not passing filter"
      ),
      Value = c(
        ncol(counts),
        n_total,
        n_pass,
        paste0(round(100 * n_pass / n_total, 1), "%"),
        n_fail,
        paste0(round(100 * n_fail / n_total, 1), "%")
      ),
      stringsAsFactors = FALSE
    )
  })
  
  #' Render scatter plot: median count vs variance
  #' Genes passing filter are colored darker, filtered-out genes are lighter
  output$var_plot <- renderPlot({
    counts   <- load_counts()
    filtered <- filter_counts()
    req(counts, filtered)
    
    plot_df <- data.frame(
      median   = filtered$gene_med,
      variance = filtered$gene_var,
      pass     = filtered$pass
    )
    
    ggplot(plot_df, aes(x = log10(median + 1),
                        y = log10(variance + 1),
                        color = pass)) +
      geom_point(size = 0.5, alpha = 0.5) +
      scale_color_manual(
        values = c("TRUE" = "#2166ac", "FALSE" = "#d9d9d9"),
        labels = c("TRUE" = "Passes filter", "FALSE" = "Filtered out")
      ) +
      labs(
        title  = "Median Count vs Variance",
        x      = "log10(Median Count + 1)",
        y      = "log10(Variance + 1)",
        color  = "Filter Status"
      ) +
      theme_bw()
  })
  
  #' Render scatter plot: median count vs number of zeros
  output$zero_plot <- renderPlot({
    counts   <- load_counts()
    filtered <- filter_counts()
    req(counts, filtered)
    
    plot_df <- data.frame(
      median = filtered$gene_med,
      zeros  = filtered$gene_zeros,
      pass   = filtered$pass
    )
    
    ggplot(plot_df, aes(x = log10(median + 1),
                        y = zeros,
                        color = pass)) +
      geom_point(size = 0.5, alpha = 0.5) +
      scale_color_manual(
        values = c("TRUE" = "#2166ac", "FALSE" = "#d9d9d9"),
        labels = c("TRUE" = "Passes filter", "FALSE" = "Filtered out")
      ) +
      labs(
        title = "Median Count vs Number of Zero Samples",
        x     = "log10(Median Count + 1)",
        y     = "Number of Zero Samples",
        color = "Filter Status"
      ) +
      theme_bw()
  })
  
  #' Render clustered heatmap of filtered counts
  #' Uses pheatmap for hierarchical clustering of both genes and samples
  output$heatmap <- renderPlot({
    filtered <- filter_counts()
    req(filtered)
    
    mat <- filtered$filtered
    
    # Heatmap needs at least 2 genes to cluster
    if (nrow(mat) < 2) {
      plot.new()
      text(0.5, 0.5, "Too few genes pass filter for heatmap.\nTry loosening the filters.",
           cex = 1.2, col = "red")
      return()
    }
    
    # Cap at 500 genes for performance - take top 500 by variance
    if (nrow(mat) > 500) {
      top_var <- order(apply(mat, 1, var), decreasing = TRUE)[1:500]
      mat <- mat[top_var, ]
    }
    
    # Optionally log-transform for visualization
    if (input$log_transform) {
      mat <- log10(mat + 1)
    }
    
    # Scale each gene (row) to show relative expression patterns
    pheatmap(
      mat,
      scale          = "row",       # z-score normalize each gene
      show_rownames  = FALSE,        # too many genes to label
      show_colnames  = TRUE,
      cluster_rows   = TRUE,
      cluster_cols   = TRUE,
      fontsize_col   = 7,
      color          = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
      main           = paste0("Clustered Heatmap (", nrow(mat), " genes)")
    )
  })
  
  #' Render PCA scatter plot
  #' Projects samples into principal component space and plots PC1 vs PC2
  output$pca_plot <- renderPlot({
    filtered <- filter_counts()
    req(filtered)
    
    mat <- filtered$filtered
    
    if (nrow(mat) < 3) {
      plot.new()
      text(0.5, 0.5, "Too few genes pass filter for PCA.\nTry loosening the filters.",
           cex = 1.2, col = "red")
      return()
    }
    
    # Log transform before PCA to stabilize variance
    mat_log <- log10(mat + 1)
    
    # Run PCA - prcomp expects samples as rows so we transpose
    pca <- prcomp(t(mat_log), scale. = TRUE)
    
    # Calculate % variance explained by each PC
    pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    
    # Build data frame for plotting - use first two PCs
    pca_df <- data.frame(
      PC1     = pca$x[, 1],
      PC2     = pca$x[, 2],
      sample  = rownames(pca$x),
      # Color by condition based on sample name prefix
      condition = ifelse(startsWith(rownames(pca$x), "H_"),
                         "Huntington's Disease", "Neurologically normal")
    )
    
    ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
      geom_point(size = 3, alpha = 0.8) +
      labs(
        title  = "PCA of Normalized Counts",
        x      = paste0("PC1 (", pct_var[1], "% variance)"),
        y      = paste0("PC2 (", pct_var[2], "% variance)"),
        color  = "Condition"
      ) +
      scale_color_manual(values = c(
        "Huntington's Disease"  = "#d6604d",
        "Neurologically normal" = "#2166ac"
      )) +
      theme_bw()
  })
}

# Run the app
shinyApp(ui = ui, server = server)