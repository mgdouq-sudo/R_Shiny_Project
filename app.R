# app.R
# Final Project: Huntington's Disease Gene Expression Analysis
# GSE64810: mRNA-Seq of human post-mortem BA9 brain tissue
# Comparison: 20 HD vs 49 Neurologically Normal controls

library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(colourpicker)
library(pheatmap)
library(ggrepel)

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("Huntington's Disease Gene Expression Explorer"),
  p("Data: GSE64810 — mRNA-Seq of post-mortem human BA9 prefrontal cortex,
     comparing 20 Huntington's Disease and 49 neurologically normal samples."),
  
  # Top-level tabs for each component
  tabsetPanel(
    
    # ── TAB 1: Sample Information ──────────────────────────────────────────
    tabPanel("Samples",
             sidebarLayout(
               sidebarPanel(
                 fileInput("sample_file", "Upload Sample Info CSV",
                           accept = ".csv"),
                 hr(),
                 selectInput("plot_col", "Variable to plot", choices = NULL),
                 selectInput("group_col", "Group by", choices = NULL),
                 actionButton("update_samples", "Update", class = "btn-primary")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", tableOutput("summary_table")),
                   tabPanel("Table",   DTOutput("sample_data_table")),
                   tabPanel("Plots",   plotOutput("violin_plot"))
                 )
               )
             )
    ),
    
    # ── TAB 2: Counts Matrix ───────────────────────────────────────────────
    tabPanel("Counts",
             sidebarLayout(
               sidebarPanel(
                 fileInput("counts_file", "Upload Normalized Counts Matrix (CSV or TSV)",
                           accept = c(".csv", ".txt", ".tsv")),
                 hr(),
                 sliderInput("var_percentile",
                             "Minimum percentile of variance",
                             min = 0, max = 100, value = 50, step = 5),
                 sliderInput("nonzero_samples",
                             "Minimum number of non-zero samples",
                             min = 0, max = 69, value = 10, step = 1),
                 checkboxInput("log_transform",
                               "Log-transform counts for heatmap", value = TRUE),
                 actionButton("update_counts", "Apply Filters", class = "btn-primary")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Filter Summary",  tableOutput("filter_summary")),
                   tabPanel("Diagnostic Plots",
                            plotOutput("var_plot"),
                            plotOutput("zero_plot")),
                   tabPanel("Heatmap", plotOutput("heatmap", height = "600px")),
                   tabPanel("PCA",     plotOutput("pca_plot"))
                 )
               )
             )
    ),
    
    # ── TAB 3: Differential Expression ────────────────────────────────────
    tabPanel("DE",
             sidebarLayout(
               sidebarPanel(
                 fileInput("de_file", "Upload DE Results (CSV or TSV)",
                           accept = c(".csv", ".txt", ".tsv")),
                 hr(),
                 radioButtons("x_axis", "X-axis",
                              choices  = c("log2FoldChange", "baseMean", "stat"),
                              selected = "log2FoldChange"),
                 radioButtons("y_axis", "Y-axis",
                              choices  = c("padj", "pvalue"),
                              selected = "padj"),
                 colourInput("color_sig", "Significant gene color",   value = "#E64B35"),
                 colourInput("color_ns",  "Non-significant gene color", value = "#CDC4B5"),
                 sliderInput("padj_slider",
                             "P-adjusted magnitude (log10)",
                             min = -300, max = -1, value = -5),
                 actionButton("update_de", "Plot", class = "btn-primary")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Table",        DTOutput("de_table")),
                   tabPanel("Volcano Plot", plotOutput("volcano_plot"))
                 )
               )
             )
    ),
    
    # ── TAB 4: GSEA ───────────────────────────────────────────────────────
    tabPanel("GSEA",
             sidebarLayout(
               sidebarPanel(
                 fileInput("gsea_file", "Upload fgsea Results CSV",
                           accept = ".csv"),
                 hr(),
                 h4("Top Results Controls"),
                 sliderInput("top_n",
                             "Number of top pathways to plot",
                             min = 5, max = 50, value = 20, step = 5),
                 radioButtons("nes_filter", "NES direction",
                              choices  = c("All", "Positive", "Negative"),
                              selected = "All"),
                 actionButton("update_barplot", "Update Barplot",
                              class = "btn-primary"),
                 hr(),
                 h4("Table Controls"),
                 sliderInput("padj_filter",
                             "Filter by adjusted p-value (max padj)",
                             min = 0, max = 1, value = 0.05, step = 0.01),
                 downloadButton("download_table", "Download Filtered Table"),
                 hr(),
                 h4("Scatter Plot Controls"),
                 sliderInput("scatter_padj",
                             "Adjusted p-value threshold for coloring",
                             min = 0, max = 1, value = 0.05, step = 0.01),
                 actionButton("update_btn", "Update", class = "btn-primary")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Top Results", plotOutput("barplot",      height = "600px")),
                   tabPanel("Table",       DTOutput("gsea_table")),
                   tabPanel("Plots",       plotOutput("scatter_plot"))
                 )
               )
             )
    )
  )
)

# ── SERVER ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # ── TAB 1: Sample Information ────────────────────────────────────────────
  
  #' load_sample_data: loads sample info CSV from upload or default path
  load_sample_data <- reactive({
    if (is.null(input$sample_file)) {
      return(read.csv("data/sample_info.csv", stringsAsFactors = FALSE))
    }
    ext <- tools::file_ext(input$sample_file$name)
    if (ext != "csv") {
      showNotification("Please upload a CSV file.", type = "error")
      return(NULL)
    }
    read.csv(input$sample_file$datapath, stringsAsFactors = FALSE)
  })
  
  # Update plot column selectors dynamically based on uploaded file
  observe({
    df <- load_sample_data()
    req(df)
    num_cols <- names(df)[sapply(df, is.numeric)]
    all_cols <- names(df)
    updateSelectInput(session, "plot_col", choices = num_cols,
                      selected = num_cols[1])
    updateSelectInput(session, "group_col", choices = all_cols,
                      selected = "diagnosis")
  })
  
  #' make_summary_table: builds column type/value summary data frame
  make_summary_table <- function(df) {
    rows <- lapply(names(df), function(col) {
      vals <- df[[col]]
      if (is.numeric(vals)) {
        m <- round(mean(vals, na.rm = TRUE), 2)
        s <- round(sd(vals,   na.rm = TRUE), 2)
        data.frame(`Column Name` = col, `Type` = "numeric",
                   `Mean (sd) or Distinct Values` = paste0(m, " (+/- ", s, ")"),
                   stringsAsFactors = FALSE, check.names = FALSE)
      } else {
        data.frame(`Column Name` = col, `Type` = "character",
                   `Mean (sd) or Distinct Values` = paste(unique(vals), collapse = ", "),
                   stringsAsFactors = FALSE, check.names = FALSE)
      }
    })
    do.call(rbind, rows)
  }
  
  output$summary_table <- renderTable({
    df <- load_sample_data(); req(df)
    make_summary_table(df)
  })
  
  output$sample_data_table <- renderDT({
    df <- load_sample_data(); req(df)
    datatable(df, options = list(pageLength = 15, scrollX = TRUE),
              rownames = FALSE)
  })
  
  #' make_violin_plot: violin + jitter plot of selected numeric column
  make_violin_plot <- function(df, col_name, group_name) {
    df <- df[!is.na(df[[col_name]]), ]
    ggplot(df, aes(x = .data[[group_name]], y = .data[[col_name]],
                   fill = .data[[group_name]])) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_jitter(width = 0.1, size = 1.5, alpha = 0.5) +
      labs(title = paste("Distribution of", col_name, "by", group_name),
           x = group_name, y = col_name, fill = group_name) +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 15, hjust = 1))
  }
  
  output$violin_plot <- renderPlot({
    input$update_samples
    df <- load_sample_data(); req(df)
    make_violin_plot(df, isolate(input$plot_col), isolate(input$group_col))
  })
  
  # ── TAB 2: Counts Matrix ─────────────────────────────────────────────────
  
  #' load_counts: loads normalized counts matrix from upload or default path
  load_counts <- reactive({
    if (is.null(input$counts_file)) {
      counts <- read.table(
        "data/GSE64810_mlhd_DESeq2_norm_counts_adjust.txt",
        header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    } else {
      ext <- tools::file_ext(input$counts_file$name)
      sep <- if (ext == "csv") "," else "\t"
      counts <- read.table(input$counts_file$datapath,
                           header = TRUE, sep = sep, row.names = 1,
                           check.names = FALSE)
    }
    as.matrix(counts)
  })
  
  #' filter_counts: filters genes by variance percentile and non-zero samples
  filter_counts <- reactive({
    input$update_counts
    counts <- load_counts(); req(counts)
    gene_var   <- apply(counts, 1, var)
    gene_med   <- apply(counts, 1, median)
    gene_zeros <- apply(counts, 1, function(x) sum(x == 0))
    var_threshold <- quantile(gene_var,
                              probs = isolate(input$var_percentile) / 100)
    nonzero_count <- ncol(counts) - gene_zeros
    pass <- gene_var >= var_threshold &
      nonzero_count >= isolate(input$nonzero_samples)
    list(filtered = counts[pass, ], pass = pass,
         gene_var = gene_var, gene_med = gene_med, gene_zeros = gene_zeros)
  })
  
  output$filter_summary <- renderTable({
    counts   <- load_counts();   req(counts)
    filtered <- filter_counts(); req(filtered)
    n_total <- nrow(counts); n_pass <- sum(filtered$pass)
    n_fail  <- n_total - n_pass
    data.frame(
      Metric = c("Number of samples", "Total number of genes",
                 "Genes passing filter", "% genes passing filter",
                 "Genes not passing filter", "% genes not passing filter"),
      Value  = c(ncol(counts), n_total, n_pass,
                 paste0(round(100 * n_pass / n_total, 1), "%"),
                 n_fail, paste0(round(100 * n_fail / n_total, 1), "%")),
      stringsAsFactors = FALSE)
  })
  
  output$var_plot <- renderPlot({
    counts <- load_counts(); filtered <- filter_counts(); req(counts, filtered)
    plot_df <- data.frame(median = filtered$gene_med,
                          variance = filtered$gene_var, pass = filtered$pass)
    ggplot(plot_df, aes(x = log10(median + 1), y = log10(variance + 1),
                        color = pass)) +
      geom_point(size = 0.5, alpha = 0.5) +
      scale_color_manual(values = c("TRUE" = "#2166ac", "FALSE" = "#d9d9d9"),
                         labels = c("TRUE" = "Passes filter",
                                    "FALSE" = "Filtered out")) +
      labs(title = "Median Count vs Variance",
           x = "log10(Median Count + 1)", y = "log10(Variance + 1)",
           color = "Filter Status") + theme_bw()
  })
  
  output$zero_plot <- renderPlot({
    counts <- load_counts(); filtered <- filter_counts(); req(counts, filtered)
    plot_df <- data.frame(median = filtered$gene_med,
                          zeros = filtered$gene_zeros, pass = filtered$pass)
    ggplot(plot_df, aes(x = log10(median + 1), y = zeros, color = pass)) +
      geom_point(size = 0.5, alpha = 0.5) +
      scale_color_manual(values = c("TRUE" = "#2166ac", "FALSE" = "#d9d9d9"),
                         labels = c("TRUE" = "Passes filter",
                                    "FALSE" = "Filtered out")) +
      labs(title = "Median Count vs Number of Zero Samples",
           x = "log10(Median Count + 1)", y = "Number of Zero Samples",
           color = "Filter Status") + theme_bw()
  })
  
  output$heatmap <- renderPlot({
    filtered <- filter_counts(); req(filtered)
    mat <- filtered$filtered
    if (nrow(mat) < 2) {
      plot.new()
      text(0.5, 0.5, "Too few genes pass filter.\nTry loosening the filters.",
           cex = 1.2, col = "red"); return()
    }
    if (nrow(mat) > 500) {
      top_var <- order(apply(mat, 1, var), decreasing = TRUE)[1:500]
      mat <- mat[top_var, ]
    }
    if (isolate(input$log_transform)) mat <- log10(mat + 1)
    pheatmap(mat, scale = "row", show_rownames = FALSE, show_colnames = TRUE,
             cluster_rows = TRUE, cluster_cols = TRUE, fontsize_col = 7,
             color = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
             main = paste0("Clustered Heatmap (", nrow(mat), " genes)"))
  })
  
  output$pca_plot <- renderPlot({
    filtered <- filter_counts(); req(filtered)
    mat <- filtered$filtered
    if (nrow(mat) < 3) {
      plot.new()
      text(0.5, 0.5, "Too few genes pass filter for PCA.",
           cex = 1.2, col = "red"); return()
    }
    mat_log <- log10(mat + 1)
    pca     <- prcomp(t(mat_log), scale. = TRUE)
    pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    pca_df  <- data.frame(
      PC1 = pca$x[, 1], PC2 = pca$x[, 2],
      sample = rownames(pca$x),
      condition = ifelse(startsWith(rownames(pca$x), "H_"),
                         "Huntington's Disease", "Neurologically normal"))
    ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(values = c("Huntington's Disease"  = "#d6604d",
                                    "Neurologically normal" = "#2166ac")) +
      labs(title = "PCA of Normalized Counts",
           x = paste0("PC1 (", pct_var[1], "% variance)"),
           y = paste0("PC2 (", pct_var[2], "% variance)"),
           color = "Condition") + theme_bw()
  })
  
  # ── TAB 3: Differential Expression ──────────────────────────────────────
  
  #' load_de: loads DE results from upload or default path
  load_de <- reactive({
    if (is.null(input$de_file)) {
      de <- read.table(
        "data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt",
        header = TRUE, sep = "\t", row.names = 1)
    } else {
      ext <- tools::file_ext(input$de_file$name)
      sep <- if (ext == "csv") "," else "\t"
      de  <- read.table(input$de_file$datapath,
                        header = TRUE, sep = sep, row.names = 1)
    }
    de[!is.na(de$padj), ]
  })
  
  #' volcano_plot: ggplot volcano with significance threshold line
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    dataf$color <- ifelse(dataf[[y_name]] < 10^slider, color1, color2)
    ggplot(dataf, aes(x = .data[[x_name]], y = -log10(.data[[y_name]]))) +
      geom_point(color = dataf$color, size = 0.8, alpha = 0.6) +
      geom_hline(yintercept = -slider, linetype = "dashed",
                 color = "black", linewidth = 0.5) +
      labs(title = "Volcano Plot: HD vs Neurologically Normal",
           x = x_name, y = paste0("-log10(", y_name, ")"),
           caption = paste0("Dashed line = padj threshold (10^", slider, ")")) +
      theme_bw()
  }
  
  #' draw_table: filters DE results and formats p-value columns
  draw_table <- function(dataf, slider) {
    filtered <- dataf[dataf$padj < 10^slider, ]
    filtered$pvalue         <- formatC(filtered$pvalue, format = "e", digits = 4)
    filtered$padj           <- formatC(filtered$padj,   format = "e", digits = 4)
    filtered$baseMean       <- round(filtered$baseMean, 2)
    filtered$log2FoldChange <- round(filtered$log2FoldChange, 3)
    filtered$lfcSE          <- round(filtered$lfcSE, 3)
    filtered$stat           <- round(filtered$stat, 3)
    filtered
  }
  
  output$de_table <- renderDT({
    input$update_de
    de <- load_de(); req(de)
    datatable(draw_table(de, isolate(input$padj_slider)),
              options = list(pageLength = 15, scrollX = TRUE),
              rownames = TRUE, filter = "top")
  })
  
  output$volcano_plot <- renderPlot({
    input$update_de
    de <- load_de(); req(de)
    volcano_plot(de, isolate(input$x_axis), isolate(input$y_axis),
                 isolate(input$padj_slider),
                 isolate(input$color_sig), isolate(input$color_ns))
  })
  
  # ── TAB 4: GSEA ─────────────────────────────────────────────────────────
  
  #' load_gsea: loads fgsea results from upload or default path
  load_gsea <- reactive({
    if (is.null(input$gsea_file)) {
      gsea <- read.csv("data/fgsea_results.csv", stringsAsFactors = FALSE)
    } else {
      gsea <- read.csv(input$gsea_file$datapath, stringsAsFactors = FALSE)
    }
    gsea[!is.na(gsea$padj), ]
  })
  
  #' filter_gsea: filters by padj threshold and NES direction
  filter_gsea <- function(gsea, padj_thresh, nes_dir) {
    filtered <- gsea[gsea$padj <= padj_thresh, ]
    if (nes_dir == "Positive") filtered <- filtered[filtered$NES > 0, ]
    if (nes_dir == "Negative") filtered <- filtered[filtered$NES < 0, ]
    filtered
  }
  
  output$barplot <- renderPlot({
    input$update_barplot
    gsea <- load_gsea(); req(gsea)
    filtered <- filter_gsea(gsea, padj_thresh = 1,
                            nes_dir = isolate(input$nes_filter))
    top <- head(filtered[order(filtered$padj), ], isolate(input$top_n))
    top$pathway_label <- gsub("HALLMARK_", "", top$pathway)
    top$pathway_label <- gsub("_", " ", top$pathway_label)
    top <- top[order(top$NES), ]
    top$pathway_label <- factor(top$pathway_label, levels = top$pathway_label)
    ggplot(top, aes(x = NES, y = pathway_label, fill = NES > 0)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(
        values = c("TRUE" = "#d6604d", "FALSE" = "#2166ac"),
        labels = c("TRUE" = "Activated in HD", "FALSE" = "Suppressed in HD")
      ) +
      labs(title = paste("Top", input$top_n, "Pathways by Adjusted P-value"),
           x = "Normalized Enrichment Score (NES)", y = "Pathway",
           fill = "Direction") +
      theme_bw() +
      theme(axis.text.y = element_text(size = 9))
  })
  
  output$gsea_table <- renderDT({
    gsea <- load_gsea(); req(gsea)
    filtered <- filter_gsea(gsea, input$padj_filter, input$nes_filter)
    filtered$NES  <- round(filtered$NES, 3)
    filtered$pval <- formatC(filtered$pval, format = "e", digits = 3)
    filtered$padj <- formatC(filtered$padj, format = "e", digits = 3)
    filtered <- filtered[, !names(filtered) %in% "leadingEdge"]
    datatable(filtered, options = list(pageLength = 15, scrollX = TRUE),
              rownames = FALSE)
  })
  
  output$download_table <- downloadHandler(
    filename = function() {
      paste0("fgsea_filtered_padj", input$padj_filter,
             "_", input$nes_filter, ".csv")
    },
    content = function(file) {
      gsea     <- load_gsea()
      filtered <- filter_gsea(gsea, input$padj_filter, input$nes_filter)
      filtered <- filtered[, !names(filtered) %in% "leadingEdge"]
      write.csv(filtered, file, row.names = FALSE)
    }
  )
  
  output$scatter_plot <- renderPlot({
    input$update_btn
    gsea <- load_gsea(); req(gsea)
    scatter_padj <- isolate(input$scatter_padj)
    gsea$pathway_label <- gsub("HALLMARK_", "", gsea$pathway)
    gsea$pathway_label <- gsub("_", " ", gsea$pathway_label)
    gsea$significant   <- gsea$padj < scatter_padj
    ggplot(gsea, aes(x = NES, y = -log10(padj), color = significant)) +
      geom_point(size = 2, alpha = 0.8) +
      ggrepel::geom_text_repel(
        data = gsea[gsea$significant, ],
        aes(label = pathway_label), size = 3, max.overlaps = 15) +
      scale_color_manual(
        values = c("TRUE" = "#d6604d", "FALSE" = "#AAAAAA"),
        labels = c("TRUE" = paste0("padj < ", scatter_padj),
                   "FALSE" = "Not significant")) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = -log10(scatter_padj),
                 linetype = "dashed", color = "red") +
      labs(title = "GSEA: NES vs Significance",
           x = "Normalized Enrichment Score (NES)",
           y = "-log10(adjusted p-value)", color = "Significance") +
      theme_bw()
  })
}

# Run the application
shinyApp(ui = ui, server = server)