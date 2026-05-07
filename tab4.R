# tab4_gsea.R
# Standalone app for Tab 4: Gene Set Enrichment Analysis
# Test this independently before combining into the full app

library(shiny)
library(bslib)
library(ggplot2)
library(DT)

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("Gene Set Enrichment Analysis Explorer"),
  sidebarLayout(
    sidebarPanel(
      # File upload for fgsea results CSV
      fileInput("gsea_file", "Upload fgsea Results CSV",
                accept = ".csv"),
      
      hr(),
      
      # Tab 1 controls: slider for top N pathways to show in barplot
      h4("Top Results Controls"),
      sliderInput("top_n",
                  "Number of top pathways to plot",
                  min = 5, max = 50, value = 20, step = 5),
      
      hr(),
      
      # Tab 2 controls: filter by padj and NES direction
      h4("Table Controls"),
      sliderInput("padj_filter",
                  "Filter by adjusted p-value (max padj)",
                  min = 0, max = 1, value = 0.05, step = 0.01),
      radioButtons("nes_filter", "NES direction",
                   choices  = c("All", "Positive", "Negative"),
                   selected = "All"),
      # Download button for filtered table
      downloadButton("download_table", "Download Filtered Table"),
      
      hr(),
      
      # Tab 3 controls: padj threshold for scatter plot coloring
      h4("Scatter Plot Controls"),
      sliderInput("scatter_padj",
                  "Adjusted p-value threshold for coloring",
                  min = 0, max = 1, value = 0.05, step = 0.01),
      
      actionButton("update_btn", "Update", class = "btn-primary")
    ),
    mainPanel(
      tabsetPanel(
        # Tab 4a: Barplot of top N pathways by padj
        tabPanel("Top Results", plotOutput("barplot", height = "600px")),
        # Tab 4b: Filterable sortable table of fgsea results
        tabPanel("Table", DTOutput("gsea_table")),
        # Tab 4c: NES vs -log10(padj) scatter plot
        tabPanel("Plots", plotOutput("scatter_plot"))
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  #' load_gsea
  #' Loads fgsea results from uploaded file or default path.
  #' Returns a data frame with columns: pathway, pval, padj, NES, etc.
  load_gsea <- reactive({
    if (is.null(input$gsea_file)) {
      # Default file for testing
      gsea <- read.csv("data/fgsea_results.csv", stringsAsFactors = FALSE)
    } else {
      gsea <- read.csv(input$gsea_file$datapath, stringsAsFactors = FALSE)
    }
    # Remove rows with NA padj
    gsea <- gsea[!is.na(gsea$padj), ]
    return(gsea)
  })
  
  #' filter_gsea
  #' Filters fgsea results by padj threshold and NES direction.
  #' @param gsea Data frame of fgsea results
  #' @param padj_thresh Maximum padj value to include
  #' @param nes_dir "All", "Positive", or "Negative"
  #' @return Filtered data frame
  filter_gsea <- function(gsea, padj_thresh, nes_dir) {
    filtered <- gsea[gsea$padj <= padj_thresh, ]
    if (nes_dir == "Positive") {
      filtered <- filtered[filtered$NES > 0, ]
    } else if (nes_dir == "Negative") {
      filtered <- filtered[filtered$NES < 0, ]
    }
    return(filtered)
  }
  
  #' Render barplot of top N pathways by adjusted p-value
  #' Bars are colored by NES: positive = red (activated in HD),
  #' negative = blue (suppressed in HD)
  output$barplot <- renderPlot({
    gsea <- load_gsea()
    req(gsea)
    
    # Apply NES direction filter before selecting top N
    filtered <- filter_gsea(gsea, padj_thresh = 1, nes_dir = input$nes_filter)
    
    # Select top N pathways by adjusted p-value
    top <- head(filtered[order(filtered$padj), ], input$top_n)
    
    # Clean up pathway names for display (remove HALLMARK_ prefix)
    top$pathway_label <- gsub("HALLMARK_", "", top$pathway)
    top$pathway_label <- gsub("_", " ", top$pathway_label)
    
    # Order bars by NES for cleaner visualization
    top <- top[order(top$NES), ]
    top$pathway_label <- factor(top$pathway_label,
                                levels = top$pathway_label)
    
    ggplot(top, aes(x = NES, y = pathway_label, fill = NES > 0)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(
        values = c("TRUE" = "#d6604d", "FALSE" = "#2166ac"),
        labels = c("TRUE" = "Activated in HD", "FALSE" = "Suppressed in HD")
      ) +
      labs(
        title = paste("Top", input$top_n, "Pathways by Adjusted P-value"),
        x     = "Normalized Enrichment Score (NES)",
        y     = "Pathway",
        fill  = "Direction"
      ) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 9))
  })
  
  #' Render filterable sortable table of fgsea results
  output$gsea_table <- renderDT({
    gsea <- load_gsea()
    req(gsea)
    
    # Apply filters from sidebar
    filtered <- filter_gsea(gsea, input$padj_filter, input$nes_filter)
    
    # Round numeric columns for display
    filtered$NES    <- round(filtered$NES, 3)
    filtered$pval   <- formatC(filtered$pval, format = "e", digits = 3)
    filtered$padj   <- formatC(filtered$padj, format = "e", digits = 3)
    
    # Drop leadingEdge column - too long to display nicely
    filtered <- filtered[, !names(filtered) %in% "leadingEdge"]
    
    datatable(
      filtered,
      options  = list(pageLength = 15, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  #' Download handler for filtered table
  output$download_table <- downloadHandler(
    filename = function() {
      paste0("fgsea_filtered_padj", input$padj_filter, "_",
             input$nes_filter, ".csv")
    },
    content = function(file) {
      gsea     <- load_gsea()
      filtered <- filter_gsea(gsea, input$padj_filter, input$nes_filter)
      # Remove leadingEdge for cleaner export
      filtered <- filtered[, !names(filtered) %in% "leadingEdge"]
      write.csv(filtered, file, row.names = FALSE)
    }
  )
  
  #' Render NES vs -log10(padj) scatter plot
  #' Gene sets below padj threshold are colored, others are grey
  output$scatter_plot <- renderPlot({
    input$update_btn  # re-run when Update is clicked
    gsea <- load_gsea()
    req(gsea)
    
    # Use isolate so it only updates on button click
    scatter_padj <- isolate(input$scatter_padj)
    
    # Clean pathway names for hover labels
    gsea$pathway_label <- gsub("HALLMARK_", "", gsea$pathway)
    gsea$pathway_label <- gsub("_", " ", gsea$pathway_label)
    
    # Color by significance threshold
    gsea$significant <- gsea$padj < scatter_padj
    
    ggplot(gsea, aes(x = NES,
                     y = -log10(padj),
                     color = significant)) +
      geom_point(size = 2, alpha = 0.8) +
      # Label the significant points with pathway names
      ggrepel::geom_text_repel(
        data = gsea[gsea$significant, ],
        aes(label = pathway_label),
        size = 3, max.overlaps = 15
      ) +
      scale_color_manual(
        values = c("TRUE" = "#d6604d", "FALSE" = "#AAAAAA"),
        labels = c("TRUE" = paste0("padj < ", scatter_padj),
                   "FALSE" = "Not significant")
      ) +
      # Add vertical line at NES = 0
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      # Add horizontal line at significance threshold
      geom_hline(yintercept = -log10(scatter_padj),
                 linetype = "dashed", color = "red") +
      labs(
        title  = "GSEA: NES vs Significance",
        x      = "Normalized Enrichment Score (NES)",
        y      = "-log10(adjusted p-value)",
        color  = "Significance"
      ) +
      theme_bw()
  })
}

# Run the app
shinyApp(ui = ui, server = server)