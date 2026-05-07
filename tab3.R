# tab3_de.R
# Standalone app for Tab 3: Differential Expression
# Builds directly on Assignment 7 volcano plot code
# Test this independently before combining into the full app

library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(colourpicker)

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("Differential Expression Explorer"),
  sidebarLayout(
    sidebarPanel(
      # File upload for DE results CSV
      fileInput("de_file", "Upload DE Results (CSV or TSV)",
                accept = c(".csv", ".txt", ".tsv")),
      hr(),
      # X and Y axis selectors for volcano plot
      radioButtons("x_axis", "X-axis",
                   choices  = c("log2FoldChange", "baseMean", "stat"),
                   selected = "log2FoldChange"),
      radioButtons("y_axis", "Y-axis",
                   choices  = c("padj", "pvalue"),
                   selected = "padj"),
      # Color pickers for significant vs non-significant genes
      colourInput("color_sig", "Significant gene color",
                  value = "#E64B35"),
      colourInput("color_ns", "Non-significant gene color",
                  value = "#CDC4B5"),
      # Slider for p-adjusted significance threshold (log10 scale)
      sliderInput("padj_slider",
                  "P-adjusted magnitude (log10)",
                  min = -300, max = -1, value = -100),
      submitButton("Plot")
    ),
    mainPanel(
      tabsetPanel(
        # Tab 3a: Sortable DE results table with gene search
        tabPanel("Table", DTOutput("de_table")),
        # Tab 3b: Volcano plot (from Assignment 7, adapted for new data)
        tabPanel("Volcano Plot", plotOutput("volcano_plot"))
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  #' load_de
  #' Loads the DE results from uploaded file or default path.
  #' The DE file has columns: symbol, baseMean, HD.mean, Control.mean,
  #' log2FoldChange, lfcSE, stat, pvalue, padj
  #' Returns a data frame.
  load_de <- reactive({
    if (is.null(input$de_file)) {
      # Default file for testing
      de <- read.table(
        "data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt",
        header = TRUE, sep = "\t", row.names = 1
      )
    } else {
      ext <- tools::file_ext(input$de_file$name)
      sep <- if (ext == "csv") "," else "\t"
      de <- read.table(
        input$de_file$datapath,
        header = TRUE, sep = sep, row.names = 1
      )
    }
    # Remove rows where padj is NA (genes that were filtered by DESeq2)
    de <- de[!is.na(de$padj), ]
    return(de)
  })
  
  #' volcano_plot
  #' Creates a volcano plot from DE results.
  #' Genes below the padj threshold are colored with color_sig,
  #' others with color_ns.
  #' @param dataf DE data frame
  #' @param x_name Column name for x-axis
  #' @param y_name Column name for y-axis (will be -log10 transformed)
  #' @param slider Negative integer for padj threshold (10^slider)
  #' @param color1 Color for significant genes
  #' @param color2 Color for non-significant genes
  #' @return ggplot object
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    # Assign color based on significance threshold
    dataf$color <- ifelse(dataf[[y_name]] < 10^slider, color1, color2)
    
    ggplot(dataf, aes(x = .data[[x_name]],
                      y = -log10(.data[[y_name]]))) +
      geom_point(color = dataf$color, size = 0.8, alpha = 0.6) +
      # Add a horizontal line at the significance threshold
      geom_hline(yintercept = -slider,
                 linetype = "dashed", color = "black", linewidth = 0.5) +
      labs(
        title = paste("Volcano Plot: HD vs Neurologically Normal"),
        x     = x_name,
        y     = paste0("-log10(", y_name, ")"),
        caption = paste0("Dashed line = padj threshold (10^", slider, ")")
      ) +
      theme_bw()
  }
  
  #' draw_table
  #' Filters DE results to genes passing the padj threshold and
  #' formats pvalue/padj columns for display.
  #' @param dataf DE data frame
  #' @param slider Negative integer for padj threshold
  #' @return Filtered and formatted data frame
  draw_table <- function(dataf, slider) {
    # Filter to significant genes only
    filtered <- dataf[dataf$padj < 10^slider, ]
    # Format small p-values in scientific notation for readability
    filtered$pvalue <- formatC(filtered$pvalue, format = "e", digits = 4)
    filtered$padj   <- formatC(filtered$padj,   format = "e", digits = 4)
    # Round other numeric columns
    filtered$baseMean       <- round(filtered$baseMean, 2)
    filtered$log2FoldChange <- round(filtered$log2FoldChange, 3)
    filtered$lfcSE          <- round(filtered$lfcSE, 3)
    filtered$stat           <- round(filtered$stat, 3)
    return(filtered)
  }
  
  #' Render the sortable DE results table
  output$de_table <- renderDT({
    de <- load_de()
    req(de)
    tbl <- draw_table(de, input$padj_slider)
    datatable(
      tbl,
      options  = list(pageLength = 15, scrollX = TRUE),
      rownames = TRUE,
      # Enable gene name search in the table
      filter   = "top"
    )
  })
  
  #' Render the volcano plot
  output$volcano_plot <- renderPlot({
    de <- load_de()
    req(de)
    volcano_plot(
      de,
      input$x_axis,
      input$y_axis,
      input$padj_slider,
      input$color_sig,
      input$color_ns
    )
  })
}

# Run the app
shinyApp(ui = ui, server = server)