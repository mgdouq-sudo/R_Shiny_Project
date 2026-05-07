# tab1_samples.R
# Standalone app for Tab 1: Sample Information Exploration
# Test this independently before combining into the full app

library(shiny)
library(bslib)
library(ggplot2)
library(DT)

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("Sample Information Explorer"),
  sidebarLayout(
    sidebarPanel(
      # File upload input - accepts CSV only
      fileInput("sample_file", "Upload Sample Info CSV",
                accept = ".csv"),
      hr(),
      # Only show these controls when on the Plots tab
      # Column selector for x-axis (which variable to plot)
      selectInput("plot_col", "Variable to plot",
                  choices = NULL),  # populated dynamically from uploaded file
      # Column selector for grouping (color by)
      selectInput("group_col", "Group by",
                  choices = NULL),  # populated dynamically from uploaded file
      submitButton("Update")
    ),
    mainPanel(
      tabsetPanel(
        # Tab 1a: Summary table of column types and values
        tabPanel("Summary", tableOutput("summary_table")),
        # Tab 1b: Raw data table with sortable columns
        tabPanel("Table", DTOutput("data_table")),
        # Tab 1c: Violin plots of continuous variables
        tabPanel("Plots", plotOutput("violin_plot"))
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  #' load_sample_data
  #' Reactive function that loads the uploaded CSV file.
  #' Falls back to the default sample_info.csv if no file is uploaded.
  #' Returns a data frame.
  load_sample_data <- reactive({
    if (is.null(input$sample_file)) {
      # Default file for testing
      return(read.csv("data/sample_info.csv", stringsAsFactors = FALSE))
    }
    # Validate that uploaded file is a CSV
    ext <- tools::file_ext(input$sample_file$name)
    if (ext != "csv") {
      showNotification("Please upload a CSV file.", type = "error")
      return(NULL)
    }
    read.csv(input$sample_file$datapath, stringsAsFactors = FALSE)
  })
  
  #' Update selectInput choices dynamically based on uploaded file columns.
  #' Numeric columns go to plot_col, all columns go to group_col.
  observe({
    df <- load_sample_data()
    req(df)
    
    # Find numeric columns for plotting
    num_cols <- names(df)[sapply(df, is.numeric)]
    # All columns for grouping
    all_cols <- names(df)
    
    updateSelectInput(session, "plot_col", choices = num_cols,
                      selected = num_cols[1])
    updateSelectInput(session, "group_col", choices = all_cols,
                      selected = "diagnosis")
  })
  
  #' make_summary_table
  #' Builds a summary data frame showing each column's type and
  #' either mean ± sd (for numeric) or distinct values (for character/factor).
  #' @param df Data frame loaded by load_sample_data()
  #' @return A data frame with columns: Column Name, Type, Mean (sd) or Distinct Values
  make_summary_table <- function(df) {
    rows <- lapply(names(df), function(col) {
      vals <- df[[col]]
      if (is.numeric(vals)) {
        # Calculate mean and sd, ignoring NAs
        m  <- round(mean(vals, na.rm = TRUE), 2)
        s  <- round(sd(vals,   na.rm = TRUE), 2)
        data.frame(
          `Column Name` = col,
          `Type`        = "numeric",
          `Mean (sd) or Distinct Values` = paste0(m, " (+/- ", s, ")"),
          stringsAsFactors = FALSE, check.names = FALSE
        )
      } else {
        # Show unique values for categorical columns
        distinct <- paste(unique(vals), collapse = ", ")
        data.frame(
          `Column Name` = col,
          `Type`        = "character",
          `Mean (sd) or Distinct Values` = distinct,
          stringsAsFactors = FALSE, check.names = FALSE
        )
      }
    })
    do.call(rbind, rows)
  }
  
  #' Render the summary table
  output$summary_table <- renderTable({
    df <- load_sample_data()
    req(df)
    
    # Add header info above the table
    summary_df <- make_summary_table(df)
    summary_df
  })
  
  # Add number of rows/columns as text above the table
  output$summary_table <- renderTable({
    df <- load_sample_data()
    req(df)
    make_summary_table(df)
  })
  
  #' Render the raw data table with sorting enabled via DT package
  output$data_table <- renderDT({
    df <- load_sample_data()
    req(df)
    datatable(df,
              options = list(pageLength = 15, scrollX = TRUE),
              rownames = FALSE)
  })
  
  #' make_violin_plot
  #' Creates a violin plot of a selected numeric column, grouped by a
  #' selected categorical column, colored by group.
  #' @param df Data frame
  #' @param col_name Name of the numeric column to plot on y-axis
  #' @param group_name Name of the column to group/color by
  #' @return A ggplot object
  make_violin_plot <- function(df, col_name, group_name) {
    # Remove rows where the selected column is NA
    df <- df[!is.na(df[[col_name]]), ]
    
    ggplot(df, aes(x = .data[[group_name]],
                   y = .data[[col_name]],
                   fill = .data[[group_name]])) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      # Add individual data points on top of violin
      geom_jitter(width = 0.1, size = 1.5, alpha = 0.5) +
      labs(
        title = paste("Distribution of", col_name, "by", group_name),
        x     = group_name,
        y     = col_name,
        fill  = group_name
      ) +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 15, hjust = 1))
  }
  
  #' Render the violin plot
  output$violin_plot <- renderPlot({
    df <- load_sample_data()
    req(df, input$plot_col, input$group_col)
    make_violin_plot(df, input$plot_col, input$group_col)
  })
}

# Run the app
shinyApp(ui = ui, server = server)