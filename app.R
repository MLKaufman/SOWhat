library(shiny)
library(Seurat)
library(bslib)
library(ggplot2)
library(magrittr)
library(gtools)

# Set max upload size to 2GB for large Seurat objects
options(shiny.maxRequestSize = 2000 * 1024^2)

ui <- page_sidebar(
  title = "Seurat Explorer High-Performance",
  sidebar = sidebar(
    fileInput("file", "Upload Seurat Object (.rds)", accept = ".rds"),
    selectInput("resolution", "Clustering Resolution", choices = NULL),
    hr(),
    card(
      card_header("Object Statistics"),
      tableOutput("stats_table")
    ),
    hr(),
    helpText("Note: Plot uses rasterization for high performance.")
  ),
  layout_column_wrap(
    width = 1,
    card(
      card_header("UMAP Projection"),
      plotly::plotlyOutput("umap_plot", height = "600px")
    ),
    card(
      card_header(
        "Cluster Annotations",
        class = "d-flex justify-content-between align-items-center",
        downloadButton("download_annotes", "Save Annotations", class = "btn-sm")
      ),
      uiOutput("cluster_table_ui")
    )
  )
)

server <- function(input, output, session) {
  # Reactive value to store the Seurat object
  raw_obj <- reactive({
    path <- if (!is.null(input$file)) {
      input$file$datapath
    } else if (file.exists("sample_seurat.rds")) {
      "sample_seurat.rds"
    } else {
      NULL
    }

    if (is.null(path)) {
      return(NULL)
    }

    withProgress(message = "Loading Seurat object...", {
      readRDS(path)
    })
  })

  # Update resolution choices when object is loaded
  observe({
    obj <- raw_obj()
    if (is.null(obj)) {
      return()
    }

    metadata_cols <- colnames(obj@meta.data)
    # Find columns that look like clustering results (e.g., RNA_snn_res.0.5)
    res_cols <- grep("res\\.", metadata_cols, value = TRUE)

    # If no specific resolution columns found, just show all factor/character columns
    if (length(res_cols) == 0) {
      res_cols <- metadata_cols
    }

    updateSelectInput(session, "resolution", choices = res_cols, selected = head(res_cols, 1))
  })

  # Render the UMAP plot
  output$umap_plot <- plotly::renderPlotly({
    obj <- raw_obj()
    req(obj, input$resolution)

    # Ensure the resolution is set as the active identity
    Idents(obj) <- input$resolution

    p <- DimPlot(obj, reduction = "umap", label = TRUE, raster = FALSE) +
      theme_minimal() +
      labs(title = paste("Resolution:", input$resolution))

    plotly::ggplotly(p) %>%
      plotly::layout(yaxis = list(scaleanchor = "x", scaleratio = 1))
  })

  # Reactive value for cluster annotations
  annotes <- reactiveValues(data = list())

  # UI for the editable table
  output$cluster_table_ui <- renderUI({
    obj <- raw_obj()
    req(obj, input$resolution)

    # Get unique clusters for the selected resolution
    clusters <- gtools::mixedsort(unique(as.character(obj@meta.data[[input$resolution]])))

    # Create HTML table
    tags$table(
      class = "table table-hover table-sm",
      tags$thead(
        tags$tr(
          tags$th("Cluster"),
          tags$th("Annotation"),
          tags$th("Reason")
        )
      ),
      tags$tbody(
        lapply(clusters, function(cl) {
          # Use namespaced IDs based on resolution and cluster
          id_ann <- paste0("ann_", input$resolution, "_", cl)
          id_reason <- paste0("reason_", input$resolution, "_", cl)

          # Retrieve existing values if any
          val_ann <- if (!is.null(annotes$data[[id_ann]])) annotes$data[[id_ann]] else ""
          val_reason <- if (!is.null(annotes$data[[id_reason]])) annotes$data[[id_reason]] else ""

          tags$tr(
            tags$td(cl),
            tags$td(textInput(id_ann, NULL, value = val_ann, width = "100%", placeholder = "Enter annotation...")),
            tags$td(textInput(id_reason, NULL, value = val_reason, width = "100%", placeholder = "Enter reason..."))
          )
        })
      )
    )
  })

  # Observer to save inputs into the reactiveValues
  observe({
    obj <- raw_obj()
    req(obj, input$resolution)
    clusters <- gtools::mixedsort(unique(as.character(obj@meta.data[[input$resolution]])))

    for (cl in clusters) {
      id_ann <- paste0("ann_", input$resolution, "_", cl)
      id_reason <- paste0("reason_", input$resolution, "_", cl)

      # Use isolate to avoid recursive invalidation
      isolate({
        if (!is.null(input[[id_ann]])) annotes$data[[id_ann]] <- input[[id_ann]]
        if (!is.null(input[[id_reason]])) annotes$data[[id_reason]] <- input[[id_reason]]
      })
    }
  })

  # Render statistics table
  output$stats_table <- renderTable({
    obj <- raw_obj()
    req(obj)

    data.frame(
      Metric = c("Number of Cells", "Number of Features", "Assays Available", "Active Resolution"),
      Value = c(
        ncol(obj),
        nrow(obj),
        paste(Assays(obj), collapse = ", "),
        input$resolution
      )
    )
  })

  # Download handler for annotations
  output$download_annotes <- downloadHandler(
    filename = function() {
      paste0("cluster_annotations_", input$resolution, ".csv")
    },
    content = function(file) {
      obj <- raw_obj()
      req(obj, input$resolution)
      clusters <- gtools::mixedsort(unique(as.character(obj@meta.data[[input$resolution]])))

      # Prepare data for export
      export_data <- data.frame(
        Cluster = clusters,
        Annotation = sapply(clusters, function(cl) {
          id <- paste0("ann_", input$resolution, "_", cl)
          if (!is.null(annotes$data[[id]])) annotes$data[[id]] else ""
        }),
        Reason = sapply(clusters, function(cl) {
          id <- paste0("reason_", input$resolution, "_", cl)
          if (!is.null(annotes$data[[id]])) annotes$data[[id]] else ""
        }),
        stringsAsFactors = FALSE
      )

      write.csv(export_data, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
