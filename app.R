library(shiny)
library(Seurat)
library(bslib)
library(ggplot2)
library(magrittr)
library(gtools)
library(DT)
library(dplyr)

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
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      width = "65%",
      navset_card_tab(
        title = "Projection & Analysis",
        nav_panel(
          "UMAP",
          plotly::plotlyOutput("umap_plot", height = "800px")
        ),
        nav_panel(
          "Dot Plot",
          card_body(
            layout_column_wrap(
              width = 1,
              textAreaInput("genes", "Enter Gene List",
                placeholder = "e.g., PTPRC, CD3E, CD19", rows = 3, width = "100%"
              ),
              div(
                class = "d-flex gap-2 mb-3",
                actionButton("preset_immune", "Immune", class = "btn-outline-primary btn-sm"),
                actionButton("preset_death", "Cell Death", class = "btn-outline-primary btn-sm"),
                actionButton("preset_prolif", "Proliferation", class = "btn-outline-primary btn-sm"),
                actionButton("preset_cycle", "Cell Cycle", class = "btn-outline-primary btn-sm")
              )
            ),
            plotOutput("dot_plot", height = "800px")
          )
        ),
        nav_panel(
          "Markers",
          layout_sidebar(
            sidebar = sidebar(
              actionButton("run_markers", "Find All Markers", class = "btn-primary w-100"),
              hr(),
              helpText("Calculates markers for all clusters at the current resolution.")
            ),
            card(
              card_header("Top 10 Markers Heatmap"),
              plotOutput("markers_heatmap")
            ),
            card(
              card_header("All Markers Table"),
              DT::DTOutput("markers_table")
            )
          )
        )
      )
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

  # Helper to format genes based on detected species
  format_genes <- function(gene_list, obj) {
    # Check if rownames are mostly uppercase (Human) or title case (Mouse)
    all_genes <- rownames(obj)
    # Sample a few genes to check
    sample_genes <- head(all_genes, 100)
    is_upper <- sum(grepl("^[A-Z0-9]+$", sample_genes)) > length(sample_genes) * 0.5

    if (is_upper) {
      return(toupper(gene_list))
    } else {
      # Title case: Uppercase first letter, lowercase the rest
      return(paste0(toupper(substr(gene_list, 1, 1)), tolower(substr(gene_list, 2, nchar(gene_list)))))
    }
  }

  # Preset observers
  observeEvent(input$preset_immune, {
    obj <- raw_obj()
    req(obj)
    genes <- format_genes(c("PTPRC", "CD3E", "CD19", "CD8A", "CD4", "MS4A1", "CD14", "FCGR3A"), obj)
    updateTextAreaInput(session, "genes", value = paste(genes, collapse = ", "))
  })

  observeEvent(input$preset_death, {
    obj <- raw_obj()
    req(obj)
    genes <- format_genes(c("BAX", "BAK1", "CASP3", "CASP8", "CASP9", "FAS"), obj)
    updateTextAreaInput(session, "genes", value = paste(genes, collapse = ", "))
  })

  observeEvent(input$preset_prolif, {
    obj <- raw_obj()
    req(obj)
    genes <- format_genes(c("MKI67", "TOP2A", "PCNA", "MCM2"), obj)
    updateTextAreaInput(session, "genes", value = paste(genes, collapse = ", "))
  })

  observeEvent(input$preset_cycle, {
    obj <- raw_obj()
    req(obj)
    genes <- format_genes(c("CCND1", "CCNE1", "CDK2", "CDK4"), obj)
    updateTextAreaInput(session, "genes", value = paste(genes, collapse = ", "))
  })

  # Render the Dot Plot
  output$dot_plot <- renderPlot({
    obj <- raw_obj()
    req(obj, input$resolution, input$genes)

    # Parse and clean gene list
    gene_list <- unlist(strsplit(input$genes, "[, \t\n\r]+"))
    gene_list <- gene_list[gene_list != ""]

    if (length(gene_list) == 0) {
      return(NULL)
    }

    # Validate genes exist in the object
    valid_genes <- intersect(gene_list, rownames(obj))

    if (length(valid_genes) == 0) {
      return(ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No valid genes found in dataset") +
        theme_void())
    }

    # Ensure the resolution is set as the active identity
    Idents(obj) <- input$resolution

    DotPlot(obj, features = valid_genes) +
      RotatedAxis() +
      labs(title = paste("Dot Plot - Resolution:", input$resolution))
  })

  # Reactive value for markers
  markers_data <- reactiveVal(NULL)

  # Run FindAllMarkers
  observeEvent(input$run_markers, {
    obj <- raw_obj()
    req(obj, input$resolution)

    withProgress(message = "Finding markers (this may take a minute)...", {
      # Ensure active identity matches selected resolution
      Idents(obj) <- input$resolution
      res <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      markers_data(res)
    })
  })

  # Render Markers Table
  output$markers_table <- DT::renderDT({
    req(markers_data())
    df <- markers_data() %>% filter(p_val_adj < 0.05)
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
  })

  # Render Markers Heatmap
  output$markers_heatmap <- renderPlot(
    {
      res <- markers_data()
      obj <- raw_obj()
      req(res, obj, input$resolution)

      # Get top 10 markers per cluster
      top10 <- res %>%
        group_by(cluster) %>%
        slice_max(n = 10, order_by = avg_log2FC)

      # Ensure active identity matches
      Idents(obj) <- input$resolution

      DoHeatmap(obj, features = top10$gene) + NoLegend()
    },
    height = function() {
      res <- markers_data()
      if (is.null(res)) {
        return(800)
      }
      top10 <- res %>%
        group_by(cluster) %>%
        slice_max(n = 10, order_by = avg_log2FC)
      # Approx 18px per gene + 200px for headers/labels
      max(800, nrow(top10) * 18 + 200)
    }
  )

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
