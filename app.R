library(shiny)
library(Seurat)
library(bslib)
library(ggplot2)
library(magrittr)
library(gtools)
library(DT)
library(dplyr)
library(clustifyr)

# Set max upload size to 2GB for large Seurat objects
options(shiny.maxRequestSize = 2000 * 1024^2)

ui <- page_sidebar(
  sidebar = sidebar(
    titlePanel("SOWhat?"),
    fileInput("file", "Upload Seurat Object (.rds)", accept = ".rds"),
    selectInput("resolution", "Clustering Resolution", choices = NULL),
    hr(),
    card(
      card_header("Object Statistics"),
      tableOutput("stats_table")
    ),
    hr()
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
              uiOutput("dot_preset_buttons")
            ),
            plotOutput("dot_plot", height = "800px")
          )
        ),
        nav_panel(
          "Module Scores",
          card_body(
            layout_column_wrap(
              width = 1,
              textAreaInput("module_genes", "Enter Gene List for Module Score",
                placeholder = "e.g., PTPRC, CD3E, CD19", rows = 3, width = "100%"
              ),
              uiOutput("module_preset_buttons")
            ),
            plotOutput("module_plot", height = "800px")
          )
        ),
        nav_panel(
          "Clustifyr",
          layout_sidebar(
            sidebar = sidebar(
              selectInput("ref_file", "Select Reference Matrix", choices = NULL),
              actionButton("run_clustifyr", "Run Clustifyr", class = "btn-primary w-100"),
              actionButton("apply_clustifyr", "Apply Top Cell Types", class = "btn-success w-100 mt-2"),
              hr(),
              helpText("Reference matrices should be .rds files in the 'refmats' folder.")
            ),
            plotOutput("clustifyr_plot", height = "800px")
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
        div(
          actionButton("reset_annotes", label = NULL, icon = icon("rotate-left"), class = "btn-sm btn-outline-danger"),
          downloadButton("download_annotes", "", class = "btn-sm")
        )
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

  # Reactive annotated object
  annotated_obj <- reactive({
    obj <- raw_obj()
    req(obj, input$resolution)

    # Set active resolution
    Idents(obj) <- input$resolution

    # Get current clusters
    clusters <- gtools::mixedsort(unique(as.character(obj@meta.data[[input$resolution]])))

    # Map old names to new names from inputs or defaults
    new_names <- sapply(clusters, function(cl) {
      id <- paste0("ann_", input$resolution, "_", cl)
      val <- input[[id]]
      if (is.null(val) || val == "") {
        # Fallback to reactiveValues if input isn't ready yet or empty
        val <- annotes$data[[id]]
      }
      if (is.null(val) || val == "") cl else val
    })

    names(new_names) <- clusters
    # RenameIdents handles named vector where names=old, values=new
    RenameIdents(obj, new_names)
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

    # Update refmats choices
    ref_files <- list.files("refmats", pattern = "\\.rds$", full.names = FALSE)
    updateSelectInput(session, "ref_file", choices = ref_files)
  })

  # Render the UMAP plot
  output$umap_plot <- plotly::renderPlotly({
    obj <- annotated_obj()
    req(obj)

    p <- DimPlot(obj, label = TRUE, raster = FALSE) +
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

  # Dynamic Preset Generation
  genefiles <- reactive({
    list.files("genelists", pattern = "\\.txt$", full.names = FALSE)
  })

  # Helper to load genes from file
  load_genes_from_file <- function(filename) {
    path <- file.path("genelists", filename)
    if (!file.exists(path)) {
      return(character(0))
    }
    content <- readLines(path, warn = FALSE)
    # Split by comma, space, or tab
    genes <- unlist(strsplit(content, "[, \t]+"))
    genes <- genes[genes != "" & genes != "."] # Clean up
    return(genes)
  }

  output$dot_preset_buttons <- renderUI({
    files <- genefiles()
    # Add hardcoded presets to the list of buttons
    hardcoded <- c("Immune", "Cell Death", "Proliferation", "Cell Cycle")

    div(
      class = "d-flex flex-wrap gap-2 mb-1",
      # Hardcoded first
      lapply(hardcoded, function(p) {
        id <- paste0("preset_", gsub(" ", "_", tolower(p)))
        actionButton(id, p, class = "btn-outline-primary btn-sm")
      }),
      # Then dynamic files
      lapply(files, function(f) {
        name <- tools::file_path_sans_ext(f)
        id <- paste0("dyn_dot_", name)
        actionButton(id, name, class = "btn-outline-secondary btn-sm")
      })
    )
  })

  output$module_preset_buttons <- renderUI({
    files <- genefiles()
    hardcoded <- c("Immune", "Cell Death", "Proliferation", "Cell Cycle")

    div(
      class = "d-flex flex-wrap gap-2 mb-1",
      # Hardcoded first
      lapply(hardcoded, function(p) {
        id <- paste0("module_preset_", gsub(" ", "_", tolower(p)))
        actionButton(id, p, class = "btn-outline-primary btn-sm")
      }),
      # Then dynamic files
      lapply(files, function(f) {
        name <- tools::file_path_sans_ext(f)
        id <- paste0("dyn_module_", name)
        actionButton(id, name, class = "btn-outline-secondary btn-sm")
      })
    )
  })

  # Observers for hardcoded presets
  observeEvent(input$preset_immune, {
    obj <- raw_obj()
    req(obj)
    genes <- format_genes(c("PTPRC", "CD3E", "CD3D", "CD4", "CD8A", "CD19", "MS4A1", "CD14", "FCGR3A", "LYZ", "HLA-DRA", "NCAM1", "CD1C", "CD68", "MZB1", "PPBP"), obj)
    updateTextAreaInput(session, "genes", value = paste(genes, collapse = ", "))
  })

  observeEvent(input$preset_cell_death, {
    obj <- raw_obj()
    req(obj)
    genes <- format_genes(c("BAX", "BAK1", "CASP3", "CASP8", "CASP9", "FAS"), obj)
    updateTextAreaInput(session, "genes", value = paste(genes, collapse = ", "))
  })

  observeEvent(input$preset_proliferation, {
    obj <- raw_obj()
    req(obj)
    genes <- format_genes(c("MKI67", "TOP2A", "PCNA", "MCM2"), obj)
    updateTextAreaInput(session, "genes", value = paste(genes, collapse = ", "))
  })

  observeEvent(input$preset_cell_cycle, {
    obj <- raw_obj()
    req(obj)
    genes <- format_genes(c("CCND1", "CCNE1", "CDK2", "CDK4"), obj)
    updateTextAreaInput(session, "genes", value = paste(genes, collapse = ", "))
  })

  observeEvent(input$module_preset_immune, {
    obj <- raw_obj()
    req(obj)
    genes <- format_genes(c("PTPRC", "CD3E", "CD3D", "CD4", "CD8A", "CD19", "MS4A1", "CD14", "FCGR3A", "LYZ", "HLA-DRA", "NCAM1", "CD1C", "CD68", "MZB1", "PPBP"), obj)
    updateTextAreaInput(session, "module_genes", value = paste(genes, collapse = ", "))
  })

  observeEvent(input$module_preset_cell_death, {
    obj <- raw_obj()
    req(obj)
    genes <- format_genes(c("BAX", "BAK1", "CASP3", "CASP8", "CASP9", "FAS"), obj)
    updateTextAreaInput(session, "module_genes", value = paste(genes, collapse = ", "))
  })

  observeEvent(input$module_preset_proliferation, {
    obj <- raw_obj()
    req(obj)
    genes <- format_genes(c("MKI67", "TOP2A", "PCNA", "MCM2"), obj)
    updateTextAreaInput(session, "module_genes", value = paste(genes, collapse = ", "))
  })

  observeEvent(input$module_preset_cell_cycle, {
    obj <- raw_obj()
    req(obj)
    genes <- format_genes(c("CCND1", "CCNE1", "CDK2", "CDK4"), obj)
    updateTextAreaInput(session, "module_genes", value = paste(genes, collapse = ", "))
  })

  # Dynamic observers
  observe({
    files <- genefiles()
    obj <- raw_obj()
    req(obj)

    for (f in files) {
      name <- tools::file_path_sans_ext(f)

      # Dot Plot Observer
      local({
        filename <- f
        id_dot <- paste0("dyn_dot_", tools::file_path_sans_ext(filename))
        observeEvent(input[[id_dot]], {
          genes_raw <- load_genes_from_file(filename)
          genes <- format_genes(genes_raw, obj)
          updateTextAreaInput(session, "genes", value = paste(genes, collapse = ", "))
        })

        # Module Score Observer
        id_mod <- paste0("dyn_module_", tools::file_path_sans_ext(filename))
        observeEvent(input[[id_mod]], {
          genes_raw <- load_genes_from_file(filename)
          genes <- format_genes(genes_raw, obj)
          updateTextAreaInput(session, "module_genes", value = paste(genes, collapse = ", "))
        })
      })
    }
  })

  # Render the Dot Plot
  output$dot_plot <- renderPlot({
    obj <- annotated_obj()
    req(obj, input$genes)

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

    DotPlot(obj, features = valid_genes) +
      RotatedAxis() +
      labs(title = paste("Dot Plot - Resolution:", input$resolution))
  })

  # Render the Module Plot
  output$module_plot <- renderPlot({
    obj <- annotated_obj()
    req(obj, input$module_genes)

    # Parse and clean gene list
    gene_list <- unlist(strsplit(input$module_genes, "[, \t\n\r]+"))
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

    # Add Module Score
    # We use a unique name to avoid collisions
    score_name <- "UserModule"
    obj <- AddModuleScore(obj, features = list(valid_genes), name = score_name)

    # VlnPlot uses the name + "1"
    VlnPlot(obj, features = paste0(score_name, "1")) +
      NoLegend() +
      labs(title = "Module Score Distribution", y = "Score")
  })

  # Reactive value for clustifyr results
  clustifyr_res <- reactiveVal(NULL)

  # Run Clustifyr
  observeEvent(input$run_clustifyr, {
    obj <- annotated_obj()
    req(obj, input$ref_file)

    withProgress(message = "Running Clustifyr analysis...", {
      ref_path <- file.path("refmats", input$ref_file)
      ref_mat <- readRDS(ref_path)

      # Clustifyr metadata formatting fix:
      # Explicitly add current annotations to the metadata
      obj$current_annotations <- as.character(Idents(obj))

      # Handle Gene Name Overlap (Case-Insensitive matching)
      # Normalize both query and reference genes to uppercase for better matching
      query_genes_orig <- rownames(obj)
      ref_genes_orig <- rownames(ref_mat)

      # Create a version of ref_mat with uppercase rownames
      ref_mat_upper <- ref_mat
      rownames(ref_mat_upper) <- toupper(ref_genes_orig)

      # Determine if we should also uppercase the query genes
      # (If the overlap is low, we'll try uppercase)
      if (length(intersect(query_genes_orig, ref_genes_orig)) < 10) {
        # Create a temporary expression matrix with uppercase gene names for clustifyr
        # Instead of modifying the Seurat object (which is heavy), we can pass the data directly
        # or just hope clustifyr handles the query_genes argument well.
        # However, the most robust way is to provide a ref_mat that matches the query nomenclature.

        # Test if query is TitleCase (Mouse)
        is_mouse <- any(grepl("^[A-Z][a-z]", head(query_genes_orig, 50)))
        if (is_mouse) {
          # Format ref genes to TitleCase
          rownames(ref_mat) <- paste0(
            toupper(substr(ref_genes_orig, 1, 1)),
            tolower(substr(ref_genes_orig, 2, nchar(ref_genes_orig)))
          )
        } else {
          # Default to everything UPPERCASE
          rownames(ref_mat) <- toupper(rownames(ref_mat))
        }
      }

      # Run clustify
      res <- clustify(
        input = obj,
        ref_mat = ref_mat,
        cluster_col = "current_annotations",
        obj_out = FALSE
      )
      clustifyr_res(res)
    })
  })

  # Render Clustifyr Heatmap
  output$clustifyr_plot <- renderPlot({
    res <- clustifyr_res()
    req(res)

    plot_cor_heatmap(
      cor_mat = res,
      cluster_rows = TRUE,
      cluster_columns = TRUE
    )
  })

  # Apply Top Cell Types from Clustifyr
  observeEvent(input$apply_clustifyr, {
    res <- clustifyr_res()
    req(res)

    # Iterate through current clusters
    clusters <- rownames(res)
    for (cl in clusters) {
      # find column with max correlation
      best_hit <- colnames(res)[which.max(res[cl, ])]

      # Determine input ID
      id <- paste0("ann_", input$resolution, "_", cl)

      # Update reactive data and UI
      annotes$data[[id]] <- best_hit
      updateTextInput(session, id, value = best_hit)
    }
  })

  # Reactive value for markers
  markers_data <- reactiveVal(NULL)

  # Run FindAllMarkers
  observeEvent(input$run_markers, {
    obj <- annotated_obj()
    req(obj)

    withProgress(message = "Finding markers (this may take a minute)...", {
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
      obj <- annotated_obj()
      req(res, obj)

      # Get top 10 markers per cluster
      top10 <- res %>%
        group_by(cluster) %>%
        slice_max(n = 10, order_by = avg_log2FC)

      DoHeatmap(obj, features = top10$gene, angle = 90) + NoLegend()
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
          tags$th("Notes")
        )
      ),
      tags$tbody(
        lapply(clusters, function(cl) {
          # Use namespaced IDs based on resolution and cluster
          id_ann <- paste0("ann_", input$resolution, "_", cl)
          id_reason <- paste0("reason_", input$resolution, "_", cl)

          # Retrieve existing values if any, default to cluster ID if set
          val_ann <- if (!is.null(annotes$data[[id_ann]]) && annotes$data[[id_ann]] != "") annotes$data[[id_ann]] else cl
          val_reason <- if (!is.null(annotes$data[[id_reason]])) annotes$data[[id_reason]] else ""

          tags$tr(
            tags$td(cl),
            tags$td(textInput(id_ann, NULL, value = val_ann, width = "100%", placeholder = "Enter annotation...")),
            tags$td(textInput(id_reason, NULL, value = val_reason, width = "100%", placeholder = "Enter notes..."))
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
      Metric = c("Object", "Cells", "Features", "Assays"),
      Value = c(
        if (is.null(input$SO) || length(input$SO) == 0) "N/A" else as.character(input$SO),
        as.character(ncol(obj)),
        as.character(nrow(obj)),
        paste(Assays(obj), collapse = ", ")
      )
    )
  })

  # Reset Annotations logic
  observeEvent(input$reset_annotes, {
    showModal(modalDialog(
      title = "Reset Cluster Annotations?",
      "This will revert all annotations for the current resolution back to their original cluster numbers. This cannot be undone.",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_reset", "Reset", class = "btn-danger")
      ),
      easyClose = TRUE
    ))
  })

  observeEvent(input$confirm_reset, {
    obj <- raw_obj()
    req(obj, input$resolution)

    clusters <- gtools::mixedsort(unique(as.character(obj@meta.data[[input$resolution]])))

    for (cl in clusters) {
      id_ann <- paste0("ann_", input$resolution, "_", cl)
      # Revert to cluster ID
      annotes$data[[id_ann]] <- cl
      updateTextInput(session, id_ann, value = cl)
    }

    removeModal()
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
          val <- if (!is.null(input[[id]])) input[[id]] else annotes$data[[id]]
          if (is.null(val) || val == "") cl else val
        }),
        Reason = sapply(clusters, function(cl) {
          id <- paste0("reason_", input$resolution, "_", cl)
          val <- if (!is.null(input[[id]])) input[[id]] else annotes$data[[id]]
          if (is.null(val)) "" else val
        }),
        stringsAsFactors = FALSE
      )

      write.csv(export_data, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
