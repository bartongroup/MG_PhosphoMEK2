### DE PHOSPHO EXPLORER

library(shiny)
library(ggbeeswarm)
library(cowplot)
library(tidyverse)
library(DT)
library(fenr)
source("func.R")

select <- dplyr::select

css <- "table{font-size: 11px; background-color: #EAF5FF}"

### Read data ###

data <- read_rds("../data_de.rds")
contrasts <- unique(data$de$contrast)
all_genes <- data$pho2gene$gene_name %>% unique()

max_points <- 500


#######################################################################

ui <- shinyUI(fluidPage(
  
  tags$style(css),
  
  titlePanel("Phosphoproteomic identification of ERK1/2 targets in human neuromesodermal progenitors"),

  p("This study used phospho-proteomics to identify phospho-dynamic proteins following acute inhibition of ERK1/2 with the MEK1/2 inhibitor PD0325901 in human embryonic stem cell-derived neuromesodermal progenitors."),

  fluidRow(
    column(12,
      fluidRow(
        column(4, 
          radioButtons("contrast", "Contrast:", choices = contrasts, inline = TRUE),
          radioButtons("plotType", "Plot type:", choices = c("Volcano", "MA"), inline = TRUE),
          plotOutput("mainPlot", height = "480px", width = "100%", brush = "plot_brush", hover = "plot_hover")
        ),
        column(3,
          textOutput("geneInfo"),
          plotOutput("prophoPlot", height = "150px", width = "100%"),
          plotOutput("peptideSeq", height = "300px", width = "100%"),
          br(),
          plotOutput("intensityPlot", height = "200px", width = "100%")
        ),
        column(5,
          p("Peptide list"),
          div(style = 'height: 200px; overflow-y: scroll', tableOutput("peptideInfo")),
          br(),
          radioButtons("enrichment", "Enrichment:", choices = c("GO-CC", "GO-BP", "GO-MF", "Reactome", "KEGG"), inline = TRUE),
          div(style = 'height: 400px; overflow-y: scroll', tableOutput("Enrichment"))
        )
      ),
      fluidRow(
        DT::dataTableOutput("allPhosphoTable")
      )
    )
  )
)
)


########################################################################################


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  get_xy_data <- function() {
    ctr <- input$contrast
    xy_data <- data$de %>% 
      filter(contrast == ctr)
    
    if (input$plotType == "Volcano") {
      xy_data$x <- xy_data$logFC
      xy_data$y <- -log10(xy_data$PValue)
    } else if (input$plotType == "MA") {
      xy_data$x <- xy_data$AveExpr
      xy_data$y <- xy_data$logFC
    }
    xy_data
  }
  
  select_phospho <- function(max_hover = 1) {
    xy_data <- get_xy_data()
    sel <- NULL
    tab_idx <- as.numeric(input$allPhosphoTable_rows_selected)
    if (!is.null(input$plot_brush)) {
      brushed <- na.omit(brushedPoints(xy_data, input$plot_brush))
      sel <- brushed %>% 
        select(id, multi)
    } else if (!is.null(input$plot_hover)) {
      near <- nearPoints(xy_data, input$plot_hover, threshold = 20, maxpoints = max_hover)
      sel <- near %>% 
        select(id, multi)
    } else if (length(tab_idx) > 0) {
      sel <- xy_data[tab_idx, ] %>% 
        select(id, multi)
    }
    return(sel)
  }
  
  select_peptides <- function(pho_sel) {
    data$pho2pep %>% 
      filter(id %in% pho_sel$id) %>% 
      pull(peptide_ids) %>% 
      unique()
  }
  
  select_proteins <- function(pho_sel) {
    data$pho2pro %>% 
      filter(id %in% pho_sel$id) %>% 
      pull(protein_ids) %>% 
      unique()
  }
  
  
  
  output$peptideInfo <- renderTable({
    xy_data <- get_xy_data()
    pho_sel <- select_phospho()
    if (is.null(pho_sel)) return(NULL)
    df <- NULL
    if (!is.null(pho_sel) && nrow(pho_sel) >= 1 && nrow(pho_sel) <= max_points) {
      peptide_ids <- select_peptides(pho_sel)
      df <- data$pep$info %>% filter(id %in% peptide_ids) %>% 
        arrange(gene_name) %>% 
        select(sequence, protein, gene_name, start_position, end_position) %>% 
        mutate_at(vars(start_position, end_position), as.integer)
    } else if (nrow(pho_sel) > max_points) {
      df <- data.frame(Error = paste0('only ',max_points,' points can be selected.'))
    }
    df
  })

  enrichmentTable <- function(fterms) {
    xy_data <- get_xy_data()
    pho_sel <- NULL
    fe <- NULL
    if (!is.null(input$plot_brush)) {
      brushed <- na.omit(brushedPoints(xy_data, input$plot_brush))
      pho_sel <- brushed %>% select(id, multi)
      sel_genes <- data$pho2gene_first %>% 
        filter(id %in% pho_sel$id) %>% 
        pull(gene_name) %>% 
        unique()
      n <- length(sel_genes)
      if (n > 2 && n <= max_points) {
        fe <- sh_functional_enrichment(all_genes, sel_genes, fterms)
      } else if (n > 2) {
        fe <- data.frame(Error = paste0('only ',max_points,' points can be selected.'))
      }
    }
    fe
  }
  
  output$Enrichment <- renderTable({
    if (input$enrichment == "GO-CC") {
      d <- data$go_cc
    } else if (input$enrichment == "GO-BP") {
      d <- data$go_bp
    } else if (input$enrichment == "GO-MF") {
      d <- data$go_mf
    } else if (input$enrichment == "Reactome") {
      d <- data$reactome
    } else if (input$enrichment == "KEGG") {
      d <- data$kegg
    }
    enrichmentTable(d)
  })
  
  
  output$geneInfo <- renderText({
    pho_sel <- select_phospho()
    if (!is.null(pho_sel) && nrow(pho_sel) == 1) {
      protein_id <- select_proteins(pho_sel)
      data$pro$info %>% 
        filter(id %in% protein_id) %>% 
        pull(protein_names)
    }
  })
  
  output$peptideSeq <- renderPlot({
    pho_sel <- select_phospho()
    if (!is.null(pho_sel) && nrow(pho_sel) == 1) {
      sh_plot_pepseq(data$pep, data$pho, data$de, pho_sel)
    }
  })
  
  output$intensityPlot <- renderPlot({
    pho_sel <- select_phospho()
    if (!is.null(pho_sel) && nrow(pho_sel) == 1) {
      protein_id <- select_proteins(pho_sel)
      plot_grid(
        sh_plot_intensities(data$pho, pho_sel, tit = "Phospho site", log_scale = FALSE),
        sh_plot_intensities(data$pro, protein_id[1], tit = "Protein", log_scale = FALSE),
        align = "h"
      )
    }
  })
  
  output$prophoPlot <- renderPlot({
    pho_sel <- select_phospho()
    if (!is.null(pho_sel) && nrow(pho_sel) == 1) {
      protein_ids <- select_proteins(pho_sel)
      #print(pho_sel)
      #print(protein_ids)
      sh_plot_full_protein(data$de, data$pro, protein_ids, pho_sel, input$contrast)
    }
  })

  output$mainPlot <- renderPlot({
    xy_data <- get_xy_data()
    tab_idx <- as.numeric(input$allPhosphoTable_rows_selected)
    
    if (input$plotType == "Volcano") {
      g <- sh_plot_volcano(xy_data)
    } else {
      g <- sh_plot_ma(xy_data)
    }
    if (length(tab_idx) >= 1) {
      g <- g + geom_point(data = xy_data[tab_idx, ], colour = "red", size = 2)
    }
    g
  })

  output$allPhosphoTable <- DT::renderDataTable({
    xy_data <- get_xy_data()
    all_table(xy_data)
  })
}

# Run the application
shinyApp(ui = ui, server = server)



# input <- list(contrast = "CTRL-LPSN", plotType = "Volcano")