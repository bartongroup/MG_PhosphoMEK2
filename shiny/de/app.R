### DE PHOSPHO EXPLORER

libDir <- "/cluster/gjb_lab/mgierlinski/R_shiny/library/4.1"
if(dir.exists(libDir)) .libPaths(libDir)

library(shiny)
library(ggbeeswarm)
library(cowplot)
library(tidyverse)
library(DT)
source("func.R")

css <- "table{font-size: 11px; background-color: #EAF5FF}"

### Read data ###

data <- read_rds("../data_de.rds")
contrasts <- unique(data$de$contrast)
all_genes <- data$pho2gene$gene_name %>% unique()

max_points <- 500


#######################################################################

ui <- shinyUI(fluidPage(
  
  tags$style(css),
  
  titlePanel("MEKi phosphoproteomics"),

  fluidRow(
    column(12,
      fluidRow(
        column(4, 
          radioButtons("contrast", "Contrast:", choices = contrasts, inline = TRUE),
          radioButtons("plotType", "Plot type:", choices = c("Volcano", "MA"), inline=TRUE),
          plotOutput("mainPlot", height="480px", width="100%", brush="plot_brush", hover="plot_hover")
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
          radioButtons("enrichment", "Enrichment:", choices = c("GO", "Reactome", "KEGG"), inline=TRUE),
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
    
    if(input$plotType == "Volcano") {
      xy_data$x <- xy_data$logFC
      xy_data$y <- -log10(xy_data$PValue)
    } else if(input$plotType == "MA") {
      xy_data$x <- xy_data$AveExpr
      xy_data$y <- xy_data$logFC
    }
    xy_data
  }
  
  select_phospho <- function(max_hover=1) {
    xy_data <- get_xy_data()
    sel <- NULL
    tab_idx <- as.numeric(input$allPhosphoTable_rows_selected)
    if(!is.null(input$plot_brush)){
      brushed <- na.omit(brushedPoints(xy_data, input$plot_brush))
      sel <- brushed$id
    } else if(!is.null(input$plot_hover)) {
      near <- nearPoints(xy_data, input$plot_hover, threshold = 20, maxpoints = max_hover)
      sel <- near$id
    } else if(length(tab_idx) > 0) {
      sel <- xy_data[tab_idx, ] %>% pull(id)
    }
    return(sel)
  }
  
  select_peptides <- function(phospho_ids) {
    data$pho2pep %>% 
      filter(id %in% phospho_ids) %>% 
      pull(peptide_ids) %>% 
      unique()
  }
  
  select_proteins <- function(phospho_ids) {
    data$pho2pro %>% 
      filter(id %in% phospho_ids) %>% 
      pull(protein_ids) %>% 
      unique()
  }
  
  output$peptideInfo <- renderTable({
    xy_data <- get_xy_data()
    phospho_ids <- select_phospho()
    df <- NULL
    if (!is.null(phospho_ids) && length(phospho_ids) >= 1 && length(phospho_ids) <= max_points) {
      peptide_ids <- select_peptides(phospho_ids)
      df <- data$pep$info %>% filter(id %in% peptide_ids) %>% 
        arrange(gene_name) %>% 
        select(sequence, protein, gene_name, start_position, end_position) %>% 
        mutate_at(vars(start_position, end_position), as.integer)
    } else if (length(phospho_ids) > max_points) {
      df <- data.frame(Error=paste0('only ',max_points,' points can be selected.'))
    }
    df
  })

  enrichmentTable <- function(terms) {
    xy_data <- get_xy_data()
    phospho_ids <- NULL
    fe <- NULL
    if(!is.null(input$plot_brush)){
      brushed <- na.omit(brushedPoints(xy_data, input$plot_brush))
      phospho_ids <- brushed$id
      sel_genes <- data$pho2gene %>% 
        filter(id %in% phospho_ids) %>% 
        pull(gene_name)
      n <- length(sel_genes)
      if(n > 2 && n <= max_points) {
        fe <- sh_functional_enrichment(all_genes, sel_genes, terms)
      } else if(n > 2) {
        fe <- data.frame(Error=paste0('only ',max_points,' points can be selected.'))
      }
    }
    fe
  }
  
  output$Enrichment <- renderTable({
    if(input$enrichment == "GO") {
      d <- data$go
    } else if(input$enrichment == "Reactome") {
      d <- data$reactome
    } else if(input$enrichment == "KEGG") {
      d <- data$kegg
    }
    enrichmentTable(d)
  })
  
  
  output$geneInfo <- renderText({
    phospho_ids <- select_phospho()
    if(!is.null(phospho_ids) && length(phospho_ids) == 1) {
      protein_id <- select_proteins(phospho_ids)
      data$pro$info %>% 
        filter(id == protein_id) %>% 
        pull(protein_names)
    }
  })
  
  output$peptideSeq <- renderPlot({
    phospho_ids <- select_phospho()
    if(!is.null(phospho_ids) && length(phospho_ids) == 1) {
      sh_plot_pepseq(data$pep, data$pho, data$de, phospho_ids)
    }
  })
  
  output$intensityPlot <- renderPlot({
    phospho_ids <- select_phospho()
    if(!is.null(phospho_ids) && length(phospho_ids) == 1) {
      protein_id <- select_proteins(phospho_ids)
      plot_grid(
        sh_plot_intensities(data$pho, phospho_ids, tit="Phospho site", log_scale=FALSE),
        sh_plot_intensities(data$pro, protein_id[1], tit="Protein", log_scale=FALSE),
        align = "h"
      )
    }
  })
  
  output$prophoPlot <- renderPlot({
    phospho_ids <- select_phospho()
    if(!is.null(phospho_ids) && length(phospho_ids) == 1) {
      protein_id <- select_proteins(phospho_ids)
      sh_plot_full_protein(data$de, data$pro, protein_id[1], phospho_ids, input$contrast)
    }
  })

  output$mainPlot <- renderPlot({
    xy_data <- get_xy_data()
    tab_idx <- as.numeric(input$allPhosphoTable_rows_selected)
    
    if(input$plotType == "Volcano") {
      g <- sh_plot_volcano(xy_data)
    } else {
      g <- sh_plot_ma(xy_data)
    }
    if(length(tab_idx) >= 1) {
      g <- g + geom_point(data=xy_data[tab_idx, ], colour="red", size=2)
    }
    g
  })

  output$allPhosphoTable <- DT::renderDataTable({
    xy_data <- get_xy_data()
    all_table(filter(xy_data, contrast == input$contrast))
  })
}

# Run the application
shinyApp(ui = ui, server = server)



# input <- list(contrast = "CTRL-LPSN", plotType = "Volcano")