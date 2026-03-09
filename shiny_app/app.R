# Load necessary libraries
library(shiny)
library(DT)
library(data.table)
library(ggplot2)
options(rsconnect.max.bundle.size = 5 * 1024^3)  # Set to 3.5 GB
#rsconnect.max.bundle.size = 6.5 * 1024^3  # Set to 3.5 GB

#getOption("rsconnect.max.bundle.size")

#setwd("/research/labs/neurology/fryer/m239830/LBD_CWOW/QTL/LBD_CWOW_QTL/shiny_app")
# Load the gene annotations
genes <- fread("genes_filtered.txt")
genes <- genes[, .(gene_id, gene_name)]

# UI
ui <- fluidPage(
  titlePanel("LBD CWOW eQTL"),
  tabsetPanel(
    tabPanel("Home",
             fluidRow(
               column(12,
                      h3(""),
                      p("Explore eQTL data for different disease conditions"),
                      p("The eQTL results are derived from first subsetting the data to only include a disease type and control samples. Model used is modelLINEAR_CROSS to  include the interaction of SNP and the last covariate in the model (0 = control, and 1 = disease) and test for signification associations."),
                      p("Trans eQTLs are limited to those with a p-value less than 1e-5. All cis eQTLs are reported, cis distance is defined as 1e6.")
               )
             )
    ),
    tabPanel("eQTL Viewer",
             sidebarLayout(
               sidebarPanel(
                 selectInput("dataset", "Select Dataset:", 
                             choices = c("Alzheimer's disease (AD) and controls" = "AD_control", 
                                         "Lewy body disease (LBD) and controls" = "LBD_control", 
                                         "Lewy body disease LBD(S) and controls" = "LBD_S_control", 
                                         "Lewy body disease LBD(AS) and controls" = "LBD_AS_control", 
                                         "Lewy body disease LBD(ATS) and controls" = "LBD_ATS_control", 
                                         "Pathological amyloid (PA) and controls" = "PA_control")),
                 radioButtons("eqtl_type", "Select eQTL Type:", 
                              choices = c("cis", "trans")),
                 textInput("search_snp", "Search SNP ID:", ""),
                 textInput("search_gene", "Search Gene:", "")
               ),
               mainPanel(
                 DTOutput("eqtl_table")
               )
             )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Reactive to load the selected dataset
  loadData <- reactive({
    file_path <- paste0(input$eqtl_type, "_eQTL_", input$dataset)
    dt <- fread(file_path)
    dt <- dt[, .(SNP, gene, beta, FDR)]  # Ensure 'gene' is the column with gene_id
    dt <- merge(dt, genes, by.x = "gene", by.y = "gene_id", all.x = TRUE)  # Merge to get gene_name
    dt <- dt[, .(SNP, gene_name, beta, FDR)]  # Keep only relevant columns
    dt <- dt[order(FDR)]  # Order by smallest FDR
    return(dt)
  })
  
  # Reactive for filtering data based on search inputs
  filteredData <- reactive({
    dt <- loadData()
    if (input$search_snp != "") {
      dt <- dt[grepl(input$search_snp, SNP)]
    }
    if (input$search_gene != "") {
      dt <- dt[grepl(input$search_gene, gene_name)]
    }
    return(dt)
  })
  
  # Render the data table with color gradients
  output$eqtl_table <- renderDT({
    dt <- filteredData()
    
    # Define breakpoints and color gradients
    brks <- seq(-1000, 1000, 100)
    clrs <- colorRampPalette(c("#6baed6", "white", "red"))(length(brks) + 1)
    
    brks2 <- seq(0, .001, .000001)
    clrs2 <- colorRampPalette(c("green", "lightgreen", "white"))(length(brks2) + 1)
    
    datatable(dt, options = list(pageLength = 10, dom = 't', autoWidth = TRUE, rownames = FALSE)) %>%
      formatStyle('FDR', 
                  backgroundColor = styleInterval(brks2, clrs2)) %>%
      formatStyle('beta', 
                  backgroundColor = styleInterval(brks, clrs))
  })
}

# Run the app
shinyApp(ui, server)
