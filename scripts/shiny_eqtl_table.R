# Load necessary libraries
library(shiny)
library(DT)
library(data.table)
library(ggplot2)

# UI
ui <- fluidPage(
  titlePanel("LBD CWOW GWAS and eQTL"),
  tabsetPanel(
    tabPanel("Home",
             fluidRow(
               column(12,
                      h3(""),
                      p("This application allows you to explore eQTL data for different disease conditions, including Alzheimer's disease, Lewy body disease, and Pathological amyloid, compared to controls. You can search for specific SNPs or genes and visualize their associations in cis or trans eQTLs. The table highlights the significance (FDR) and effect size (beta). 
                        Created by Dr. Kimberly Olney in the lab of Dr. John Fryer at Mayo Clinic in Scottsdale Arizona")
               )
             )
    ),
    tabPanel("eQTL Viewer",
             sidebarLayout(
               sidebarPanel(
                 selectInput("dataset", "Select Dataset:", 
                             choices = c("Alzheimer's disease (AD) and controls" = "AD_control", 
                                         "Lewy body disease (LBD) and controls" = "LBD_control", 
                                         "Pathological amyloid and controls" = "PA_control")),
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
    dt <- dt[, .(SNP, gene, beta, FDR)]  # Remove t-test and p-value columns
    return(dt)
  })
  
  # Reactive for filtering data based on search inputs
  filteredData <- reactive({
    dt <- loadData()
    if (input$search_snp != "") {
      dt <- dt[grepl(input$search_snp, SNP)]
    }
    if (input$search_gene != "") {
      dt <- dt[grepl(input$search_gene, gene)]
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
