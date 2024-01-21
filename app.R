library(shiny)
library(Seurat)
library(ggplot2)

#load(file="../GSE148826_intergrated_hemocytes.Robj")
#integrated_hemocytes <- DietSeurat(integrated_hemocytes, counts = F,misc = F)
#integrated_hemocytes <- integrated_hemocytes[,sample(1:19344, 15000, replace=F)]
#save(integrated_hemocytes,file="integrated_hemocytes.Robj")

load(file="integrated_hemocytes.Robj")

ui <- fluidPage(
  headerPanel('Single cell RNA-seq of Drosophila larval hemocytes'),
  sidebarPanel(
    textAreaInput('gene', 'Input gene names one at a time', value = "", width = NULL, placeholder = 'e.g. FBgn0032422'),
    actionButton("do", "Evaluate!")
  ),
  mainPanel(
    p("Goal is to identify genes that are constitutively activated following long-term parasite exposure and those that are transiently induced by parasite infection"),
    p("Input gene name as Flybase gene ID, which can be obtained from https://flybase.org/"),
    p("PLASM1-2: plasmatocytes, LAM1-2: immature lamellocytes, LAM3: mature lamellocytes, CC: crystal cells, MET: metabolism, AMP: anti-microbial peptide"),
    p("Selection: 26 generations of selection under high parasite (Leptopilina boulardi) pressure"),
    p("No Selection: 26 generations under standard conditions without parasite exposure"),
    p("Infection: Infection with parasitic wasp"),
    p("No Infection: No infection control"),
    p("Expression level: Log-normalised read counts from Seurat V3"),
    p("Full data available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148826"),
    p("Please post questions at https://github.com/arunkumarramesh/dmel_scRNA_hemocyte"),
    plotOutput("plot1")
  )
)

server <- function(input, output,session) {
  gene <- eventReactive(input$do, {
    unlist(strsplit(as.character(input$gene), ',', fixed=TRUE))
  }, ignoreNULL= T)
  p <- reactive({
    VlnPlot(integrated_hemocytes,features=c(gene()),split.by = "poptreat",pt=0,adjust=0)+
      geom_boxplot()
  })
  observeEvent(input$do, {
    output$plot1 <- renderPlot({
      print(p())
    })
  })
}

shinyApp(ui, server)
