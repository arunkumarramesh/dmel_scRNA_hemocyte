library(shiny)
library(Seurat)
library(ggplot2)

load(file="integrated_hemocytes.Robj")

ui <- fluidPage(
  headerPanel('Single cell hemocyte D.mel'),
  sidebarPanel(
    textAreaInput('gene', 'Input gene names one at a time', value = "", width = NULL, placeholder = 'e.g. FBgn0032422'),
    actionButton("do", "Evaluate!")
  ),
  mainPanel(plotOutput("plot1")
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



