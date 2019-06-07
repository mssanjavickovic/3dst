## ui.R ##
# Make functions
library(shiny)
library(shinydashboard)
options(rgl.useNULL=TRUE)
library(shinyRGL)
library(rgl)
library(plot3Drgl)
library(akima)

### Read in data and functions
source('global.R')

### Sidebar Design
sidebar <- ## Sidebar content
  dashboardSidebar(
    sidebarLayout(
      sidebarMenu(
        menuItem("RA1", tabName = "RA1", icon = icon("cube")),
        menuItem("RA2", tabName = "RA2", icon = icon("cube")),
        mainPanel(selectInput('min','Min value for Log(normalized expression)', choices = c(-2:12), selected = "-2"),
                  selectInput('max','Max value for Log(normalized expression)', choices = c(12:-2), selected = "6"))
      ),
      uiOutput("conds")
    )
  )

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "RA1", fluidRow(column(1),
                                      column(8,box(selectInput("Gene_1",label = "Choose 1st gene to display", multiple = F, choices = unique(t1_names, t2_names), selected = "CCL19"),
                                                   rglwidgetOutput("visData.1", height = 400, width = "90%")),
                                             fluidRow(box(selectInput("Gene_2",label = "Choose 2nd gene to display", multiple = F, choices = unique(t1_names, t2_names), selected = "CD52"),
                                                          rglwidgetOutput("visData.2", height = 400, width = "90%"))))),
            fluidRow(column(1),column(8,(box(selectInput("Gene_3",label = "Choose 3rd gene to display", multiple = F, choices = unique(t1_names, t2_names), selected = "MS4A1"),
                                             rglwidgetOutput("visData.3", height = 400, width = "90%"))),
                                      fluidRow(box(selectInput("Gene_4",label = "Choose 4th gene to display", multiple = F, choices = unique(t1_names, t2_names), selected = "SSR4"),
                                                   rglwidgetOutput("visData.4", height = 400, width = "90%")))))),
    tabItem(tabName = "RA2", fluidRow(column(1),
                                      column(8,box(selectInput("GeneRA2_1",label = "Choose 1st gene to display", multiple = F, choices = unique(t1_names, t2_names), selected = "CCL19"),
                                                   rglwidgetOutput("visData.5", height = 400, width = "90%")),
                                             fluidRow(box(selectInput("GeneRA2_2",label = "Choose 2nd gene to display", multiple = F, choices = unique(t1_names, t2_names), selected = "CD52"),
                                                          rglwidgetOutput("visData.6", height = 400, width = "90%"))))),
            fluidRow(column(1),
                     column(8,(box(selectInput("GeneRA2_3",label = "Choose 3rd gene to display", multiple = F, choices = unique(t1_names, t2_names), selected = "MS4A1"),
                                   rglwidgetOutput("visData.7", height = 400, width = "90%"))),
                            fluidRow(box(selectInput("GeneRA2_4",label = "Choose 4th gene to display", multiple = F, choices = unique(t1_names, t2_names), selected = "SSR4"),
                                         rglwidgetOutput("visData.8", height = 400, width = "90%")))))))
)

# Put them together into a dashboardPage
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "3DSeTH: Spatial Transcriptomics for Heatmap gene expression visualization",titleWidth = 850),
  sidebar,
  body
)

server <- function(input, output, session){
  output$conds <- renderUI({
    con <- c(input$RA1,input$RA2)
  })
while (rgl.cur() > 0) { rgl.close() }
  output$visData.1 <- renderRglwidget({
    validate(
      need(!input$Gene_1 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.4(sample = "RA1", gene=input$Gene_1, t1.1, t1.2, t1.3, t1.4, s1.1,s1.2,s1.3,s1.4, x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA1_",input$Gene_1,".png", sep=""))
    rglwidget()})
  output$visData.2 <- renderRglwidget({
    validate(
      need(!input$Gene_2 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.4(sample = "RA1",gene=input$Gene_2,t1.1, t1.2, t1.3, t1.4, s1.1,s1.2,s1.3,s1.4, x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA1_",input$Gene_2,".png", sep=""))
    rglwidget()})
  output$visData.3 <- renderRglwidget({
    validate(
      need(!input$Gene_3 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.4(sample = "RA1",gene=input$Gene_3,t1.1, t1.2, t1.3, t1.4, s1.1,s1.2,s1.3,s1.4, x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA1_",input$Gene_3,".png", sep=""))
    rglwidget()})
  output$visData.4 <- renderRglwidget({
    validate(
      need(!input$Gene_4 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.4(sample = "RA1",gene=input$Gene_4,t1.1, t1.2, t1.3, t1.4, s1.1,s1.2,s1.3,s1.4, x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA1_",input$Gene_4,".png", sep=""))
    rglwidget()})
  output$visData.5 <- renderRglwidget({
    validate(
      need(!input$GeneRA2_1 %in% t2_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.9(sample = "RA2",gene=input$GeneRA2_1, t2.1, t2.2, t2.3, t2.4, t2.5, t2.6, t2.7,  s2.1,s2.2,s2.3,s2.4,s2.5,s2.6,s2.7, x = 65, y = 50, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA2_",input$Gene_1,".png", sep=""))
    rglwidget()})
  output$visData.6 <- renderRglwidget({
    validate(
      need(!input$GeneRA2_2 %in% t2_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.9(sample = "RA2",gene=input$GeneRA2_2,t2.1, t2.2, t2.3, t2.4, t2.5, t2.6, t2.7,  s2.1,s2.2,s2.3,s2.4,s2.5,s2.6,s2.7, x = 65, y = 50, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA2_",input$Gene_2,".png", sep=""))
    rglwidget()})
  output$visData.7 <- renderRglwidget({
    validate(
      need(!input$GeneRA2_3 %in% t2_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.9(sample = "RA2",gene=input$GeneRA2_3,t2.1, t2.2, t2.3, t2.4, t2.5, t2.6, t2.7,  s2.1,s2.2,s2.3,s2.4,s2.5,s2.6,s2.7, x = 65, y = 50, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA2_",input$Gene_3,".png", sep=""))
    rglwidget()})
  output$visData.8 <- renderRglwidget({
    validate(
      need(!input$GeneRA2_4 %in% t2_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.9(sample = "RA2",gene=input$GeneRA2_4, t2.1, t2.2, t2.3, t2.4, t2.5, t2.6, t2.7,  s2.1,s2.2,s2.3,s2.4,s2.5,s2.6,s2.7, x = 65, y = 50, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA2_",input$Gene_4,".png", sep=""))
    rglwidget()})
  session$onSessionEnded(function() { while (rgl.cur() > 0) { rgl.close() }})
}

shinyApp(ui, server)

