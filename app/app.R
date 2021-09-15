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
        menuItem("RA3", tabName = "RA3", icon = icon("cube")),
        menuItem("RA4", tabName = "RA4", icon = icon("cube")),
        menuItem("RA5", tabName = "RA5", icon = icon("cube")),
        menuItem("RA6", tabName = "RA6", icon = icon("cube")),
        mainPanel(selectInput('min','Min value for Log(normalized expression)', choices = c(-2:12), selected = "0"),
                  selectInput('max','Max value for Log(normalized expression)', choices = c(12:-2), selected = "6"))
      ),
      uiOutput("conds")
    )
  )

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "RA1", fluidRow(column(1),
                                      column(8, box(selectizeInput("Gene_1", label = "Choose 1st gene to display", multiple = F, choices = NULL),
                                                   rglwidgetOutput("visData.1", height = 400, width = "90%")),
                                             fluidRow(box(selectizeInput("Gene_2",label = "Choose 2nd gene to display", multiple = F, choices = NULL),
                                                          rglwidgetOutput("visData.2", height = 400, width = "90%"))))),
            fluidRow(column(1), column(8, (box(selectizeInput("Gene_3",label = "Choose 3rd gene to display", multiple = F, choices = NULL),
                                             rglwidgetOutput("visData.3", height = 400, width = "90%"))),
                                      fluidRow(box(selectizeInput("Gene_4",label = "Choose 4th gene to display", multiple = F, choices = NULL),
                                                   rglwidgetOutput("visData.4", height = 400, width = "90%")))))),
    tabItem(tabName = "RA2", fluidRow(column(1),
                                     column(8,box(selectizeInput("GeneRA2_1",label = "Choose 1st gene to display", multiple = F, choices = NULL),
                                                  rglwidgetOutput("visData.5", height = 400, width = "90%")),
                                            fluidRow(box(selectizeInput("GeneRA2_2",label = "Choose 2nd gene to display", multiple = F, choices = NULL),
                                                         rglwidgetOutput("visData.6", height = 400, width = "90%"))))),
           fluidRow(column(1),
                    column(8,(box(selectizeInput("GeneRA2_3",label = "Choose 3rd gene to display", multiple = F, choices = NULL),
                                  rglwidgetOutput("visData.7", height = 400, width = "90%"))),
                           fluidRow(box(selectizeInput("GeneRA2_4",label = "Choose 4th gene to display", multiple = F, choices = NULL),
                                        rglwidgetOutput("visData.8", height = 400, width = "90%")))))),
    tabItem(tabName = "RA3", fluidRow(column(1),
                                     column(8,box(selectizeInput("GeneRA3_1",label = "Choose 1st gene to display", multiple = F, choices = NULL),
                                                  rglwidgetOutput("visData.9", height = 400, width = "90%")),
                                            fluidRow(box(selectizeInput("GeneRA3_2",label = "Choose 2nd gene to display", multiple = F, choices = NULL),
                                                         rglwidgetOutput("visData.10", height = 400, width = "90%"))))),
           fluidRow(column(1),column(8,(box(selectizeInput("GeneRA3_3",label = "Choose 3rd gene to display", multiple = F, choices = NULL),
                                            rglwidgetOutput("visData.11", height = 400, width = "90%"))),
                                     fluidRow(box(selectizeInput("GeneRA3_4",label = "Choose 4th gene to display", multiple = F, choices = NULL),
                                                  rglwidgetOutput("visData.12", height = 400, width = "90%")))))),
    tabItem(tabName = "RA4", fluidRow(column(1),
                                     column(8,box(selectizeInput("GeneRA4_1",label = "Choose 1st gene to display", multiple = F, choices = NULL),
                                                  rglwidgetOutput("visData.13", height = 400, width = "90%")),
                                            fluidRow(box(selectizeInput("GeneRA4_2",label = "Choose 2nd gene to display", multiple = F, choices = NULL),
                                                         rglwidgetOutput("visData.14", height = 400, width = "90%"))))),
           fluidRow(column(1),column(8,(box(selectizeInput("GeneRA4_3",label = "Choose 3rd gene to display", multiple = F, choices = NULL),
                                            rglwidgetOutput("visData.15", height = 400, width = "90%"))),
                                     fluidRow(box(selectizeInput("GeneRA4_4",label = "Choose 4th gene to display", multiple = F, choices = NULL),
                                                  rglwidgetOutput("visData.16", height = 400, width = "90%")))))),
    tabItem(tabName = "RA5", fluidRow(column(1),
                                     column(8,box(selectizeInput("GeneRA5_1",label = "Choose 1st gene to display", multiple = F, choices = NULL),
                                                  rglwidgetOutput("visData.17", height = 400, width = "90%")),
                                            fluidRow(box(selectizeInput("GeneRA5_2",label = "Choose 2nd gene to display", multiple = F, choices = NULL),
                                                         rglwidgetOutput("visData.18", height = 400, width = "90%"))))),
           fluidRow(column(1),column(8,(box(selectizeInput("GeneRA5_3",label = "Choose 3rd gene to display", multiple = F, choices = NULL),
                                            rglwidgetOutput("visData.19", height = 400, width = "90%"))),
                                     fluidRow(box(selectizeInput("GeneRA5_4",label = "Choose 4th gene to display", multiple = F, choices = NULL),
                                                  rglwidgetOutput("visData.20", height = 400, width = "90%")))))),
    tabItem(tabName = "RA6", fluidRow(column(1),
                                      column(8,box(selectizeInput("GeneRA6_1",label = "Choose 1st gene to display", multiple = F, choices = NULL),
                                                   rglwidgetOutput("visData.21", height = 400, width = "90%")),
                                             fluidRow(box(selectizeInput("GeneRA6_2",label = "Choose 2nd gene to display", multiple = F, choices = NULL),
                                                          rglwidgetOutput("visData.22", height = 400, width = "90%"))))),
            fluidRow(column(1),column(8,(box(selectizeInput("GeneRA6_3",label = "Choose 3rd gene to display", multiple = F, choices = NULL),
                                             rglwidgetOutput("visData.23", height = 400, width = "90%"))),
                                      fluidRow(box(selectizeInput("GeneRA6_4",label = "Choose 4th gene to display", multiple = F, choices = NULL),
                                                   rglwidgetOutput("visData.24", height = 400, width = "90%"))))))
  )
)

# Put them together into a dashboardPage
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "3DSeTH: Spatial Transcriptomics for Heatmap gene expression visualization",titleWidth = 850),
  sidebar,
  body
)

server <- function(input, output, session){
  
  # Server side selectize for RA1
  updateSelectizeInput(session, "Gene_1", choices = unique(t1_names, t2_names), selected = "CCL19", server = TRUE)
  updateSelectizeInput(session, "Gene_2", choices = unique(t1_names, t2_names), selected = "CD52", server = TRUE)
  updateSelectizeInput(session, "Gene_3", choices = unique(t1_names, t2_names), selected = "MS4A1", server = TRUE)
  updateSelectizeInput(session, "Gene_4", choices = unique(t1_names, t2_names), selected = "SSR4", server = TRUE)
  
  # Server side selectize for RA2
  updateSelectizeInput(session, "GeneRA2_1", choices = unique(t1_names, t2_names), selected = "CCL19", server = TRUE)
  updateSelectizeInput(session, "GeneRA2_2", choices = unique(t1_names, t2_names), selected = "CD52", server = TRUE)
  updateSelectizeInput(session, "GeneRA2_3", choices = unique(t1_names, t2_names), selected = "MS4A1", server = TRUE)
  updateSelectizeInput(session, "GeneRA2_4", choices = unique(t1_names, t2_names), selected = "SSR4", server = TRUE)
  
  # Server side selectize for RA3
  updateSelectizeInput(session, "GeneRA3_1", choices = unique(t1_names, t2_names), selected = "CCL19", server = TRUE)
  updateSelectizeInput(session, "GeneRA3_2", choices = unique(t1_names, t2_names), selected = "CD52", server = TRUE)
  updateSelectizeInput(session, "GeneRA3_3", choices = unique(t1_names, t2_names), selected = "MS4A1", server = TRUE)
  updateSelectizeInput(session, "GeneRA3_4", choices = unique(t1_names, t2_names), selected = "SSR4", server = TRUE)
  
  # Server side selectize for RA4
  updateSelectizeInput(session, "GeneRA4_1", choices = unique(t1_names, t2_names), selected = "CCL19", server = TRUE)
  updateSelectizeInput(session, "GeneRA4_2", choices = unique(t1_names, t2_names), selected = "CD52", server = TRUE)
  updateSelectizeInput(session, "GeneRA4_3", choices = unique(t1_names, t2_names), selected = "MS4A1", server = TRUE)
  updateSelectizeInput(session, "GeneRA4_4", choices = unique(t1_names, t2_names), selected = "SSR4", server = TRUE)
  
  # Server side selectize for RA5
  updateSelectizeInput(session, "GeneRA5_1", choices = unique(t1_names, t2_names), selected = "CCL19", server = TRUE)
  updateSelectizeInput(session, "GeneRA5_2", choices = unique(t1_names, t2_names), selected = "CD52", server = TRUE)
  updateSelectizeInput(session, "GeneRA5_3", choices = unique(t1_names, t2_names), selected = "MS4A1", server = TRUE)
  updateSelectizeInput(session, "GeneRA5_4", choices = unique(t1_names, t2_names), selected = "SSR4", server = TRUE)
  
  # Server side selectize for RA6
  updateSelectizeInput(session, "GeneRA6_1", choices = unique(t1_names, t2_names), selected = "CCL19", server = TRUE)
  updateSelectizeInput(session, "GeneRA6_2", choices = unique(t1_names, t2_names), selected = "CD52", server = TRUE)
  updateSelectizeInput(session, "GeneRA6_3", choices = unique(t1_names, t2_names), selected = "MS4A1", server = TRUE)
  updateSelectizeInput(session, "GeneRA6_4", choices = unique(t1_names, t2_names), selected = "SSR4", server = TRUE)
  
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
  output$visData.9 <- renderRglwidget({
    validate(
      need(!input$GeneRA3_1 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.4(sample = "RA3", gene=input$GeneRA3_1, t3.1, t3.2, t3.3, t3.4, s3.1,s3.2,s3.3,s3.4, x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA3_",input$Gene_1,".png", sep=""))
    rglwidget()})
  output$visData.10 <- renderRglwidget({
    validate(
      need(!input$GeneRA3_2 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.4(sample = "RA3",gene=input$GeneRA3_2,t3.1, t3.2, t3.3, t3.4, s3.1,s3.2,s3.3,s3.4,  x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA3_",input$Gene_2,".png", sep=""))
    rglwidget()})
  output$visData.11 <- renderRglwidget({
    validate(
      need(!input$GeneRA3_3 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.4(sample = "RA3",gene=input$GeneRA3_3,t3.1, t3.2, t3.3, t3.4, s3.1,s3.2,s3.3,s3.4,  x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA3_",input$Gene_3,".png", sep=""))
    rglwidget()})
  output$visData.12 <- renderRglwidget({
    validate(
      need(!input$GeneRA3_4 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.4(sample = "RA3",gene=input$GeneRA3_4,t3.1, t3.2, t3.3, t3.4, s3.1,s3.2,s3.3,s3.4,  x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA3_",input$Gene_4,".png", sep=""))
    rglwidget()})
  output$visData.13 <- renderRglwidget({
    validate(
      need(!input$GeneRA4_1 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.5(sample = "RA4", gene=input$GeneRA4_1, t4.1, t4.2, t4.3, t4.4, t4.5, s4.1,s4.2,s4.3,s4.4, s4.5, x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA4_",input$Gene_1,".png", sep=""))
    rglwidget()})
  output$visData.14 <- renderRglwidget({
    validate(
      need(!input$GeneRA4_2 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.5(sample = "RA4",gene=input$GeneRA4_2,t4.1, t4.2, t4.3, t4.4, t4.5, s4.1,s4.2,s4.3,s4.4, s4.5, x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA4_",input$Gene_2,".png", sep=""))
    rglwidget()})
  output$visData.15 <- renderRglwidget({
    validate(
      need(!input$GeneRA4_3 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.5(sample = "RA4",gene=input$GeneRA4_3,t4.1, t4.2, t4.3, t4.4, t4.5, s4.1,s4.2,s4.3,s4.4, s4.5,  x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA4_",input$Gene_3,".png", sep=""))
    rglwidget()})
  output$visData.16 <- renderRglwidget({
    validate(
      need(!input$GeneRA4_4 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.5(sample = "RA4",gene=input$GeneRA4_4,t4.1, t4.2, t4.3, t4.4, t4.5, s4.1,s4.2,s4.3,s4.4, s4.5, x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA4_",input$Gene_4,".png", sep=""))
    rglwidget()})
  output$visData.17 <- renderRglwidget({
    validate(
      need(!input$GeneRA5_1 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.3(sample = "RA5", gene=input$GeneRA5_1, t5.1, t5.2, t5.3, s5.1, s5.2, s5.3, x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA5_",input$Gene_1,".png", sep=""))
    rglwidget()})
  output$visData.18 <- renderRglwidget({
    validate(
      need(!input$GeneRA5_2 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.3(sample = "RA5",gene=input$GeneRA5_2,t5.1, t5.2, t5.3, s5.1, s5.2, s5.3,  x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA5_",input$Gene_2,".png", sep=""))
    rglwidget()})
  output$visData.19 <- renderRglwidget({
    validate(
      need(!input$GeneRA5_3 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.3(sample = "RA5",gene=input$GeneRA5_3,t5.1, t5.2, t5.3, s5.1, s5.2, s5.3,   x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA5_",input$Gene_3,".png", sep=""))
    rglwidget()})
  output$visData.20 <- renderRglwidget({
    validate(
      need(!input$GeneRA5_4 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.3(sample = "RA5",gene=input$GeneRA5_4,t5.1, t5.2, t5.3, s5.1,s5.2,s5.3,  x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA5_",input$Gene_4,".png", sep=""))
    rglwidget()})
  output$visData.21 <- renderRglwidget({
    validate(
      need(!input$GeneRA6_1 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.4(sample = "RA6",gene=input$GeneRA6_1,t6.1, t6.2, t6.3, t6.4, s6.1,s6.2,s6.3,s6.4,  x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA4_",input$Gene_3,".png", sep=""))
    rglwidget()})
  output$visData.22 <- renderRglwidget({
    validate(
      need(!input$GeneRA6_2 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.4(sample = "RA6",gene=input$GeneRA6_2,t6.1, t6.2, t6.3, t6.4, s6.1,s6.2,s6.3,s6.4,  x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA4_",input$Gene_3,".png", sep=""))
    rglwidget()})
  output$visData.23 <- renderRglwidget({
    validate(
      need(!input$GeneRA6_3 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.4(sample = "RA6",gene=input$GeneRA6_3,t6.1, t6.2, t6.3, t6.4, s6.1,s6.2,s6.3,s6.4,  x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA4_",input$Gene_3,".png", sep=""))
    rglwidget()})
  output$visData.24 <- renderRglwidget({
    validate(
      need(!input$GeneRA6_4 %in% t1_names == FALSE, "The selected gene is not expressed in this dataset.")
    )
    plot.gene.3d.4(sample = "RA6",gene=input$GeneRA6_4,t6.1, t6.2, t6.3, t6.4, s6.1,s6.2,s6.3,s6.4,  x = 59, y = 40, transparency = 0.1, min = as.numeric(input$min), max = as.numeric(input$max))
    #snapshot3d(filename = paste("tmp/RA4_",input$Gene_3,".png", sep=""))
    rglwidget()})
  session$onSessionEnded(function() { while (rgl.cur() > 0) { rgl.close() }})
}

shinyApp(ui, server)

