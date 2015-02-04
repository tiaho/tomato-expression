# ui.R

library(shiny)

shinyUI(fluidPage(
  titlePanel("Tomato Expression"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("gene", label = h5("Enter desired gene name(s) separated by commas"), value = "Solyc02g081130.1.1"),
      br(),
      br(),
      checkboxInput("logscale", "Plots in a log scale", value = TRUE),
      br(),
      radioButtons("table_options", label = h5("Display values in a table"),
                   choices = list("Normalized CPM" = 1,
                                  "log2(Normalized CPM)" = 2,
                                  "FDR Corrected p-values for Overall Significance" = 3,
                                  "FDR Corrected p-values for Pairwise Significance" = 4),
                   selected = 2),
      br(),
      downloadButton('download_table', 'Download Table'),
      br(),
      br(),
      textOutput("overall_significance"),
      br(),
      textOutput("pairwise_significance")
    ),
    
    mainPanel(
      plotOutput("graph"),
      br(),
      br(),
      tableOutput("table")
    )
  )
    
))