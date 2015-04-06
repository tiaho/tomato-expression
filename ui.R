# ui.R

library(shiny)

seedlings <- read.csv("data/seedling_fitted_20Jun11.csv",as.is=T) #character mathcing works better if we don't convert to factors
adjusted <- read.csv("data/adjusted.csv",as.is=T)
WMadjusted <- read.csv("data/WMadjusted.csv",as.is=T)

gene_list <- c(seedlings$X, adjusted$X, WMadjusted$X)
gene_list <- unique(gene_list)

shinyUI(fluidPage(
  titlePanel("Tomato Expression"),
  
  sidebarLayout(
    sidebarPanel(
      selectizeInput("gene", label = h5("Enter desired gene name(s) separated by commas"),
                     multiple = TRUE, selected = "Solyc02g081130.1.1",
                     choices = gene_list, 
      ),
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