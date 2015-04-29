# ui.R

library(shiny)

seedlings <- read.csv("data/seedling_fitted_20Jun11.csv",as.is=T) #character macthing works better if we don't convert to factors
adjusted <- read.csv("data/adjusted.csv",as.is=T)
WMadjusted <- read.csv("data/WMadjusted.csv",as.is=T)

# gene_list <- c(seedlings$X, adjusted$X, WMadjusted$X)
# gene_list <- sort(unique(gene_list))

gene_list <- sort(seedlings$X)

shinyUI(fluidPage(
  titlePanel("Tomato Expression"),
  
  sidebarLayout(
    sidebarPanel(
      selectizeInput("gene", label = h5("Enter desired gene name(s) separated by commas"),
                     multiple = TRUE, selected = "Solyc02g081130.1.1",
                     choices = gene_list
      ),
      br(),
      br(),
      radioButtons("logscale", label = h5("Plot/Display Options"),
                   choices = list("Normalized CPM" = 1,
                                  "log2(Normalized CPM)" = 2),
                   selected = 1),
      br(),
      downloadButton('download_table_cpm', 'Download the CPM Table'),
      br(),
      br(),
      radioButtons("table_options", label = h5("FDR Options"),
                   choices = list("FDR Corrected p-values for Overall Significance" = 1,
                                  "FDR Corrected p-values for Pairwise Significance" = 2),
                   selected = 1),
      br(),
      downloadButton('download_table_fdr', 'Download the FDR Table'),
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
      textOutput("title1"),
      tableOutput("table_cpm"),
      br(),
      br(),
      textOutput("title2"),
      tableOutput("table_fdr")
    )
  )
    
))