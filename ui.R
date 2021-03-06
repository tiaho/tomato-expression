# ui.R

library(shiny)

shinyUI(fluidPage(
  titlePanel("Tomato Expression"),

  sidebarLayout(
    sidebarPanel(
      selectizeInput("gene",
                     label = h5("Enter desired gene name(s)"),
                     multiple = TRUE,
                     selected = "Solyc02g081130.1.1",
                     choices = c("Solyc02g081130.1.1")
      ),
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
      downloadButton('download_table_fdr', 'Download the FDR Table')
    ),

    mainPanel(
      plotOutput("graph", width = "100%", height = "auto"),
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
