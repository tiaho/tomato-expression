# server.R

library(ggplot2)
seedlings <- read.csv("data/seedling_fitted_20Jun11.csv",as.is=T) #character mathcing works better if we don't convert to factors
adjusted <- read.csv("data/adjusted.csv",as.is=T)
WMadjusted <- read.csv("data/WMadjusted.csv",as.is=T)

names(seedlings)[names(seedlings) == 'X'] <- 'gene'
names(adjusted)[names(adjusted) == 'X'] <- 'gene'
names(WMadjusted)[names(WMadjusted) == 'X'] <- 'gene'

gene_list <- sort(seedlings$gene)
# seedlings$gene.shortID <- sub("(\\.[0-9]+)+$", "",seedlings$gene)

shinyServer(function(input, output, session) {

  updateSelectizeInput(session, 'gene', choices  = gene_list,
                                        server   = TRUE,
                                        selected = "Solyc02g081130.1.1")
  
  # sets the height for the plot
#   my_height <- reactive({
#     800
# #     200 + (150 * length(input$gene))
#   })
#   #   plot_height <- reactive({
#   #     300 + (150 * (length(input$gene) / 4))
#   #   })
  
  # graph
  output$graph <- renderPlot({
    data <- vector()
    for (i in 1:length(input$gene)){
      place_seedlings <- match(input$gene[i], seedlings$gene, nomatch = 0)
      if (place_seedlings == 0){
        stop(paste(input$gene, "does not exist"))
      }else{
        Species <- names(seedlings)[2:5]
        CPM <- as.numeric(seedlings[place_seedlings, 2:5])
        if (input$logscale == 2) {CPM <- log2(CPM)}
        gene <- rep(input$gene[i], 4)
        tmp_data <- cbind(Species, CPM, gene)
        data <- rbind(data, tmp_data)
      } 
    }

      df <- as.data.frame(data)
      df$CPM <- as.numeric(as.character((df$CPM)))
  
      ggplot(df, aes(x = Species, y = CPM)) + 
        geom_bar(stat = "identity", position = "identity", aes(fill = Species)) + 
        facet_wrap( ~ gene, ncol = 4) +
        ggtitle(paste("Expression Level of Selected Genes"))
  })
  
  
  
  # produces the data for the CPM table
  table_data_cpm <- reactive({
        if (input$logscale == 1){ # Normalized CPM
          data <- data.frame()
          for (i in 1:length(input$gene)){
            place_seedlings <- match(input$gene[i], seedlings$gene, nomatch = 0)
            if (place_seedlings == 0){
              stop(paste(input$gene, "does not exist"))
            }else{
                data[i,1] <- seedlings[place_seedlings, 1]
                data[i,2] <- seedlings[place_seedlings, 4]
                data[i,3] <- seedlings[place_seedlings, 2]
                data[i,4] <- seedlings[place_seedlings, 5]
                data[i,5] <- seedlings[place_seedlings, 3]
            }
          }
          colnames(data) <- c("gene", "SHA", "SLY", "SPE", "SPI")
          data
        }else if (input$logscale == 2) { # log2(Normalized CPM
          place_seedlings <- match(input$gene, seedlings$gene, nomatch = 0)
          if (place_seedlings == 0){
            stop(paste(input$gene, "does not exist"))
          }else{
            data <- seedlings[place_seedlings, 1:5]
            data[1,1] <- seedlings[place_seedlings, 1]
            data[1,2] <- seedlings[place_seedlings, 4]
            data[1,3] <- seedlings[place_seedlings, 2]
            data[1,4] <- seedlings[place_seedlings, 5]
            data[1,5] <- seedlings[place_seedlings, 3]
            data[,2:5] <- log2(data[,2:5])
            colnames(data) <- c("gene", "SHA", "SLY", "SPE", "SPI")
            data
          }
        }
  })
#   }, height = my_height)
  
  # creates the table for CPM
  output$table_cpm <- renderTable({
    table_data_cpm()
  }, digits = 4, include.rownames = FALSE)
  
  # makes a title for the CPM table
  output$title1 <- renderText({
    if (input$logscale == 1){
      type = "Normalized CPM"
    } else{
      type = "log2(Normalized CPM)"
    }
    paste(type)
  })
  
  # download the CPM table
  output$download_table_cpm <- downloadHandler(
    filename = function() {
      if (input$logscale == 1){
        type = "Normalized CPM"
      } else{
        type = "log2(Normalized CPM)"
      }
      paste(input$gene, '_CPM_', type, '.csv', sep='')
    },
    content = function(file) {
      write.csv(table_data_cpm(), file)
    })
  
  
  # produces the data for the FDR table
  table_data_fdr <- reactive({
  if (input$table_options == 1) { # FDR Corrected p-values for Overall Significance
      place_WMadjusted <- match(input$gene, WMadjusted$gene, nomatch = 0)
      if (place_WMadjusted == 0){
        stop(paste(input$gene, "does not exist"))
      }else{
        data <- WMadjusted[place_WMadjusted, 1:2]
        colnames(data) <- c("gene", "spe")
        data
      }
    }else if (input$table_options == 2) { # FDR Corrected p-values for Pairwise Significance
      place_adjusted <- match(input$gene, adjusted$gene, nomatch = 0)
      if (place_adjusted == 0){
        stop(paste("Pairwise species comparison data does not exist for", input$gene))
      }else{
        data <- adjusted[place_adjusted, 1:7]
        colnames(data) <- c("gene", "SLY_SPI", "SLY_SHA", "SLY_SPE", "SPI_SHA", "SPI_SPE", "SHA_SPE")
        data
      }
    }
  })
  
  # creates the table for FDR
  output$table_fdr <- renderTable({
    table_data_fdr()
  }, digits = 4, include.rownames = FALSE)
  
  # makes a title for the FDR table
  output$title2 <- renderText({
    if (input$table_options == 1){
      type = "Overall Significance"
    } else{
      type = "Pairwise Significance"
    }
    paste(type)
  })
  
  # download the FDR table
  output$download_table_fdr <- downloadHandler(
    filename = function() {
      if (input$table_options == 1){
        type = "Overall Significance"
      } else{
        type = "Pairwise Significance"
      }
      paste(input$gene, '_FDR_', type, '.csv', sep='')
      },
    content = function(file) {
      write.csv(table_data_fdr(), file)
    })
  
})