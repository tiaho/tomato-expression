# server.R

library(ggplot2)
seedlings <- read.csv("data/seedling_fitted_20Jun11.csv",as.is=T) #character mathcing works better if we don't convert to factors
adjusted <- read.csv("data/adjusted.csv",as.is=T)
WMadjusted <- read.csv("data/WMadjusted.csv",as.is=T)

names(seedlings)[names(seedlings) == 'X'] <- 'gene'
names(adjusted)[names(adjusted) == 'X'] <- 'gene'
names(WMadjusted)[names(WMadjusted) == 'X'] <- 'gene'

seedlings$gene.shortID <- sub("(\\.[0-9]+)+$", "",seedlings$gene)

shinyServer(function(input, output) {
  
  # plots the graph
  output$graph <- renderPlot({
    genes.input <- sub(" ","",input$gene) #strip white spaces
    genes.input <- unlist(strsplit(genes.input,split=","))
    genes.input <- sub("(\\.[0-9]+)+$", "", genes.input)

    names <- seedlings$gene[match(genes.input, seedlings$gene.shortID)]

    if (any(is.na(names)) ){
      missing <- genes.input[!genes.input %in% seedlings$gene.shortID]
      stop(paste(missing, "does not exist"))
    } else {
      # populates a group of vectors with the needed data for the genes found
      combined <- vector() 
      for (genes in names){
        place_seedlings <- match(genes, seedlings$gene, nomatch = 0)
        Species <- names(seedlings)[2:5]
        CPM <- as.numeric(seedlings[place_seedlings, 2:5])
        gene_name <- rep(genes, 4)
        if (input$logscale == TRUE) {CPM <- log2(CPM)}
        tmp <- cbind(Species, CPM, gene_name)
        combined <- rbind(combined, tmp)
      }
      df <- as.data.frame(combined)      
      
      # makes a string of the gene names
      gene_names_string <- character()
      for (genes in names){
        gene_names_string <- paste(gene_names_string, genes)
      }

      # plots
      ggplot(df, aes(x = Species, y = as.numeric(as.character(CPM)))) + 
        geom_bar(stat = "identity", position = "identity", aes(fill = Species)) + 
        facet_wrap(~ gene_name, ncol = 4) + 
        ggtitle(paste("Expression level of", gene_names_string)) +
        ylab("CPM")
      
    }
  })
  output$graph <- renderPlot({
    place_seedlings <- match(input$gene, seedlings$gene, nomatch = 0)
    if (place_seedlings == 0){
      stop(paste(input$gene, "does not exist"))
    }else{
      Species <- names(seedlings)[2:5]
      CPM <- as.numeric(seedlings[place_seedlings, 2:5])
      if (input$logscale == TRUE) {CPM <- log2(CPM)}
      combined <- cbind(Species, CPM)
      df <- as.data.frame(combined)      
      ggplot(df, aes(x = Species, y = CPM)) + 
        geom_bar(stat = "identity", position = "identity", aes(fill = Species)) + 
        ggtitle(paste("Expression level of", input$gene))
    }
  })
  
  # produces the data for 1 table, use radio buttons to determine which to plot
  table_data <- reactive({
    if (input$table_options == 1){ # Normalized CPM
      place_seedlings <- match(input$gene, seedlings$gene, nomatch = 0)
      if (place_seedlings == 0){
        stop(paste(input$gene, "does not exist"))
      }else{
        data <- data.frame()
        data[1,1] <- seedlings[place_seedlings, 1]
        data[1,2] <- seedlings[place_seedlings, 4]
        data[1,3] <- seedlings[place_seedlings, 2]
        data[1,4] <- seedlings[place_seedlings, 5]
        data[1,5] <- seedlings[place_seedlings, 3]
        colnames(data) <- c("gene", "SHA", "SLY", "SPE", "SPI")
        data
      }
    }else if (input$table_options == 2) { # log2(Normalized CPM
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
    }else if (input$table_options == 3) { # FDR Corrected p-values for Overall Significance
      place_WMadjusted <- match(input$gene, WMadjusted$gene, nomatch = 0)
      if (place_WMadjusted == 0){
        stop(paste(input$gene, "does not exist"))
      }else{
        data <- WMadjusted[place_WMadjusted, 1:2]
        colnames(data) <- c("gene", "spe")
        data
      }
    }else if (input$table_options == 4) { # FDR Corrected p-values for Pairwise Significance
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
  
  # creates the table
  output$table <- renderTable({
    table_data()
  }, digits = 4, include.rownames = FALSE)
  
  # determines if there is overall significance or not
  output$overall_significance <- renderText({
    place_WMadjusted <- match(input$gene, WMadjusted$gene, nomatch = 0)
    if (place_WMadjusted == 0) {
      stop(paste(input$gene, "does not exist"))
    }else{
      if (WMadjusted[place_WMadjusted, 2] <= 0.05) {paste("There are differences across the species for gene", input$gene)}
      else {paste("There are no differences across the species for gene", input$gene)}
    }
  })
  
  # determines if there is pairwise significance or not
  output$pairwise_significance <- renderText({
    places <- grep(input$gene, seedlings$gene)
    names <- vector()
    for (i in 1:length(places)){
      names[i] <- as.character(seedlings$gene[places[i]])
    }
    names
    place_adjusted <- match(input$gene, adjusted$gene, nomatch = 0)
    if (place_adjusted == 0){
      stop(paste("Pairwise species comparison data does not exist for", input$gene))
    }else{
      sig <- as.numeric(adjusted[place_adjusted, 2:7]) <= 0.05
      if (all(sig == rep(FALSE, 6))) {
        "There are no significant pairwise comparisons."
      }else {
        sig_places <- which(sig)
        pairs <- vector()
        for (i in 1:length(sig_places)){
          colnum <- sig_places[i] + 1
          if (i == 1){
            switch(names(adjusted)[colnum],
                   SLY_SPI = {pairs <- "S. lycopersicum and S. pimpinellifolium"},
                   SLY_SHA = {pairs <- "S. lycopersicum and S. habrochaites"},
                   SLY_SPE = {pairs <- "S. lycopersicum and S. pennellii"},
                   SPI_SHA = {pairs <- "S. pimpinellifolium and S. habrochaites"},
                   SPI_SPE = {pairs <- "S. pimpinellifolium and S. pennellii"},
                   SHA_SPE = {pairs <- "S. habrochaites and S. pennellii"}
            )}else {
              switch(names(adjusted)[colnum],
                     SLY_SPI = {pairs <- paste(pairs, "S. lycopersicum and S. pimpinellifolium", sep = "; ")},
                     SLY_SHA = {pairs <- paste(pairs, "S. lycopersicum and S. habrochaites", sep = "; ")},
                     SLY_SPE = {pairs <- paste(pairs, "S. lycopersicum and S. pennellii", sep = "; ")},
                     SPI_SHA = {pairs <- paste(pairs, "S. pimpinellifolium and S. habrochaites", sep = "; ")},
                     SPI_SPE = {pairs <- paste(pairs, "S. pimpinellifolium and S. pennellii", sep = "; ")},
                     SHA_SPE = {pairs <- paste(pairs, "S. habrochaites and S. pennellii", sep = "; ")}
              )
            }
        }
        paste("There are significant pairwise comparisons for:", pairs)
      }
    }
  })
  
  output$download_table <- downloadHandler(
    filename = function() { paste(input$gene, '.csv', sep='') },
    content = function(file) {
      write.csv(table_data(), file)
    }
  )
})