library(shiny)
library(ggplot2)
library(dplyr)
library(bslib)
library(Rtsne)
library(tibble)
library(tidyr)

# Load aggregated dataset 
aggregated_data <- read.csv("aggregated_data.csv")
colnames(aggregated_data) <- c("Mouse_Strain", "Phenotype_Group", "Phenotype", "Gene", "Aggregated_PValue")

# Define UI
ui <- page_navbar(
  title = "Genotype-Phenotype Viewer: Knockout Mice Analysis",
  bg = "#2D89C8",
  inverse = TRUE,
  tags$head(
    tags$style(HTML("
      .navbar-nav {
        margin-left: auto;  /* Moves the navbar tabs to the right */
      }
    "))
  ),
  # First Tab: Phenotype Visualization
  nav_panel(
    title = "Phenotype Visualization",
    sidebarLayout(
      sidebarPanel(
        h4("Filters"),
        # Dropdowns for filtering data
        selectInput("mouse_strain", "Select Mouse Strain:", 
                    choices = c("All", unique(aggregated_data$Mouse_Strain)), 
                    selected = unique(aggregated_data$Mouse_Strain)[1]),
        selectInput("phenotype_group", "Select Phenotype Group:",
                    choices = c("All", unique(aggregated_data$Phenotype_Group)),
                    selected = "All"),
        selectInput("phenotype", "Select Phenotype:",
                    choices = c("All"),
                    selected = "All"),
        selectInput("gene", "Select Gene:", 
                    choices = c("All", unique(aggregated_data$Gene)), 
                    selected = "All"),
        hr(),
        # Checkbox to filter by significance
        checkboxInput("significant_only", "Show Only Significant Results (P < 0.05)", TRUE),
        hr(),
        # Chart type selection
        radioButtons("chart_type", "Select Chart Type:",
                     choices = c("Scatter Plot" = "scatter", 
                                 "Bar Plot" = "bar")),
        # Export buttons
        downloadButton("download_plot", "Download Plot"),
        downloadButton("download_table", "Download Table")
      ),
      mainPanel(
        # Main display: Data plot and table outputs
        tabsetPanel(
          tabPanel("Data Plot", 
                   plotOutput("main_plot", height = "700px")),
          tabPanel("Data Table", 
                   tableOutput("filtered_table")),
        )
      )
    )
  ),
  # Second Tab: Gene Clustering
  nav_panel(
    title = "Gene Clustering",
    sidebarLayout(
      sidebarPanel(
        h4("Filters"),
        
        # Select mouse strain
        selectInput("mouse_strain_clustering", "Select Mouse Strain:", 
                    choices = c(unique(aggregated_data$Mouse_Strain)), 
                    selected = unique(aggregated_data$Mouse_Strain)[1]),
        
        # Selection between Phenotype_Group and Phenotype
        radioButtons("selection_type_clustering", "Choose Filter Type:",
                     choices = c("Phenotype Group" = "group", "Phenotype" = "phenotype"),
                     selected = "group"),
        
        # Conditional UI for Phenotype Group
        conditionalPanel(
          condition = "input.selection_type_clustering == 'group'",
          selectInput("phenotype_group_clustering", "Select Phenotype Group:",
                      choices = c(unique(aggregated_data$Phenotype_Group)),
                      selected = unique(aggregated_data$Phenotype_Group)[1])
        ),
        
        # Conditional UI for Phenotype
        conditionalPanel(
          condition = "input.selection_type_clustering == 'phenotype'",
          selectInput("phenotype_clustering", "Select Phenotype:",
                      choices = c(unique(aggregated_data$Phenotype)),
                      selected = unique(aggregated_data$Phenotype)[1])
      ),
      # Export button for PCA and t-SNE plots
      downloadButton("download_pca_plot", "Download PCA"),
      downloadButton("download_tsne_plot", "Download t-SNE")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("PCA Plot", plotOutput("pca_plot", height = "700px")),
          tabPanel("t-SNE Plot", plotOutput("tsne_plot", height = "700px"))),
        textOutput("errorMessage_clustering")
      )
    )
  )
)


# Define Server
server <- function(input, output, session) {
  
  # Update Phenotype Dropdown based on Group Selection
  observeEvent(input$phenotype_group, {
    phenotypes_in_group <- if (input$phenotype_group == "All") {
      unique(aggregated_data$Phenotype)
    } else {
      aggregated_data %>%
        filter(Phenotype_Group == input$phenotype_group) %>%
        pull(Phenotype) %>%
        unique()
    }
    
    # Update phenotype choices
    updateSelectInput(session, "phenotype",
                      choices = c("All", phenotypes_in_group),
                      selected = "All")
  })
  
  # Reactive dataset for main filtering
  filtered_data <- reactive({
    data <- aggregated_data
    
    if (input$mouse_strain != "All") {
      data <- data %>% filter(Mouse_Strain == input$mouse_strain)
    }
    if (input$phenotype_group != "All") {
      data <- data %>% filter(Phenotype_Group == input$phenotype_group)
    }
    if (input$phenotype != "All") {
      data <- data %>% filter(Phenotype == input$phenotype)
    }
    if (input$gene != "All") {
      data <- data %>% filter(Gene == input$gene)
    }
    if (input$significant_only) {
      data <- data %>% filter(Aggregated_PValue < 0.05)
    }
    
    data
  })
  
  # Generate Main Plot
  
  # Define the plot reactively
  plot_data <- reactive({
    data <- filtered_data()
    
    if (nrow(data) == 0) {
      validate(need(FALSE, "No data available based on the current filter selections. Please adjust your filter criteria."))
    }
    
    data <- data %>%
      mutate(Significance = ifelse(Aggregated_PValue < 0.05, "Significant", "Not Significant"))
    
    # Create the plot based on the chart type selected
    p <- NULL
    if (input$chart_type == "scatter") {
      p <- ggplot(data, aes(x = Phenotype, y = Aggregated_PValue, color = Significance)) +
        geom_point() +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue")) +
        labs(title = paste("Mouse Strain:", input$mouse_strain, " - P-Values by Phenotype"),
             x = "Phenotype", y = "Aggregated P-Value",
             color = "Significance")
    } else if (input$chart_type == "bar") {
      p <- ggplot(data, aes(x = Phenotype, y = Aggregated_PValue, fill = Gene)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Bar Plot of Phenotypes", x = "Phenotype", y = "Aggregated P-Value")
    }
    
    p  # Return the plot object
  })
  
  # Render the plot in the UI
  output$main_plot <- renderPlot({
    plot_data()  # Call the reactive expression that generates the plot
  })
  
  # Download handler for the plot
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Get the plot object from the reactive expression
      p <- plot_data()
      
      # Save the plot to the specified file
      ggsave(file, plot = p, width = 8, height = 6)
    }
  )
  
  
  # Data table output
  output$filtered_table <- renderTable({
    data <- filtered_data()
    
    if (nrow(data) == 0) {
      return("No data available based on the current filter selections. Please adjust your filter criteria.")
    }
    
    data
  })
  
  # Reactive dataset for table download
  output$download_table <- downloadHandler(
    filename = function() {
      paste("filtered_data_", Sys.Date(), ".csv", sep = "")  # Create a filename with the current date
    },
    content = function(file) {
      # Get the filtered data
      data <- filtered_data()
      
      # Write the filtered data to the CSV file
      write.csv(data, file, row.names = FALSE)
    }
  )
  
  
  # Filter data for gene clustering
  filtered_data_pca <- reactive({
    # Start with filtering by Mouse Strain
    data <- aggregated_data[aggregated_data$Mouse_Strain == input$mouse_strain_clustering,]
    
    # Apply additional filtering based on selection type
    if (input$selection_type_clustering == "group") {
      data <- data[data$Phenotype_Group == input$phenotype_group_clustering,]
    } else if (input$selection_type_clustering == "phenotype") {
      data <- data[data$Phenotype == input$phenotype_clustering,]
    }
    
    return(data)
  })
  
  # Reactive PCA plot definition
  pca_plot_reactive <- reactive({
    tryCatch({
      data <- filtered_data_pca()  # Assuming filtered_data_pca() gives the required data
      
      # Prepare data for PCA dynamically based on selection type
      pivot_column <- if (input$selection_type_clustering == "group") "Phenotype_Group" else "Phenotype"
      
      # Prepare data for PCA
      data <- data %>%
        mutate(transformed_pvalue = -log10(Aggregated_PValue)) %>%
        pivot_wider(names_from = all_of(pivot_column), values_from = transformed_pvalue) %>%
        {
          if (input$selection_type_clustering == "group") {
            mutate(., Gene = paste0(Gene, "_", row_number()))  # Add unique row number for group selection
          } else {
            .  # Return the dataset unchanged
          }
        } %>%
        column_to_rownames("Gene")  # Name the rownames with Gene
      
      # Extract numeric columns for PCA
      numeric_data <- data %>%
        select(where(is.numeric))
      
      # Impute missing values with column means and remove zero-variance columns
      numeric_data <- numeric_data %>%
        mutate(across(
          everything(),
          ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)  # Impute missing values
        )) %>%
        select(where(~ !any(is.na(.)) && sd(.x, na.rm = TRUE) > 0))  # Remove columns with NA or zero variance
      
      # Standardize data
      standardized_data <- scale(numeric_data)
      
      # K-means clustering
      set.seed(123)  # For reproducibility
      kmeans_result <- kmeans(standardized_data, centers = 3)  # Choose appropriate k
      
      # Add cluster information
      numeric_data$cluster <- as.factor(kmeans_result$cluster)
      
      # PCA for visualization with error handling
      pca_result <- prcomp(standardized_data, center = TRUE, scale. = TRUE)
      
      # Convert PCA results to data frame
      pca_data <- as.data.frame(pca_result$x) %>%
        cbind(cluster = numeric_data$cluster)
      
      # Plot
      pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
        geom_point() +
        geom_text(aes(label = rownames(pca_data), vjust = -0.5, size = 3)) +
        ggtitle("PCA of Genes with Similar Phenotype Score") +
        theme_minimal()
      
      return(pca_plot)
      
    }, error = function(e) {
      output$errorMessage_clustering <- renderText({
        print("Error: The data is not suitable or sufficient for PCA")
      })
      return(NULL)
    })
  })
  
  # Render the PCA plot in the UI
  output$pca_plot <- renderPlot({
    pca_plot_reactive()  # Call the reactive function to generate the plot
  })
  
  # Download handler for PCA plot
  output$download_pca_plot <- downloadHandler(
    filename = function() {
      paste("pca_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Get the PCA plot from the reactive expression
      pca_plot <- pca_plot_reactive()
      
      # Save the plot to the file
      ggsave(file, plot = pca_plot, width = 8, height = 6, units = "in")
    }
  )
  # Reactive t-SNE plot definition
  tsne_plot_reactive <- reactive({
    tryCatch({
      data <- filtered_data_pca()  # Assuming filtered_data_pca() gives the required data
      
      # Prepare data for t-SNE dynamically based on selection type
      pivot_column <- if (input$selection_type_clustering == "group") "Phenotype_Group" else "Phenotype"
      
      # Prepare data for t-SNE
      data <- data %>%
        mutate(transformed_pvalue = -log10(Aggregated_PValue)) %>%
        pivot_wider(names_from = all_of(pivot_column), values_from = transformed_pvalue) %>%
        {
          if (input$selection_type_clustering == "group") {
            mutate(., Gene = paste0(Gene, "_", row_number())) # Add unique row number for group selection
          } else {
            . # Return the dataset unchanged
          }
        } %>%
        column_to_rownames("Gene") # Name the rownames with Gene
      
      # Extract numeric columns for t-SNE
      numeric_data <- data %>%
        select(where(is.numeric))
      
      # Impute missing values with column means and remove zero-variance columns
      numeric_data <- numeric_data %>%
        mutate(across(
          everything(),
          ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)  # Impute missing values
        )) %>%
        select(where(~ !any(is.na(.)) && sd(.x, na.rm = TRUE) > 0))  # Remove columns with NA or zero variance
      
      # Standardize data
      standardized_data <- scale(numeric_data)
      
      # K-means clustering
      set.seed(123)  # For reproducibility
      kmeans_result <- kmeans(standardized_data, centers = 3)  # Choose appropriate k
      
      # Add cluster information
      numeric_data$cluster <- as.factor(kmeans_result$cluster)
      
      # Perform t-SNE
      tsne_result <- Rtsne(standardized_data, dims = 2, perplexity = 2, verbose = FALSE)
      tsne_data <- as.data.frame(tsne_result$Y)
      colnames(tsne_data) <- c("Dim1", "Dim2")
      
      # Add row names as a column for labeling
      tsne_data$Gene <- rownames(data)
      tsne_data$cluster <- numeric_data$cluster
      
      # Plot
      ggplot(tsne_data, aes(x = Dim1, y = Dim2, color = cluster)) +
        geom_point() +
        geom_text(aes(label = Gene), vjust = -0.5, size = 4) +
        ggtitle("t-SNE of Genes with Similar Phenotype Score") +
        theme_minimal()
      
    }, error = function(e){
      output$errorMessage_clustering <- renderText({
        print("Error: The data is not suitable or sufficient for t-SNE")
      })
      return(NULL)
    })
  })
  
  # Render the t-SNE plot in the UI
  output$tsne_plot <- renderPlot({
    tsne_plot_reactive()  # Call the reactive function to generate the plot
  })
  # Download handler for t-SNE plot
  output$download_tsne_plot <- downloadHandler(
    filename = function() {
      paste("tsne_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Get the t-SNE plot from the reactive expression
      tsne_plot <- tsne_plot_reactive()
      
      # Save the plot to the file
      ggsave(file, plot = tsne_plot, width = 8, height = 6, units = "in")
    }
  )
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
