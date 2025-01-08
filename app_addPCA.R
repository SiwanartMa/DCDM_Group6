library(shiny)
library(ggplot2)
library(dplyr)
library(bslib)
library(Rtsne)
library(tibble)
library(tidyr)

# Load aggregated dataset (replace with your actual file path)
aggregated_data <- read.csv("aggregated_data.csv")
colnames(aggregated_data) <- c("Mouse_Strain", "Phenotype_Group", "Phenotype", "Gene", "Aggregated_PValue")

# Define UI
ui <- page_navbar(
  title = "IMPC Visualizations",
  bg = "#2D89C8",
  inverse = TRUE,
  nav_panel(
    title = "Visualization",
    sidebarLayout(
      sidebarPanel(
        h4("Filters"),
        selectInput("mouse_strain", "Select Mouse Strain:", 
                    choices = c("All", unique(aggregated_data$Mouse_Strain)), 
                    selected = unique(aggregated_data$Mouse_Strain)[1]),
        selectInput("phenotype", "Select Phenotype:", 
                    choices = c("All", unique(aggregated_data$Phenotype)), 
                    selected = "All"),
        selectInput("gene", "Select Gene:", 
                    choices = c("All", unique(aggregated_data$Gene)), 
                    selected = "All"),
        hr(),
        checkboxInput("significant_only", "Show Only Significant Results (P < 0.05)", TRUE),
        hr(),
        h4("Visualization Options"),
        radioButtons("chart_type", "Select Chart Type:",
                     choices = c("Scatter Plot" = "scatter", 
                                 "Bar Plot" = "bar"))
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Visualization", plotOutput("main_plot")),
          tabPanel("Data Table", tableOutput("filtered_table"))
        )
      )
    )
  ),
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
      )
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("PCA Plot", plotOutput("pca_plot")),
          tabPanel("t-SNE Plot", plotOutput("tsne_plot"))),
        textOutput("errorMessage_clustering")
      )
    )
  )
)



# Define Server
server <- function(input, output) {
  # Reactive dataset based on filters
  filtered_data <- reactive({
    data <- aggregated_data
    
    if (input$gene != "All") {
      data <- data %>% filter(Gene == input$gene)
    }
    if (input$phenotype != "All") {
      data <- data %>% filter(Phenotype == input$phenotype)
    }
    if (input$mouse_strain != "All") {
      data <- data %>% filter(Mouse_Strain == input$mouse_strain)
    }
    if (input$significant_only) {
      data <- data %>% filter(Aggregated_PValue < 0.05)
    }
    data
  })
  
  # Main plot output
  output$main_plot <- renderPlot({
    data <- filtered_data()
    
    data <- data %>%
      mutate(Significance = ifelse(Aggregated_PValue < 0.05, "Significant", "Not Significant"))
    
    if (input$chart_type == "scatter") {
      ggplot(data, aes(x = Phenotype, y = Aggregated_PValue, color = Significance)) +
        geom_point(size = 3) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue")) +
        labs(title = paste("Mouse Strain:", input$mouse_strain, " - P-Values by Phenotype"),
             x = "Phenotype", y = "Aggregated P-Value",
             color = "Significance")
    } else if (input$chart_type == "bar") {
      ggplot(data, aes(x = Phenotype, y = Aggregated_PValue, fill = Gene)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Bar Plot of Phenotypes", x = "Phenotype", y = "Aggregated P-Value")
    }
  })
  
  # Data table output
  output$filtered_table <- renderTable({
    filtered_data()
  })
  
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
  
  # PCA plot output
  output$pca_plot <- renderPlot({
    
  tryCatch({
    data <- filtered_data_pca()
    
    # Prepare data for PCA dynamically based on selection type
    pivot_column <- if (input$selection_type_clustering == "group") "Phenotype_Group" else "Phenotype"
    
    # Prepare data for PCA
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
    pca_result <-prcomp(standardized_data, center = TRUE, scale. = TRUE)
    
    # Convert PCA results to data frame and check for valid results
    pca_data <-as.data.frame(pca_result$x) %>%
        cbind(cluster = numeric_data$cluster)

    # Plot
    ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
      geom_point() +
      geom_text(aes(label = rownames(pca_data), vjust = -0.5, size = 3)) +
      ggtitle("PCA of Genes with Similar Phenotype Score") +
      theme_minimal()
    
  }, error = function(e){
    output$errorMessage_clustering <- renderText({
      print("Error: The data is not suitable or sufficient for PCA")
    })
    return(NULL)
  })
  
  })
  
  # t-SNE plot output
  output$tsne_plot <- renderPlot({
  tryCatch({
      data <- filtered_data_pca()
      
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
}
   

# Run the application
shinyApp(ui = ui, server = server)
