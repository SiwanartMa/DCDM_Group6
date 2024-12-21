#Setting the working directory
setwd("/users/k24034333/DCDM_Group6/")

#Setting the data directory
data_path = "/scratch_tmp/grp/msc_appbio/DCDM_group6/originals/Group6/data"

#Reading the SOP
sop <- read.csv("/users/k24034333/DCDM_Group6/metadata/IMPC_SOP.csv")

#Listing files in the data directory
data_files <- list.files(path = data_path, full.names = TRUE)


#Examining the contents of one of the data files
sample_data <- read.csv(data_files[1], header = FALSE, 
                        col.names = c("dataField", "Values"))

#Function to Validate and Convert Data files based on the SOP
validate_and_convert <- function(value, expected_type) {
  tryCatch({
    if (expected_type == "String") {
      # Convert to string
      return(as.character(value))  
    } else if (expected_type == "Float") {
      # Convert to float
      return(as.numeric(value))  
    } else {
      # Return NA for unsupported types
      return(NA)  
    }
  }, error = function(e) {
    # Return NA if conversion fails
    return(NA)  
  })
}

# Initialize an empty dataframe for the merged data
merged_df <- data.frame(matrix(ncol = nrow(sop), nrow = 0))
colnames(merged_df) <- sop$dataField

# Loop through each data file
for (file in data_files) {
  
  # Read the File
  sample_file <- read.csv(file, header = FALSE, 
                          col.names = c("fieldnames", "fieldvalues"))
  
  # Reorder the File based on SOP dataField order
  reordered_data <- sample_file[match(sop$dataField, sample_file$fieldnames), ]
  
  # Initialize a vector to store validation results
  validation_results <- c()
  
  # Loop through each field in the reordered data to validate
  for (i in 1:nrow(reordered_data)) {
    field_name <- reordered_data$fieldnames[i]
    field_value <- reordered_data$fieldvalues[i]
    
    # Find the SOP rule for the field
    sop_row <- sop[sop$dataField == field_name, ]
    
    # If SOP rule exists, validate the field value
    if (nrow(sop_row) > 0) {
      datatype <- sop_row$dataType
      valid_value <- validate_and_convert(field_value, datatype)
      
      # Append the validation result to the vector
      validation_results <- c(validation_results, valid_value)
    } else {
      # If no SOP rule exists, mark as invalid
      validation_results <- c(validation_results, "No SOP Rule")
    }
  }
  
  # Add the validated and reordered data as a row to the final data frame
  merged_df <- rbind(merged_df, validation_results)
}

# Adding column names
colnames(merged_df) <- sop$dataField

# saved the merged dataframe
write.csv(merged_df, "/users/k24034333/DCDM_Group6/outputs/tables/Merged_df_final.csv")
