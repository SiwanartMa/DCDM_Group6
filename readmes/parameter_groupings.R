# Author: Payel Sardar

# Setting the working directory
setwd("D:/KCL2024/Courses/7BBG1003_Data_Cleaning_Mgmt/Project_group6")

# Load necessary libraries
library(dplyr)
library(stringr)

# Load the dataset
parameter_data <- read.csv("impc_parameter.csv", stringsAsFactors = FALSE)

# Preview the data
head(parameter_data)

# Add a new column for grouping
parameter_data$parameter_group <- NA

grouping_rules <- list(
  Brain = c("brain", "cerebral cortex", "cerebellum", "hippocampus", 
            "hypothalamus", "olfactory lobe", "peripheral nervous system", 
            "neuro", "cognitive", "memory"),
  Liver = c("liver", "hepatic"),
  Kidney = c("kidney", "renal"),
  Cardiovascular = c("heart", "cardio", "vascular", "blood pressure"),
  Skin = c("skin", "epidermis", "dermal"),
  Blood = c("blood", "hemoglobin", "hematology"),
  Metabolic = c("metabolic", "glucose", "insulin"),
  Behavior = c("behavior", "activity", "locomotion", "movement"),
  Immune = c("immune", "immunity", "lymphocyte", "infection", "antibody"),
  Respiratory = c("lung", "respiratory", "pulmonary", "breathing"),
  Digestive = c("digestive", "stomach", "intestinal", "colon", "gut"),
  Reproductive = c("reproductive", "fertility", "sperm", "ovary", 
                   "uterus", "testis", "gonad"),
  Bone = c("bone", "skeletal", "joint"),
  Motor_Function = c("limb", "strength", "grip", "motor"),
  Weight_and_Body_Composition = c("weight", "body mass", "obesity"),
  Imaging = c("image", "imaging", "picture", "scan"),
  General_and_Equipment = c("comments", "general", "equipment"),
  Others = c("other")
)

# Assign groups
for (group_name in names(grouping_rules)) {
  keywords <- grouping_rules[[group_name]]
  pattern <- paste(keywords, collapse = "|")
  parameter_data$parameter_group <- ifelse(
    is.na(parameter_data$parameter_group) & str_detect(tolower(parameter_data$name), pattern),
    group_name,
    parameter_data$parameter_group
  )
}

parameter_data$parameter_group[is.na(parameter_data$parameter_group)] <- "Others"

# Save grouped data to a new file
write.csv(parameter_data, "Grouped_IMPC_Parameters.csv", 
          row.names = FALSE,
          fileEncoding =  "UTF-8")

# Summary of group counts
group_summary <- parameter_data %>% group_by(parameter_group) %>% summarise(count = n())
print(group_summary)
