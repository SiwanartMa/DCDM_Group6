-- Diseases table
CREATE TABLE `Disease_information` (
  `id` INT AUTO_INCREMENT PRIMARY KEY,       -- Surrogate primary key
  `disease_id` VARCHAR(50) NOT NULL,         -- Non-unique disease identifier
  `disease_term` TEXT NOT NULL,      
  `gene_accession_id` VARCHAR(15) NOT NULL,  
  `phenodigm_score` DECIMAL(10, 5),
  INDEX `idx_disease_gene` (`gene_accession_id`), -- Index for filtering/joining by gene_accession_id
  INDEX `idx_disease_id` (`disease_id`)           -- Index for disease_id lookups
);

-- Sample_Records table
CREATE TABLE `Sample_Records` (
  `analysis_id` VARCHAR(15) PRIMARY KEY,     -- Unique identifier for the analysis/sample
  `gene_accession_id` VARCHAR(15) NOT NULL, 
  `gene_symbol` VARCHAR(15) NOT NULL,       
  `mouse_strain` VARCHAR(10) NOT NULL,      
  `mouse_life_stage` VARCHAR(25),           
  `parameter_id` VARCHAR(25) NOT NULL,      
  `parameter_name` VARCHAR(255),            
  `pvalue` DECIMAL(10, 5),
  INDEX `idx_sample_gene` (`gene_accession_id`), -- Index for filtering/joining by gene_accession_id
  INDEX `idx_sample_param` (`parameter_id`)      -- Index for filtering by parameter_id
);

-- Diseases_Sample_Records table, a mapping table between Disease_information and Sample Records
CREATE TABLE `Diseases_Sample_Records` (
  `disease_ref_id` INT NOT NULL,             -- References the `id` column in Diseases
  `analysis_id` VARCHAR(30) NOT NULL,        
  `gene_accession_id` VARCHAR(50) NOT NULL,  -- links Diseases and Samples
  PRIMARY KEY (`disease_ref_id`, `analysis_id`, `gene_accession_id`), -- Composite key
  FOREIGN KEY (`disease_ref_id`) REFERENCES `Disease_information` (`id`) ON DELETE CASCADE,
  FOREIGN KEY (`analysis_id`) REFERENCES `Sample_Records` (`analysis_id`) ON DELETE CASCADE,
  INDEX `idx_dsr_disease_ref` (`disease_ref_id`), -- Index for filtering/joining by disease_ref_id
  INDEX `idx_dsr_analysis` (`analysis_id`),       -- Index for filtering by analysis_id
  INDEX `idx_dsr_gene` (`gene_accession_id`)      -- Index for filtering by gene_accession_id
);

-- IMPC_procedure table
CREATE TABLE `IMPC_Procedure` (
  `impcParameterOrigId` bigint NOT NULL PRIMARY KEY,
  `name` varchar(255) NOT NULL,
  `description` text,
  `isMandatory` boolean NOT NULL
);

-- Creating IMPC_Parameter table
CREATE TABLE `IMPC_Parameter` (
  `impcParameterOrigId` BIGINT NOT NULL PRIMARY KEY,
  `name` TEXT NOT NULL,
  `description` TEXT,
  `parameterId` TEXT NOT NULL,
  `parameter_group` TEXT NOT NULL,
  CONSTRAINT `fk_impc_parameter_impc_procedure` 
    FOREIGN KEY (`impcParameterOrigId`) 
    REFERENCES `IMPC_Procedure` (`impcParameterOrigId`)
    ON UPDATE CASCADE -- to propagate updates to IMPC_Procedure
);

-- Sample_Records_IMPC_Parameter table (junction table to link Sample_Records and IMPC_Parameter)
CREATE TABLE `Sample_Records_IMPC_Parameter` (
  `analysis_id` VARCHAR(15) NOT NULL,                    -- Foreign key to Sample_Records
  `parameterId` VARCHAR(50) NOT NULL,                    -- link between Sample Records and IMPC_Parameter
  `impcParameterOrigId` BIGINT NOT NULL,                 -- Foreign key to IMPC_Parameter
  PRIMARY KEY (`analysis_id`, `parameterId`, `impcParameterOrigId`), -- Composite primary key
  FOREIGN KEY (`analysis_id`) REFERENCES `Sample_Records` (`analysis_id`) ON DELETE CASCADE, -- Link to Sample_Records
  FOREIGN KEY (`impcParameterOrigId`) REFERENCES `IMPC_Parameter` (`impcParameterOrigId`) ON DELETE CASCADE, -- Link to IMPC_Parameter
  INDEX `idx_srip_analysis_id` (`analysis_id`),           -- Index for filtering by analysis_id
  INDEX `idx_srip_parameter_id` (`parameterId`),           -- Index for filtering by parameterId
  INDEX `idx_srip_impcParameterOrigId` (`impcParameterOrigId`) -- Index for filtering by impcParameterOrigId
);
