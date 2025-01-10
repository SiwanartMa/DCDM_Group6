-- Populating the Disease_Sample_Records table mapping Disease_information and Sample_Records
INSERT INTO Diseases_Sample_Records (disease_ref_id, analysis_id, gene_accession_id)
SELECT d.id, s.analysis_id, s.gene_accession_id
FROM Disease_information d
JOIN Sample_Records s ON d.gene_accession_id = s.gene_accession_id;


-- Populating the Sample_Records_IMPC_Parameter table mapping Sample_Records and IMPC_Parameter
INSERT INTO Sample_Records_IMPC_Parameter (analysis_id, parameterId, impcParameterOrigId)
SELECT s.analysis_id , ip.parameterId , ip.impcParameterOrigId 
FROM IMPC_Parameter ip
JOIN Sample_Records s ON ip.parameterId = s.parameter_id ;
