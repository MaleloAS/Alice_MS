library(data.table) 
 
vep <- fread(snakemake@input[['tsv']]) # Read VEP file
 
# Example of some filters
vep_short <- vep[library_id != 'PbANKA_820' & library_id != 'GNP_PbANKA_820' & library_id != 'm7_GNP_PbANKA_820']
 
# Write-out f_f_filtered table
write.table(x= vep_short, file= snakemake@output[['tsv']], sep= '\t', row.names= FALSE, quote= FALSE)
