library(data.table) 
 
all_lines_gene_description <- fread(snakemake@input[['tsv']]) # Read  file
 
# Example of some filters
all_lines_filtered <- all_lines_gene_description[qual > 10 & sum_alt_freq > 0.8 & biotype != 'pseudogenic_transcript' & impact == 'HIGH']
 
# Write-out filtered table
write.table(x= all_lines_filtered, file= snakemake@output[['tsv']], sep= '\t', row.names= FALSE, quote= FALSE)
