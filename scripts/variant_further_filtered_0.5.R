library(data.table) 
 
all_lines_filtered <- fread(snakemake@input[['tsv']]) # Read VEP file
 
# Example of some filters
all_lines_further_filtered_0.5 <- all_lines_filtered[ref_depth > 0 & alt_freq > 0.5 & description != 'PIR+protein' & description != 'fam-b+protein' & description != 'fam-c+protein' & description != 'conserved+Plasmodium+protein,+unknown+function' & description != 'conserved+protein,+unknown+function' & description != 'serine/threonine+protein+kinase+VPS15,+putative' & description != 'regulator+of+chromosome+condensation-PP1-interacting+protein']
 
# Write-out further filtered table
write.table(x= all_lines_further_filtered_0.5, file= snakemake@output[['tsv']], sep= '\t', row.names= FALSE, quote= FALSE)