library(data.table)

all_lines <- fread(snakemake@input[['tsv']]) # Read VEP file
gene_names <- fread(snakemake@input[['csv']]) # Read gene_description file
gene_names <- gene_names[, list(gene_id, description)]

# Merge tables
all_lines_gene_description <- merge(all_lines, gene_names, by.x = 'gene', by.y = 'gene_id', all.x = TRUE, all.y = FALSE)

# Write-out filtered table
write.table(x= all_lines_gene_description, file= snakemake@output[['tsv']], sep= '\t', row.names= FALSE, quote= FALSE)
