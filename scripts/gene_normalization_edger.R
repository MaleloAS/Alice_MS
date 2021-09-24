library(data.table) 
library(edgeR)
 
counts <- fread(snakemake@input[['tsv']])
 
# Change table to wide-format with one column per library
counts <- dcast(data= counts, gene_id + gene_length ~ library_id, value.var= 'count')
 
# Put the gene lengths in a separate table and remove gene_length from counts table
gl <- counts[, list(gene_id, gene_length)]
counts[, gene_length := NULL]
 
# Turn the data.frame into a matrix
mat <- as.matrix(counts, rownames= 'gene_id')
 
# Estimate library sizes and compute expression values by normalising for gene
# length and library size
y <- DGEList(mat)
y <- calcNormFactors(y)
 
rpkm <- rpkm(y, log= TRUE, gene.length= gl$gene_length)
 
# Change RPKM table from "wide" format (one column per library) to "long"
# format (a single column for RPKM).
xrpkm <- as.data.table(rpkm, keep.rownames= 'gene_id')
xrpkm <- melt(data= xrpkm, id.vars= 'gene_id', variable.name= 'library_id', value.name= 'logrpkm')
 
write.table(xrpkm, snakemake@output[['tsv']], sep= '\t', row.names= FALSE, quote= FALSE)
