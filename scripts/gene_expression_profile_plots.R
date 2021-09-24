library(data.table) 
library(edgeR)
library(ggplot2)

# Reading tables
ss= read.table("sample_sheet.tsv", header=FALSE, sep="\t")
gene_names= read.table("gene_names.tsv", header=FALSE, sep="\t")

# Reload tables using header and row names
ss= read.table("sample_sheet.tsv", header=TRUE, sep="\t")
gene_names= read.table("gene_names.tsv", header=TRUE, sep="\t")

# Load only the columns we need
ss <- fread('sample_sheet.tsv', select= c('library_id', 'time_point', 'rapamycin', 'library_type', 'genome','stage'))

# Not strictly necessary but let's do it:
ss <- ss[library_type %in% 'RNAseq' & genome == 'PbergheiANKA']

# Read the all_lines table
counts <- fread('all_lines.tsv')

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

# Write table
logrpkm <- xrpkm

# list of selected mutated genes for K173
mutated_genes <- c('PBANKA_0615700', 'PBANKA_0704400', 'PBANKA_0827100', 'PBANKA_1020500', 'PBANKA_1407700', 'PBANKA_0609500', 'PBANKA_1446500') 

# Now dat should have all the information we need for plotting in a suitable format:
dat <- merge(logrpkm[gene_id %in% mutated_genes], ss, by= 'library_id')

# Add gene description information. `gene_description` is a table with columns "gene_id" and "description"
dat <- merge(dat, gene_names, by= 'gene_id')

#Join gene_id and description
dat[, gene_desc := paste(gene_id, description)]

# Replace ugly "+" signs
dat[, gene_desc := gsub('+', ' ', gene_desc, fixed= TRUE)]

# Plot individual libraries, colored by Rapamycin variable
all_k173 <- ggplot(data= dat[!time_point %in% c(-1, -2, -3, -4)], aes(x= time_point, y= logrpkm, color= rapamycin, group= rapamycin)) + 
  geom_point() + 
  facet_wrap(~gene_id) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
all_k173

# Add line connecting averages of each group
stage_avg <- dat[, list(logrpkm= mean(logrpkm)), by= list(gene_id, time_point, rapamycin)]

all_k173 <- all_k173 + geom_line(data= stage_avg[!time_point %in% c(-1, -2, -3, -4)])
all_k173

# Add males and females. Individual libraries and averages
all_k173 <-  all_k173 +
  geom_col(data= stage_avg[time_point == -1], colour= 'dark blue', fill= 'grey80') +
  geom_col(data= stage_avg[time_point == -2], colour= 'dark blue', fill= 'grey80') +
  geom_col(data= stage_avg[time_point == -3], colour= 'red', fill= 'grey80') +
  geom_col(data= stage_avg[time_point == -4], colour= 'red', fill= 'grey80') +
  geom_point(data= dat[time_point == -1], colour= 'dark blue') +
  geom_point(data= dat[time_point == -2], colour= 'dark blue') +
  geom_point(data= dat[time_point == -3], colour= 'red') +
  geom_point(data= dat[time_point == -4], colour= 'red') 
all_k173