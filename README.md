# Description
Malaria is a parasitic disease of great global concern. Understanding its transmission, through the identification of genes involved in this process might contribute to its complete eradication; in that, these genes could be targets for drug and/or vaccine development. This study developed a pipeline using bioinformatics programs for the analysis of genomes of gametocyte non-producer lines (GNPs), as well as gametocyte producer lines of _Plasmodium berghei and knowlesi_. The pipeline was also able to analyse RNA-seq datasets of these species. Through this analysis we were able to identify variants present in the genomes of these species, and their expression profiles at the different malaria life-cycle stages was determined based on the analysed RNA-seq datasets. The analysis was performed using the workflow manager, **snakemake**. The programs, commands and codes, and scripts used throughout the study are shown in the requirement.txt, Snakefile, and scripts files, respectively.

# Set up
Create and activate an environment dedicated to the running of programs using commands below:
  conda create -n my-env
  conda activate my-env

# Run pipeline
To run the pipeline, use command below:
  To know jobs to be run by snakemake use command:
    snakemake --dry-run --printshellcmds --jobs 10 --directory output --config sample_sheet=$PWD/sample_sheet.csv

  To have the jobs executed and done use command:
    snakemake --printshellcmds --jobs 10 --directory output --config sample_sheet=$PWD/sample_sheet.csv


# Results
**Fastq folder**
  contains details of fastq files analysed.

**Fastqc and multiqc folder**
  contains files generated after performing quality assessment on the fastq files. A summary of the quality assessment report for all the fastq files is shown in an html files.

**cutadapt folder**
  contains files that have had their reads trimmed and filtered after quality assessment, for an accurate and robust downstream analysis.

**ref folder**
  contains reference genomes and their annotation files downloaded from PlasmoDB, and indexed before read mapping.

**bwa and hisat2 folder**
contains BAM files, and their indexed bai files of the WGS and RNA-seq datasets generated after read mapping to the reference genome.

**mosdepth folder**
  contains files showing the calculated genome coverage of the libraries

**freebayes folder**
  contains VCF files showing the variant calling for each library

**vep folder**
  contains annotated VCF files of each library, and a tsv table - single concatenated- for all the annotated VCF files.

**gene_description**
  contains a tsv file showing gene_id and their gene_description.

**featureCounts folder**
  contains tsv files showing details of the exon of the gene - chromosome details, start and end position, and length. A summary of the assigned and unassigned reads is also provided.

**edger folder**
  contains a tsv file showing the information of normalisation of the reads.

**Pberghei/Pknowlesi-tables**
  contains tables generated at the end of the analysis, to which several thresholds were applied in order to obtain a list of variants.

**scripts**
  contains written scripts used for the study

**requirement.txt**
  contains bioinformatics programs, in their latest version used for the study

**sample_sheet.csv**
  contains information of the WGS and RNA-seq datasets used in the study

**Snakefile**
  contains the commands and codes used for the study, run using snakemake
