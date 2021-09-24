Date 02-06-2021 Setting-Up

Set up a cluster account, install VPN client, download mobaxterm, install Atom and configure it to cluster. install package- ftp-remote-edit.
Use cluster account details to login, always connect to the ssh username@headnode03.cent.gla.ac.uk, then connect to username@bioinf03.hpc.gla.ac.uk

All programs will managed using conda package manager. To install bioconda repositories use command
  'curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh, sh Miniconda3-latest-Linux-x86_64.sh'
To use bioconda, configuration of conda is done using commands
  conda config --add channels defaults,
  conda config --add channels bioconda,
  conda config --add channels conda-forge
A faster version of conda- mamba- is installed using command
  'conda install mamba -n base -c conda-forge'

Conda environment
  set up a conda environment where all needed programs will be installed. To initialize a conda environment use command 'conda create -n my-env'. my-env can be replaced with preferred name that only contains letters, digits and underscore.
  To activate conda use command
    'conda activate my-env'

Set up a project directory.
  It is important to set up a directory specifically for your project in which you will keep all files for you project.
  In this directory have a:
    README.md file which acts as the lab manual,
    requirements.txt file that will contain programs with latest version details.


Date 22-06-2021 Data collection from ENA, and Sample_sheet table

Collect data for plasmodium berghei from ENA (https://www.ebi.ac.uk/ena/browser/home)  under accession numbers ERP000253, ERS004442, ERS002990.
Download supplementary table of Sinha et al 2014 from (https://static-content.springer.com/esm/art%3A10.1038%2Fnature12970/MediaObjects/41586_2014_BFnature12970_MOESM44_ESM.xlsx) to use the specific names for the parasite lines.
Create a table with columns:
  library_name
  secondary_sample_accession
  line (i.e. K173, 820 etc.)
  fastq_ftp_r1 (link to raw data file for first-in-pair reads)
  fastq_ftp_r2 (link to raw data file for second-in-pair reads)
  sample_alias, and
  links from which sample lines were obtained.
Save table with appropriate name in a folder on your PC.


Date 23-06-2021 Adding programs to requirements.txt and Installing using mamba
log into the cluster, activate conda environment, and change to your project directory

To the file requirements.txt, add the following programs for a start :
  bwa =0.7.17
  freebayes =1.3.5
  snakemake =6.4.1
  hisat2 =2.2.1
  fastqc =0.11.9
  samtools
  bioconductor-edger =3.34
  ensembl-vep =104.3
  multiqc =1.10.1
Install programs using command
  'mamba install --freeze-installed -n my-env --yes --file requirements.txt'

Note: above programs and versions can be obtains from google search 'my-program bioconda'. for example for  bwa link is 'https://bioconda.github.io/recipes/bwa/README.html'


Date 24-06-2021 Downloading Reference genome, gene annotation file, and Fastq files.

Go to plasmodb.org website to download the Reference genome in fasta format and gene annotation list in GFF format.
  On the homepage click on the 'Data' tab and select 'Download data files'.
  In the new window click on 'current release' folder and select 'Pberghei ANKA' folder
  Then click on the 'DATA' folder and select fasta format and download.
  For the gene annotation list use the same steps but in 'DATA' folder select gff format and download.
  Save files appropriate.

Download the fastq_ftp_r1 and fastq_ftp_r2 in the sample_sheet excel sheet table created on 22-06-2021.
  Save them in your PC. No need to extract the zipped files.


Date 25-06-2021 Copying and Running Sample_sheet.tsv and Snakefile files.
log into the cluster, activate conda environment, and change to your project directory

On the cluster under '/export/home/db291g/Tritume/alice/' sample_sheet.tsv and Snakefile were copied to the project directory.
  Note: Sample_sheet.tsv is excel sheet with few edits from Dario. Snakefile is a file given to snakemake program to enable execution of various jobs for the project.

To show which commands will be executed in snakemake,
  execute command 'snakemake --dry-run --printshellcmds --jobs 1 --directory output --config sample_sheet=$PWD/sample_sheet.tsv'
  To run the jobs, execute command 'snakemake --printshellcmds --jobs 1 --directory output --config sample_sheet=$PWD/sample_sheet.tsv'
The command generates an output sub_directory in you project directory which contains .snakemake folder, fastq folder, fastqc folder which contains zipped and html(shows generated report from the fastqc run on the fastq files in sample_sheet.tsv- download these and view in web browser), and multiqc folder.


Date 28-06-2021 Analysis of multiqc/fastqc_report.html, Reference and annotation genome file download, and installation of cutadapt
log into the cluster, activate conda environment, and change to your project directory

Transfer the generated 'multiqc/fastqc_report.html'from the cluster to PC and view it using a web browser.
  The file might contain three signs seem to these of 'traffic lights' -good/warning/fail.
  few bad quality reads were seen in ERR019300, which require exclusion from the sequencing run.
To trim bad read, use cutadapt. Get the latest version from 'my-program bioconda'
  Add it to the requirements.txt file and install using 'mamba install --freeze-installed --yes --file requirements.txt'

To download the reference genome in fasta format, and the gene annotation in GFF format using snakemake, add to the Snakefile a rule to download.
  Refer to rules in the Snakefile.


Date 29-06-2021 Editing sample_sheet.tsv and Snakefile, running snakemake (Trimming, reference indexing, read mapping, sorting, and indexing)
log into the cluster, activate conda environment, and change to your project directory

To exclude fastq file with poor quality reads, edit the sample_sheet.tsv, by adding the prefix "#" at the beginning of the row. in this case row with ERR019300.
A command is added to the Snakefile to recognize the exclusion of this line. Refer to Snakefile for command added.
To the sample_sheet.tsv file, another column is added 'run_id' this is useful for line with multiple runs, for example line K173 in this case.
The run_id is obtained from the fastq URL, for example: "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR264/ERR264579/ERR264579_2.fastq.gz" the run ID would be ERR264579
To recognize the added column, a command line is added just at the top of the Snakefile after the 'import' statements. Refer to the Snakefile for command.

To the Snakefile, a line is added to the rule all, just after the statement 'multiqc/fastqc_report.html'
  Refer to Snakefile for the command - 'expand('bwa/{run_id}.bam', run_id= sample_sheet.run_id),'

A cutadapt rule is added to the Snakefile.
  The command specifies the Phred score to be used, the which reads to be trimmed, and that all reads below a length to be cut off.
  Refer to Snakefile for the rules.

A bwa-index rule for the reference genome is also added to the Snakefile.
  The command is bwa index, and the output of the indexing is a *.bwt file.
  This produces five types of files (amb, ann, bwt, pac, and sa)

A mapping, sorting, and indexing rule for the reads is added also to the Snakefile
  The read mapping rule(bwa-mem-align) is piped to the sort and index rule (samtools)
  The fastq files are aligned to the reference genome, using bwa-mem, resulting into a bam file.
  The output.bam file will then be sorted and indexed using samtools, for easy access to their genomic locations.

Execution of these added rules to Snakefile, are run first by:
  checking jobs to be run by using command 'snakemake --dry-run --printshellcmds --jobs 10 --directory output --config sample_sheet=$PWD/sample_sheet.tsv',
  then to actually run these jobs used 'snakemake --printshellcmds --jobs 10 --directory output --config sample_sheet=$PWD/sample_sheet.tsv'

To view files generated after jobs are done running:
  go to the project directory, check the output folder that will contain files according to your output commands names instructed in the Snakefile.


Date 30-06-2021 Merging bam files and Variant calling using freebayes
log into the cluster, activate conda environment, and change to your project directory

Note: Merging is done for bam files belonging to the same library, in our case K173, and those within the same job, marking duplicate reads is done. Duplicate reads are those mapping exactly the same genomic locations, and these are mostly PCR duplicates.

Adding programs to directory
  Create a sub_directory (scripts), to it a program copied from '/export/home/db291g/Tritume/add_allele_frequency_vcf.py' is added. This program adds to the output of freebayes information that facilitates the filtering of variants of interest.
  To the requirements.txt file add programs and run 'mamba install ...':
    bcftools without its version
    pysam with latest version

Editing Snakefile
  Add another function immediately after the 'get_fastq' function. Refer to Snakefile for function
    The command returns bam files for each run_id.
  Before the 'rule all':
    add a function -refer to Snakefile- that will tell snakemake to take the values for run_id and library_id exactly as they are in the sample_sheet.
  To the 'rule all':
    delete the 'expand('bwa/{run_id}.bam', run_id= sample_sheet.run_id),'
    add 'expand('freebayes/{library_id}.vcf', library_id= sample_sheet.library_id),', which will give a list of output files for freebayes
    Note: Our interest is the final output file.
  Add the rule merge_bam
    refer to Snakefile for rule.
  Add rule freebayes
      this rule is piped to 'bcftools norms' which applies some reformatting to the output of freebayes. It is also piped to program in scripts directory.

Execution of snakemake commands
  use the 'snakemake --dry-run.....' command to see jobs to be done.
  without '--dry-run', execute the actual jobs.


Date 01-07-2021 Variant annotation using vep
  log into the cluster, activate conda environment, and change to your project directory

  The freebayes will detect the variants but will not show what type of variants they are. To be able to know the type of variants detected and their biological context, variant annotation is done using ensembl-vep.

  To carry out the variant annotation of the files generated by freebayes, we use:
    Input files:
      vcf generated by freebayes (Refer to 30/06/2021).
      reference 52_PbergheiANKA_Genome fasta file (refer to rule download on 28/06/2021).
      52_PbergheiANKA_genome_annotations gff file (refer to rule download on 28/06/2021).
    output file:
      this will be a vcf file similar to that of freebayes. Refer to Snakefile for output function.

  Editing Snakefile
    To the 'rule all':
      delete the 'expand('freebayes/{library_id}.vcf', library_id= sample_sheet.library_id),'
      add 'expand('vep/{library_id}.vcf', library_id= sample_sheet.library_id),', which will give a list of output files for vep.
    Add rule vep
        the vep command is:
          'vep -i {input.vcf} --vcf --gff {input.gff} --fasta {input.fasta} --distance 1000 --force_overwrite -o {output.vcf}'

  Execution of snakemake commands
    use the 'snakemake --dry-run.....' command to see jobs to be done.
    without '--dry-run', execute the actual jobs.


Date 02-07-2021 Extraction of annotation information into a readable table.
log into the cluster, activate conda environment, and change to your project directory

The vcf files generated from the variant annotation step (01-07-2021) can be difficult to read and make sense out of them because it contains a lot of details. To make it readable, parsable, and easy to understand, extraction of the annotation information is done and put in a table with the following columns:
  chromosome and position of the variant
  reference allele (REF)
  alternate allele (ALT)
  qual
  type of variant (snp or indel)
  reference allele ([%AD{{0}}])
  alternate allele ([%AD{{1}}])
  alternate allele frequency ([%ALT_AF])
  sum of all alternate allele frequencies at this locus ([%SUM_ALT_AF]). This equals ALT_AF if there is only one alternate allele
  name of feature (typically the PBANKA transcript ID)
  type of feature-biotype-(typically mRNA)
  Consequence, which is the type of modification that the variant causes on the transcript (if it is in or near a transcript)
  Impact, this is the subjective assessment of the severity of the variant (for more details check https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html).
  Amino_acids and codons which are the changes in sequence in the protein and transcript

To carry out the annotation_extract of the files generated by vep, we use:
  Input file:
    vcf generated by vep (Refer to 01-07-2021).
  output file:
    this will be a tsv file since we need a table. Refer to Snakefile for output function.

Editing Snakefile
  To the 'rule all':
    delete the 'expand('vep/{library_id}.vcf', library_id= sample_sheet.library_id),'
    add 'expand('vep/{library_id}.tsv', library_id= sample_sheet.library_id),', which will give a list of output files in the vep folder as *.tsv.
  Add rule annotation_extract
      the annotation_extract command is:
        echo "chrom pos ref alt type ref_depth alt_depth alt_freq sum_alt_freq feature biotype consequence impact amino_acids codons library_id" | tr ' ' '\t' > {output.tsv}

        bcftools +split-vep -d -f '%CHROM %POS %REF %ALT %TYPE [%AD{0}] [%AD{1}] [%ALT_AF] [%SUM_ALT_AF] %Feature %BIOTYPE %Consequence %IMPACT %Amino_acids %Codons [%SAMPLE]\n' {input.vcf} | tr ' ' '\t' >> {output.tsv}

  Note:
    The first command prepares the header line of the output table and writes to file {output.tsv}.
    The second command appends to {output.tsv} the actual data. We use the program bcftools to extract from each line the information we are interested in.
  Refer to Snakefile for rule annotation_extract.

Execution of snakemake commands
  use the 'snakemake --dry-run.....' command to see jobs to be done.
  without '--dry-run', execute the actual jobs.
