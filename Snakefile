import pandas
import os

def get_fastq(sample_sheet, run_id):
    """Query the sample sheet to return a list of fastq file(s) for the
    corresponding run_id
    """
    r1 = sample_sheet[sample_sheet.run_id == run_id].fastq_r1.iloc[0]
    fq = [r1]
    r2 = sample_sheet[sample_sheet.run_id == run_id].fastq_r2.iloc[0]
    if r2 is not None:
        fq.append(r2)
    return fq


def get_bam_for_library(sample_sheet, library_id, genome):
    """Query the sample sheet and returns the path to bam files for library_id
    """
    run_id = sample_sheet[(sample_sheet.library_id == library_id) & (sample_sheet.genome == genome)].run_id.unique()
    bam = [f'{genome}/bwa/{x}.bam' for x in run_id]
    return bam


def is_pe_library(sample_sheet, library_id):
    r2 = sample_sheet[sample_sheet.library_id == library_id].fastq_ftp_r2.iloc[0]
    if pandas.isna(r2) or r2 is None or r2 == '':
        return False
    else:
        return True


sample_sheet = pandas.read_csv(config['sample_sheet'], sep= ',', comment= '#')

if len(sample_sheet[sample_sheet.library_id == sample_sheet.run_id]) != 0:
    raise Exception('Error in sample sheet: library_id and run_id cannot be the same')


sample_sheet['fastq_r1'] = None
sample_sheet['fastq_r2'] = None

fastq_ftp = {}
for i,row in sample_sheet.iterrows():
    fq = 'fastq/' + os.path.basename(row.fastq_ftp_r1)
    sample_sheet.at[i, 'fastq_r1'] = fq
    fastq_ftp[fq] = row.fastq_ftp_r1

    if not pandas.isna(row.fastq_ftp_r2):
        fq = 'fastq/' + os.path.basename(row.fastq_ftp_r2)
        sample_sheet.at[i, 'fastq_r2'] = fq
        fastq_ftp[fq] = row.fastq_ftp_r2


fastqc_reports = []
for fq in fastq_ftp:
    fq = re.sub('.fastq.gz', '', os.path.basename(fq))
    fastqc_reports.append(fq + '_fastqc.zip')

wildcard_constraints:
    run_id= '|'.join([re.escape(x) for x in sample_sheet.run_id]),
    library_id= '|'.join([re.escape(x) for x in sample_sheet.library_id]),


rule all:
    input:
        'multiqc/fastqc_report.html',
        expand('{genome}/vep/all_lines.tsv', genome= sample_sheet[sample_sheet.library_type == 'WGS'].genome),
        expand('{genome}/mosdepth/dist.html', genome= sample_sheet[sample_sheet.library_type == 'WGS'].genome),
        expand('{genome}/edger/logrpkm.tsv', genome= sample_sheet[sample_sheet.library_type == 'RNAseq'].genome),
        expand('{genome}/gene_description/gene_names.csv', genome= sample_sheet[sample_sheet.library_type == 'WGS'].genome),
        expand('{genome}/tables/all_lines_gene_description.tsv', genome= sample_sheet[sample_sheet.library_type == 'WGS'].genome),
        expand('{genome}/tables/all_lines_filtered.tsv', genome= sample_sheet[sample_sheet.library_type == 'WGS'].genome),
        'PbergheiANKA/tables/all_lines_further_filtered_0.5.tsv',
        'PbergheiANKA/tables/all_lines_further_filtered_0.35.tsv',


rule download_fastq:
    output:
        fastq= 'fastq/{fq_name}.fastq.gz',
    params:
        ftp= lambda wc: fastq_ftp[f'fastq/{wc.fq_name}.fastq.gz']
    shell:
        r"""
        curl -L --output {output.fastq} {params.ftp}
        """


rule fastqc:
    input:
        fastq= 'fastq/{fq_name}.fastq.gz',
    output:
        qc= 'fastqc/{fq_name}_fastqc.zip',
    shell:
        r"""
        fastqc --noextract --outdir fastqc {input.fastq}
        """


rule multiqc:
    input:
        fastqc_reports= [f'fastqc/{x}' for x in fastqc_reports],
    output:
        'multiqc/fastqc_report.html',
    shell:
        r"""
        multiqc --force --outdir multiqc --filename fastqc_report.html {input.fastqc_reports}
        """


rule cutadapt:
    input:
        fastq= lambda wc: get_fastq(sample_sheet, wc.run_id),
    output:
        fastq_r1= temp('cutadapt/{run_id}_1.fastq.gz'),
        fastq_r2= temp(touch('cutadapt/{run_id}_2.fastq.gz')),
        report= 'cutadapt/{run_id}.log',
    run:
        if len(input.fastq) == 2:
            shell(r"""
            cutadapt --quality-cutoff 15 --minimum-length 10 --cores 4 \
                -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o {output.fastq_r1} -p {output.fastq_r2} {input.fastq} > {output.report}
            """)
        else:
            shell(r"""
            cutadapt --quality-cutoff 15 --minimum-length 10 --cores 4 \
                -a AGATCGGAAGAGC -o {output.fastq_r1} {input.fastq} > {output.report}
            """)


rule download_genome:
    output:
        fasta= 'ref/PlasmoDB-52_{genome}_Genome.fasta',
    shell:
        r"""
        curl https://plasmodb.org/common/downloads/release-52/{wildcards.genome}/fasta/data/PlasmoDB-52_{wildcards.genome}_Genome.fasta > {output.fasta}
        """


rule download_genome_annotations:
    output:
        gff= 'ref/PlasmoDB-52_{genome}.2.gff.gz',
    shell:
        r"""
        curl --silent -L https://plasmodb.org/common/downloads/release-52/{wildcards.genome}/gff/data/PlasmoDB-52_{wildcards.genome}.gff \
        | awk -v FS='\t' -v OFS='\t' '$1 !~ "^#" {{
                                        biotype = ";biotype="$3;
                                        if($3 == "protein_coding_gene") {{$3 = "gene"}}
                                        if($3 == "ncRNA_gene") {{$3 = "gene"}}
                                        print $0 biotype}}' \
        | sort -k1,1 -k4,4n -k5,5n -t$'\t' \
        | bgzip > ref/PlasmoDB-52_{wildcards.genome}.2.gff.gz

        tabix --force --preset gff ref/PlasmoDB-52_{wildcards.genome}.2.gff.gz
        """

rule bwa_index:
    input:
        fasta= 'ref/PlasmoDB-52_{genome}_Genome.fasta',
    output:
        index= 'ref/PlasmoDB-52_{genome}_Genome.fasta.bwt',
    shell:
        r"""
        bwa index {input.fasta}
        """


rule bwa_mem_align:
    input:
        fasta= 'ref/PlasmoDB-52_{genome}_Genome.fasta',
        index= 'ref/PlasmoDB-52_{genome}_Genome.fasta.bwt',
        fastq_r1= 'cutadapt/{run_id}_1.fastq.gz',
        fastq_r2= 'cutadapt/{run_id}_2.fastq.gz',
    output:
        bam= temp('{genome}/bwa/{run_id}.bam'),
        bai= temp('{genome}/bwa/{run_id}.bam.bai'),

    shell:
        r"""
        bwa mem -t 8 {input.fasta} {input.fastq_r1} {input.fastq_r2} \
        | samtools sort -@ 4 > {output.bam}

        samtools index {output.bam}
        """


rule merge_bam:
    input:
        bam= lambda wc: get_bam_for_library(sample_sheet, wc.library_id, wc.genome),
    output:
        bam= '{genome}/bwa/{library_id}.bam',
        md= '{genome}/bwa/{library_id}.md.txt',
        stats= '{genome}/bwa/{library_id}.stats'
    shell:
        r"""
        samtools merge -@ 4 - {input.bam} \
        | samtools sort -@ 4 -n \
        | samtools fixmate -m -@ 4 - - \
        | samtools sort -@ 4 \
        | samtools markdup -@ 4 -f {output.md} - {output.bam}

        samtools index -@ 4 {output.bam}
        samtools stats -@ 4 {output.bam} > {output.stats}
        """


rule mosdepth:
    input:
        bam= '{genome}/bwa/{library_id}.bam',
    output:
        txt= '{genome}/mosdepth/{library_id}.mosdepth.global.dist.txt',
    shell:
        r"""
        outdir=$(dirname {output.txt})
        mosdepth -x -n -t 4 $outdir/{wildcards.library_id} {input.bam}
        """


rule mosdepth_plot:
    input:
        txt= lambda wc: expand('{{genome}}/mosdepth/{library_id}.mosdepth.global.dist.txt',
            library_id= sample_sheet[(sample_sheet.library_type == 'WGS') & (sample_sheet.genome == wc.genome)].library_id.unique()),
    output:
        html= '{genome}/mosdepth/dist.html',
    shell:
        r"""
        python {workflow.basedir}/scripts/mosdepth_script_plot_dist.py --output {output.html} {input.txt}
        """


rule freebayes:
    input:
        bam= '{genome}/bwa/{library_id}.bam',
        fasta= 'ref/PlasmoDB-52_{genome}_Genome.fasta',
    output:
        vcf= '{genome}/freebayes/{library_id}.vcf',
    shell:
        r"""
        freebayes -f {input.fasta} --min-alternate-count 5 --ploidy 1 \
            --min-base-quality 5 \
            --min-alternate-fraction 0.1 --pooled-continuous {input.bam} \
        | awk -v OFS='\t' -v FS='\t' -v n=1 '{{if($0 !~ "^#"){{$3=n; n++}} print $0}}' \
        | {workflow.basedir}/scripts/add_allele_frequency_vcf.py - \
        | bcftools norm -m- -f {input.fasta} \
        | awk -v OFS='\t' -v name='{wildcards.library_id}' '{{if($1 == "#CHROM" && $10 == "unknown") {{$10 = name}} print $0}}' > {output.vcf}
        """


rule vep:
    input:
        vcf= '{genome}/freebayes/{library_id}.vcf',
        gff= 'ref/PlasmoDB-52_{genome}.2.gff.gz',
        fasta= 'ref/PlasmoDB-52_{genome}_Genome.fasta',
    output:
        vcf= '{genome}/vep/{library_id}.vcf',
    shell:
        r"""
        vep -i {input.vcf} --vcf --gff {input.gff} --fasta {input.fasta} --distance 1000 --force_overwrite -o {output.vcf}
        """


rule annotation_extract:
    input:
        vcf= '{genome}/vep/{library_id}.vcf',
    output:
        tsv= temp('{genome}/vep/{library_id}.tsv'),
    shell:
        r"""
        echo "chrom pos locus_id ref alt qual type ref_depth alt_depth alt_freq sum_alt_freq symbol gene feature biotype consequence impact amino_acids codons library_id" | tr ' ' '\t' > {output.tsv}

        bcftools +split-vep -d -f '%CHROM %POS %ID %REF %ALT %QUAL %TYPE [%AD{{0}}] [%AD{{1}}] [%ALT_AF] [%SUM_ALT_AF] %SYMBOL %Gene %Feature %BIOTYPE %Consequence %IMPACT %Amino_acids %Codons [%SAMPLE]\n' {input.vcf} | tr ' ' '\t' >> {output.tsv}
        """


rule concatenate_annotation:
    input:
        tsv= lambda wc: expand('{{genome}}/vep/{library_id}.tsv',
            library_id= sorted(sample_sheet[(sample_sheet.library_type == 'WGS') & (sample_sheet.genome == wc.genome)].library_id.unique())),
    output:
        tsv= '{genome}/vep/all_lines.tsv',
    shell:
        r"""
        awk 'NR == 1 || $1 != "chrom"' {input.tsv} > {output.tsv}
        """

rule gene_description:
    input:
        gff= 'ref/PlasmoDB-52_{genome}.2.gff.gz',
    output:
        csv= '{genome}/gene_description/gene_names.csv',
    shell:
        r"""
        zcat {input.gff} \
        | {workflow.basedir}/scripts/getGeneDescriptionFromGFF.py > {output.csv}
        """


rule merge_tables:
    input:
        tsv= '{genome}/vep/all_lines.tsv',
        csv= '{genome}/gene_description/gene_names.csv',
    output:
        tsv= '{genome}/tables/all_lines_gene_description.tsv',
    script:
        os.path.join(workflow.basedir, 'scripts/merge_tables.R')


rule variant_filtered:
    input:
        tsv= '{genome}/tables/all_lines_gene_description.tsv',
    output:
        tsv= '{genome}/tables/all_lines_filtered.tsv',
    script:
        os.path.join(workflow.basedir, 'scripts/variant_filtered.R')


rule variant_further_filtered_Pberghei:
    input:
        tsv= 'PbergheiANKA/tables/all_lines_filtered.tsv',
    output:
        tsv= 'PbergheiANKA/tables/all_lines_further_filtered_0.5.tsv',
    script:
        os.path.join(workflow.basedir, 'scripts/variant_further_filtered_0.5.R')


rule variant_f_filtered_Pberghei:
    input:
        tsv= 'PbergheiANKA/tables/all_lines_filtered.tsv',
    output:
        tsv= 'PbergheiANKA/tables/all_lines_further_filtered_0.35.tsv',
    script:
        os.path.join(workflow.basedir, 'scripts/variant_further_filtered_0.35.R')


rule hisat2_index:
    input:
        fasta= 'ref/PlasmoDB-52_{genome}_Genome.fasta',
    output:
        index= 'ref/PlasmoDB-52_{genome}_Genome.fasta.1.ht2',
    shell:
        r"""
        hisat2-build -p 7 {input.fasta} {input.fasta}
        """


rule hisat2_align:
    input:
        fasta= 'ref/PlasmoDB-52_{genome}_Genome.fasta',
        index= 'ref/PlasmoDB-52_{genome}_Genome.fasta.1.ht2',
        fastq_r1= 'cutadapt/{run_id}_1.fastq.gz',
        fastq_r2= 'cutadapt/{run_id}_2.fastq.gz',
        n_fastq= lambda wc: get_fastq(sample_sheet, wc.run_id),
    output:
        bam= temp('{genome}/hisat2/{run_id}.bam'),
        bai= temp('{genome}/hisat2/{run_id}.bam.bai'),
        log= '{genome}/hisat2/{run_id}.log',
    run:
        if len(input.n_fastq) == 2:
            shell(r"""
            hisat2 -p 7 --summary-file {output.log} --new-summary --max-intronlen 5000 \
                -x {input.fasta} -1 {input.fastq_r1} -2 {input.fastq_r2} \
            | samtools sort -@ 4 > {output.bam}

            samtools index {output.bam}
            """)
        else:
            shell(r"""
            hisat2 -p 7 --summary-file {output.log} --new-summary --max-intronlen 5000 \
                -x {input.fasta} -U {input.fastq_r1} \
            | samtools sort -@ 4 > {output.bam}

            samtools index {output.bam}
            """)


rule rename_hisat2:
    input:
        bam= lambda wc: '{genome}/hisat2/%s.bam' % sample_sheet[sample_sheet.library_id == wc.library_id].run_id.iloc[0],
        bai= lambda wc: '{genome}/hisat2/%s.bam.bai' % sample_sheet[sample_sheet.library_id == wc.library_id].run_id.iloc[0],
        log= lambda wc: '{genome}/hisat2/%s.log' % sample_sheet[sample_sheet.library_id == wc.library_id].run_id.iloc[0],
    output:
        bam= '{genome}/hisat2/{library_id}.bam',
        bai= '{genome}/hisat2/{library_id}.bam.bai',
        log= '{genome}/hisat2/{library_id}.log',
    run:
        # sanity check we have only one run per library
        assert len(sample_sheet[sample_sheet.library_id == wildcards.library_id]) == 1

        shell(r"""
        mv {input.bam} {output.bam}
        mv {input.bai} {output.bai}
        mv {input.log} {output.log}
        """)


rule featureCounts:
    input:
        gff= 'ref/PlasmoDB-52_{genome}.2.gff.gz',
        bam= '{genome}/hisat2/{library_id}.bam',
    output:
        tsv= '{genome}/featureCounts/{library_id}.tsv',
    params:
        is_pe= lambda wc: '-p' if is_pe_library(sample_sheet, wc.library_id) else '',
        strand= lambda wc: int(sample_sheet[sample_sheet.library_id == wc.library_id].featureCountsStrand.iloc[0]),
    shell:
        r"""
        featureCounts {params.is_pe} -T 10 -Q 10 -s {params.strand} -t exon -g gene_id -a {input.gff} -o {output.tsv} {input.bam}
        """


rule concatenate_featureCounts:
    input:
        tsv= lambda wc: expand('{{genome}}/featureCounts/{library_id}.tsv',
            library_id= sorted(sample_sheet[(sample_sheet.library_type == 'RNAseq') & (sample_sheet.genome == wc.genome)].library_id.unique())),
    output:
        tsv= '{genome}/featureCounts/all_lines.tsv',
    shell:
        r"""
        # First line of output is the column header
        echo "gene_id gene_length count library_id" | tr ' ' '\t' > {output.tsv}
        for x in {input.tsv}
	    do
	    library_id=`basename $x`     # Remove directory path from filename
    	library_id=${{library_id%.tsv}} # Remove .tsv suffix
        awk -v library_id="$library_id" -v OFS='\t' 'NR > 2 {{print $1, $6, $7, library_id}}' $x >> {output.tsv}
        done
        """


rule gene_normalization_edger:
    input:
        tsv= '{genome}/featureCounts/all_lines.tsv',
    output:
        tsv= '{genome}/edger/logrpkm.tsv',
    script:
        os.path.join(workflow.basedir, 'scripts/gene_normalization_edger.R')
