# EXPRESSION USING STRINGTIE
# Abundance estimation from a BAM alignment
# Outputs a gtf and an abundance matrix (no raw counts)

rule stringtie_expression:
    input:
        bam="results/alignment/{sample}.bam"
    output:
        gtf="results/expression/{sample}/{sample}.gtf",
        counts="results/expression/{sample}/{sample}.gene_abund.anno.tab"
    params:
        refgen=config["genome"]["annotation_file"]
    threads: 8
    log:
        "workflow/logs/stringtie_expression/{sample}.log"
    benchmark:
        repeat("workflow/benchmarks/stringtie_expression/{sample}.tsv", 3)
    conda:
        "../envs/stringtie.yaml"
    shell:
        """
        stringtie {input.bam} -p {threads} --rf -l {wildcards.sample} \
        -o {output.gtf} -G {params.refgen} -A {output.counts} -eB &>> {log}
        """

# Measure stringtie accuracy

rule stringtie_quality:
    input:
        gtf="results/expression/{sample}/{sample}.gtf"
    output:
        stats="results/expression/{sample}/compare.stats"
    params:
        refgen=config["genome"]["annotation_file"],
        prefix="results/expression/{sample}/compare"
    threads: 2
    log:
        "workflow/logs/stringtie_quality/{sample}.log"
    benchmark:
        repeat("workflow/benchmarks/stringtie_quality/{sample}.tsv", 3)
    conda:
        "../envs/stringtie.yaml"
    shell:
        """
        gffcompare -R -r {params.refgen} {input.gtf} -o {params.prefix} &>> {log}
        """

# EXPRESSION USING HTSEQ
# Abundance estimation from a BAM alignment
# Outputs a raw count matrix

rule htseq_expression:
    input:
        bam="results/alignment/star/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        counts="results/expression/htseq/{sample}_gene_counts.txt"
    params:
        refgen=config["genome"]["annotation_file"]
    threads: 1
    log:
        "workflow/logs/htseq_expression/{sample}.log"
    benchmark:
        repeat("workflow/benchmarks/htseq_expression/{sample}.tsv", 3)
    conda:
        "../envs/htseq.yaml"
    shell:
        """
        htseq-count -r pos -s reverse -t exon -i gene_id \
        {input.bam} {params.refgen} > {output.counts} 2>> {log}
        """

# EXPRESSION USING FEATURECOUNTS
# Abundance estimation from a BAM alignment
# Outputs a raw count matrix

rule featurecounts_expression:
    input:
        bam="results/alignment/star/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        counts="results/expression/featurecounts/{sample}_gene_counts.tsv"
    params:
        refgen=config["genome"]["annotation_file"]
    threads: 1
    log:
        "workflow/logs/featurecounts_expression/{sample}.log"
    benchmark:
        repeat("workflow/benchmarks/featurecounts_expression/{sample}.tsv", 3)
    conda:
        "../envs/featurecounts.yaml"
    shell:
        """
        featureCounts -t exon -g gene_id -s 2 -p -T {threads}\
        -a {params.refgen} -o {output.counts} {input.bam} &>> {log}
        """

# EXPRESSION USING SALMON
# Index the reference transcriptome

rule salmon_index:
    input:
        transcriptome=config["genome"]["transcriptome"]
    output:
        tindex=directory("results/expression/salmon/{genome}_index")
    threads: 1
    log:
        "workflow/logs/salmon_index/{genome}.log"
    benchmark:
        repeat("workflow/benchmarks/salmon_index/{genome}.tsv", 3)
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon index -t {input.transcriptome} -i {output.tindex} -p {threads} &>> {log}
        """
# Abundance estimation via quasi-mapping
# Outputs a raw count matrix

rule salmon_expression:
    input:
        fwfasta="results/preprocessed/{sample}_R1.trimmed.fastq.gz",
        rvfasta="results/preprocessed/{sample}_R2.trimmed.fastq.gz",
        index=f"results/expression/salmon/{config['genome']['genome_name']}_index"
    output:
        counts="results/expression/salmon/{sample}/quant.sf"
    params:
        outdir="results/expression/salmon/{sample}",
        library="ISR"
    threads: 8
    log:
        "workflow/logs/salmon_expression/{sample}.log"
    benchmark:
        repeat("workflow/benchmarks/salmon_expression/{sample}.tsv", 3)
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon quant -i {input[2]} -l {params.library} \
        -1 {input[0]} -2 {input[1]} \
        -p {threads} -o {params.outdir} &>> {log}
        """

# EXPRESSION USING RSEM
# Create reference from genome and annotation

rule rsem_reference:
    input:
        genome=config["genome"]["genome_file"],
        gtf=config["genome"]["annotation_file"]
    output:
        expand("results/expression/rsem/{{genome}}_rsem_reference.{file}", file=["grp", "ti", "transcripts.fa", "seq", "chrlist", "idx.fa", "n2g.idx.fa"])
    params:
        prefix="results/expression/rsem/{genome}_rsem_reference"
    threads: 4
    log:
        "workflow/logs/rsem_reference/{genome}_rsem_reference.log"
    benchmark:
        repeat("workflow/benchmarks/rsem_reference/{genome}_rsem_reference.tsv", 3)
    conda:
        "../envs/rsem.yaml"
    shell:
        """
        rsem-prepare-reference --gtf {input.gtf} \
            {input.genome} {params.prefix} &>> {log}
        """

# Abundance estimation from a BAM alignment
# Outputs a raw count matrix of genes and transcripts

rule rsem_expression:
    input:
        trbam="results/alignment/star/{sample}_Aligned.toTranscriptome.out.bam",
        ref=f"results/expression/rsem/{config['genome']['genome_name']}_rsem_reference.transcripts.fa"
    output:
        countsgen="results/expression/rsem/{sample}.genes.results",
        countst="results/expression/rsem/{sample}.isoforms.results"
    params:
        reference=f"results/expression/rsem/{config['genome']['genome_name']}_rsem_reference",
        prefix="results/expression/rsem/{sample}"
    threads: 8
    log:
        "workflow/logs/rsem_expression/{sample}.log"
    benchmark:
        repeat("workflow/benchmarks/rsem_expression/{sample}.tsv", 3)
    conda:
        "../envs/rsem.yaml"
    shell:
        """
        rsem-calculate-expression --paired-end --strandedness reverse -p {threads} \
            --no-bam-output --alignments {input.trbam} {params.reference} {params.prefix} &>> {log}
        """