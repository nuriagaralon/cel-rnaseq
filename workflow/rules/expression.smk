rule stringtie_expression:
    input:
        "results/alignment/{sample}.bam"
    output:
        "results/expression/{sample}/{sample}.gtf",
        "results/expression/{sample}/{sample}.gene_abund.anno.tab"
    params:
        refgen=config["genome"]["annotation_file"]
    threads: 8
    log:
        "workflow/logs/stringtie_expression/{sample}.log"
    benchmark:
        "workflow/benchmarks/stringtie_expression/{sample}.tsv"
    conda:
        "../envs/stringtie.yaml"
    shell:
        """
        stringtie {input} -p {threads} --rf -l {wildcards.sample} \
        -o {output[0]} -G {params.refgen} -A {output[1]} -eB &>> {log}
        """

rule stringtie_quality:
    input:
        "results/expression/{sample}/{sample}.gtf"
    output:
        "results/expression/{sample}/compare.stats"
    params:
        refgen=config["genome"]["annotation_file"],
        prefix="results/expression/{sample}/compare"
    threads: 2
    log:
        "workflow/logs/stringtie_quality/{sample}.log"
    benchmark:
        "workflow/benchmarks/stringtie_quality/{sample}.tsv"
    conda:
        "../envs/stringtie.yaml"
    shell:
        """
        gffcompare -R -r {params.refgen} {input} -o {params.prefix} &>> {log}
        """

rule htseq_expression:
    input:
        "results/alignment/star/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "results/expression/htseq/{sample}_gene_counts.txt"
    params:
        refgen=config["genome"]["annotation_file"]
    threads: 1
    log:
        "workflow/logs/htseq_expression/{sample}.log"
    benchmark:
        "workflow/benchmarks/htseq_expression/{sample}.tsv"
    conda:
        "../envs/htseq.yaml"
    shell:
        """
        htseq-count -r pos -s reverse -t exon -i gene_id \
        {input} {params.refgen} > {output} 2>> {log}
        """
# --additional-attr gene_name

rule featurecounts_expression:
    input:
        "results/alignment/star/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "results/expression/featurecounts/{sample}_gene_counts.tsv"
    params:
        refgen=config["genome"]["annotation_file"]
    threads: 1
    log:
        "workflow/logs/featurecounts_expression/{sample}.log"
    benchmark:
        "workflow/benchmarks/featurecounts_expression/{sample}.tsv"
    conda:
        "../envs/featurecounts.yaml"
    shell:
        """
        featureCounts -t exon -g gene_id -s 2 -p -T {threads}\
        -a {params.refgen} -o {output} {input} &>> {log}
        """

rule salmon_index:
    input:
        config["genome"]["transcriptome"]
    output:
        "results/expression/salmon/index"
    params:
        genome=config["genome"]["genome_name"]
    threads: 1
    log:
        "workflow/logs/salmon_index/{params.genome}.log"
    benchmark:
        "workflow/benchmarks/salmon_index/{params.genome}.tsv"
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon index -t {input} -i {output} -p {threads}
        """

rule salmon_expression:
    input:
        "results/preprocessed/{sample}_R1.trimmed.fastq.gz",
        "results/preprocessed/{sample}_R2.trimmed.fastq.gz",
        "results/expression/salmon/index"
    output:
        "results/expression/salmon/{sample}/quant.sf"
    params:
        outdir="results/expression/salmon/{sample}",
        library="ISR"
    threads: 8
    log:
        "workflow/logs/salmon_expression/{sample}.log"
    benchmark:
        "workflow/benchmarks/salmon_expression/{sample}.tsv"
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon quant -i {input[2]} -l {params.library} \
        -1 {input[0]} -2 {input[1]} \
        -p {threads} -o {params.outdir}
        """