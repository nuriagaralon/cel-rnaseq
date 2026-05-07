# EXPRESSION USING FEATURECOUNTS
# Abundance estimation from a BAM alignment
# Outputs a raw count matrix

rule featurecounts_expression:
    input:
        bam="results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        counts="results/expression/{sample}_gene_counts.tsv"
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
        -a {params.refgen} -o {output.counts} {input.bam} &>> {log}
        """