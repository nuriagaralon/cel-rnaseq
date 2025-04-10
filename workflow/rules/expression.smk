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
        "results/alignment/star/{sample}_Aligned.toTranscriptome.out.bam",
        "results/alignment/{sample}.bam"
    output:
        "results/expression/{sample}/{sample}.gtf",
        "results/expression/{sample}/{sample}.gene_abund.anno.tab"
    params:
        refgen=config["genome"]["annotation_file"]
    threads: 8
    log:
        "workflow/logs/htseq_expression/{sample}.log"
    benchmark:
        "workflow/benchmarks/htseq_expression/{sample}.tsv"
    conda:
        "../envs/htseq.yaml"
    shell:
        """
        htseq-count -r pos -s yes -t exon -i gene_id --additional-attr gene_name \
        /path/to/output/{sample}_Aligned.toTranscriptome.out.bam \
        /path/to/genome/annotation.gtf > gene_counts.txt 2>> {log}
        """