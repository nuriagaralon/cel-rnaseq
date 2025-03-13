#Besides stringtie, check htseq and stuff helena gave to me
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
        "envs/stringtie.yaml"
    shell:
        """
        stringtie {input} -p {threads} --rf -l {wildcards.sample} \
        -o {output[0]} -G {params.refgen} -A {output[1]} -eB
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
        "envs/stringtie.yaml"
    shell:
        """
        gffcompare -R -r {params.refgen} {input} -o {params.prefix}
        """
