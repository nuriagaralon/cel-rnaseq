rule hisat2_index:
    input:
        config["genome"]["genome_file"]
    output:
        expand("results/alignment/index/{{genome}}.{i}.ht2", i=range(1,9))
    params:
        outindex="results/alignment/index/{genome}"
    threads: 16
    log:
        "workflow/logs/hisat2_index/{genome}.log"
    benchmark:
        "workflow/benchmarks/hisat2_index/{genome}.tsv"
    conda:
        "../envs/hisat.yaml"
    shell:
        """
        hisat2-build -p {threads} {input} {params.outindex} &>> {log}
        """

rule hisat2_align:
    input:
        "results/preprocessed/{sample}_R1.trimmed.fastq.gz",
        "results/preprocessed/{sample}_R2.trimmed.fastq.gz",
        expand("results/alignment/index/{genome}.{i}.ht2", genome=config["genome"]["genome_name"], i=range(1,9))
    output:
        bam="results/alignment/{sample}.bam",
        bambai="results/alignment/{sample}.bam.bai"
    params:
        strandness="RF",
        index=f"results/alignment/index/{config['genome']["genome_name"]}"
    threads: 8
    log:
        "workflow/logs/hisat2_align/{sample}.log"
    benchmark:
        "workflow/benchmarks/hisat2_align/{sample}.tsv"
    conda:
        "../envs/hisat.yaml"
    shell:
        """
        (hisat2 -p {threads} --dta --rna-strandness {params.strandness} \
        -x {params.index} -1 {input[0]} -2 {input[1]} | samtools view -bhS | samtools sort -o {output.bam}
        sambamba index {output.bam}) &>> {log}
        """
#--summary-file {output.sum} --met-file {output.met}
