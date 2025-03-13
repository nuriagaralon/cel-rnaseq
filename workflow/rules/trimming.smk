rule trimmomatic_trim:
    input:
        expand("raw_data/samples/{{sample}}_{pr}.fastq.gz", pr=config['pairedreads'])
    output:
        "results/preprocessed/{sample}_R1.trimmed.fastq.gz",
        temp("results/preprocessed/{sample}_SE1.trimmed.fastq.gz"),
        "results/preprocessed/{sample}_R2.trimmed.fastq.gz",
        temp("results/preprocessed/{sample}_SE2.trimmed.fastq.gz")
    params:
        adapter=config['adapters'],
        seed_mismatch=2,
        palindrome_clip=30,
        simple_clip=10,
        adapter_length=5,
        crop_lowq="true",
        leading=3,
        trailing=3,
        window_size=4,
        window_quality=15,
        minlen=36
    threads: 4
    log:
        "workflow/logs/trimmomatic_trim/{sample}.log"
    benchmark:
        "workflow/benchmarks/trimmomatic_trim/{sample}.tsv"
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        """
        trimmomatic PE -threads {threads} {input} {output} \
        ILLUMINACLIP:{params.adapter}:{params.seed_mismatch}:{params.palindrome_clip}:{params.simple_clip}:{params.adapter_length}:{params.crop_lowq} \
        LEADING:{params.leading} \
        TRAILING:{params.trailing} \
        SLIDINGWINDOW:{params.window_size}:{params.window_quality} \
        MINLEN:{params.minlen} &>> {log}
        """

rule trim_join_SE:
    input:
        "results/preprocessed/{sample}_SE1.trimmed.fastq.gz",
        "results/preprocessed/{sample}_SE2.trimmed.fastq.gz"
    output:
        "results/preprocessed/{sample}_SE.trimmed.fastq.gz"
    threads: 1
    log:
        "workflow/logs/trim_join_SE/{sample}.log"
    benchmark:
        "workflow/benchmarks/trim_join_SE/{sample}.tsv"    
    shell:
        "cat {input[0]} {input[1]} > {output} 2>> {log}"
