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
        adapter_length=config["trimming"]["adaptlen"],
        crop_lowq="true",
        leading=3,
        trailing=3,
        window_size=4,
        window_quality=config["trimming"]["quality"],
        minlen=config["trimming"]["minlen"]
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

def get_join_input(wildcards):
    if config["tools"]["trim"] == "trimmomatic":
        return ["results/preprocessed/{wildcards.sample}_SE1.trimmed.fastq.gz",
                "results/preprocessed/{wildcards.sample}_SE2.trimmed.fastq.gz"]
    elif config["tools"]["trim"] == "trimgalore":
        return ["results/preprocessed/{wildcards.sample}_unpaired_1.fq.gz",
                "results/preprocessed/{wildcards.sample}_unpaired_2.fq.gz"]
    else:
        raise ValueError("Unknown trimming tool specified in the config file.")

rule trim_join_SE:
    input:
        get_join_input
    output:
        "results/preprocessed/{sample}_SE.trimmed.fastq.gz"
    threads: 1
    log:
        "workflow/logs/trim_join_SE/{sample}.log"
    benchmark:
        "workflow/benchmarks/trim_join_SE/{sample}.tsv"    
    shell:
        "cat {input[0]} {input[1]} > {output} 2>> {log}"

rule trimgalore_trim:
    input:
        expand("raw_data/samples/{{sample}}_{pr}.fastq.gz", pr=config['pairedreads'])
    output:
        "results/preprocessed/{sample}_val_1.fq.gz",
        temp("results/preprocessed/{sample}_unpaired_1.fq.gz"),
        "results/preprocessed/{sample}_val_2.fq.gz",
        temp("results/preprocessed/{sample}_unpaired_2.fq.gz")
    params:
        adapters="--illumina",
        paired="--paired --retain_unpaired",
        quality=config["trimming"]["quality"],
        stringency=config["trimming"]["adaptlen"],
        minlen=config["trimming"]["minlen"],
        outdir="results/preprocessed"
    threads: 1
# Cores are odd for trimgalore: check --help -j/--cores
    log:
        "workflow/logs/trimgalore_trim/{sample}.log"
    benchmark:
        "workflow/benchmarks/trimgalore_trim/{sample}.tsv"
    conda:
        "../envs/trimgalore.yaml"
    shell:
        """
        trim_galore {params.adapters} {params.paired}\
        -q {params.quality} --stringency {params.stringency} --length {params.minlen} \
        --cores {threads} -o {params.outdir} {input[0]} {input[1]} &>> {log}
        """
