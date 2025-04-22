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
        repeat("workflow/benchmarks/hisat2_index/{genome}.tsv", 3)
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
        repeat("workflow/benchmarks/hisat2_align/{sample}.tsv", 3)
    conda:
        "../envs/hisat.yaml"
    shell:
        """
        (hisat2 -p {threads} --dta --rna-strandness {params.strandness} \
        -x {params.index} -1 {input[0]} -2 {input[1]} | samtools view -bhS | samtools sort -o {output.bam}
        sambamba index {output.bam}) &>> {log}
        """
#--summary-file {output.sum} --met-file {output.met}


rule star_index:
    input:
        config["genome"]["genome_file"],
        config["genome"]["annotation_file"]
    output:
        "results/alignment/star_index/Genome",
        "results/alignment/star_index/"
    params:
        genome=config["genome"]["genome_name"],
        overhang=100,
        outdir="results/alignment/star_index",
        SAindexNbases=12 #for a 100Mb genome, 2.2.5 manual
    threads: 16
    log:
        "workflow/logs/star_index/{params.genome}.log"
    benchmark:
        repeat("workflow/benchmarks/star_index/{params.genome}.tsv", 3)
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR --runMode genomeGenerate \
            --genomeFastaFiles {input[0]} \
            --sjdbGTFfile  {input[1]} \
            --sjdbOverhang {params.overhang} \
            --runThreadN {threads} \
            --genomeSAindexNbases {params.SAindexNbases} \
            --genomeDir {params.outdir} &>> {log}
        """


rule star_align:
    input:
        "results/preprocessed/{sample}_R1.trimmed.fastq.gz",
        "results/preprocessed/{sample}_R2.trimmed.fastq.gz",
        "results/alignment/star_index/Genome"
    output:
        "results/alignment/star/{sample}_Aligned.toTranscriptome.out.bam",
        "results/alignment/star/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        index="results/alignment/star_index/",
        outpref="results/alignment/star/{sample}_"
    threads: 8
    log:
        "workflow/logs/star_align/{sample}.log"
    benchmark:
        repeat("workflow/benchmarks/star_align/{sample}.tsv", 3)
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR --genomeDir {params.index} \
            --readFilesIn {input[0]} {input[1]} \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFileNamePrefix {params.outpref} \
            --outSAMunmapped Within \
            --outSAMattributes NH HI AS nM \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode TranscriptomeSAM \
            --genomeLoad NoSharedMemory &>> {log}
        """