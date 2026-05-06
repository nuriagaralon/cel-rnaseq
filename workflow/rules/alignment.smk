rule star_index:
    input:
        config["genome"]["genome_file"],
        config["genome"]["annotation_file"]
    output:
        expand("results/alignment/star_index/{{genome}}/{file}", file=["Genome", "SA", "SAindex"])
    params:
        overhang=100,
        outdir="results/alignment/star_index/{genome}",
        SAindexNbases=12 #for a 100Mb genome, 2.2.5 manual
    threads: 16
    log:
        "workflow/logs/star_index/{genome}.log"
    benchmark:
        "workflow/benchmarks/star_index/{genome}.tsv"
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
        f"results/alignment/star_index/{config['genome']['genome_name']}/Genome"
    output:
        "results/alignment/{sample}_Aligned.toTranscriptome.out.bam",
        "results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        index=f"results/alignment/star_index/{config['genome']['genome_name']}",
        outpref="results/alignment/{sample}_"
    threads: 8
    log:
        "workflow/logs/star_align/{sample}.log"
    benchmark:
        "workflow/benchmarks/star_align/{sample}.tsv"
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