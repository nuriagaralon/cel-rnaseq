# ALIGNMENT USING STAR
# Index the reference genome

rule star_index:
    input:
        genome=config["genome"]["genome_file"],
        gtf=config["genome"]["annotation_file"]
    output:
        index=expand("results/alignment/index/{{genome}}/{file}", file=["Genome", "SA", "SAindex"])
    params:
        overhang=100,
        outdir="results/alignment/index/{genome}",
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
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile  {input.gtf} \
            --sjdbOverhang {params.overhang} \
            --runThreadN {threads} \
            --genomeSAindexNbases {params.SAindexNbases} \
            --genomeDir {params.outdir} &>> {log}
        """
# Align sample reads to indexed genome
# Outputs a sorted BAM

rule star_align:
    input:
        fwfasta="results/preprocessed/{sample}_R1.trimmed.fastq.gz",
        rvfasta="results/preprocessed/{sample}_R2.trimmed.fastq.gz",
        index=f"results/alignment/index/{config['genome']['genome_name']}/Genome"
    output:
        bam="results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        index=f"results/alignment/index/{config['genome']['genome_name']}",
        outprefix="results/alignment/{sample}_"
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
            --readFilesIn {input.fwfasta} {input.rvfasta} \
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
            --outFileNamePrefix {params.outprefix} \
            --outSAMunmapped Within \
            --outSAMattributes NH HI AS nM \
            --outSAMtype BAM SortedByCoordinate \
            --genomeLoad NoSharedMemory &>> {log}
        """