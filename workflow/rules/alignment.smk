# ALIGNMENT USING HISAT2
# Index the reference genome

rule hisat2_index:
    input:
        genome=config["genome"]["genome_file"]
    output:
        index=expand("results/alignment/index/{{genome}}.{i}.ht2", i=range(1,9))
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
        hisat2-build -p {threads} {input.genome} {params.outindex} &>> {log}
        """

# Align sample reads to indexed genome
# Outputs a BAM and BAM.BAI index file

rule hisat2_align:
    input:
        fwfasta="results/preprocessed/{sample}_R1.trimmed.fastq.gz",
        rvfasta="results/preprocessed/{sample}_R2.trimmed.fastq.gz",
        index=expand("results/alignment/index/{genome}.{i}.ht2", genome=config["genome"]["genome_name"], i=range(1,9))
    output:
        bam="results/alignment/{sample}.bam",
        bambai="results/alignment/{sample}.bam.bai"
    params:
        strandness="RF",
        index=f"results/alignment/index/{config['genome']['genome_name']}"
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
        -x {params.index} -1 {input.fwfasta} -2 {input.rvfasta} | samtools view -bhS | samtools sort -o {output.bam}
        sambamba index {output.bam}) &>> {log}
        """
# ALIGNMENT USING STAR
# Index the reference genome

rule star_index:
    input:
        genome=config["genome"]["genome_file"],
        gtf=config["genome"]["annotation_file"]
    output:
        index=expand("results/alignment/star_index/{{genome}}/{file}", file=["Genome", "SA", "SAindex"])
    params:
        overhang=100,
        outdir="results/alignment/star_index/{genome}",
        SAindexNbases=12 #for a 100Mb genome, 2.2.5 manual
    threads: 16
    log:
        "workflow/logs/star_index/{genome}.log"
    benchmark:
        repeat("workflow/benchmarks/star_index/{genome}.tsv", 3)
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
# Outputs a transcriptome BAM and sorted BAM

rule star_align:
    input:
        fwfasta="results/preprocessed/{sample}_R1.trimmed.fastq.gz",
        rvfasta="results/preprocessed/{sample}_R2.trimmed.fastq.gz",
        index=f"results/alignment/star_index/{config['genome']['genome_name']}/Genome"
    output:
        trbam="results/alignment/star/{sample}_Aligned.toTranscriptome.out.bam",
        coordbam="results/alignment/star/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        index=f"results/alignment/star_index/{config['genome']['genome_name']}",
        outprefix="results/alignment/star/{sample}_"
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
            --quantMode TranscriptomeSAM \
            --genomeLoad NoSharedMemory &>> {log}
        """