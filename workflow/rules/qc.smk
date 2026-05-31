# QC USING FASTQC AND MULTIQC
# Run FastQC on all raw fasta files
# Aggregate stats results from alignment and expression

rule fastqc_raw:
    input:
        "raw_data/samples/{sample_pr}.fastq.gz"
    output:
        "results/read_quality/QC_raw/fastqc/{sample_pr}_fastqc.zip"
    params:
        outdir="results/read_quality/QC_raw/fastqc"
    threads: 1
    log:
        "workflow/logs/fastqc_raw/{sample_pr}.log"
    benchmark: 
        repeat("workflow/benchmarks/fastqc_raw/{sample_pr}.tsv", 3)
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastqc -t {threads} -o {params.outdir} --dir {params.outdir} {input} &>> {log}
        """

# Aggregate FastQC reports in one file

rule multiqc_raw:
    input:
        expand("results/read_quality/QC_raw/fastqc/{sample}_{pr}_fastqc.zip", sample=sample_names, pr=config['pairedreads'])
    output:
        outfile="results/read_quality/QC_raw/qcreport_raw.html",
        outdata=directory("results/read_quality/QC_trimmed/qcreport_raw_data")
    params:
        outdir="results/read_quality/QC_raw"
    threads: 1
    log:
        "workflow/logs/multiqc_raw/multiqc.log"
    benchmark:
        repeat("workflow/benchmarks/multiqc_raw/multiqc.tsv", 3)
    conda:
        "../envs/qc.yaml"
    shell:
        """
        FILENAME=$(basename {output.outfile})
        multiqc {params.outdir} -n $FILENAME -o {params.outdir} &>> {log}
        """

# Run FastQC on all trimmed fasta files

rule fastqc_trimmed:
    input:
        "results/preprocessed/{sample}_{type}.trimmed.fastq.gz"
    output:
        "results/read_quality/QC_trimmed/fastqc/{sample}_{type}.trimmed_fastqc.zip"
    params:
        outdir="results/read_quality/QC_trimmed/fastqc"
    threads: 1
    log:
        "workflow/logs/fastqc_trimmed/{sample}_{type}.log"
    benchmark: 
        repeat("workflow/benchmarks/fastqc_trimmed/{sample}_{type}.tsv", 3)
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastqc -t {threads} -o {params.outdir} --dir {params.outdir} {input} &>> {log}
        """

# Aggregate FastQC reports in one file (include Forward, reverse and unpaired reads)

rule multiqc_trimmed:
    input:
        expand("results/read_quality/QC_trimmed/fastqc/{sample}_{pr}.trimmed_fastqc.zip", sample=sample_names, pr=config['trimqc']) 
    output:
        outfile="results/read_quality/QC_trimmed/qcreport_trimmed.html",
        outdata=directory("results/read_quality/QC_trimmed/qcreport_trimmed_data")
    params:
        outdir="results/read_quality/QC_trimmed"
    threads: 1
    log:
        "workflow/logs/multiqc_trimmed/multiqc.log"
    benchmark:
        repeat("workflow/benchmarks/multiqc_trimmed/multiqc.tsv", 3)
    conda:
        "../envs/qc.yaml"
    shell:
        """
        FILENAME=$(basename {output.outfile})
        multiqc {params.outdir} -n $FILENAME -o {params.outdir} &>> {log}
        """

# Aggregate STAR quality logs

rule aggregate_star:
    input:
        stats=expand("results/alignment/star/{sample}_Log.final.out", sample=config["samples"])
    output:
        stardata="results/read_quality/star_qc.tsv"
    threads: 1
    log:
        "workflow/logs/aggregate_star.log"
    benchmark:
        "workflow/benchmarks/aggregate_star.tsv"  
    conda:
        "../envs/Raggregate.yaml"
    script:
        "../scripts/star_qc.R"

# Aggregate featureCounts summaries

rule aggregate_fc:
    input:
        stats=expand("results/expression/featurecounts/{sample}_gene_counts.tsv.summary", sample=config["samples"])
    output:
        fcdata="results/read_quality/featurecounts_qc.tsv"
    threads: 1
    log:
        "workflow/logs/aggregate_fc.log"
    benchmark:
        "workflow/benchmarks/aggregate_fc.tsv"  
    conda:
        "../envs/Raggregate.yaml"
    script:
        "../scripts/featurecounts_qc.R"

