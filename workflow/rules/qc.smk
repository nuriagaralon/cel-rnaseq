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


rule multiqc_raw:
    input:
        expand("results/read_quality/QC_raw/fastqc/{sample}_{pr}_fastqc.zip", sample=config['samples'], pr=config['pairedreads'])
    output:
        "results/read_quality/QC_raw/qcreport_raw.html"
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
        FILENAME=$(basename {output})
        multiqc {params.outdir} -n $FILENAME -o {params.outdir} &>> {log}
        """

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

MULTIFR = ["R1", "R2", "SE"]

rule multiqc_trimmed:
    input:
        expand("results/read_quality/QC_trimmed/fastqc/{sample}_{pr}.trimmed_fastqc.zip", sample=config['samples'], pr=MULTIFR) 
    output:
        "results/read_quality/QC_trimmed/qcreport_trimmed.html"
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
        FILENAME=$(basename {output})
        multiqc {params.outdir} -n $FILENAME -o {params.outdir} &>> {log}
        """