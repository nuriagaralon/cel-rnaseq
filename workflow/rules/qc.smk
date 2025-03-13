rule fastqc_raw:
    input:
        "raw_data/samples/{sample_pr}.fastq.gz"
    output:
        "results/read_quality/QC_raw/fastqc/{sample_pr}_fastqc.zip"
    params:
        outdir="results/read_quality/QC_raw/fastqc"
    threads: 4
    log:
        "workflow/logs/fastqc_raw/{sample_pr}.log"
    benchmark: 
        "workflow/benchmarks/fastqc_raw/{sample_pr}.tsv"
    conda:
        "envs/qc.yaml"
    shell:
        """
        fastqc -t {threads} -o {params.outdir} --dir {params.outdir} {input}
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
        "workflow/benchmarks/multiqc_raw/multiqc.tsv"
    conda:
        "envs/qc.yaml"
    shell:
        """
        FILENAME=$(basename {output})
        multiqc {params.outdir} -n $FILENAME -o {params.outdir}
        """

rule fastqc_trimmed:
    input:
        "results/preprocessed/{sample}_{type}.trimmed.fastq.gz"
    output:
        "results/read_quality/QC_trimmed/fastqc/{sample}_{type}_fastqc.zip"
    params:
        outdir="results/read_quality/QC_trimmed/fastqc"
    threads: 4
    log:
        "workflow/logs/fastqc_trimmed/{sample}_{type}.log"
    benchmark: 
        "workflow/benchmarks/fastqc_trimmed/{sample}_{type}.tsv"
    conda:
        "envs/qc.yaml"
    shell:
        """
        fastqc -t {threads} -o {params.outdir} --dir {params.outdir} {input}
        """

MULTIFR = ["R1", "R2", "SE"]

rule multiqc_trimmed:
    input:
        expand("results/read_quality/QC_trimmed/fastqc/{sample}_{pr}_fastqc.zip", sample=config['samples'], pr=MULTIFR) 
    output:
        "results/read_quality/QC_trimmed/qcreport_trimmed.html"
    params:
        outdir="results/read_quality/QC_trimmed"
    threads: 1
    log:
        "workflow/logs/multiqc_trimmed/multiqc.log"
    benchmark:
        "workflow/benchmarks/multiqc_trimmed/multiqc.tsv"
    conda:
        "envs/qc.yaml"
    shell:
        """
        FILENAME=$(basename {output})
        multiqc {params.outdir} -n $FILENAME -o {params.outdir}
        """