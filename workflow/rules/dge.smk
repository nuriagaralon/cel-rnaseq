# DGE DATA AGGREGATION
# Gathers count tables to one R data file

rule aggregate_dge:
    input:
        counts=expand("results/expression/{sample}_gene_counts.tsv", sample=config["samples"])
        features=config["genome"]["features_file"]
        metadata=config["metadata"]
    output:
        rdata="results/dge/gene_counts_data.RData"
    threads: 1
    log:
        "workflow/logs/aggregate_dge.log"
    benchmark:
        "workflow/benchmarks/aggregate_dge.tsv"  
    conda:
        "../envs/Raggregate.yaml"
    script:
        "workflow/scripts/dge_aggregate.R"
        