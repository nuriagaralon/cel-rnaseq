#!/bin/bash

# SCRIPT TO HELP RUN THE PIPELINE 
# Activate conda
conda activate snakemake

# Create DAG and rulegraphs:
snakemake all --dag | dot -T svg > images/all_dag.svg
snakemake qc_bench --dag | dot -T svg > images/qc_dag.svg
snakemake trim_bench --dag | dot -T svg > images/trim_dag.svg
snakemake alignment_bench --dag | dot -T svg > images/alignment_dag.svg
snakemake expression_bench --dag | dot -T svg > images/expression_dag.svg

snakemake --rulegraph | dot -T svg > images/bench_rulegraph.svg

# Dry run
snakemake -np all

# Run pipeline
snakemake -p --use-conda all
#snakemake -p --use-conda --cores 8 all