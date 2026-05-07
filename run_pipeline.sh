#!/bin/bash

# SCRIPT TO HELP RUN THE PIPELINE 
# Activate conda
conda activate snakemake

# Create DAG and rulegraph:
snakemake all --dag | dot -T svg > images/dag.svg

snakemake --rulegraph | dot -T svg > images/rulegraph.svg

# Dry run
snakemake -np all

# Run pipeline
snakemake -p --use-conda all
#snakemake -p --use-conda --cores 8 all