#!/bin/bash

snakemake --dag all | dot -T svg > images/general_dag.svg
snakemake --dag all_qc | dot -T svg > images/qc_dag.svg
snakemake --dag all_align | dot -T svg > images/align_dag.svg
snakemake --dag all_expression | dot -T svg > images/expression_dag.svg

snakemake --rulegraph | dot -T svg > images/rulegraph.svg

snakemake -np all
