rule fastqc_raw:
    input:
    "raw_data/samples/Exp_2_3_1.fastq"

    output:
    "QC_raw/Exp_2_3_1.zip"



rule multiqc_raw:




# Variables
INPUT="raw_data/RNAseq_repeat_fastq/*.fastq.gz"
OUTPUT="read_quality/rawReads"

# running FastQC
fastqc -t 4\
    -o $OUTPUT \
    -d $OUTPUT \
$INPUT

multiqc $OUTPUT -o $OUTPUT

mamba deactivate