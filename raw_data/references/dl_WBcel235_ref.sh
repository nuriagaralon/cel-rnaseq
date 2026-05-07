#!/bin/bash
# SCRIPT TO DOWNLOAD, CHECK AND CLEAN REFERENCE DATA
set -e

# Download C. elegans WBcel235 reference files
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_rna.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/md5checksums.txt

# Checksums
grep "WBcel235_genomic.fna\|WBcel235_genomic.gtf\|WBcel235_rna.fna"  md5checksums.txt > selected_checksums.txt
md5sum -c selected_checksums.txt || {
    echo "ERROR: Checksum verification failed. Please download again."
    rm -f GCF_000002985.6_WBcel235_genomic.fna.gz
    rm -f GCF_000002985.6_WBcel235_genomic.gtf.gz
    rm -f GCF_000002985.6_WBcel235_rna.fna.gz
    rm -f md5checksums.txt 
    rm -f selected_checksums.txt
    exit 1
}

# Remove checksum files
rm -f md5checksums.txt 
rm -f selected_checksums.txt

# Unzip
gunzip *.gz

# Update gtf file (Code generated with ChatGPT)
awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} !($3 ~ /^(gene|start_codon|stop_codon)$/) {
n=$9
gsub(/db_xref "[^"]*"; ?/,"",n)
gsub(/locus_tag "[^"]*"; ?/,"",n)
gsub(/product "[^"]*"; ?/,"",n)
gsub(/standard_name "[^"]*"; ?/,"",n)
gsub(/transcript_biotype "[^"]*"; ?/,"",n)
gsub(/exon_number "[^"]*"; ?/,"",n)
gsub(/note "[^"]*"; ?/,"",n)

match(n,/transcript_id "([^"]+)"/,t)
match(n,/gene_id "([^"]+)"/,g)
match(n,/gene "([^"]+)"/,gn)

$9="transcript_id \""t[1]"\"; gene_id \""g[1]"\"; gene_name \""gn[1]"\";"

print
}' GCF_000002985.6_WBcel235_genomic.gtf > temp.gtf


# Remove unassigned transcripts (required for RSEM)
awk 'BEGIN{FS=OFS="\t"}
/^#/ {print; next}

$9 !~ /transcript_id "unassigned_transcript_/ {
    print
}' temp.gtf > GCF_000002985.6_WBcel235_clean.gtf

rm temp.gtf