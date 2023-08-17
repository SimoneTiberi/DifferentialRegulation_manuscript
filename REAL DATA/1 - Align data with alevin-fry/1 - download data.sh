# Tutorial: How to download raw sequence data from GEO/SRA
# https://www.biostars.org/p/111040/

# SRA IDS: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA531650&o=acc_s%3Aa
# click on “Accession List”

# https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR8869267

# wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR8869267/SRR8869267.1

# Libraries from different samples were pooled based on molar concentrations and sequenced on a NextSeq 500 instrument (Illumina) with 26 bases for read 1, 57 bases for read 2 and 8 bases for index 1. 

cd /home/Shared_sherborne/simone/Diff_Velo/discovery/01_data

fastq_dump="/home/stiber/software/sratoolkit.3.0.6-ubuntu64/bin/fastq-dump"

$fastq_dump -I --split-files SRR8869247
$fastq_dump -I --split-files SRR8869248
$fastq_dump -I --split-files SRR8869249
$fastq_dump -I --split-files SRR8869262
$fastq_dump -I --split-files SRR8869263
$fastq_dump -I --split-files SRR8869264
