cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation

#!/bin/bash
salmon="../software/salmon-1.5.2/bin/salmon"

#############################################################################
# Build index
#############################################################################

#!/bin/bash
grep ">" 01_annotation/GRCm38.primary_assembly.genome.fa | cut -d ">" -f 2 | cut -d " " -f 1 > 01_annotation/GRCm38.primary_assembly.genome.chrnames.txt

$salmon index \
-t <(cat 01_annotation/gencode.vM24.annotation.expanded.fa 01_annotation/GRCm38.primary_assembly.genome.fa) \
-i 01_annotation/gencode.vM24.annotation.expanded.sidx --gencode -p 32 \
-d 01_annotation/GRCm38.primary_assembly.genome.chrnames.txt

alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

R -e 'tximeta::makeLinkedTxome(
  indexDir = "01_annotation/gencode.vM24.annotation.expanded.sidx", 
  source = "GENCODE", genome = "GRCm38", 
  organism = "Mus musculus", release = "M24", 
  fasta = "01_annotation/gencode.vM24.annotation.expanded.fa", 
  gtf = "01_annotation/gencode.vM24.annotation.expanded.gtf", 
  write = TRUE, jsonFile = "01_annotation/gencode.vM24.annotation.expanded.json"
)'

#############################################################################
# Align and quantify
#############################################################################

## NORMAL1
$salmon alevin -l ISR -i 02_alevin/index \
-1 fastq/SRR6337197_1.fastq fastq/SRR6337198_1.fastq fastq/SRR6337199_1.fastq \
fastq/SRR6337200_1.fastq fastq/SRR6337201_1.fastq fastq/SRR6337202_1.fastq \
fastq/SRR6337203_1.fastq fastq/SRR6337204_1.fastq \
-2 fastq/SRR6337197_2.fastq fastq/SRR6337198_2.fastq fastq/SRR6337199_2.fastq \
fastq/SRR6337200_2.fastq fastq/SRR6337201_2.fastq fastq/SRR6337202_2.fastq \
fastq/SRR6337203_2.fastq fastq/SRR6337204_2.fastq \
-o 02_alevin/normal1/ -p 32 --chromium --dumpFeatures --dumpBfh \
--tgMap 01_annotation/gencode.vM24.annotation.expanded.tx2gene.tsv

## NORMAL2
$salmon alevin -l ISR -i 02_alevin/index \
-1 fastq/SRR6337205_1.fastq fastq/SRR6337206_1.fastq fastq/SRR6337207_1.fastq \
fastq/SRR6337208_1.fastq fastq/SRR6337209_1.fastq fastq/SRR6337210_1.fastq \
fastq/SRR6337211_1.fastq fastq/SRR6337212_1.fastq \
-2 fastq/SRR6337205_2.fastq fastq/SRR6337206_2.fastq fastq/SRR6337207_2.fastq \
fastq/SRR6337208_2.fastq fastq/SRR6337209_2.fastq fastq/SRR6337210_2.fastq \
fastq/SRR6337211_2.fastq fastq/SRR6337212_2.fastq \
-o 02_alevin/normal2/ -p 32 --chromium --dumpFeatures --dumpBfh \
--tgMap 01_annotation/gencode.vM24.annotation.expanded.tx2gene.tsv

## NORMAL3
$salmon alevin -l ISR -i 02_alevin/index \
-1 fastq/SRR6337213_1.fastq fastq/SRR6337214_1.fastq fastq/SRR6337215_1.fastq \
fastq/SRR6337216_1.fastq fastq/SRR6337217_1.fastq fastq/SRR6337218_1.fastq \
fastq/SRR6337219_1.fastq fastq/SRR6337220_1.fastq \
-2 fastq/SRR6337213_2.fastq fastq/SRR6337214_2.fastq fastq/SRR6337215_2.fastq \
fastq/SRR6337216_2.fastq fastq/SRR6337217_2.fastq fastq/SRR6337218_2.fastq \
fastq/SRR6337219_2.fastq fastq/SRR6337220_2.fastq \
-o 02_alevin/normal3/ -p 32 --chromium --dumpFeatures --dumpBfh \
--tgMap 01_annotation/gencode.vM24.annotation.expanded.tx2gene.tsv

## NORMAL4
$salmon alevin -l ISR -i 02_alevin/index \
-1 fastq/SRR6337221_1.fastq fastq/SRR6337222_1.fastq fastq/SRR6337223_1.fastq \
fastq/SRR6337224_1.fastq fastq/SRR6337225_1.fastq fastq/SRR6337226_1.fastq \
fastq/SRR6337227_1.fastq fastq/SRR6337228_1.fastq \
-2 fastq/SRR6337221_2.fastq fastq/SRR6337222_2.fastq fastq/SRR6337223_2.fastq \
fastq/SRR6337224_2.fastq fastq/SRR6337225_2.fastq fastq/SRR6337226_2.fastq \
fastq/SRR6337227_2.fastq fastq/SRR6337228_2.fastq \
-o 02_alevin/normal4/ -p 32 --chromium --dumpFeatures --dumpBfh \
--tgMap 01_annotation/gencode.vM24.annotation.expanded.tx2gene.tsv
