alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation

R

library(eisaR)
library(Biostrings)
library(BSgenome)
library(stringr)
library(GenomicFeatures)
library(roe)

# extra_spliced and extra_unspliced ??
# read_length and flank_trim_length -> default values in minnow repo.

filename_prefix = "transcriptome_splici"
output_dir = tempdir()

gtf_path <- file.path("01_annotation/gencode.vM24.annotation.gtf.gz")
genome_path = file.path("01_annotation/GRCm38.primary_assembly.genome.fa")
read_length = 101
flank_trim_length = 5 # in theory not in minnow data
output_dir = file.path("05_alevin")

make_splici_txome(
  genome_path = genome_path,
  gtf_path = gtf_path,
  read_length = read_length,
  output_dir = output_dir,
  flank_trim_length = flank_trim_length,
  filename_prefix = filename_prefix)

grep("transcriptome_splici", dir(output_dir), value = TRUE)
