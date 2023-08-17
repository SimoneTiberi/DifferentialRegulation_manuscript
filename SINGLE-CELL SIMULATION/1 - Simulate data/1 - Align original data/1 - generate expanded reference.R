alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation

R
# header: explain code

# load libraries
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biostrings)
  library(BSgenome)
  library(eisaR)
  library(GenomicFeatures)
  library(SummarizedExperiment)
  library(tximeta)
  library(rjson)
  library(reticulate)
  library(SingleCellExperiment)
  library(scater)
})

# set global variables
options(timeout = 300)
path <- "01_annotation/"
fa_file <- "GRCm38.primary_assembly.genome.fa"
gtf_file <- "gencode.vM24.annotation.gtf.gz"

# download data
if (!file.exists(paste0(path, fa_file))) {
  #download.file(url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/GRCm38.primary_assembly.genome.fa.gz",
  #              destfile = paste0(path, fa_file, ".gz"))
  #gunzip(paste0(path, fa_file, ".gz"))
  print("NOT THERE!")
}

if (!file.exists(paste0(path, gtf_file))) {
  #download.file(url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gtf.gz",
  #              destfile = paste0(path, gtf_file))
  print("NOT THERE!")
}

# prepare data
gtf <- paste0(path, gtf_file)
grl <- eisaR::getFeatureRanges(
  gtf = gtf,
  featureType = c("spliced", "intron"),
  intronType = "separate",
  flankLength = 90L, # what read length to use?
  joinOverlappingIntrons = FALSE,
  verbose = TRUE
)
#Warning message:
#  In .get_cds_IDX(mcols0$type, mcols0$phase) :
#  The "phase" metadata column contains non-NA values for features of type
#stop_codon. This information was ignored.

genome <- Biostrings::readDNAStringSet(
  paste0(path, fa_file)
)
#Error in x$.self$finalize() : attempt to apply non-function
#Warning message:
#  call dbDisconnect() when finished working with a connection

names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = grl
)

Biostrings::writeXStringSet(
  seqs, filepath = paste0(path, "gencode.vM24.annotation.expanded.fa")
)

eisaR::exportToGtf(
  grl, 
  filepath = paste0(path, "gencode.vM24.annotation.expanded.gtf")
)

write.table(
  metadata(grl)$corrgene, 
  file = paste0(path, "gencode.vM24.annotation.expanded.features.tsv"),
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)

df <- eisaR::getTx2Gene(
  grl, filepath = paste0(path, "gencode.vM24.annotation.expanded.tx2gene.tsv")
)

