# get reference files:
cd /home/Shared_sherborne/simone/Diff_Velo/discovery_bulk

# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf.gz
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/GRCm39.primary_assembly.genome.fa.gz

alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'
R

genome = "reference/GRCm39.primary_assembly.genome.fa.gz"
gtf = "reference/gencode.vM33.annotation.gtf.gz"
outdir = "reference"

suppressPackageStartupMessages({
  library(eisaR)
  library(Biostrings)
  library(GenomicFeatures)
  library(BSgenome)
})

bnm <- sub("\\.gtf", "", basename(gtf))

## Extract a GRanges object with spliced and unspliced transcripts
grl <- eisaR::getFeatureRanges(
  gtf = gtf,
  featureType = c("spliced", "unspliced"), 
  verbose = TRUE
)

## Extract spliced and unspliced transcript sequences
genome <- Biostrings::readDNAStringSet(genome)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = grl
)

## Remove duplicated sequences (mostly transcripts without introns), from both 
## transcript sequences and GRangesList object
tx_to_keep <- names(seqs)[which(!duplicated(seqs))]
grl <- grl[tx_to_keep]
seqs <- seqs[tx_to_keep]

## Additionally remove the few unspliced transcripts where the spliced 
## variant has been removed
(only_unspliced <- paste0(setdiff(gsub("-U", "", names(seqs)), names(seqs)), "-U"))
tx_to_keep <- names(seqs)[!(names(seqs) %in% only_unspliced)]
grl <- grl[tx_to_keep]
seqs <- seqs[tx_to_keep]

stopifnot(all(!duplicated(seqs)))

## We don't want to quantify spliced and unspliced sequences as different genes - 
## remove the "-U" suffix from the 'unspliced' gene IDs
## Use the information in the "corrgene" data frame
head(metadata(grl)$corrgene, 3)
## Unlist to make the whole thing run faster
g <- unlist(grl, use.names = FALSE)
## First check that relisting gives the same thing as before (minus the metadata)
tmp <- grl
metadata(tmp) <- list()
identical(relist(g, grl), tmp)
## Then replace the gene IDs and relist
idx <- which(g$gene_id %in% metadata(grl)$corrgene$unspliced)
g$gene_id[idx] <- metadata(grl)$corrgene$spliced[match(g$gene_id[idx], metadata(grl)$corrgene$unspliced)]
grlr <- relist(g, grl)  ## metadata is removed

stopifnot(all(names(seqs) == names(grlr)))

## Export
df <- eisaR::getTx2Gene(
  grlr, filepath = file.path(outdir, paste0(bnm, ".expanded.tx2gene.tsv"))
)
Biostrings::writeXStringSet(
  seqs, filepath = file.path(outdir, paste0(bnm, ".expanded.transcripts.fa"))
)
eisaR::exportToGtf(
  grlr, filepath = file.path(outdir, paste0(bnm, ".expanded.gtf"))
)
