alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation

R

fa="01_annotation/GRCm38.primary_assembly.genome.fa"
gtf="01_annotation/gencode.vM24.annotation.gtf.gz"
outdir="01_annotation_for_minnow/"

#-----------------------------------------------------------------#
# Load packages
#-----------------------------------------------------------------#
suppressPackageStartupMessages({
  library(eisaR)
  library(tximeta)
  library(BSgenome)
})

#-----------------------------------------------------------------#
# main step
#-----------------------------------------------------------------#
process <- function(gtf, fa, outdir) {
  # setup output file name
  dir.create(file.path(outdir), showWarnings = FALSE)
  out_fa=file.path(outdir,"annotation.expanded.fa")
  out_gtf=file.path(outdir,"annotation.expanded.gtf")
  out_tx2gene=file.path(outdir,"annotation.expanded.tx2gene.tsv")
  
  # Initialization
  featureType = c("spliced", "unspliced")
  suffixes <- c(unspliced = "-U")
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")
  grlfull <- GenomicRanges::GRangesList()
  featurelist <- list()
  ebt <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  t2g <- AnnotationDbi::select(txdb, keys = names(ebt), keytype = "TXNAME", 
                               columns = "GENEID")
  e2 <- BiocGenerics::unlist(ebt)
  e2$transcript_id <- names(e2)
  e2$gene_id = t2g$GENEID[match(e2$transcript_id, t2g$TXNAME)]
  e2$exon_id <- e2$exon_name
  e2$exon_name <- NULL
  e2$type <- "exon"
  names(e2) <- NULL
  S4Vectors::mcols(e2) <- S4Vectors::mcols(e2)[, c("exon_id", "exon_rank", 
                                        "transcript_id", "gene_id", "type")]
  ebt <- BiocGenerics::relist(e2, ebt)
  ebg <- GenomicFeatures::exonsBy(txdb, by = "gene")
  corrtx <- data.frame(spliced = unique(names(ebt)), stringsAsFactors = FALSE)
  corrtx$unspliced <- paste0(corrtx$spliced, suffixes["unspliced"])
  corrgene <- data.frame(spliced = unique(unlist(ebt)$gene_id),
                         stringsAsFactors = FALSE)
  corrgene$unspliced <- paste0(corrgene$spliced, suffixes["unspliced"])
  
  # process spliced txps
  message("Extracting spliced transcript features")
  featurelist$spliced <- names(ebt)
  grlfull <- c(grlfull, ebt)
  
  #process unspliced txps
  message("Extracting unspliced transcript features")
  ebtr <- range(ebg)
  e2 <- BiocGenerics::unlist(ebtr)
  
  e2$exon_rank <- 1L
  e2$transcript_id <- names(e2)
  e2$gene_id = names(e2)
  
  e2$type <- "exon"
  e2$transcript_id <- paste0(e2$transcript_id, suffixes["unspliced"])
  e2$exon_id <- e2$transcript_id
  names(e2) <- NULL
  S4Vectors::mcols(e2) <- S4Vectors::mcols(e2)[, c("exon_id", "exon_rank", 
                                        "transcript_id", "gene_id", "type")]
  ebtr <- BiocGenerics::relist(e2, ebtr)
  names(ebtr) <- paste0(names(ebtr), suffixes["unspliced"])
  featurelist$unspliced <- names(ebtr)
  grlfull <- c(grlfull, ebtr)
  
  # post processing
  S4Vectors::metadata(grlfull)$corrtx <- corrtx[, colnames(corrtx) %in%
                                                  featureType, drop = FALSE]
  S4Vectors::metadata(grlfull)$corrgene <- corrgene[, colnames(corrgene) %in%
                                                      featureType, drop = FALSE]
  S4Vectors::metadata(grlfull)$featurelist <- featurelist
  grl = grlfull
  
  # rm(corrgene,corrtx,e2,ebg,ebt,ebtr,featurelist,grlfull,spliced_txp_length,t2g,unspliced_txp_length)
  
  # writing outputs
  message("Reading genome...")
  
  genome <- Biostrings::readDNAStringSet(
    fa
  )
  names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
  
  message("Extracting transcript sequences...")
  seqs <- GenomicFeatures::extractTranscriptSeqs(
    x = genome, 
    transcripts = grl
  )
  
  message("Writing transcript sequences...")
  Biostrings::writeXStringSet(
    seqs, out_fa
  )

  message("Writing gtf file...")
  eisaR::exportToGtf(
    grl, 
    filepath = out_gtf
  )
  
  message("Writing tx2gene file...")
  df <- eisaR::getTx2Gene(
    grl, filepath = out_tx2gene
    )
  
}

process(gtf=gtf, fa=fa, outdir=outdir)
