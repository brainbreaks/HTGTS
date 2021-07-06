library(readr)
library(stringr)
library(dplyr)
library(GenomeInfoDb)


options(shiny.maxRequestSize=5*1024^3)
options(shiny.sanitize.errors = TRUE)

# baits_df = readr::read_tsv("data/baits.tsv")
# samples_df = readr::read_tsv("data/samples.tsv")
# offtargets_df = readr::read_tsv("data/offtargets.tsv")
repeatmasker_df = readr::read_tsv("data/ucsc_hg19_repeatmasker.tsv", col_names=names(repeatmasker_cols$cols), col_types=repeatmasker_cols, skip=1) %>%
  dplyr::mutate(repeatmasker_id=1:n())
repeatmasker_ranges = GenomicRanges::makeGRangesFromDataFrame(repeatmasker_df %>% dplyr::mutate(seqnames=repeatmasker_chrom, start=repeatmasker_start, end=repeatmasker_end), keep.extra.columns=T)

#
#
# default_genome <- "mm9"
# default_chrom <- "chr1"
# chromsizes_cols <- readr::cols(seqnames=col_character(), seqlengths=col_double())
# bed_cols <- cols(seqnames=readr::col_character(), start=readr::col_double(), end=readr::col_double(), name=readr::col_character(), score=readr::col_character(), strand=readr::col_character())
# genome_info_df <- readr::read_tsv(stringr::str_glue("data/{id}/{id}.chrom.sizes", id=default_genome), col_names=names(chromsizes_cols$cols), col_types=chromsizes_cols) %>%
#     dplyr::filter(!grepl(".*_.*", seqnames))
# genome_info <- with(genome_info_df, GenomeInfoDb::Seqinfo(seqnames, seqlengths, isCircular=rep(F, length(seqnames)), genome=rep(default_genome, length(seqnames))))
# default_chrom_choices <- genome_info_df$seqnames
# default_range_max <- genome_info_df %>% dplyr::filter(seqnames==default_chrom) %>% .$seqlengths
# genome_txdb <- load_txdb(stringr::str_glue('data/{id}/{id}.refGene.gtf.gz', id=default_genome))
# genome_genes <- txdb2genes(genome_txdb$longest_transcript)
