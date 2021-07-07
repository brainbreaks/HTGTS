library(stringr)
library(dplyr)
library(readr)

macs_cols = cols(
  macs_chrom=col_character(), macs_start=col_double(), macs_end=col_double(), macs_length=col_character(), macs_summit_abs=col_double(),
  macs_pileup=col_double(), macs_pvalue=col_double(), macs_fc=col_double(), macs_qvalue=col_double(), macs_name=col_character(), macs_comment=col_character()
)

tlx_cols = cols(
  Qname=readr::col_character(), JuncID=readr::col_character(), Rname=readr::col_character(), Junction=readr::col_double(),
  Strand=readr::col_character(), Rstart=readr::col_double(), Rend=readr::col_double(),
  B_Rname=readr::col_character(), B_Rstart=readr::col_double(), B_Rend=readr::col_double(), B_Strand=readr::col_double(),
  B_Qstart=readr::col_double(), B_Qend=readr::col_double(), Qstart=readr::col_double(), Qend=readr::col_double(), Qlen=readr::col_double(),
  B_Cigar=readr::col_character(), Cigar=readr::col_character(), Seq=readr::col_character(), J_Seq=readr::col_character(), Barcode=readr::col_logical(),
  unaligned=readr::col_double(), baitonly=readr::col_double(), uncut=readr::col_double(), misprimed=readr::col_double(), freqcut=readr::col_double(),
  largegap=readr::col_double(), mapqual=readr::col_double(), breaksite=readr::col_double(), sequential=readr::col_double(), repeatseq=readr::col_double(), duplicate=readr::col_double()
)

repeatmasker_cols = cols(
  repeatmasker_bin=readr::col_double(),
  repeatmasker_score=readr::col_double(),
  repeatmasker_mismatches_per_kb=readr::col_double(),
  repeatmasker_deletions_per_kb=readr::col_double(),
  repeatmasker_insertions_per_kb=readr::col_double(),
  repeatmasker_chrom=readr::col_character(),
  repeatmasker_start=readr::col_double(),
  repeatmasker_end=readr::col_double(),
  repeatmasker_genoLeft=readr::col_double(),
  repeatmasker_strand=readr::col_character(),
  repeatmasker_name=readr::col_character(),
  repeatmasker_class=readr::col_character(),
  repeatmasker_family=readr::col_character(),
  repeatmasker_repStart=readr::col_character(),
  repeatmasker_repEnd=readr::col_character(),
  repeatmasker_repLeft=readr::col_character(),
  repeatmasker_id=readr::col_character()
)

macs2 = function(name, sample, control=NULL, qvalue=0.01, extsize=200, slocal=1000, output_dir="data/macs2", llocal=10000000) {
  bed_sample = paste("-t", sample)
  bed_control = ifelse(is.null(control), "", paste("-c", control))

  cmd = stringr::str_glue("macs2 callpeak {bed_sample} {bed_control} --seed 123 -f BED -g hs --keep-dup all -n {name} --outdir {output_dir} --nomodel --slocal {slocal} --extsize {extsize} -q {qvalue} --llocal {llocal} --bdg --trackline", bed_sample=bed_sample, bed_control=bed_control, name=name, output_dir=output_dir, extsize=extsize, qvalue=qvalue, llocal=sprintf("%0.0f", llocal), slocal=sprintf("%0.0f", slocal))
  print(cmd)
  system(cmd)

  readr::read_tsv(paste0(output_dir, "/", name, "_peaks.xls"), comment="#", col_names=names(macs_cols$cols), col_types=macs_cols) %>%
    dplyr::slice(-1) %>%
    dplyr::select(-macs_comment)
}