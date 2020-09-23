#' Parse PAVfinder 0.2.0 output
#'
#' @param events_path Path to the PAVfinder events file
#'
#' @return A data frame with parsed PAVfincer results
#' @export
parse_pavfinder_0.2.0 <- function (events_path) {
  col_names = c("ID", "event", "chrom1", "pos1", "orient1", "chrom2",
                "pos2", "orient2", "size", "contigs", "contig_breaks",
                "contig_support_span", "homol_seq", "homol_coords",
                "homol_len", "novel_sequence", "gene1", "transcript1",
                "exon1",  "exon_bound1", "gene2", "transcript2",
                "exon2", "exon_bound2", "sense_fusion", "5'gene",
                "3'gene", "support_reads", "jn_depth", "ref5_jn_coord",
                "ref5_jn_depth", "ref3_jn_coord", "ref3_jn_depth",
                "transcript1_tpm", "transcript2_tpm", "filter")
  # Start with everything as a character - alter as necessary
  col_types = c("iccccccccccccccccccccccccccddcdcdddc")
  # TODO - rework this so I'm not doing as much counting!
  # E.g.:
  #    struct <- data.frame(
  #      ID = 'i', etc...
  # Read the TSV data, specifying column types
  events.dat <- readr::read_tsv(events_path,
                                skip=1,
                                col_names = col_names,
                                col_types = col_types,
                                na = c("", "NA", "-"))
  return(events.dat)
}


#' Parse PAVfinder 0.3.0 output
#'
#' @param events_path Path to the PAVfinder events file
#'
#' @return A data frame with parsed PAVfincer results
#' @export
parse_pavfinder_0.3.0 <- function (events_path) {
  # Set column names
  col_names = c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2',
                'name', 'score', 'strand1', 'strand2', 'orient1', 'orient2',
                'event', 'size', 'gene1', 'transcript1', 'transcript_break1',
                'exon1', 'exon_bound1', 'gene2', 'transcript2', 'transcript_break2',
                'exon2', 'exon_bound2', 'gene_5prime', 'gene_3prime', 'exon_5prime',
                'exon_3prime', 'feature', 'seq_id', 'seq_breaks', 'ins_seq',
                'homol_seq', 'homol_seq_coords', 'novel_seq', 'novel_seq_coords',
                'copy_number_change', 'repeat_seq', 'in_frame', 'probe',
                'support_span', 'spanning_reads', 'flanking_pairs', 'transcript1_tpm',
                'transcript2_tpm', 'filter')
  # Set column types
  col_types = readr::cols(
    chrom1 = readr::col_character(),
    start1 = readr::col_integer(),
    end1 = readr::col_integer(),
    chrom2 = readr::col_character(),
    start2 = readr::col_integer(),
    end2 = readr::col_integer(),
    name = readr::col_integer(),
    score = readr::col_character(),
    strand1 = readr::col_character(),
    strand2 = readr::col_character(),
    orient1 = readr::col_character(),
    orient2 = readr::col_character(),
    event = readr::col_character(),
    size = readr::col_integer(),
    gene1 = readr::col_character(),
    transcript1 = readr::col_character(),
    transcript_break1 = readr::col_character(),
    exon1 = readr::col_integer(),
    exon_bound1 = readr::col_character(),
    gene2 = readr::col_character(),
    transcript2 = readr::col_character(),
    transcript_break2 = readr::col_character(),
    exon2 = readr::col_integer(),
    exon_bound2 = readr::col_character(),
    gene_5prime = readr::col_character(),
    gene_3prime = readr::col_character(),
    exon_5prime = readr::col_integer(),
    exon_3prime = readr::col_integer(),
    feature = readr::col_character(),
    seq_id = readr::col_character(),
    seq_breaks = readr::col_character(),
    ins_seq = readr::col_character(),
    homol_seq = readr::col_character(),
    homol_seq_coords = readr::col_character(),
    novel_seq = readr::col_character(),
    novel_seq_coords = readr::col_character(),
    copy_number_change = readr::col_character(),
    repeat_seq = readr::col_character(),
    in_frame = readr::col_character(),
    probe = readr::col_character(),
    support_span = readr::col_character(),
    spanning_reads = readr::col_integer(),
    flanking_pairs = readr::col_integer(),
    transcript1_tpm = readr::col_double(),
    transcript2_tpm = readr::col_double(),
    filter = readr::col_character()
  )
  # Read the TSV data
  events.dat <- readr::read_tsv(events_path,
                                skip=2,
                                col_names = col_names,
                                col_types = col_types,
                                na = c("", "NA", "-", "na"))
  return(events.dat)
}

#' Parse output from RNA-SeQC v1.1.8
#'
#' @param rna_seqc_path file
#'
#' @return df
#' @export
#'
#' @examples
#'#' qc.dat <- rna_seqc_parser(rna_seqc_path)
rna_seqc_parser <- function ( rna_seqc_path ) {
  # Read the TSV data, using readr::read_tsv
  rna_seqc.dat <- readr::read_tsv(rna_seqc_path)
  return(rna_seqc.dat)
}

#' Parse gene-level output from Sailfish >=0.7.0
#'
#' @param sailfish_path file
#'
#' @return df
#' @export
#'
#' @examples
#'#' sailfish.genes.df <- sailfish_gene_parser_post_0.7.0(sailfish_path)
sailfish_gene_parser_post_0.7.0 <- function ( sailfish_path ) {
  # The sailfish files in the DB are just the directory name
  sailfish_file <- stringr::str_c(sailfish_path, '/quant.genes.sf.gz', '/')
  # Remove the trailing slash
  # http://stackoverflow.com/questions/23413331/
  #  how-to-remove-last-n-characters-from-every-element-in-the-r-vector
  sailfish_file <- gsub('.{1}$', '', sailfish_file)
  # Set column names directly
  sailfish_col_names <- c("Name",
                          "Length",
                          "EffectiveLength",
                          "TPM",
                          "NumReads")
  # Read the TSV data, specifying column types
  sailfish.dat <- readr::read_tsv(sailfish_file,
                                  skip=1,        # Skip header rows
                                  col_names = sailfish_col_names,
                                  col_types = "cdddd")
  return(sailfish.dat)
}

#' Parse isoform-level output from Sailfish <0.7.0
#'
#' @param sailfish_path file
#'
#' @return df
#' @export
#'
#' @examples
#'#' sailfish.df <- sailfish_isoform_parser_pre_0.7.0(sailfish_path)
sailfish_isoform_parser_pre_0.7.0 <- function ( sailfish_path ) {
  # The sailfish files in the DB are just the directory name
  sailfish_file <- stringr::str_c(sailfish_path, '/quant_bias_corrected.sf', '/')
  # Remove the trailing slash
  # http://stackoverflow.com/questions/23413331/
  #  how-to-remove-last-n-characters-from-every-element-in-the-r-vector
  sailfish_file <- gsub('.{1}$', '', sailfish_file)
  # Set column names directly
  sailfish_col_names <- c("Transcript",
                          "Length",
                          "TPM",
                          "RPKM",
                          "KPKM",
                          "EstimatedNumKmers",
                          "EstimatedNumReads")
  # Read the TSV data, specifying column types
  sailfish.dat <- readr::read_tsv(sailfish_file,
                                  skip=5,        # Skip header rows
                                  col_names = sailfish_col_names,
                                  col_types = "cdddddd")
  return(sailfish.dat)
}

#' Parse isoform-level output from Sailfish >=0.7.0
#'
#' @param sailfish_path file
#'
#' @return df
#' @export
#'
#' @examples
#'#' sailfish.df <- sailfish_isoform_parser_post_0.7.0(sailfish_path)
sailfish_isoform_parser_post_0.7.0 <- function ( sailfish_path ) {
  # The sailfish files in the DB are just the directory name
  sailfish_file <- stringr::str_c(sailfish_path, '/quant.sf.gz', '/')
  # Remove the trailing slash
  # http://stackoverflow.com/questions/23413331/
  #  how-to-remove-last-n-characters-from-every-element-in-the-r-vector
  sailfish_file <- gsub('.{1}$', '', sailfish_file)
  # Set column names directly
  sailfish_col_names <- c("Name",
                          "Length",
                          "EffectiveLength",
                          "TPM",
                          "NumReads")
  # Read the TSV data, specifying column types
  sailfish.dat <- readr::read_tsv(sailfish_file,
                                  skip=1,        # Skip header rows
                                  col_names = sailfish_col_names,
                                  col_types = "cdddd")
  return(sailfish.dat)
}


