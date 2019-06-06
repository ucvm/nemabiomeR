# Count primers



#' Count primers
#'
#' This function counts the number of times it finds a match to the given primer in a fastq file
#'
#' @param r1_files list of fastq files with the forward read
#' @param fwd_primer forward primer sequence as a string
#' @param rev_primer reverse primer sequence as a string
#' @param r2_files list of fastq files with the reverse read, optional
#'
#' @return a tibble with the counts of the primers hits in different orientations, one row per file
#'
#' @export
count_primers = function(r1_files, fwd_primer, rev_primer, r2_files = NULL) {

  out = r1_files %>% purrr::set_names(basename(.)) %>%
    purrr::map_dfr(function(f) {
      list(fwd = fwd_primer, fwd_rev = rev_comp(fwd_primer), rev = rev_primer, rev_rev = rev_comp(rev_primer)) %>%
        purrr::map_int(~countp(.x, f)) %>%
        tibble::enframe() %>%
        tidyr::spread(name, value)
    }, .id = "filename")

  if (!is.null(r2_files)) {
    out2 = r2_files %>% purrr::set_names(basename(.)) %>%
      purrr::map_dfr(function(f) {
        list(fwd = fwd_primer, fwd_rev = rev_comp(fwd_primer), rev = rev_primer, rev_rev = rev_comp(rev_primer)) %>%
          purrr::map_int(~countp(.x, f)) %>%
          tibble::enframe() %>%
          tidyr::spread(name, value)
      }, .id = "filename")

    out = dplyr::bind_rows(out, out2)
  }

  return(out)

}



countp = function(primer, file) {
  num_hits = Biostrings::vcountPattern(primer, ShortRead::sread(ShortRead::readFastq(file)), fixed = FALSE)
  return(sum(num_hits > 0))
}
