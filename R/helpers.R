# ----------------------
# Helpers
# ---------------------



countp = function(primer, file) {
  num_hits = Biostrings::vcountPattern(primer, ShortRead::sread(ShortRead::readFastq(file)), fixed = FALSE)
  return(sum(num_hits > 0))
}

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


#' Reverse complement
#'
#' Helper to reverse complement a DNA sequence that's a string
#'
#' @param seq string, DNA sequence
#'
#' @export
rev_comp = function(seq) {
  if (length(seq) < 1) {
    return(character(0))
  }
  else if (length(seq) == 1) {
    as(Biostrings::reverseComplement(Biostrings::DNAString(seq)), "character")
  }
  else {
    as(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq)), "character")
  }
}



#' Clip primers using cutadapt
#'
#' Uses the external cutadapt progrom to trim off primer sequences
#'
#' @param fwd_files list of file names of the forward reads
#' @param rev_files list of file names of the reverse reads
#' @param fwd_primer forward primer sequnce as a string
#' @param rev_primer reverse primer sequnce as a string
#' @param cutadapt_bin location of the cutadapt binary. If you're not sure, run
#'   \code{which cutadapt} on the command line to find out
#' @param samples character vector of samples names, if not provided the the forward file names will be used
#' @param max_n discard reads with more than this many Ns
#' @param q quality score for 3' end-trimming
#' @param min_len minimum length of sequence to retain
#'
#'
clip_primers = function(fwd_files, rev_files, fwd_primer, rev_primer,
                        out_dir, cutadapt_bin, samples = NULL,
                        max_n = 0, q = 20, min_len = 50) {

  cutadapt = path.expand(cutadapt_bin)

  # Make sure it works
  check = sys::exec_internal(cutadapt, args = "--version")

  attempt::stop_if(check$status == 1,
          msg = stringr::str_c(sys::as_text(check$stderr), collapse = "\n")
  )

  message(glue::glue("Found cutadapt {sys::as_text(check$stdout)}"))

  fwd_cut = file.path(path.expand(out_dir), basename(fwd_files))
  rev_cut = file.path(path.expand(out_dir), basename(rev_files))

  if (is.null(samples)) {
    samples = basename(fwd_files)
  }

  names(fwd_cut) = samples
  names(rev_cut) = samples


  cutadapt_args = list(r1_out = fwd_cut, r2_out = rev_cut, r1_in = fwd_files, r2_in = rev_files) %>%
    pmap(function(r1_out, r2_out, r1_in, r2_in) {
      c("-m", min_len, "-q", q,
        "--max-n", max_n,
        "-g", fwd_primer, "-G", rev_primer,
        "-a", rev_comp(rev_primer), "-A", rev_comp(fwd_primer),
        "-o", r1_out, "-p", r2_out,
        "--trimmed-only", "--report", "full",
        r1_in, r2_in)
    }
    )

  codes = cutadapt_args %>% purrr::set_names(samples) %>%
    furrr::future_imap(~sys::exec_wait(cutadapt, .x,
                                std_out = path.expand(file.path(out_dir, str_c(.y, ".log"))),
                                std_err = path.expand(file.path(out_dir, str_c(.y, ".log")))
    )
    )


  attempt::warn_if_any(codes != 0,
      msg = glue::glue("Sample(s) {names(codes)[codes != 0]} did not complete successfully.  Check log file"))

  message("cutadapt run completed")

}


#' Get number of uniques sequences from a dada2 object
#'
#' @export
getN = function(x) {
  sum(getUniques(x))
}
