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
