% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cutadapt.R
\name{clip_primers}
\alias{clip_primers}
\title{Clip primers using cutadapt}
\usage{
clip_primers(fwd_files, rev_files, fwd_primer, rev_primer, out_dir,
  cutadapt_bin, samples = NULL, max_n = 0, q = 20, min_len = 50)
}
\arguments{
\item{fwd_files}{list of file names of the forward reads}

\item{rev_files}{list of file names of the reverse reads}

\item{fwd_primer}{forward primer sequnce as a string}

\item{rev_primer}{reverse primer sequnce as a string}

\item{out_dir}{path to directory to write the clipped files}

\item{cutadapt_bin}{location of the cutadapt binary. If you're not sure, run
\code{which cutadapt} on the command line to find out}

\item{samples}{character vector of samples names, if not provided the the forward file names will be used}

\item{max_n}{discard reads with more than this many Ns}

\item{q}{quality score for 3' end-trimming}

\item{min_len}{minimum length of sequence to retain}
}
\description{
Uses the external cutadapt progrom to trim off primer sequences
}
