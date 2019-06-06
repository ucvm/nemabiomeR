# ----------------------------------------------
# Script to run the nemabiome analysis pipeline
# ----------------------------------------------

# Version 0.1.0

library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(stringr)
library(purrr)
library(nemabiomeR)


# Parameters to set -------------------------------------------------------

path = "~/MiSeq_Run"
cutadapt_bin = "~/miniconda3/bin/cutadapt"
fwd_primer = "ACGTCTGGTTCAGGGTTGTT"
rev_primer = "TTAGTTTCTTTTCCTCCGCT"

# if you've got multiple cores
plan(multiprocess)

# Setup -------------------------------------------------------------------


fwd_files = sort(list.files(path, pattern = "R1", full.names = TRUE))
rev_files = sort(list.files(path, pattern = "R2", full.names = TRUE))

stop_if_any(list(is_empty(fwd_files), is_empty(rev_files)),
            msg = "Couldn't find files")


samples = str_extract(basename(fwd_files), "^[^_]+")

stop_if_not(length(samples) > 0,
            msg = "Couldn't detect sample names")


names(fwd_files) = samples
names(rev_files) = samples


count_primers(fwd_files[1:4], fwd_primer, rev_primer, rev_files[1:4])



# Primer removal ----------------------------------------------------------

# Create an output directory to store the clipped files
cut_dir =file.path(path, "cutadapt")
if (!dir.exists(cut_dir)) dir.create(cut_dir)

clip_primers(fwd_files, rev_files, fwd_primer, rev_primer, cut_dir, cutadapt_bin = cutadapt_bin)

# quick check that we got something
head(list.files(cut_dir))



# Check quality scores ----------------------------------------------------

plotQualityProfile(fwd_cut[1:2]) + ggtitle("Forward")
plotQualityProfile(rev_cut[1:2]) + ggtitle("Reverse")

# Same as for the clippling we create an output directory to store the filtered files
filt_dir = file.path(path, "filtered")
if (!dir.exists(filt_dir)) dir.create(filt_dir)

fwd_filt = file.path(filt_dir, basename(fwd_files))
rev_filt = file.path(filt_dir, basename(rev_files))

names(fwd_filt) = samples
names(rev_filt) = samples

filtered_out = filterAndTrim(
  fwd = fwd_cut,
  filt = fwd_filt,
  rev = rev_cut,
  filt.rev = rev_filt,
  maxEE = c(2, 5),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
  )

head(filtered_out)



# run dada2 core pipeline -------------------------------------------------

# -- Errors
err_fwd = learnErrors(fwd_filt, multithread = TRUE)
err_rev = learnErrors(rev_filt, multithread = TRUE)

plotErrors(err_fwd, nominalQ = TRUE)

# --- Denoise
dada_fwd = dada(fwd_filt, err = err_fwd, multithread = TRUE)
dada_rev = dada(rev_filt, err = err_rev, multithread = TRUE)


# --- Merge pairs
mergers = mergePairs(
  dadaF = dada_fwd,
  dadaR = dada_rev,
  derepF = fwd_filt,
  derepR = rev_filt,
  maxMismatch = 1,
  verbose = TRUE
)


# --- Sequence table

seqtab = makeSequenceTable(mergers)
dim(seqtab)


# --- Remove chimeras
seqtab_nochim = removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab_nochim)



# How did we do? ----------------------------------------------------------

# --- Sequence length distribution
table(nchar(getSequences(seqtab_nochim)))


# --- Number of reads at each step
track = cbind(
  filtered_out,
  sapply(dada_fwd, getN),
  sapply(dada_rev, getN),
  sapply(mergers, getN),
  rowSums(seqtab_nochim)
)

colnames(track) = c("raw", "filtered", "denoised_fwd", "denoised_rev", "merged", "no_chim")
rownames(track) = samples
head(track)

