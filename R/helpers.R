# ----------------------
# Helpers
# ---------------------



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



#' Get number of uniques sequences from a dada2 object
#'
#' @param x object to get unique sequences from.  See \code{\link[dada2]{getUniques}} for details.
#' @export
getN = function(x) {
  sum(getUniques(x))
}
