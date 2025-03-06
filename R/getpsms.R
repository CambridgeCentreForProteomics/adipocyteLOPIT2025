## Count number of PSMs with the same sequence
getPSMsCount <- function(object) {
  fData(object)$Annotated.Sequence <- toupper(fData(object)$Annotated.Sequence)
  seqs <- fData(object)$Annotated.Sequence
  counts <- sapply(seqs, function(z) length(which(z == seqs)))
  fData(object)$PSM.count <- counts
  return(object)
} 

## record the number of NA's across the fractions per PSM
recordPSMsNA <- function(object) {
  r <- apply(exprs(object), 1, function(z) sum(is.na(z)))
  fData(object)$PSM.count.na <- r
  return(object)
}

## record the maximum number of NA's in a PSM for a given peptide
## e.g. we have 7 PSMs, 1 of which had 2 missing values, 3 of which 
## had 1 missing value, we record the maximum number of NA's here as 3
maxPSMsNA <- function(object) {
  fData(object)$Annotated.Sequence <- toupper(fData(object)$Annotated.Sequence)
  seqs <- fData(object)$Annotated.Sequence
  ind <- lapply(seqs, function(z) 
    which(z == fData(object)$Annotated.Sequence))
  mymax <- sapply(ind, function(z) max(fData(object)[z, "PSM.count.na"]))
  fData(object)$PSM.count.na.max <- mymax
  return(object)
} 