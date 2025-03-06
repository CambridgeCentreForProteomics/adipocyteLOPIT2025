## add protein desciption info
addProteinInfo <- function(msnset, fcol) {
  
  .info <- fData(msnset)[[fcol]]
  
  ## Gene Name
  .protInfo <- sapply(.info, getGN, USE.NAMES = FALSE)
  if (any(duplicated(.protInfo))) {
    ## add an the uniprot id as an index for duplicated Gene Names
    ## we require unique gene names for plotting
    .ind <- which(duplicated(.protInfo))
    .protInfo[.ind] <- paste0(.protInfo[.ind], "_", featureNames(msnset)[.ind])
  }
  fData(msnset)[["Gene.Name"]] <- .protInfo
  
  ## Gene description
  .protInfo <- sapply(.info, getDesc, USE.NAMES = FALSE)
  if (any(duplicated(.protInfo))) {
    .ind <- which(duplicated(.protInfo))
    .protInfo[.ind] <- paste0(.protInfo[.ind], "_", featureNames(msnset)[.ind])
  }
  fData(msnset)[["Gene.Description"]] <- .protInfo

  
  ## Gene OS
  .protInfo <- sapply(.info, getOS, USE.NAMES = FALSE)
  fData(msnset)[["OS"]] <- .protInfo
  
  
  ## Gene OX
  .protInfo <- sapply(.info, getOX, USE.NAMES = FALSE)
  fData(msnset)[["OX"]] <- .protInfo
  
  
  ## Gene PE
  .protInfo <- sapply(.info, getPE, USE.NAMES = FALSE)
  fData(msnset)[["PE"]] <- .protInfo
  
  ## Gene SV
  .protInfo <- sapply(.info, getSV, USE.NAMES = FALSE)
  fData(msnset)[["SV"]] <- .protInfo
  
  return(msnset)
}


# ## Extract gene name from protein description
# addGeneName <- function(msnset, fcol, col.name = "Gene.Name") {
#   .descr <- fData(msnset)[[fcol]]
#   .gn <- sapply(.descr, getGN, USE.NAMES = FALSE)
#   if (any(duplicated(.gn))) {
#     ## add an the uniprot id as an index for duplicated Gene Names
#     ## we require unique gene names for plotting
#     .ind <- which(duplicated(.gn))
#     .gn[.ind] <- paste0(.gn[.ind], "_", featureNames(msnset)[.ind])
#   }
#   fData(msnset)[[col.name]] <- .gn
#   return(msnset)
# }
# 
# addDesc <- function(msnset, fcol, col.name = "Gene.Description") {
#   .descr <- fData(msnset)[[fcol]]
#   .gn <- sapply(.descr, getDesc, USE.NAMES = FALSE)
#   if (any(duplicated(.gn))) {
#     ## add an the uniprot id as an index for duplicated Gene Names
#     ## we require unique gene names for plotting
#     .ind <- which(duplicated(.gn))
#     .gn[.ind] <- paste0(.gn[.ind], "_", featureNames(msnset)[.ind])
#   }
#   fData(msnset)[[col.name]] <- .gn
#   return(msnset)
# }

## helpers to extract gene name for description in msnset
getGN <- function(x) {
  x <- as.character(x)
  pn2 <- strsplit(x, split = "GN=")
  ind <- grep(pattern = "PE=", pn2[[1]])
  if (length(ind) > 0)
    gn <- strsplit(pn2[[1]][ind], split = " PE=")[[1]][1]
  else
    gn <- pn2[[1]][length(pn2[[1]])]
  return(gn)
}

getDesc <- function(x) {
  x <- as.character(x)
  pn2 <- strsplit(x, split = " OS=")
  if (length(pn2) > 0)
    de <- pn2[[1]][1]
  else
    de <- pn2[[1]][length(pn2[[1]])]
  return(de)
}

getOS <- function(x) {
  x <- as.character(x)
  pn2 <- strsplit(x, split = " OX=")
  if (length(pn2) > 0) {
    pn2 <- strsplit(pn2[[1]][1], split = " OS=")
    os <- pn2[[1]][2]
  } else {
    os <- pn2[[1]][length(pn2[[1]])]
  }
  return(os)
}

getOX <- function(x) {
  x <- as.character(x)
  pn2 <- strsplit(x, split = " GN=")
  if (length(pn2) > 0) {
    pn2 <- strsplit(pn2[[1]][1], split = " OX=")
    ox <- pn2[[1]][2]
  } else {
    ox <- pn2[[1]][length(pn2[[1]])]
  }
  return(ox)
}

getPE <- function(x) {
  x <- as.character(x)
  pn2 <- strsplit(x, split = " SV=")
  if (length(pn2) > 0) {
    pn2 <- strsplit(pn2[[1]][1], split = " PE=")
    pe <- pn2[[1]][2]
  } else {
    pe <- pn2[[1]][length(pn2[[1]])]
  }
  return(pe)
}

getSV <- function(x) {
  x <- as.character(x)
  pn2 <- strsplit(x, split = " SV=")
  if (length(pn2) > 0) {
    sv <- pn2[[1]][2]
  } else {
    sv <- pn2[[1]][length(pn2[[1]])]
  }
  return(sv)
}
