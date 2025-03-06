# Utility functions useful across all analyses
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(grid))
suppressMessages(library(VennDiagram))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(MSnbase))
suppressMessages(library(robustbase))
suppressMessages(library(biobroom))
suppressMessages(library(archive))



print_n_feature <- function(features_df, message){
  cat(sprintf("%s\t%s\n", length(rownames(features_df)), message))
}

print_n_prot <- function(features_df, master_protein_col){
  cat(sprintf("These features are associated with %s master proteins\n",
              length(unique(features_df[[master_protein_col]]))))}

print_summaries <- function(features_df, master_protein_col, message){
  print_n_feature(features_df, message)
  print_n_prot(features_df, master_protein_col)
}

removeCrap <- function(obj, protein_col="Protein.Accessions"){
  cat(sprintf("Input data: %s rows\n", nrow(obj)))
  obj <- obj %>% filter(!grepl("cRAP", !!as.symbol(protein_col), ignore.case=FALSE))
  cat(sprintf("Output data: %s rows\n", nrow(obj)))
  return(obj)
}

# Steps are:
# 1. Filter modifications
# 2. Center-median normalise PSM
# 3. Aggregate PSM level data to unique peptide sequence + modification
# 4. Sum normalise
# 5. Impute missing values
# 6. (Optional) Aggregate to unique peptide sequence, sum normalise and impute missing values (starts from step 4)
# 7. (Optional) Aggregate to unique protein ID, sum normalise and impute missing values (starts from step 6 prior to sum normalisation)
#
# Input :
#   raw_psm: [REQUIRED] The PSM level dataframe. This can be obtained by parsing PD output with `parse_features`. If PSMs contain PTMs,
#             you may want to further parse with `parsePTMScores`, `addPTMPositions` and `addSiteSequence` first.
#   sample_infile: [REQUIRED] Filename for table mapping TMT tags to sample names. e.g:
#                   Tag     Sample_name
#                   126     100
#                   127N    400
#                   [...]   [...]

# ----------------------------------------------------------------------------------------------------------------------
# Function	: parse_features
# Aim		: Parse output from PD and filter
#
# Steps are:
# 1. Read in PD data
# 2. Exclude features without a master protein
# 3. (Optional) Exclude features without a unique master protein (Number.of.Protein.Groups==1)
# 4. (Optional) Exclude features matching a cRAP protein
# 5. (Optional) Exclude features matching a proteins associated with a cRAP protein (see below)
# 6. Filter out features without quantification values (only if (TMT or SILAC) & peptide level input)
#
# "Associated cRAP"
# I've observed that the cRAP database does not contain all the possible cRAP proteins.
# E.g some features can be assigned to a keratin which is not in the cRAP database.
# In most cases, these "associated cRAP" proteins will have at least one feature shared with a cRAP
# protein. Thus, use this to remove them: see simplified example below
#
# feature Protein.Accessions          Master.Protein.Accessions
# f1      protein1, protein2, cRAP1   protein 1
# f2      protein1, protein3          protein 3
# f3      protein2                    protein 2
#
# Here, f1 indicates that protein1 and protein2 are associated with a cRAP protein.
# f2 and f3 are therefore filtered out as well, regardless of the Master.Protein.Accession column value
#
#
# Input:
#    infile: [REQUIRED] file containing output from Proteome Discoverer
#    master_protein_col: [REQUIRED; DEFUALT="Master.Protein.Accessions"]
#                        The column containing the master protein.
#    protein_col: [REQUIRED; DEFUALT="Protein.Accessions"]
#                 The column containing all the protein matches.
#    unique_master: [DEFAULT=TRUE] Filter out features without a unique master protein
#    silac: [DEFAULT=FALSE] SILAC experiment
#    TMT: [DEFAULT=FALSE] TMT experiment
#    level: [REQUIRED; DEFAULT="peptide"] Input level. Must be "peptide" or "PSM
#    filter_crap: [DEFAULT=TRUE] Filter out the features which match a cRAP protein
#    crap_fasta: [DEFAULT=NULL] Fasta file containing the cRAP proteins. Expects fasta header format thusly
#    >sp|cRAP002|P02768|ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2, e.g Uniprot in 3rd position ('|' delimited)
#    filter_associated_crap: [DEFAULT=TRUE] Filter out the features which match a cRAP associated protein
#
# Output: Filtered PD results
# ----------------------------------------------------------------------------------------------------------------------
parse_features <- function(infile,
                           sep="\t",
                           master_protein_col="Master.Protein.Accessions",
                           protein_col="Protein.Accessions",
                           unique_master=TRUE,
                           silac=FALSE,
                           TMT=FALSE,
                           level="peptide",
                           filter_crap=TRUE,
                           crap_fasta=NULL,
                           filter_associated_crap=TRUE,
                           protein_group_col = "Number.of.Protein.Groups") {

  
  if(!level %in% c("PSM", "peptide")){
    stop("level must be PSM or peptide")
  }
  
  features_df <- read.delim(archive_read(infile),  sep=sep, header=T, stringsAsFactors=FALSE)

  # add this line as protein_group_col is not always called "Number of protein groups"
  colnames(features_df)[which(colnames(features_df) == protein_group_col)] <- "Number.of.Protein.Groups"
  
  cat("Tally of features at each stage:\n")
  
  print_summaries(features_df, master_protein_col, "All features")
  
  features_df <- features_df %>% filter(UQ(as.name(master_protein_col))!="")
  print_summaries(features_df, master_protein_col, "Excluding features without a master protein")
  
  if(filter_crap){
    
    if(is.null(crap_fasta)){
      stop('must supply the crap fasta argument to filter cRAP proteins')
    }
    
    con <- file(crap_fasta, open="r")
    lines <- readLines(con)
    
    #print(lines %>% strsplit(split="\\|"))
    crap_proteins <- lines %>% strsplit("\\|") %>% lapply(function(x){
      if(substr(x[[1]],1,1)!=">"){
        return()
      }
      else{
        return(x[[3]])
      }
    }) %>% unlist()
    
    close(con)
    
    if(filter_associated_crap){
      associated_crap <- features_df %>%
        filter((UQ(as.name(master_protein_col)) %in% crap_proteins)|
                 grepl("cRAP", !!as.symbol(protein_col), ignore.case=FALSE)) %>%
        pull(UQ(as.name(protein_col))) %>%
        strsplit("; ") %>%
        unlist()
      associated_crap <- associated_crap[!grepl("cRAP", associated_crap)]
    
    }

    features_df <- features_df %>% filter(!UQ(as.name(master_protein_col)) %in% crap_proteins,
                                          !grepl("cRAP", !!as.symbol(protein_col), ignore.case=FALSE))
    print_summaries(features_df, master_protein_col, "Excluding features matching a cRAP protein")
    
    if(filter_associated_crap){
      cat(sprintf("Identified an additional %s proteins as 'cRAP associated'\n", length(associated_crap)))

      if(length(associated_crap)>0){
        # remove isoforms
        associated_crap_no_isoform <- unique(sapply(strsplit(associated_crap, "-"), "[[", 1))
        associated_crap_regex <- paste(associated_crap_no_isoform, collapse="|")
        features_df <- features_df %>% filter(!grepl(associated_crap_regex, UQ(as.name(protein_col))))
        print_summaries(features_df, master_protein_col, "Excluding features associated with a cRAP protein")
      }
    }
    
  }

  if(unique_master){
    features_df <- features_df %>% filter(Number.of.Protein.Groups==1)
    print_summaries(
      features_df, master_protein_col, "Excluding features without a unique master protein")
  }
  
  
  if(silac|TMT & level=="peptide"){
    features_df <- features_df %>% filter(Quan.Info=="")
    print_summaries(features_df, master_protein_col, "Excluding features without quantification")
  }
  return(features_df)
}

