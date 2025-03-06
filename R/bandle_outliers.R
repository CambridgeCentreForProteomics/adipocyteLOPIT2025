## Packages called here for convenience but should be moved when this code
## is transferred to bandle
library(coda)
library(gtools)

## ================================================
## Modified from plotConvergence in bandle package
## ================================================
# This function generates a rank plot for each chain for each condition. 
# For example, the barplot for chain 1 condition 1 tells us how many times 
# chain 1 appeared as having the first, second, third or fourth highest 
# number of outliers.


plotConvergenceRank <- function(params, seed = 1) {
  set.seed(seed)
  n <- as.numeric(seq(params@chains@chains[[1]]@n))
  nchains <- length(params@chains@chains)
  out <- vector("list", 2)
  for (j in 1:2) {
    out[[j]] <- vapply(params@chains@chains, 
                       function(x) colSums(x@outlier[[j]]), n)
    toplot <- apply(out[[j]], 1, rank, ties.method = "random")
    for (i in seq.int(nrow(toplot))) {
      chain_rank <- factor(toplot[i, ], levels = 1:nchains)
      barplot(table(chain_rank),
              main = paste("Condition ", j, "- Chain", i),
              xlab = "Rank",
              col =  brewer.pal(nchains, "Set1")[i])}
  }
}

## ================================================
## Taken and modified from mcmc_get_outliers
## compute the number of outliers for each iteration for each chain
## ================================================
bandle_get_outliers <- function(params) {
  ## check for validity
  stopifnot(inherits(params, "bandleParams"))
  ## get outliers
  out_cond1 <- lapply(params@chains@chains, function(mc) 
    coda::mcmc(colSums(1 - mc@outlier$cond1)))
  out_cond2 <- lapply(params@chains@chains, function(mc) 
    coda::mcmc(colSums(1 - mc@outlier$cond2)))
  out <- list(cond1 = out_cond1, 
              cond2 = out_cond2)
  return(out)
}

## ================================================
## Taken and modified from pRoloc TAGM MCMC code
## ================================================

plotOutliers <- function(params, 
                         auto.layout = TRUE) {
  

  ## check for validity
  stopifnot(inherits(params, "bandleParams")) 
  
  ## get outliers
  out <- bandle_get_outliers(params)
  
  ## set plotting parameters
  oldpar <- NULL
  on.exit(par(oldpar))
  if (auto.layout) {
    mfrow <- c(length(out$cond1), 4)
    oldpar <- par(mfrow = mfrow)
  }
  
  ## set names for plots
  chain_id <- seq_along(out$cond1)
  t_names1 <- paste0("Trace - condition ", 1, ", chain ", chain_id)
  d_names1 <- paste0("Density - condition ", 1, ", chain ", chain_id)
  t_names2 <- paste0("Trace - condition ", 2, ", chain ", chain_id)
  d_names2 <- paste0("Density - condition ", 2, ", chain ", chain_id)
  
  ## colours
  cols <- brewer.pal(n = 9, "Set1")[seq_along(chain_id)]
  
  ## use coda for plots
  for (i in seq_along(chain_id)) {
    traceplot(out$cond1[[i]], main = t_names1[i], col = cols[i])
    densplot(out$cond1[[i]], main = d_names1[i], col = cols[i])
    traceplot(out$cond2[[i]], main = t_names2[i], col = cols[i])
    densplot(out$cond2[[i]], main = d_names2[i], col = cols[i])
  }
}

## ================================================
## Calculate the Gelman diagnostics for all pairwise 
## combinations of chains in the bandle run
## ================================================

calculateGelman <- function(params) {
  ## get number of conditions (default is 2 but add here in case of pkg developemnt to allow more conditions)
  n_conds <- length(params@chains@chains[[1]]@outlier)
  ## initiate list
  gelmanStats <- vector("list", n_conds)
  ## Look over both conditions (n_conds = 2)
  for (j in seq(n_conds)) {
    ## get MCMC outputs
    outliers_cond <- lapply(params@chains@chains, function(mc)
      coda::mcmc(colSums(1 - mc@outlier[[j]])))
    ## get combinations
    p <- combinations(length(params@chains), 2)
    m <- vector("list", nrow(p))
    for (i in seq(nrow(p))) {
      m[[i]] <- as.vector(p[i, ])
    }
    ## calculate gelmans
    gelmanStats[[j]] <- sapply(m, function(z) 
      gelman.diag(outliers_cond[z])$psrf)
    rownames(gelmanStats[[j]]) <- c("Point_Est", "Upper_CI") 
    colnames(gelmanStats[[j]]) <- sapply(m, function(x) 
      paste0("comb_", paste0(x, collapse = "")))
  }
  return(gelmanStats)
}