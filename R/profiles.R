library(ggplot2)
library(dplyr)


# make data frame for ggplot
plotFacetProfiles <- function(data, 
                              fcol,
                              replicate.column.name, 
                              col,
                              reorderByRep = TRUE,
                              plotReps = TRUE,
                              ribbons = TRUE,
                              ...) {
  
  if (missing(fcol)) 
    stop("fcol missing, please specify fcol")
  if (!(fcol %in% fvarLabels(data)))
    stop("fcol not found in data")
  
  if (missing(col)) {
    n <- getMarkerClasses(data, fcol)
    col <- getStockcol()[seq(n)]
  }
  
  
  if (reorderByRep) {
    pd <- pData(data)
    if (missing(replicate.column.name))
      stop(message(paste("replicate.column.name must be provided")))
    new_order <- order(pd[[replicate.column.name]])
    data <- data[, new_order]
    message(paste("Reordering data according to replicate levels"))
  }
  fd <- fData(data)
  pd <- pData(data)
  names(col) <- getMarkerClasses(data, fcol)
  data <- exprs(data)
  intensities = NULL
  mrk = NULL
  if (missing(replicate.column.name)) {
    # message(paste("Replicate information not provided, assuming 1 replicate only"))
    repInfo <- rep(1, ncol(data))
    reps <- FALSE
  } else {
    repInfo <- pd[, replicate.column.name]
    reps <- TRUE
  }
  
  ## Re-organise data structure for ggplot 
  .rn <- rownames(data)
  .cn <- colnames(data)
  plot_data <- data.frame("id" = rep(.rn, ncol(data)),  
                          "fraction" = rep(.cn, each = nrow(data)), # variable
                          "intensities" = as.vector(data),  # value
                          "rep" = factor(rep(repInfo, each = nrow(data))),
                          "mrk" = rep(fd[, fcol], ncol(data)))
  ## convert the fraction names to factors using the within command
  plot_data <- within(plot_data, fraction <- factor(fraction, levels = colnames(data)))
  
  df <- plot_data %>% group_by(mrk, fraction, rep) %>%
    dplyr::summarise(min = min(intensities, na.rm = TRUE),
                     quant_05 = quantile(intensities, 0.25, na.rm = TRUE),
                     mean = mean(intensities, na.rm = TRUE),
                     quant_95 = quantile(intensities, 0.75, na.rm = TRUE),
                     max = max(intensities, na.rm = TRUE), .groups = "keep",
                     na.rm = TRUE)
  
  fracLev <- levels(df$fraction)
  if (plotReps) {
    repLev <- levels(df$rep)
    p <- ggplot()
    for(i in seq(repLev)){
      if (ribbons) {
        p <- p +
          geom_ribbon(data = subset(df, rep == repLev[i]),
                      mapping = aes(fraction, ymin=min, ymax=max, group = mrk, 
                                    colour = mrk, fill = mrk), 
                      linetype = 0, # remove upper and lower ribbon lines
                      alpha=0.5) +
          geom_line(data = subset(df, rep == repLev[i]),
                    mapping = aes(fraction, mean, group = mrk, color = mrk))
      } else {
        p <- p + 
          geom_line(data = subset(plot_data, rep == repLev[i]), 
                    aes(x = fraction, y = intensities, group = id, color = mrk),
                    alpha=0.3)
      }
    }
  } else {
    if (ribbons) {
      p <- 
        ggplot() + geom_ribbon(data = df,
                               mapping = aes(fraction, ymin=min, ymax=max, group = mrk, 
                                             colour = mrk, fill = mrk), 
                               linetype = 0, # remove upper and lower ribbon lines
                               alpha=0.5) +
        geom_line(data = df,
                  mapping = aes(fraction, mean, group = mrk, color = mrk))
    } else {
      p <- ggplot() + 
        geom_line(data = plot_data, 
                  aes(x = fraction, y = intensities, group = id, color = mrk),
                  alpha=0.3) 
      
    }
  }
  
  ## extract colours for organelles in the data 
  col <- c(col, "unknown" = "darkgrey")
  if (is.factor(df$mrk)) 
    col <- col[levels(df$mrk)]
  else
    col <- col[unique(df$mrk)]
  
  ## only show every other label on the x-axis for clarity
  lab_x <- fracLev
  lab_x[seq(2,length(lab_x),by=2)] <- ""
  
  ## plot data and customise
  p <- p + 
    scale_x_discrete(limits=fracLev, 
                     breaks = fracLev, 
                     labels = lab_x) +  # show only every other label on x-axis
    ylab("Normalised intensities") + xlab("") +
    scale_fill_manual(values = col, aesthetics = c("fill","colour")) +
    scale_color_manual(values = col, aesthetics = c("fill, colour")) +
    theme_bw() +
    theme(panel.spacing = unit(1, "lines"),
          legend.position = "none", 
          axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1),
          strip.text.x = element_text(size = 12, face="bold"),
          panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12 + 2, colour = rgb(0, 0, 0)),
          axis.text.y = element_text(size = 12, colour = rgb(0, 0, 0))) +
    facet_wrap(~ mrk, scales = "fixed", ... ) 
  return(p)
}


plotProtsOfInterest <- function(mylist,
                         foi = "P14142",
                         col = c("black", "blue", "cyan"),
                         replicate.column.name = "Replicate",
                         listNames = c("Basal", "Insulin", "CL")) {
  data <- lapply(mylist, exprs)
  data <- do.call(rbind, lapply(data, function(z) z[foi, ]))
  rownames(data) <- listNames
  listNames = c("Basal", "Insulin", "CL")
  p <- strsplit(colnames(data), "_Basal_")
  foo <- function(i) paste0(p[[i]][1], "_", p[[i]][2])
  colnames(data) <- sapply(seq(ncol(data)), foo)

  ## prep data for ggplot
  repInfo <- pData(mylist[[1]])[, "Replicate"]
  .rn <- rownames(data)
  .cn <- colnames(data)
  plot_data <- data.frame(id = rep(.rn, ncol(data)),
                          fraction = rep(.cn, each = nrow(data)), # variable
                          intensities = as.vector(data),  # value
                          rep = factor(rep(repInfo, each = nrow(data))))
  plot_data <- within(plot_data, fraction <- factor(fraction, levels = colnames(data)))

  ## line plot
  p <- ggplot()
  p <- p +
    geom_line(data = subset(plot_data, id=="Basal"),
              aes(x = fraction, y = intensities, group = id, color = col[1]),
              alpha=0.9) 
  p <- p +
    geom_line(data = subset(plot_data, id=="Insulin"),
              aes(x = fraction, y = intensities, group = id, color = col[2]),
              alpha=0.9)
  p <- p +
    geom_line(data = subset(plot_data, id=="CL"),
              aes(x = fraction, y = intensities, group = id, color = col[3]),
              alpha=0.9) 
  p <- p +
    # scale_x_discrete(limits=fracLev, breaks = fracLev[seq(1,length(fracLev),by=2)]) +  # show only every other label on x-axis
    ylab("Normalised intensities") +
    xlab("Fraction/replicate") +
    scale_fill_manual(values = col, aesthetics = c("fill","colour"),
                      labels = c("Basal", "Insulin", "CL")) +
    scale_color_manual(values = col, aesthetics = c("fill, colour"),
                       labels = c("Basal", "Insulin", "CL")) +
    theme_minimal() +
    theme(panel.spacing = unit(1, "lines"),
          # legend.position = "none",
          axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1),
          # strip.text.x = element_text(size = 8, face="bold"),
          # panel.background = element_rect(fill = "gray95"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10 + 2, colour = rgb(0, 0, 0)),
          axis.text.y = element_text(size = 10, colour = rgb(0, 0, 0)),
          legend.position = "top") +
    guides(color=guide_legend("")) +
    ggtitle(label = paste("Protein profiles for", foi))
  p 

}

