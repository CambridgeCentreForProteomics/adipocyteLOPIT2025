library(colorspace)
## ===========pretty t-SNE===========
prettyMap <- function(map_matrix, object, 
                      fcol = "markers",
                       main = "", 
                       mainCol = getStockcol(), 
                       outlineCol = darken(getStockcol()), 
                       pchUn = 21, lwd = 1,
                       ...) {
  fData(object)[, "tmp"] <- fData(object)[, fcol]
  setUnknowncol(NULL)
  setUnknowncol(paste0(getUnknowncol(), 90))
  # setStockcol(NULL)
  plot2D(map_matrix, method = "none",  methargs = list(object),
         fcol = NULL, cex = 1, pch = pchUn, grid = FALSE, main = main, 
         cex.axis = 1.5,
         cex.lab = 1.5, lwd = lwd,
         bg = getUnknowncol(), col = darken(getUnknowncol()), ...)
  cl <- getMarkerClasses(object, fcol = "tmp")
  names(cl) <- getMarkerClasses(object, fcol = "tmp")
  col1 <- outlineCol
  col2 <- mainCol
  names(col1) <- cl
  names(col2) <- cl
  for (i in seq(cl)) {
    ind <- which(fData(object)[, "tmp"] == cl[i])
    points(map_matrix[ind, ], col = col1[i], pch = 21, 
           bg = col2[i], cex = 1.5, lwd = lwd)
  }
}

prettyMap_overlay <- function(map_matrix, object, fcol = "markers",
                               main = "", 
                               mainCol = paste0(getStockcol(), 70), 
                               outlineCol = darken(getStockcol()), ...) {
  fData(object)[, "tmp"] <- fData(object)[, fcol]
  setUnknowncol(NULL)
  setUnknowncol(paste0(getUnknowncol(), 70))
  plot2D(map_matrix, method = "none",  methargs = list(object),
         fcol = NULL, pch = 19, grid = FALSE, main = main, 
         cex.axis = 1.5,
         cex.lab = 1.5, 
         col = getUnknowncol(), ...)
  cl <- getMarkerClasses(object, fcol = "tmp")
  names(cl) <- getMarkerClasses(object, fcol = "tmp")
  col1 <- outlineCol
  col2 <- mainCol
  names(col1) <- cl
  names(col2) <- cl  
  for (i in seq(cl)) {
    ind <- which(fData(object)[, "tmp"] == cl[i])
    points(map_matrix[ind, ], col = col1[i], pch = 21, 
           bg = col2[i], cex = 1.5)
  }
}