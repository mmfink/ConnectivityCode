##############################################################################
# Data prep and explore functions, to be run before creating models
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Created ages ago, last updated 01/10/2024
#
# -----------------------------------------------------------------
# Code licensed under the GNU General Public License version 3.
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see https://www.gnu.org/licenses/
#############################################################################
require(dplyr)
require(data.table)
library(ggplot2)
library(terra)

myextract <- function(inpts, inrow){ 
  # Extract raster values at a set of points
  # used within a foreach statement
  rasname <- inrow$raster
  inlbl <- inrow$label
  inras <- rast(rasname)
  x <- extract(inras, inpts, ID=FALSE)
  names(x) <- inlbl
  return(x)
}

comb <- function(x, ...) {
  #https://stackoverflow.com/questions/19791609/saving-multiple-outputs-of-foreach-dopar-loop
  # used within a foreach loop
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

numCorr <- function(dat){
  # create a matrix of variable covariance
  # dat; data.table, just the numeric input vars, nothing else
  cmat <- cor(dat, use = "pairwise.complete.obs")
  smat <- cor(dat, method = "spearman", use = "pairwise.complete.obs")
  
  if (dim(dat)[1] < 2000) {
    kmat <- cor(dat, method = "kendall", use = "pairwise.complete.obs")
  } else {
    s <- sample(seq(1:dim(dat)[1]), size = 2000, replace = FALSE)
    kmat <- cor(dat[s, ], method = "kendall", use = "pairwise.complete.obs")
  }
  
  cmat = pmax(abs(cmat), abs(smat), abs(kmat), na.rm = TRUE)
  cmat[is.na(cmat)] <- 0
  
  return(cmat)
}

getConfusionMatrix <- function(tbl) {
  #https://stats.stackexchange.com/questions/35609/why-do-i-need-bag-composition-to-calculate-oob-error-of-combined-random-forest-m/35613#35613
  class.error = vector()
  
  for (i in 1:nrow(tbl)) {
    rowSum = sum(tbl[i,])
    accurate = diag(tbl)[i]
    error = rowSum - accurate
    
    class.error[i] = error / rowSum
  }
  return(cbind(tbl, class.error))
}

myblocksize <- function(totrow){
  minrows <- floor(totrow/60)
  xtra <- totrow - (minrows * 60)
  if(xtra > 0){
    totb <- 61
    numrow <- c(rep(minrows, 60), xtra)
  } else {
    totb <- 60
    numrow <- rep(minrows, 60)
  }
  blocktbl <- data.table(row=seq(from=1, to=totrow, by=minrows),
                         nrows=numrow,
                         n=1:totb)
  return(blocktbl)
}

# FIXME # Prediction functions to feed to workers
procbloc_class <- function(iBlock, randfor){
  predValues <- predict(randfor, iBlock, type='response')
  classValues <- as.numeric(levels(predValues))[predValues]
  chunk <- matrix(classValues, nrow = 1, byrow = T)
  return(chunk)
}

classpal <- c("#C2523C", "#F2B60E", "#77ED00", "#1BAA7D", "#0B2C7A")

plotRF <- function(rfOutput){
  # https://stackoverflow.com/questions/39330728/plot-legend-random-forest-r
  # Get OOB data from plot and coerce to data.table
  oobData = as.data.table(plot(rfOutput))
  
  # Define trees as 1:ntree
  oobData[, trees := .I]
  
  # Cast to long format
  oobData2 = melt(oobData, id.vars = "trees")
  setnames(oobData2, "value", "error")
  
  # Plot using ggplot
  out <- ggplot(data = oobData2, aes(x = trees, y = error, color = variable)) + 
    geom_line(linewidth=0.8) + theme_light() + labs(title="Out of Bag Error") +
    scale_color_manual(values=c("#000000", classpal)) +
    guides(color=guide_legend(title="Classes"))
  
  return(out)
}