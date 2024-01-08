##############################################################################
# Supervised Classification of landscape permeability
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Created 11/7/2023, last updated 01/08/2024
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

library(data.table)
library(dplyr)
library(randomForest)
library(terra)
library(foreach)
library(doSNOW)

source("model_func.R")

ncores <- 6
ntrees <- 4200
projdir <- "D:/GIS/Projects/CPW_Statewide_Habt"
spp <- "Elk"
run_name <- "Run3_Jan24"
classSampNums <- 5000
training <- rast(file.path(projdir, "Elk_perm_training3.tif"))
inputs <- fread("EnvInputs.csv")

# Create sampling points ----------------------------------------------------
trainpts <- spatSample(training, classSampNums, method="stratified", values=TRUE, 
                       as.points=TRUE)

# Make response table -------------------------------------------------------
outvals <- foreach(r=1:nrow(inputs),
                   .combine="cbind") %do%
  myextract(trainpts, inputs[r, 1:2])

outdt <- data.table(outvals)
respdt <- cbind(trainpts$Elk_perm_training3, outdt)
fwrite(respdt, file=file.path(projdir, paste0(spp, run_name, "_response.csv")))

## if already done:
respdt <- fread(file=file.path(projdir, paste0(spp, run_name, "_response.csv")))
## - ##

respdt <- setnames(respdt, c("V1"), c("Response"))
varstouse <- as.character(inputs[use==1, label])

chkdat <- numCorr(respdt[, ..varstouse])
chkdat

cdat <- c("Response", varstouse)
dat <- respdt[, ..cdat]
vnames <- varstouse

# Tune RF -------------------------------------------------------------------
sub1 <- dat %>% dplyr::filter(Response == 1) %>% sample_frac(0.1)
sub2 <- dat %>% dplyr::filter(Response == 2) %>% sample_frac(0.1)
sub3 <- dat %>% dplyr::filter(Response == 3) %>% sample_frac(0.1)
sub4 <- dat %>% dplyr::filter(Response == 4) %>% sample_frac(0.1)
sub5 <- dat %>% dplyr::filter(Response == 5) %>% sample_frac(0.1)

subtune <- as.data.table(bind_rows(sub1, sub2, sub3, sub4, sub5))

x <- tuneRF(dat[, c(2:7)], y=as.factor(dat[, Response]),
            ntreeTry = 151,
            stepFactor = 1.5,
            improve=0.01)
mt <- x[x[,2] == min(x[,2]),1]

# Run model -----------------------------------------------------------------

treeSubs <- ntrees/ncores

cl <- snow::makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)
rf.fit <- foreach(tree = rep(treeSubs,ncores),
                  .combine = randomForest::combine,
                  .packages = c("randomForest", "data.table"),
                  .multicombine = TRUE) %dopar% {
                     randomForest(as.factor(Response) ~.,
                                  data=dat,
                                  importance=TRUE,
                                  ntree=tree,
                                  mtry=mt,
                                  replace = TRUE,
                                  norm.votes = TRUE)
                  }
snow::stopCluster(cl)

ptab_rf.mod <- table(rf.fit$predicted, rf.fit$y)
err_rate <- 1-sum(diag(ptab_rf.mod))/sum(ptab_rf.mod)
rf.fit$confusion <- getConfusionMatrix(ptab_rf.mod)
rf.fit$err.rate <- err_rate
saveRDS(rf.fit, file=file.path(projdir, paste0(spp, run_name, ".rds")))
varImpPlot(rf.fit)
rf.fit$confusion

## Skip all that
rf.fit <- readRDS(file.path(projdir, paste0(spp, run_name, ".rds")))
vnames <- names(rf.fit$forest$ncat)
## - ##

# Test Area Spatial Prediction ---------------------------------------------------
testarea <- rast("D:/GIS/Projects/CPW_Statewide_Habt/LDIbase1_testarea.tif")
testext <- ext(testarea)
rasfiles <- inputs[label %in% vnames, raster]
spatlist <- lapply(rasfiles, rast)
names(spatlist) <- inputs[raster %in% rasfiles, label]
spatstack <- rast(spatlist)
teststack <- crop(spatstack, testext)
testout <- predict(teststack, rf.fit, type='response',
                   filename=file.path(projdir, paste0("test_", run_name, ".tif")),
                   wopt=list(gdal=c("COMPRESS=LZW", "TFW=YES", "BIGTIFF=YES")))


# Full Spatial Prediction --------------------------------------------------------
rasfiles <- inputs[label %in% vnames, raster]
spatlist <- lapply(rasfiles, rast)
names(spatlist) <- inputs[raster %in% rasfiles, label]
spatstack <- rast(spatlist)
allrows <- nrow(spatstack)
allcols <- ncol(spatstack)
#bs <- myblocksize(allrows)
wopt=list(gdal=c("COMPRESS=LZW", "TFW=YES", "BIGTIFF=YES"))

outx <- rast(nrows=allrows,
             ncols=ncol(spatstack),
             extent=ext(spatstack))

bs <- writeStart(outx, filename=file.path(projdir, paste0(spp, run_name, "_perm.tif")),
                 datatype='INT1U', wopt=wopt)

for(i in 1:bs$n) {
  iRow <- values(spatstack, row=bs$row[i], nrows=bs$nrows[i])
  subx <- procbloc_class(iRow, rf.fit)
  writeValues(outx, subx, bs$row[i], bs$nrows[i])
}

writeStop(outx)

# Nice try but no
# cl <- makeCluster(10, type = 'SOCK')
# registerDoSNOW(cl)
# 
# x <- foreach(i=1:max(bs$n), 
#              .combine = rbind) %do% {
#                startrow <- bs$row[i]
#                numrows <- bs$nrows[i]
#                colchunks <- seq(1, (numrows * allcols), by=allcols)
#                imageBlock <- values(spatstack, row=startrow, nrows=numrows)
#                subx <- foreach(r=1:numrows,
#                                .combine=rbind,
#                                .packages = c("terra", "randomForest")) %dopar% {
#                                  chnk <- imageBlock[colchunks[r]:(colchunks[r + 1] - 1), 1:5]
#                                  procbloc_class(chnk, rfmod)
#                                }
#                writeValues(outx, subx, startrow, numrows)
#              }
# 
# writeStop(outx)
# snow::stopCluster(cl)
