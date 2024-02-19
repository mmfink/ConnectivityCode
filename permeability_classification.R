##############################################################################
# Supervised Classification of landscape permeability
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Created 11/7/2023, last updated 02/07/2024
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

source("model_func.R")

ncores <- 4
ntrees <- 2001
projdir <- "H:/CPW_Statewide_Habt"
spp <- "Elk"
run_name <- "Run2_Feb24"
classSampNums <- 6000
mindist <- 120
training <- rast(file.path(projdir, "Elk_perm_training9.tif"))
inputs <- fread("EnvInputs.csv")

# Create sampling points ----------------------------------------------------
trainpts <- spatSample(training, classSampNums, method="stratified", values=TRUE, 
                       as.points=TRUE)
# remove points too close together, if any
trainpts <- filterPts(trainpts, mindist)
trainpts$idx <- seq.int(1, nrow(trainpts))

# Make response table -------------------------------------------------------
outvals <- foreach(r=1:nrow(inputs),
                   .combine="cbind") %do%
  myextract(trainpts, inputs[r, 1:2])

outdt <- data.table(outvals)
respdt <- cbind(as.data.frame(trainpts, geom="XY"), outdt)
setDT(respdt, key="idx")
respdt <- setnames(respdt, c("Band_1"), c("Response"))
fwrite(respdt, file=file.path(projdir, paste0(spp, run_name, "_response.csv")))

## if already done:
respdt <- fread(file=file.path(projdir, paste0(spp, run_name, "_response.csv")))
## - ##

varstouse <- as.character(inputs[use==1, label])

chkdat <- numCorr(respdt[, ..varstouse])
chkdat

cdat <- c("Response", varstouse)
dat <- respdt[, ..cdat]

# Tune RF -------------------------------------------------------------------
sub1 <- dat %>% dplyr::filter(Response == 1) %>% sample_frac(0.2)
sub2 <- dat %>% dplyr::filter(Response == 2) %>% sample_frac(0.2)
sub3 <- dat %>% dplyr::filter(Response == 3) %>% sample_frac(0.2)
sub4 <- dat %>% dplyr::filter(Response == 4) %>% sample_frac(0.2)
sub5 <- dat %>% dplyr::filter(Response == 5) %>% sample_frac(0.2)

subtune <- as.data.table(bind_rows(sub1, sub2, sub3, sub4, sub5))

x <- tuneRF(subtune[, ..varstouse], y=as.factor(subtune[, Response]),
            ntreeTry = 501,
            stepFactor = 1.5,
            improve=0.01)
mt <- x[x[,2] == min(x[,2]),1]

# Run model -----------------------------------------------------------------
rf.fit <- randomForest(as.factor(Response) ~.,
                       data=dat,
                       importance=TRUE,
                       ntree=ntrees,
                       mtry=mt,
                       replace = TRUE,
                       sampsize = c(2000,2000,2000,2000,2000),
                       norm.votes = TRUE)

plotRF(rf.fit)
varImpPlot(rf.fit)
rf.fit$confusion
saveRDS(rf.fit, file=file.path(projdir, paste0(spp, run_name, ".rds")))

## if already done:
rf.fit <- readRDS(file.path(projdir, paste0(spp, run_name, ".rds")))
varstouse <- names(rf.fit$forest$ncat)

# Test Spatial Prediction --------------------------------------------------------
testarea <- rast("H:/CPW_Statewide_Habt/LDIbase1_testarea.tif")
testext <- ext(testarea)
rasfiles <- inputs[label %in% varstouse, raster]
spatlist <- lapply(rasfiles, rast)
names(spatlist) <- inputs[raster %in% rasfiles, label]
spatstack <- rast(spatlist)
teststack <- crop(spatstack, testext)
testout <- predict(teststack, rf.fit, type='response', cores=ncores,
                   cpkgs=c("randomForest"),
                   filename=file.path(projdir, paste0("test_", run_name, ".tif")),
                   wopt=list(gdal=c("COMPRESS=LZW", "TFW=YES", "BIGTIFF=YES")))


# Full Spatial Prediction --------------------------------------------------------
rasfiles <- inputs[label %in% varstouse, raster]
spatlist <- lapply(rasfiles, rast)
names(spatlist) <- inputs[raster %in% rasfiles, label]
spatstack <- rast(spatlist)
wopt=list(gdal=c("COMPRESS=LZW", "TFW=YES", "BIGTIFF=YES"))
permout <- predict(spatstack, rf.fit, type="response", cores=ncores,
                   cpkgs=c("randomForest"),
                   filename=file.path(projdir, paste0(spp, run_name, ".tif")),
                   wopt=wopt)
## Took 3.5 hours


