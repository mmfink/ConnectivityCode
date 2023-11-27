# Updated for 11/29 presentation
library(dplyr)
library(terra)

outpth <- "D:/GIS/Projects/CPW_Statewide_Habt/Omniscape"
wopt <- list(gdal=c("COMPRESS=LZW", "TFW=YES", "BIGTIFF=YES"))

# Mean ntiles of Raw Current ------------------------------------------------
r133 <- "C:/Users/finkm/Documents/SyncroSim/redo_tests/Results5mi/cum_currmap.tif"
r268 <- "C:/Users/finkm/Documents/SyncroSim/redo_tests/Results10mi/cum_currmap.tif"
r536 <- "C:/Users/finkm/Documents/SyncroSim/redo_tests/Results20mi/cum_currmap.tif"
r1340 <- "C:/Users/finkm/Documents/SyncroSim/redo_tests/Results50mi/cum_currmap.tif"

ras133 <- rast(r133)
ras268 <- rast(r268)
ras536 <- rast(r536)
ras1340 <- rast(r1340)

testarea <- ext(ras133)
allrows <- nrow(ras133)
allcols <- ncol(ras133)

v133 <- values(ras133, mat=FALSE)
v268 <- values(ras268, mat=FALSE)
v536 <- values(ras536, mat=FALSE)
v1340 <- values(ras1340, mat=FALSE)

ptile133 <- ntile(v133, 100)
ptile268 <- ntile(v268, 100)
ptile536 <- ntile(v536, 100)
ptile1340 <- ntile(v1340, 100)

plist <- list(rast(nrows=allrows, ncols=allcols, crs=crs(ras133), extent=testarea, vals=ptile133),
              rast(nrows=allrows, ncols=allcols, crs=crs(ras133), extent=testarea, vals=ptile268),
              rast(nrows=allrows, ncols=allcols, crs=crs(ras133), extent=testarea, vals=ptile536),
              rast(nrows=allrows, ncols=allcols, crs=crs(ras133), extent=testarea, vals=ptile1340))
names(plist) <- c("radius133", "radius268", "radius536", "radius1340")

ptstack <- rast(plist)
testout <- app(ptstack, fun=mean, filename=file.path(outpth, "omnitest3.tif"), wopt=wopt)

# Mean of Rescaled Normalized Current ---------------------------------------
n133 <- "C:/Users/finkm/Documents/SyncroSim/redo_tests/Results5mi/normalized_cum_currmap.tif"
n268 <- "C:/Users/finkm/Documents/SyncroSim/redo_tests/Results10mi/normalized_cum_currmap.tif"
n536 <- "C:/Users/finkm/Documents/SyncroSim/redo_tests/Results20mi/normalized_cum_currmap.tif"
n1340 <- "C:/Users/finkm/Documents/SyncroSim/redo_tests/Results50mi/normalized_cum_currmap.tif"

nras133 <- rast(n133)
nras268 <- rast(n268)
nras536 <- rast(n536)
nras1340 <- rast(n1340)

raslist <- list(nras133, nras268, nras536, nras1340)
outlist <- list()

normrescale <- function(inras){
  setMinMax(inras, force=TRUE)
  inrange <- minmax(inras)
  orange1 <- 0.6 - inrange[1,1]
  orange2 <- inrange[2,1] - 1.4
  outras <- ifel(inras < 0.6, ((inras - inrange[1,1])*9)/orange1,
                 ifel(inras > 1.4, (((inras - 1.4)*10)/orange2) + 10, 10))
  return(outras)
}

for(ras in raslist){
  rscale <- normrescale(ras)
  outlist <- append(outlist, rscale)
}

names(outlist) <- c("radius133", "radius268", "radius536", "radius1340")

normout <- app(outlist, fun=mean, 
               filename=file.path(outpth, "normalizedtest3.tif"), 
               wopt=wopt)
