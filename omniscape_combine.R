# Combining multiple window size outputs into a single output.
library(dplyr)
library(terra)

r67 <- "C:/Users/finkm/Documents/SyncroSim/Scenario-4/omniscape_Results/cum_currmap.tif"
r133 <- "C:/Users/finkm/Documents/SyncroSim/Scenario-5/Results_noSource/cum_currmap.tif"
r188 <- "C:/Users/finkm/Documents/SyncroSim/Scenario_6/Results7mi/cum_currmap.tif"
r268 <- "C:/Users/finkm/Documents/SyncroSim/Scenario_6/Results10mi/cum_currmap.tif"
outpth <- "D:/GIS/Projects/CPW_Statewide_Habt/Omniscape"
wopt <- list(gdal=c("COMPRESS=LZW", "TFW=YES", "BIGTIFF=YES"))

ras67 <- rast(r67)
ras133 <- rast(r133)
ras188 <- rast(r188)
ras268 <- rast(r268)

testarea <- ext(ras188)
window(ras67) <- testarea
allrows <- nrow(ras188)
allcols <- ncol(ras188)

v67 <- values(ras67, mat=FALSE)
v133 <- values(ras133, mat=FALSE)
v188 <- values(ras188, mat=FALSE)
v268 <- values(ras268, mat=FALSE)

ptile67 <- ntile(v67, 100)
ptile133 <- ntile(v133, 100)
ptile188 <- ntile(v188, 100)
ptile268 <- ntile(v268, 100)

plist <- list(rast(nrows=allrows, ncols=allcols, crs=crs(ras188), extent=testarea, vals=ptile67),
              rast(nrows=allrows, ncols=allcols, crs=crs(ras188), extent=testarea, vals=ptile133),
              rast(nrows=allrows, ncols=allcols, crs=crs(ras188), extent=testarea, vals=ptile188),
              rast(nrows=allrows, ncols=allcols, crs=crs(ras188), extent=testarea, vals=ptile268))
names(plist) <- c("radius67", "radius133", "radius188", "radius268")

ptstack <- rast(plist)
testout <- app(ptstack, fun=mean, filename=file.path(outpth, "omnitest2.tif"), wopt=wopt)
