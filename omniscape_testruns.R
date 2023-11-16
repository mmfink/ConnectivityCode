# Note this just rescales and crops the resistance layer.
# TODO - would be nice to batch-run julia commands using different window sizes.

library(terra)

linearscl <- function(x, oldmin, oldmax, newmin, newmax) {
  oldrange <- oldmax - oldmin
  newrange <- newmax - newmin
  newx <- (((x - oldmin) * newrange) / oldrange) + newmin
  return(newx)
}

outpth <- "D:/GIS/Projects/CPW_Statewide_Habt/Omniscape"
wopt <- list(gdal=c("COMPRESS=LZW", "TFW=YES", "BIGTIFF=YES"))
testarea <- rast(file.path(outpth, "omnitest2.tif"))
ldi <- rast("D:/GIS/Projects/Scorecard/ldi_2023/LDI_8_NoAg.tif")
ldirange <- minmax(ldi)
rscldi <- linearscl(ldi, ldirange[1,1], ldirange[2,1], 1, 100)
testldi <- crop(rscldi, ext(testarea), 
                filename = file.path(outpth, "ldiNoAg_omnitest.tif"),
                wopt = wopt)
