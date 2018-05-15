require(bfastSpatial)
require(bfast)
require(bayts)
require(rgdal)
require(raster)
require(sp)  


#load s1
s1vv <- brick("/media/DATA4/balli001/Sumatra_test_site/S1VV_stack_mt55_test.gri")
s1vv_date <- substr(getZ(s1vv),1,10)
s1vv_date


s1vh <- brick("/media/DATA4/balli001/Sumatra_test_site/S1VH_stack_mt55_test.gri")
s1vh_date <- substr(getZ(s1vh),1,10)
s1vh_date

#load landsat ndvi
ndvi <- brick("/media/DATA4/balli001/Sumatra_test_site/NDVI_res_stack_test.gri")
ndvi_date <- substr(getZ(ndvi),1,10)
ndvi_date

plot(s1vv,1)
plot(ndvi,30)

#hansen 2015 forest mask, 1 = forest, NA=non-forest
fm <- raster("/media/DATA4/balli001/Sumatra_test_site/Hansen_GFC-2015_10N_100E_2015_forest_mask_test.tif")
fm <- projectRaster(fm,s1vv)
plot(fm)

#Subsetting
sub <- TRUE
if (sub==TRUE){
  e_new <- c(101.85,101.9,0.90,1)
  s1vv <- crop(s1vv,e_new)
  s1vh <- crop(s1vh,e_new)
  ndvi <- crop(ndvi,e_new)
  fm <- crop(fm,e_new)
  } 


plot(s1vh,1)
plot(s1vh,88)
cell <- click(s1vv, n=1, cell=TRUE)[,1]
cell <- 10424

tsSvv <- bfastts(as.vector(s1vv[cell]),as.Date(s1vv_date),type=c("irregular"))
tsSvh <- bfastts(as.vector(s1vh[cell]),as.Date(s1vh_date),type=c("irregular"))
tsndvi <- bfastts(as.vector(ndvi[cell])/1000,as.Date(ndvi_date),type=c("irregular"))

plotts(list(tsSvv,tsSvh,tsndvi),ylimL=list(c(-25,0),c(-25,0),c(-1,1)),labL=list("S1 VV [dB]","S1 VH [dB]","Landsat NDVI"))



############
#test bayts with standard pdfs
s1vv_pdf <- c(c("gaussian","gaussian"),c(-10,1.5),c(-14,3))   
s1vh_pdf <- c(c("gaussian","gaussian"),c(-15,1.5),c(-20,3))  
lndvi_pdf <- c(c("gaussian","gaussian"),c(0.7,0.075),c(0.4,0.125))

# (1b) Theshold of deforestation probability at which flagged change is confirmed (chi)
chi = 0.9
# (1c) Start date of monitoring
start = 2016
# (1d)
bwf=c(0.1,0.9)

# (2) apply bayts (combine original time series into a time series of NF probabilities and detect deforestation)
# (2a) apply bayts
bts <- bayts(tsL=list(tsndvi,tsSvv,tsSvh),pdfL=list(lndvi_pdf,s1vv_pdf,s1vh_pdf),chi=chi,bwf=bwf,start=start)


# (2b) plot original time series; including flagged and detected changes
plotBayts(bts$bayts,start=start)
#plotBayts(bts$bayts,start=start,ylimL=list(c(0,1),c(-15,0),c(-25,0)),labL=list("S1 VV [dB]","S1 VH [dB]","Landsat NDVI"))


# Test the baytsSpatial on a subset
btss <- baytsSpatial(bL = list(ndvi,s1vv,s1vh), datesL = list(ndvi_date,s1vv_date,s1vh_date),
                     pdfL = list(lndvi_pdf,s1vv_pdf,s1vh_pdf), bwf = c(0.1, 0.9), chi = 0.9, PNFmin = 0.5,
                     start = 2016, end = NULL, pptype = "irregular", out_file = all_input_test,
                     mc.cores = 8)



