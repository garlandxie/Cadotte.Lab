#### Metadata ####

## metric.250: landscape metrics for each nestbox site 
# field_1 (chr): nestbox site
# class (int): landscape cover type; 1-tree, 2-grass/shrubs, 3-bare earth, 4-water, 5-buildings, 6-roads, 7-agriculture
# n.patches (num): number of patches of a landscape cover type for each site
# prop.lanscape (num): proportion of landscape cover type for each site
# edge.density (num): landscape heterogeneity of a cover type for each site 
# mean.patch.area (num): mean patch area of a cover type for each site
# max.patch.area (num): maximum patch area of a cover type for each site 


### Libraries ####
library(ade4)
library(sp)
library(raster)
library(rgeos)
library(spatialEco)
library(rgdal)
library(reshape2)
library(data.table)
library(PerformanceAnalytics)
library(maptools)
library(spdep)
library(pscl)
library(ggplot2)
library(ggrepel)

### Load Libraries ####
load("/Users/garlandxie/Google Drive/UTSC/UTSC 2016-2017/BIOD99/D99_Data/D99_RData/D99_Phylo_Land%.RData")

# Import #

# 2011-2013 community matrix (raw abundances)
comm <- read.csv("/Users/garlandxie/Google Drive/UTSC/UTSC 2016-2017/BIOD99/D99_Comm_Matrices/D99 _SiteXSpecies_2011-2013.csv")

# MFD
funcMPD <- read.csv("/Users/garlandxie/Google Drive/UTSC/UTSC 2016-2017/BIOD99/functMPD_BEEFD_Feb2017.csv")
funcMPD <- data.table(funcMPD)

# pairwise distance matrix
func.dist <- read.csv("/Users/garlandxie/Google Drive/UTSC/UTSC 2016-2017/BIOD99/D99_Data/D99_Traits/beeFunctDistMatrix.csv")
rownames(func.dist) <- func.dist$X
func.dist <- func.dist[,-c(1)] # remove X
func.dist <- as.matrix(func.dist)

# SppXtraits matrix
spp_tr <- read.csv("/Users/garlandxie/Google Drive/UTSC/UTSC 2016-2017/BIOD99/D99_Traits/D99_Spp_Traits.csv", 
                   stringsAsFactors = F)
rownames(spp_tr) <- spp_tr$Spp

### Landscape Metric Calculations ####

# Check to see if you have the proper drivers
# Should have ESRI Shapefiles!
ogrDrivers()$name

# Set working directory 
setwd("~/Google Drive/UTSC/UTSC 2016-2017/BIOD99/D99_Spatial_Data/JSM_LatLongs/D99_Nestbox_Buffers")

# Read 250 buffers of nestboxes as a SpatialPolygonsDataFrame obj
nest_250 <- readOGR(dsn = "D99_Buffers_250",
                    layer = "D99_Buffers_250_Removed")

colnames(nest_250@data) <- c("Site", "Lat", "Long")

# Read 500 buffers of nestboxes as a SpatialPolygonsDataFrame obj
nest_500 <- readOGR(dsn = "D99_Buffers_500", 
                    layer = "D99_Buffer_500_Updated")

colnames(nest_500@data) <- c("Site", "Lat", "Long")

# Read 2007 land cover data in raster class
lc_2007 <- raster("/Users/garlandxie/Google Drive/UTSC/UTSC 2016-2017/BIOD99/D99_Spatial_Data/D99_TO_Land_Cover/toronto_2007_landcover.img")

# Plot to make sure it looks reasonable
plot(nest_250)
plot(nest_500) # buffers are larger than nest_250, makes sense!

# Redefine projection for vector shapefiles; reason: rasterfile reprojection takes too long. 
# Caution: redefining CRS for vectors may be risky!!
nest_250 <- spTransform(nest_250, 
                        CRS("+proj=tmerc +lat_0=0 +lon_0=-79.5 +k=0.9999 +x_0=304800 +y_0=0 
                            +datum=NAD27 +units=m +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"))

nest_500 <- spTransform(nest_500, 
                        CRS("+proj=tmerc +lat_0=0 +lon_0=-79.5 +k=0.9999 +x_0=304800 +y_0=0 
                            +datum=NAD27 +units=m +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"))

# Check to make sure shapefile overlaps raster
# They overlap! CRS for both vector and raster files are consistent
plot(lc_2007); plot(nest_250, add = TRUE)
plot(lc_2007); plot(nest_500, add = TRUE)

# Calculate Landscape metrics
# metrics: 2 = n.patches, 4 = prop.landscapes, 7 = edge.density, 10 = mean.patch.area, 13 = mean.shape.index, 
l_met_250 <- land.metrics(x = nest_250, y = lc_2007, metrics = c(2, 4, 7, 10, 13))
l_met_500 <- land.metrics(x = nest_500, y = lc_2007, metrics = c(2, 4, 7, 10, 13))



##### Functional diversity metrics ######

ses.mfd.taxa <- ses.mpd(samp = comm.total[ ,-c(1,2,3,4)], 
                        # remove site, community type, lat, long
                        # pairwise functional distance (Gower's distance)
                        dis = func.dist,
                        # use tip-label swapping to preserve phenotypes
                        null.model = "taxa", 
                        # use abundance-weighed ses.MFD values
                        abundance.weighted = T, 
                        runs = 999, 
                        iterations = 1000)

# change rownames from ses.mpd to ses.mfd 
colnames(ses.mfd.taxa) <- c("ntaxa", "mfd.obs", "mfd.rand.mean", "mdf.rand.sd", "mfd.obs.rank", "mfd.obs.z", "mpd.obs.p", "runs")

# add in site, lat, longs from comm.total
ses.mfd.taxa$Site <- comm.total$Site
ses.mfd.taxa$Lat <- comm.total$Lat
ses.mfd.taxa$Long <- comm.total$Long

# merge with environmental variables using lat-longs as common key
final <- merge(ses.mfd.taxa, final, by = c("Lat", "Long"))

### Exploratory Data Analysis ####

## Open RData files 

load("/Users/garlandxie/Google Drive/UTSC/UTSC 2016-2017/BIOD99/D99_Spatial_Data/D99_LandMetrics/D99_R_Land%/metric.250.RData")

## Convert into data.table structure

metric.250 <- data.table(metric.250)

## Reclassify data from for landscape metrics
metric.250$class <- as.integer(metric.250$class)
metric.250$n.patches <- as.numeric(metric.250$n.patches)
metric.250$prop.landscape <- as.numeric(metric.250$prop.landscape)
metric.250$edge.density <- as.numeric(metric.250$edge.density)
metric.250$mean.patch.area <- as.numeric(metric.250$mean.patch.area)
metric.250$max.patch.area <- as.numeric(metric.250$max.patch.area)

## Change column names
colnames(metric.250) <- c("Site", "class", "n.patches", "prop.landscape", 
                          "edge.density", "mean.patch.area", "max.patch.area")

## 250 spatial scale: histograms for landscape composition
# Variables of interests are: trees(1), grasses(2), buildings(5), roads(6), agriculture(7)

hist(as.numeric(metric.250[metric.250$class == 1,]$prop.landscape), 
     main = "250m Spatial Scale",
     xlab = "% Tree Canopy")

hist(as.numeric(metric.250[metric.250$class == 2,]$prop.landscape), 
     main = "250m Spatial Scale",
     xlab = "% Grasslands")

hist(as.numeric(metric.250[metric.250$class == 5,]$prop.landscape), 
     main = "250m Spatial Scale",
     xlab = "% Buildings")

hist(as.numeric(metric.250[metric.250$class == 6,]$prop.landscape), 
     main = "250m Spatial Scale",
     xlab = "% Roads")

hist(as.numeric(metric.250[metric.250$class == 7,]$prop.landscape), 
     main = "250m Spatial Scale",
     xlab = "% Other Paved Surfaces")

## 250 spatial scale: histograms for landscape heterogeneity

hist(metric.250[metric.250$class == 1,]$edge.density, 
     main = "250m Spatial Scale",
     xlab = "Edge Density (Trees)")

hist(metric.250[metric.250$class == 2,]$edge.density, 
     main = "250m Spatial Scale",
     xlab = "Edge Density (Grasses)")

hist(metric.250[metric.250$class == 5,]$edge.density, 
     main = "250m Spatial Scale",
     xlab = "Edge Density (Buildings)")

hist(as.numeric(metric.250[metric.250$class == 6,]$edge.density), 
     main = "250m Spatial Scale",
     xlab = "% Roads")

hist(as.numeric(metric.250[metric.250$class == 7,]$prop.landscape), 
     main = "250m Spatial Scale",
     xlab = "% Other Paved Surfaces")

### Create database ####

## Aggregate prop.landscape and sum.ED by site to determine any errors
metric.250.agg <- metric.250[ ,.(sum_250 = sum(prop.landscape, na.rm = T), 
                                ED_250 = sum(edge.density, na.rm = T)), 
                                by = Site]

# Ideally, every aggregate measure (sum) should be equal to 1
# var() is an approximate way to check the equality of all rows 
var(metric.250.agg$sum) # 2.793186e-31, close enough to zero

# Note: compared lowest sum.ED and highest sum.ED values, visualised landscape patterns in QGIS, and data makes sense
# histogram to see any landscape heterogeneity
hist(metric.250.agg$sum.ED, 
     main = "250m Spatial Scale",
     xlab = "Edge Density")

# Reshape metric.250: class changes rows to columns with prop of landscape
wide_250 <- dcast(metric.250, Site ~ class, value.var = "prop.landscape") 
colnames(wide_250) <- c("Site", "Tree_250", "Grass_250", 
                        "Earth_250", "Water_250", "Buildings_250", 
                        "Roads_250", "OtherGrey_250", "Agriculture_250")
final_250 <- merge(metric.250.agg, wide_250, by = "Site")

## Add spatial coordinates into final_250

# Import March latlongs
jsm_latlongs <- read.csv("/Users/garlandxie/Google Drive/UTSC/UTSC 2016-2017/BIOD99/D99_Spatial_Data/JSM_LatLongs/JSM_LatLongs_Apr_14_2015.csv",
                         header = F)

# Rename columns
colnames(jsm_latlongs) <- c("Site", "Long", "Long")

# Merge final_250 and jsm_latlongs using the common key "Site"
final_250 <- merge(jsm_latlongs, final_250, by = "Site")

# Convert NA's into 0's for prop.landscape classes
final_250[is.na(final_250)] <- 0

# Create new column that represents proportion of urban area
final_250$Urb_250 <- apply(final_250[, c("Buildings_250", "Roads_250", "OtherGrey_250")], 1, FUN = sum)

### Covariance Matrix: 250m spatial scale #####

## Perform covariance matrix using Pearson's product-moment correlation coefficient for 250m spatial scale
cor(final_250[,c("ED_250","Tree_250", "Grass_250", "Buildings_250", "Roads_250", "OtherGrey_250", "Urb_250")], 
    method = "spearman")


#                   ED_250   Tree_250  Grass_250   Buildings_250  Roads_250  OtherGrey_250  Urb_250
#ED_250         1.00000000 -0.4803650  0.04255109     0.7149427  0.7506984    0.14808246  0.6340149
#Tree_250      -0.48036501  1.0000000 -0.15335307    -0.7412585 -0.6875961   -0.74610210 -0.8734628
#Grass_250      0.04255109 -0.1533531  1.00000000    -0.2314527 -0.1683065    0.06821473 -0.1699437
#Buildings_250  0.71494268 -0.7412585 -0.23145271     1.0000000  0.8101168    0.47738631  0.9244671
#Roads_250      0.75069843 -0.6875961 -0.16830647     0.8101168  1.0000000    0.41262111  0.8475750
#OtherGrey_250  0.14808246 -0.7461021  0.06821473     0.4773863  0.4126211    1.00000000  0.7138082
#Urb_250        0.63401491 -0.8734628 -0.16994369     0.9244671  0.8475750    0.71380817  1.0000000

# Comments: 
# ED_250 and Urb_250 are moderately correlated (0.63); possible set of candidate predictor variables
# ED_250, Tree_250, and Grass_250 are moderately correlated with each other; possible set of candidate predictor variables
# histograms of ED_250 and Urb_250 shows a large variance in landscape metrics; some evidence of a gradient 

# Plot covariance matrix with histograms, scatterpolots, coefficient values, and p-values (0.5, 0.1*, 0.05*, 0.01** etc.)
(chart.Correlation(final_250[,c("ED_250","Tree_250", "Grass_250", "Buildings_250", "Roads_250", "OtherGrey_250", "Urb_250")], 
                   method = "spearman", cex.label = 4.0))


### Covariance Matrix: 500m spatial scale #####

## Perform covariance matrix using Pearson's product-moment correlation coefficient for 500m spatial scale
cor(final_500[,c("ED_500","Tree_500", "Grass_500", "Buildings_500", "Roads_500", "OtherGrey_500", "Urb_500")], 
    method = "spearman")

#                  ED_500   Tree_500   Grass_500  Buildings_500  Roads_500 OtherGrey_500    Urb_500
#ED_500         1.00000000 -0.4859612  0.01382977     0.7626842  0.7724547     0.1816956  0.6497225
#Tree_500      -0.48596117  1.0000000 -0.12957894    -0.7651227 -0.6981510    -0.7898522 -0.8865659
#Grass_500      0.01382977 -0.1295789  1.00000000    -0.2286502 -0.1632329     0.1186228 -0.1529913
#Buildings_500  0.76268420 -0.7651227 -0.22865025     1.0000000  0.8587741     0.5086522  0.9347781
#Roads_500      0.77245467 -0.6981510 -0.16323289     0.8587741  1.0000000     0.4015620  0.8663279
#OtherGrey_500  0.18169564 -0.7898522  0.11862279     0.5086522  0.4015620     1.0000000  0.7234642
#Urb_500        0.64972248 -0.8865659 -0.15299134     0.9347781  0.8663279     0.7234642  1.0000000


# Comments: 
# ED_250 and Urb_250 are moderately correlated (0.65); possible set of candidate predictor variables
# ED_250, Tree_250, and Grass_250 are moderately correlated with each other; possible set of candidate predictor variables
# histograms of ED_250 and Urb_250 shows a large variance in landscape metrics; some evidence of a gradient 

# Plot covariance matrix with histograms, scatterpolots, coefficient values, and p-values (0.5, 0.1*, 0.05*, 0.01** etc.)
(chart.Correlation(final_500[,c("ED_500","Tree_500", "Grass_500", "Buildings_500", "Roads_500", "OtherGrey_500", "Urb_500")], 
                   method = "kendall", cex.label = 4.0))

### Spatial Autocorrelation ####

# merge both 500m and 250m landscape metrics together 
final <- merge(final_500, final_250, by = c("Lat", "Long"))

# columns sum_250 and sum_500 can be removed using data.table syntax
final <- data.table(final)
final <- final[ ,":="(sum_250 = NULL, sum_500 = NULL)] 

# merge ses.mpd null models (taxa.labels) to landscape metrics
final <- merge(ses.mpd.taxa, final, by = c("Lat", "Long"))

# remove any unncessary columns from final using data.table syntax
ses.mpd.taxa <- data.table(ses.mpd.taxa)
final <- final[ ,":="(mpd.rand.mean = NULL, mpd.rand.sd = NULL, mpd.obs.rank = NULL, 
                      mpd.obs.p = NULL, runs = NULL, Site.x = NULL)] 

# Calculate distance-based neighbours: use k-nearest neighbour algorithm
# Binary spatial weigh matrices don't account for complex intersite relationships
# But how to account for 
k1_nb <- knn2nb(knearneigh(cbind(final$Lat,final$Long), k = 8), row.names = final$Site)

# Calculate spatial weight matrix
swm <- nb2listw(k1_nb, glist = NULL, style="W", zero.policy = NULL)

# Calculate global Moran's I index using a Monte Carlo simulation (permutation)
# replace first argument with a dependent variable
# number of simulations is 999
# set zero.policy to TRUE to accound for no-neighbour observations
moran.mc(final$ntaxa, swm, na.action = na.omit, nsim = 999, zero.policy = TRUE)


### Monte-Carlo simulation of Moran I

## data:  final$mpd.obs 
# weights: swm 
# omitted: 1, 5, 6, 7, 9, 10, 11, 12, 14, 23, 28, 30, 31, 32, 33, 38, 44, 52, 53, 54, 60, 63, 66, 67, 68, 69, 72, 75, 76, 77, 79, 81, 84, 86, 87, 88, 103, 105, 106, 110, 112, 115, 116, 117, 120, 122, 123, 128, 129, 131, 133, 138, 140, 143, 144, 149, 150, 151, 153, 156, 165 
# number of simulations + 1: 1000 
# statistic = 0.039788, observed rank = 809, p-value = 0.191
# alternative hypothesis: greater

## data:  final$mpd.obs.z 
# weights: swm 
# omitted: 1, 5, 6, 7, 9, 10, 11, 12, 14, 23, 28, 30, 31, 32, 33, 38, 44, 52, 53, 54, 60, 63, 66, 67, 68, 69, 72, 75, 76, 77, 79, 81, 84, 86, 87, 88, 103, 105, 106, 110, 112, 115, 116, 117, 120, 122, 123, 128, 129, 131, 133, 138, 140, 143, 144, 149, 150, 151, 153, 156, 165 
# number of simulations + 1: 1000 
# statistic = 0.041466, observed rank = 790, p-value = 0.21
# alternative hypothesis: greater

## Monte-Carlo simulation of Moran I

# data:  final$ntaxa 
# weights: swm  
# number of simulations + 1: 1000 
# statistic = 0.0091759, observed rank = 701, p-value = 0.299
# alternative hypothesis: greater


### GLM - Phylogenetic Structure #####
mpd_250_gls <- gls(mpd.obs.z ~ Urb_250 + ED_250 , data = final, na.action = na.omit)

# semivariogram: high nugget effect, weak or no spatial autocorrelation
plot(Variogram(mpd_250_0, form=~ final[!is.na(final$mpd.obs.z), ]$Lat + final[!is.na(final$mpd.obs.z), ]$Long))

mpdz_250_lm0 <- lm(mpd.obs.z ~ Urb_250, data = final)
mpdz_250_lm <- lm(mpd.obs.z ~ Urb_250 + ED_250, data = final)

# Call:
#  lm(formula = mpd.obs.z ~ Urb_250 + ED_250, data = final)

# Residuals:
#  Min       1Q   Median       3Q      Max 
# -1.94022 -0.65474 -0.03254  0.48716  1.74968 

# Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.9682     0.2248  -4.307 3.78e-05 ***
#  Urb_250       1.1850     0.5260   2.253   0.0264 *  
#  ED_250       -0.3903     0.7650  -0.510   0.6110    
---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.8053 on 103 degrees of freedom
# (63 observations deleted due to missingness)
# Multiple R-squared:  0.06275,	Adjusted R-squared:  0.04455 
# F-statistic: 3.448 on 2 and 103 DF,  p-value: 0.03553


mpd_500_gls <- gls(mpd.obs.z ~ Urb_500 + ED_500, data = final, na.action = na.omit)

# testing for spatial autocorrelation using semivariogram
# high nugget indicates a weak or no spatial autocorrelation
plot(Variogram(mpd_500_gls, form =~ final[!is.na(final$mpd.obs.z), ]$Lat + final[!is.na(final$mpd.obs.z), ]$Long))

mpdz_500_lm0 <- lm(mpd.obs.z ~ Urb_500, data = final)
mpdz_lm_500 <- lm(mpd.obs.z ~ Urb_500 + ED_500, data = final)

# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     -1.0328     0.1987  -5.198 1.04e-06 ***
#    Urb_500       0.9603     0.4153   2.313   0.0228 *  
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Null deviance: 69.968  on 103  degrees of freedom
# Residual deviance: 66.482  on 102  degrees of freedom
# (61 observations deleted due to missingness)
# AIC: 254.6

mpd_lm_250 <- lm(mpd.obs ~ Urb_250 + ED_250, data = final)

# Call:
#  lm(formula = mpd.obs ~ Urb_250 + ED_250, data = final)

# Residuals:
#  Min      1Q  Median      3Q     Max 
# -0.4390 -0.1563  0.0063  0.1327  0.5384 

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.3332     0.0620   5.373 4.83e-07 ***
#  Urb_250       0.4231     0.1451   2.916  0.00435 ** 
#  ED_250       -0.1946     0.2110  -0.922  0.35873    
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.2221 on 103 degrees of freedom
# (63 observations deleted due to missingness)
# Multiple R-squared:  0.09208,	Adjusted R-squared:  0.07445 
# F-statistic: 5.223 on 2 and 103 DF,  p-value: 0.006909


mpd_lm_500 <- lm(mpd.obs ~ Urb_500 + ED_500, data = final)

# Call:
#  lm(formula = mpd.obs ~ Urb_500 + ED_500, data = final)

# Residuals:
#  Min       1Q   Median       3Q      Max 
# -0.45531 -0.18100  0.01118  0.13630  0.53462 

# Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.29504    0.06893   4.280 4.19e-05 ***
#  Urb_500      0.31894    0.15388   2.073   0.0407 *  
#  ED_500       0.01381    0.24085   0.057   0.9544    
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.2246 on 103 degrees of freedom
# (63 observations deleted due to missingness)
# Multiple R-squared:  0.07184,	Adjusted R-squared:  0.05382 
# F-statistic: 3.986 on 2 and 103 DF,  p-value: 0.02151



#### GLM - Species Richness #####

# GLS can account for some correlation between residuals 
SR_250_gls <- gls(ntaxa ~ ED_250 + Urb_250, data = final)

# semivariogram: high nugget effect, weak or no spatial autocorrelation
plot(Variogram(SR_250_lm, form = ~final[!is.na(final$mpd.obs.z), ]$Lat + final[!is.na(final$mpd.obs.z), ]$Long))

SR_250_glm0 <- glm.nb(ntaxa ~ ED_250, data = final)
SR_250_glm <- glm.nb(ntaxa ~ ED_250 + Urb_250, data = final)

# Coefficients:
#                Estimate Std. Error z value Pr(>|z|)   
#  (Intercept)  0.52431    0.17974   2.917  0.00353 **
#  ED_250       1.16025    0.56243   2.063  0.03912 * 
#  Urb_250     -0.09745    0.36412  -0.268  0.78899   
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Null deviance: 201.01  on 164  degrees of freedom
# Residual deviance: 195.46  on 162  degrees of freedom
# AIC: 662.22

SR_500_gls <-  gls(ntaxa ~ ED_500 + Urb_500, data = final, na.action = na.omit)

# semivariogram: high nugget effect, weak or no spatial autocorrelation
plot(Variogram(SR_500_gls, form = ~final[!is.na(final$mpd.obs.z), ]$Lat + final[!is.na(final$mpd.obs.z), ]$Long))

SR_500_glm0 <- glm.nb(ntaxa ~ ED_500, data = final)
SR_500_glm <- glm.nb(ntaxa ~ ED_500 + Urb_500, data = final)

# Call:
#   glm.nb(formula = ntaxa ~ ED_500, data = final, init.theta = 3.721727164, link = log)

# Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
# -2.0722  -0.9732  -0.1050   0.5957   2.3998  

#  Coefficients:
#  Estimate Std. Error z value Pr(>|z|)  
#  (Intercept)   0.4969     0.1998   2.487   0.0129 *
#  ED_500        1.1365     0.5304   2.143   0.0321 *
#   ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# (Dispersion parameter for Negative Binomial(3.7217) family taken to be 1)

# Null deviance: 200.33  on 164  degrees of freedom
# Residual deviance: 195.68  on 163  degrees of freedom
# AIC: 661.1

# Number of Fisher Scoring iterations: 1


### GLM - Functional Structure #####

mfdz_lm0_500 <- lm(mfd.obs.z ~ Urb_500, data = final)

# Call:
#  lm(formula = mfd.obs.z ~ Urb_500, data = final)

# Residuals:
#  Min      1Q  Median      3Q     Max 
# -2.4070 -0.7329  0.1308  0.6396  1.9792 

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)  -0.9467     0.2228  -4.249 4.69e-05 ***
#  Urb_500       0.6837     0.4639   1.474    0.144    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.9056 on 104 degrees of freedom
# (63 observations deleted due to missingness)
# Multiple R-squared:  0.02046,	Adjusted R-squared:  0.01104 
# F-statistic: 2.172 on 1 and 104 DF,  p-value: 0.1436

mfdz_lm_500 <- lm(mfd.obs.z ~ Urb_500 + ED_500, data = final)

# Call:
#  lm(formula = mfd.obs.z ~ Urb_500 + ED_500, data = final)

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.31270 -0.74468  0.07519  0.54486  2.08426 

# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)   
#   (Intercept)  -0.7343     0.2771  -2.650  0.00931 **
#   Urb_500       1.2104     0.6185   1.957  0.05306 . 
#   ED_500       -1.2415     0.9681  -1.282  0.20257   
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.9028 on 103 degrees of freedom
# (63 observations deleted due to missingness)
# Multiple R-squared:  0.03585,	Adjusted R-squared:  0.01713 
# F-statistic: 1.915 on 2 and 103 DF,  p-value: 0.1525

mfdz_lm0_250 <- lm(mfd.obs.z ~ Urb_250, data = final) 

# Call:
#  lm(formula = mfd.obs.z ~ Urb_250, data = final)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.3969 -0.7970  0.1182  0.6037  1.9048 

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -1.0025     0.2037  -4.921 3.24e-06 ***
#   Urb_250       0.8459     0.4355   1.942   0.0548 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.8988 on 104 degrees of freedom
# (63 observations deleted due to missingness)
# Multiple R-squared:  0.03501,	Adjusted R-squared:  0.02573 
# F-statistic: 3.773 on 1 and 104 DF,  p-value: 0.05478

mfdz_lm_250 <- lm(mfd.obs.z ~ ED_250 + Urb_250, data = final)  

# Call:
# lm(formula = mfd.obs.z ~ ED_250 + Urb_250, data = final)

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.25496 -0.78860  0.04125  0.57952  2.00045 

# Coefficients:
#               Estimate    Std. Error  t value Pr(>|t|)   
# (Intercept)   -0.8089      0.2500  -3.236  0.00163 **
#  ED_250       -1.1290     0.8508  -1.327  0.18745   
#  Urb_250       1.3665     0.5849   2.336  0.02142 * 
#  ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.8956 on 103 degrees of freedom
# (63 observations deleted due to missingness)
# Multiple R-squared:  0.05123,	Adjusted R-squared:  0.03281 
# F-statistic: 2.781 on 2 and 103 DF,  p-value: 0.06664

mfd_lm_250 <- lm(mfd.obs ~ ED_250 + Urb_250, data = final)  

# Call:
#  lm(formula = mfd.obs ~ ED_250 + Urb_250, data = final)

# Residuals:
#  Min        1Q    Median        3Q       Max 
# -0.171756 -0.051289 -0.002831  0.058697  0.163767 

# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)  0.13023    0.02016   6.461  3.5e-09 ***
#   ED_250      -0.12056    0.06860  -1.757   0.0818 .  
#   Urb_250      0.14683    0.04716   3.113   0.0024 ** 
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.07221 on 103 degrees of freedom
# (63 observations deleted due to missingness)
# Multiple R-squared:  0.08761,	Adjusted R-squared:  0.06989 
# F-statistic: 4.945 on 2 and 103 DF,  p-value: 0.008898

mfd_lm0_500 <- lm(mfd.obs ~ Urb_500, data = final)  


mfd_lm_500 <- lm(mfd.obs ~ ED_500 + Urb_500, data = final)  

# Call:
#  lm(formula = mfd.obs ~ ED_500 + Urb_500, data = final)

# Residuals:
#  Min        1Q    Median        3Q       Max 
# -0.161273 -0.049099  0.000927  0.052296  0.181238 

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.12188    0.02270   5.370  4.9e-07 ***
#  ED_500      -0.04817    0.07931  -0.607   0.5449    
#  Urb_500      0.09858    0.05067   1.946   0.0544 .  
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.07395 on 103 degrees of freedom
# (63 observations deleted due to missingness)
# Multiple R-squared:  0.04296,	Adjusted R-squared:  0.02437 
# F-statistic: 2.312 on 2 and 103 DF,  p-value: 0.1042



### Phylogenetic Signals ####

# phylogenetic signal: strength of covariance between differences in trait values and phylo-distance separating species
# strong signal: closely-related species share similar traits 

# pairwise correlation cofficient matrix: determine independent or dependent traits 
# remove these columns; 1: Spp, 2: Group, 3: Native Status

cor(spp_tr[, -c(1, 2, 3)], use = "complete.obs", method = "spearman")
chart.Correlation(spp_tr[, -c(1, 2, 3)], use = "complete.obs", method = "spearman")

#             Emer.X     Nest.Mat   Pol.Col     Lecty       Volt        M.Bl       IT
# Emer.X    1.0000000 -0.3863054  0.42529344  0.28178446 -0.2483116  0.38001838  0.2946802
# Nest.Mat -0.3863054  1.0000000 -0.81012928 -0.10040323  0.4897804 -0.64578170 -0.4949972
# Pol.Col   0.4252934 -0.8101293  1.00000000  0.09295112 -0.7313355  0.60532133  0.3638485
# Lecty     0.2817845 -0.1004032  0.09295112  1.00000000 -0.1149932 -0.05897173 -0.3340716
# Volt     -0.2483116  0.4897804 -0.73133551 -0.11499323  1.0000000 -0.70371373 -0.4329231
# M.Bl      0.3800184 -0.6457817  0.60532133 -0.05897173 -0.7037137  1.00000000  0.8338565
# IT        0.2946802 -0.4949972  0.36384850 -0.33407164 -0.4329231  0.83385647  1.0000000

# high correlation: M.Bl + IT, Volt + Pol.Col (|r| > 0.7), 
# moderate correlation: M.Bl + Nest.Mat,  (0.5 < |r| < 0.7)
# low correlation: Emer.X + Nest.Mat, Emer.X + Pol.Col, Emer.X + Volt, Emer.X + M.Bl, Emer.X + IT (|r| < 0.5)

# Conclusion: remove mean body length, and keep the raw data of the variables for phylogenetic tests. 

## Pagel's lambda: Measuring and testing phylogenetic signal
# Munkemuller et al. 2012 shows that Pagel's lambda is the most robust measure in terms of size, polytomies and branch length 

## star phylogeny of current tree 
lambda0 <-rescale(tree.abv, "lambda", 0)

## IT 
spp_tr <- spp_tr[tree.abv$tip.label, ]
IT <- spp_tr[complete.cases(spp_tr), c("IT")] # drop Coel.alt, Col.moe, Coel.say, Stel.var
names(IT) <- rownames(spp_tr[complete.cases(spp_tr), ])

# use MLE approach to estimate lambda for IT traits under a resolved ultrametric tree 
# fitContinuous is used for continuous functional traits
IT_lambda <- fitContinuous(dat = IT, 
                           phy = drop.tip(tree.abv, tip = c("Hyl.vert", "Hyl.mes", "Meg.rel", "Meg.iner")), 
                           model = "lambda")

# estimate lambda for IT traits under a star phylogeny (lambda0)
IT_lambda0 <- fitContinuous(dat = IT, 
                            phy = drop.tip(lambda0, tip = c("Hyl.vert", "Hyl.mes", "Meg.rel", "Meg.iner")), 
                            model = "lambda")

# significance test: p < 0.01, lambda: 1, lambda0: 
LR.IT <- 2*(IT_lambda$opt$lnL - IT_lambda0$opt$lnL)
pchisq(LR.IT, df = 1, lower.tail = F)


## Nest.Mat 
# drop polytomies to form a fully bifurcating tree 
spp_tr <- spp_tr[tree.abv$tip.label, ]
Nest.Mat <- spp_tr[ , c("Nest.Mat")]
names(Nest.Mat) <- rownames(spp_tr)

Nest_lambda <- fitDiscrete(dat = Nest.Mat, 
                            phy = drop.tip(tree.abv, tip = c("Hyl.vert", "Hyl.mes", "Stel.var", "Meg.rel", "Meg.iner")), 
                            transform = "lambda")

Nest_lambda0 <- fitDiscrete(dat = Nest.Mat, 
                            phy = drop.tip(lambda0, tip = c("Hyl.vert", "Hyl.mes", "Meg.rel", "Meg.iner")), 
                            transform = "lambda")

# significant p < 0.01; lambda0: 0, lambda: 1 (caution: restricted to 1 b/c of bound limits)
LR.Nest <- 2*(Nest_lambda$opt$lnL - Nest_lambda0$opt$lnL)
pchisq(LR.Nest, df = 1, lower.tail = F)

# Lecty ()
# drop polytomies to form a fully bifurcating tree 
spp_tr <- spp_tr[tree.abv$tip.label, ]
Lecty <- spp_tr[, c("Lecty")]
names(Lecty) <- rownames(spp_tr)
Lecty_lambda <- fitDiscrete(dat = Lecty, 
                            phy = drop.tip(tree.abv, tip = c("Hyl.vert", "Hyl.mes", "Meg.rel", "Meg.iner")), 
                            transform = "lambda")

Lecty_lambda0 <- fitDiscrete(dat = Lecty, 
                            phy = drop.tip(lambda0, tip = c("Hyl.vert", "Hyl.mes", "Meg.rel", "Meg.iner")), 
                            transform = "lambda")

# significance test: p < 0.05, lambda: 0.911, lambda0: 0.25
LR.Lecty <- 2*(Lecty_lambda$opt$lnL - Lecty_lambda0$opt$lnL)
pchisq(LR.Lecty, df = 1, lower.tail = F)

# Pol.Col 
spp_tr <- spp_tr[tree.abv$tip.label, ]
Polcol <- spp_tr[, c("Pol.Col")]
names(Polcol) <- rownames(spp_tr)

Pol_lambda <- fitDiscrete(dat = Polcol, 
                            phy = drop.tip(tree.abv, tip = c("Hyl.vert", "Hyl.mes", "Meg.rel", "Meg.iner")), 
                            transform = "lambda")

Pol_lambda0 <- fitDiscrete(dat = Polcol, 
                             phy = drop.tip(lambda0, tip = c("Hyl.vert", "Hyl.mes", "Meg.rel", "Meg.iner")), 
                             transform = "lambda")

# significance test: p < 0.01, lambda: 1, lambda0: 0.26
LR.pol <- 2*(Pol_lambda$opt$lnL - Pol_lambda0$opt$lnL)
pchisq(LR.pol, df = 1, lower.tail = F)

# Emer.X 
spp_tr <- spp_tr[tree.abv$tip.label, ]
Emer <- spp_tr[ , c("Emer.X")]
names(Emer) <- rownames(spp_tr)

Emer_lambda <- fitDiscrete(dat = Emer, 
                          phy = drop.tip(tree.abv, tip = c("Hyl.vert", "Hyl.mes", "Meg.rel", "Meg.iner")), 
                          transform = "lambda")

Emer_lambda0 <- fitDiscrete(dat = Emer, 
                           phy = drop.tip(lambda0, tip = c("Hyl.vert", "Hyl.mes", "Meg.rel", "Meg.iner")), 
                           transform = "lambda")

# significance test: p < 0.01, lambda: 0.84, lambda_null: 0
LR.Emer <- 2*(Emer_lambda$opt$lnL - Emer_lambda0$opt$lnL)
pchisq(LR.Emer, df = 1, lower.tail = F)

# Volt 
spp_tr <- spp_tr[tree.abv$tip.label, ]
Volt <- spp_tr[ , c("Volt")]
names(Volt) <- rownames(spp_tr)

Volt_lambda <- fitDiscrete(dat = Volt, 
                           phy = drop.tip(tree.abv, tip = c("Hyl.vert", "Hyl.mes", "Meg.rel", "Meg.iner")), 
                           transform = "lambda")

Volt_lambda0 <- fitDiscrete(dat = Volt, 
                            phy = drop.tip(lambda0, tip = c("Hyl.vert", "Hyl.mes", "Meg.rel", "Meg.iner")), 
                            transform = "lambda")

# significance test: p < 0.01, lambda: 0.83, lambda0: 0.15
LR.Volt <- 2*(Volt_lambda$opt$lnL - Volt_lambda0$opt$lnL)
pchisq(LR.Volt, df = 1, lower.tail = F)

# Female Body Length
spp_tr <- spp_tr[tree.abv$tip.label, ]
M.Bl <- spp_tr[ , c("M.Bl")]
names(M.Bl) <- rownames(spp_tr)

Mbl_lambda <- fitDiscrete(dat = M.Bl, 
                           phy = drop.tip(tree.abv, tip = c("Hyl.vert", "Hyl.mes", "Meg.rel", "Meg.iner")), 
                           transform = "lambda")

Mbl_lambda0 <- fitDiscrete(dat = M.Bl, 
                            phy = drop.tip(lambda0, tip = c("Hyl.vert", "Hyl.mes", "Meg.rel", "Meg.iner")), 
                            transform = "lambda")

# significance test: p < 0.01, lambda: 1, lambda0: 0
LR.Mbl <- 2*(Mbl_lambda$opt$lnL - Mbl_lambda0$opt$lnL)
pchisq(LR.Mbl, df = 1, lower.tail = F)

## Extended RLQ (Pavoine et al. 2011) #####

## Spatial Matrix 

# geographical coordinates of samples sites 
mxy <- data.frame(Site = final$Site, Long = final$Long, Lat = final$Lat)

# create neighbourhood matrix using k-nearest neighbour algorithm
neigh <- nb2neig(graph2nb(gabrielneigh(cbind(mxy$Lat, mxy$Long))))

# plot samples & neighbourhood graph
# get rid of Site column (id: 1)
# make sure the map looks like the shape of Toronto 
s.label(dfxy = mxy[,c(-1)], neig = neigh,
        include.origin = FALSE, addaxes = FALSE,
        cpoint = 2, 
        cneig = 0.1)

## Species Abundance

# Make sure that comm.total has the same sites used as final
# use which() to determine how many sites in comm are in final 
# Remove columns 1 to 4 
comm <- comm[which(comm$Site %in% sort(final$Site)), -c(1:4, 36:39)]
rownames(comm) <- final$Site

coa_spp <- dudi.coa(comm, scan = FALSE, nf = 55)

vecspa <- scores.neig(neigh)
pcaspa <- dudi.pca(vecspa, coa_spp$lw, scan = F, nf = 96)

## Environmental variables 

# extract relevant columns
# perform log transformation on quantitative variables 
final <- final[order(final$Site), ]
env <- final[ , c("ED_500", "Urb_500", "ED_250", "Urb_250")]

# perform PCA 
pcaenv <- dudi.pca(env, coa_spp$lw, scan = FALSE, nf = 8)

## Test for spatial autocorrelation in environmental variables 
gm <- gearymoran(neig2mat(neigh), env) 

# class: krandtest 
# Monte-Carlo tests
# Call: as.krandtest(sim = matrix(res$result, ncol = nvar, byrow = TRUE), 
#                                 obs = res$obs, alter = alter, names = test.names)

# Number of tests:   4 

# Adjustment method for multiple comparisons:   none 
# Permutation number:   999 
# Test       Obs   Std.Obs   Alter Pvalue
# 1  ED_500 0.6157128 10.688177 greater  0.001
# 2 Urb_500 0.6889907 11.338082 greater  0.001
# 3  ED_250 0.4884916  8.477963 greater  0.001
# 4 Urb_250 0.5555339  9.275421 greater  0.001

## Traits

# Remove coelioxys and stelis from spp_tr
spp_tr <- spp_tr[-c(32:35), ]

rlq_Volt <- data.frame(spp_tr$Volt)
rlq_NestMat <- data.frame(spp_tr$Nest.Mat)
rlq_PolCol <- data.frame(spp_tr$Pol.Col)
rlq_Emer <- data.frame(spp_tr$Emer.X)
rlq_Lecty <- data.frame(spp_tr$Lecty)
rlq_M_Bl <- data.frame(spp_tr$M.Bl)
rlq_IT <- data.frame(spp_tr$IT)

# Assign names to each trait

names(rlq_Volt) <- "Voltinism"
names(rlq_NestMat) <- "Nesting Material"
names(rlq_PolCol) <- "Pollen Transportation"
names(rlq_Emer) <- "Emergence Time"
names(rlq_Lecty) <- "Diet Breadth"
names(rlq_M_Bl) <- "Female Body Length"
names(rlq_IT) <- "IT"

# metric for ordinal variables (1), and quantitative variables (1)
listdis <- ldist.ktab(ktab.list.df(list(rlq_Volt,
                                        rlq_NestMat,
                                        rlq_PolCol,
                                        rlq_Lecty,
                                        rlq_Emer,
                                        rlq_M_Bl,
                                        rlq_IT)),
                                        c("O", "O", "O", "O", "O", "Q", "Q"), 
                                        scan = TRUE)

## fourth-corner analysis

funtest <- function(dis){
  trlqO <- rlq(pcaenv, coa_spp, dudi.pco(dis, coa_spp$cw, full=T),
               scan = FALSE, nf = 8)
  funi <- function(i){
    e <- sample(1:ncol(comm.total[, abbrev]))
    return(sum(rlq(pcaenv, coa_spp, dudi.pco(as.dist(as.matrix(dis)[e, e]),
                                            coa_spp$cw, full = TRUE), scan = FALSE, nf = 8)$eig))
  }
  trlqS <- lapply(1:999, funi)
  return(as.randtest(unlist(trlqS), sum(trlqO$eig), alter="greater"))
} 




## RLQ and 4th Corner (Dray et al. 2014) #####

# correspondence analysis for species composition 
# remove coelioxys and stelis from analysis since they have missing IT values
comm <- comm[which(comm$Site %in% final$Site), -c(1:4 ,36:39)]
rownames(comm) <- final$Site
afcL.comm <- dudi.coa(comm, scannf = F)

# principal coordinate analysis for environmental variables 
final <- final[order(final$Site), ]
env <-  final[ , c("ED_500", "ED_250", 
                  "Urb_500", "Urb_250",
                  "Tree_500", "Tree_250",
                  "Grass_500", "Grass_250")]
rownames(env) <- final$Site
acpR.env <- dudi.pca(env, row.w = afcL.comm$lw, scannf = F)

# hill-smith analysis for traits (quantitative and qualitative)

# remove first three columns since they aren't considered traits
spp_tr <- spp_tr[complete.cases(spp_tr), -c(1:3)]

spp_tr$Emer.X <- as.numeric(spp_tr$Emer.X) 

spp_tr$Nest.Mat <- factor(spp_tr$Nest.Mat, 
                          levels = c(1, 2, 3, 4, 5, 7, 9), 
                          labels = c("Leaves", "Mud", "Leaves + Mud", 
                                     "Resin", "Leaves Chewed", 
                                     "None", "Wood Pulp"))

spp_tr$Pol.Col <- factor(spp_tr$Pol.Col, 
                        levels = c(1, 2, 3),
                         labels = c("Corbica", "Gut", "Scopa"))

spp_tr$Lecty <- as.numeric(spp_tr$Lecty)
spp_tr$Volt <- as.numeric(spp_tr$Volt)

acpQ <- dudi.hillsmith(spp_tr, row.w  = afcL.comm$cw, scannf = F)

# combining R, Q and L
rlq.bees <- rlq(acpR.env, afcL.comm, acpQ, scannf = F)

# Class: rlq dudi
# Call: rlq(dudiR = acpR.env, dudiL = afcL.comm, dudiQ = acpQ, scannf = F)

# Total inertia: 1.482

# Eigenvalues:
#  Ax1       Ax2       Ax3       Ax4       Ax5 
# 1.1461335 0.2836772 0.0458660 0.0056480 0.0004999 

# Projected inertia (%):
#  Ax1      Ax2      Ax3      Ax4      Ax5 
# 77.32725 19.13911  3.09448  0.38106  0.03373 

# Cumulative projected inertia (%):
#  Ax1   Ax1:2   Ax1:3   Ax1:4   Ax1:5 
# 77.33   96.47   99.56   99.94   99.98 

# (Only 5 dimensions (out of 8) are shown)

# Eigenvalues decomposition:
#  eig     covar      sdR      sdQ      corr
# 1 1.1461335 1.0705762 1.935624 1.511049 0.3660313
# 2 0.2836772 0.5326136 1.637784 1.516449 0.2144509

# Inertia & coinertia R (acpR.env):
#  inertia      max     ratio
# 1  3.746639 4.490630 0.8343236
# 12 6.428974 6.860909 0.9370441

# Inertia & coinertia Q (acpR.env):
#  inertia      max     ratio
# 1  2.283268 2.816649 0.8106329
# 12 4.582886 5.572057 0.8224766

# Correlation L (afcL.comm):
#   corr       max     ratio
# 1 0.3660313 0.7969012 0.4593183
# 2 0.2144509 0.7649655 0.2803406

## Correlations between axe scores and environmental variables

#                RS1        RS2
# ED_500    -0.52869088  0.1672816
# ED_250    -0.50866536  0.1831467
# Urb_500   -0.35072303 -0.2867327
# Urb_250   -0.25941590 -0.4477556
# Tree_500   0.13979954  0.4789756
# Tree_250   0.06758709  0.5872258
# Grass_500  0.38055785 -0.2285512
# Grass_250  0.32016566 -0.1711319



### Plots ####

# SR_250
plot(ntaxa ~ ED_250, data = final, 
     xlab = 'Edge Density (250m scale)', 
     ylab = "Species Richness", 
     col = "grey", lwd = 3, lty = 2, pch = 16)
abline(SR_250_glm, lwd = 2)

# SR_500 
plot(ntaxa ~ ED_500, data = final, 
     xlab = 'Edge Density (500m scale)', 
     ylab = "Species Richness", 
     col = "grey", lwd = 3, lty = 2, pch = 16)
abline(SR_500_glm, lwd = 2)

# sesMPD  vs Urb_250
termplot(mpdz_250_lm, 
         # include Urb_250 as the indepedent variable from a multi-regression model
         terms = "Urb_250", 
         # include partial residuals
         partial.resid = T, se = T, 
         # x and y labels, 
         xlab = "Proportion of Urban Area (250m)",
         ylab = "ses.MPD", 
         # use filled-in grey points
         pch = 16)
abline(h = 0, lty = 2)

# sesMPD  vs Urb_500
termplot(mpdz_lm_500, 
         # include Urb_250 as the indepedent variable from a multi-regression model
         terms = "Urb_500", 
         # include partial residuals
         partial.resid = T, se = T, 
         # x and y labels, 
         xlab = "Proportion of Urban Area (500m)",
         ylab = "ses.MPD", 
         # use filled-in grey points
         pch = 16)
abline(h = 0, lty = 2)

# sesMFD  vs Urb_250
termplot(mfdz_lm_250, 
         # include Urb_250 as the indepedent variable from a multi-regression model
         terms = "Urb_250", 
         # include partial residuals and standard error curve
         partial.resid = T, se = T,
         # x and y labels, 
         xlab = "Proportion of Urban Area (250m)",
         ylab = "ses.MFD", 
         # use filled-in grey points
         pch = 16)
abline(h = 0, lty = 2)

# sesMFD vs Urb_500 
plot(mfd.obs.z ~ Urb_500, data = final, 
     xlab = 'Proportion of Urban Area (500m)', 
     ylab = "ses.MFD", 
     col = "grey", lwd = 3, lty = 2, pch = 16)
abline(h = 0, lty = 2)

# MFD vs Urb_250
termplot(mfd_lm_250, 
         # include Urb_250 as the indepedent variable from a multi-regression model
         terms = "Urb_250", 
         # include partial residuals
         partial.resid = T, se = T,
         # x and y labels, 
         xlab = "Proportion of Urban Area (250m)",
         ylab = "MFD", 
         # use filled-in grey points
         pch = 16)

# MFD vs Urb_500 
plot(mfd.obs ~ Urb_500, data = final, 
     xlab = 'Proportion of Urban Area (500m)', 
     ylab = "MFD", 
     col = "grey", lwd = 3, lty = 2, pch = 16)

# MPD vs Urb_250
termplot(mpd_lm_250, 
         # include Urb_250 as the indepedent variable from a multi-regression model
         terms = "Urb_250", 
         # include partial residuals and standard error curve
         partial.resid = T, se = T,
         # x and y labels, 
         xlab = "Proportion of Urban Area (250m)",
         ylab = "MPD", 
         # use filled-in grey points
         pch = 16)

# MPD vs Urb_500
termplot(mpd_lm_500, 
         # include Urb_250 as the indepedent variable from a multi-regression model
         terms = "Urb_500", 
         # include partial residuals and standard error curve
         partial.resid = T, se = T,
         # x and y labels, 
         xlab = "Proportion of Urban Area (500m)",
         ylab = "MPD", 
         # use filled-in grey points
         pch = 16)

# phylogenetic tree 
plot.phylo(tree, label.offset = 0.05)

# rlq biplots: community composition along RLQ axes
ggplot(rlq.bees$lQ, aes(x=AxcQ1, y=AxcQ2)) +
       # insert vertical line at (0, 0)
       geom_vline(xintercept = 0, lwd = 0.5, col = "grey") + 
       # insert horizonal line at (0, 0)
       geom_hline(yintercept = 0, lwd = 0.5, col = "grey") +
       # use points to indicate species composition on biplots
       geom_point(shape = 1, pch = 0.5) +
       # create non-overlapping labels 
       geom_text_repel(aes(x=AxcQ1, y=AxcQ2), label = rownames(rlq.bees$lQ)) + 
       # increase space for cover type symbols
       ylim(-10, 10) + 
       # x-axis and y-axis labels
       xlab("RLQ Axis 1") + 
       ylab("RLQ Axis 2") +
       # remove gridlines and grey background
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       panel.background = element_blank(), axis.line = element_line(colour = "black")) 

#rlq bilots: trait distribution along RLQ axes 
ggplot(rlq.bees$c1, aes(x = CS1, y = CS2)) +
       # insert vertical line at (0, 0)
       geom_vline(xintercept = 0, lwd = 0.5, col = "grey") + 
       # insert horizonal line at (0, 0)
       geom_hline(yintercept = 0, lwd = 0.5, col = "grey") +
       # use points to indicate position of traits on biplots
       geom_point(shape = 1 ,pch = 0.5) + 
       # create non-overlapping labels 
       geom_text_repel(aes(x = CS1, y = CS2), 
                       label = rownames(rlq.bees$c1),
                       box.padding = unit(0.45, "lines")) + 
       # increase space for cover type symbols
       ylim(-4, 4) + 
       # remove gridlines and grey background 
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       panel.background = element_blank(), axis.line = element_line(colour = "black")) 

# rlq biplots: environmental variables along RLQ axes
s.arrow(rlq.bees$l1, boxes = F, grid = F) # remove boxes and gridlines
