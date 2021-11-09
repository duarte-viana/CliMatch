###########################################
# Supporting Information
# "Ecological traits linking climate with bird abundances and distributions"
# Duarte S. Viana, Jonathan M. Chase
###########################################

# Data preparation

library(tidyr)
library(sp)
library(rgdal)
library(raster)
library(rgeos)
library(geosphere)
library("mapdata")

# Load BBS data
#bbs.clim <- read.table("BBS_BIOCLIM.txt",sep="\t",header=T)
#save(bbs.clim,file="BBS_BIOCLIM.Rdata")
load("BBS_BIOCLIM.Rdata")
bbs.wide <- spread(bbs.clim, key=species_id, value=abundance)

# All sites surveyed across all years
sites <- bbs.clim[!duplicated(bbs.clim$site_id),c("long","lat","site_id")]
row.names(sites) <- sites$site_id
sites.sp<-SpatialPointsDataFrame(sites[,1:2],data=data.frame(sites$site_id),
                                 proj4string=CRS("+proj=longlat +datum=WGS84"))

#----------------------------------------------------

# Visualisation and exploration

# Map sites
map('worldHires')
points(sites.sp, col='red')
(e<-bbox(extent(sites.sp)*1.2))

quartz(height=4,width=6)
par(mar=c(0,0,0,0))
map('worldHires', xlim = e[1, ], ylim = e[2, ])
points(sites.sp, col='darkgreen', pch=16, cex=0.3)

# Visualise cross-validation procedure
sites.in.id <- unique(bbs.clim$site_id[bbs.clim$year==2018 & !is.na(bbs.clim$bio1)])
sites.in <- sites.sp[sites.sp$sites.site_id %in% sites.in.id,]
lonlat <- coordinates(sites.in)
centers <-  regularCoordinates(16)
dists <-  pointDistance(centers, lonlat, longlat = TRUE)
radius <-  1.8e5
in.test <- sites.in[apply(dists, 2, function(x) min(x) <= radius),]
quartz(height=4,width=6)
par(mar=c(0,0,0,0))
e<-bbox(extent(sites.in)*1.2)
map('worldHires', xlim = e[1, ], ylim = e[2, ])
#box()
points(sites.in, col="dodgerblue3", pch=16, cex=0.5)
points(in.test, col="brown3", pch=16, cex=0.5)
#legend("bottomleft",legend=c("train data","test data"),pch=c(16,16),col=c("dodgerblue3","brown3"),bty="n")



################################################################################
################################################################################

# Data selection for subsequent analyses

# Temporal partitioning: all available sites for each year
years <- 1983:2018 # exclude 1982 because it has a smaller spatial extent
bbs.wide <- spread(bbs.clim, key=species_id, value=abundance)

data.id <- data.frame(year=numeric(0), aou=numeric(0), Npres=numeric(0), Npres.train=numeric(0), 
                      Npres.test=numeric(0), abund.mean=numeric(0), abund.cv=numeric(0), abund.total=numeric(0))

# All the points within inner.radius of a center point will be in the test set.
# Everything more than outer.radius away will be in the training set.
# radii are in meters
# Adapted from Harris 2015
l.centers <- list()
jk <- 0
for(j in 1:5){
  centers <-  regularCoordinates(j+14)
  for(k in 1:10){
    jk <- jk+1
    l.centers[[jk]] <- jitter(centers, 70)
  }
}
radius <-  1.8e5

itrain <- list()
itest <- list()
id.sp <- 0

for(y in years){
  # All sites surveyed in year y
  bbs.y <- bbs.wide[bbs.wide$year==y,]
  bbs.y <- bbs.y[!duplicated(bbs.y$site_id),]
  bbs.y <- bbs.y[!is.na(bbs.y$bio1),]
  sites.in <- unique(bbs.y$site_id)
  
  lonlat <- bbs.y[,c("long","lat")]
  # split training and test data (for 50-fold CV)
  in.train <- list()
  in.test <- list()
  for(j in 1:50){
    centers <-  l.centers[[j]]
    dists <-  pointDistance(centers, lonlat, longlat = TRUE)
    in.train[[j]] <- sites.in[apply(dists, 2,function(x) min(x) > radius)]
    in.test[[j]] <- sites.in[apply(dists, 2, function(x) min(x) <= radius)]
  }
  itrain[[y-1982]] <- in.train
  itest[[y-1982]] <- in.test
  
  # species information
  for(i in 1:424){
    id.sp <- id.sp+1
    bbs.sp.i <- cbind(bbs.y[,c("long","lat","site_id")], abund=bbs.y[,23+i])
    bbs.sp.i$abund[is.na(bbs.sp.i$abund)] <- 0
    Npres.train.j <- c()
    Npres.test.j <- c()
    for(j in 1:50){
      Npres.train.j[j] <- nrow(bbs.sp.i[bbs.sp.i$site_id %in% in.train[[j]] & bbs.sp.i$abund>0,])
      Npres.test.j[j] <- nrow(bbs.sp.i[bbs.sp.i$site_id %in% in.test[[j]] & bbs.sp.i$abund>0,])
    }
    Npres.train <- min(Npres.train.j)
    Npres.test <- min(Npres.test.j)
    
    bbs.sp.i <- bbs.sp.i[bbs.sp.i$abund>0,]
    Npres <- nrow(bbs.sp.i)
    if(Npres<3) data.id[id.sp,] <- c(y, names(bbs.y)[23+i], Npres, 
                                     Npres.train, Npres.test, NA,NA,NA)
    if(Npres>=3){
      # abundance metrics
      abund.mean <- mean(bbs.sp.i$abund)
      abund.cv <- sd(bbs.sp.i$abund)/mean(bbs.sp.i$abund)
      abund.total <- sum(bbs.sp.i$abund)
      
      data.id[id.sp,] <- c(y, names(bbs.y)[23+i], Npres, Npres.train, Npres.test, 
                           abund.mean, abund.cv, abund.total)
    }
  }
}

# Sites for each year
idata <- list()
for(y in years){
  # All sites surveyed in year y
  bbs.y <- bbs.wide[bbs.wide$year==y,]
  bbs.y <- bbs.y[!duplicated(bbs.y$site_id),]
  bbs.y <- bbs.y[!is.na(bbs.y$bio1),]
  idata[[y-1982]] <- unique(bbs.y$site_id)
}
names(idata) <- 1983:2018

# cases for climatic models
cases <- data.frame(year=1983:2018)
# Replicate dataset to apply other methods (BRT and GAM)
ni <- nrow(cases)
cases <- cases[rep(1:ni,2),,drop=F]
cases$data.type <- rep(c("abund","pres"),each=ni)
ni <- nrow(cases)
cases <- cases[rep(1:ni,5),]
cases$method <- rep(c("BRT_lr.01","BRT_lr.001","GAM_kest","GAM_k3","BRT_Gau"),each=ni)
cases <- cases[!(cases$method=="BRT_Gau" & cases$data.type=="pres"),]

save(idata,itrain,itest,data.id,cases,file="Data_id.Rdata")



################################################################################
################################################################################

# Get species names and AOU codes

library(DBI)
#library(RSQLite)
library(dplyr)
#library(dbdplyr)

setwd("~/Documents/Papers/iDiv/Paper_sDiv/Data/BBS")

bbs_db <- dbConnect(RSQLite::SQLite(), 'bbs_sqlite.db')
sps <- tbl(bbs_db, "breed_bird_survey_species") %>% data.frame()

sci.names <- paste(sps$genus,sps$species,sep="_")
sps$sci <- sci.names

# Add BirdTree names
birdtree.sp <- read.csv("BLIOCPhyloMasterTax.csv")
sps <- left_join(sps, birdtree.sp[,c("TipLabel","English")], by=c("english_common_name" = "English"))

save(sps,file="species_names.Rdata")


################################################################################
################################################################################

# Extract specialization index (Martin and Fahrig 2018)

library(tabulizer)
options(java.parameters = "-Xmx8000m")

tab <- extract_tables("https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.2428&file=ecy2428-sup-0002-AppendixS2.pdf",
                      output = "data.frame",pages=1:10)
tab2 <- lapply(tab,function(x) x[-1,1:4])
for(i in 1:length(tab2)){
  tabi <- tab2[[i]]
  names(tabi) <- c("sci","common.name","N","SSI")
  tab2[[i]] <- tabi
}
ssi <- do.call("rbind",tab2)
ssi$sci[ssi$sci=="Campylorhynchus"] <- "Campylorhynchus brunneicapillus"
ssi <- ssi[ssi$sci!="brunneicapillus",]
ssi$common.name[ssi$common.name=="black-throated green"] <- "black-throated green warbler"
ssi <- ssi[ssi$common.name!="warbler",]
ssi$common.name[ssi$common.name=="northern rough-winged"] <- "northern rough-winged swallow"
ssi <- ssi[ssi$common.name!="swallow",]
ssi$sci[ssi$sci=="Xanthocephalus"] <- "Xanthocephalus xanthocephalus"
ssi <- ssi[ssi$sci!="xanthocephalus",]
ssi$sci <- gsub(' ', '_', ssi$sci)
setwd("~/Documents/Papers/iDiv/Paper_BBS_R2-traits")
save(ssi,file="specialization_index.Rdata")



