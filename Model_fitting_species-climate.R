###########################################
# Supporting Information
# "Ecological traits linking climate with bird abundances and distributions"
# Duarte S. Viana, Jonathan M. Chase
###########################################

# Fit BRT and GAM models

# Passing arguments to R from command lines
args = commandArgs(trailingOnly=TRUE)
output.file <- args[1]

# try to get SGE_TASK_ID from submit script, otherwise fall back to 1
task.id = as.integer(Sys.getenv("SGE_TASK_ID", "1"))

library(tidyr)
library(gbm)
library(mgcv)
source("R_Functions.R")
# load data info for the full dataset
#load("Data_id.Rdata")
# load data info for the subsampled datasets
load("Data_id.Rdata")
for(i in 1:ncol(data.id)) data.id[,i] <- as.numeric(data.id[,i])
# load BBS data
load("BBS_BIOCLIM.Rdata")

# Prepare data

# subset data
wsites <- idata[[paste(cases$year[task.id])]]
dat <- bbs.clim[bbs.clim$year==cases$year[task.id] & bbs.clim$site_id %in% wsites,]
dat.wide <- spread(dat, key=species_id, value=abundance)
names(itrain) <- 1983:2018
in.train <- itrain[[paste(cases$year[task.id])]]
names(itest) <- 1983:2018
in.test <- itest[[paste(cases$year[task.id])]]
sp.data <- data.id[data.id$year==cases$year[task.id],]
in.sp <- sp.data$aou[sp.data$Npres.train>=5 & sp.data$Npres.test>=5]
site.id <- dat.wide[dat.wide$site_id %in% wsites,"site_id"]

# Set data matrices
mat.sp0 <- dat.wide[dat.wide$site_id %in% wsites,24:ncol(dat.wide)]
mat.sp0[is.na(mat.sp0)] <- 0
mat.sp <- mat.sp0[,names(mat.sp0) %in% in.sp] # only include species with at least 5 presences in train data and 5 in test data
if(cases$data.type[task.id]=="pres") mat.sp[mat.sp>0] <- 1
if(cases$method[task.id]=="BRT_Gau") {
  for(j in 1:ncol(mat.sp)){
    spj <- mat.sp[,j]
    mat.sp[,j] <- sqrt(spj)/max(sqrt(spj))
  } 
}
# Climatic variables
# Group 1: isothermality (BIO3), mean temperature of the warmest quarter (BIO10), 
# precipitation seasonality (BIO15), precipitation of the wettest quarter (BIO16) 
# Group 2: max temperature of warmest month (BIO5), min temperature of coldest month (BIO6),
# precipitation of wettest month (BIO13), precipitation of driest month (BIO14)
mat.env <- dat.wide[dat.wide$site_id %in% wsites,c("bio2","bio3","bio5","bio8","bio9","bio15","bio16","bio18")]


# Prepare output object
sp.ids <- unique(bbs.clim$species_id)
res.vp <- array(NA, c(1, length(sp.ids), 4),
                dimnames=list(case=task.id,sp=paste(sp.ids),
                              stat=c("D2","R2","D2.CV","R2.CV")))


# Modelling
try({
  if(cases$method[task.id]=="BRT_lr.01"){
    # BRT
    if(cases$data.type[task.id]=="abund"){
      r2.env <- R2.BRT(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      distr="poisson", inter.depth=2, lr=0.01)
    } 
    if(cases$data.type[task.id]=="pres"){
      r2.env <- R2.BRT(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      distr="bernoulli", inter.depth=2, lr=0.01)
    } 
  }
  
  if(cases$method[task.id]=="BRT_lr.001"){
    # BRT
    if(cases$data.type[task.id]=="abund"){
      r2.env <- R2.BRT(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      distr="poisson", inter.depth=2, lr=0.001)
    } 
    if(cases$data.type[task.id]=="pres"){
      r2.env <- R2.BRT(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      distr="bernoulli", inter.depth=2, lr=0.001)
    } 
  }
  
  if(cases$method[task.id]=="GAM_kest"){
    # GAM with optimised splines 
    if(cases$data.type[task.id]=="abund"){
      r2.env <- R2.GAM(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      family="poisson", env.eff = "splines")
    }
    if(cases$data.type[task.id]=="pres"){
      r2.env <- R2.GAM(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      family="binomial", env.eff = "splines")
    }
  }
  
  if(cases$method[task.id]=="GAM_k3"){
    # GAM with fixed number of knots
    if(cases$data.type[task.id]=="abund"){
      r2.env <- R2.GAM(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      family="poisson", env.eff = "splines", k=3)
    }
    if(cases$data.type[task.id]=="pres"){
      r2.env <- R2.GAM(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      family="binomial", env.eff = "splines", k=3)
    }
  }
  
  if(cases$method[task.id]=="BRT_Gau"){
    r2.env <- R2.BRT(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                     distr="gaussian", inter.depth=2, lr=0.01)
  } 
  
  # Store results
  res.vp[1,names(mat.sp),1] <- r2.env[[1]]
  res.vp[1,names(mat.sp),2] <- r2.env[[2]]
  res.vp[1,names(mat.sp),3] <- r2.env[[3]]
  res.vp[1,names(mat.sp),4] <- r2.env[[4]]
})

# Save output to a .Rdata file
save(res.vp,file=output.file)



