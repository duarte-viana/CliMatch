###########################################
# Supporting Information
# "Ecological traits linking climate with bird abundances and distributions"
# Duarte S. Viana, Jonathan M. Chase
###########################################

# Prepare data 
# Take the climate matching estimates (species-climate r2) and add species traits
library(dplyr)
library(tidyr)

# Load data
setwd("~/Documents/Papers/iDiv/Paper_BBS_R2-traits")
load("results.Rdata")
load("Data_id.Rdata")
for(i in 1:ncol(data.id)) data.id[,i] <- as.numeric(data.id[,i])

id.sp <- unique(data.id$aou[data.id$Npres.train>=5 & data.id$Npres.test>=5])
id <- data.id[data.id$aou %in% id.sp,]

# species-climate r2
r2s <- as.data.frame(res[,,4])

# raw data
lres <- id[rep(1:nrow(id),4),]
lres$method <- rep(c("BRT_lr.01","BRT_lr.001","GAM_kest","GAM_k3"),each=nrow(id))
lres$aou <- as.character(lres$aou)
lres$r2 <- NA
data.type <- "abund"
#data.type <- "pres"
for(j in unique(cases$method)){
  for(i in unique(lres$aou)){
    lres$r2[lres$method==j & lres$aou==i] <- r2s[which(cases$data.type==data.type & cases$method==j),names(r2s)==i]
  }
}
lres <- na.exclude(lres)

load("species_names.Rdata")
sps <- sps[!duplicated(sps$aou),]
sps$aou <- as.character(sps$aou)
lresi <- left_join(lres,sps[,c("aou","sci","sporder","family")],by="aou")

sp.names <- unique(lresi$sci)
birdtree.sp <- read.csv("BLIOCPhyloMasterTax.csv")
all(sp.names %in% birdtree.sp$TipLabel)

# Correct unmatched species names
mis.sp <- data.frame(sci=sp.names[!(sp.names %in% birdtree.sp$TipLabel)])
mis.sp <- left_join(mis.sp, sps[,c("TipLabel","sci")], by="sci")
mis.sp$TipLabel <- as.character(mis.sp$TipLabel)
mis.sp$TipLabel[mis.sp$sci=="Circus_hudsonius"] <- "Circus_hudsonius"
#mis.sp$TipLabel[mis.sp$sci=="Ammospiza_leconteii"] <- "Ammodramus_leconteii"
mis.sp$TipLabel[mis.sp$sci=="Setophaga_coronata"] <- "Dendroica_coronata"
mis.sp$TipLabel[mis.sp$sci=="Setophaga_nigrescens"] <- "Dendroica_nigrescens"
mis.sp$TipLabel[mis.sp$sci=="Artemisiospiza_nevadensis"] <- "Amphispiza_belli"
mis.sp$TipLabel[mis.sp$sci=="Aphelocoma_woodhouseii"] <- "Aphelocoma_californica"
lresi$phylo <- lresi$sci
for(i in mis.sp$sci) lresi$phylo[lresi$sci==i] <- mis.sp$TipLabel[mis.sp$sci==i]
lresi <- na.exclude(lresi)
lresi <- left_join(lresi,birdtree.sp[,c("TipLabel","BLFamilyLatin","BLFamilyEnglish")],by=c("phylo"="TipLabel"))

# Add species traits
traits <- read.table("traits.txt",header=T,sep="\t")
names(traits) <- c("species","lifespan","log.mass","body.length","wingspan",names(traits[,6:ncol(traits)]))
lresi1 <- left_join(lresi,traits,by=c("phylo"="species"))
lresi2 <- left_join(lresi,traits,by=c("sci"="species"))
lresi1[is.na(lresi1$log.mass),] <- lresi2[is.na(lresi1$log.mass),]
lresi <- as.data.frame(lresi1)

# Add habitat specialization index
load("specialization_index.Rdata")
lresi1 <- left_join(lresi,ssi,by=c("phylo"="sci"))
lresi2 <- left_join(lresi,ssi,by="sci")
lresi1[is.na(lresi1$SSI),] <- lresi2[is.na(lresi1$SSI),]
lresi <- as.data.frame(lresi1)
lresi$SSI <- as.numeric(lresi$SSI)

# Add competition index
load("competition_index.Rdata")
lresi <- left_join(lresi,ci,by="aou")

# Add threat status
iucn <- read.csv("IUCN_status.csv")
iucn$sci <- paste(iucn$genusName,iucn$speciesName,sep="_")
lresi <- left_join(lresi,iucn[,c("sci","redlistCategory")],by="sci")
lresi$redlistCategory <- as.character(lresi$redlistCategory)
lresi$redlistCategory[lresi$phylo=="Picoides_villosus"] <- "Least Concern"
lresi$redlistCategory[lresi$phylo=="Icterus_bullockii"] <- "Least Concern"
lresi$redlistCategory[lresi$phylo=="Vermivora_ruficapilla"] <- "Least Concern"
lresi$redlistCategory[lresi$phylo=="Coccothraustes_vespertinus"] <- "Vulnerable"
lresi$redlistCategory[lresi$phylo=="Vermivora_celata"] <- "Least Concern"
lresi$redlistCategory[lresi$phylo=="Aphelocoma_californica"] <- "Least Concern"
lresi$redlistCategory[lresi$phylo=="Dryocopus_pileatus"] <- "Least Concern"

# save(lresi,file="lresi.Rdata")
# save(lresi,file="lresi.pres.Rdata")


#----------------------------------------------------------

library(dplyr)
library(tidyr)
setwd("~/Documents/Papers/iDiv/Paper_BBS_R2-traits")
load("lresi.Rdata")

lres2 <- lresi[lresi$method=="BRT_lr.01",]

# Descriptive statistics of climate matching
mean(lres2$r2)
sd(lres2$r2)
quantile(lres2$r2,c(0.025,0.975))
range(lres2$r2)
length(unique(lres2$aou))


#######################################################################################
#######################################################################################

# Phylogenetic signal

# Get species names (for the filtered data) and remove underscore
sp.names <- gsub('_', ' ', unique(lresi$phylo))
lres2 <- lresi[lresi$method=="BRT_lr.01",]
lres2 <- lres2[,c("r2","year","phylo","abund.mean","abund.cv","log.mass","migration","main.habitat","SSI","c.max")]
lres2 <- na.exclude(lres2)
sp.names2 <- gsub('_', ' ', unique(lres2$phylo))

#Export the vector as a one-column .csv file 
#write.csv(sp.names, "sp.names.csv", row.names = FALSE)
#write.csv(sp.names2, "sp.names2.csv", row.names = FALSE)

# Once exported, copy the species names and paste into the `Select species` form 
# on the [Phylogeny Subsets](http://birdtree.org/subsets/) page at Birdtree.org. 
# We downloaded 100 trees from the *Ericsson All Species* dataset for our analyses. 
# Once processed, save the resulting .nex file and read in the multiphylo object 
# using functions in the `ape` package
setwd("/Users/Viana/Documents/Papers/iDiv/Paper_BBS_R2-traits/Phylo_100_BirdTree")
sp.trees <- ape::read.nexus("output.nex")
setwd("/Users/Viana/Documents/Papers/iDiv/Paper_BBS_R2-traits/Phylo_100_BirdTree2")
sp.trees2 <- ape::read.nexus("output.nex")
setwd("/Users/Viana/Documents/Papers/iDiv/Paper_BBS_R2-traits")
# Now check to make sure that all of the species names are represented in the tree. 
# Here, we have to include the underscore once again, as this is included in Birdtree.org phylogeny subsets. 
# This call should return `TRUE` if there are no unmatched names
sp.names.underscore <- gsub(' ', '_', sp.names)
all(sp.names.underscore %in% sp.trees[[1]]$tip.label) # Circus hudsonius is not in the tree!

# Phylogenetic covariance
cons.tree <- sp.trees[[1]] # 1 phylogenetic tree
# From https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html
# The phylo object contains information on the relationship between species. 
# Using this information, we can construct a covariance matrix of species (Hadfield & Nakagawa, 2010).

# Build a consensus tree
library(phytools)
cons.tree <- consensus.edges(sp.trees,method="mean.edge")
cons.tree <- di2multi(cons.tree)
A <- ape::vcv.phylo(cons.tree)
cons.tree2 <- consensus.edges(sp.trees2,method="mean.edge")
cons.tree2 <- di2multi(cons.tree2)
A2 <- ape::vcv.phylo(cons.tree2)

# Covariance matrix for a nested design
# A block-diagonal matrix with phylogenetic effect within years
# Following an example in phyr package: https://rdrr.io/cran/phyr/man/pglmm.html
# Response values are assumed to be phylogenetically uncorrelated between years (i.e. covariance=0)
# Use a subsample of years, otherwise the model gets too heavy
library(Matrix)
y.in <- seq(1983,2018,5)
ny <- length(y.in)
lA <- list()
for(i in 1:ny) lA[[i]] <- A
At <- as.matrix(bdiag(lA))
spl <- paste(rep(colnames(A),ny))
colnames(At) <- paste(spl,rep(y.in,each=ncol(A)),sep="_")
rownames(At) <- paste(spl,rep(y.in,each=ncol(A)),sep="_")
A <- At
diag(A) <- diag(A)+0.00000002 # to ensure positive definitiveness
lA <- list()
for(i in 1:ny) lA[[i]] <- A2
At <- as.matrix(bdiag(lA))
spl <- paste(rep(colnames(A2),ny))
colnames(At) <- paste(spl,rep(y.in,each=ncol(A2)),sep="_")
rownames(At) <- paste(spl,rep(y.in,each=ncol(A2)),sep="_")
A2 <- At
diag(A2) <- diag(A2)+0.00000002 # to ensure positive definitiveness
lresi <- lresi[lresi$abund.cv<7,]
# save(lresi,A,A2,file="data_for_brms.Rdata")


# Plot phylogeny with trait
# http://www.phytools.org/Cordoba2017/ex/15/Plotting-methods.html

# Mean of R2 across years
lres2 <- lresi[lresi$method=="BRT_lr.01",]
lres2 <- lres2[lres2$phylo!="Circus_hudsonius",]
lres.y <- lres2 %>%
  group_by(phylo) %>%
  summarize(r2 = mean(r2))
lres.y <- as.data.frame(lres.y)
trait.y<-setNames(lres.y$r2,lres.y$phylo)
tax <- lres2[!duplicated(lres2$phylo),c("phylo","BLFamilyLatin","family","sporder")]

# Plot
library(RColorBrewer)
pal <-  colorRampPalette(c("blue","yellow","red"))
obj<-contMap(cons.tree,trait.y,plot=FALSE)
obj<-setMap(obj,colors=pal(20)) # colors=rev(brewer.pal(11,"RdYlBu"))
quartz(height=4,width=5)
plot(obj,fsize=c(0.5,0.7),outline=FALSE,lwd=c(2,7),leg.txt="",legend=30,
     type="fan",ftype="off",xlim=c(-140,145),ylim=c(-120,120))
for(i in unique(tax$family)){
  spi <- tax$phylo[tax$family==i]
  if(length(spi)>1) arc.cladelabels(text=i,node=findMRCA(cons.tree,spi),orientation="horizontal",
                                    ln.offset=1.02,lab.offset=1.07,mark.node=FALSE,cex=0.7)
  if(length(spi)==1) arc.cladelabels(text=i,node=which(cons.tree$tip.label==spi),orientation="horizontal",
                                     ln.offset=1.02,lab.offset=1.07,mark.node=FALSE,cex=0.7)
}

# Phylogentic autocorrelation
library(phylosignal)
library(phylobase)
corg <- phyloCorrelogram(phylo4d(x = cons.tree,tip.data = trait.y), n.points = 20)
plot.phylocorrelogram(corg)

# Phylogenetic signal
library(pez)
phylosig(cons.tree, trait.y, method = "lambda", test = TRUE)

#######################################################################################
#######################################################################################

# Model climate matching as a function of species traits

library(brms)
library(performance)
library(car)

load("data_for_brms.Rdata")

#----------------------------------------------------------

# Models with phylogenetic autocorrelation

if(task.id==1){
  # Model traits with phylogenetic effect
  lres2 <- lresi[lresi$method=="BRT_lr.01",]
  lres2 <- lres2[,c("r2","year","phylo","abund.mean","abund.cv","log.mass","migration","main.habitat","SSI","c.max")]
  lres2 <- na.exclude(lres2)
  lres2 <- lres2[lresi$phylo!="Circus_hudsonius",]
  lres2 <- lres2[lres2$year %in% seq(1983,2018,5),]
  lres2$year <- as.factor(lres2$year)
  lres2$log.mass <- scale(lres2$log.mass,scale=FALSE) 
  lres2$abund.mean <- scale(log(lres2$abund.mean), scale=FALSE)
  lres2$abund.cv <- scale(log(lres2$abund.cv+1), scale=FALSE)
  lres2$main.habitat <- factor(lres2$main.habitat)
  lres2$SSI <- scale(as.numeric(lres2$SSI,scale=FALSE))
  lres2$c.max <- scale(as.numeric(lres2$c.max,scale=FALSE))
  lres2$sp.year <- paste(lres2$phylo,lres2$year,sep="_")
  A <- A2
  
  mm <- brm(r2 ~ poly(log.mass,2) + migration + main.habitat + SSI + c.max + abund.mean + poly(abund.cv,2) +
            (poly(log.mass,2) + migration + main.habitat + SSI + c.max + abund.mean + poly(abund.cv,2) ||year) +
            (1|gr(sp.year, cov = A)),
            data=lres2, data2 = list(A = A), family=Beta, chains=4, thin=1, iter=3000, warmup=1000, 
            control = list(adapt_delta = .99, max_treedepth = 15), cores=4)
}

if(task.id==2){
  # Model IUCN threat status
  lres2 <- lresi[lresi$method=="BRT_lr.01",]
  lres2 <- lres2[,c("r2","year","phylo","abund.mean","abund.cv","redlistCategory")]
  lres2 <- na.exclude(lres2)
  lres2 <- lres2[lresi$phylo!="Circus_hudsonius",]
  lres2 <- lres2[lres2$year %in% seq(1983,2018,5),]
  lres2$year <- as.factor(lres2$year)
  lres2$abund.mean <- scale(log(lres2$abund.mean))
  lres2$abund.cv <- scale(log(lres2$abund.cv+1))
  #lres2$redlistCategory[lres2$redlistCategory %in% c("Near Threatened","Vulnerable")] <- "Threatened"
  lres2$redlistCategory <- factor(lres2$redlistCategory)
  lres2$sp.year <- paste(lres2$phylo,lres2$year,sep="_")
  
  mm <- brm(r2 ~ redlistCategory + abund.mean + poly(abund.cv,2) +
            (redlistCategory + abund.mean + poly(abund.cv,2) ||year) +
            (1|gr(sp.year, cov = A)),
            data=lres2, data2 = list(A = A), family=Beta, chains=4, thin=1, iter=3000, warmup=1000, 
            control = list(adapt_delta = .99, max_treedepth = 15), cores=4)
}

if(task.id==3){
  # Null model with phylogenetic effect
  lres2 <- lresi[lresi$method=="BRT_lr.01",]
  lres2 <- lres2[,c("r2","year","phylo")]
  lres2 <- lres2[lresi$phylo!="Circus_hudsonius",]
  lres2 <- na.exclude(lres2)
  lres2 <- lres2[lres2$year %in% seq(1983,2018,5),]
  lres2$year <- as.factor(lres2$year)
  lres2$sp.year <- paste(lres2$phylo,lres2$year,sep="_")
  
  mm <- brm(r2 ~ 1 +
            (1|gr(sp.year, cov = A)),
            data=lres2, data2 = list(A = A), family=Beta, chains=4, thin=1, iter=3000, warmup=1000, 
            control = list(adapt_delta = .99, max_treedepth = 15), cores=4)
}



#----------------------------------------------------------

# Single predictor models and nested models (without phylogenetic autocorrelation)

#load("lresi.Rdata")
load("lresi.pres.Rdata")
lres2 <- lresi[lresi$method=="BRT_lr.01",]
lres2 <- lres2[lres2$abund.cv<7,]
lres2$year <- factor(lres2$year)
lres2$migration <- factor(lres2$migration)
lres2$main.habitat <- factor(lres2$main.habitat)
lres2$abund.mean <- log(lres2$abund.mean)
lres2$abund.cv <- log(lres2$abund.cv+1)

xvars <- c("poly(log.mass,2)","migration","main.habitat","SSI","c.max","abund.mean","poly(abund.cv,2)")
mod.id <- data.frame(xvar=c("all",rep(xvars,each=2)))
mod.id$type <- c("full",rep(c("uni","drop"),7))
mod.id$xvar <- as.character(mod.id$xvar)


if(mod.id$xvar[task.id]=="all"){
  lres2 <- lres2[,c("r2","year","abund.mean","abund.cv","log.mass","migration","main.habitat","SSI","c.max")]
  lres2 <- na.exclude(lres2)
  lres2$log.mass <- scale(lres2$log.mass,scale=FALSE) 
  lres2$abund.mean <- scale(lres2$abund.mean, scale=FALSE)
  lres2$abund.cv <- scale(lres2$abund.cv, scale=FALSE)
  lres2$main.habitat <- factor(lres2$main.habitat)
  lres2$SSI <- scale(lres2$SSI,scale=FALSE)
  lres2$c.max <- scale(lres2$c.max,scale=FALSE)
  mm <- brm(r2 ~ poly(log.mass,2) + migration + main.habitat + SSI + c.max + abund.mean + poly(abund.cv,2) +
                     (poly(log.mass,2) + migration + main.habitat + SSI + c.max + abund.mean + poly(abund.cv,2) ||year),
                   data=lres2, family=Beta, chains=4, thin=10, iter=3000, warmup=1000, 
                   control = list(adapt_delta = .99), cores=4)
}

if(mod.id$xvar[task.id]!="all"){
  if(mod.id$type[task.id]=="uni"){
    if(mod.id$xvar[task.id]=="poly(log.mass,2)") lres2 <- lres2[,c("r2","year","log.mass")]
    if(mod.id$xvar[task.id]=="poly(abund.cv,2)") lres2 <- lres2[,c("r2","year","abund.cv")]
    if(!(mod.id$xvar[task.id] %in% c("poly(log.mass,2)","poly(abund.cv,2)"))) lres2 <- lres2[,c("r2","year",mod.id$xvar[task.id])]
    lres2 <- na.exclude(lres2)
    if(is.numeric(lres2[,3])) lres2[,3] <- scale(lres2[,3],scale=FALSE)
    if(is.factor(lres2[,3])) lres2[,3] <- factor(lres2[,3])
    rand.part <- paste("+ (", mod.id$xvar[task.id], " ||year)", sep="")
    formula.mod <- as.formula(paste("r2 ~ ", mod.id$xvar[task.id], rand.part, sep=""))
    mm <- brm(formula.mod, data=lres2, family=Beta, chains=4, thin=10, iter=3000, warmup=1000, 
                  control = list(adapt_delta = .99), cores=4)
  }
  if(mod.id$type[task.id]=="drop"){
    lres2 <- lres2[,c("r2","year","abund.mean","abund.cv","log.mass","migration","main.habitat","SSI","c.max")]
    lres2 <- na.exclude(lres2)
    lres2$log.mass <- scale(lres2$log.mass,scale=FALSE) 
    lres2$abund.mean <- scale(lres2$abund.mean, scale=FALSE)
    lres2$abund.cv <- scale(lres2$abund.cv, scale=FALSE)
    lres2$main.habitat <- factor(lres2$main.habitat)
    lres2$SSI <- scale(lres2$SSI,scale=FALSE)
    lres2$c.max <- scale(lres2$c.max,scale=FALSE)
    mod.terms <- paste(xvars[xvars!=mod.id$xvar[task.id]], collapse = "+")
    rand.part <- paste("+ (", mod.terms, " ||year)", sep="")
    formula.mod <- as.formula(paste("r2 ~ ", mod.terms, rand.part, sep=""))
    mm <- brm(formula.mod, data=lres2, family=Beta, chains=4, thin=10, iter=3000, warmup=1000, 
              control = list(adapt_delta = .99), cores=4)
  }
}

#----------------------------------------------------------

# Model IUCN threat status
lres2 <- lresi[lresi$method=="BRT_lr.01",]
lres2 <- lres2[,c("r2","year","phylo","abund.mean","abund.cv","redlistCategory")]
lres2 <- lres2[lres2$abund.cv<7,]
lres2 <- na.exclude(lres2)
lres2$year <- as.factor(lres2$year)
lres2$abund.mean <- scale(log(lres2$abund.mean), scale=FALSE)
lres2$abund.cv <- scale(log(lres2$abund.cv+1), scale=FALSE)
#lres2$redlistCategory[lres2$redlistCategory %in% c("Near Threatened","Vulnerable")] <- "Threatened"
lres2$redlistCategory <- factor(lres2$redlistCategory)

vif(lm(r2 ~ redlistCategory + abund.mean + abund.cv, data=lres2))

m.iucn <- brm(r2 ~ redlistCategory + abund.mean + abund.cv +
                (redlistCategory + abund.mean + abund.cv ||year),
              data=lres2, family=Beta, chains=4, thin=10, iter=3000, warmup=1000, 
              control = list(adapt_delta = .99), cores=4)






