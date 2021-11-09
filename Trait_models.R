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

# Traits dataset
traits <- lresi[!duplicated(lresi$aou),c("sci","phylo")]

# Add habitat specialization index (from Martin and Fahrig 2018)
load("specialization_index.Rdata")
traits1.1 <- left_join(traits,ssi[,c("sci","SSI")],by=c("phylo"="sci"))
traits1.2 <- left_join(traits,ssi[,c("sci","SSI")],by="sci")
traits1.1[is.na(traits1.1$SSI),] <- traits1.2[is.na(traits1.1$SSI),]
traits <- as.data.frame(traits1.1)
traits$SSI <- as.numeric(traits$SSI)
a <- which(is.na(traits$SSI))
unique(traits$sci[a])

# Add diet and body mass
# from Wilman et al.; EltonTraits 1.0; https://doi.org/10.6084/m9.figshare.c.3306933.v1
traits2 <- read.csv2("Traits_Wilman.csv")
traits2$Scientific <- gsub(' ', '_', traits2$Scientific)
traits2$Scientific[traits2$Scientific=="Circus_cyaneus"] <- "Circus_hudsonius"
traits2 <- traits2[,c(8,36,10:19)]
#traits2$diet.breadth <- 1/apply((traits2[,-c(1,2)]/100)^2,1,sum)
traits <- left_join(traits,traits2,by=c("phylo"="Scientific"))
# Species with missing traits
a <- which(is.na(traits$BodyMass.Value))
unique(traits$sci[a])
traits$log.mass <- log(as.numeric(traits$BodyMass.Value))

# Calculate a continuous variable for diet 
# from Cooke et al. (https://doi.org/10.1038/s41467-019-10284-z)
library(FD)
diet_all <- traits %>% 
  # select diet data
  dplyr::select(sci, contains("Diet")) %>% 
  tibble::remove_rownames() %>% 
  # add species names to rownames (needed for gowdis function)
  tibble::column_to_rownames(var = "sci") %>% 
  as.data.frame()
diet_all <- diet_all[,-ncol(diet_all)]
# calculate species x species gower distance matrix based on traits
diet_gd <- FD::gowdis(diet_all)
# perform principal coordinates analysis (PCoA)
diet_pco <- ade4::dudi.pco(diet_gd, scannf = FALSE)
pc_diet <- diet_pco$tab
summary(diet_pco)
# Projected inertia Ax1 = 27.733
# Projected inertia Ax2 = 17.059
plot(cbind(pc_diet[,1],diet_all))
plot(cbind(pc_diet[,2],diet_all))
# Add diet
diet <- data.frame(sci=row.names(pc_diet),diet.pc1=pc_diet[,1])
traits <- left_join(traits,diet,by="sci")

# Pairwise correlation plots
library(psych)
diet.cors <- data.frame(cbind(PC1=pc_diet[,1],diet_all))
quartz(height=8,width=8)
par(mar=c(1,1,1,1))
pairs.panels(diet.cors, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = FALSE # show correlation ellipses
)

# Add more traits
# from Sheard et al. 2020; https://doi.org/10.1038/s41467-020-16313-6 
traits3 <- read.table("Traits_Sheard_HWI.txt",header=T,sep="\t")
traits3$phylo <- gsub(' ', '_', traits3$"Tree.name")
traits3$phylo[traits3$phylo=="Circus_cyaneus"] <- "Circus_hudsonius"
traits <- left_join(traits,traits3[,c("phylo","HWI","Territoriality","Habitat","Migration.1")],by="phylo")
a <- which(is.na(traits$HWI))
unique(traits$phylo[a])
traits$Migration.1[traits$sci=="Pica_hudsonia"] <- "Not a Migrant"
traits$migration <- NA
traits$migration[traits$Migration.1=="Full Migrant"] <- "migrant"
traits$migration[traits$Migration.1!="Full Migrant"] <- "non-migrant"
traits$main.habitat <- factor(traits$Habitat)
levels(traits$main.habitat) <- c("dense","semi-open","open")

# Add threat status
iucn <- read.csv("IUCN_status.csv")
iucn$sci <- paste(iucn$genusName,iucn$speciesName,sep="_")
traits <- left_join(traits,iucn[,c("sci","redlistCategory")],by="sci")
traits$redlistCategory <- as.character(traits$redlistCategory)
traits$redlistCategory[traits$phylo=="Picoides_villosus"] <- "Least Concern"
traits$redlistCategory[traits$phylo=="Icterus_bullockii"] <- "Least Concern"
traits$redlistCategory[traits$phylo=="Vermivora_ruficapilla"] <- "Least Concern"
traits$redlistCategory[traits$phylo=="Coccothraustes_vespertinus"] <- "Vulnerable"
traits$redlistCategory[traits$phylo=="Vermivora_celata"] <- "Least Concern"
traits$redlistCategory[traits$phylo=="Aphelocoma_californica"] <- "Least Concern"
traits$redlistCategory[traits$phylo=="Dryocopus_pileatus"] <- "Least Concern"

# Traits table
write.table(traits,row.names=FALSE,sep="\t",file="Trait_data.txt")

# All data
lresi <- left_join(lresi,traits[,c("sci","log.mass","migration","HWI","main.habitat","diet.pc1","SSI","Territoriality","redlistCategory")],by="sci")
lresi$phylo[lresi$phylo=="Circus_hudsonius"] <- "Circus_cyaneus"
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

# From https://github.com/nicholasjclark/BBS.occurrences/blob/master/Clark_etal_analysis/Appendix_S4_PhyloTraitData.Rmd
# Get species names (for the filtered data) and remove underscore
sp.names <- gsub('_', ' ', unique(lresi$phylo))
lres2 <- lresi[lresi$method=="BRT_lr.01",]
lres2 <- lres2[,c("r2","year","phylo","abund.mean","abund.cv","log.mass","migration","HWI","main.habitat","diet.pc1","SSI","Territoriality")]
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
lres.y <- lres2 %>%
  group_by(phylo) %>%
  summarize(r2 = mean(r2))
lres.y <- as.data.frame(lres.y)
trait.y<-setNames(lres.y$r2,lres.y$phylo)
tax <- lres2[!duplicated(lres2$phylo),c("phylo","BLFamilyLatin","family","sporder")]

# Plot
library(RColorBrewer)
pal <-  colorRampPalette(c("blue","yellow","red"))
cons.tree <- consensus.edges(sp.trees,method="mean.edge")
cons.tree <- multi2di(cons.tree)
obj<-contMap(cons.tree,trait.y,plot=FALSE)
obj<-setMap(obj,colors=pal(20)) # colors=rev(brewer.pal(11,"RdYlBu"))
quartz(height=4,width=5.5)
plot(obj,fsize=c(0.5,0.7),outline=FALSE,lwd=c(2,7),leg.txt="",legend=30,
     type="fan",ftype="off",xlim=c(-140,145),ylim=c(-120,120))
for(i in unique(tax$family[!tax$family %in% c("Ptiliogonatidae","Bombycillidae","Sturnidae","Regulidae","Aegithalidae","Alaudidae")])){
  spi <- tax$phylo[tax$family==i]
  if(length(spi)>1) arc.cladelabels(text=i,node=findMRCA(cons.tree,spi),orientation="horizontal",
                                    ln.offset=1.02,lab.offset=1.07,mark.node=FALSE,cex=0.7)
  if(length(spi)==1) arc.cladelabels(text=i,node=which(cons.tree$tip.label==spi),orientation="horizontal",
                                     ln.offset=1.02,lab.offset=1.07,mark.node=FALSE,cex=0.7)
}
for(i in unique(tax$family[tax$family %in% c("Ptiliogonatidae","Bombycillidae","Sturnidae","Regulidae","Aegithalidae","Alaudidae")])){
  spi <- tax$phylo[tax$family==i]
  if(length(spi)>1) arc.cladelabels(text=NULL,node=findMRCA(cons.tree,spi),orientation="horizontal",
                                    ln.offset=1.02,lab.offset=1.07,mark.node=FALSE,cex=0.7)
  if(length(spi)==1) arc.cladelabels(text=NULL,node=which(cons.tree$tip.label==spi),orientation="horizontal",
                                     ln.offset=1.02,lab.offset=1.07,mark.node=FALSE,cex=0.7)
}


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
  lres2 <- lres2[,c("r2","phylo","year","abund.mean","abund.cv","log.mass","migration","HWI","main.habitat","diet.pc1","SSI","Territoriality")]
  lres2 <- na.exclude(lres2)
  lres2 <- lres2[lres2$year %in% seq(1983,2018,5),]
  lres2$year <- as.factor(lres2$year)
  lres2$log.mass <- scale(lres2$log.mass,scale=FALSE) 
  lres2$abund.mean <- scale(log(lres2$abund.mean), scale=FALSE)
  lres2$abund.cv <- scale(log(lres2$abund.cv+1), scale=FALSE)
  lres2$main.habitat <- factor(lres2$main.habitat)
  lres2$Territoriality <- factor(lres2$Territoriality)
  lres2$SSI <- scale(as.numeric(lres2$SSI,scale=FALSE))
  lres2$HWI <- scale(as.numeric(lres2$HWI,scale=FALSE))
  lres2$sp.year <- paste(lres2$phylo,lres2$year,sep="_")
  A <- A2
  
  mm <- brm(r2 ~ poly(log.mass,2) + migration + poly(HWI,2) + main.habitat + diet.pc1 + SSI + Territoriality + abund.mean + poly(abund.cv,2) +
              (poly(log.mass,2) + migration + poly(HWI,2) + main.habitat + diet.pc1 + SSI + Territoriality + abund.mean + poly(abund.cv,2) ||year) +
              (1|gr(sp.year, cov = A)),
            data=lres2, data2 = list(A = A), family=Beta, chains=4, thin=1, iter=3000, warmup=1000, 
            control = list(adapt_delta = .99, max_treedepth = 15), cores=4)
}

if(task.id==2){
  # Model traits with phylogenetic effect, but without SSI (to include all species)
  lres2 <- lresi[lresi$method=="BRT_lr.01",]
  lres2 <- lres2[,c("r2","phylo","year","abund.mean","abund.cv","log.mass","migration","HWI","main.habitat","diet.pc1","Territoriality")]
  lres2 <- na.exclude(lres2)
  lres2 <- lres2[lres2$year %in% seq(1983,2018,5),]
  lres2$year <- as.factor(lres2$year)
  lres2$log.mass <- scale(lres2$log.mass,scale=FALSE) 
  lres2$abund.mean <- scale(log(lres2$abund.mean), scale=FALSE)
  lres2$abund.cv <- scale(log(lres2$abund.cv+1), scale=FALSE)
  lres2$main.habitat <- factor(lres2$main.habitat)
  lres2$Territoriality <- factor(lres2$Territoriality)
  lres2$HWI <- scale(as.numeric(lres2$HWI,scale=FALSE))
  lres2$sp.year <- paste(lres2$phylo,lres2$year,sep="_")
  
  mm <- brm(r2 ~ poly(log.mass,2) + migration + poly(HWI,2) + main.habitat + diet.pc1 + Territoriality + abund.mean + poly(abund.cv,2) +
              (poly(log.mass,2) + migration + poly(HWI,2) + main.habitat + diet.pc1 + Territoriality + abund.mean + poly(abund.cv,2) ||year) +
              (1|gr(sp.year, cov = A)),
            data=lres2, data2 = list(A = A), family=Beta, chains=4, thin=1, iter=3000, warmup=1000, 
            control = list(adapt_delta = .99, max_treedepth = 15), cores=4)
}

if(task.id==3){
  # Model IUCN threat status
  lres2 <- lresi[lresi$method=="BRT_lr.01",]
  lres2 <- lres2[,c("r2","year","phylo","abund.mean","abund.cv","redlistCategory")]
  lres2 <- na.exclude(lres2)
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

if(task.id==4){
  # Null model with phylogenetic effect
  lres2 <- lresi[lresi$method=="BRT_lr.01",]
  lres2 <- lres2[,c("r2","year","phylo")]
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

# data variations
xvars <- c("poly(log.mass,2)","migration","poly(HWI,2)","main.habitat","diet.pc1","SSI","Territoriality","abund.mean","poly(abund.cv,2)")
mod.id0 <- data.frame(xvar=c("all",rep(xvars,each=2)))
mod.id0$type <- c("full",rep(c("uni","drop"),9))
mod.id0$xvar <- as.character(mod.id0$xvar)
data.id <- c("all.counts","all.pres","nonmig.counts","mig.counts")
mod.id <- mod.id0[rep(1:nrow(mod.id0),length(data.id)),]
mod.id$data.id <- rep(data.id,each=nrow(mod.id0))
mod.id <- mod.id[!(mod.id$data.id %in% c("nonmig.counts","mig.counts") & mod.id$xvar=="migration"),]

# run models
for(task.id in 69:nrow(mod.id)){
  # load data
  if(mod.id$data.id[task.id]=="all.pres") load("lresi.pres.Rdata")
  if(mod.id$data.id[task.id]!="all.pres"){
    load("lresi.Rdata")
    if(mod.id$data.id[task.id]=="nonmig.counts") lresi <- lresi[lresi$migration=="non-migrant",]
    if(mod.id$data.id[task.id]=="mig.counts") lresi <- lresi[lresi$migration=="migrant",]
  } 
  
  # Prepare data
  lres2 <- lresi[lresi$method=="BRT_lr.01",]
  lres2 <- lres2[lres2$abund.cv<7,]
  lres2$year <- factor(lres2$year)
  lres2$migration <- factor(lres2$migration)
  lres2$main.habitat <- factor(lres2$main.habitat)
  lres2$Territoriality <- factor(lres2$Territoriality)
  lres2$abund.mean <- log(lres2$abund.mean)
  lres2$abund.cv <- log(lres2$abund.cv+1)
  
  # Run models
  
  if(mod.id$xvar[task.id]=="all"){
    lres2 <- lres2[,c("r2","year","abund.mean","abund.cv","log.mass","migration","HWI","main.habitat","diet.pc1","SSI","Territoriality")]
    lres2 <- na.exclude(lres2)
    lres2$log.mass <- scale(lres2$log.mass,scale=FALSE) 
    lres2$abund.mean <- scale(lres2$abund.mean, scale=FALSE)
    lres2$abund.cv <- scale(lres2$abund.cv, scale=FALSE)
    lres2$main.habitat <- factor(lres2$main.habitat)
    lres2$Territoriality <- factor(lres2$Territoriality)
    lres2$SSI <- scale(lres2$SSI,scale=FALSE)
    lres2$HWI <- scale(lres2$HWI,scale=FALSE)
    if(!(mod.id$data.id[task.id] %in% c("nonmig.counts","mig.counts"))){
      mm <- brm(r2 ~ poly(log.mass,2) + migration + poly(HWI,2) + main.habitat + diet.pc1 + SSI + Territoriality + abund.mean + poly(abund.cv,2) +
                  (poly(log.mass,2) + migration + poly(HWI,2) + main.habitat + diet.pc1 + SSI + Territoriality + abund.mean + poly(abund.cv,2) ||year),
                data=lres2, family=Beta, chains=4, thin=10, iter=3000, warmup=1000, 
                control = list(adapt_delta = .99), cores=4)
    } 
    if(mod.id$data.id[task.id] %in% c("nonmig.counts","mig.counts")){
      mm <- brm(r2 ~ poly(log.mass,2) + poly(HWI,2) + main.habitat + diet.pc1 + SSI + Territoriality + abund.mean + poly(abund.cv,2) +
                  (poly(log.mass,2) + poly(HWI,2) + main.habitat + diet.pc1 + SSI + Territoriality + abund.mean + poly(abund.cv,2) ||year),
                data=lres2, family=Beta, chains=4, thin=10, iter=3000, warmup=1000, 
                control = list(adapt_delta = .99), cores=4)
    } 
  }
  
  if(mod.id$xvar[task.id]!="all"){
    if(mod.id$type[task.id]=="uni"){
      if(mod.id$xvar[task.id]=="poly(log.mass,2)") lres2 <- lres2[,c("r2","year","log.mass")]
      if(mod.id$xvar[task.id]=="poly(abund.cv,2)") lres2 <- lres2[,c("r2","year","abund.cv")]
      if(mod.id$xvar[task.id]=="poly(HWI,2)") lres2 <- lres2[,c("r2","year","HWI")]
      if(!(mod.id$xvar[task.id] %in% c("poly(log.mass,2)","poly(abund.cv,2)","poly(HWI,2)"))) lres2 <- lres2[,c("r2","year",mod.id$xvar[task.id])]
      lres2 <- na.exclude(lres2)
      if(is.numeric(lres2[,3])) lres2[,3] <- scale(lres2[,3],scale=FALSE)
      if(is.factor(lres2[,3])) lres2[,3] <- factor(lres2[,3])
      rand.part <- paste("+ (", mod.id$xvar[task.id], " ||year)", sep="")
      formula.mod <- as.formula(paste("r2 ~ ", mod.id$xvar[task.id], rand.part, sep=""))
      mm <- brm(formula.mod, data=lres2, family=Beta, chains=4, thin=10, iter=3000, warmup=1000, 
                control = list(adapt_delta = .99), cores=4)
    }
    if(mod.id$type[task.id]=="drop"){
      lres2 <- lres2[,c("r2","year","abund.mean","abund.cv","log.mass","migration","HWI","main.habitat","diet.pc1","SSI","Territoriality")]
      lres2 <- na.exclude(lres2)
      lres2$log.mass <- scale(lres2$log.mass,scale=FALSE) 
      lres2$abund.mean <- scale(lres2$abund.mean, scale=FALSE)
      lres2$abund.cv <- scale(lres2$abund.cv, scale=FALSE)
      lres2$main.habitat <- factor(lres2$main.habitat)
      lres2$Territoriality <- factor(lres2$Territoriality)
      lres2$SSI <- scale(lres2$SSI,scale=FALSE)
      lres2$HWI <- scale(lres2$HWI,scale=FALSE)
      if(!(mod.id$data.id[task.id] %in% c("nonmig.counts","mig.counts"))){
        mod.terms <- paste(xvars[xvars!=mod.id$xvar[task.id]], collapse = "+")
        rand.part <- paste("+ (", mod.terms, " ||year)", sep="")
        formula.mod <- as.formula(paste("r2 ~ ", mod.terms, rand.part, sep=""))
        mm <- brm(formula.mod, data=lres2, family=Beta, chains=4, thin=10, iter=3000, warmup=1000, 
                  control = list(adapt_delta = .99), cores=4)
      } 
      if(mod.id$data.id[task.id] %in% c("nonmig.counts","mig.counts")){
        mod.terms <- paste(xvars[xvars!=mod.id$xvar[task.id] & xvars!="migration"], collapse = "+")
        rand.part <- paste("+ (", mod.terms, " ||year)", sep="")
        formula.mod <- as.formula(paste("r2 ~ ", mod.terms, rand.part, sep=""))
        mm <- brm(formula.mod, data=lres2, family=Beta, chains=4, thin=10, iter=3000, warmup=1000, 
                  control = list(adapt_delta = .99), cores=4)
      } 
    }
  }
  save(mm,file=paste("brms_all_",task.id,".Rdata",sep=""))
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
lres2$redlistCategory <- factor(lres2$redlistCategory)

m.iucn <- brm(r2 ~ redlistCategory + abund.mean + abund.cv +
                (redlistCategory + abund.mean + abund.cv ||year),
              data=lres2, family=Beta, chains=4, thin=10, iter=3000, warmup=1000, 
              control = list(adapt_delta = .99), cores=4)





