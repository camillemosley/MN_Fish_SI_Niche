### Stable Isotope Data merge 
# CLM 9-2-2024

#clear environment and plots 
rm(list=ls())
graphics.off()

#libraries 
library(tidyverse)
library(SIBER)

##### SIBER analysis of ellipses/trophic niche for lake communities ####

#working directory in Gdrive 
setwd("/Users/cammosley/Google Drive/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/") #iMac
setwd("/Users/cammosley/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM")#macbook

#bring in data
BigLakes<-read.csv("/Users/cammosley/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/Large lakes 2017 2018/2017 and 2018 Organism and Isotope Data_Canon as of 08012019_flagged for baselines 10292019.csv")
glimpse(BigLakes)
BigLakes<-BigLakes[,c(3,4,7,13,17,18,25,30,31,33,34)] # "Year","Lake","Group","Site1","Taxa_general","Taxa_specific1","YOY","Percent.N","N15.corrected","Percent.C","C13.corrected" 

Lakes<-read.csv("~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/Cross lakes 2019 2021/Copy of WAE_FDWBS_duplicate_correction_record_NB.xlsx - Sheet1.csv")
glimpse(Lakes)
Lakes<-Lakes[,c(4,5,8,13,17,18,25,29,30,32,33)]#"Year","Lake" ,"Group","Site1","Taxa_general","Taxa_specific1","YOY","Percent.N","N15.corrected","Percent.C","C13.corrected" 
LakesAll<-rbind(BigLakes,Lakes)#join lake data 

FWLakes<-read.csv("~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/zm_foodwebs_iso_data_22Jul24.csv")

#getting columns of interest - group, lake, year, taxa, %C and %Nvalues 
#looking at fish samples
#BigLakes<-BigLakes[BigLakes$Group=="Fish",]

LakesAll<-LakesAll[!is.na(LakesAll$Percent.C),]
LakesAll<-LakesAll[!is.na(LakesAll$Percent.N),]

a<-ggplot(LakesAll, aes(Percent.C,Percent.N,))+
  geom_point(aes(colour=Taxa_general))+
  xlab("σ C (%)")+
  ylab("σ N (%)")+
  theme_bw()+
  facet_wrap(~Lake)

addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
# Apply on original plot
addSmallLegend(a)  

#format data for SIBER 
# iso 1 =x, iso 2 = y
#iso 1 = %C , iso 2 = %N

LakesAll$iso1<-LakesAll$Percent.C
LakesAll$iso2<-LakesAll$Percent.N

#Making index for community (lake)
LakesAll<- LakesAll %>%
     mutate(index = as.integer(factor(Lake, levels = unique(Lake))))
LakesAll$community<-LakesAll$index #1-11 note lakes are labeled 1-11 alphabetically 

##Making index for group (organism type)
LakesAll<- LakesAll %>%
  mutate(index = as.integer(factor(Group, levels = unique(Group))))
LakesAll$group<-LakesAll$index #1=Inverts, 2=Zoops, 3=Fish, 4=Bird
#LakesAll<-LakesAll[LakesAll$group==c(1,2,3),]#removing bird groups from the dataset 
#LakesAll<-LakesAll[LakesAll$community==c(1,3,5:6),]#subsetting for half the lakes 

Test<-LakesAll[,c(12,13,16,15)]
#make dataset a SIBER object for analysis

Test<-createSiberObject(Test)

#Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")



par(mfrow=c(1,1))
plotSiberObject(Test,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille')
)

par(mfrow=c(1,1))

community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

# this time we will make the points a bit smaller by 
# cex = 0.5
plotSiberObject(Test,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'permille'),
                ylab=expression({delta}^15*N~'permille'),
                cex = 0.5
)

group.ML <- groupMetricsML(Test)
print(group.ML)

#plot the group ellipses 
plotGroupEllipses(Test, n = 100, p.interval = 0.95,
                  lty = 1, lwd = 2)

community.ML <- communityMetricsML(Test) 
print(community.ML)

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(Test, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior) 

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('permille' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group",
                 ylims = c(0,60)
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

# extract the posterior means
mu.post <- extractPosteriorMeans(Test, ellipses.posterior)

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)


# Visualise the first community
siberDensityPlot(layman.B[[1]], xticklabels = colnames(layman.B[[1]]), 
                 bty="L", ylim = c(0,10))

# add the ML estimates (if you want). Extract the correct means 
# from the appropriate array held within the overall array of means.
comm1.layman.ml <- laymanMetrics(Test$ML.mu[[1]][1,1,],
                                 Test$ML.mu[[1]][1,2,]
)
points(1:6, comm1.layman.ml$metrics, col = "red", pch = "x", lwd = 2)


# Visualise the second community
siberDensityPlot(layman.B[[2]], xticklabels = colnames(layman.B[[2]]), 
                 bty="L", ylim = c(0,10))

# add the ML estimates. (if you want) Extract the correct means 
# from the appropriate array held within the overall array of means.
comm2.layman.ml <- laymanMetrics(Test$ML.mu[[2]][1,1,],
                                 Test$ML.mu[[2]][1,2,]
)
points(1:6, comm2.layman.ml$metrics, col = "red", pch = "x", lwd = 2)

# Alternatively, pull out TA from both and aggregate them into a 
# single matrix using cbind() and plot them together on one graph.

# go back to a 1x1 panel plot
par(mfrow=c(1,1))

siberDensityPlot(cbind(layman.B[[1]][,"TA"], layman.B[[2]][,"TA"]),
                 xticklabels = c("Community 1", "Community 2"), 
                 bty="L", ylim = c(0,10),
                 las = 1,
                 ylab = "TA - Convex Hull Area",
                 xlab = "")

##### Trophic position ####
library(tRophicPosition)
library(stringr)
library(dplyr)

####### Baseline Corrections ######

#C and N ratios of fish relative to their baselines within in each lake 
#looking at data in package to compare to my data formating for TP package 
data("Bilagay")
Bilagay <- Bilagay %>% mutate(Community = paste(Study,"-", Location, sep = ""))

#bring in data *BIG LAKES HAS BASELINES wait until the end* 
Lakes<-read.csv("~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/Cross lakes 2019 2021/Copy of WAE_FDWBS_duplicate_correction_record_NB.xlsx - Sheet1.csv")
FWLakes<-read.csv("~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/zm_foodwebs_iso_data_22Jul24.csv")

#check baselines flag when present 
#lakes should have both baselines 

FWLakes<-FWLakes %>%
  mutate(Community=paste(lake,"-",year,sep="")) %>%
  rename(Spp=taxa_general,d13C=d.13.C,d15N=d.15.N,Group=organism_type,Taxa_specific1=taxa_specific_1,Length=length, Baseline.Flag=baseline_flag)

FWLakes<-FWLakes[,c(4,8,9,13,14,15:17,18,25)] #"Group","Spp" ,Taxa_specific1","Length","YOY","FFG","Baseline.Flag","d13C" ,"d15N","Community"   

Lakes<-Lakes %>%
  mutate(Community=paste(Lake,"-",Year,sep="")) %>%
  rename(Spp=Taxa_general,d13C=d.13.C,d15N=d.15.N)

Lakes<-Lakes[,c(8,17,18,23,25,22,43,46,56,61)]#"Group","Spp","Taxa_specific1","Length","YOY","FFG","d13C","d15N","Baseline.Flag","Community"  
Lakes<-Lakes[c(1:38),]#remove empty rows from gsheets

#Big Sand 2020 lakes have all in spp colunmn for zoop samples?

#join datasets together with columns needed for TP function 
All<-rbind(Lakes,FWLakes) 

#label baseline samples 
All$FG<-ifelse(All$Spp %in% c("Zooplankton","Deepwter Zooplankton","Copepods & cladocerans",
                    "Bivalvia","Dreissena","Unionidae"), "Pelagic_BL","Littoral_BL" )
#All$FG<-ifelse(All$FFG %in% c("CG","SC","SH"), "Littoral_BL" ,"Pelagic_BL")
#label consumer samples
Frows<-which(All$Group=="Fish") #fish rows 
All[c(Frows),11]="Consumer"

#fix species labels for fish 
Fish<-All[c(Frows),]

#make subset with baseline data to rejoin after fixing species labels
BL<-which(All$FG=="Littoral_BL")

BP<-which(All$FG=="Pelagic_BL" )
#remove omnivores 
BLs<-All[c(BL,BP),]
BLsFix<-BLs[BLs$Community%in%c("Buffalo-2021","Big Sandy-2021","Alexander-2021"),]
#making a spp label for zooplankton samples for pelagic baseline sort 
BLsFix$Spp[BLsFix$Group=="Zooplankton"]="Zooplankton"
BLsFix$Taxa_specific1[BLsFix$Group=="Zooplankton"]="Zooplankton"
BLsFix$FG[BLsFix$Group=="Zooplankton"]="Pelagic_BL"
BLsFix<-BLsFix[!BLsFix$Taxa_specific1=="unionidae",]

BLs<-BLs[!is.na(BLs$Baseline.Flag),]

#Knife 2019 needs gastropod samples
BLKnife<-BLs[BLs$Community=="Knife-2019",]
BLs<-BLs[!BLs$Spp=="Gastropoda",]#removing gastropods 
BLs<-BLs[!BLs$Taxa_specific1%in%c("sphaeriidae","ephemeridae"),]
BLs<-BLs[!BLs$Community=="Alexander-2021",]
BLs<-rbind(BLs,BLsFix)
BLs<-BLs[!is.na(BLs$Community),]

BLsFix<-BLs[BLs$Community=="Bemidji-2022",]
BLs<-BLs[!BLs$Community=="Bemidji-2022",]
BLsFix<-BLsFix[!BLsFix$Spp=="Bivalvia",]
BLs<-rbind(BLs,BLsFix)
BLs<-rbind(BLs,BLKnife)#adding knife samples back in

#get aplhabetical sort of names 
sort(unique(Fish$Taxa_specific1)) #35values
sort(unique(Fish$Spp)) #23 values

#matching species names with codes from isotope_data_entry_2018_08062018.xlsx (gdrive)
Fish$Spp<-case_when(
  Fish$Spp == "Amploplites" ~ "RKB",
  Fish$Spp == "Catostomus" ~ "WTS",
  Fish$Spp == "Culaea" ~ "BST",
  Fish$Spp == "Cyprinella" ~ "SFS",
  Fish$Spp == "Esox" ~ "NOP",
  Fish$Spp == "Fundulus" ~ "BKF",
  Fish$Spp == "Labidesthes" ~ "BKS",
  Fish$Spp == "Lota" ~ "BUB",
  Fish$Spp == "Luxilus" ~ "CSH",
  Fish$Spp == "Moxostoma" ~ "SHR",
  Fish$Spp == "Notemigonus" ~ "GOS",
  Fish$Spp == "Perca" ~ "YEP",
  Fish$Spp == "Percina" ~ "LGP",
  Fish$Spp == "Percopsis" ~ "TRP",
  Fish$Spp == "Pomoxis" ~ "BLC",
  Fish$Spp == "Sander" ~ 'WAE',
  Fish$Taxa_specific1 == "artedi" ~ "CIS",
  Fish$Taxa_specific1== "clupeaformis" ~ "LKW",
  Fish$Taxa_specific1 == "notatus" ~ "BNM",
  Fish$Taxa_specific1 == "nigrum" ~ "JND",
  Fish$Taxa_specific1 == "flabellare" ~ "FTD",
  Fish$Taxa_specific1 == "macrochirus" ~ "BLG",
  Fish$Taxa_specific1 == "gibbosus" ~ "PMK",
  Fish$Taxa_specific1 == "cyanellus" ~ "GSF",
  Fish$Taxa_specific1 == "salmoides" ~ "LMB",
  Fish$Taxa_specific1 == "dolomieu" ~ "SMB",
  Fish$Taxa_specific1 == "volucellus" ~ "MMS",
  Fish$Taxa_specific1 == "heterodon" ~ "BCS",
  Fish$Taxa_specific1 == "hudsonius" ~ "SPO",
  Fish$Taxa_specific1 == "promelas" ~ "FHM",
  Fish$Taxa_specific1 == "heterolepis" ~ "BNS",
  Fish$Taxa_specific1 == "exile" ~ "IOD",
  Fish$Taxa_specific1 == "rupestris" ~ "RKB",
  Fish$Taxa_specific1 == "dorsalis" ~ "BMS",
  Fish$Taxa_specific1 == "Hybrid" ~ "HSF",
  Fish$Taxa_specific1 == "microlophus" ~ "RSF") #redear sunfish

check<-Fish[is.na(Fish$Spp),]#0 observations that need fixing 
unique(check$Taxa_specific1)


RKB<-"Ambloplites"
WTS<-"Catostomus"
CIS<-"Coregonus"
BST<-"Culaea"
SFS<-"Cyprinella"
NOP <-"Esox"
# JND<-"Etheostoma nigrum"
# FTD<-"Etheostoma" flabellare
BKF<-"fundulus" 
BKS<-"Labidesthes"
# BLG<-"Lepomis" macrochirus
# PMK<-"Lepomis"gibbosus
# GSF<-"Lepomis" cyanellus
BUB<-"Lota"
CSH<-"Luxilus"
# LMB<-micropterus salmoides
# SMB<-Micropterus dolomieu
SHR<-"Moxostoma"
GOS<-"Notemigonus"
# MMS<-Notropis volucellus
# BCS<-Notropis heterodon
# SPO<-Notropis hudsonius
# BNS <-Notropis heterolepsis 
YEP<-"Perca"
LGP<-"Percina"
TRP<-"Percopsis"
# BNM<-Pimephales notatus
# FHM<-Pimephales promelas
BLC<-"Pomoxis"
WAE<-"Sander"
#microlophus - redear sunfish -> RSF


#Pelagic-Bivalves, zoops, anything eats plankton 

#cg collector gatherers, SH shredder, SC scrapper 

#littoral grazers and scrapers

#primary consumers filter feeders and zoops 

#predators not good for baselines

#get YOY labels
for (i in 1:nrow(Fish)){
  if (Fish$YOY[i]%in%c("Y")){
    Fish$Spp[i]=paste(Fish$Spp[i],"Y",sep = "_")
  }
  else{
    Fish$Spp[i]=Fish$Spp[i]
  }
}
#join data for spp corrected dataset
All<-rbind(BLs,Fish)
#remove Knife lake from samples
All<-All[!All$Community=="Knife-2019",]

#subset for life stages 
Fish_adult<-Fish[Fish$YOY=="N",]
All_adult<-rbind(BLs,Fish_adult)# 3080 observations 
Fish_yoy<-Fish[Fish$YOY=="Y",]
All_YOY<-rbind(BLs,Fish_yoy)#1901 observations 

#arrange by community ie lake year samples + remove na values/empty rows
All<- All %>% arrange(Community)

All<-All[,c(1:5,7,8,10,11)] #"Group","Spp",Taxa_specific1","Length","YOY",FFG","d13C","d15N","Community","FG"  

#making edits on baselines given data exploration plots

#make Isotope data object in list format for TP package 
#changed the Consumner column from SPP to FG to get plots for each lake -year
AllList<-extractIsotopeData(All,b1="Pelagic_BL", b2="Littoral_BL", baselineColumn = "FG", consumersColumn = "FG",
                            groupsColumn = "Community", d13C = "d13C",d15N = "d15N")
str(AllList)

save(All,file = "MNLakesIsotopesJoined_010925.csv")#saving output of cleaned SI samples

#get summary stats and plots for each community 
for (community in AllList) {
  print(summary(community))
  plot(community)
}


#each community has to be input seperately into this function for it to run
#for example Alexander Lake
screenIsotopeData(isotopeData = AllList$`Alexander-2021-Consumer`,
                  density = "both",
                  consumer = "Consumer",
                  b1 = "Pelagic_BL",
                  b2 = "Littoral_BL",
                  legend = c(1.15, 1.15),
                  title = "Alexander Lake Stable Isotope Samples",
                  xylim = NULL,)
screenIsotopeData(isotopeData = AllList$`Rush-2022-Consumer`,
                  density = "both",
                  consumer = "Consumer",
                  b1 = "Pelagic_BL",
                  b2 = "Littoral_BL",
                  legend = c(1.15, 1.15),
                  title = "Rush 2022",
                  xylim = NULL,)
#make subset of data to check values
Sub<-All[,c(2,7,8,10)]
Sub$Species<-Sub$Spp
Sub<-Sub[,2:5]
screenFoodWeb(Sub)
extractIsotopeData()

###### WALLEYE #####
#trophic discrimantion factors
TDF_values <- TDF()

Fish<-loadIsotopeData(All,consumer = "WAE",consumersColumn = "Spp",b1="Pelagic_BL",b2="Littoral_BL", baselineColumn = "FG",
                      groupsColumn = "Community", deltaC = TDF_values$deltaC,
                      deltaN = TDF_values$deltaN)
#Default multi species TP uses lamda =2 (BLs), n. chains = number of MCMC simulations, n.iter etc. are iterations
Lake_models <- multiSpeciesTP(Fish, model = "twoBaselinesFull",
                                 n.adapt = 10000, n.iter = 10000,
                                 burnin = 10000, n.chains = 5, print = FALSE)

#The function multiSpeciesTP() returns 4 objects: 1. a list named multiSpeciesTP, which includes the raw data returned by posteriorTP() for each consumer/species per group/community/sampling location; 
#2. a data frame named df with the mode, median and credibility confidence interval for trophic position and alpha (if a two baselines model was chosen), grouped by model, consumer/species and group/community/sampling location; 
#3. a list named TPs with the posterior samples of trophic position for each group/community/location and consumer/species; and 
#4. a list named Alphas which includes the posterior samples of alpha (relative contribution of baseline 1) for each species by group/community/location.

# By default the mode is used in both trophic position and alpha plots
credibilityIntervals(Lake_models$df, x = "group", xlab ="Community")

# If you want to use the median instead of the mode,
# just add y1 and y2 as arguments
credibilityIntervals(Lake_models$df, x = "group", xlab ="Community", 
                     y1 = "median", y2 = "alpha.median")

# To get a numerical summary
sapply(Lake_models$"TPs", quantile, probs = c(0.025, 0.5, 0.975)) %>% round(3)

# To get the mode
getPosteriorMode(Lake_models$"TPs")
#note these outputs group all consumers, need to make edit to loop through each community specific iso data object

##### TP for inividial species in lakes ####

IsoData_Sp<-extractIsotopeData(All,b1="Pelagic_BL", b2="Littoral_BL", baselineColumn = "FG", consumersColumn = "FG",
                            groupsColumn = "Community", d13C = "d13C",d15N = "d15N")
str(IsoData_Sp)

#remove NAs
IsoData <-lapply(IsoData_Sp,na.omit)

#using Alexander lake for the momement will need to create for loop to look at consumers over all samples

  #trying to look at mutliple models to calcualte tp for a species across mutliple locations , see TP appendix 3

#multiple models for alexander lake
#By default multiModelsTP() defines a lambda = 2 for the baselines, uses 2 chains
#(n.chains = 2) to do the Bayesian calculation with 20,000 adaptive iterations (n.adapt = 20000), 20,000
#actual iterations (n.iter = 20000), 20,000 iterations as burnin (burnin = 20000) and a thinning of 10 (thin = 10).
                                                                                                       
LM<-multiModelTP(AllList$`Alexander-2021-Consumer`, model = "twoBaselinesFull",
             n.adapt = 10000, n.iter = 10000,
             burnin = 10000, n.chains = 5, print = FALSE)

str(LM)
credibilityIntervals(LM$gg,x="model") 
#code works but makes a interval for all consumers together, grouping not correct for inference

#try subsetting for 2 species of interest
Fish_BW<-Fish[Fish$Spp%in%c("WAE"),]
BW<-rbind(BLs,Fish_BW)
BWList<-extractIsotopeData(BW,b1="Pelagic_BL", b2="Littoral_BL", baselineColumn = "FG", consumersColumn = "FG",
                            groupsColumn = "Community", d13C = "d13C",d15N = "d15N")
str(BWList)

#remove NAs
ID <-lapply(BWList,na.omit)

#get summary stats and plots for each community 
for (community in ID) {
  print(summary(community))
  plot(community)
}

####### DATA CHECKING BY LAKE RAW SAMPS #########
#lakes that need fixing 
#washington 2020, Tenmile 2020, Steamboat 2021, Rush 2022, Round 2022, potato 2022
#Plantagenet-2022, Melissa 2019, Leech 2022, koronis 2019, knife 2019
#Island 2020,Horeshoe 2019, clear water 2019, chippewa, buffalo 2021, big sandy
#bemidjii 2022, Belle 2020, alexander 2021

#Mayb check: north lida 2021,  Pelican 2020, Gull 2019, green 2020, Big sand 2020

library(ggplot2)
library(readr)
library(dplyr)
library(ggrepel)
library(grDevices)

#subset data using ALL df from line 400
Check=subset(All,d15N>0)
Check$FFG[Check$Group=="Fish"]="Fish"
#summarizing data based on speciic taxa within the lake
gd<- Check %>%
  group_by(Community,FG,FFG, Group,Spp) %>%
  summarise(mean.C=mean(d13C, na.rm=T), sd.C=sd(d13C, na.rm=T), mean.N=mean(d15N, na.rm=T), sd.N=sd(d15N, na.rm=T), count=length(d13C))
#creating factor variables for plotting 
gd$FG=factor(gd$FG)
gd$FFG=factor(gd$FFG)
gd$Spp=factor(gd$Spp)
#remove NAs from sd columns (values where there is only one observation)
gd<-gd[!is.na(gd$sd.C),]
#standardize colors and shapes across plots 
library(RColorBrewer)
my_colors = data.frame("color.list"=brewer.pal(length(levels(gd$FFG)), "Paired"), "FFG"=levels(gd$FFG))
my_shapes = data.frame("shape.list"=c(1,19,2), "Habitat"=levels(gd$FG))
#plot by group, taxa, and lake 
lake.list=unique(gd$Community)
for(i in 1:length(lake.list)){
  current.data=subset(gd, Community==lake.list[i])
  current.colors=as.character(my_colors[my_colors$FFG%in%unique(current.data$FFG),1])
  current.shapes=my_shapes[my_shapes$Habitat%in%unique(current.data$FG),1]
  ggplot(current.data, aes(mean.C, mean.N, colour=FFG, shape=FG, ymax=mean.N+sd.N, ymin=mean.N-sd.N, xmax=mean.C+sd.C, xmin=mean.C-sd.C))+geom_point(size=2)+theme_bw()+geom_text_repel(aes(label=Spp), size=4)+geom_errorbar(alpha=.5)+geom_errorbarh(alpha=0.5)+theme(panel.grid = element_blank(), legend.position="bottom", axis.title = element_text(size=12), axis.text = element_text(size=12))+ggtitle(lake.list[i])+scale_colour_manual(values=current.colors)+scale_shape_manual(values=current.shapes)
  ggsave(paste(lake.list[i],"_isotope_data_by_taxa_includes_compromised.png", sep=""), height=10, width=8, units="in")
}


###### ALL LAKE DATA ######
#change format of datasets 
#BigLakes * does not have d15 and d13 columns using corrected C and N columns for now double chwck with DNR
BigLakes<-read.csv("~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/Large lakes 2017 2018/2017 and 2018 Organism and Isotope Data_Canon as of 08012019_flagged for baselines 10292019.csv")
#make community each lake, project ID combination 
BigLakes<-BigLakes %>%
mutate(Community=paste(Lake,"-",Year,sep="")) %>%
rename(Spp=Taxa_general,d13C=C13.corrected,d15N=N15.corrected)
BigLakes<-BigLakes[,c(8,7,17,18,22,23,25,31,34,36:39)] #getting columns of interest
#adult baseline df
BigFish_BL<-BigLakes[BigLakes$Adult.Baseline.Flag=="Y",]
BigFish_BL<-BigLakes[!is.na(BigLakes$Adult.Baseline.Flag),]
#yoy baseline df
YBigFish_BL<-BigLakes[BigLakes$YOY.Baseline.Flag=="Y",]
YBigFish_BL<-YBigFish_BL[!is.na(YBigFish_BL$YOY.Baseline.Flag),]
#fish dfs Y-yoy samples
BigFish<-BigLakes[BigLakes$Group=="Fish",]
BigF<-BigFish[!BigFish$YOY=="FALL",]
YBigFish<-BigF[!is.na(BigF$Spp),]
#sort for adults by getting na in yoy column 
BigFish<-BigFish[is.na(BigFish$YOY),]
BigFish<-BigFish[!is.na(BigFish$Spp),]
#cfix labels on adult dataset
BigFish$FFG<-"Consumer"
BL<-which(BigFish_BL$Habitat=="Littoral")
BP<-which(BigFish_BL$Habitat=="Offshore")
BigFish_L<-BigFish_BL[c(BL),]
BigFish_P<-BigFish_BL[c(BP),]
BigFish_P$FFG<-"Pelagic_BL"
BigFish_L$FFG<-"Littoral_BL"
BigFish_BL<-rbind(BigFish_L,BigFish_P)
#Fix labels on yoy data set 
YBigFish$FFG<-"Consumer"
YBigFish$Spp=paste(YBigFish$Spp,"Y",sep = "_")
BL<-which(YBigFish_BL$Habitat=="Littoral")
BP<-which(YBigFish_BL$Habitat=="Offshore")
YBigFish_L<-YBigFish_BL[c(BL),]
YBigFish_P<-YBigFish_BL[c(BP),]
YBigFish_P$FFG<-"Pelagic_BL"
YBigFish_L$FFG<-"Littoral_BL"
YBigFish_BL<-rbind(YBigFish_L,YBigFish_P)
#removing columns to match other SI df
BigFish_BL<-BigFish_BL[,c(2:9,12,13)]
BigFish<-BigFish[,c(2:9,12,13)]
BigLakes<-rbind(BigFish_BL,BigFish)
#do the same for YOY samples 
YBigFish_BL<-YBigFish_BL[,c(2:9,12,13)]
YBigFish<-YBigFish[,c(2:9,12,13)]
BigLakes_Y<-rbind(YBigFish_BL,YBigFish)

#combining the yoy and adult samples from the big lakes 
BigLakes<-rbind(BigLakes,BigLakes_Y)
BigLakes<-BigLakes[,c(1:8,10)]#remove baseline flag column

#removing FFG from ALL dataset to join big lakes
All$FFG<-All$FG
All<-All[,c(1:8,10)]#removing FFG column
Samples<-rbind(BigLakes,All)#7835
write.csv(Samples,"~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Samples_All_2025_FEB.csv")

Samples<-read.csv("Samples_All_2025_FEB.csv")
### checking N and C values ####
check<-Samples[Samples$d15N<0,]
Samples<-Samples[!Samples$d15N<0,]#7835

#using adult subset for moment need to fix species options for yoys 
#All_adult<-All_adult[!is.na(All_adult$Community),]
#All_adult<-All_adult[,c(1:5,7:ncol(All_adult))]
#All_adult$FFG<-All_adult$FG
#All_adult<-All_adult[,c(1:9,11)]
#Samples_A<-rbind(BigLakes,All_adult)
#check<-which(Samples_A$Spp=="TLC")
#Samples_A$Spp[c(check)]="CIS"

#write_csv(Samples_A,"~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Samples_reformatted.csv")

Clist<-unique(Samples$Community) #39 unique lake year combinations 
FishSamps<-Samples[Samples$Group=="Fish",]#6609 total observations 
df<-data.frame(rbind(table(FishSamps$Community,FishSamps$Spp)))
write.csv(df,"~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Species_Fish_all.csv")

#figure out which spp are represented across all lakes
FishSamps %>%
  group_by(Community) %>%
  summarise(count_unique = n_distinct(Spp)) %>%
  print(n=39) %>%
  order() #each lake year has min 9 different species-life stages present 


# *check SCU = scuplin?
Spp_detail<-read_csv(file = "~/Google Drive/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Species_Fish.csv")
df<-data.frame(Names=c(Spp_detail[39,]))

####### GO to tp Calcs script #######

##### littoral reliance  #####
#necessary packages 
library(MixSIAR)
library(tidyverse)
library(R2jags)

#bring in joined dataframe 
test<-read.csv(file = "~/Google Drive/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Samples_All_2025.csv")
test<-test[,c(2:ncol(test))]

#make separate df for consumer and baselines (source)
consumer<-test[test$Group=="Fish",]
baseline<-test[test$Group%in%c("Invertebrate","Zooplankton"),] 

#Leech 2018 not in consumer data remove community when BIG lakes is joined 
#baseline<-baseline[!baseline$Community=="Leech-2018",]

#change baseline labels to match trophic discrimination factor info 
baseline$Group <- gsub("Zooplankton", "Pelagic", baseline$Group)
baseline$Group <- gsub("Invertebrate", "Littoral", baseline$Group)

#store file for mix siar functions
write.csv(baseline,"~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/Base_Isotope_2025.csv")

#attach temp info (not used in this calculation- step can happen later)
#temp<-read.csv(file = "~/Google Drive/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/fish_tempCorrected.csv")
#temp <- temp %>% rename(Spp=species_code)
#temp <-temp[,c(1,2,4:7,18)] #get columns of interest 
#fish_temp<-full_join(temp,consumer,by="Spp")
#fish_temp<-fish_temp[!is.na(fish_temp$Spp),]
#fish_temp<-fish_temp[!is.na(fish_temp$Community),]

#store file for mix siar functions
write.csv(consumer,"~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/FishIso_2025.csv")
FishIso<-read.csv(file = "~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/FishIso_2025.csv")

#species check 
#Test<-FishIso[,c(7,16)]

#df<- Test %>% 
  #group_by(Community) %>% 
  #summarise(Spp = paste(unique(Spp), collapse = ', ')) %>%
  #print(n=39)

#SPP<-read.csv("~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Spp.csv")

#Loading mixture data for MixSIAR analysis 
mix_A = load_mix_data(filename="~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/FishIso_2025.csv", #file path for consumer 
                    iso_names=c("d13C","d15N"), #Identifies carbon isotope for analysis 
                    factors=c("Community"), #Lists factors
                    fac_random=c(TRUE), #Thermal guild is fixed factor, Lake is a random factor
                    fac_nested=c(FALSE), 
                    cont_effects=NULL) #Analysis does not contain continuous effects

mix_A$data<-mix_A$data[!is.na(mix_A$data$Community),]
#species subsets
mix_A$data<-mix_A$data[mix_A$data$Spp=="WAE",] #Walleye
mix_A$data<-mix_A$data[mix_A$data$Spp=="LMB",] #Largemouth Bass
mix_A$data<-mix_A$data[mix_A$data$Spp=="CIS",] #Cisco
mix_A$data<-mix_A$data[mix_A$data$Spp=="YEP",] #yellow perch
mix_A$data<-mix_A$data[mix_A$data$Spp=="BLG",] #bluegill
mix_A$data<-mix_A$data[mix_A$data$Spp=="NOP",] #pike
mix_A$data<-mix_A$data[mix_A$data$Spp=="SMB",] #smallmouth bass
mix_A$data<-mix_A$data[mix_A$data$Spp=="BLC",] #black crappie


#Loading source data for MixSIAR analysis 
source = load_source_data(filename="~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/Base_Isotope.csv",
                          conc_dep=FALSE, #No concentration dependence 
                          data_type = "raw",
                          mix = mix_A) #Telling the function which object is the fish data

model_file = "MixSIAR_model_blc.txt"   #Writes the JAGS model for exploration

#Loading in Trophic Enrichment Factors (ignore warning message)
discr = load_discr_data("~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/TDF.csv", mix_A)

process_err <- FALSE
resid_err <- TRUE
write_JAGS_model(model_file, #Supply the model name as an object
                 process_err = FALSE,
                 resid_err = TRUE,
                 mix_A, source) #Supplies mixture and source data to the model


##Model output is provided below as a .rds
###Run function gives the length of the model run, it is set to test. Model run of "extreme" was used for analysis
model.run = run_model(run="normal", #length of model run
                      mix_A, source, #mixture and source objects
                      discr = discr, #no trophic discrimination factors
                      model_file = "MixSIAR_model_blc.txt", #Model object
                      alpha.prior = c(1,1), #1 = uninformative prior for each baseline source
                      process_err,
                      resid_err) #Residual error included in analysis 
              

# Long model from analysis saved as RDS
#model.run <- readRDS("~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/AllLakesModel_Output.rds")

setwd("~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Plots and figures/")
                          
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics_blc",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density_blc",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot_blc",
                       sup_xy = FALSE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "xy_plot_blc",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics_blc",
                       indiv_effect = FALSE,
                       plot_post_save_png = TRUE,
                       plot_pairs_save_png = TRUE,
                       plot_xy_save_png = TRUE,
                       diag_save_ggmcmc = TRUE)
#Viewing output
output_diagnostics(model.run, mix_A, source, output_options)
df<-output_stats(model.run, mix_A, source, output_options)#returns summary statistics from a fit MixSIAR model
output_JAGS(model.run, mix_A, source, output_options)

plot_data(filename = 'isospace_plot_blc',plot_save_pdf = T,plot_save_png = F, mix_A, source, discr, return_obj = T)

BASELine<-read.csv("Base_Isotope_2025.csv")
L_BL<-BASELine[BASELine$FFG=="Littoral_BL",]
P_BL<-BASELine[BASELine$FFG=="Pelagic_BL",]

