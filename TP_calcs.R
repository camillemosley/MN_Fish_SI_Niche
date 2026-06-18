#Getting clean TP df to join with littoral reliance 
#CLM 11-21-2024, adapted from leech tp script from DL
#clear environment and plots 
rm(list=ls())
graphics.off()


##### Trophic position ####
library(tRophicPosition)
library(stringr)
library(dplyr)

setwd("~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data")

Samples<-read.csv("~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/Samples_reformatted.csv")
Samples<-Samples[!is.na(Samples$Community),]
Samples<-Samples[!Samples$Community=="Knife-2019",]#knife only baselines collected no fish samples to use 
Fish<-Samples[Samples$Group=="Fish",]
Fish<-Fish[!is.na(Fish$Spp),]
Fish$FFG<-"Consumer"
BL<-Samples[!Samples$Group=="Fish",]
Samples<-rbind(BL,Fish)
Samples<-Samples[!is.na(Samples$FFG),]
Samples$FFG<-gsub("Litoral_BL","Littoral_BL",Samples$FFG)

Cisco<-Samples[Samples$Spp=="CIS",] #make sure all spp have FFG and double check spp codes for small vs big lakes 

#fix Littoral BL spelling 

#make Isotope data object in list format for TP package 
#changed the Consumer column from SPP to FG to get plots for each lake -year
AllList<-extractIsotopeData(Samples,b1="Pelagic_BL", b2="Littoral_BL", baselineColumn = "FFG", consumersColumn = "Spp",
                            groupsColumn = "Community", d13C = "d13C",d15N = "d15N")
str(AllList)

write.csv(Samples,file = "MNLakesIsotopesJoined_09252025.csv")#saving output of cleaned SI samples

#get summary stats and plots for each community 
for (community in AllList) {
  print(summary(community))
  plot(community)
  screenIsotopeData(isotopeData = AllList$community,
                    density = "both",
                    consumer = "Consumer",
                    b1 = "Pelagic_BL",
                    b2 = "Littoral_BL",
                    legend = c(1.15, 1.15),
                    title = community,
                    xylim = NULL,)
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


###### WALLEYE #####
#trophic discrimantion factors
TDF_values <- TDF()

##### TP for individual species in lakes ####

#one species 
BL<- Samples %>%
  filter(FFG %in%c("Littoral_BL","Pelagic_BL"))
Consumer<-Samples %>%
  filter(FFG %in%c("Consumer"))

WAE<-Consumer[Consumer$Spp=="WAE",]
WAE_samps<-rbind(WAE,BL)

IsoData_Sp<-extractIsotopeData(WAE_samps,b1="Pelagic_BL", b2="Littoral_BL", baselineColumn = "FFG", consumersColumn = "FFG",
                               groupsColumn = "Community", d13C = "d13C",d15N = "d15N")
str(IsoData_Sp)

#remove NAs
IsoData <-lapply(IsoData_Sp,na.omit)

#define baysiean model
model.string <-jagsTwoBaselinesFull()
#intialize the model
model.WAE<-TPmodel(data=IsoData$`Alexander-2021-Consumer`, model.string = model.string)
#sampling and plotting posterior distribution 
samples.wae<-posteriorTP(model.WAE)
plot(samples.wae)
#save the posterior samples of trophic position (only first chain)
WAE_Alex.TP <- as.data.frame(samples.wae[[1]][,"TP"])$var1
# And then we combine it with the posterior samples from the second chain
WAE_Alex.TP <- c(WAE_Alex.TP,as.data.frame(samples.wae[[2]][,"TP"])$var1)
# And we check the length of the variable we created
length(WAE_Alex.TP)

# With this code we make the Species variable, including the species name
# "Orestias chungarensis"
Species <- c(rep("Walleye", length(WAE_Alex.TP)))
# and then we combine it with the posterior samples of Orestias trophic position
df <- data.frame(WAE_Alex.TP, Species)
# Changing the variable names of the dataframe "df"
colnames(df) <- c("TP", "Species")
# looking at output summary stats and plotting 
summary(df)
trophicDensityPlot(df)

#multiple species in one lake

All<-rbind(BL,Consumer)
#example Alexander lake
#Alex<- All %>% 
  #filter(Community %in%c("Alexander-2021"))
#make ISO object 
#IsoAlex<-extractIsotopeData(Alex,b1="Pelagic_BL", b2="Littoral_BL", baselineColumn = "FFG", consumersColumn = "Spp",
                               groupsColumn = "Community", d13C = "d13C",d15N = "d15N")
#check data format for model function
#str(IsoAlex)
#remove NAs
#IsoAlex <-lapply(IsoAlex,na.omit)
#IsoAlex<-IsoAlex[c(1:4,6:13)]#NA spp
#intialize the model for mutliple species (12)
#Tp.ALex<-multiSpeciesTP(IsoAlex, model = "twoBaselinesFull",
               #n.adapt = 20000, n.iter = 100000,
               #burnin = 50000, n.chains = 4, print = TRUE)
#save model output 
#saveRDS(Tp.ALex, "final_tpAlex.rds")

# numerical summary
#sapply(Tp.ALex$"TPs", quantile, probs = c(0.05, 0.5, 0.95)) %>% 
  #round(3)
#sapply(Tp.ALex$"Alphas", quantile, probs = c(0.05, 0.5, 0.95)) %>% 
  #round(3)

#str(Tp.ALex)

####### All samples TP calc #####

setwd("~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/TP models/")
lakes<-unique(Samples$Community)

for (lake in lakes) {
  tryCatch({ #try Catch to get errors and continue loop for lakes that have enough data points for each species
    
    df<- Samples %>% 
    filter(Community==lake)
    
  Iso<-extractIsotopeData(df,b1="Pelagic_BL", b2="Littoral_BL", baselineColumn = "FFG", consumersColumn = "Spp",
                              groupsColumn = "Community", d13C = "d13C",d15N = "d15N")
  #remove NAs
  Iso <-lapply(Iso,na.omit)
  
  #get TP and alpha (pelagic reliance) for all species in the lake
  Tp<-multiSpeciesTP(Iso, model = "twoBaselinesFull",
                          n.adapt = 20000, n.iter = 100000,
                          burnin = 50000, n.chains = 4, print = TRUE)
  #save model output 
  saveRDS(Tp,file= paste(lake,".rds",sep = "_"))},error=function(e){})
  
}

#Note pelagic is baseline 1 so aplha is pelagic reliance and littoral is 1-alpha

