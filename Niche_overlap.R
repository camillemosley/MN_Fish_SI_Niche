### Stable Isotope Data merge 
# CLM 9-2-2024

#clear environment and plots 
remove(list = ls())
graphics.off()

#libraries 
library(ggplot2)
library(SIBER)
library(ggpubr)
library(readxl)
library(dplyr)
library(stringr)
library(pheatmap)
library(plotly)
library(RColorBrewer)
library(tidyverse)
library(viridis)
library(ggrepel)
library(ggcorrplot)


##### SIBER analysis of ellipses/trophic niche for lake communities ####

#working directory in Gdrive 
setwd("/Users/cammosley/Google Drive/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/") #iMac
setwd("/Users/cammosley/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/TP models/")#macbook

#read in lake TP model outputs
TPFiles <- lapply(Sys.glob("*_.rds"), readRDS)#reading in all rds files 

#remove empty outputs that didn't have enough data
TPFiles<-TPFiles[c(1:14,16,17,19:32,34,36,37)] # remove 15, 18, 33, 35

TPFiles<-TPFiles[c(1:16,18:32)]
#TPFiles<-TPFiles[c(1:6,8,9,11:14,16,19:21,23:26,28:32,34,37)]#remove 7, 10, 15, 17,18, 22, 27, 33, 35,36
#make empty dataframe to store values from loop
TP_final<-data.frame(group=factor(),consumer=factor(),median=integer(),upper=integer(),
                   lower=integer(),alpha.median=integer(),alpha.upper=integer(),
                     alpha.lower=integer())

#get df values from all the lakes (need to check why some lakes do not have accurate rds files)
for (i in 1:length(TPFiles)){
  tryCatch({
  Data<-TPFiles[[i]][["df"]] #select the df of TP estimates from model output file
 # if ()#make loop go to i +1 if ncol < 12
  Data<-Data[,c(2:6,8:10)] #subset dataframe to have group=lake, consumer, median TP, median pelagic reliance 
  TP_final <- rbind(TP_final,Data)
  print(i)
  error=function(e){}})
}

TP_final$littoral_reliance<-1-TP_final$alpha.median
TP_final<-TP_final[!is.na(TP_final$consumer),]
TP_final<-TP_final[!is.na(TP_final$littoral_reliance),]
TP_final<-TP_final[TP_final$littoral_reliance>0,]
plot(TP_final$littoral_reliance,TP_final$median)
write.csv(TP_final,file="TrophicP_LittorialR_median04302026.csv")

TP_final<-read.csv("TrophicP_LittorialR_median04302026.csv")

#ggplot 
a<-ggplot(data = TP_final, aes(x=littoral_reliance,y=median, colour = group))+
     geom_point()+
     facet_wrap(~consumer)+ylab("Trophic Position")+xlab("Littoral Reliance")
b<-ggplot(data = TP_final, aes(x=littoral_reliance,y=median, colour = consumer, label=consumer))+
  geom_point()+
  facet_wrap(~group)+ylab("Trophic Position")+xlab("Littoral Reliance")+geom_text(hjust=0, vjust=0)

#list all species 
ALLFish<-merge(Physi,TP_final)

#dividing observations by waterbody 
split_fish <- split(TP_final, TP_final$group)

my_plot <- function(dat) {
  ggplot(dat, aes(x=littoral_reliance,y=median, colour = consumer, label=consumer)) +
    geom_point() +
    ylab("Trophic Position")+xlab("Littoral Reliance")+geom_text(hjust=0, vjust=0)
}

lapply(split_fish, my_plot)

#plots for species of interest 
FishList<-TP_final[TP_final$consumer%in%c("WAE","YEP","LMB","NOP","BLG","BLC","CIS","SMB"),]
ggplot(FishList, aes(x=littoral_reliance,y=median, colour = consumer, label=consumer)) +
  geom_point() +
  ylab("Trophic Position")+xlab("Littoral Reliance")+geom_text(hjust=0, vjust=0)

#color labels by ct max 
setwd("/Users/cammosley/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/")#macbook
Physi<-read.csv("fish_tempCorrected.csv")
Physi<-Physi[,c(4,5:7,18)] # get columns of interest 
Physi$consumer<-Physi$species_code
FishCT<-merge(FishList,Physi)
FishCT$Temperature.preference.class<-as.factor(FishCT$Temperature.preference.class)

ggplot(FishCT, aes(x=littoral_reliance,y=median, colour = FTP, label=group)) +
  geom_point() +
  ylab("Trophic Position")+xlab("Littoral Reliance")+geom_text(hjust=0, vjust=0)+
  scale_color_viridis_c(option = "rocket") +facet_wrap(~consumer)+theme_dark()

#dividing observations by waterbody 
split_fish <- split(FishList, FishList$group)
pdf("All_Fish_Plots.pdf", width = 8, height = 6)   # open PDF device                                         

my_plot <- function(dat) {
  lakename<-dat$group
  ggplot(dat, aes(x=littoral_reliance,y=median, colour = consumer, label=consumer)) +
    geom_jitter() +
    ylim(2,5)+
    xlim(0,1)+
    labs(title=lakename)+
    ylab("Trophic Position")+xlab("Littoral Reliance")+geom_text(hjust=0, vjust=0)
}

plots <- lapply(split_fish, my_plot)               # create all plots
for(p in plots) print(p)                           # write each plot to the PDF

dev.off() 

#b<-ggplot(data = FishList, aes(x=littoral_reliance,y=median, colour = consumer, label=consumer))+
  #geom_point()+ylim(0,5.5)+
  #facet_wrap(~group)+ylab("Trophic Position")+xlab("Littoral Reliance")+geom_text(hjust=0, vjust=0)

##### Figure 2 in paper #####
#remove year from the lake names
FishCT$group<-substr(FishCT$group,1,nchar(FishCT$group)-5)

b<-ggplot(data = FishCT, aes(x=littoral_reliance,y=median, colour = group))+
  geom_jitter(size=3)+
  theme_bw()+
  theme(axis.title = element_text(face = "bold"))+
  facet_wrap(~consumer)+
  ylab("Trophic Position")+
  xlab("Littoral Reliance")+
  labs(color= "Lake Community")

#look at littoral reliance values 
FishList_NOCis<-FishList[!FishList$consumer=="CIS",]
CIS<-FishList[FishList$consumer=="CIS",]
sd_LR<-sd(FishList$littoral_reliance)
sd_TP<-sd(FishList$median)

### order plot by thermal metric then plots spceis info across lakes
# 1) compute one FTP value per consumer (use median)
consumer_order <- FishCT %>%
  group_by(consumer) %>%
  summarise(FTP_consumer = median(FTP, na.rm = TRUE)) %>%
  arrange(FTP_consumer) %>%
  pull(consumer)

# 2) set factor order
FishCT <- FishCT %>%
  mutate(consumer = factor(consumer, levels = consumer_order))

# 3) Determine global axis limits
x_min <- min(FishCT$littoral_reliance)
x_max <- max(FishCT$littoral_reliance)
y_min <- min(FishCT$median)
y_max <- max(FishCT$median)

# 4) Plot with common limits
b <- ggplot(FishCT, aes(x = littoral_reliance, y = median, colour = group))+
  theme_bw() +geom_point()+
  facet_wrap(~ consumer) +                  # ← same axis limits
  xlab("Littoral Reliance") +
  ylab("Trophic Position") +
  coord_cartesian(xlim = c(x_min, x_max),
                  ylim = c(y_min, y_max))   # ← apply global axis limits

b


LMB<-TP_final[TP_final$consumer=="LMB",]
WAE<-TP_final[TP_final$consumer=="WAE" & TP_final$group%in%(LMB$group),]
fig1<-plot_ly(WAE, x=WAE$littoral_reliance,y=WAE$median, type="scatter",color = ~group, mode="markers", colors = brewer.pal(n = min(length(unique(WAE$group)), 12), "Set3"))%>%
  layout(xaxis = list(title = 'Littoral Reliance'), yaxis = list(title = 'Trophic Position'))

fig2<-plot_ly(LMB, x=LMB$littoral_reliance,y=LMB$median, type="scatter",color = ~group,mode="markers", colors = brewer.pal(n = min(length(unique(WAE$group)), 12), "Set3"))%>%
  layout(xaxis = list(title = 'Littoral Reliance'), yaxis = list(title = 'Trophic Position'))

fig<-subplot(fig1,fig2) %>% layout(title = 'Trophic position and Littoral Reliance',
                               plot_bgcolor='#e5ecf6', 
                               xaxis = list( 
                                 zerolinecolor = '#ffff', 
                                 zerolinewidth = 2, 
                                 gridcolor = 'ffff'), 
                               yaxis = list( 
                                 zerolinecolor = '#ffff', 
                                 zerolinewidth = 2, 
                                 gridcolor = 'ffff'))
annotations = list( 
  list( 
    x = 0.2,  
    y = 1.0,  
    text = "Walleye",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ),  
  list( 
    x = 0.8,  
    y = 1,  
    text = "Largemouth Bass",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ))  

fig<-fig %>% layout(annotations = annotations)

w<-ggplot(data=WAE,aes(x=littoral_reliance,y=median,color = group))+geom_point()+ylab("Trophic Position")+xlab("Littoral Reliance")+title(main = "Walleye")+theme(legend.position="none")

l<-ggplot(data=LMB,aes(x=littoral_reliance,y=median,color = group))+geom_point()+ylab("Trophic Position")+xlab("Littoral Reliance")+title(main = "Largemouth Bass")+theme(legend.position="none")
ggarrange(w,l,
          labels = c("A", "B"))

#### TP and LR plots clean ####

TP_species<-TP_final[TP_final$consumer%in%c("WAE","YEP","LMB","NOP","BLG","BLC","CIS"),]

# Summarize to one row per group (median values)
TP_summary <- TP_species %>%
  group_by(group) %>%
  summarise(
    median_TP = median(median, na.rm = TRUE),               # median trophic position
    median_LR = median(littoral_reliance, na.rm = TRUE),    # median littoral reliance
    n_species = n_distinct(consumer)                        # optional: number of consumers per group
  )

# Check summary
print(head(TP_summary))

# Biplot: median trophic position vs. littoral reliance
TP_biplot <- ggplot(TP_summary, aes(x = median_LR, y = median_TP)) +
  geom_point(
    aes(size = n_species, color = median_TP),  # color = trophic position gradient
    alpha = 0.8
  ) +
  ggrepel::geom_text_repel(
    aes(label = group),
    hjust = -0.1,
    vjust = 0.5,
    size = 4,
    check_overlap = TRUE
  ) +
  scale_color_viridis_c(option = "plasma") +
  scale_size_continuous(range = c(2, 6)) +
  labs(
    title = "Median Trophic Position vs. Littoral Reliance (per species)",
    x = "Median Littoral Reliance",
    y = "Median Trophic Position",
    color = "Trophic Position",
    size  = "No. of Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

TP_biplot

#make biplots for each species 
TP_summary <- TP_species %>%
  group_by(group, consumer) %>%
  summarise(
    median_TP = median(median, na.rm = TRUE),
    median_LR = median(littoral_reliance, na.rm = TRUE),
    n_species = n_distinct(consumer),
    .groups = "drop"
  )

xlim_vals <- range(TP_summary$median_LR, na.rm = TRUE)
ylim_vals <- range(TP_summary$median_TP, na.rm = TRUE)

species_list <- unique(TP_summary$consumer)

pdf("Species_Biplots.pdf", width = 6, height = 5)

for (sp in species_list) {
  
  df_sp <- TP_summary %>% filter(consumer == sp)
  
  p <- ggplot(df_sp, aes(x = median_LR, 
                         y = median_TP, 
                         color = median_TP)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(aes(label = group), size = 3, max.overlaps = 20) +
    scale_color_viridis(option = "plasma") +
    xlim(xlim_vals) +
    ylim(ylim_vals) +
    labs(
      title = paste("Trophic Position vs. Littoral Reliance for", sp),
      x = "Median Littoral Reliance",
      y = "Median Trophic Position",
      color = "TP"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  print(p)
}

dev.off()

#change tp column name to not confuse median call in plotting code
FishList<- FishList %>%
  rename(TP=median)
#look at the median values by group to make sure input to ggplot is correct 
FishList %>% 
  group_by(consumer) %>%
  summarise(medianLR=median(littoral_reliance),medianTP=median(TP)) %>%
  arrange(desc(medianLR))

# Compute median littoral reliance by consumer
median_order <- FishList %>%
  group_by(consumer) %>%
  summarise(med_littoral = median(littoral_reliance, na.rm = TRUE)) %>%
  arrange(med_littoral) %>%
  pull(consumer)

# Apply that ordering to the factor
FishList <- FishList %>%
  mutate(consumer = factor(consumer, levels = median_order))

# Identify outliers per consumer
FishList <- FishList %>%
  group_by(consumer) %>%
  mutate(
    Q1 = quantile(littoral_reliance, 0.25, na.rm = TRUE),
    Q3 = quantile(littoral_reliance, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = littoral_reliance < (Q1 - 1.5 * IQR) |
      littoral_reliance > (Q3 + 1.5 * IQR)
  ) %>%
  ungroup()

# Figure X: Manuscript  ####
Lr_plot <- ggplot(FishList, aes(x = consumer, y = littoral_reliance)) +
  # gray points for non-outliers only
  geom_jitter(
    data = subset(FishList, !is_outlier),
    color = "gray60",
    width = 0.2,
    alpha = 0.5,
    size = 1.5
  ) +
  # boxplots with red outliers on top
  geom_boxplot(
    aes(fill = consumer),
    alpha = 0.7,
    outlier.colour = "red",
    outlier.shape = 21,
    outlier.size = 2,
    outlier.stroke = 2
  ) +
  # Add group labels for outliers only
  geom_text_repel(
    data = subset(FishList, is_outlier),
    aes(label = group),
    color = "black",
    size = 3
  ) +
  labs(
    #title = "Littoral Reliance values across study lakes",
    x = "Species (ordered by median Littoral Reliance)",
    y = "Littoral Reliance"
  ) +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0, unit = "pt")), # Increase space for the x-axis title (add margin to the top)
    legend.position = "none"
  ) +
  scale_fill_viridis_d(option = "viridis")

Lr_plot

FishList <- FishList %>%
  group_by(consumer) %>%
  mutate(mean_litt = mean(littoral_reliance, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(consumer = reorder(consumer, mean_litt))

tukey_letters <- multcompLetters4(anova_model, tukey)

letters_vec <- tukey_letters$consumer$Letters

letters_df <- data.frame(
  consumer = names(letters_vec),
  letters = letters_vec,
  stringsAsFactors = FALSE
)

summary_df <- FishList %>%
  group_by(consumer) %>%
  summarise(
    mean_litt = mean(littoral_reliance, na.rm = TRUE),
    lower = mean_litt - qt(0.975, df = n()-1) * sd(TP)/sqrt(n()),
    upper = mean_litt + qt(0.975, df = n()-1) * sd(TP)/sqrt(n()),
    .groups = "drop"
  ) %>%
  left_join(letters_df, by = "consumer")

LR_plot <- ggplot(summary_df, aes(x = consumer, y = mean_litt)) +
  
  # CI error bars
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    size = 0.8
  ) +
  
  # mean points
  geom_point(size = 4, aes(color = consumer)) +
  
  # significance letters
  geom_text(
    aes(y = upper + 0.1, label = letters),
    size = 4
  ) +
  
  # optional raw data (faded)
  geom_jitter(
    data = FishList,
    aes(x = consumer, y = littoral_reliance),
    inherit.aes = FALSE,
    width = 0.2,
    alpha = 0.3,
    color = "gray60"
  ) +
  
  labs(
    x = "Species (ordered by mean Littoral Reliance)",
    y = "Mean Littoral Reliance ± 95% CI"
  ) +
  
  theme_minimal() +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0, unit = "pt")),
        legend.position = "none"
  ) +
  
  scale_color_viridis_d(option = "viridis")

LR_plot

# Compute median trophic position by consumer
median_order <- FishList %>%
  group_by(consumer) %>%
  summarise(med_TP = median(TP, na.rm = TRUE)) %>%
  arrange(med_TP) %>%
  pull(consumer)

# Apply that ordering to the factor
FishList <- FishList %>%
  mutate(consumer = factor(consumer, levels = median_order))

# Identify outliers per consumer group
FishList <- FishList %>%
  group_by(consumer) %>%
  mutate(
    Q1 = quantile(TP, 0.25, na.rm = TRUE),
    Q3 = quantile(TP, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = TP < (Q1 - 1.5 * IQR) |
      TP > (Q3 + 1.5 * IQR)
  ) %>%
  ungroup()

# Figure X: Manuscript  ####
TP_plot <- ggplot(FishList, aes(x = consumer, y = TP)) +
  # Gray jitter points for non-outliers
  geom_jitter(
    data = subset(FishList, !is_outlier),
    color = "gray60",
    width = 0.2,
    alpha = 0.5,
    size = 1.5
  ) +
  # Boxplots with red outliers
  geom_boxplot(
    aes(fill = consumer),
    alpha = 0.7,
    outlier.colour = "red",
    outlier.shape = 21,
    outlier.size = 2,
    outlier.stroke = 2
  ) +
  stat_summary(
    fun.data = median_hilow,
    fun.args = list(conf.int = 0.95),
    geom = "errorbar",
    width = 0.2
  ) +
  # Add group labels for outliers only
  geom_text_repel(
    data = subset(FishList, is_outlier),
    aes(label = group),
    color = "black",
    size = 3,
    nudge_y = 0.05,    # vertical offset
    max.overlaps = Inf # allow multiple labels if needed
  ) +
  labs(
       x = "Species (ordered by median Trophic Position)",
       y = "Trophic Postion") +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0, unit = "pt")),
    legend.position = "none"
  ) +
  scale_fill_viridis_d(option = "viridis")

TP_plot

FishList <- FishList %>%
  group_by(consumer) %>%
  mutate(mean_TP = mean(TP, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(consumer = reorder(consumer, mean_TP))

tukey_letters <- multcompLetters4(anova_model, tukey)

letters_vec <- tukey_letters$consumer$Letters

letters_df <- data.frame(
  consumer = names(letters_vec),
  letters = letters_vec,
  stringsAsFactors = FALSE
)

summary_df <- FishList %>%
  group_by(consumer) %>%
  summarise(
    mean_TP = mean(TP, na.rm = TRUE),
    lower = mean_TP - qt(0.975, df = n()-1) * sd(TP)/sqrt(n()),
    upper = mean_TP + qt(0.975, df = n()-1) * sd(TP)/sqrt(n()),
    .groups = "drop"
  ) %>%
  left_join(letters_df, by = "consumer")

TP_plot <- ggplot(summary_df, aes(x = consumer, y = mean_TP)) +
  
  # CI error bars
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    size = 0.8
  ) +
  
  # mean points
  geom_point(size = 4, aes(color = consumer)) +
  
  # significance letters
  geom_text(
    aes(y = upper + 0.1, label = letters),
    size = 4
  ) +
  
  # optional raw data (faded)
  geom_jitter(
    data = FishList,
    aes(x = consumer, y = TP),
    inherit.aes = FALSE,
    width = 0.2,
    alpha = 0.3,
    color = "gray60"
  ) +
  
  labs(
    x = "Species (ordered by mean Trophic Position)",
    y = "Mean Trophic Position ± 95% CI"
  ) +
  
  theme_minimal() +
  theme(axis.title = element_text(face = "bold"),
               axis.text.x = element_text(angle = 45, hjust = 1),
               axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0, unit = "pt")),
               legend.position = "none"
  ) +
  
  scale_color_viridis_d(option = "viridis")

TP_plot


#format data for SIBER 
# iso 1 =x, iso 2 = y
#iso 1 = %C , iso 2 = %N

# load in the included demonstration dataset
data("demo.siber.data")

#read in  joined in isotope data #####
IsoData<-read.csv("/Users/cammosley/Google Drive/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/LAKE_IDs_Samples_All_2025_FEB.csv")
#subset for columns needed for siber and reformat df to make siber object 
IsoData<-IsoData[ ,c(3,8,9,10)]
IsoData<-IsoData[!IsoData$YOY=="Y",]#no YOY
###### several lakes #####
IsoData<-IsoData[IsoData$Spp%in%c("WAE","YEP","LMB","BLG","SMB","BLC","NOP","RKB"),] 
IsoData<-IsoData[IsoData$Spp%in%c("WAE","YEP","LMB","BLG","SMB","BLC","NOP","RKB"),] 
IsoData<-IsoData[IsoData$Community%in%c("Belle-2020"),]#
#The header names in your data file must match 
#c("iso1", "iso2", "group", "community") exactly
IsoData<-IsoData %>% 
  rename(iso1=d13C, iso2=d15N, group=Spp, community=Community)

IsoData<-data.frame(iso1=IsoData$iso1,iso2=IsoData$iso2,group=IsoData$group,community=IsoData$community)
#create SIBER object
SI<-createSiberObject(IsoData)
#t least one of your groups has less than 5 observations.
#The absolute minimum sample size for each group is 3

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)

# plot the raw data
par(mfrow=c(1,1))
plotSiberObject(SI,
                ax.pad = 2, 
                hulls = T, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille')
)

legend("topleft",
        as.character(paste("Group ",(SI$original.data$group))),
       pch=19,
       col=1:length(SI$original.data$group))
legend("topright",
       as.character(paste("Lake ",(SI$original.data$community))),
       pch=19,
       col=1:length((SI$original.data$community)))


# add the confidence interval of the means to help locate
# the centre of each data cluster
plotGroupEllipses(SI, n = 100, p.interval = 0.95,
                  ci.mean = T, lty = 1, lwd = 2)

# Fit the Bayesian models

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
ellipses.posterior <- siberMVN(SI, parms, priors)

# extract the posterior means
mu.post <- extractPosteriorMeans(SI, ellipses.posterior)

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)


# Visualise the first community
# drop the 3rd column of the posterior which is TA using -3.

#Communites are labled 1-6 in the output metrics tp get info for community 2 put layman.B[[2][ , -3], 
#for X community input is layman.B[[X][ , -3]
siberDensityPlot(layman.B[[1]][ , -3], 
                 xticklabels = colnames(layman.B[[1]][ , -3]), 
                 bty="L", ylim = c(0,20))

# add the ML estimates (if you want). Extract the correct means 
# from the appropriate array held within the overall array of means.
comm1.layman.ml <- laymanMetrics(SI$ML.mu[[1]][1,1,],
                                 SI$ML.mu[[1]][1,2,]
)

# again drop the 3rd entry which relates to TA
points(1:5, comm1.layman.ml$metrics[-3], 
       col = "red", pch = "x", lwd = 2)

siberDensityPlot(layman.B[[2]][ , -3], 
                 xticklabels = colnames(layman.B[[2]][ , -3]), 
                 bty="L", ylim = c(0,20))

# add the ML estimates. (if you want) Extract the correct means 
# from the appropriate array held within the overall array of means.
comm2.layman.ml <- laymanMetrics(SI$ML.mu[[2]][1,1,],
                                 SI$ML.mu[[2]][1,2,]
)
points(1:5, comm2.layman.ml$metrics[-3], 
       col = "red", pch = "x", lwd = 2)

# Alternatively, pull out TA from both and aggregate them into a 
# single matrix using cbind() and plot them together on one graph.

# go back to a 1x1 panel plot
par(mfrow=c(1,1))

# Now we only plot the TA data. We could address this as either
# layman.B[[1]][, "TA"]
# or
# layman.B[[1]][, 3]
siberDensityPlot(cbind(layman.B[[1]][ , "TA"], 
                       layman.B[[2]][ , "TA"],
                       layman.B[[3]][ , "TA"],
                       layman.B[[4]][ , "TA"],
                       layman.B[[5]][ , "TA"]),
                 xticklabels = c("Community 1", "Community 2","Community 3",
                                 "Community 4",
                                 "Community 5"), 
                 bty="L", ylim = c(0, 10),
                 las = 1,
                 ylab = "TA - Convex Hull Area",
                 xlab = "")

TA1_lt_TA2 <- sum(layman.B[[1]][,"TA"] < 
                    layman.B[[2]][,"TA"]) / 
  length(layman.B[[1]][,"TA"])

print(TA1_lt_TA2)

#### 2 lakes ######
IsoData<-read.csv("/Users/cammosley/Google Drive/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/MNLakesIsotopesJoined_112124.csv")
#subset for columns needed for siber and reformat df to make siber object 
IsoData<-IsoData[,c(3,8,9,11)]
IsoData<-IsoData[IsoData$Spp%in%c("WAE","YEP","LMB","NOP","BLG","SMB","BLC"),] #7 groups 
#The header names in your data file must match 
#c("iso1", "iso2", "group", "community") exactly
IsoData<-IsoData %>% 
  rename(iso1=d13C, iso2=d15N, group=Spp, community=Community)
IsoData<-data.frame(iso1=IsoData$iso1,iso2=IsoData$iso2,group=IsoData$group,community=IsoData$community)
#create SIBER object
SI<-createSiberObject(IsoData)
#t least one of your groups has less than 5 observations.
#The absolute minimum sample size for each group is 3

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

# plot the raw data
par(mfrow=c(1,1))
plotSiberObject(SI,
                ax.pad = 2, 
                hulls = T, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille')
)

# add the confidence interval of the means to help locate
# the centre of each data cluster
plotGroupEllipses(SI, n = 100, p.interval = 0.95,
                  ci.mean = T, lty = 1, lwd = 2)

# Fit the Bayesian models

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
ellipses.posterior <- siberMVN(SI, parms, priors)

# extract the posterior means
mu.post <- extractPosteriorMeans(SI, ellipses.posterior)

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)

# drop the 3rd column of the posterior which is TA using -3.

#Communites are labled 1-6 in the output metrics tp get info for community 2 put layman.B[[2][ , -3], 
#for X community input is layman.B[[X][ , -3]
siberDensityPlot(layman.B[[1]][ , -3], 
                 xticklabels = colnames(layman.B[[1]][ , -3]), 
                 bty="L", ylim = c(0,20))

# add the ML estimates (if you want). Extract the correct means 
# from the appropriate array held within the overall array of means.
comm1.layman.ml <- laymanMetrics(SI$ML.mu[[1]][1,1,],
                                 SI$ML.mu[[1]][1,2,]
)

# again drop the 3rd entry which relates to TA
points(1:5, comm1.layman.ml$metrics[-3], 
       col = "red", pch = "x", lwd = 2)

siberDensityPlot(layman.B[[2]][ , -3], 
                 xticklabels = colnames(layman.B[[2]][ , -3]), 
                 bty="L", ylim = c(0,20))

# add the ML estimates. (if you want) Extract the correct means 
# from the appropriate array held within the overall array of means.
comm2.layman.ml <- laymanMetrics(SI$ML.mu[[2]][1,1,],
                                 SI$ML.mu[[2]][1,2,]
)
points(1:5, comm2.layman.ml$metrics[-3], 
       col = "red", pch = "x", lwd = 2)


# go back to a 1x1 panel plot
par(mfrow=c(1,1))

# Now we only plot the TA data. We could address this as either
# layman.B[[1]][, "TA"]
# or
# layman.B[[1]][, 3]
siberDensityPlot(cbind(layman.B[[1]][ , "TA"], 
                       layman.B[[2]][ , "TA"]),
                 xticklabels = c("Community 1", "Community 2"), 
                 bty="L", ylim = c(0, 10),
                 las = 1,
                 ylab = "TA - Convex Hull Area",
                 xlab = "")

TA1_lt_TA2 <- sum(layman.B[[1]][,"TA"] < 
                    layman.B[[2]][,"TA"]) / 
  length(layman.B[[1]][,"TA"])

print(TA1_lt_TA2)

##### All Lakes ######

IsoData<-read.csv("/Users/cammosley/Google Drive/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/Samples_All_2025_FEB.csv")
#subset for columns needed for siber and reformat df to make siber object 
IsoData<-IsoData[!IsoData$YOY=="Y",]#no YOY
IsoData<-IsoData[,c(3,8:10)]
IsoData<-IsoData[IsoData$Spp%in%c("WAE","YEP","LMB","BLG","SMB","BLC","NOP","CIS"),] 
#IsoData<-IsoData[IsoData$Spp%in%c("WAE","YEP","LMB","BLG","SMB","BLC","NOP"),]
#IsoData<-IsoData[IsoData$Spp%in%c("WAE","YEP","LMB","NOP","BLG","SMB","CIS","LAT","LKW"),]
#The header names in your data file must match 
#c("iso1", "iso2", "group", "community") exactly
IsoData<-IsoData %>% 
  dplyr::rename(iso1=d13C, iso2=d15N, group=Spp, community=Community)
IsoData<-data.frame(iso1=IsoData$iso1,iso2=IsoData$iso2,group=IsoData$group,community=IsoData$community)

#removing samples sizes that are too small
#probably remove LKW 
#IsoData<-IsoData[IsoData$community%in%c("Washington-2020","Alexander-2021"),]
IsoData <-IsoData[ order(IsoData$group),]#grouping the species codes and sorting them alphabetically 

sample_size <- IsoData %>% 
  group_by(group, community) %>%
  filter(n() > 4) 

A<-sample_size[sample_size$community=="Alexander-2021",]
#BM<-sample_size[sample_size$community=="Bemidji-2022",]
BE<-sample_size[sample_size$community=="Belle-2020",]     
BC<-sample_size[sample_size$community=="Big Cormorant-2019",]  
BS<-sample_size[sample_size$community=="Big Sandy-2021",]
#EV<-sample_size[sample_size$community=="East Vermilion-2018",]  
Leech22<-sample_size[sample_size$community=="Leech-2022",]  
WA<-sample_size[sample_size$community=="Washington-2020",]  
CA<-sample_size[sample_size$community=="Cass-2018",]
#KA<-sample_size[sample_size$community=="Kabetogama-2018",]
#Leech17<-sample_size[sample_size$community=="Leech-2017",]
LOTW<-sample_size[sample_size$community=="Lake of the Woods-2018",]
ML<-sample_size[sample_size$community=="Mille Lacs-2017",]
#LR<-sample_size[sample_size$community=="Lower Red-2017",]
#UR<-sample_size[sample_size$community=="Upper Red-2017",]
RA<-sample_size[sample_size$community=="Rainy-2018",]
#WV<-sample_size[sample_size$community=="West Vermilion-2018",]
WI<-sample_size[sample_size$community=="Winnie-2018",]
BSA<-sample_size[sample_size$community=="Big Sand-2020",]
BU<-sample_size[sample_size$community=="Buffalo-2021",]
CH<-sample_size[sample_size$community=="Chippewa-2021",]
CL<-sample_size[sample_size$community=="Clearwater-2019",]
GR<-sample_size[sample_size$community=="Green-2020",]
GU<-sample_size[sample_size$community=="Gull-2019",]
HO<-sample_size[sample_size$community=="Horseshoe-2019",]
IS<-sample_size[sample_size$community=="Island-2020",]
KO<-sample_size[sample_size$community=="Koronis-2019",]
LB<-sample_size[sample_size$community=="Little Boy-2021",]
ME<-sample_size[sample_size$community=="Melissa-2019",]
NL<-sample_size[sample_size$community=="North Lida-2021",]
PE<-sample_size[sample_size$community=="Pelican-2020",]
PL<-sample_size[sample_size$community=="Plantagenet-2022",]
PO<-sample_size[sample_size$community=="Potato-2022",]
RO<-sample_size[sample_size$community=="Round-2022",]
RU<-sample_size[sample_size$community=="Rush-2022",]
SH<-sample_size[sample_size$community=="Shamineau-2021",]
ST<-sample_size[sample_size$community=="Steamboat-2021",]
TE<-sample_size[sample_size$community=="Tenmile-2020",]

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
ellipses.posterior <- siberMVN(SI, parms, priors)

# extract the posterior means
mu.post <- extractPosteriorMeans(SI, ellipses.posterior)

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)

#overlap between largemouth bass and walleye 
#lmb group 3 and wae group 6 
ellipse1<-"Alexander-2021.SMB"
ellipse2<-"Alexander-2021.WAE"
# The overlap of the maximum likelihood fitted standard ellipses are 
# estimated using
sea.overlap <- maxLikOverlap(ellipse1, ellipse2, SI, 
                             p.interval = NULL, n = 100)
# the overlap betweeen the corresponding 95% prediction ellipses is given by:
ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse2, SI, 
                                   p.interval = 0.95, n = 100)
# so in this case, the overlap as a proportion of the non-overlapping area of 
# the two ellipses, would be
prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] + 
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])
#and the corresponding Bayesian estimates for the overlap between the 
# 95% ellipses is given by:
bayes95.overlap <- bayesianOverlap(ellipse1, ellipse2, ellipses.posterior,
                                   draws = 100, p.interval = 0.95, n = 100)
# and we can calculate the corresponding credible intervals using (Andrew L Jackson)
overlap.credibles <- lapply(
  as.data.frame(bayes95.overlap), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(overlap.credibles)

# the proportion of ellipse 2 that overlaps with ellipse 1 
prop.of.second <- as.numeric(ellipse95.overlap["overlap"] / ellipse95.overlap["area.2"])
print(prop.of.second)

#the proportion of 1 and 2 that overlap with each other
prop.of.both <- as.numeric(ellipse95.overlap["overlap"] / (ellipse95.overlap["area.1"] + ellipse95.overlap["area.2"]))
print(prop.of.both)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = c("BLC","BLG","NOP","SMB","WAE","YEP"), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('permille' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "Ellipse Area of fish in Lake Alexander - 2021", 
                 ylims = c(0,5))

# Calculate the various Layman metrics on each of the communities.
community.ML <- communityMetricsML(SI) 
print(community.ML)

# Calculate sumamry statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(SI)
print(group.ML)


###### GGplot ####
write.csv(IsoData,"IsoPlotting_data.csv")
# Get unique lakes
communities <- unique(IsoData$community)

# Create plots for each community
for (comm in communities) {
  # Subset data for the community
  df_subset <- IsoData[IsoData$community == comm, ]
  p.ell <- 0.95 #95% confidence interval around ellipse
  # Create ggplot for each data point in the community
  p <- ggplot(df_subset, aes(x = iso1, y = iso2, color = group)) +
    geom_point() +
    stat_ellipse(aes(group = interaction(group, community), 
                     color = group), 
                 alpha = 0.25, 
                 level = p.ell,
                 type = "norm",
                 geom = "polygon") +
    labs(
      x = expression(delta^13*C~'\u2030'),
      y = expression(delta^15*N~'\u2030'),
      color = "Species"
    ) +
    theme_minimal()+
    ggtitle(paste0("Plot of ",comm))
  
  # Print the plot
  print(p)
  
  # Optional: Save each plot as an image file
  ggsave(filename = paste0("Ellipse_", comm, ".png"), plot = p, width = 6, height = 4, dpi = 300)
}


first.plot<-ggplot(IsoData, aes(x = iso1, y = iso2, color = group)) +
  geom_point()+ facet_wrap(~community)

ggplot(IsoData, aes(x = iso1, y = iso2, color = group)) +
  geom_point() +
  stat_ellipse() +
  labs(
    x = expression(delta^13*C~'\u2030'),
    y = expression(delta^15*N~'\u2030'),
    color = "Species"
  ) +
  theme_minimal()


##### ellipse area calcs ######

#make empty dfs to store SIBER output 
SIBER_comm<-data.frame(community=character(),TA_C=numeric(),MNND=numeric(), SDNND=numeric())
SIBER_group<-data.frame(matrix(ncol = 0, nrow = 3, dimnames = list(c("TA","SEA", "SEAc"), NULL)))

df_list<-list(A,BE,BC, BS,Leech22,
              WA,CA,ML,RA,WI,BSA,BU,
              CH,CL,GR,GU,HO,IS,KO,LB,ME,NL,PE,PL,PO,
              RO,RU,SH,ST,TE)

#loop through SI lake dataframes to calculate community and population ellipse metrics 
for (i in seq_along(df_list)){
  df <- df_list[[i]]#get individual lake
  
  # Check if dataframe is empty
  if (nrow(df) == 0) {
    cat("Skipping empty dataframe", i, "\n")
    next  # Skip to the next dataframe
  }
  
  SI<-createSiberObject(df)
  
  # Calculate the various Layman metrics on each of the communities.
  community.ML <- communityMetricsML(SI) 
  
  # Calculate sumamry statistics for each group: TA, SEA and SEAc
  group.ML <- groupMetricsML(SI)
  
  #add results to dataframe
  results_comm<-data.frame(community=colnames(community.ML),TA_C=community.ML[3,],MNND=community.ML[5,], SDNND=community.ML[6,])
  results_group<-data.frame(group.ML)
  
  #append to master df 
  SIBER_comm<-rbind(results_comm,SIBER_comm)
  SIBER_group<-cbind(results_group,SIBER_group)
  
  #ploting ellipses using jags model
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
  ellipses.posterior <- siberMVN(SI, parms, priors)
  
  SEA <- siberEllipses(ellipses.posterior)
  
  siberDensityPlot(SEA, xticklabels = c(unique(df$group)), 
                   xlab = c("Community | Group"),
                   ylab = expression("Standard Ellipse Area " ('permille' ^2) ),
                   bty = "L",
                   las = 1,
                   main = paste("Ellipse Area by group in",df$community[1]), 
                   ylims = c(0,(max(group.ML[2,])+4.5))) #getting range based of max and min values for each lake
  
  community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
  group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
  group.hull.args      <- list(lty = 2, col = "grey20")
  plotSiberObject(SI,
                  ax.pad = 2, 
                  hulls = T, community.hulls.args, 
                  ellipses = F, group.ellipses.args,
                  group.hulls = F, group.hull.args,
                  bty = "L",
                  iso.order = c(1,2),
                  xlab = expression({delta}^13*C~'permille'),
                  ylab = expression({delta}^15*N~'permille')
  )
  legend("topleft",
         as.character(paste("Group ",unique(SI$original.data$group))),
         pch=19,
         col=1:length(unique(SI$original.data$group)))
  legend("topright",
         as.character(paste("Lake ",unique(SI$original.data$community))),
         pch=19,
         col=1:length(unique(SI$original.data$community)))
  plotGroupEllipses(SI, n = 100, p.interval = 0.95,
                    ci.mean = T, lty = 1, lwd = 2)
  
  print(paste("DataFrame", i))
}

SIBER_spp<-as.data.frame(t(SIBER_group))
SIBER_spp$group <- substr(rownames(SIBER_spp), nchar(rownames(SIBER_spp)) - 3 + 1, nchar(rownames(SIBER_spp)))

write.csv(SIBER_comm,file = "CommunitySiberResults_11172025.csv")
write.csv(SIBER_spp,file = "GroupSiberResults11172025.csv")

ggplot(SIBER_spp, aes(x = SEAc)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black") +
  facet_wrap(~ group) +
  labs(title = "Histogram of Ellipse Area by Species", x = "Standard Ellipse Area corrected for sample size", y = "Count")

median_order <- SIBER_spp %>%
  group_by(group) %>%
  summarise(med_SEA = median(SEAc, na.rm = TRUE)) %>%
  arrange(med_SEA) %>%
  pull(group)

SIBER_spp<-SIBER_spp[!SIBER_spp$group=="RKB",]

SIBER_spp <- SIBER_spp %>%
  mutate(group = factor(group, levels = median_order))
ggplot(SIBER_spp, aes(y=SEAc, x=group))+
  geom_jitter(color = "gray60",
              width = 0.2,
              alpha = 0.5,
              size = 1.5)+
  geom_boxplot(
    aes(fill = group),
    alpha = 0.7)+
  labs(
    title = "Standard Ellipse Area across all lakes",
    x = "Species (ordered by median SEAc)",
    y = "Standard Ellipse Area (corrected for sample size)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_viridis_d(option = "viridis")
  

#loop through niche dfs to calculate overlap of cool and warm water species of interest
NicheArea<-read.csv("GroupSiberResults11172025.csv")#ellipse area for each spp lake combination 
NicheArea$community<-substr(NicheArea$X, 1, nchar(NicheArea$X) - 4)#getting col of lake years
#### final species list  ---
species_list<-c("BLC","BLG","CIS","LMB","NOP","SMB","WAE","YEP")

#function to calcualte overlapping region for all species combination that overlap in a given lake
run_overlap_wide <- function(SI, species_list, community, p.interval = 0.95, n = 100) {
  overlap_values <- list() #make empty list to store results
  
  for (i in 1:(length(species_list) - 1)) { #loop through by using unique species combinations
    for (j in (i + 1):length(species_list)) {
      sp1 <- species_list[i]
      sp2 <- species_list[j]
      
      sp1_full <- paste0(community, ".", sp1)
      sp2_full <- paste0(community, ".", sp2)
      
      result <- tryCatch({ #make it skip over null results and continue to loop through other combinations 
        maxLikOverlap(sp1_full, sp2_full, SI, p.interval = p.interval, n = n)
      }, error = function(e) {
        cat("Error in overlap for", sp1, "vs", sp2, ":", e$message, "\n")
        return(NULL) 
      })
      
      if (!is.null(result)) { #if it is not a null result capture results and store them to append to df
        # result is a named numeric vector
        overlap_value <- result["overlap"]
        pair_name <- paste(sp1, sp2, sep = "_vs_") #record the species and the order in when they are input into the function
        overlap_values[[pair_name]] <- overlap_value
        cat("Computed overlap for", sp1, "vs", sp2, ":", overlap_value, "\n")
      }
    }
  }
  
  # Convert to dataframe row (if anything was computed)
  if (length(overlap_values) == 0) {
    cat("No overlaps computed for", community, "\n")
    return(NULL)
  }
  
  df_row <- as.data.frame(overlap_values) #store overlap results for each combination 
  df_row$community <- community
  df_row <- df_row[, c("community", setdiff(names(df_row), "community"))]
  
  return(df_row)
}


final_results<-data.frame()#store overlap results
#looping through all lakes to calculate niche overlap metrics 
for(i in seq_along(df_list)){
  df <- df_list[[i]]#get individual lake
  
  # Check if dataframe is empty
  if (nrow(df) == 0) {
    cat("Skipping empty dataframe", i, "\n")
    next  # Skip to the next dataframe
  }
  
  SI<-createSiberObject(df)
  check<-run_overlap_wide(SI,species_list=unique(SI$original.data$group),community=unique(SI$original.data$community),p.interval = 0.95, n = 100)
  
  #append to master df 
  final_results<-bind_rows(final_results,check)
  
}
# the overlap betweeen the corresponding 95% prediction ellipses is given by:
write.csv(final_results,"NicheOverlapCalcs_092425.csv")

#### function for overlap and adding in 0s where the species don't overlap but coexist and NAs where they dont co-occur 
run_overlap_wide_O <- function(SI, species_list, community, p.interval = 0.95, n = 100) {
  # Make empty named vector to hold all pair results (including 0s)
  overlap_values <- setNames(rep(NA, choose(length(species_list), 2)), character(0))
  
  pair_index <- 1
  
  for (i in 1:(length(species_list) - 1)) {
    for (j in (i + 1):length(species_list)) {
      sp1 <- species_list[i]
      sp2 <- species_list[j]
      
      sp1_full <- paste0(community, ".", sp1)
      sp2_full <- paste0(community, ".", sp2)
      pair_name <- paste(sp1, sp2, sep = "_vs_")
      
      # Check if both species are present in this community
      sp_present <- sp1 %in% dimnames(SI$ML.mu[[community]])[[3]] &&
        sp2 %in% dimnames(SI$ML.mu[[community]])[[3]]
      
      if (sp_present) {
        result <- tryCatch({
          maxLikOverlap(sp1_full, sp2_full, SI, p.interval = p.interval, n = n)
        }, error = function(e) {
          cat("Error in overlap for", sp1, "vs", sp2, ":", e$message, "\n")
          return(NA_real_)  # set to NA instead of 0 if there's an error
        })
        
        if (is.numeric(result) && "overlap" %in% names(result)) {
          overlap_val <- result["overlap"]
          if (is.na(overlap_val)) {
            overlap_values[pair_name] <- NA
          } else if (overlap_val < 0.001) { #putting 0s for extremely small values 
            overlap_values[pair_name] <- 0
            cat("Overlap for", sp1, "vs", sp2, ": 0 (zero overlap)\n")
          } else {
            overlap_values[pair_name] <- overlap_val #if its greater than ~ 0 place value in column 
            cat("Overlap for", sp1, "vs", sp2, ":", overlap_val, "\n")
          }
        } else {
          # Result was NA due to an error in computation
          overlap_values[pair_name] <- NA #if model estimation does not converge 
        }
        
      } else {
        cat("Skipping", sp1, "vs", sp2, "— not both present in", community, "\n")
        overlap_values[pair_name] <- NA  # species not co-occurring place NA to differ from 0 
      }
      
      pair_index <- pair_index + 1
    }
  }
  
  # Always return a row with the community and full set of overlap values
  df_row <- as.data.frame(as.list(overlap_values))
  df_row$community <- community
  df_row <- df_row[, c("community", setdiff(names(df_row), "community"))]
  
  return(df_row)
}

final_results<-data.frame()#store overlap results
#looping through all lakes to calculate niche overlap metrics but placing 0s for pairs that co-occur but don't overlap 
for(i in seq_along(df_list)){
  df <- df_list[[i]]#get individual lake
  
  # Check if dataframe is empty
  if (nrow(df) == 0) {
    cat("Skipping empty dataframe", i, "\n") #some datasets don't have enough sample size to evaluate. loop will skip 
    next  # Skip to the next dataframe
  }
  
  SI<-createSiberObject(df)
  check<-run_overlap_wide_O(SI,species_list=unique(SI$original.data$group),community=unique(SI$original.data$community),p.interval = 0.95, n = 100)
  
  #append to master df 
  final_results<-bind_rows(final_results,check)
  
}
# the overlap betweeen the corresponding 95% prediction ellipses is given by:
final_results<-final_results[, grep("^(NA)", names(final_results), value = TRUE, invert = TRUE)]#removing NA columns created by model errors
write.csv(final_results,"NicheOverlapCalcs_11172025.csv")

final_results<-read.csv("NicheOverlapCalcs_11172025.csv")
final_results<-final_results[,-1]
#average overlap acorss communities 
# Remove df_name column and calculate column means
overlap_means <- colMeans(final_results[ , -1], na.rm = TRUE)
#overlap_means<-overlap_means[2:34]

# Convert to dataframe for easier display
average_summary <- data.frame(
  species_pair = names(overlap_means),
  average_overlap = overlap_means,
  row.names = NULL
)

# View summary
average_summary<-average_summary[!is.na(average_summary$average_overlap),]
print(average_summary)

#plot
ggplot(average_summary, aes(x=average_overlap,y=species_pair))+labs(x="Average overlap (niche area)",title="Average niche overlap area for species pairs")+geom_point()

# Reshape to long format
long_data <- final_results %>%
  pivot_longer(
    cols = -community,  # keep 'community' as is, pivot everything else
    names_to = "comparison",
    values_to = "value"
  )

# Plot heatmap
ggplot(long_data, aes(x = comparison, y = community, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(na.value = "grey90") +  # or use scale_fill_gradient
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Comparison",
    y = "Community",
    fill = "Value",
    title = "Heatmap of Niche Overlap by Lake and Species Comparisons"
  )
# Plot heatmap

ggplot(long_data, aes(x = comparison, y = community, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  theme_ipsum() +
  labs(
    x = "Comparison",
    y = "Community",
    fill = "Value",
    title = "Heatmap of Niche Overlap by Lake and Species Comparisons",
    )+theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### overlap switching species denominator in equation to calculate different overlap combos###
run_overlap_wide_SP2 <- function(SI, species_list, community, p.interval = 0.95, n = 100) {
  overlap_values <- list()       # Store overlap values
  proportions_list <- list()     # Store proportion metrics
  
  for (i in 1:(length(species_list) - 1)) {
    for (j in (i + 1):length(species_list)) {
      sp1 <- species_list[i]
      sp2 <- species_list[j]
      
      sp1_full <- paste0(community, ".", sp1)
      sp2_full <- paste0(community, ".", sp2)
      pair_name <- paste(sp2, sp1, sep = "_vs_")
      
      # Check for co-occurrence in lake
      sp_present <- sp1 %in% dimnames(SI$ML.mu[[community]])[[3]] &&
        sp2 %in% dimnames(SI$ML.mu[[community]])[[3]]
      
      if (!sp_present) {
        overlap_values[[pair_name]] <- NA
        cat("Skipping", sp2, "vs", sp1, "— not both present in", community, "\n")
        next
      }
      
      # Try to compute overlap
      result <- tryCatch({
        maxLikOverlap(sp2_full, sp1_full, SI, p.interval = p.interval, n = n)
      }, error = function(e) {
        cat("Error in overlap for", sp2, "vs", sp1, ":", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(result)) {
        raw_overlap <- as.numeric(result["overlap"])
        
        if (is.na(raw_overlap) || raw_overlap < 0.001) {
          overlap_values[[pair_name]] <- 0
          cat("Overlap for", sp2, "vs", sp1, ":", raw_overlap, " — recorded as 0\n")
        } else {
          overlap_values[[pair_name]] <- raw_overlap
          cat("Overlap for", sp2, "vs", sp1, ":", raw_overlap, "\n")
        }
        
        # Compute proportions each pair overlap ellipse given species are both present 
        prop.of.first  <- as.numeric(result["overlap"] / result["area.1"])
        prop.of.second <- as.numeric(result["overlap"] / result["area.2"])
        prop.of.both   <- as.numeric(result["overlap"] / (result["area.1"] + result["area.2"]))
        
        proportions_list[[pair_name]] <- data.frame( #making empty df to store output 
          community = community,
          pair = pair_name,
          prop.of.first = prop.of.first,
          prop.of.second = prop.of.second,
          prop.of.both = prop.of.both,
          stringsAsFactors = FALSE
        )
      } else {
        overlap_values[[pair_name]] <- 0  # optional: you can leave this out to only handle NULL with NA or 0
        cat("Overlap result NULL for", sp2, "vs", sp1, "— recorded as 0\n")
      }
    }
  }
  
  # Convert overlap standard ellipse area values to wide-format data.frame
  overlap_df <- as.data.frame(t(as.data.frame(overlap_values)))
  overlap_df$community <- community
  
  # Reorder columns
  overlap_df <- overlap_df[, c("community", setdiff(names(overlap_df), "community"))]
  
  # Count species present in this community
  species_present <- species_list[species_list %in% dimnames(SI$ML.mu[[community]])[[3]]]
  overlap_df$n_species_present <- length(species_present)
  
  # Combine proportions
  proportions_df <- if (length(proportions_list) > 0) {
    do.call(rbind, proportions_list)
  } else {
    data.frame()
  }
  
  return(list(
    overlap_wide = overlap_df,
    overlap_props = proportions_df
  ))
}


#running function to get overlap across communities 
# Create empty master dataframes to store results
final_overlap_wide <- data.frame()
final_overlap_props <- data.frame()

for (i in seq_along(df_list)) {
  df <- df_list[[i]]
  
  if (nrow(df) == 0) {
    cat("Skipping empty dataframe", i, "\n")
    next
  }
  
  SI <- createSiberObject(df)
  species_list <- unique(SI$original.data$group)
  community <- unique(SI$original.data$community)
  
  # Run overlap function
  check <- run_overlap_wide_SP2(SI, species_list, community, p.interval = 0.95, n = 100)
  
  if (!is.null(check)) {
    # Append both outputs separately
    final_overlap_wide <- bind_rows(final_overlap_wide, check$overlap_wide)
    final_overlap_props <- bind_rows(final_overlap_props, check$overlap_props)
  }
}

write.csv(final_overlap_wide, "overlap_wide_results_11172025.csv", row.names = FALSE)
write.csv(final_overlap_props, "overlap_proportions_long_11172025.csv", row.names = FALSE)

#read in data now 
final_overlap_props<-read.csv("/Users/cammosley/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/overlap_proportions_long_092425.csv")
final_overlap_wide<-read.csv("/Users/cammosley/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/overlap_wide_results_092425.csv")

#proportionality in overlap plots 
a<-ggplot(data = final_overlap_props, aes(prop.of.first,prop.of.second,colour = pair))+geom_point()
a+geom_abline()+facet_wrap(~pair)#adding line to look at which species are dominating niche space in pairs across lakes 

final_overlap_props  <- final_overlap_props %>% 
  mutate(across(where(is.numeric), ~round(., 3))) #round the values to get 0s for extremely small values 
  
#get columns for each species in the pair and attach the lake data 
final_overlap_props$sp1<-substr(final_overlap_props$pair,1,3)
final_overlap_props$sp2<-substr(final_overlap_props$pair,8,10)
final_overlap_props<-final_overlap_props[!final_overlap_props$sp1=="RKB",]
final_overlap_props<-final_overlap_props[!final_overlap_props$sp2=="RKB",]

#read in lake data 
DATA<-read.csv("/Users/cammosley/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/TP_Lake_LR_Data_04282025.csv")
final_overlap_props$Community<-final_overlap_props$community#change column name to match for joining 
test<-full_join(final_overlap_props,DATA,by="Community")

#test<-full_join(test,final_overlap_wide,by="community")
test$nhdhr<-substr(test$nhdhr_id,7,15) 
#add in invasive species 
#read in MN dnr lakes list and fix column names 
infested_waters_excel <- read_excel("~/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/lake covariate data/infested-waters.xlsx", sheet = 1, cell_rows(2:1326)) # read it in as an excel file
# The March 21st 2023 is the most up to date as of April 26th 23

#clean the data
IW <- infested_waters_excel %>% 
  rename(dow = 6) %>% # these renames basically make it easier to work with the columns
  rename(Water_Body_Name = 1) %>%
  rename (County = 2) %>%
  rename(AIS_Species = 3) %>%
  rename(Year_Infested = 4)%>%
  rename (Year_Confirmed = 5)%>%
  select(dow,Water_Body_Name,County,AIS_Species,Year_Infested,Year_Confirmed)%>% # this just selected the columns that I wanted, you can select as many or few as you want
  filter(!dow == "none")%>% #removes rows where DOW is listed as "none"
  filter(!dow == "none, part of Winnibigoshish") # not sure why this isn't listed as Winnibigoshish...

#### Must reformat DOW (currently 6 digits with dash, we want it to be 8 digits no dashes)
IW_fixed <- IW %>% mutate(dow = gsub("-","",IW$dow)) %>% # removes all dashes from dow IW_nodash$dow
  mutate(DOW = str_pad(dow,8, side="right", pad= "0")) %>% # This adds zeroes so that all DOWs are 8 digits
  rename(parent_dow=dow) # sometimes DNR calls the first 6 digits of the dow the parent dow, this ignores the last two digits that identify sub-basins

IW_parent_dow <- IW_fixed %>% # ensures everything is the same length with no additional spaces
  mutate(PARENT_DOW2 = str_trunc(IW_fixed$DOW, ellipsis= "", side = "right",6))%>%
  mutate(PARENT_DOW = str_trim(PARENT_DOW2, side = "right"))%>%
  select(-PARENT_DOW2, - parent_dow, -Year_Confirmed)

# This is my code to just get zebra mussel data
ZM <- IW_parent_dow %>% 
  filter(AIS_Species== "zebra mussel") %>%
  filter(!is.na(DOW))%>% #filters out rivers and creeks (anything without a DOW)
  select(AIS_Species, Year_Infested, DOW)

#getting lake ID fron nldas dataset 
MN_lake_all_temp_metrics_all_years<-read.csv("/Users/cammosley/Google Drive/Shared drives/Hansen Lab/Data resources/Lake temp metrics/MN_lake_all_temp_metrics_all_years.csv")
Check <- MN_lake_all_temp_metrics_all_years[MN_lake_all_temp_metrics_all_years$site_id%in%c(test$nhdhr_id),]
Check<-Check[Check$year > 2020,]
Check<-distinct(Check)

ZM<-ZM[ZM$DOW%in%c(Check$DOW),]#ZM lakes that match study lakes dows
ZM$DOW<-as.integer(ZM$DOW)
ZM_check<-full_join(ZM,Check)
ZM_check<-full_join(ZM_check,test)
ZM_check<-ZM_check[!is.na(ZM_check$pair),]
test<-ZM_check
invaded<-c("Alexander-2021","Big Cormorant-2019","Chippewa-2021","Clearwater-2019","Green-2020","Gull-2019",
           "Leech-2022","Melissa-2019","North Lida-2021","Pelican-2020","Rush-2022",
           "Steamboat-2021","Tenmile-2020")
test$ZM<-ifelse(test$community %in%c(invaded), 1, 0)
test$ZM<-as.factor(test$ZM)

#coldwater fish present 
CF<-read.csv("/Users/cammosley/Google Drive/Shared drives/Hansen Lab/RESEARCH PROJECTS/Fish Survey Data/MN_Data/Data for relative abundance modeling/coldwater_fish_PA.csv")
CF<-CF[CF$DOWNumber%in%c(Check$DOW),]
coldfish<-c("Alexander-2021","Chippewa-2021","Clearwater-2019","Gull-2019",
           "Leech-2022","Green-2020" ,"North Lida-2021","Pelican-2020","Rush-2022",
           "Steamboat-2021","Plantagenet-2022","Big Sand-2020","Potato-2022","Koronis-2019")
test$cold<-ifelse(test$community %in%c(coldfish), 1, 0)
test$cold<-as.factor(test$cold)

#thermal habitat data 
#VHOD<-read.csv("/Users/cammosley/Downloads/VHOD.csv")
#test<-left_join(final_overlap_props,test)
test<-test[!is.na(test$pair),]
test<-test[ ,c(1:187,189:217)]
test<-distinct(test)

test_unique <- test %>% distinct(community, pair, .keep_all = TRUE)

merged <- final_overlap_props %>%
  left_join(test_unique, by = c("community", "pair"))

#merged<-merged[!is.na(merged$prop.of.first.y),]

# Standardization function
standardize_name <- function(name) {
  name <- tolower(name)
  name <- gsub("\\.", "-", name)  # Replace dots with hyphens
  name <- gsub(" ", "-", name)    # Replace spaces with hyphens
  return(name)
}

# Find common communities
common_communities <- intersect(
  standardize_name(NicheArea$community),
  standardize_name(merged$community)
)

# Filter merged to only include matching communities
merged_filtered <- merged[standardize_name(merged$community) %in% common_communities, ]

test<-merged[merged$community%in%c(NicheArea$community),]


write.csv(merged_filtered,"RF_Data_092625.csv")
Merge_Data<-read.csv("RF_Data_092425.csv")
#look at bass and walleye interactions 
#proporition of bass in walleye niches 
WAE_test<- test %>% filter(sp1.x == "LMB" | sp2.x == "WAE")
a<-lm(prop.of.first.x ~ lakearea.z, data = WAE_test)
summary(a)

b<-ggplot(data=WAE_test,aes(lakearea.z,prop.of.first))+geom_point()+geom_smooth()
c<-ggplot(data=WAE_test,aes(mean_gdd_0c.z,prop.of.first))+geom_point()+geom_smooth()
d<-ggplot(data=WAE_test,aes(IceOff,prop.of.first))+geom_point()+geom_smooth()
e<-ggplot(data=WAE_test,aes(clarity.z,prop.of.first))+geom_point()+geom_smooth()
f<-ggplot(data=WAE_test,aes(lakeperimeter.z,prop.of.first))+geom_point()+geom_smooth()
g<-ggplot(data=WAE_test,aes(secchi.z,prop.of.first))+geom_point()+geom_smooth()
h<-ggplot(data=WAE_test,aes(lakeshorelinefactor.z,prop.of.first))+geom_point()+geom_smooth()
i<-ggplot(data=WAE_test,aes(elevation.z,prop.of.first))+geom_point()+geom_smooth()
j<-ggplot(data=WAE_test,aes(depth.z,prop.of.first))+geom_point()+geom_smooth()
k<-ggplot(data=WAE_test,aes(total.dev.z,prop.of.first))+geom_point()+geom_smooth()
l<-ggplot(data=WAE_test,aes(total.ag.z,prop.of.first))+geom_point()+geom_smooth()
m<-ggplot(data=WAE_test,aes(total.for.z,prop.of.first))+geom_point()+geom_smooth()
n<-ggplot(data=WAE_test,aes(mean_surf.z,prop.of.first))+geom_point()+geom_smooth()
o<-ggplot(data=WAE_test,aes(x=n_species_present))+geom_bar(stat="count")
p<-ggplot(WAE_test, aes(x = ZM, y = prop.of.first))+geom_boxplot()
q<-ggplot(WAE_test, aes(x = cold, y = prop.of.first))+geom_boxplot()
r<-ggplot(data=WAE_test,aes(,prop.of.first))+geom_point()+geom_smooth()

#look at averages across communities using tidy verse
avg_overlap <- final_overlap_props %>%
  group_by(pair) %>%
  summarise(
    avg_prop_first = mean(prop.of.first, na.rm = TRUE),
    avg_prop_second = mean(prop.of.second, na.rm = TRUE),
    avg_prop_both = mean(prop.of.both, na.rm = TRUE),
    n = n()
  ) %>%
  arrange(desc(avg_prop_both))

#reformatting data results to look at overall overlap with walleye
heatmap_long <- final_overlap_props %>%
  select(community, pair, prop.of.both) 

# Heatmap with ggplot2 
#the proportion of species 1 and 2 that overlap with each other
ggplot(heatmap_long, aes(x = pair, y = community, fill = prop.of.both)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(na.value = "grey90", name = "Overlap\n(prop. of both)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Ellipse Overlap",
    x = "Species Pair",
    y = "Community"
  )

# Extract species names from pair
avg_matrix_data <- final_overlap_props %>%
  separate(pair, into = c("sp1", "sp2"), sep = "_vs_") %>%
  group_by(sp1, sp2) %>%
  summarise(mean_overlap = mean(prop.of.both, na.rm = TRUE), .groups = "drop")

# Make symmetrical by duplicating
avg_matrix_sym <- bind_rows(
  avg_matrix_data,
  avg_matrix_data %>% rename(sp1 = sp2, sp2 = sp1)
)

# Convert to matrix format
avg_heatmap_data <- avg_matrix_sym %>%
  pivot_wider(names_from = sp2, values_from = mean_overlap) %>%
  column_to_rownames("sp1") %>%
  as.matrix()

avg_heatmap_data<-avg_matrix_data

# Save plot as high-resolution PNG
png("average_overlap_heatmap.png", width = 2000, height = 1800, res = 300)
pheatmap(avg_heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         number_format = "%.2f",
         color = viridis::viridis(100),
         fontsize = 10,
         main = "Average Overlap Across Communities (Prop. of Both)")
dev.off()
library(tidyverse)

# 1. Reshape long data to wide matrix
plot_data <- avg_matrix_data %>%
  pivot_wider(names_from = sp2, values_from = mean_overlap) %>%
  column_to_rownames("sp1") %>%
  as.matrix()

# 2. Force it to be a full square matrix (filling missing species pairs with 0 or NA)
all_species <- sort(unique(c(avg_matrix_data$sp1, avg_matrix_data$sp2)))
final_matrix <- matrix(NA, nrow=length(all_species), ncol=length(all_species), 
                       dimnames=list(all_species, all_species))
final_matrix[as.matrix(avg_matrix_data[1:2])] <- avg_matrix_data$mean_overlap

# 3. Create the half-heatmap (Lower Triangle)
# Masking the upper triangle
final_matrix[upper.tri(final_matrix)] <- NA

# 4. Plot using ggplot2
df_melted <- as.data.frame(final_matrix) %>%
  rownames_to_column("sp1") %>%
  pivot_longer(-sp1, names_to = "sp2", values_to = "overlap") %>%
  filter(!is.na(overlap)) # Remove the masked upper triangle

library(viridis)

ggplot(df_melted, aes(sp1, sp2, fill = overlap)) +
  scale_fill_viridis(option    = "magma",
                     limits    = c(0, 0.5),      # ── scale to actual data range
                     breaks    = c(0, 0.05, 0.10, 0.15, 0.20),  # clean legend ticks
                     labels    = c("0", "0.05", "0.10", "0.15", "0.20"),
                     name      = "Proportion\nof Overlap") +
  geom_tile() +
  geom_text(aes(label  = round(overlap, 3),
                colour = overlap > 0.11),         # midpoint of actual range
            size = 5, show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE"  = "black",
                                 "FALSE" = "white")) +
  theme_classic() +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    axis.text    = element_text(size = 10)
  ) +xlab("")+ylab("")

ggplot(df_melted, aes(sp1, sp2, fill = overlap)) +
  
  # ── swap to fish species palette ─────────────────────────────────────────
  scale_fill_gradientn(
    colours = c("#313695","#4575b4","#74add1","#e6f598",
                "#fee090","#f46d43","#d73027","#a50026"),
    limits  = c(0, 0.3),
    breaks  = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
    labels  = c("0", "0.10", "0.20", "0.30", "0.40", "0.50"),
    name    = "Proportion\nof Overlap"
  ) +
  
  geom_tile() +
  
  geom_text(aes(label  = round(overlap, 3),
                colour = overlap > 0.25),         # midpoint of 0-0.5 range
            size = 7, show.legend = FALSE) +
  
  scale_colour_manual(values = c( "black")) +
  theme_classic() +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    axis.text    = element_text(size = 13)
  ) +
  xlab("") + ylab("")


# Assuming 'final_matrix' is a symmetric matrix with NA in the upper triangle

df_melted <- as.data.frame(final_matrix) %>%
  rownames_to_column("sp1") %>%
  pivot_longer(-sp1, names_to = "sp2", values_to = "overlap") %>%
  filter(!is.na(overlap)) %>%
  # Ensure both sp1 and sp2 are factors with the same order
  mutate(sp1 = factor(sp1, levels = unique(sp1)),
         sp2 = factor(sp2, levels = unique(sp1)))

ggplot(df_melted, aes(sp1, sp2, fill = overlap)) +
  geom_tile(color = "white") +
  scale_fill_gradient2() +
  geom_text(aes(label = round(overlap, 3)), size = 5) +
  # --- Key modification: reverse the y-axis ---
  scale_y_reverse() + 
  theme_classic2() +
  labs(title = "Average Niche Overlap", fill = "Proportion of overlap") +
  xlab("") + ylab("")


#Community niche area by lake plots 
library(data.table)

get_hull_df <- function(data, x_var, y_var, group_var) {
  
  dt <- as.data.table(data)
  
  # Compute hull for each group
  dt_hull <- dt[, .SD[chull(get(x_var), get(y_var))], by = c(group_var)]
  
  # Ensure only one group column exists
  names(dt_hull)[names(dt_hull) == group_var] <- "community"
  
  return(as.data.frame(dt_hull))
}
hulls <- get_hull_df(IsoData, "iso1", "iso2", "community")
#getting year off lake names 
hulls$community<-substr(hulls$community, 1, nchar(hulls$community) -5)
hull_plot <- ggplot() +
  geom_polygon(
    data = hulls,
    aes(x = iso1, y = iso2, fill = community, group = community),
    alpha = 0.3,
    color = "black"
  ) +
  labs(
    x = bquote(Carbon ~ delta^{13} * C),
    y = bquote(Nitrogen ~ delta^{15} * N)
  ) +
  scale_x_continuous(
    limits = range(IsoData$iso1),
    breaks = pretty(IsoData$iso1, n = 5)
  ) +
  scale_y_continuous(
    limits = range(IsoData$iso2),
    breaks = pretty(IsoData$iso2, n = 5)
  ) +
  theme_minimal()+
  theme(legend.position = "none")

a<-hull_plot + facet_wrap(~community)

hull_all <- ggplot() +
  geom_polygon(
    data = hulls,
    aes(x = iso1, y = iso2, fill = community, group = community),
    alpha = 0.3,
    color = "black"
  ) +
  #geom_point(
    #data = IsoData,
    #aes(x = iso1, y = iso2, color = community),
   # alpha = 0.6
 # ) +
  labs(
    title = "Convex Hulls Areas of Fish Communities by Lake",
    x = bquote(Carbon ~ delta^{13} * C),
    y = bquote(Nitrogen ~ delta^{15} * N)
  ) +
  theme_minimal()+
  theme(legend.position = "none")
#join the plots 
ggarrange(hull_all,a)


hist(SIBER_comm$MNND)#mean nearest neighbour distance of the means
hist(SIBER_comm$SDNND)#the standard deviation of the nearest neighbour distance

ggplot(SIBER_comm, aes(y=MNND,x=community))+geom_point()

# Add original data points
#geom_point(data = iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) 
#geom text option for area value labels label = paste("Area:", round(hull_areas$area, 2)

hist(hull_areas$area)
abline(v=mean(hull_areas$area),col="red",lwd=3)
lines(density(hull_areas$area),col="green",lwd=3)

# Add a legend
legend("topright", legend = c("Kernel Density", "Normal Fit"), col = c("darkgreen", "darkblue"), lwd = 2, lty = c(1, 2))

#look at grouped ellipse area by species 
hull_areas_S <- get_hull_area(IsoData, "iso1", "iso2", "group")
ggplot() +
  geom_sf(data = hull_areas_S, aes(fill = group), alpha = 0.3, inherit.aes = FALSE) +
  geom_text(data = as.data.frame(st_coordinates(st_centroid(hull_areas_S))), 
            aes(x = X, y = Y, label = " "), 
            color = "black", size = 4) +
  labs(title = "Species Convex Hull Areas of Stable Isotope Data",
       x = bquote( Carbon ~delta^{13} * C) ,
       y = bquote( Nitrogen ~delta^{15} * N)) +
  theme_minimal()

#stats test of overalp proportions 

# removing NAs to make data calculations without errors
df_test <- final_overlap_props %>%
  filter(!is.na(prop.of.both)) %>%
  separate(pair, into = c("sp1", "sp2"), sep = "_vs_")

library(lme4)
library(lmerTest)  # gives p-values

# Random effect for community, fixed effects for species pairs, check if certain species have higher overlap and random lake effect 
model <- lmer(prop.of.both ~ sp1 + sp2 + (1 | community), data = df_test)
summary(model)

#look at lake differences using lake as fixed effect this time for diff species pairs
model_fixed <- lm(prop.of.both ~ sp1 + sp2 + community, data = df_test)
anova(model_fixed)

#test if overlap is predicted by species specific comboinations 
model_inter <- lm(prop.of.both ~ sp1 * sp2 + community, data = df_test)
anova(model_inter)

library(emmeans)
emmeans(model_inter, pairwise ~ sp1:sp2)#calcuate stats for each species pairs 
ggplot(df_test, aes(x = sp1, y = prop.of.both, fill = sp2)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#box plot of overlap 
ggplot(final_overlap_props, aes(x = pair, y = prop.of.both)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Overlap Differences by Species Pair", y = "Proportion of both species ellipse overlap", x = "Species Pair")

#data test using overlap sea c values
final_overlap_wide$pair=substr(rownames(final_overlap_wide),1,10)#pair column
final_overlap_wide$sp1=substr(final_overlap_wide$pair,1,3)
final_overlap_wide$sp2=substr(final_overlap_wide$pair,8,11)

#remove NAs and make species factor variables instead of characters 
cols<-c("sp1","sp2")
df_test_SEA <- final_overlap_wide %>%
  filter(!is.na(V1)) %>%
  mutate_at(cols,factor)

# Random effect for community, fixed effects for species pairs, check if certain species have higher overlap and random lake effect 
model <- lmer(V1 ~ sp1 + sp2 + (1 | community), data = df_test_SEA)
summary(model)

#look at lake differences using lake as fixed effect this time for diff species pairs
model_fixed <- lm(V1 ~ sp1 + sp2 + community, data = df_test_SEA)
anova(model_fixed)

#test if overlap is predicted by species specific comboinations 
model_inter <- lm(V1 ~ sp1 * sp2 + community, data = df_test_SEA)
anova(model_inter)
emmeans(model_inter, ~ sp1:sp2) #post hoc comparison to look at combinations 

#more plotting 
em_results <- as.data.frame(emmeans(model_inter, ~ sp1:sp2))
em_clean <- em_results %>%
  filter(!is.na(emmean) & lower.CL > 0) %>%
  mutate(
    pair = paste(sp2, "vs", sp1, sep = "_"),  # label for x-axis
    pair = factor(pair, levels = unique(pair))  # optional: preserve original order
  )
ggplot(em_clean, aes(x = pair, y = emmean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(
    x = "Species Pair",
    y = "Estimated Mean Overlap (SEAc)",
    title = "Pairwise Niche Overlap with 95% confidence intervals"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

#look at walleye 
wae_pairs <- em_clean %>%
  filter(sp1 == "WAE" | sp2 == "WAE")

ggplot(wae_pairs, aes(x = pair, y = emmean, fill=sp2)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(
    x = "Species Pair",
    y = "Estimated Mean Overlap (SEAc)",
    title = "Pairwise Niche Overlap with 95% confidence intervals"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

#### linear model comparison SEA ####
# Random effect for community, fixed effects for species , check if certain species have higher standard ellipse area/niche breadth and random lake effect 
NicheArea_test<-NicheArea[!NicheArea$group=="RKB",]
model <- aov( SEAc ~ community, data = NicheArea_test) 
model2<-aov( SEAc ~ community + group, data = NicheArea_test)
anova(model,model2)

ggplot(data = NicheArea_test, aes(x=SEAc, fill= group))+geom_histogram(binwidth = 1)+ theme_bw()+labs(title = "Frequency of Standard Ellipse Area across Fish Species")

NicheAreaMeans <- NicheArea_test %>% 
  group_by(group) %>%
  summarise(SEAc=mean(SEAc))

ggplot(NicheAreaMeans, aes(y=SEAc,x=group, colour = group))+geom_point(size=3)+theme_bw()+labs(title = "Mean Standard Ellipse Area", x="Species",y="Standard Ellipse Area (sample size corrected)")

#niche area summary

NicheArea$lake<-substr(NicheArea$community, 1, nchar(NicheArea$community)-5)
TrophicP_LittorialR_median111420205$Lake<-substr(TrophicP_LittorialR_median111420205$group, 1, nchar(TrophicP_LittorialR_median111420205$group)-5)
TP<-TrophicP_LittorialR_median111420205[,c(3:7)]

#raw iso samples for species in all lakes 
FishIsoCHECK<-FishIso[FishIso$Spp%in%c(TPFish$consumer),]
Fish_samps_table<-table(FishIsoCHECK$Spp,FishIsoCHECK$Community)
Fish_samps_margins <-addmargins(Fish_samps_table)
print(Fish_samps_margins)
Fish_samps_df<-as.data.frame(Fish_samps_margins)
write.csv(Fish_samps_margins,file="Fish_iso_samps_sum.csv")

TEST$lake<-substr(TEST$Community, 1, nchar(TEST$Community)-5)

#lake size check 
MNLakeInfo<-read.csv("~/Google Drive/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/lake covariate data/Copy of mn_lake_list.csv")


