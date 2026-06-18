#clear environment and plots 
remove(list = ls())
graphics.off()

#loading in packages 
library(randomForest)
library(earth)
library(doParallel)
library(caret)
library(vip)
library(pdp)
library(iml)
library(dplyr)
library(ceterisParibus)
library(DALEX)
library(ggiraph)
library(hstats)
library(patchwork)
library(ggplot2)
library(grid)
library(reprtree)
library(lme4)
library(Hmisc)
#library(brms)

####### USE REAL DATA ######
#working directory in Gdrive 
setwd("/Users/cammosley/Google Drive/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/") #iMac
setwd("/Users/cammosley/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/TP models/")#macbook

#load in niche data 
#Dat<-read.csv("RF_Data.csv")
Dat<-read.csv("RF_Data_092625.csv") #correct data for subset of lakes with TP and LR

#add in physiological metrics 
setwd("/Users/cammosley/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/")#macbook
Physi<-read.csv("fish_tempCorrected.csv")
Physi<-Physi[,c(4:7,9,13,14,18)] # get columns of interest 
Physi$sp1.x=Physi$species_code# make columns for sp1 and sp2 to join physio data for each niche pair 
Dat<-full_join(Dat,Physi,by="sp1.x")

#ADD in physio data for species 2 in the pair 
Physi$sp2.x=Physi$species_code
Dat<-full_join(Dat,Physi,by="sp2.x")

#make a metric that takes the difference of CT max values,OGT,and FTP
Dat$CTdiff<-abs(Dat$Ctmax.x - Dat$Ctmax.y)#Difference in CT max of overlapping individuals 
Dat$OGTdiff<-abs(Dat$OGT.x - Dat$OGT.y)#Difference in CT max of overlapping individuals 
Dat$FTPdiff<-abs(Dat$FTP.x - Dat$FTP.y)#Difference in CT max of overlapping individuals 
Dat$thermalguild<-paste(Dat$Temperature.preference.class.x,Dat$Temperature.preference.class.y) #make categorical thermal guild variable
Dat<-Dat[!is.na(Dat$Temperature.preference.class.x),]
Dat<-Dat[!is.na(Dat$Temperature.preference.class.y),]

#standardize thermal guild categorical variable
# Function to sort words within each value
sort_words <- function(x) {
  sapply(strsplit(as.character(x), " "), function(words) {
    paste(sort(words), collapse = " ")
  })
}

# Apply to the thermalguild column
Dat$thermalguild <- sort_words(Dat$thermalguild)

# Check the results and make a factor variable 
unique(Dat$thermalguild)
Dat$thermalguild<-as.factor(Dat$thermalguild)

#Walleye sort 
#Dat <- Dat %>%
  #filter(
   # grepl("WAE", pair) & 
   #   (grepl("LMB", pair) | grepl("SMB",pair))
 # )

#making sure to remove NAs from dataframe to predict certain metric
Dat<-Dat[!is.na(Dat$community),]
Dat$community<-as.factor(Dat$community)

#clean df
write.csv(Dat,"RF_model_clean_df_11182025.csv")

Dat<-read.csv("RF_model_clean_df_11182025.csv") #requires data formatting if you read in csv file

#remove rkb values 
Dat<-Dat[!Dat$sp1.x.x=="RKB",]
Dat<-Dat[!Dat$sp2.x=="RKB",]

#get columns of interest
#metric testing prop.of.both 
Dat <- Dat %>%
  dplyr::select(community,
         prop.of.both.x,
         pair,
         lakeperimeter.z,
         mean_gdd_0c.z,
         secchi.z,
         IceOff,
         post_ice_warm_rate,
         winter_dur_0_4,
         mean_epi_hypo_ratio,
         lakeshorelinefactor.z,
         elevation.z,
         depth.z,
         max_surf_jul,
         mean_surf.z,
         lakearea.z,
         total.for.z,
         total.dev.z,
         total.dev.z,
         stratification_duration,
         height_19.3_23.3,
         height_27_32,
         height_10.6_11.2,
         days_26_28,
         ZM,
         OGTdiff,
         CTdiff,
         FTPdiff,
         cold,
         thermalguild)

DatNA<-Dat[is.na(Dat$mean_epi_hypo_ratio),]
Dat<-distinct(Dat)
#corecing factor variables to integers for impute function and rf model 

#Dat$thermalguild<-unclass(factor(Dat$thermalguild))
#Dat$community<-unclass(factor(Dat$community))
#Dat$sp1.x.x<-unclass(factor(Dat$sp1.x.x))
#Dat$sp2.x<-unclass(factor(Dat$sp2.x))

Dat$thermalguild<-factor(Dat$thermalguild)
Dat$community<-factor(Dat$community)
Dat$pair<-factor(Dat$pair)
#Dat$sp1.x.x<-factor(Dat$sp1.x.x)
#Dat$sp2.x<-factor(Dat$sp2.x)

# Impute missing values using rfImpute()
# Note: rfImpute requires a response variable without missing values
Dat <- rfImpute(prop.of.both.x ~ ., data = Dat, iter = 5)

set.seed(222)
ind <- sample(2, nrow(Dat), replace = TRUE, prob = c(0.7, 0.3))
train <- Dat[c(1:314),] #218 observations
test <- Dat[315,] #97 observations

rf <- randomForest(prop.of.both.x~., data=train, ntree=500, proximity=TRUE, importance=TRUE,) 
print(rf)

#Call:
#Call:
  #randomForest(formula = prop.of.both.x ~ ., data = train, ntree = 500,      proximity = TRUE) 
#Type of random forest: regression
#Number of trees: 500
#No. of variables tried at each split: 9

#Mean of squared residuals: 0.00635712
#% Var explained: 33.44

#prediction and confusion matrix 
p1 <- predict(rf, train)

#? confusionMatrix(p1,train$prop.of.both.x) #Error: `data` and `reference` should be factors with the same levels.

plot(rf) #looking at error rate, 500 trees bring error to about zero
#plotting number of tree nodes
hist(treesize(rf),
     main = "No. of Nodes for the Trees",
     col = "green")

#Variable Importance
varImpPlot(rf,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance for Random Forest Model (niche overlap)")
importance_matrix<-importance(rf,scale = TRUE)
importance_df <- as.data.frame(importance_matrix)
importance_df$Feature<-rownames(importance_df)
importance_matrix <- data.frame(
  Rank = 1:nrow(importance_df),
  Feature = rownames(importance_df),
  IncNodePurity = importance_df$IncNodePurity,
  IncNodePurity_pct = round(importance_df$IncNodePurity * 100, 2)
)

var_IMP<-importance(rf)
print(var_IMP)
var_IMP<-as.data.frame(var_IMP)

#sort by highest value
IMP_sorted <- var_IMP[order(-var_IMP$`%IncMSE`),]
print(IMP_sorted)
IMP_sorted$feature<-rownames(IMP_sorted)
IMP_sorted$Rank<-1:28

#visualize rf model decison tee 
#plotting represetative model tree 
reprtree::plot.getTree(rf)

### fix error code, 20 fold cross validation 


# Get top 10 variables for %IncMSE
top10_mse <- importance_df %>%
  arrange(desc(`%IncMSE`)) %>%
  head(10)

# Get top 10 for IncNodePurity (no error bars typically shown)
top10_purity <- importance_df %>%
  arrange(desc(IncNodePurity)) %>%
  head(10)

# Plot 1: %IncMSE with error bars
p1 <- ggplot(top10_mse, aes(x = `%IncMSE`, y = reorder(Feature, `%IncMSE`))) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = `%IncMSE` - IncMSE_SE, 
                     xmax = `%IncMSE` + IncMSE_SE),
                 height = 0.2, alpha = 0.6) +
  labs(x = "%IncMSE", y = "", 
       title = "Top 10 - Variable Importance") +
  theme_bw() +
  theme(panel.grid.major.y = element_line(linetype = "dotted", color = "gray80"),
        panel.grid.major.x = element_blank())


# Sort by IncNodePurity
importance_matrix <- importance_matrix[order(-importance_matrix$IncNodePurity), ]

# Reassign ranks after sorting
importance_matrix$Rank <- 1:nrow(importance_matrix)
print(importance_matrix)

# Save as CSV 
write.csv(IMP_sorted, "feature_importance_ranked_12032025.csv", row.names = FALSE)

#error 
a<-ggplot(data.frame(Predicted = p1[1:97], Actual = test$prop.of.both.x), 
          aes(x = Actual, y = Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Random Forest Predictions on Niche Overlap Dataset prop of both", 
       x = "Actual Values", y = "Predicted Values")
a

#plot model features together
X=Dat[which(names(Dat) !="prop.of.both.x")]
model=Predictor$new(rf,data=X,y=Dat$prop.of.both.x)
effect=FeatureEffects$new(model)
# Create the feature effects plot - keep original x-axis labels
A <- effect$plot(features = c("OGTdiff", "FTPdiff", "CTdiff")) + #top 3 variables 
  plot_annotation(
    title = "Random forest model feature effects on niche overlap",
    subtitle = "Top 3 numerical model features",
    caption=" " #add in space for predictor universal x axis lab
  )

# Print the plot
print(A)

# Add universal x-axis label below everything using grid
grid.text("Predictor Variables", 
          x = 0.5, 
          y = 0.03,  # Adjust this value to position it below the individual labels
          gp = gpar(fontsize = 10, fontface = "bold"))
# Add universal y-axis label
grid.text("Predicted Proportion of Niche Overlap", 
          x = 0.007,  # Position on left side
          y = 0.5,   # Centered vertically
          rot = 90,  # Rotate 90 degrees
          gp = gpar(fontsize = 10, fontface = "bold"))

effect$plot(features=c("mean_epi_hypo_ratio","lakearea.z","stratification_duration"))
# Create the feature effects plot - keep original x-axis labels
B <- effect$plot(features=c("mean_epi_hypo_ratio","lakearea.z","stratification_duration")) + #top 6 (next top 3)
  plot_annotation(
    title = "Random forest model feature effects on niche overlap",
    subtitle = "Numerical model features",
    caption=" " #add in space for predictor universal x axis lab
  )

# Print the plot
print(B)

# Add universal x-axis label below everything using grid
grid.text("Predictor Variables", 
          x = 0.5, 
          y = 0.03,  # Adjust this value to position it below the individual labels
          gp = gpar(fontsize = 10, fontface = "bold"))
# Add universal y-axis label
grid.text("Predicted Proportion of Niche Overlap", 
          x = 0.007,  # Position on left side
          y = 0.5,   # Centered vertically
          rot = 90,  # Rotate 90 degrees
          gp = gpar(fontsize = 10, fontface = "bold"))
### Ceteris Paribus (Individual Conditional Expectation) Plots by Community
# Ceteris Paribus (ICE) Plots with Different Curves for Each Community
  
  # Now create function for all features
  create_ice_plot <- function(feature_name) {
    tryCatch({
      print(paste("\nCreating plot for:", feature_name))
      
      ice_effect <- FeatureEffect$new(
        model,
        feature = feature_name,
        method = "ice",
        grid.size = 30
      )
      
      ice_data <- ice_effect$results
      ice_data$row_index <- as.numeric(as.character(ice_data$.id))
      ice_data$community <- Dat$community[ice_data$row_index]
      
      # Find feature column
      feature_col <- names(ice_data)[!names(ice_data) %in% 
                                       c(".type", ".id", ".value", "community", "row_index")]
      
      if(length(feature_col) > 0) {
        # Use the actual feature column name
        print(paste("Using feature column:", feature_col[1]))
        
        # Add feature_value first
        ice_data$feature_value <- ice_data[[feature_col[1]]]
        
        # Calculate PDP more simply
        pdp_data <- aggregate(.value ~ feature_value, data = ice_data, FUN = mean)
        colnames(pdp_data)[colnames(pdp_data) == ".value"] <- "pdp_value"
        
        # Create plot
        p <- ggplot() +
          geom_line(data = ice_data, 
                    aes(x = feature_value, y = .value, color = community, group = .id), 
                    alpha = 0.4, linewidth = 0.5) +
          geom_line(data = pdp_data, 
                    aes(x = feature_value, y = pdp_value),
                    color = "black", linewidth = 2, 
                    linetype = "dashed") +
          labs(
            title = paste("ICE Plot:", feature_name),
            subtitle = "Individual curves by community (black dashed = average PDP)",
            x = feature_name,
            y = "Predicted niche overlap (proportion of ellipse area)",
            color = "Community"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 10),
            legend.position = "right",
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 9, face = "bold"),
            legend.key.height = unit(0.3, "cm")
          )
        
        return(list(plot = p, data = ice_data, pdp = pdp_data))
      } else {
        return(NULL)
      }
      
    }, error = function(e) {
      print(paste("Error:", e$message))
      return(NULL)
    })
  }
  
  # Create plots for all features
  all_features <- c("CTdiff","FTPdiff", "OGTdiff", "mean_epi_hypo_ratio", 
                    "lakearea.z","stratification_duration","lakeperimeter.z",
                    "winter_dur_0_4", "secchi.z",
                    "depth.z","lakeshorelinefactor.z",
                    "total.dev.z","post_ice_warm_rate",
                    "height_19.3_23.3",
                    "max_surf_JulAugSep",
                    "total.ag.z",
                    "height_10.6_11.2",
                    "elevation.z",
                    "mean_gg_0c.z",
                    "IceOff",
                    "mean_surf.z",
                    "total.for.z",
                    "days_26_28",
                    "height_27_32")
  
  cp_plots_by_feature <- list()
  
  for(feat in all_features) {
    result <- create_ice_plot(feat)
    if(!is.null(result)) {
      print(result$plot)
      cp_plots_by_feature[[feat]] <- result
    }
  }
  

# Return results
results_cp <- list(
  plots = if(exists("cp_plots_by_feature")) cp_plots_by_feature else NULL,
  model = model
)

results_cp


##### random forest leave one out ####

n_obs <- nrow(Dat) # Number of observations in your dataset
predictions_loo <- numeric(n_obs) # To store predictions for each left-out observation

for (i in 1:n_obs) {
  # Create training and validation sets for the current fold
  train_data <- Dat[-i, ] # All data except the i-th observation
  validation_data <- Dat[i, ] # The i-th observation
  
  # Train the Random Forest model on the training data
  rf_model_loo <- randomForest(prop.of.both.x ~ ., data = train_data, ntree = 500)
  
  # Predict for the left-out observation
  predictions_loo[i] <- predict(rf_model_loo, newdata = validation_data)
}

#look at model output
print(rf_model_loo)
importance(rf_model_loo)
plot(rf_model_loo) #looking at error rate, 500 trees bring error to about zero
#plotting number of tree nodes
hist(treesize(rf_model_loo),
     main = "No. of Nodes for the Trees",
     col = "green")

#error 
a<-ggplot(data.frame(Predicted = rf_model_loo$predicted[1:314], Actual = Dat$prop.of.both.x[1:314]), 
          aes(x = Actual, y = Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Random Forest Predictions on Niche Overlap Dataset prop of both", 
       x = "Actual Values", y = "Predicted Values")
a
#Variable Importance
varImpPlot(rf_model_loo,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance for Random Forest Model (niche overlap)")
importance_matrix<-importance(rf_model_loo,scale = TRUE)
importance_df <- as.data.frame(importance_matrix)
importance_df$Feature<-rownames(importance_df)
importance_matrix <- data.frame(
  Rank = 1:nrow(importance_df),
  Feature = rownames(importance_df),
  IncNodePurity = importance_df$IncNodePurity,
  IncNodePurity_pct = round(importance_df$IncNodePurity * 100, 2)
)
randomForest::importance(rf_model_loo, type=2)
print(var_IMP)
varImpPlot(rf_model_loo)

#visualize rf model decison tee 
#plotting represetative model tree 
reprtree::plot.getTree(rf_model_loo)

##### interactions in Niche overlap model #####
# IML Interaction Analysis using original data structure

# Use training data that the model was actually trained on
X_train <- train[, !names(train) %in% "prop.of.both.x"]
y_train <- train$prop.of.both.x

# Create predictor with explicit predict function
predictor <- Predictor$new(
  model = rf,
  data = X_train,
  y = y_train,
  predict.function = function(model, newdata) {
    as.numeric(predict(model, newdata, type = "response"))
  }
)

# Test that predictor works
print("Testing predictor...")
test_pred <- predictor$predict(X_train[1:5, ])
print(test_pred)

# Now try interactions
print("Calculating interactions...")
interact <- Interaction$new(predictor)
plot(interact)


####Interaction statistics using hstats 
#code from Liv Nyffeler
#system.time(
#  s <- hstats(rf, X = Dat[, -which(names(Dat) == "prop.of.both.x")], 
#              approx = TRUE)  # Set approx = TRUE for faster computation
#)
# Print the results
#print(s)

# Plot the interaction statistics 
#plot(s)  # Or summary(s) for numeric output
#summary(s)

# Print the results
#print(s)

#Make pdps based on interactions in the data 
#make mars model 
propboth.mars <- earth(prop.of.both.x ~. , data = Dat) 
summary(propboth.mars)
# Get variable importance using evimp
importance <- evimp(propboth.mars, trim = TRUE)
# Print the importance scores
print(importance)
# Plot the importance scores
plot(importance)

library(doParallel) # load the parallel backend
cl <- makeCluster(4) # use 4 workers
registerDoParallel(cl) # register the parallel backend

#Comb OGT influence on niche overlap
pdp::partial(propboth.mars, pred.var = c("thermalguild","total.dev.z"), plot = TRUE,
        chull = TRUE, parallel = TRUE, paropts = list(.packages = "earth")) # Figure showing most important variables 
stopCluster(cl) # good practice

#mean_epi_hypo_ratio _Comb CT influence on niche overlap
cl <- makeCluster(4) # use 4 workers
registerDoParallel(cl) # register the parallel backend
pdp::partial(propboth.mars, pred.var = c("CTdiff","depth.z"), plot = TRUE,
        chull = TRUE, parallel = TRUE, paropts = list(.packages = "earth")) # Figure showing most important variables 
stopCluster(cl) # good practice

#CombCT influence on niche overlap
cl <- makeCluster(4) # use 4 workers
registerDoParallel(cl) # register the parallel backend
pdp::partial(propboth.mars, pred.var = c("CTdiff","lakearea.z"), plot = TRUE,
        chull = TRUE, parallel = TRUE, paropts = list(.packages = "earth")) # Figure showing most important variables 
stopCluster(cl) # good practice


#Lake area influence on niche overlap
cl <- makeCluster(4) # use 4 workers
registerDoParallel(cl) # register the parallel backend

pdp::partial(propboth.mars, pred.var = c("total.dev.z"), plot = TRUE,
        chull = TRUE, parallel = TRUE, paropts = list(.packages = "earth")) # Figure showing most important variables 
stopCluster(cl) # good practice

#Lake area:depth influence on niche overlap
cl <- makeCluster(4) # use 4 workers
registerDoParallel(cl) # register the parallel backend

pdp::partial(propboth.mars, pred.var = c("CombCT","depth.z"), plot = TRUE,
        chull = TRUE, parallel = TRUE, paropts = list(.packages = "earth")) # Figure showing most important variables 
stopCluster(cl) # good practice

#Lake area:Comb OGT influence on niche overlap
cl <- makeCluster(4) # use 4 workers
registerDoParallel(cl) # register the parallel backend

pdp::partial(propboth.mars, pred.var = c("mean_epi_hypo_ratio","CombCT"), plot = TRUE,
        chull = TRUE, parallel = TRUE, paropts = list(.packages = "earth")) # Figure showing most important variables 
stopCluster(cl) # good practice


# Add a label to the colorkey
lattice::trellis.focus("legend", side = "right", clipp.off = TRUE, highlight = FALSE)
grid::grid.text("proportion of overlap of both species", x = 0.2, y = 1.05, hjust = 0.85, vjust = 1)
lattice::trellis.unfocus()

cl <- makeCluster(4) # use 4 workers
registerDoParallel(cl) # register the parallel backend

pdp::partial(propboth.mars, pred.var = c("lakeperimeter.z"), plot = TRUE,
        chull = TRUE, parallel = TRUE, paropts = list(.packages = "earth")) # Figure showing most important variables 
stopCluster(cl) # good practice


#cl <- makeCluster(4) # use 4 workers
#registerDoParallel(cl) # register the parallel backend

#partial(propboth.mars, pred.var = c("lakearea.z", "total.dev.z","post_ice_warm_rate"), plot = TRUE,
        #chull = TRUE, parallel = TRUE, paropts = list(.packages = "earth")) # Figure showing most important variables 
#stopCluster(cl) # good practice

##### ICE curves ####

#individual conditional expectation plots dispay one line per instance that shows how the instance's prediction
#changes when a feature changes, see the dependence of the prediction on a feature for each instance seperately, resutling in oe line per
#instance of a dataset

#use partial to obtain ICE curves 
pred.fun <- function(object, newdata) {
  mean(predict(object, newdata), na.rm = TRUE)
}
pred.ice<-function(propboth.mars, test) predict(propboth.mars,test)
lakeperim.ice<-pdp::partial(propboth.mars,pred.var = "lakeperimeter.z")

pdp::plotPartial(lakeperim.ice,rug=T,train=train,alpha=0.3)

#check range of veriable to fix arguement before plotting
pdp::partial(propboth.mars,pred.var ="lakearea.z",pred.grid = data.frame(lakearea.z=0:5.5),plot = TRUE, chull=TRUE)
pdp::partial(propboth.mars,pred.var ="total.dev.z",pred.grid = data.frame(total.dev.z=-1:2),plot = TRUE, chull=TRUE)

#plotting thermal guild overlap pair densities 
legend_title<-"Thermal guild of species pairs"
#get pairs with overlapping area
DAT_percent <-Dat[Dat$prop.of.both.x>0.05,]
a<-ggplot(data=DAT_percent,aes(thermalguild, fill=thermalguild))+geom_bar()+scale_fill_manual(legend_title,values = c("lightblue", "royalblue", "blue", "blue3" ,"darkblue"))
#create wrapper percent formatting function from the scales package 
pct_format = scales::percent_format(accuracy = .1)
#adding in text to show percent and raw values of thermal guild pairs observed in the niche overlap calcs 
a <- a + geom_text(
  aes(fontface="bold",
    label = sprintf(
      '%d (%s)',
      ..count..,
      pct_format(..count.. / sum(..count..))
    )
  ),
  stat = 'count',
  nudge_y = 5,
  colour = 'black',
  size = 5, 
)
a + theme_bw()


### TP and LR ######
#examine the baseline niche metric using all parameters of interest 

#### littoral reliance 
Dat<-read.csv("RF_model_clean_df.csv")
Dat<-Dat[,c(3:ncol(Dat))]

#get columns of interest
#metric testing littoral reliance 
Dat<-Dat[,c(1,18,22,30,82,92,109,133,142,171,207,209:219,221:222)]

TP<-read.csv("TrophicP_LittorialR_median.csv")
TP$community<-TP$group

New<-merge(TP,Dat, all.x = TRUE, by="community")
medianNew<-merge(TP,Dat, all.x = TRUE, by="community")
Test<-New[,c(1,4,7,8:ncol(New))]
Test<-distinct(Test)
Test<-Test[Test$consumer%in%c("WAE","YEP","LMB","BLG","SMB","BLC","NOP","CIS"),]
#Test<-Test[Test$consumer%in%c("WAE","LMB"),]
Test<-Test[!is.na(Test$littoral_reliance),]
Test$species_code<-Test$consumer
Test<-Test[,c(1,3:ncol(Test))]

#add in physiological metrics 
Physi<-read.csv("fish_tempCorrected.csv")
Physi<-Physi[,c(4,5:7,18)] # get columns of interest 
# Apply to the thermalguild column
Physi$thermalguild<-as.factor(Physi$Temperature.preference.class)

#Join LR and physio data
Test<-inner_join(Test,Physi,by="species_code")
Test<-Test[,c(1:28,30,31)]#get individual species thermal metrics and dropping combined phsiology metrics
Test<-distinct(Test)
#use lakes only from overlap analysis and make sure its a factor variable
Test<-Test[Test$community%in%c(Dat$community),]
Test$community<-as.factor(Test$community)
Test$species<-as.factor(Test$species_code)
Test<-Test[,c(1:25,27:31)]

# Impute missing values using rfImpute() note* can only do this for numeric variables 
Test <- Test %>%
  mutate_at(c(11,25),as.numeric) #changing integer vectors to numeric
Test <- rfImpute(littoral_reliance ~ ., data = Test, iter = 5)

# Split the dataset into training and testing sets, leave one out 
n_obs <- nrow(Test) # Number of observations in your dataset
predictions_loo <- numeric(n_obs) # To store predictions for each left-out observation

for (i in 1:n_obs) {
  # Create training and validation sets for the current fold
  train_data <- Test[-i, ] # All data except the i-th observation
  validation_data <- Test[i, ] # The i-th observation
  
  # Train the Random Forest model on the training data
  rf_model_loo <- randomForest(littoral_reliance ~ ., data = train_data, ntree = 500, importance=T)
  
  # Predict for the left-out observation
  predictions_loo[i] <- predict(rf_model_loo, newdata = validation_data)
}

#look at model output
print(rf_model_loo)
imp<-importance(rf_model_loo, type = 2)
plot(rf_model_loo, main="Littoral reliance rf model error") #looking at error rate, 500 trees bring error to about zero
#plotting number of tree nodes
hist(treesize(rf_model_loo),
     main = "No. of Nodes for the Trees (littoral reliance)",
     col = "green")

#error 
a<-ggplot(data.frame(Predicted = predictions_loo[1:162], Actual = Test$littoral_reliance[1:162]), 
          aes(x = Actual, y = Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Random Forest Predictions on Littoral Reliance", 
       x = "Actual Values", y = "Predicted Values")
a
#Variable Importance
varImpPlot(rf_model_loo,
           sort = T,
           n.var = 28,
           main = "Top 10 - Variable Importance for Random Forest Model (littoral reliance)", 
           type = 1)
importance_matrix<-importance(rf_model_loo,scale = TRUE)
importance_df <- as.data.frame(importance_matrix)
importance_df$Feature<-rownames(importance_df)
importance_matrix <- data.frame(
  Rank = 1:nrow(importance_df),
  Feature = rownames(importance_df),
  IncNodePurity = importance_df$IncNodePurity,
  IncNodePurity_pct = round(importance_df$IncNodePurity * 100, 2)
)
# Sort by IncNodePurity
importance_matrix <- importance_matrix[order(-importance_matrix$IncNodePurity), ]

# Reassign ranks after sorting
importance_matrix$Rank <- 1:nrow(importance_matrix)
print(importance_matrix)

# Save as CSV 
write.csv(importance_df, "feature_importance_ranked_05212026_LR.csv", row.names = FALSE)

# Visualize the predicted vs actual values
d<-ggplot(data.frame(Predicted = predictions_loo, Actual = Test$littoral_reliance[1:162]), 
          aes(x = Actual, y = Predicted))+
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Random Forest Predictions on Niche Overlap Dataset Littoral Reliance", 
       x = "Actual Values", y = "Predicted Values")
d


####Interaction statistics using hstats code from Liv Nyffeler

system.time(
  s <- hstats(rf_model_loo, X = Test[, -which(names(Test) == "littoral_reliance")], 
              approx = TRUE)  # Set approx = TRUE for faster computation
)
# Print the results
print(s)

# Plot the interaction statistics 
plot(s)  # Or summary(s) for numeric output
summary(s)

#plot model features together
X=Test[which(names(Test) !="littoral_reliance")]
model=Predictor$new(rf_model_loo,data=X,y=Test$littoral_reliance)
effect=FeatureEffects$new(model)
# Create the feature effects plot - keep original x-axis labels
A <- effect$plot(features = c("OGT", "FTP", "Ctmax")) + #top 3 variables 
  plot_annotation(
    title = "Random forest model feature effects on littoral reliance",
    subtitle = "Physiology parameter model features",
    caption=" " #add in space for predictor universal x axis lab
  )

# Print the plot
print(A)

# Add universal x-axis label below everything using grid
grid.text("Predictor Variables", 
          x = 0.5, 
          y = 0.03,  # Adjust this value to position it below the individual labels
          gp = gpar(fontsize = 10, fontface = "bold"))
# Add universal y-axis label
grid.text("Predicted Littoral Reliance", 
          x = 0.007,  # Position on left side
          y = 0.5,   # Centered vertically
          rot = 90,  # Rotate 90 degrees
          gp = gpar(fontsize = 10, fontface = "bold"))

# Create the feature effects plot - keep original x-axis labels
B <- effect$plot(features=c("depth.z","mean_epi_hypo_ratio","stratification_duration")) + #top 6 (next top 3)
  plot_annotation(
    title = "Random forest model feature effects on littoral reliance",
    subtitle = "Numerical model features",
    caption=" " #add in space for predictor universal x axis lab
  )

# Print the plot
print(B)

# Add universal x-axis label below everything using grid
grid.text("Predictor Variables", 
          x = 0.5, 
          y = 0.03,  # Adjust this value to position it below the individual labels
          gp = gpar(fontsize = 10, fontface = "bold"))
# Add universal y-axis label
grid.text("Predicted Littoral Reliance", 
          x = 0.007,  # Position on left side
          y = 0.5,   # Centered vertically
          rot = 90,  # Rotate 90 degrees
          gp = gpar(fontsize = 10, fontface = "bold"))
### Ceteris Paribus (Individual Conditional Expectation) Plots by Community
# Ceteris Paribus (ICE) Plots with Different Curves for Each Community

# Now create function for all features
create_ice_plot <- function(feature_name) {
  tryCatch({
    print(paste("\nCreating plot for:", feature_name))
    
    ice_effect <- FeatureEffect$new(
      model,
      feature = feature_name,
      method = "ice",
      grid.size = 30
    )
    
    ice_data <- ice_effect$results
    ice_data$row_index <- as.numeric(as.character(ice_data$.id))
    ice_data$community <- Test$community[ice_data$row_index]
    
    # Find feature column
    feature_col <- names(ice_data)[!names(ice_data) %in% 
                                     c(".type", ".id", ".value", "community", "row_index")]
    
    if(length(feature_col) > 0) {
      # Use the actual feature column name
      print(paste("Using feature column:", feature_col[1]))
      
      # Add feature_value first
      ice_data$feature_value <- ice_data[[feature_col[1]]]
      
      # Calculate PDP more simply
      pdp_data <- aggregate(.value ~ feature_value, data = ice_data, FUN = mean)
      colnames(pdp_data)[colnames(pdp_data) == ".value"] <- "pdp_value"
      
      # Create plot
      p <- ggplot() +
        geom_line(data = ice_data, 
                  aes(x = feature_value, y = .value, color = community, group = .id), 
                  alpha = 0.4, linewidth = 0.5) +
        geom_line(data = pdp_data, 
                  aes(x = feature_value, y = pdp_value),
                  color = "black", linewidth = 2, 
                  linetype = "dashed") +
        labs(
          title = paste("ICE Plot:", feature_name),
          subtitle = "Individual curves by community (black dashed = average PDP)",
          x = feature_name,
          y = "Predicted littoral reliance",
          color = "Community"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          legend.position = "right",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 9, face = "bold"),
          legend.key.height = unit(0.3, "cm")
        )
      
      return(list(plot = p, data = ice_data, pdp = pdp_data))
    } else {
      return(NULL)
    }
    
  }, error = function(e) {
    print(paste("Error:", e$message))
    return(NULL)
  })
}

# Create plots for all features
all_features <- c("CTmax","FTP", "OGT", "mean_epi_hypo_ratio", 
                  "lakearea.z","stratification_duration","lakeperimeter.z",
                  "winter_dur_0_4", "secchi.z",
                  "depth.z","lakeshorelinefactor.z",
                  "total.dev.z","post_ice_warm_rate",
                  "height_19.3_23.3",
                  "max_surf_JulAugSep",
                  "total.ag.z",
                  "height_10.6_11.2",
                  "elevation.z",
                  "mean_gg_0c.z",
                  "IceOff",
                  "mean_surf.z",
                  "total.for.z",
                  "days_26_28",
                  "height_27_32")

cp_plots_by_feature <- list()

for(feat in all_features) {
  result <- create_ice_plot(feat)
  if(!is.null(result)) {
    print(result$plot)
    cp_plots_by_feature[[feat]] <- result
  }
}


# Return results
results_cp <- list(
  plots = if(exists("cp_plots_by_feature")) cp_plots_by_feature else NULL,
  model = model
)

results_cp


### trophic position
#### looking at littoral reliance 
Dat<-read.csv("RF_model_clean_df.csv")
Dat<-Dat[,c(3:ncol(Dat))]

#get columns of interest
#metric testing littoral reliance 
Dat<-Dat[,c(1,18,22,30,82,92,109,133,142,171,207,209:219,221:222)]

TP<-read.csv("TrophicP_LittorialR_median.csv")
TP$community<-TP$group

New<-merge(TP,Dat, all.x = TRUE, by="community")
Test<-New[,c(1,4,5,8:ncol(New))]
Test<-distinct(Test)
Test<-Test[Test$consumer%in%c("WAE","YEP","LMB","BLG","SMB","BLC","NOP","CIS"),]
#Test<-Test[Test$consumer%in%c("WAE","LMB"),]
Test<-Test[!is.na(Test$median),]
Test$species_code<-Test$consumer
Test<-Test[,c(1,3:ncol(Test))]

#add in physiological metrics 
Physi<-read.csv("fish_tempCorrected.csv")
Physi<-Physi[,c(4,5:7,18)] # get columns of interest 
# Apply to the thermalguild column
Physi$thermalguild<-as.factor(Physi$Temperature.preference.class)

#Join LR and physio data
Test<-inner_join(Test,Physi,by="species_code")
Test<-Test[,c(1:28,30,31)]#get individual species thermal metrics and dropping combined phsiology metrics
Test<-distinct(Test)
#use lakes only from overlap analysis and make sure its a factor variable
Test<-Test[Test$community%in%c(Dat$community),]
Test$community<-as.factor(Test$community)
Test$species<-as.factor(Test$species_code)
Test<-Test[,c(1:25,27:31)]

# Impute missing values using rfImpute() note* can only do this for numeric variables 
Test <- Test %>%
  mutate_at(c(3,4,11,24,25),as.numeric) #changing integer vectors to numeric
Test <- rfImpute(median ~ ., data = Test, iter = 5)

# Split the dataset into training and testing sets, leave one out 
n_obs <- nrow(Test) # Number of observations in your dataset
predictions_loo <- numeric(n_obs) # To store predictions for each left-out observation

for (i in 1:n_obs) {
  # Create training and validation sets for the current fold
  train_data <- Test[-i, ] # All data except the i-th observation
  validation_data <- Test[i, ] # The i-th observation
  
  # Train the Random Forest model on the training data
  rf_model_loo <- randomForest(median ~ ., data = train_data, ntree = 500, importance=T)
  
  # Predict for the left-out observation
  predictions_loo[i] <- predict(rf_model_loo, newdata = validation_data)
}

#look at model output
print(rf_model_loo)
imp<-importance(rf_model_loo, type = 2)
plot(rf_model_loo, main="Trophic Position rf model error") #looking at error rate, 500 trees bring error to about zero
#plotting number of tree nodes
hist(treesize(rf_model_loo),
     main = "No. of Nodes for the Trees (trophic position)",
     col = "green")

#error 
a<-ggplot(data.frame(Predicted = predictions_loo[1:162], Actual = Test$median[1:162]), 
          aes(x = Actual, y = Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Random Forest Predictions on Trophic Position", 
       x = "Actual Values", y = "Predicted Values")
a
#Variable Importance
varImpPlot(rf_model_loo,
           sort = T,
           n.var = 28,
           main = "Top 10 - Variable Importance for Random Forest Model (trophic position)", 
           type = 1)
importance_matrix<-importance(rf_model_loo,scale = TRUE)
importance_df <- as.data.frame(importance_matrix)
importance_df$Feature<-rownames(importance_df)
importance_matrix <- data.frame(
  Rank = 1:nrow(importance_df),
  Feature = rownames(importance_df),
  IncNodePurity = importance_df$IncNodePurity,
  IncNodePurity_pct = round(importance_df$IncNodePurity * 100, 2)
)
# Sort by IncNodePurity
importance_matrix <- importance_matrix[order(-importance_matrix$IncNodePurity), ]

# Reassign ranks after sorting
importance_matrix$Rank <- 1:nrow(importance_matrix)
print(importance_matrix)

# Save as CSV 
write.csv(importance_df, "feature_importance_ranked_05212026_TP.csv", row.names = FALSE)

#plot model features together
X=Test[which(names(Test) !="trophic position")]
model=Predictor$new(rf_model,data=X,y=Test$alpha.median)
effect=FeatureEffects$new(model)
# Create the feature effects plot - keep original x-axis labels
A <- effect$plot(features = c("OGT", "FTP", "Ctmax")) + #top 3 variables 
  plot_annotation(
    title = "Random forest model feature effects on trophic position",
    subtitle = "Top 3 numerical model features",
    caption=" " #add in space for predictor universal x axis lab
  )

# Print the plot
print(A)

# Add universal x-axis label below everything using grid
grid.text("Predictor Variables", 
          x = 0.5, 
          y = 0.03,  # Adjust this value to position it below the individual labels
          gp = gpar(fontsize = 10, fontface = "bold"))
# Add universal y-axis label
grid.text("Predicted Trophic Position", 
          x = 0.007,  # Position on left side
          y = 0.5,   # Centered vertically
          rot = 90,  # Rotate 90 degrees
          gp = gpar(fontsize = 10, fontface = "bold"))

# Create the feature effects plot - keep original x-axis labels
B <- effect$plot(features=c("depth.z","mean_epi_hypo_ratio","stratification_duration")) + #top 6 (next top 3)
  plot_annotation(
    title = "Random forest model feature effects on trophic position",
    subtitle = "Numerical model features",
    caption=" " #add in space for predictor universal x axis lab
  )

# Print the plot
print(B)

# Add universal x-axis label below everything using grid
grid.text("Predictor Variables", 
          x = 0.5, 
          y = 0.03,  # Adjust this value to position it below the individual labels
          gp = gpar(fontsize = 10, fontface = "bold"))
# Add universal y-axis label
grid.text("Predicted Trophic Position", 
          x = 0.007,  # Position on left side
          y = 0.5,   # Centered vertically
          rot = 90,  # Rotate 90 degrees
          gp = gpar(fontsize = 10, fontface = "bold"))
### Ceteris Paribus (Individual Conditional Expectation) Plots by Community
# Ceteris Paribus (ICE) Plots with Different Curves for Each Community

# Now create function for all features
create_ice_plot <- function(feature_name) {
  tryCatch({
    print(paste("\nCreating plot for:", feature_name))
    
    ice_effect <- FeatureEffect$new(
      model,
      feature = feature_name,
      method = "ice",
      grid.size = 30
    )
    
    ice_data <- ice_effect$results
    ice_data$row_index <- as.numeric(as.character(ice_data$.id))
    ice_data$community <- Test$community[ice_data$row_index]
    
    # Find feature column
    feature_col <- names(ice_data)[!names(ice_data) %in% 
                                     c(".type", ".id", ".value", "community", "row_index")]
    
    if(length(feature_col) > 0) {
      # Use the actual feature column name
      print(paste("Using feature column:", feature_col[1]))
      
      # Add feature_value first
      ice_data$feature_value <- ice_data[[feature_col[1]]]
      
      # Calculate PDP more simply
      pdp_data <- aggregate(.value ~ feature_value, data = ice_data, FUN = mean)
      colnames(pdp_data)[colnames(pdp_data) == ".value"] <- "pdp_value"
      
      # Create plot
      p <- ggplot() +
        geom_line(data = ice_data, 
                  aes(x = feature_value, y = .value, color = community, group = .id), 
                  alpha = 0.4, linewidth = 0.5) +
        geom_line(data = pdp_data, 
                  aes(x = feature_value, y = pdp_value),
                  color = "black", linewidth = 2, 
                  linetype = "dashed") +
        labs(
          title = paste("ICE Plot:", feature_name),
          subtitle = "Individual curves by community (black dashed = average PDP)",
          x = feature_name,
          y = "Predicted trophic position",
          color = "Community"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          legend.position = "right",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 9, face = "bold"),
          legend.key.height = unit(0.3, "cm")
        )
      
      return(list(plot = p, data = ice_data, pdp = pdp_data))
    } else {
      return(NULL)
    }
    
  }, error = function(e) {
    print(paste("Error:", e$message))
    return(NULL)
  })
}

# Create plots for all features
all_features <- c("Ctmax","FTP", "OGT", "mean_epi_hypo_ratio", 
                  "lakearea.z","stratification_duration","lakeperimeter.z",
                  "winter_dur_0_4", "secchi.z",
                  "depth.z","lakeshorelinefactor.z",
                  "total.dev.z","post_ice_warm_rate",
                  "height_19.3_23.3",
                  "max_surf_JulAugSep",
                  "total.ag.z",
                  "height_10.6_11.2",
                  "elevation.z",
                  "mean_gg_0c.z",
                  "IceOff",
                  "mean_surf.z",
                  "total.for.z",
                  "days_26_28",
                  "height_27_32")

cp_plots_by_feature <- list()

for(feat in all_features) {
  result <- create_ice_plot(feat)
  if(!is.null(result)) {
    print(result$plot)
    cp_plots_by_feature[[feat]] <- result
  }
}


# Return results
results_cp <- list(
  plots = if(exists("cp_plots_by_feature")) cp_plots_by_feature else NULL,
  model = model
)

results_cp


#plotting summary statisctics of TP and LR values 

# Lollipop Charts for Range of Median and Littoral Reliance by Group

# Calculate range of median values per group
median_range <- TP %>%
  group_by(group) %>%
  summarise(
    min_median = min(median, na.rm = TRUE),
    max_median = max(median, na.rm = TRUE),
    range_median = max_median - min_median,
    .groups = "drop"
  ) %>%
  arrange(desc(range_median))

print(median_range)

# Calculate range of littoral reliance values per group
littoral_col <- grep("littoral", names(TP), ignore.case = TRUE, value = TRUE)[1]

if(is.na(littoral_col) || length(littoral_col) == 0) {
  print("Warning: Could not find littoral reliance column. Using placeholder.")
  print("Available columns:")
  print(names(TP))
  littoral_col <- "littoral_reliance"  # Adjust this to your actual column name
}

print(paste("Using column for littoral reliance:", littoral_col))

littoral_range <- TP %>%
  group_by(group) %>%
  summarise(
    min_littoral = min(littoral_reliance, na.rm = TRUE),
    max_littoral = max(littoral_reliance, na.rm = TRUE),
    range_littoral = max_littoral - min_littoral,
    .groups = "drop"
  ) %>%
  arrange(desc(range_littoral))

print(littoral_range)

# PLOT 1: Lollipop chart for Median Range
p1 <- ggplot(median_range, aes(x = reorder(group, range_median), y = range_median)) +
  geom_segment(aes(x = reorder(group, range_median), xend = reorder(group, range_median),
                   y = 0, yend = range_median),
               color = "steelblue", linewidth = 1.5) +
  geom_point(size = 4, color = "steelblue") +
  coord_flip() +
  labs(
    title = "Range of Trophic Position values",
    subtitle = "Mixing Model output",
    x = "Lake-Survey Year",
    y = "Trophic Position samples range [all species] (Max - Min)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p1)

# PLOT 2: Lollipop chart for Littoral Reliance Range
p2 <- ggplot(littoral_range, aes(x = reorder(group, range_littoral), y = range_littoral)) +
  geom_segment(aes(x = reorder(group, range_littoral), xend = reorder(group, range_littoral),
                   y = 0, yend = range_littoral),
               color = "coral", linewidth = 1.5) +
  geom_point(size = 4, color = "coral") +
  coord_flip() +
  labs(
    title = "Range of Littoral Reliance Values",
    subtitle = "Mixing Model output",
    x = "Lake-Survey Year",
    y = "Littoral Reliance samples range [all species] (Max - Min)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p2)

#summary plots TP
ggplot(data=TP, aes(x=consumer,y=median, fill=consumer))+geom_boxplot()+labs(title="Range of Trophic position values")+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()+theme(legend.position = "none")

# Plot with custom 95% CI error bars
ggplot(TP, aes(x = consumer, y = median, fill = consumer)) +
  # Adds the standard box (IQR and median)
  geom_boxplot(coef=0, width = 0.5) +
  # Overrides the error bars (whiskers) to use 95% CI of the mean
  stat_summary(
    fun.data = "mean_cl_normal", # Calculates mean and 95% CI (uses t-distribution)
    geom = "errorbar",            # Uses error bar geometry
    width = 0.2,                  # Controls the width of the error bar tips
    color = "black"                 # Optional: change color for emphasis
  ) +
  # Optional: add a point for the mean
  stat_summary(
    fun = mean,
    geom = "point",
    color = "black"
  ) + labs(title="Range of Trophic position values wth 95% confidence intervals", y="Trophic position")+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()+theme(legend.position = "none")
#summary plots LR
ggplot(data=TP, aes(x=consumer,y=littoral_reliance, fill=39.7consumer))+geom_boxplot()+labs(title="Range of Littoral reliance values")+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()+theme(legend.position = "none")


# Store results
lollipop_results <- list(
  plot_median = p1,
  plot_littoral = p2,
  median_range_data = median_range,
  littoral_range_data = littoral_range
)

lollipop_results

model <- aov(median ~ group, data = TP) 
model2<-aov( median ~ group + consumer, data = TP)
anova(model,model2)

#### map of lakes ####
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(ggspatial)
library(cowplot)
library(RColorBrewer)
library(tigris)
library(rnaturalearthhires)

#get study lake coordinates
setwd("/Users/cammosley/Library/CloudStorage/GoogleDrive-cmosley@umn.edu/Shared drives/Hansen Lab/RESEARCH PROJECTS/Isotope cross lake analysis - CM/Data/TP models/")#macbook
Dat<-read.csv("RF_Data_092625.csv") #correct data for subset of lakes with TP and LR
Coord<-Dat[ ,c(2,206:207)]
Coord<-distinct(Coord)
Coord<-Coord[!is.na(Coord$lon),]
Coord$community<-substr(Coord$community,1,nchar(Coord$community)-5)

# Get MN boundary as sf - DL
mn_map <- map_data("state", region = "minnesota")

#make study lake coordinates to shape object 
my_sf <-st_as_sf(Coord, coords = c('lon','lat'))
my_sf <-st_set_crs(my_sf, value = 4326)

#first try 
ggplot(my_sf)+geom_sf(aes(color=community))

ggplot(data = world)+geom_sf()+
  geom_point(data=Coord,aes(x=lon,y=lat),size = 4)+coord_sf(xlim = c(-105,-85),ylim = c(40,50), expand = F)

#getting color gradient for points
# Choose a base palette (e.g., "Paired" from RColorBrewer)
base_colors <- brewer.pal(n = 12, name = "Paired")

# Expand the palette to 30 colors
my_palette <- colorRampPalette(base_colors)(26)

# View the hexadecimal color codes
print(my_palette)
library(ggplot2)
library(ggspatial)
library(rnaturalearth)
library(sf)
library(cowplot)

# Base Minnesota map
lakes <- ggplot() +
  # Plot the Minnesota map
  geom_polygon(data = mn_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black") +
  labs(x = "Longitude", y = "Latitude", 
       title = "Stable Isotope Food Web Study Lakes", 
       colour = "Lake Community") +
  geom_point(data = Coord, aes(x = lon, y = lat, colour = community), 
             size = 5.5, shape = 18) + 
  scale_color_manual(values = my_palette) +
  theme_minimal() +
  # Add the north arrow
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(.5, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  # Add the scale bar
  annotation_scale(location = "br", width_hint = 0.2, 
                   pad_x = unit(.5, "in"), pad_y = unit(0.2, "in")) +
  coord_fixed()

# Create US inset map
us_sf <- ne_states(country = "United States of America", returnclass = "sf") %>%
  st_transform(26915) %>% 
  filter(!name %in% c("Alaska", "Hawaii"))

mn_boundary <- us_sf %>% filter(name == "Minnesota")

us_inset <- ggplot() +
  geom_sf(data = us_sf, fill = "gray90", color = "white", size = 0.3) +
  geom_sf(data = mn_boundary, fill = "orange", color = "black", size = 0.5) +
  theme_void() +
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))

# Combine with inset
final_map <- ggdraw() +
  draw_plot(lakes) +
  draw_plot(us_inset, x = 0.6, y = 0.75, width = 0.25, height = 0.25)

final_map

###### BRMS ######
Dat<-read.csv("RF_model_clean_df_11182025.csv") #requires data formatting if you read in csv file

#remove rkb values 
Dat<-Dat[!Dat$sp1.x.x=="RKB",]
Dat<-Dat[!Dat$sp2.x=="RKB",]

#get columns of interest
#metric testing prop.of.both 
Dat <- Dat %>%
  dplyr::select(community,
                prop.of.both.x,
                pair,
                lakeperimeter.z,
                mean_gdd_0c.z,
                secchi.z,
                IceOff,
                post_ice_warm_rate,
                winter_dur_0_4,
                mean_epi_hypo_ratio,
                lakeshorelinefactor.z,
                elevation.z,
                depth.z,
                max_surf_jul,
                mean_surf.z,
                lakearea.z,
                total.for.z,
                total.dev.z,
                total.dev.z,
                stratification_duration,
                height_19.3_23.3,
                height_27_32,
                height_10.6_11.2,
                days_26_28,
                ZM,
                OGTdiff,
                CTdiff,
                FTPdiff,
                cold,
                thermalguild)

DatNA<-Dat[is.na(Dat$mean_epi_hypo_ratio),]
Dat<-distinct(Dat)
#corecing factor variables to integers for impute function and rf model 

#Dat$thermalguild<-unclass(factor(Dat$thermalguild))
#Dat$community<-unclass(factor(Dat$community))
#Dat$sp1.x.x<-unclass(factor(Dat$sp1.x.x))
#Dat$sp2.x<-unclass(factor(Dat$sp2.x))

Dat$thermalguild<-factor(Dat$thermalguild)
Dat$community<-factor(Dat$community)
Dat$pair<-factor(Dat$pair)
#Dat$sp1.x.x<-factor(Dat$sp1.x.x)
#Dat$sp2.x<-factor(Dat$sp2.x)

# Impute missing values using rfImpute()
# Note: rfImpute requires a response variable without missing values
Dat <- rfImpute(prop.of.both.x ~ ., data = Dat, iter = 5)

#model formula 
mod1<-brm(prop.of.both.x ~ pair +
          lakeperimeter.z +
          mean_gdd_0c.z +
          secchi.z +
          IceOff+
          post_ice_warm_rate+
          winter_dur_0_4+
          mean_epi_hypo_ratio+
          lakeshorelinefactor.z+
          elevation.z+
          depth.z+
          max_surf_jul+
          mean_surf.z+
          lakearea.z+
          total.for.z+
          total.dev.z+
          total.dev.z+
          stratification_duration+
          height_19.3_23.3+
          height_27_32+
          height_10.6_11.2+
          days_26_28+
          ZM+
          OGTdiff+
          CTdiff+
          FTPdiff+
          cold+
          thermalguild +
            (1 | community), data=Dat, 
          iter = 50000,
          warmup = 20000,
          family = gaussian(), 
          chains = 5, 
          save_pars = save_pars(all = TRUE), 
          control = list(max_treedepth = 15, adapt_delta = 0.99))
mod1fit<-add_criterion(mod1, "loo", moment_match = TRUE)
summary(mod1)
plot(mod1)
plot(conditional_effects(mod1))
ranef_mod1<-ranef(mod1)

pp_check(mod1)

mod2

#note run time is 1 day and rhats for most parameters are 2+ 

#### glmm #####

ScaledDat<- Dat %>% 
  mutate_if(is.numeric, scale)


lm1<-lmer(prop.of.both.x ~ pair +
            lakeperimeter.z +
            mean_gdd_0c.z +
            secchi.z +
            IceOff+
            post_ice_warm_rate+
            winter_dur_0_4+
            mean_epi_hypo_ratio+
            lakeshorelinefactor.z+
            elevation.z+
            depth.z+
            max_surf_jul+
            mean_surf.z+
            lakearea.z+
            total.for.z+
            total.dev.z+
            total.dev.z+
            stratification_duration+
            height_19.3_23.3+
            height_27_32+
            height_10.6_11.2+
            days_26_28+
            ZM+
            OGTdiff+
            CTdiff+
            FTPdiff+
            cold+
            thermalguild +
            (1 | community), 
          data=ScaledDat, REML = FALSE)
summary(lm1)
plot(lm1)

lm2<-lmer(prop.of.both.x ~ pair +
            lakeperimeter.z +
            mean_gdd_0c.z +
            secchi.z +
            IceOff+
            post_ice_warm_rate+
            winter_dur_0_4+
            mean_epi_hypo_ratio+
            lakeshorelinefactor.z+
            elevation.z+
            depth.z+
            max_surf_jul+
            mean_surf.z+
            lakearea.z+
            total.for.z+
            total.dev.z+
            total.dev.z+
            stratification_duration+
            height_19.3_23.3+
            height_27_32+
            height_10.6_11.2+
            days_26_28+
            ZM+
            OGTdiff+
            CTdiff+
            FTPdiff+
            cold+
            thermalguild + 
            (1 | pair) +
            (1 | community), 
          data=ScaledDat, REML = FALSE)

summary(lm2)
plot(lm2)

#visualize model predicted vs actual values for model fit assessment 
ScaledDat$pred<-predict(lm1)
ggplot(ScaledDat,aes(x=prop.of.both.x,y=pred,colour=community, group=community)) + geom_point() + geom_line() + theme(legend.position="bottom", legend.direction = "horizontal")

