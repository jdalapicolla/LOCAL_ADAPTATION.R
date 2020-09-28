####################### VALE INSTITUTE OF TECHNOLOGY ##########################

########################      GEA-GPA TUTORIAL       ##########################
######################   STEP 02: PILOCARPINE VALUES   ########################

library(dplyr)
library(caret)
library(mgcv)
library(glmnet)
library(ranger)
library(gbm)
library(kernlab)

#### Load data ----
## Pilocarpine
piloc = read.table("./data_raw/jaborandi_database_piloc.csv", h=T, sep=";", dec=".") %>%
  dplyr::filter(!is.na(pilocarpina)) %>%
  dplyr::filter(!parcela %in% c("PCMA","LAMA", "PVNMA", "CPI", "PGPI"))

head(piloc)
nrow(piloc)

## Leaf data
leaf = read.table("./data_raw/jaborandi_database_folha.csv", h=T, sep=";", dec=".") %>%
  rename(cod_folhas = amostra) ## change "amostra" column name as in piloc dataframe
head(leaf)
nrow(leaf)

#merge pilocarpine and leaf tables
nrow(piloc); nrow(leaf)
df = left_join(piloc, leaf, by="cod_folhas") 
head(df)
nrow(df) ## OK

dfF <- df %>%
  dplyr::select(-c("cod_folhas", "cod_solo", "UTM1", "UTM2", 
            "Numero.Amostra", "Data.da.OS", "Data.Prog.Fim")) %>%
  dplyr::select(site, everything())

head(dfF)
nrow(dfF)

## Soil data
soil = read.table("./data_raw/jaborandi_database_solo.csv", h=T, sep=";", dec=".")
nrow(soil) ## Lower N - exclude!

#### Explore data----
str(dfF)
names(dfF)
## Check for outliers
names <- names(dfF[, 3:ncol(dfF)])
for(i in 1:length(names)){
  dotchart(dfF[, names[i]], main=paste(names[i]))
} 
## Outliers: in Cu_I, Fe_I, Zn_I
dfF[dfF$Cu_l>=80,] #53
dfF[dfF$Fe_l>=15000,] #53
dfF[dfF$Zn_l>=400,] #53

dfF <- dfF %>% dplyr::filter(Cu_l<=80) # Removing sample 53

##Normality
hist(dfF$pilocarpina) # OK

##### Model selection using machine learning and the caret package ----
## 10-fold validation with 5 repeats
## Pre-processing data by centering, scaling and using PCA axes to avoid collinearity
summary(dfF)

## Random Forest (no need to worry abour colinearity)
modelLookup("ranger")

set.seed(42)
model_rf <- train(
  pilocarpina ~ N_l + P_l + K_l + Ca_l + Mg_l + S_l + B_l + Zn_l + Fe_l
  + Mn_l + Cu_l,
  dfF,
  method = "ranger",
  importance = "impurity",
  trControl = trainControl(
    method = "repeatedcv", p=0.8, number = 10, repeats=5, verbose = TRUE
  ),
  tuneLength = 5
)

ggplot(model_rf)
plot(varImp(model_rf)) # Exclude B_l, Mn_l, Fe-L, S_l in subsequent models

## lm
modelLookup("lm")

set.seed(42)
model_lm <- train(
  pilocarpina ~ N_l + P_l + K_l + Ca_l + Mg_l + Zn_l + Cu_l,
  dfF,
  method = "lm",
  trControl = trainControl(
    method = "repeatedcv", p=0.8, number = 10, repeats=5, verbose = TRUE
  ),
  preProcess = c( "center", "scale", "pca")
)

min(model_lm$results$RMSE) ## Try to minimize RMSE
print(model_lm)

## glmnet
modelLookup("glmnet")

set.seed(42)
model_glmnet <- train(
  pilocarpina ~ N_l + P_l + K_l + Ca_l + Mg_l + Zn_l + Cu_l,
  dfF,
  method = "glmnet",
  #tuneLength=10,
  tuneGrid = expand.grid(
    lambda = 0:10 / 10,
    alpha = seq(0, 1, 0.1)
  ),
  trControl = trainControl(
    method = "repeatedcv", p=0.8, number = 10, repeats=5, verbose = TRUE
  ),
  preProcess = c( "center", "scale", "pca")
)

min(model_glmnet$results$RMSE) ## Try to minimize RMSE 
ggplot(model_glmnet)
print(model_glmnet)

## gam
modelLookup("gam")

set.seed(42)
model_gam <- train(
  pilocarpina ~ N_l + P_l + K_l + Ca_l + Mg_l + Zn_l + Cu_l,
  dfF,
  method = "gam",
  trControl = trainControl(
    method = "repeatedcv", p=0.8, number = 10, repeats=5, verbose = TRUE
  ),
  preProcess = c( "center", "scale", "pca"), 
  tuneLength = 5
)

min(model_gam$results$RMSE) ## Try to minimize RMSE
ggplot(model_gam)
print(model_gam)

## random forest
modelLookup("ranger")

set.seed(42)
model_rf <- train(
  pilocarpina ~ N_l + P_l + K_l + Ca_l + Mg_l + S_l + B_l + Zn_l + Fe_l
  + Mn_l + Cu_l,
  dfF,
  method = "ranger",
  importance = "impurity",
  trControl = trainControl(
    method = "repeatedcv", p=0.8, number = 10, repeats=5, verbose = TRUE
  ),
  #preProcess = c( "center", "scale"),
  tuneLength = 5
)

min(model_rf$results$RMSE) ## Try to minimize RMSE
ggplot(model_rf)
plot(varImp(model_rf))

## gradient boosting machines
modelLookup("gbm")

set.seed(42)
model_gbm <- train(
  pilocarpina ~ N_l + P_l + K_l + Ca_l + Mg_l + S_l + B_l + Zn_l + Fe_l
  + Mn_l + Cu_l,
  dfF,
  method = "gbm",
  trControl = trainControl(
    method = "repeatedcv", number = 10, repeats=5, verbose = TRUE
  ),
  #preProcess = c( "center", "scale", "pca"),
  tuneLength = 10
)

min(model_rf$results$RMSE) ## Try to minimize RMSE
ggplot(model_gbm)

## Support vector machines
modelLookup("svmPoly")

set.seed(42)
model_svm <- train(
  pilocarpina ~ N_l + P_l + K_l + Ca_l + Mg_l + S_l + B_l + Zn_l + Fe_l
  + Mn_l + Cu_l,
  dfF,
  method = "svmPoly",
  trControl = trainControl(
    method = "repeatedcv", p=0.8, number = 10, repeats=5, verbose = TRUE
  ),
  #preProcess = c( "center", "scale", "pca"),
  #tuneLength = 5
)

min(model_rf$results$RMSE) ## Try to minimize RMSE
ggplot(model_svm)

#### Compare models
model_list <- list(
  lm = model_lm,
  glmnet = model_glmnet,
  gam = model_gam,
  rf = model_rf, 
  gbm = model_gbm,
  svm = model_svm
)

## Collect resamples from the CV folds
resamps <- resamples(model_list)
resamps
summary(resamps)
dotplot(resamps, metric = "RMSE") ## svm is best model
bwplot(resamps, metric = "RMSE") 
xyplot(resamps, metric = "RMSE")


#### Predict new data using the best model ----
model_svm ## best model

#Load the genotyped samples:
new.data1 = read.csv(file = "./Results/Outputs/imputed_Leaf.csv", row.names = 1)

#Compare col names:
names(new.data1)
model_svm$call

#organized cols:
new.data = new.data1[,c(3:9,13,11,12,10)]
colnames(new.data) = c("N_l", "P_l", "K_l", "Ca_l", "Mg_l", "S_l", "B_l", "Zn_l", "Fe_l", "Mn_l", "Cu_l")

#verify
names(new.data)
model_svm$call

predicted.pilocarpine <- predict(model_svm, new.data) ## Predicted pilocarpina for new dataset using svm
predicted.pilocarpine ## These are the predicted values

#save results:
pilo_pred = cbind(new.data1[,1:2], predicted.pilocarpine)
head(pilo_pred)
write.csv(pilo_pred, file = "./Results/Outputs/GPA_pilo_pred.csv")

#END