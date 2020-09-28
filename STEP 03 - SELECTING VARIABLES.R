####################### VALE INSTITUTE OF TECHNOLOGY ##########################

########################     GEA-GPA TUTORIAL       ###########################
######################## STEP 03: SELECTING VARIABLES #########################

### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###



# PRE-ANALYSIS ---- 

#AIMS:
#A. EXTRACT VALUES FROM ENVIRONMENTAL VARIABLES (RASTERS LAYERS)
#B. SELECT UNCORRELATED VARIABLES FOR FOLLOWING STEPS


# LOAD PACKAGES:
library(r2vcftools)
library(adegenet)
library(ggfortify)
library(raster)
library(rgdal)
library(sp)
library(rgeos)
library(RStoolbox)
library(usdm)
library(tidyverse)


# LOAD FUNCTIONS:
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}





## 1. LOADING GEOGRAPHICAL INFORMATION -------------------------------------

#load adapative .vcf file filtered in step 1, with geographical information:
snps = vcfLink("vcf/pilocarpus_adap_filtered.vcf", overwriteID=T)
VCFsummary(snps) #264 individuals and 19025 SNPs.

# Transform Coordinates into Spatial Points Dataframe
coords_gen = snps@meta[,c(7,4:5)]
coordinates(coords_gen) = coords_gen[ ,2:3]
nrow(coords_gen) #264 individuals
class(coords_gen) #"SpatialPointsDataFrame"
projection(coords_gen) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m")
#verify
plot(coords_gen, axes=TRUE)








#2. EXTRACTING ENVIRONMENTAL VALUES -------------------------------------------

#------------------------------------------- CLIMATIC VARIABLES - ENVIREM:
# loading rasters
current.list = list.files(path="./climatic_Envir/", pattern =".tif", full.names=TRUE)
Bioclimate = stack(current.list)
names(Bioclimate)
projection(Bioclimate) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m"

# Extracting climatic data
climatic_E = as.data.frame(extract(Bioclimate, coords_gen)) #extract can bug with tidyverse package
class(climatic_E)
head(climatic_E)
summary(climatic_E)
#Remove the variable "monthCountByTemp10" because all values are 12.
climatic_E = climatic_E[-10]
row.names(climatic_E) = snps@meta$ind_ID
head(climatic_E)

#save raw extracted values:
write.csv(climatic_E, file = "./Results/Outputs/extracted_values_climaticE.csv")


#------------------------------------------- CLIMATIC VARAIBLES - WorldClim 2.1:
# loading rasters
current.list = list.files(path="./climatic_WC/2.1/wc0.5", pattern =".tif", full.names=TRUE)
Bioclimate_WC2 = stack(current.list)
names(Bioclimate_WC2)
projection(Bioclimate_WC2) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m"

# Extracting climatic data
climatic_WC2 = as.data.frame(extract(Bioclimate_WC2, coords_gen))
class(climatic_WC2)
head(climatic_WC2)
summary(climatic_WC2)
row.names(climatic_WC2) = snps@meta$ind_ID
head(climatic_WC2)

#save raw extracted values:
write.csv(climatic_WC2, file = "./Results/Outputs/extracted_values_climaticW2.csv")








#3. VERIFY ENVIRONMENTAL HETEROGENEITY BY RASTERS PCA - OPTIONAL --------
# Choose a extension for the study area:
min(coords_gen@coords[,1])
max(coords_gen@coords[,1])
min(coords_gen@coords[,2])
max(coords_gen@coords[,2])

#define the extension:
ext = extent(-51, -49.5, -6.5, -5.5)

# Crop variables:
Bioclimate_crop = crop(Bioclimate, ext) #change the rasters

#load any shapefile or raster to visualization:
cangas = shapefile("maps/shapefiles/Cangas.shp") #if your shapefile is smaller than study area you do not need to crop. If you do it, the object will be NULL
#cangas = crop(cangas, ext) 
crs(cangas)

#if you need to fix the CRS base on a model, also a raster:
cangas = spTransform(cangas, proj4string(Bioclimate_crop))
crs(cangas)


#Perform a climate PCA
#spca = standard; nSamples = pixels number used, maskCheck = avoid NA; nComp = number of PC to use
Bioclimate_crop = Bioclimate_crop[[-10]] #remove "monthCountByTemp10" because all values are 12.
PCA_climate_E = rasterPCA(Bioclimate_crop, nSamples = NULL, nComp = 3, spca = F, maskCheck = TRUE)

#Save rasterPCA as raster multiband:
writeRaster(PCA_climate_E$map, filename="./Results/PCA/PCA_climate_E.tif", format="GTiff", overwrite=TRUE, options=c('TFW=YES'), bylayer=F)


#Save as a figure:
pdf ("./Results/PCA/Map_climate_E_pilocarpus.pdf", family="Helvetica", onefile = F)
#plot the map with the 3 first PC
plotRGB(PCA_climate_E$map, r = 1, g = 2, b = 3, alpha =255, stretch="hist", axes = F)
plot(cangas, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lwd=1, col = "beige")
points(coords_gen@coords[,1], coords_gen@coords[,2], col = 'white', bg = NULL, cex = 1.5, pch=21)

dev.off()









#4. MERGING DATASETS FOR GEA ANALYSES -----------------------------------------
#If you need to reload extracted values:
climatic_E = read.csv("./Results/Outputs/extracted_values_climaticE.csv", row.names = 1)
head(climatic_E)
climatic_WC2 = read.csv("./Results/Outputs/extracted_values_climaticW2.csv", row.names = 1)
head(climatic_WC2)
soil = read.csv("./Results/Outputs/imputed_Soil.csv", row.names = 1)
head(soil)

#verify environmental heterogeneity by histograms
pdf("./Results/Histograms/Histograms_ClimaticE.pdf", onefile = T)
ggplot(gather(climatic_E), aes(value)) + 
  geom_histogram(bins = 30) + 
  facet_wrap(~key, scales = 'free_x')+
  theme_bw()
dev.off()

#Remove variables with <5 classes: continentality
names(climatic_E)
climatic_E = climatic_E[,-4]
names(climatic_E)


#verify environmental heterogeneity by histograms
pdf("./Results/Histograms/Histograms_ClimaticWC2.pdf", onefile = T)
ggplot(gather(climatic_WC2), aes(value)) + 
  geom_histogram(bins = 30) + 
  facet_wrap(~key, scales = 'free_x')+
  theme_bw()
dev.off()

#Remove variables with <5 classes: all wind and vapor variables
names(climatic_WC2)
climatic_WC2 = climatic_WC2[,-c(20:55)]
names(climatic_WC2)

#verify environmental heterogeneity by histograms
pdf("./Results/Histograms/Histograms_Soil.pdf", onefile = T)
ggplot(gather(soil[,3:23]), aes(value)) + 
  geom_histogram(bins = 30) + 
  facet_wrap(~key, scales = 'free_x')+
  theme_bw()
dev.off()

#Remove variables with <5 classes: Keep all varaibles

#verify order of samples:
identical(rownames(climatic_E), rownames(climatic_WC2))
identical(rownames(climatic_E), as.character(soil$id))

#merge df:
GEA_vars = cbind(soil, climatic_E, climatic_WC2)
head(GEA_vars)

#save GEA variables:
write.csv(GEA_vars, file = "./Results/Outputs/GEA_vars_raw.csv")








#5. SELECTING GEA VARIABLES BY Carvalho et al. 2020 ------------------

#----------------------------------- SELECTING ENVIREM VARIABLES:
# Run PCA:
pca_E_BS = prcomp(climatic_E, center=TRUE, scale=TRUE)
summary(pca_E_BS) #verify if same PC% 

#Selecting PC by Broken Stick: 2PCs
screeplot(pca_E_BS, bstick=TRUE, type="lines")
screeplot(pca_E_BS, bstick=TRUE, type="barplot")
#Broken Stick Rule: Principal components should be retained as long as observed eigenvalues are higher than corresponding random broken stick components. Jackson 1993 & Legendre & Legendre 2012

summary(pca_E_BS)
#2PCs are 89.64% of variance

# Selecting variables with more contribution by PC
PC1_var = names(which.max(abs(pca_E_BS$rotation[,1])))
PC1_var
PC2_var = names(which.max(abs(pca_E_BS$rotation[,2])))
PC2_var

#save results
climatic_data_E = cbind(GEA_vars[,c(1:2)], climatic_E[PC1_var], climatic_E[PC2_var])
str(climatic_data_E)  
write.csv(climatic_data_E, "Results/Outputs/climatic_data_E_pilocarpus.csv")



#----------------------------------- SELECTING WorldClim 2.1 VARIABLES:
# Run PCA:
pca_WC2_BS = prcomp(climatic_WC2, center=TRUE, scale=TRUE)
summary(pca_WC2_BS) #verify if same PC% was recovered by dudi.pca function

#Selecting PC by Broken Stick:
screeplot(pca_WC2_BS, bstick=TRUE, type="lines")
screeplot(pca_WC2_BS, bstick=TRUE, type="barplot")
#Broken Stick Rule: Principal components should be retained as long as observed eigenvalues are higher than corresponding random broken stick components. Jackson 1993 & Legendre & Legendre 2012

summary(pca_WC2_BS)
#2PCs are 87.22% of variance


# Selecting variables with more contribution by PC
PC1_var = names(which.max(abs(pca_WC2_BS$rotation[,1])))
PC1_var
PC2_var = names(which.max(abs(pca_WC2_BS$rotation[,2])))
PC2_var


#save results
climatic_data_WC2 = cbind(GEA_vars[,c(1:2)], climatic_WC2[PC1_var], climatic_WC2[PC2_var])
str(climatic_data_WC2)  

write.csv(climatic_data_WC2, "Results/Outputs/climatic_data_WC2_pilocarpus.csv")


#----------------------------------- SELECTING SOIL VARIABLES:
# Run PCA:
pca_soil_BS = prcomp(soil[3:23], center=TRUE, scale=TRUE)
summary(pca_soil_BS) #verify if same PC%

#Selecting PC by Broken Stick: 2PCs
screeplot(pca_soil_BS, bstick=TRUE, type="lines")
screeplot(pca_soil_BS, bstick=TRUE, type="barplot")
#Broken Stick Rule: Principal components should be retained as long as observed eigenvalues are higher than corresponding random broken stick components. Jackson 1993 & Legendre & Legendre 2012

summary(pca_soil_BS)
#2PCs are 63.06% of variance

# Selecting variables with more contribution by PC
PC1_var = names(which.max(abs(pca_soil_BS$rotation[,1])))
PC1_var
PC2_var = names(which.max(abs(pca_soil_BS$rotation[,2])))
PC2_var

#save results
soil_data = cbind(GEA_vars[,c(1:2)], soil[PC1_var], soil[PC2_var])
str(soil_data)  

write.csv(soil_data, "Results/Outputs/soil_data_pilocarpus.csv")


#----------------------------------- VIF/CORRELATION AMONG VARIABLES:
# load data for GEA analysis
climatic_E_filtered = read.csv("Results/Outputs/climatic_data_E_pilocarpus.csv", row.names = 1)
head(climatic_E_filtered)

climatic_WC2_filtered = read.csv("Results/Outputs/climatic_data_WC2_pilocarpus.csv", row.names = 1)
head(climatic_WC2_filtered)

soil_filtered = read.csv("Results/Outputs/soil_data_pilocarpus.csv", row.names = 1)
head(soil_filtered)

# creating data 
data_cor = cbind(climatic_E_filtered[3:length(climatic_E_filtered)], climatic_WC2_filtered[3:length(climatic_WC2_filtered)], soil_filtered[3:length(soil_filtered)])
names(data_cor)

#Calculate VIF for all varaibles using correlation 0.7:
vif_cal = vifcor(data_cor, th = 0.7)
vif_cal #all variables with VIF < 2, max cor = -0.613, keep 5 variables; remove 1 variable.

#remove variables:
keepvars_cal = vif_cal@results$Variables
GEA_vars_cal = GEA_vars[colnames(GEA_vars) %in% keepvars_cal]

#merge individual and locality cols
GEA_vars_filtered = cbind(GEA_vars[,1:2], GEA_vars_cal)
#verify
length(GEA_vars_filtered)
head(GEA_vars_filtered)

#save GEA variables dataset:
write.csv(GEA_vars_filtered, file = "./GEA/Inputs/GEA_vars_filtered.csv")








#6. SELECTING GEA VARIABLES BY VIF --------------------------------------------
#If you need to load the file:
GEA_vars = read.csv("./Results/Outputs/GEA_vars_raw.csv", row.names = 1)
head(GEA_vars)

#Calculate VIF for all varaibles using threshold 0.7:
vif_vif = vifcor(GEA_vars[, -c(1:2)], th = 0.7)
vif_vif #all variables with VIF < 8, max cor = 0.691, keep 22 variables

#remove variables:
keepvars_vif = vif_vif@results$Variables
GEA_vars_vif = GEA_vars[colnames(GEA_vars) %in% keepvars_vif]

#merge individual and locality cols
GEA_vars_vif = cbind(GEA_vars[,1:2], GEA_vars_vif)
#verify
length(GEA_vars_vif)
head(GEA_vars_vif)

#save GEA variables dataset:
write.csv(GEA_vars_vif, file = "./GEA/Inputs/GEA_vars_vif.csv")








#7. GPA VARIABLES BY PILOCARPINE MODEL -----------------------------
#load leaf dataset:
leaf = read.csv("./Results/Outputs/imputed_Leaf.csv", row.names = 1)
#verify
head(leaf)

#Load predicted values for pilocarpine. From STEP2:
pilo = read.csv("./Results/Outputs/GPA_pilo_pred.csv", row.names = 1)
head(pilo)

#verify order of samples:
identical(as.character(pilo$id), as.character(leaf$id))

#merge df:
GPA_vars = cbind(leaf, pilo[,3])
names(GPA_vars)[14] = c("pilocarpine")
head(GPA_vars)


#verify environmental heterogeneity by histograms
pdf("./Results/Histograms/Histograms_Leaf.pdf", onefile = T)
ggplot(gather(GPA_vars[, -c(1:2)]), aes(value)) + 
  geom_histogram(bins = 30) + 
  facet_wrap(~key, scales = 'free_x')+
  theme_bw()
dev.off()

#keep all variables

#save raw extracted values:
write.csv(GPA_vars, file = "./Results/Outputs/phenotype_leaf.csv")





#9. SELECTING GPA VARAIBLES MANUALLY ---------------------------------------
#if you need to load:
GPA_vars = read.csv("./Results/Outputs/phenotype_leaf.csv", row.names = 1)
head(GPA_vars)
#More important variables: Zn, P, N, Mg, Ca, K, Cu in STEP2 + pilocarpine

#selecting variables:
names(GPA_vars)
leaf_vars = GPA_vars[,c(13, 4, 3, 7, 6, 5, 10, 14)]
names(leaf_vars)

#Calculate VIF for all varaibles using correlation 0.7:
vif_leaf_vars = vifcor(leaf_vars, th = 0.7)
vif_leaf_vars #all variables with VIF < 6, max cor = 0.6166, keep 8 variables

#remove variables:
keepvars_leaf_vars = vif_leaf_vars@results$Variables
GPA_leaf_vars = GPA_vars[colnames(GPA_vars) %in% keepvars_leaf_vars]

#merge individual and locality cols
GPA_leaf_vars = cbind(GPA_vars[,1:2], GPA_leaf_vars)
#verify
length(GPA_leaf_vars)
head(GPA_leaf_vars)

#save GEA variables dataset:
write.csv(GPA_leaf_vars, file = "./GPA/Inputs/GPA_vars_manual.csv")


#END