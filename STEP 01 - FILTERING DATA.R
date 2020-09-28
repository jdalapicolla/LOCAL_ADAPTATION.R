####################### VALE INSTITUTE OF TECHNOLOGY ##########################

########################      GEA-GPA TUTORIAL       ##########################
########################   STEP 01: FILTERING DATA   ##########################

### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###



# PRE-ANALYSIS ---- 

#AIMS:
#A. CONTROL OUTLIERS AND NA's FROM PHENOTYPE (LEAF DATASET) AND ENVIRONMENTAL VARIABLES (SOIL DATASET)
#B. IMPUTE MISSING DATA IN BOTH DATASETS
#C. KEEP FOR NEXT STEP ONLY INDIVIDUALS THAT CONTAIN GENETIC, PHENOTYPE, AND ENVIRONMENTAL DATA 

#LOAD PACKAGES:
library(r2vcftools)
library(tidyverse)
library(missMDA)

#LOAD FUNCTIONS:
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}


#1. LOADING DATASETS ----
# Genetic data:
# load snps data set to take the snps and individuals ID
snps =  vcfLink("vcf/pilocarpus_filtered_adaptative.vcf", overwriteID=T)
VCFsummary(snps) ## 277 individuals and 19025 SNPs.

# verify name of samples
snps@meta$ind_ID

# estimating coverage depth by sample: 
coverage_ind = c()
for (p in 1:length(snps@sample_id)){
  beta = Query(Subset(snps, samples = p), type="site-mean-depth")
  coverage_ind[p] = mean(beta$MEAN_DEPTH, na.rm = T)}
#verify
coverage_ind
#save as metafile
snps@meta$coverage = coverage_ind

#save coverage depth by sample:
write.csv(as.data.frame(cbind(snps@meta$ind_ID,coverage_ind)), "./Results/Outputs/Depth_coverage_bysamples.csv")


## SOIL AND LEAF DATA:
#loading
leaf = read.csv("./data_raw/leaf_data.csv", h=T)
head(leaf)
soil = read.csv("./data_raw/soil_data.csv", h=T)
head(soil)

#merge soil and leaf tables
nrow(soil); nrow(leaf)
df = merge(soil, leaf, by.x="id", by.y="id")
head(df)
nrow(df) ## OK

#sort the Leaf data file according to the genetic file order:
df_sorted = df %>%
  slice(match(snps@meta$ind_ID, id))

head(df_sorted)
names(df_sorted)

#removing cols 'site' that are factor/class and not a numeric variable 
df_sorted = df_sorted[,-27]
names(df_sorted)

#checking if there are some difference between files. If they are corrected, function will return "character (0)". In our case a whole locality didn't soil variables estimated.
identical(snps@meta$ind_ID, as.character(df_sorted$id))

setdiff(snps@meta$ind_ID, df_sorted$id)
setdiff(df_sorted$id, snps@meta$ind_ID)

#There are some samples with genetic information but without phenotype information!

#save data merged
write.csv(df_sorted, "./Results/Outputs/Data_merged_pilocarpus.csv")


#2. CHECKING FOR OUTLIERS ----
#verify data:
str(df_sorted)
names(df_sorted)

#check for outliers
names = names(df_sorted[, 3:ncol(df_sorted)])

for(i in 1:length(names)){
  pdf(paste0("./Results/Outliers/Dotchart_",names[i],".pdf"), onefile = T)
  dotchart(as.numeric(df_sorted[, names[i]]), main=paste(names[i]))
  text(df_sorted[, names[i]], 1:length(as.vector(df_sorted[, names[i]])), labels = rownames(df_sorted), pos=2, cex=0.8, col="red")
  dev.off()
} 

## Replace the extreme points, replacing them with NAs, and load the file again for the next steps:

#Soil dataset:
#K: Ind 244
#Mg: Ind 14
#Na: Ind 188

#Leaf dataset:
#Cu: Ind 26 0
#Cu: Ind 63



#3. CHECKING FOR NA's ----
#load file without outliers:
df_clean = read.csv("./Results/Outputs/Data_merged_pilocarpus.csv", header = T, row.names = 1)
head(df_clean)

#counting NA by row (individuals) and columns (variables):
row = apply(df_clean, 1, function(x)sum(is.na(x)))
length(names(df_clean))-2 #35 variables
max(row) # maximum number of NA per individual was 4

col = apply(df_clean, 2, function(x)sum(is.na(x)))
length(df_clean[,1]) #264 individual
max(col) # 31 maximum number of NA per variable

#Verify variables position in df with too many NA's 
remove_vars = which(col>10) #I used 10 but you can change it 
remove_vars

#deleting cols with too many NA, in my case: argila, silte and areia:
names(df_clean)
df_clean = df_clean[-remove_vars]
names(df_clean)

#checking total number of NA
df_clean %>%
  summarise(sum(is.na(.)))
#9 missings values:

col = apply(df_clean, 2, function(x)sum(is.na(x)))
col
#7 Na's for soil
#2 Na's for leaf    



#4. IMPUTING NA's VALUES ----
# using the imputePCA function from package missMDA, choosing numeric cols
names(df_clean)
res.comp = imputePCA(df_clean[,3:length(df_clean)], ncp=2, method="Regularized")
# replace values
df_clean[,3:length(df_clean)] = res.comp$completeObs
head(df_clean)

#Checking for missing data
df_clean %>%
  summarise(sum(is.na(.)))
#0 NA

#save imputated data:
write.csv(df_clean, "./Results/Outputs/imputed_SoilLeaf.csv")



#5. SPLITING SOIL AND LEAF VALUES ----
# if you need to load the imputed data:
imputed_data = read.csv("./Results/Outputs/imputed_SoilLeaf.csv", row.names = 1)
#verify
names(imputed_data)

#split variable for leaf:
soil_imp = imputed_data[c(1:23)]
names(soil_imp)

#split variable for leaf:
leaf_imp = imputed_data[c(1:2,24:34)]
names(leaf_imp)

#save datasets:
write.csv(soil_imp, "./Results/Outputs/imputed_Soil.csv")
write.csv(leaf_imp, "./Results/Outputs/imputed_Leaf.csv")

#6. FILTERING GENETIC DATASET ----
# selecting individuals:
keep_ind = df_clean$id

# subsetting VCF:         
UNIND = snps@sample_id[snps@meta$ind_ID %in% keep_ind]
snps_adap_filtered = Subset(snps, samples=UNIND)

# saving vcf
Save(snps_adap_filtered, paste0("vcf/pilocarpus_adap_filtered.vcf"))

  
#END