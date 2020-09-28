####################### VALE INSTITUTE OF TECHNOLOGY ##########################


########################      GEA-GPA TUTORIAL       #########################
########################## STEP 05: LFMM FOR GEA #############################

### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###



# PRE-ANALYSIS-----------------------------------------------------------------


#AIMS:
#A. SELECTING SNPS UNDER SELECTION USING ENVIRONMENTAL VARIABLES
#B. EXTRACTING REFERENCE SEQUENCES AND LOCI THAT CONTAIN THESE LOCI 

## LOAD THE PACKAGES
library(r2vcftools)
library(LEA)
library(lfmm)
library(usdm)
library(vegan)
library(seqinr)
library(qvalue)

##5. LOAD FUNCTIONS TO BE USED ON THIS STEP.
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}





#1. LOAD GENOMIC DATASET-------------------------------------------------------
## Load snps data set to take the snps and individuals ID
snps =  vcfLink("vcf/pilocarpus_adap_filtered.vcf", overwriteID=T)
snps@meta
VCFsummary(snps) ## 264 individuals and 19025 SNPs.

## Read geno file
pilocarpus_geno = read.geno("vcf/pilocarpus_filtered_imputed.geno")
pilocarpus_geno[1:10,1:10]
colnames(pilocarpus_geno) = snps@site_id
rownames(pilocarpus_geno) = snps@sample_id
dim(pilocarpus_geno)

## Number of missing data
sum(is.na(pilocarpus_geno)) #0








#2. LOAD ENVIRONMENTAL DATASET------------------------------------------------ 
#loading extracted values per variable in Step 3 - Variables following Carvalho et al 2020:
GEA_variables = read.csv("GPA/Inputs/GPA_vars_manual.csv", row.names = 1)
head(GEA_variables)

#selecting only variables used in RDA to comparison
#names(GEA_variables)
#GEA_variables = GEA_variables[, c(3:6)]

#checking VIF:
vif(GEA_variables[, -c(1:2)]) #VIF  <2

#convert variables to PCs
#Selecting PC by Broken Stick: 3PCs
pca_BS = prcomp(GEA_variables[, -c(1:2)], center=TRUE, scale=TRUE)
summary(pca_BS)
stats::screeplot(pca_BS, bstick=TRUE, type="lines")
stats::screeplot(pca_BS, bstick=TRUE, type="barplot")
#Broken Stick Rule: Principal components should be retained as long as observed eigenvalues are higher than corresponding random broken stick components. Jackson 1993 & Legendre & Legendre 2012

summary(pca_BS)
#3PCs are 65.16% of variance

# Verify which varaibles are correlated to choosen PCs:
GEA_variables$PC1 = pca_BS$x[,1]
GEA_variables$PC2 = pca_BS$x[,2]
GEA_variables$PC3 = pca_BS$x[,3]

### Identify the most relevant environemntal variables, with >0.7 of correlation to the PC
variables_pca_cor = as.data.frame(cor(GEA_variables[, -c(1:2)]))
variables_pca_cor

PC1_var = row.names(variables_pca_cor[which(abs(variables_pca_cor[variables_pca_cor$PC1!=1, "PC1"]) > 0.7), ])
PC1_var
PC2_var = row.names(variables_pca_cor[which(abs(variables_pca_cor[variables_pca_cor$PC2!=1, "PC2"]) > 0.7), ])
PC2_var
PC3_var = row.names(variables_pca_cor[which(abs(variables_pca_cor[variables_pca_cor$PC3!=1, "PC3"]) > 0.7), ])
PC3_var

# PC1_var = "N_GPA"       "pilocarpine"
# PC2_var = "Cu_GPA"
# PC3_var = "Zn_GPA"


## selecting the PCs for analysis
pred_scale = as.data.frame (pca_BS$x[,1:3])






#3. LFMM2 ANALYSIS------------------------------------------------------------- 
## Carry out the lfmm using the package lfmm
mod.lfmm = lfmm_lasso( Y = pilocarpus_geno, #Sparse LFMM
                         X = pred_scale, 
                         K = 1)  ## K is the number of populations
  
## Performs association testing using the fitted model:
pv = lfmm_test(Y = pilocarpus_geno, 
                  X = pred_scale, 
                  lfmm = mod.lfmm,
                  calibrate="gif")



  

#4. OUTLIER DETECTIONS ------------------------------------------------------
## Estimate adjusted p-values based on gif
pv$gif# the GIFs for the predictors
#      ORIINAL      CALIBRATED   
#PC1 = 2.737822     1.8
#PC2 = 2.396284     1.8
#PC3 = 4.566041     2.85

# Reminder:
# GIF of 1=well calibrated, >1=liberal (too many small p-values), <1=conservative (too few small p-values)
# Note: GIFs > 2 indicate poor test calibration; try increasing K...and be skeptical!

## Estimate adjusted p-values
pvalues <- pv$pvalue 
zs <- pv$score

p.PC1 <- pchisq(zs[,1]^2/pv$gif[1], df = 1, lower = FALSE)  ## too conservative
adj.PC1 <- pchisq(zs[,1]^2/1.8, df = 1, lower = FALSE)
hist(p.PC1, main="Histogram of p-values PC1")
hist(adj.PC1, main="Histogram of adjusted p-values PC1")

qv_adj.PC1 <- which(qvalue::qvalue(adj.PC1, fdr=0.05)$signif)
length(qv_adj.PC1)

p.PC2 <- pchisq(zs[,2]^2/pv$gif[2], df = 1, lower = FALSE) 
adj.PC2 <- pchisq(zs[,2]^2/1.8, df = 1, lower = FALSE) 
hist(p.PC2, main="Histogram of p-values PC2")
hist(adj.PC2, main="Histogram of adjusted p-values PC2")

qv_adj.PC2 <- which(qvalue::qvalue(adj.PC2, fdr=0.05)$signif)
length(qv_adj.PC2)

p.PC3 <- pchisq(zs[,3]^2/pv$gif[3], df = 1, lower = FALSE) 
adj.PC3 <- pchisq(zs[,3]^2/2.85, df = 1, lower = FALSE) 
hist(p.PC3, main="Histogram of p-values PC3")
hist(adj.PC3, main="Histogram of adjusted p-values PC3")

qv_adj.PC3 <- which(qvalue::qvalue(adj.PC3, fdr=0.05)$signif)
length(qv_adj.PC3)


#Subset SNPs names:
snps_fil_lfmm_candidate1 <- Subset(snps, sites=qv_adj.PC1)
snps_fil_lfmm_candidate1@site_id ##These are the candidate SNPs associated with PC1. 

snps_fil_lfmm_candidate2 <- Subset(snps, sites=qv_adj.PC2)
snps_fil_lfmm_candidate2@site_id ##These are the candidate SNPs associated with PC2.

snps_fil_lfmm_candidate3 <- Subset(snps, sites=qv_adj.PC3)
snps_fil_lfmm_candidate3@site_id ##These are the candidate SNPs associated with PC3.


## Remove duplicated snps names (some SNPs might be associated with more than on PC)
snps_candidate <- unique(c(snps_fil_lfmm_candidate1@site_id, snps_fil_lfmm_candidate2@site_id, snps_fil_lfmm_candidate3@site_id))

length(snps_candidate) ## 204
snps_fil_lfmm_candidate_all <- Subset(snps, sites=snps_candidate)
snps_fil_lfmm_candidate_all@site_id ##These are all the candidate SNPs. From here to FASTA files and BLAST


###Save Results:
pdf("GPA/LFMM/GPA_Outliers_GIF_Original_PC1.pdf", onefile = T)
hist(p.PC1, main="Histogram of p-values PC1")
dev.off()

pdf("GPA/LFMM/GPA_Outliers_GIF_Calibrated_PC1.pdf", onefile = T)
hist(adj.PC1, main="Histogram of adjusted p-values PC1")
dev.off()

pdf("GPA/LFMM/GPA_Outliers_GIF_Original_PC2.pdf", onefile = T)
hist(p.PC2, main="Histogram of p-values PC2")
dev.off()

pdf("GPA/LFMM/GPA_Outliers_GIF_Calibrated_PC2.pdf", onefile = T)
hist(adj.PC2, main="Histogram of adjusted p-values PC2")
dev.off()

pdf("GPA/LFMM/GPA_Outliers_GIF_Original_PC3.pdf", onefile = T)
hist(p.PC3, main="Histogram of p-values PC3")
dev.off()

pdf("GPA/LFMM/GPA_Outliers_GIF_Calibrated_PC3.pdf", onefile = T)
hist(adj.PC2, main="Histogram of adjusted p-values PC3")
dev.off()

#### Save filtered vcf
Save(snps_fil_lfmm_candidate_all, "GPA/LFMM/snps_candidate_lfmm_GPA_pilocarpus.vcf")








#5. SAVE CANDIDATES LOCI AND THEIR FASTA---------------------------------------
## Retrieve chromosome IDS
CANDIDATES_d1 <- Chrom(snps_fil_lfmm_candidate1)
CH1 <- unique(CANDIDATES_d1[, 1])

CANDIDATES_d2 <- Chrom(snps_fil_lfmm_candidate2)
CH2 <- unique(CANDIDATES_d2[, 1])

CANDIDATES_d3 <- Chrom(snps_fil_lfmm_candidate3)
CH3 <- unique(CANDIDATES_d3[, 1])

contigs <- unique(c(CH1,CH2,CH3))
length(contigs) # 202

## Retrieve sequences from fasta file
fastafile <- read.fasta(file = "fasta/ref_pilocarpusredo.fasta", seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
head(fastafile)
length(fastafile)

SEQ1 <- fastafile[names(fastafile) %in% CH1]
SEQ2 <- fastafile[names(fastafile) %in% CH2]
SEQ3 <- fastafile[names(fastafile) %in% CH3]

total <- fastafile[names(fastafile) %in% contigs]

## Save new fasta files containing candidate sequences
write.fasta(sequences = SEQ1, names = names(SEQ1), nbchar = 150, file.out = "GPA/LFMM/SEQ1_lfmm_pilocarpus_GPA.fasta")
write.fasta(sequences = SEQ2, names = names(SEQ2), nbchar = 150, file.out = "GPA/LFMM/SEQ2_lfmm_pilocarpus_GPA.fasta")
write.fasta(sequences = SEQ3, names = names(SEQ3), nbchar = 150, file.out = "GPA/LFMM/SEQ3_lfmm_pilocarpus_GPA.fasta")
write.fasta(sequences = total, names = names(total), nbchar = 150, file.out = "GPA/LFMM/SEQ_total_lfmm_pilocarpus_GPA.fasta")

#END
