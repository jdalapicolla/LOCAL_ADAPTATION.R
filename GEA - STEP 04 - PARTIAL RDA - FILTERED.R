####################### VALE INSTITUTE OF TECHNOLOGY ##########################

########################      GEA-GPA TUTORIAL       #########################
#####################   STEP 04: PARTIAL RDA FOR GEA   #######################

### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###



# PRE-ANALYSIS-----------------------------------------------------------------

#AIMS:
#A. SELECTING SNPS UNDER SELECTION USING ENVIRONMENTAL VARIABLES
#B. EXTRACTING REFERENCE SEQUENCES AND LOCI THAT CONTAIN THESE LOCI 


##LOAD PACKAGES
library(vegan)
library(usdm)
library(ade4)
library(psych)
library(seqinr)
library(robust)  
library(qvalue)
library(r2vcftools)
library(LEA)
library(ggplot2)
library(ggrepel)


## LOAD FUNCTIONS TO BE USED ON THIS STEP.
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}







#1. CONVERT VCF TO GENO FORMAT FOR RDA AND LFMM ANALYSES-------------------------------
#Some analysis as RDA and LFMM2 does  not allow missing data, so we used LFMM to
#impute genetic missing data, based on the population that each individual
#belongs, using the package LEA

snps = vcfLink("vcf/pilocarpus_adap_filtered.vcf", overwriteID=T)
snps@meta
VCFsummary(snps) #264 individuals and 19025 SNPs.

gen = GenotypeMatrix(snps)
gen[gen==-1] = NA
sum(is.na(gen)) #643628
dim(gen) #264 19025

###Create LFMM file
#in the case we know the number of populations (k = 4)
snps_fil_lfmm = Lfmm(snps, "vcf/pilocarpus_filtered.lfmm")
project.snmf = snmf("vcf/pilocarpus_filtered.lfmm", K = 4, 
                    entropy = TRUE, repetitions = 10,
                    project = "new")

# select the run with the lowest cross-entropy value
best = which.min(LEA::cross.entropy(project.snmf, K = 4))

# Impute the missing genotypes
impute(project.snmf, "vcf/pilocarpus_filtered.lfmm", method = 'mode', K = 4, run = best)

# Convert lfmm to geno and save
lfmm2geno("vcf/pilocarpus_filtered.lfmm_imputed.lfmm", output.file = "vcf/pilocarpus_filtered_imputed.geno", force = TRUE)





#2. LOAD SNPs DATASET---------------------------------------------------------- 
## This script was developed by Brenna Forester and adapted to our data

## Load snps data set to take the snps and individuals ID
snps =  vcfLink("vcf/pilocarpus_adap_filtered.vcf", overwriteID=T)
snps@meta
VCFsummary(snps) ##264 individuals and 19025 SNPs.

## Read geno file
gen.imp <- read.geno("vcf/pilocarpus_filtered_imputed.geno")
gen.imp[1:10,1:10]
colnames(gen.imp) <- snps@site_id
rownames(gen.imp) <- snps@sample_id
dim(gen.imp) #264 19025

## Number of missing data
sum(is.na(gen.imp)) #0






#3. LOAD ENVIRONMENTAL DATASET------------------------------------------------
#loading extracted values per variable in Step 3 - Variables following Carvalho et al 2020:
GEA_variables = read.csv("GEA/Inputs/GEA_vars_filtered.csv", row.names = 1)
head(GEA_variables)

#Scale the variables:
pred_scale = scale(GEA_variables[,3:length(GEA_variables)])
#Renames cols and rows to graphics:
names(GEA_variables)
rownames(pred_scale) = GEA_variables[,1]
colnames(pred_scale) = c("H+AL", "SB", "ARIDITY", "BIO_02", "BIO_05")
head(pred_scale)
class(pred_scale)





#4. LOAD QMATRIX BASED ON sNMF------------------------------------------------- 
# Load files
pop = read.csv("data_raw/pop_LEA_DAPC_TESS_pilocarpus.csv", row.names = 1)
metafiles = as.data.frame(snps@meta)

# selecting the individuals:
pop_filtered = pop[pop$ind_ID %in% metafiles$ind_ID,]

# Check if there are some difference between files. If they are corrected, function will return "character (0)"
setdiff(as.character(pop_filtered$ind_ID), as.character(GEA_variables$id))
setdiff(as.character(GEA_variables$id), as.character(pop_filtered$ind_ID))
# Check if identical samples were selected in the same order. If they are corrected, function will return "TRUE"
identical(as.character(pop_filtered$ind_ID), as.character(GEA_variables$id))

#saving the Qmatrix of sNMF to applied in the RDA
pilo_qmatrix = as.matrix(pop_filtered[, 8:11])
head(pilo_qmatrix)







#5. PERFORM PARTIAL RDA-------------------------------------------------------- 
#If you have more than one population, use the population assigned for each individual 
#in the Condition argument of RDA to control for population structure
# m.rda = rda(gen.imp ~ pred_scale + Condition(as.factor(pop_lea))) 

## RDA with mahalanobis distance
m.rda = rda(gen.imp ~ pred_scale + Condition(pilo_qmatrix))
m.rda

options(scipen = 999)
control_vif = as.data.frame(vif.cca(m.rda)) #at least VIF<10  (PRED)
control_vif
#pred_scaleBIO_05 46.898596

#remover
head(pred_scale)
pred_scale = pred_scale[, -5]
head(pred_scale)


## RDA with mahalanobis distance - 2nd round
m.rda = rda(gen.imp ~ pred_scale + Condition(pilo_qmatrix))
m.rda

control_vif = as.data.frame(vif.cca(m.rda)) #at least <10 (PRED)
control_vif
#pred_scaleARIDITY =  5.470459

#Save the VIF results
write.csv(control_vif, file="GEA/partial_RDA/control_vif_partialrda.csv")


#Quality
RsquareAdj(m.rda)
#$r.squared 0.03639314
#$adj.r.squared 0.02537703


#Axis
plot(m.rda, scaling=3)
screeplot(m.rda)

summary(eigenvals(m.rda, model= "constrained"))
#for graphics:
AX1 = "RDA1 (44.79%)"
AX2 = "RDA2 (30.13%)"
AX3 = "RDA3 (16.51%)"
AX4 = "RDA4 (8.57%)"


#check our RDA model for significance using formal tests.
signif.full = anova.cca(m.rda, permutations = 100)
signif.full
#           Df Variance      F   Pr(>F)   
# Model      4    196.0 3.1865 0.009901 **
# Residual 255   3921.4  

signif.terms = anova.cca(m.rda, by="terms", permutations = 100)
signif.terms
#           Df Variance      F   Pr(>F)   
# Model      4    196.0 3.1865 0.009901 **
# Residual 255   3921.4

## Check axis significance
axis_sig = anova.cca(m.rda, by= "axis", permutations = 100) 
axis_sig ## 3 Axes!!
#            Df Variance      F   Pr(>F)   
#  RDA1       1     87.8 5.7089 0.009901 **
#  RDA2       1     59.1 3.8404 0.009901 **
#  RDA3       1     32.4 2.1042 0.009901 **
#  RDA4       1     16.8 1.0927 0.069307 . 
#  Residual 255   3921.4

## we chose the first 2 axes, which explained 74.92% of the variation

#Save RDA results:
sum_rda = summary(m.rda)
write.csv(as.data.frame(sum_rda$species), file="GEA/partial_RDA/RDA_values_SNPs.csv") #snps
write.csv(as.data.frame(sum_rda$sites), file="GEA/partial_RDA/RDA_values_individuals.csv") #individuals
write.csv(as.data.frame(sum_rda$biplot), file="GEA/partial_RDA/RDA_values_variables.csv") #variables






#6. OUTLIER DETECTION --------------------------------------------------------- 
## Estimate adjusted p-values based on gif
## Save loading of the first and second pcs
load.rda = sum_rda$species[,1:2]
K = 2  # the number of RDA axes you're looking at
zscale = apply(load.rda, 2, scale)        # I'm not sure scaling is the best approach here
mscale = covRob(zscale, distance=TRUE, na.action=na.omit, estim="pairwiseGK")$dist
      (gif = median(mscale)/qchisq(0.5, df=K)) #1.577679

## Look at histogram
rda.pval <- pchisq(mscale/gif, df=K, lower.tail=FALSE)   # Remember: you can always change the GIF!!
hist(rda.pval)

## Try different gif values if necessary
rda.pval2 <- pchisq(mscale/1.1, df=K, lower.tail=FALSE)  # Less conservative GIF (smaller)
hist(rda.pval2)

## Select snps with FDR < 0.05
rda.qval <- qvalue(rda.pval2)$qvalues
m.snps <- (rda.FDR <- colnames(gen.imp)[which(rda.qval < 0.05)])
length(m.snps)

#Results for Calibrated GIF:
  # GIF     SNPs
  # 1.0     3,267
  # 1.1     2,816 ***
  # 1.2     2,432
  # 1.3     2,095
  # 1.4     1,770
  # 1.5     1,513
  # ORI     1,309


## Save the snps names to compare with other approaches
pdf("GEA/partial_RDA/GEA_Outliers_GIF_Original.pdf", onefile = T)
hist(rda.pval)
dev.off()

pdf("GEA/partial_RDA/GEA_Outliers_GIF_Calibrated.pdf", onefile = T)
hist(rda.pval2)
dev.off()

write.csv(m.snps, "GEA/partial_RDA/Pilocarpus_snps_candidates_GEA_partialRDA.csv")







#7. GRAPHICS FOR PARTIAL RDA - CANDIDATES-------------------------------
## Create Graphic objects:
# Let's look at where these candidates are in the ordination space:
snp.color = as.data.frame(colnames(gen.imp), stringsAsFactors=F)
snp.color[,2] = 'non-candidate'
head(snp.color)

#positions for candidates SNPs
positions = which(snp.color$`colnames(gen.imp)` %in% rda.FDR == TRUE)
snp.color$V2[positions] = 'candidate'

#Create objects graphs:

sp = as.data.frame(sum_rda$species[,1:2]) #Selecting 2 axis for RDA #Depending on the drawing result, the drawing data can be enlarged or reduced to a certain extent, as follows
head(sp) #snps

st = as.data.frame(sum_rda$sites[,1:2])
head(st) #individuals

yz = as.data.frame(sum_rda$biplot[,1:2])
head(yz)
row.names(yz) = c("H+AL", "SB", "ARIDITY", "BIO02") #Correcting variable names
head(yz) #biplot

grp = as.data.frame(snp.color$V2) #Grouping snps by outliers snps
colnames(grp)="group" #
head(grp)
colors = c('red', 'white') # set colors
legend = c("Candidate", "Non-candidate") #legend itens - alphabetichal order

#RDA results. It will be added as axes labels:
AX1
AX2


## Plot Axes 1 & 2 for all candidates
pdf("GEA/partial_RDA/GEA_partialRDA_canditates_all.pdf", onefile = T)
ggplot() +
  geom_point(data = sp,
             mapping = aes(x= RDA1, y = RDA2, fill=factor(grp$group)),
             size=2.5, shape=21)+
  scale_fill_manual(values =colors, labels = legend)+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.5,"cm"),
                             type = "closed"),linetype=1, size=0.9,colour = "black")+
  geom_label_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x=AX1, 
       y=AX2)+
  guides(fill=guide_legend(title="SNPs", override.aes = list(size=5)))+
  theme_bw()
dev.off()





#8. CORRELATION AMONG CANDIDATE SNPS AND PREDICTORS --------------------------- 
## Add in the correlations of each candidate SNP with the environmental predictors:
foo = matrix(nrow=length(m.snps), ncol=4)  # columns for all predictors
colnames(foo) = c("H+AL", "SB","ARIDITY", "BIO02") #order following pred_scale object

for (i in 1:length(m.snps)) {
  nam = m.snps[i]
  snp.gen = gen.imp[,nam]
  foo[i,] = apply(pred_scale,2,function(x) cor(x,snp.gen, method="spearman"))
}

cand = cbind.data.frame(m.snps,foo)  
head(cand)

## Looking for duplicate detections
length(cand$m.snps[duplicated(cand$m.snps)])
cand <- cand[!duplicated(cand$m.snps),] # remove duplicate detections
head(cand)

## See which of the predictors each candidate SNP is most strongly correlated with. Edit cols number according to the number of variables
for (i in 1:length(cand$m.snps)) {
  bar = cand[i,]
  cand[i,6] = names(which.max(abs(bar[2:5]))) # gives the variable
  cand[i,7] = max(abs(bar[2:5]))              # gives the absolute correlation value
  cand[i,8] = ifelse (abs(max((bar[2:5])))>abs(min((bar[2:5]))),max((bar[2:5])), min((bar[2:5]))) # gives the raw correlation value
  
}

colnames(cand)[6] <- "predictor"
colnames(cand)[7] <- "correlation"
colnames(cand)[8] <- "correlation_real"


#Calculate Adaptative Scores:
cand_predictors = as.data.frame(table(cand$predictor))
cand_predictors$Adapt_Scores = cand_predictors$Freq/length(cand$m.snps)
cand_predictors

#Save the cor results:
write.csv(cand, file="GEA/partial_RDA/Correlation_candidates_snps_partialrda_all.csv")
write.csv(cand_predictors, file="GEA/partial_RDA/Adapt_Scores_candidates_snps_partialrda_all.csv")


###Selecting only SNPs Candidate with high r² values:
## Analyzing correlation values
min(cand$correlation) #0.03
max(cand$correlation) #0.63

#by decil:
Min_Correlation = c("All", "0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")
subsets_cand = length(cand$m.snps)
subsets_cand[2] = length(which(cand$correlation > 0.1)) #2729
subsets_cand[3] = length(which(cand$correlation > 0.2)) #1923
subsets_cand[4] = length(which(cand$correlation > 0.3)) #901
subsets_cand[5] = length(which(cand$correlation > 0.4)) #300
subsets_cand[6] = length(which(cand$correlation > 0.5)) #72
subsets_cand[7] = length(which(cand$correlation > 0.6)) #14
subsets_cand[8] = length(which(cand$correlation > 0.7)) #0
subsets_cand[9] = length(which(cand$correlation > 0.8)) #0
subsets_cand[10] = length(which(cand$correlation > 0.9)) #0

subsets_cand = as.data.frame(cbind(Min_Correlation, subsets_cand))
subsets_cand



#Selected SNPs with more than 40% of correlation
cand2 = cand[which(cand$correlation > 0.4), ]

#Calculate Adaptative Scores:
cand_predictors2 = as.data.frame(table(cand2$predictor))
cand_predictors2$Adapt_Scores = cand_predictors2$Freq/length(cand2$m.snps)
cand_predictors2

#Save the cor results:
write.csv(cand2, file="GEA/partial_RDA/Correlation_candidates_snps_partialrda_subset.csv")
write.csv(subsets_cand, file="GEA/partial_RDA/Subsets_candidates_snps_partialrda.csv")
write.csv(cand_predictors2, file="GEA/partial_RDA/Adapt_Scores_candidates_snps_partialrda_subset.csv")


#Selected SNPs from subset by variables:

cand2_ari = cand2[which(cand2$predictor == "ARIDITY"), ]
cand2_bio2 = cand2[which(cand2$predictor == "BIO02"), ]
cand2_h_al = cand2[which(cand2$predictor == "H+AL"), ]

#Save the cor results:
write.csv(cand2_ari, file="GEA/partial_RDA/Correlation_candidates_snps_partialrda_aridity.csv")
write.csv(cand2_bio2, file="GEA/partial_RDA/Correlation_candidates_snps_partialrda_bio2.csv")
write.csv(cand2_h_al, file="GEA/partial_RDA/Correlation_candidates_snps_partialrda_h_al.csv")



#9. GRAPHICS FOR PARTIAL RDA - PREDICTORS------------------------------------------
#Create objects graphs:
## Creating Groupiing variables per SNPs and Variable
sel = cand2$m.snps
env = cand2$predictor

## Color by predictor:
col.pred = rownames(m.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) { # color code candidate SNPs
  foo = match(sel[i],col.pred)
  col.pred[foo] = env[i]
}

col.pred[grep("snp",col.pred)] = 'non-candidate' # non-candidate SNPs
col.pred[1:20]
length(col.pred)

#Grouping variable
grp=as.data.frame(col.pred) #Grouping snps by variable
colnames(grp)="group"
colors = c('yellow', 'red', 'blue', 'white')# set colors
legend = c("ARIDITY", "BIO02", "H+AL", "Non-candidate") #legend itens - alphabetichal order, no candidate fro SB


# Axes 1 & 2
pdf("GEA/partial_RDA/GEA_partialRDA_preditors_all_subset.pdf", onefile = T)
ggplot() +
  geom_point(data = sp,
             mapping = aes(x= RDA1, y = RDA2, fill=factor(grp$group)),
             size=2.5, shape=21)+
  scale_fill_manual(values =colors, labels = legend)+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.5,"cm"),
                             type = "closed"),linetype=1, size=0.9,colour = "black")+
  geom_label_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x=AX1, 
       y=AX2)+
  guides(fill=guide_legend(title="Variables", override.aes = list(size=5)))+
  theme_bw()
dev.off()


# Let's look at where these subset candidates are in the ordination space:
snp.color = as.data.frame(colnames(gen.imp), stringsAsFactors=F)
snp.color[,2] = 'non-candidate'
head(snp.color)

#positions for candidates SNPs
positions = which(snp.color$`colnames(gen.imp)` %in% cand2$m.snps == TRUE)
snp.color$V2[positions] = 'candidate'

grp = as.data.frame(snp.color$V2) #Grouping snps by outliers snps
colnames(grp)="group_subset" #
head(grp)
colors = c('red', 'white') # set colors
legend = c("Candidate", "Non-candidate") #legend itens - alphabetichal order


## Plot Axes 1 & 2
pdf("GEA/partial_RDA/GEA_partialRDA_canditates_subset.pdf", onefile = T)
ggplot() +
  geom_point(data = sp,
             mapping = aes(x= RDA1, y = RDA2, fill=factor(grp$group)),
             size=2.5, shape=21)+
  scale_fill_manual(values =colors, labels = legend)+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.5,"cm"),
                             type = "closed"),linetype=1, size=0.9,colour = "black")+
  geom_label_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x=AX1, 
       y=AX2)+
  guides(fill=guide_legend(title="SNPs", override.aes = list(size=5)))+
  theme_bw()
dev.off()







#10. SAVE CANDIDATES LOCI AND THEIR FASTA--------------------------------------
# Load snps data set 
snps =  vcfLink("vcf/pilocarpus_adap_filtered.vcf", overwriteID=T)
snps@meta
VCFsummary(snps) ## 264 individuals and 19025 SNPs.

sel = cand2$m.snps ##These are all the candidate SNPs. From here to FASTA files and BLAST
sel

snps_fil_rda_candidate <- Subset(snps, sites=as.character(sel))
snps_fil_rda_candidate@site_id ##These are all the candidate SNPs. From here to FASTA files and BLAST

## Save filtered vcf
Save(snps_fil_rda_candidate, "GEA/partial_RDA/snps_candidate_gea_RDA_mahalanobis_pilocarpus_subset.vcf")

## Retrieve chromosome IDS
snp_chrom <- Chrom(snps_fil_rda_candidate)
CH1 <- unique(snp_chrom[, 1]) #293

## Retrieve sequences from fasta file 
fastafile <- read.fasta(file = "fasta/ref_pilocarpusredo.fasta", seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
head(fastafile)
length(fastafile)

SEQ1 <- fastafile[names(fastafile) %in% CH1]
length(SEQ1) #293

## Save new fasta files containing candidate sequences
seqinr::write.fasta(sequences = SEQ1, names = names(SEQ1), nbchar = 150, file.out = "GEA/partial_RDA/SEQ_gea_RDA_mahalanobis_pilocarpus_subset.fasta")

##END