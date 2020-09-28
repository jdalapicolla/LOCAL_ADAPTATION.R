####################### VALE INSTITUTE OF TECHNOLOGY ##########################

########################      GEA-GPA TUTORIAL       ##########################
######################### STEP 06: COMPARING RESULTS ##########################

### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###



#PRE-ANALYSIS ------------------------------------------------------------------


#AIMS:
#A. COMPARING CANDIDATES SNPS IDENTIFIED WITH RDA AND LFMM2


##LOAD THE PACKAGES
library(LEA)
library(r2vcftools)
library(vegan)    
library(usdm)
library(spdep)
library(adespatial)
library(spacemakeR) 
library(ade4)
library(psych)
library(seqinr)
library(robust) 
library(ggplot2)
library(ggrepel)
library(VennDiagram)


##5. LOAD FUNCTIONS TO BE USED ON THIS STEP.
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}



#1. LOAD GENOMIC DATASET-------------------------------------------------------
#Load snps data set to take the snps and individuals ID
snps =  vcfLink("vcf/pilocarpus_adap_filtered.vcf", overwriteID=T)
snps@meta
VCFsummary(snps) ## 264 individuals and 19025 SNPs.

## Read geno file
gen.imp = read.geno("vcf/pilocarpus_filtered_imputed.geno")
gen.imp[1:10,1:10]
colnames(gen.imp) = snps@site_id
rownames(gen.imp) = snps@sample_id
dim(gen.imp)

## Number of missing data
sum(is.na(gen.imp)) #0









#2. LOAD GPA RESULTS-----------------------------------------------------------
## Load RDA and LFMM filtered vcf - Steps 4 and 5:

rda_snps = vcfLink("GPA/partial_RDA/snps_candidate_gea_RDA_mahalanobis_pilocarpus_subset.vcf", overwriteID=F) # here overwrite need to be F, otherwise we will loss the snps name

lfmm_snps = vcfLink("GPA/LFMM/snps_candidate_lfmm_GPA_pilocarpus.vcf", overwriteID=F)

all_candidate = unique(c(rda_snps@site_id,lfmm_snps@site_id))
all_candidate_gen = gen.imp[,all_candidate]
dim(all_candidate_gen) #335

inter_candidate = intersect(c(rda_snps@site_id), c(lfmm_snps@site_id))
inter_candidate_gen = gen.imp[,inter_candidate]
dim(inter_candidate_gen) #2











#3. VENN DIAGRAM--------------------------------------------------------------- 
#select SNPs names:
RDA = Chrom(rda_snps)[,1]
LFMM2 = Chrom(lfmm_snps)[,1]

## Creating a Venn Diagram and save it as pdf
pdf("GPA/VennDiagram/VennDiagram_SNPscandidates_GPA.pdf", onefile = T)
plot.new()
grid.draw(w1 <- venn.diagram(list(RDA=RDA, LFMM2=LFMM2),
                             lty = c("blank", "blank"),
                             fill = c("red", "blue"),
                             alpha = c(0.5, 0.5), cat.cex = 1.2, cex= 1.5,  cat.pos = 0, 
                             filename=NULL ))

dev.off()








#4. LOAD ENVIRONMENTAL DATASET------------------------------------------------- 
#loading extracted values per variable in Step 3 - Variables following Carvalho et al 2020:
GEA_variables = read.csv("GPA/Inputs/GPA_vars_manual.csv", row.names = 1)
head(GEA_variables)

#selecting only variables used in RDA to comparison
#names(GEA_variables)
#GEA_variables = GEA_variables[, c(1:6)]

#checking VIF:
vif(GEA_variables[, -c(1:2)]) #VIF  <2

#Scale the variables:
pred_scale = scale(GEA_variables[,3:10])
#Renames cols and rows to graphics:
names(GEA_variables)
rownames(pred_scale) = GEA_variables[,1]
colnames(pred_scale) = c("N","P", "K", "Ca", "Mg", "Cu", "Zn", "Pilo")
head(pred_scale)
class(pred_scale)









#5. POPULATION QMATRIX DATASET-------------------------------------------------
# Load files
pop = read.csv("data_raw/pop_LEA_DAPC_TESS_pilocarpus.csv", row.names = 1)
metafiles = as.data.frame(snps@meta)

# selecting the individuals:
pop_filtered = pop[pop$ind_ID %in% metafiles$ind_ID,]

#E. Check if there are some difference between files. If they are corrected, function will return "character (0)"
setdiff(as.character(pop_filtered$ind_ID), as.character(GEA_variables$id))
setdiff(as.character(GEA_variables$id), as.character(pop_filtered$ind_ID))
#F. Check if identical samples were selected in the same order. If they are corrected, function will return "TRUE"
identical(as.character(pop_filtered$ind_ID), as.character(GEA_variables$id))


#saving the Qmatrix of sNMF to applied in the RDA
pilo_qmatrix = as.matrix(pop_filtered[, 8:11])
head(pilo_qmatrix)












#6. PARTIAL RDA ANALYSES - ALL CANDIDATES-------------------------------------- 
## RDA with mahalanobis distance
m.rda = rda(all_candidate_gen ~ pred_scale + Condition(pilo_qmatrix)) #change for intersected candidates
m.rda

options(scipen = 999)
control_vif = as.data.frame(vif.cca(m.rda)) #at least <10  (PRED)
control_vif

#Save the VIF results
write.csv(control_vif, file="GPA/RDA2/control_vif_rda2.csv")


#Quality
RsquareAdj(m.rda) 
#r.squared 0.04444495
#adj.r.squared  0.02640762

#Axis
plot(m.rda, scaling=3)
screeplot(m.rda)

summary(eigenvals(m.rda, model= "constrained"))
#for graphics:
AX1 = "RDA1 (32.70%)"
AX2 = "RDA2 (18.15%)"
AX3 = "RDA3 (14.51%)"
AX4 = "RDA4 (10.22%)"
AX5 = "RDA5 (7.11%)"
AX6 = "RDA6 (6.22%)"
AX7 = "RDA7 (5.81%)"
AX8 = "RDA8 (4.87%)"


#check our RDA model for significance using formal tests.
signif.full = anova.cca(m.rda, permutations = 1000)
signif.full
#            Df Variance     F   Pr(>F)   
#  Model      8    6.194 2.4104 0.000999 ***
#  Residual 251   80.629  

signif.terms = anova.cca(m.rda, by="terms", permutations = 1000)
signif.terms
#            Df Variance     F   Pr(>F)   
#  Model      8    6.194 2.4104 0.000999 ***
#  Residual 251   80.629  


## Check axis significance
axis_sig = anova.cca(m.rda, by= "axis", permutations = 1000) 
axis_sig ## 4 Axes!!
#            Df Variance       F   Pr(>F)    
#  RDA1       1    2.025 6.3048 0.000999 ***
#  RDA2       1    1.124 3.5003 0.000999 ***
#  RDA3       1    0.899 2.7984 0.000999 ***
#  RDA4       1    0.633 1.9710 0.000999 ***
#  RDA5       1    0.440 1.3701 0.123876    
#  RDA6       1    0.410 1.2769 0.159840    
#  RDA7       1    0.360 1.1212 0.388611    
#  RDA8       1    0.302 0.9403 0.669331    
#  Residual 251   80.629 

## we chose the first 4 axes, which explained 75.58% of the variation

#Save RDA results:
sum_rda = summary(m.rda)
write.csv(as.data.frame(sum_rda$species), file="GPA/RDA2/RDA2_values_SNPs.csv") #snps
write.csv(as.data.frame(sum_rda$sites), file="GPA/RDA2/RDA2_values_individuals.csv") #individuals
write.csv(as.data.frame(sum_rda$biplot), file="GPA/RDA2/RDA2_values_variables.csv") #variables






#7. PLOT RESULTS - ALL CANDIDATES ------------------------------------------------------ 
#Create objects graphs:

sp = as.data.frame(sum_rda$species[,1:4])*2 #Selecting 5 axis for RDA #Depending on the drawing result, the drawing data can be enlarged or reduced to a certain extent, as follows
head(sp) #snps

st = as.data.frame(sum_rda$sites[,1:4])
head(st) #individuals

yz = as.data.frame(sum_rda$biplot[,1:4])*5
head(yz)
row.names(yz) = c("N","P", "K", "Ca", "Mg", "Cu", "Zn", "Pilocarpine") #Correcting variable names
head(yz) #biplot

grp = as.data.frame(pop_filtered$PopID_snmf) #Grouping snps by variable/selecion
colnames(grp)="group" #
head(grp)
colors = c('red', 'blue', 'yellow', "pink") #
legend = c("POP1", "POP2", "POP3", "POP4") #Alphabetichal

#RDA results. It will be added as axes labels:
AX1
AX2
AX3
AX4

## Axes 1 TO 4
pdf("GPA/RDA2/GEA_RDA_IND_ALL_1_4.pdf", onefile = T)
ggplot() +
  geom_point(data = sp,
             mapping = aes(x= RDA1, y = RDA4),
             size=1.5, shape=21, fill ="grey", color="grey")+
  geom_point(data = st,
             mapping = aes(x= RDA1, y = RDA4, fill=factor(grp$group)),
            size=2, shape=21)+
  scale_fill_manual(values =colors, labels = legend)+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA4), 
               arrow = arrow(angle=22.5,length = unit(0.5,"cm"),
                             type = "closed"),linetype=1, size=0.9,colour = "black")+
  geom_label_repel(data = yz,aes(RDA1,RDA4,label=row.names(yz)))+
  labs(x=AX1, 
       y=AX4)+
  guides(fill=guide_legend(title="Populations",override.aes = list(size=5)))+
  theme_bw()
dev.off()

#Checking localities of clustering individuals in RDA space:
indi = as.data.frame(sum_rda$sites[,1:4])

#RDA1 and RDA2 <-2
name_clusters = which(indi$RDA2 < -2)

local_label = snps@meta$local_ID[name_clusters]
local_label


## Axes 1 & 2 with localities labels
pdf("GPA/RDA2/GPA_RDA_IND_ALL_labels_1_2.pdf", onefile = T)
ggplot() +
  geom_point(data = sp,
             mapping = aes(x= RDA1, y = RDA2),
             size=1.5, shape=21, fill ="grey", color="grey")+
  geom_point(data = st,
             mapping = aes(x= RDA1, y = RDA2, fill=factor(grp$group)),
             size=2, shape=21)+
  scale_fill_manual(values =colors, labels = legend)+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.5,"cm"),
                             type = "closed"),linetype=1, size=0.9,colour = "black")+
  geom_label_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  geom_label_repel(data = st[name_clusters,],aes(RDA1,RDA2,label=local_label))+
  labs(x=AX1, 
       y=AX2)+
  guides(fill=guide_legend(title="Populations",override.aes = list(size=5)))+
  theme_bw()
dev.off()



## END
