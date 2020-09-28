#This code serves to estimated polygenic scores of candidate loci and model them with 
#environmental variables.
#First of all, we will estimate polygenic scores. For that we will analyse candidate loci detected from RDA approach. 
#We will assess the correlation between each candidate loci (0,1,2 values) with the environmental
#variables that they were most correlated. If this correlation is negative, we will invert the genotype values (0,1,2 to 2,1,0).
#Them we will sum the genotype values for each individual and it represent the polygenic scores. 
#Finally we will model polygenic scores in relation to environmental variable.

rm(list=ls())

# Load packages and functions -----------------------------------------------------------

library(r2vcftools)
library(ggplot2)
library(MuMIn)

### Summary stats
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}


# Load candidate snps -----------------------------------------------------

##Here we will load candidate loci selected from RDA approach. 
#Them we will separate the snps that were correlated with each enviromental variables, 
#and then insert genotype values for those that present negative correlation.

## load rda candidate loci

rda_snps = vcfLink("GPA/partial_RDA/snps_candidate_gea_RDA_mahalanobis_pilocarpus_subset.vcf", overwriteID=F)
VCFsummary(rda_snps) #264 individuals and  133 SNPs.

# load dataframe containing correlation value of each snps with environmental variables.
# colunm correlation- absolute value of correlation between the snps and predictor that is most correlated with such snps
# colunm correlation_real- value of correlation between the snps and predictor that is most correlated with such snps

# this dataframe was built in RDA_select_candidate_snps.R

cor_data <- read.csv("GPA/partial_RDA/Correlation_candidates_snps_partialrda_subset.csv", row.names = 1)
head(cor_data)

## Then we will select the snps with negative correlations
cor_neg <- subset(cor_data, correlation_real < 0)
snps_cor_neg <- as.vector(cor_neg$m.snps)
snps_neg_vcf <- Subset(rda_snps, sites=snps_cor_neg)
VCFsummary(snps_neg_vcf) #264 individuals and 64 SNPs.
snps_neg_vcf@site_id

## And then we invert the genoype values
snps_neg_gen <- GenotypeMatrix(snps_neg_vcf)
dim(snps_neg_gen)
snps_neg_gen[1:10, 1:10]
summary(snps_neg_gen)

# With negative correlation we will invert the genotype values:
#(0,1,2 to 2,1,0).
#0 = homo for 1st allele
#1 = hetero
#2 = homo for second allele

for (i in 1:ncol(snps_neg_gen)){
  snps_neg_gen[,i][which(snps_neg_gen[,i]==0)] = 2
  snps_neg_gen[,i][which(snps_neg_gen[,i]==2)] = 0
  snps_neg_gen[,i][which(snps_neg_gen[,i]==-1)] = NA
}

## Finally we will combine the genotype of snps with positive correlation with environmental variables
#and the new snps dataset with genotype inverted.

cor_pos <- subset(cor_data, correlation_real > 0)
snps_cor_pos <- as.vector(cor_pos$m.snps)
snps_pos_vcf <- Subset(rda_snps, sites=snps_cor_pos)
VCFsummary(snps_pos_vcf) #264 individuals and 69 SNPs.
snps_pos_vcf@site_id

snps_pos_gen <- GenotypeMatrix(snps_pos_vcf)
for (i in 1:ncol(snps_pos_gen)){
  snps_pos_gen[,i][which(snps_pos_gen[,i]==-1)]<-NA
}
dim(snps_pos_gen)
dim(snps_neg_gen)

snps_recorde <- cbind(snps_neg_gen, snps_pos_gen)
dim(snps_recorde)
head(snps_recorde)

## Now we will select the name of the snps that are most related to each variable

head(cor_data)

snps_n <- as.vector(subset(cor_data, predictor=="N")$m.snp)
snps_p <- as.vector(subset(cor_data, predictor=="P")$m.snp) 
snps_k <- as.vector(subset(cor_data, predictor=="K")$m.snp)
snps_ca <- as.vector(subset(cor_data, predictor=="Ca")$m.snp) #0 snps so we removed it!
snps_mg <- as.vector(subset(cor_data, predictor=="Mg")$m.snp)
snps_cu <- as.vector(subset(cor_data, predictor=="Cu")$m.snp)
snps_zn <- as.vector(subset(cor_data, predictor=="Zn")$m.snp) #0 snps so we removed it!
snps_pilo <- as.vector(subset(cor_data, predictor=="Pilocarpine")$m.snp)


## Then we will select the snps genotypes that are most related to each variable

snps_n_gen <- snps_recorde[,snps_n]
snps_p_gen <- snps_recorde[,snps_p]
snps_k_gen <- snps_recorde[,snps_k]
snps_mg_gen <- snps_recorde[,snps_mg]
snps_cu_gen <- snps_recorde[,snps_cu]
snps_pilo_gen <- snps_recorde[,snps_pilo]

###Now we will calculate the polygenic score by summing across the rows

polygenic_score <- as.data.frame(rowSums(snps_recorde, na.rm=T))
colnames (polygenic_score) <- "SCORE"

polygenic_score_n <- as.data.frame(rowSums(snps_n_gen, na.rm=T))
colnames (polygenic_score_n) <- "SCORE"

polygenic_score_p <- as.data.frame(rowSums(snps_p_gen, na.rm=T))
colnames (polygenic_score_p) <- "SCORE"

polygenic_score_k <- as.data.frame(rowSums(snps_k_gen, na.rm=T))
colnames (polygenic_score_k) <- "SCORE"

polygenic_score_mg <- as.data.frame(rowSums(snps_mg_gen, na.rm=T))
colnames (polygenic_score_mg) <- "SCORE"

polygenic_score_cu <- as.data.frame(rowSums(snps_cu_gen, na.rm=T))
colnames (polygenic_score_cu) <- "SCORE"

polygenic_score_pilo <- as.data.frame(rowSums(snps_pilo_gen, na.rm=T))
colnames (polygenic_score_pilo) <- "SCORE"

#Save Results:
polygenic_score_total = cbind(polygenic_score, polygenic_score_n, polygenic_score_p, polygenic_score_k,polygenic_score_mg, polygenic_score_cu, polygenic_score_pilo)
colnames(polygenic_score_total) = c("TOTAL", "N", "P", "K", "Mg", "Cu", "Pilocarpine")
head(polygenic_score_total)

write.csv(polygenic_score_total, file="PolygenicScores/Phenotype/Polygenic_scores_total.csv")



### Now we will load the environmental variables

env =read.csv("GPA/Inputs/GPA_vars_manual.csv", row.names = 1)
head(env)
str(env)

pred <- env[,-c(1:2,6,9)]
rownames(pred) <- env$id
head(pred)

### Test relationship between polygenic scores and predictors variables using linear and quadract models.

#######N
mod1_n <- lm(polygenic_score_n$SCORE ~ pred$N_GPA) #linear
mod2_n <- lm(polygenic_score_n$SCORE ~ pred$N_GPA + I(pred$N_GPA^2))#quadratic
mod3_n <- lm(polygenic_score_n$SCORE ~ 1)#Null

model_sel_n = model.sel(mod1_n,mod2_n,mod3_n)
model_sel_n #model2
summary(mod2_n)
bestmodel_n = capture.output(summary(mod2_n))

#Save results:
write.csv(model_sel_n, file="PolygenicScores/Phenotype/Model_Selection_N.csv")
write.csv(bestmodel_n, file="PolygenicScores/Phenotype/Best_Model_N.csv")


prd <- data.frame(hp=seq(from=range(pred$N_GPA)[1], to=range(pred$N_GPA)[2], length.out = nrow(pred)))
err <- predict(mod2_n, newdata=prd, se.fit=T)
prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

pdf("PolygenicScores/Phenotype/Graphic_N.pdf", onefile = T)
ggplot(prd,aes(x=pred$N_GPA, y=fit))+
  labs(x= "N", y="Polygenic Score")+
  geom_point(aes(x=pred$N_GPA, y=polygenic_score_n$SCORE), alpha=0.5) +
  geom_smooth(method = "loess", color ="red", size=1.5)+
  geom_smooth(aes(y = lci), color = "red", linetype = 2 , size=0.5) +
  geom_smooth(aes(y = uci), color = "red", linetype = 2 , size=0.5) +
  theme_bw() +
  theme(axis.title.y = element_text(size=12, color = "black", face = "bold"),
        axis.title.x = element_text(size=12, color = "black", face = "bold"))
dev.off()




#######P
mod1_p <- lm(polygenic_score_p$SCORE ~ pred$P_GPA) #linear
mod2_p <- lm(polygenic_score_p$SCORE ~ pred$P_GPA + I(pred$P_GPA^2))#quadratic
mod3_p <- lm(polygenic_score_p$SCORE ~ 1)#Null

model_sel_p = model.sel(mod1_p,mod2_p,mod3_p)
model_sel_p #model2
summary(mod2_p)
bestmodel_p = capture.output(summary(mod2_p))

#Save results:
write.csv(model_sel_p, file="PolygenicScores/Phenotype/Model_Selection_P.csv")
write.csv(bestmodel_p, file="PolygenicScores/Phenotype/Best_Model_P.csv")


prd <- data.frame(hp=seq(from=range(pred$P_GPA)[1], to=range(pred$P_GPA)[2], length.out = nrow(pred)))
err <- predict(mod2_p, newdata=prd, se.fit=T)
prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

pdf("PolygenicScores/Phenotype/Graphic_P.pdf", onefile = T)
ggplot(prd,aes(x=pred$P_GPA, y=fit))+
  labs(x= "P", y="Polygenic Score")+
  geom_point(aes(x=pred$P_GPA, y=polygenic_score_p$SCORE), alpha=0.5) +
  geom_smooth(method = "loess", color ="red", size=1.5)+
  geom_smooth(aes(y = lci), color = "red", linetype = 2 , size=0.5) +
  geom_smooth(aes(y = uci), color = "red", linetype = 2 , size=0.5) +
  theme_bw() +
  theme(axis.title.y = element_text(size=12, color = "black", face = "bold"),
        axis.title.x = element_text(size=12, color = "black", face = "bold"))
dev.off()


#######K
mod1_k <- lm(polygenic_score_k$SCORE ~ pred$K_GPA) #linear
mod2_k <- lm(polygenic_score_k$SCORE ~ pred$K_GPA + I(pred$K_GPA^2))#quadratic
mod3_k <- lm(polygenic_score_k$SCORE ~ 1)#Null

model_sel_k = model.sel(mod1_k,mod2_k,mod3_k)
model_sel_k #model2
summary(mod2_k)
bestmodel_k = capture.output(summary(mod2_k))

#Save results:
write.csv(model_sel_k, file="PolygenicScores/Phenotype/Model_Selection_K.csv")
write.csv(bestmodel_k, file="PolygenicScores/Phenotype/Best_Model_K.csv")


prd <- data.frame(hp=seq(from=range(pred$K_GPA)[1], to=range(pred$K_GPA)[2], length.out = nrow(pred)))
err <- predict(mod2_k, newdata=prd, se.fit=T)
prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

pdf("PolygenicScores/Phenotype/Graphic_K.pdf", onefile = T)
ggplot(prd,aes(x=pred$K_GPA, y=fit))+
  labs(x= "K", y="Polygenic Score")+
  geom_point(aes(x=pred$K_GPA, y=polygenic_score_k$SCORE), alpha=0.5) +
  geom_smooth(method = "loess", color ="red", size=1.5)+
  geom_smooth(aes(y = lci), color = "red", linetype = 2 , size=0.5) +
  geom_smooth(aes(y = uci), color = "red", linetype = 2 , size=0.5) +
  theme_bw() +
  theme(axis.title.y = element_text(size=12, color = "black", face = "bold"),
        axis.title.x = element_text(size=12, color = "black", face = "bold"))
dev.off()



#######Mg
mod1_mg <- lm(polygenic_score_mg$SCORE ~ pred$Mg_GPA) #linear
mod2_mg <- lm(polygenic_score_mg$SCORE ~ pred$Mg_GPA + I(pred$Mg_GPA^2))#quadratic
mod3_mg <- lm(polygenic_score_mg$SCORE ~ 1)#Null

model_sel_mg = model.sel(mod1_mg,mod2_mg,mod3_mg)
model_sel_mg #model2
summary(mod1_mg)
bestmodel_mg = capture.output(summary(mod1_mg))

#Save results:
write.csv(model_sel_mg, file="PolygenicScores/Phenotype/Model_Selection_Mg.csv")
write.csv(bestmodel_mg, file="PolygenicScores/Phenotype/Best_Model_Mg.csv")


prd <- data.frame(hp=seq(from=range(pred$Mg_GPA)[1], to=range(pred$Mg_GPA)[2], length.out = nrow(pred)))
err <- predict(mod1_mg, newdata=prd, se.fit=T)
prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

pdf("PolygenicScores/Phenotype/Graphic_Mg.pdf", onefile = T)
ggplot(prd,aes(x=pred$Mg_GPA, y=fit))+
  labs(x= "Mg", y="Polygenic Score")+
  geom_point(aes(x=pred$Mg_GPA, y=polygenic_score_mg$SCORE), alpha=0.5) +
  geom_smooth(method = "loess", color ="red", size=1.5)+
  geom_smooth(aes(y = lci), color = "red", linetype = 2 , size=0.5) +
  geom_smooth(aes(y = uci), color = "red", linetype = 2 , size=0.5) +
  theme_bw() +
  theme(axis.title.y = element_text(size=12, color = "black", face = "bold"),
        axis.title.x = element_text(size=12, color = "black", face = "bold"))
dev.off()


#######Cu
mod1_cu <- lm(polygenic_score_cu$SCORE ~ pred$Cu_GPA) #linear
mod2_cu <- lm(polygenic_score_cu$SCORE ~ pred$Cu_GPA + I(pred$Cu_GPA^2))#quadratic
mod3_cu <- lm(polygenic_score_cu$SCORE ~ 1)#Null

model_sel_cu = model.sel(mod1_cu, mod2_cu, mod3_cu)
model_sel_cu #model2
summary(mod1_cu)
bestmodel_cu = capture.output(summary(mod1_cu))

#Save results:
write.csv(model_sel_cu, file="PolygenicScores/Phenotype/Model_Selection_Cu.csv")
write.csv(bestmodel_cu, file="PolygenicScores/Phenotype/Best_Model_Cu.csv")


prd <- data.frame(hp=seq(from=range(pred$Cu_GPA)[1], to=range(pred$Cu_GPA)[2], length.out = nrow(pred)))
err <- predict(mod1_cu, newdata=prd, se.fit=T)
prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

pdf("PolygenicScores/Phenotype/Graphic_Cu.pdf", onefile = T)
ggplot(prd,aes(x=pred$Cu_GPA, y=fit))+
  labs(x= "Cu", y="Polygenic Score")+
  geom_point(aes(x=pred$Cu_GPA, y=polygenic_score_cu$SCORE), alpha=0.5) +
  geom_smooth(method = "loess", color ="red", size=1.5)+
  geom_smooth(aes(y = lci), color = "red", linetype = 2 , size=0.5) +
  geom_smooth(aes(y = uci), color = "red", linetype = 2 , size=0.5) +
  theme_bw() +
  theme(axis.title.y = element_text(size=12, color = "black", face = "bold"),
        axis.title.x = element_text(size=12, color = "black", face = "bold"))
dev.off()


#######Pilocarpine
mod1_pilo <- lm(polygenic_score_pilo$SCORE ~ pred$pilocarpine) #linear
mod2_pilo <- lm(polygenic_score_pilo$SCORE ~ pred$pilocarpine + I(pred$pilocarpine^2))#quadratic
mod3_pilo <- lm(polygenic_score_pilo$SCORE ~ 1)#Null

model_sel_pilo = model.sel(mod1_pilo, mod2_pilo, mod3_pilo)
model_sel_pilo #model2
summary(mod2_pilo)
bestmodel_cu = capture.output(summary(mod2_pilo))

#Save results:
write.csv(model_sel_pilo, file="PolygenicScores/Phenotype/Model_Selection_Pilocarpine.csv")
write.csv(bestmodel_pilo, file="PolygenicScores/Phenotype/Best_Model_Pilocarpine.csv")


prd <- data.frame(hp=seq(from=range(pred$pilocarpine)[1], to=range(pred$pilocarpine)[2], length.out = nrow(pred)))
err <- predict(mod2_pilo, newdata=prd, se.fit=T)
prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

pdf("PolygenicScores/Phenotype/Graphic_Pilocarpine.pdf", onefile = T)
ggplot(prd,aes(x=pred$pilocarpine, y=fit))+
  labs(x= "Pilocarpine", y="Polygenic Score")+
  geom_point(aes(x=pred$pilocarpine, y=polygenic_score_pilo$SCORE), alpha=0.5) +
  geom_smooth(method = "loess", color ="red", size=1.5)+
  geom_smooth(aes(y = lci), color = "red", linetype = 2 , size=0.5) +
  geom_smooth(aes(y = uci), color = "red", linetype = 2 , size=0.5) +
  theme_bw() +
  theme(axis.title.y = element_text(size=12, color = "black", face = "bold"),
        axis.title.x = element_text(size=12, color = "black", face = "bold"))
dev.off()


###Graphics by population:
#####6. PILOCARPINE VALUES BY POPULATIONS ---------------------------------
#A.Load Pilocarpine Values:
pilo_pred = read.csv("./Results/Outputs/GPA_pilo_pred.csv")
head(pilo_pred)


#B. LOAD POP INFORMATION BASED ON sNMF
# Load files
pop = read.csv("data_raw/pop_LEA_DAPC_TESS_pilocarpus.csv", row.names = 1)
metafiles = as.data.frame(snps@meta)
# selecting the individuals:
pop_filtered = pop[pop$ind_ID %in% metafiles$ind_ID,]
# Check if there are some difference between files. If they are corrected, function will return "character (0)"
setdiff(as.character(pop_filtered$ind_ID), as.character(pilo_pred$id))
setdiff(as.character(pilo_pred$id), as.character(pop_filtered$ind_ID))
# Check if identical samples were selected in the same order. If they are corrected, function will return "TRUE"
identical(as.character(pop_filtered$ind_ID), as.character(pilo_pred$id))
#verify
pop_filtered$PopID_snmf


#Graphic:
pdf("PolygenicScores/Phenotype/Graphic_Pilocarpine_by_Population.pdf", onefile = T)
ggplot(prd,aes(x=pred$pilocarpine, y=fit))+
  labs(x= "Pilocarpine", y="Polygenic Score")+
  geom_point(aes(x=pred$pilocarpine, y=polygenic_score_pilo$SCORE, fill=factor(pop_filtered$PopID_snmf)), shape=21, size=2.5) +
  scale_fill_manual(values=c("red", "blue", "yellow", "pink"), labels= c("POP1", "POP2", "POP3", "POP4")) +
  geom_smooth(method = "loess", color ="black", size=1.5)+
  geom_smooth(aes(y = lci), color = "black", linetype = 2 , size=0.5) +
  geom_smooth(aes(y = uci), color = "black", linetype = 2 , size=0.5) +
  theme_bw() +
  guides(fill=guide_legend(title="Populations", override.aes = list(size=5)))+
  theme(axis.title.y = element_text(size=12, color = "black", face = "bold"),
        axis.title.x = element_text(size=12, color = "black", face = "bold"))
dev.off()

##END
