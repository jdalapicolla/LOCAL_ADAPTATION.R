####################### VALE INSTITUTE OF TECHNOLOGY ##########################

########################      GEA-GPA TUTORIAL       ###########################
######################## STEP 07: FINAL VENNDIAGRAMM  ##########################

### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###



#PRE-ANALYSIS ------------------------------------------------------------------


#AIMS:
#A. COMPARING CANDIDATES SNPS IDENTIFIED WITH RDA AND LFMM2 USING A VENNDIAGRAM


##LOAD THE PACKAGES
library(VennDiagram)


## LOAD FUNCTIONS TO BE USED ON THIS STEP.
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}


#1. LOAD GEA AND GPA RESULTS-----------------------------------------------------------
rda_snps_gea = vcfLink("GEA/partial_RDA/snps_candidate_gea_RDA_mahalanobis_pilocarpus_subset.vcf", overwriteID=F) # here overwrite need to be F, otherwise we will loss the snps name

lfmm_snps_gea = vcfLink("GEA/LFMM/snps_candidate_lfmm_GEA_pilocarpus.vcf", overwriteID=F)

rda_snps_gpa = vcfLink("GPA/partial_RDA/snps_candidate_gea_RDA_mahalanobis_pilocarpus_subset.vcf", overwriteID=F) # here overwrite need to be F, otherwise we will loss the snps name

lfmm_snps_gpa = vcfLink("GPA/LFMM/snps_candidate_lfmm_GPA_pilocarpus.vcf", overwriteID=F)


#3. VENN DIAGRAM--------------------------------------------------------------- 
#select SNPs names:
RDA_GEA = unique(Chrom(rda_snps_gea)[,1])
LFMM2_GEA = unique(Chrom(lfmm_snps_gea)[,1])
RDA_GPA = unique(Chrom(rda_snps_gpa)[,1])
LFMM2_GPA = unique(Chrom(lfmm_snps_gpa)[,1])

## Creating a Venn Diagram and save it as pdf
pdf("Results/VennDiagram/VennDiagram_SNPscandidates_GPA_GEA.pdf", onefile = T)
plot.new()
grid.draw(w1 <- venn.diagram(list(RDA_GEA=RDA_GEA, LFMM2_GEA=LFMM2_GEA, RDA_GPA=RDA_GPA, LFMM2_GPA=LFMM2_GPA),
                             lty = c("blank", "blank", "blank", "blank"),
                             fill = c("red", "blue", "yellow", "black"),
                             alpha = c(0.5, 0.5, 0.5, 0.5), cat.cex = 1.2, cex= 1.5,  cat.pos = 0, 
                             filename=NULL ))

dev.off()

  ##END   


