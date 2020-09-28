
library(vcfR)
library(PopGenReport)
library(RColorBrewer)
library(lattice)

#F. Create a Genind object based on VCFR
vcf_cand = read.vcfR("./sPCA/Pilocarpine/Pilocarpus_Pilo_SNPs_Candidates.vcf", verbose = FALSE)
vcf_cand
genind = vcfR2genind(vcf_cand)
genind

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

#C. Adding POP information to genind:
pop(genind) = pop_filtered$PopID_snmf
genind


#Calculated Frequency 
allele_freq = allele.dist(genind, mk.figures = F)

#Rename lists and reorganized the df:
names(allele_freq$frequency) = vcf_cand@fix[,3]
df = as.data.frame(dplyr::bind_rows(allele_freq$frequency, .id = NULL)) 
df = t(df)
head(df)

#Created a matrix with REF allele values for each population
#POP1
pop1 = df[which(grepl(pattern = "(.1)$", rownames(df))),]
rownames(pop1) = vcf_cand@fix[,3]
pop1
pop1 = pop1[,1]

#POP2
pop2 = df[which(grepl(pattern = "(.2)$", rownames(df))),]
rownames(pop2) = vcf_cand@fix[,3]
pop2
pop2 = pop2[,1]

#POP3
pop3 = df[which(grepl(pattern = "(.3)$", rownames(df))),]
rownames(pop3) = vcf_cand@fix[,3]
pop3
pop3 = pop3[,1]

#POP4
pop4 = df[which(grepl(pattern = "(.4)$", rownames(df))),]
rownames(pop4) = vcf_cand@fix[,3]
pop4
pop4 = pop4[,1]

# Create a merged matrix 
mt = cbind(pop1, pop2, pop3, pop4)
class(mt)
head(mt)


#Count number of SNPs with a specific frenquency in intervals of 0.1 by population:
breaks_count = seq(0.0, 1.0, 0.1)

pop1_count = hist(pop1, breaks_count, plot=F)$counts
pop2_count = hist(pop2, breaks_count, plot=F)$counts
pop3_count = hist(pop3, breaks_count, plot=F)$counts
pop4_count = hist(pop4, breaks_count, plot=F)$counts


# Create a merged df
df_count = as.data.frame(cbind(pop1_count, pop2_count, pop3_count, pop4_count))
rownames(df_count) = breaks_count[-1]
colnames(df_count) = c("POP1", "POP2", "POP3", "POP4")
head(df_count)


# Define the number of colors you want
n_max = max(df_count) #7 counts is the maximum value
nb.cols = 10*n_max
mycolors = colorRampPalette(brewer.pal(9, "Blues"))(nb.cols)


# plot ASF
pdf("./AFS/Pilocarpine_AFS_Population.pdf", onefile =F)
levelplot(as.matrix(df_count), col.regions=mycolors,
          xlab=NULL, ylab=NULL,  sub="Allele Frequency",
          main = "Allele Spectrum Frequency\n28 SNPs Candidates\nPilocarpine")
dev.off()


#END