####################### VALE INSTITUTE OF TECHNOLOGY ##########################

########################      GEA-GPA TUTORIAL       #########################
########################   STEP 11: sPCA for GPA     #########################

### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###



# PRE-ANALYSIS-----------------------------------------------------------------

#Load Packages:
library(r2vcftools)
library(adespatial)
library(spdep)
library(raster)
library(rgdal)
library(adegenet)


### Summary stats
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}

fgraph <- function(obj){
  # use multispati summary
  sum.obj <- summary(obj)
  # compute Imin and Imax
  L <- listw2mat(eval(as.list(obj$call)$listw))
  Imin <- min(eigen(0.5*(L+t(L)))$values)
  Imax <- max(eigen(0.5*(L+t(L)))$values)
  I0 <- -1/(nrow(obj$li)-1)
  # create labels
  labels <- lapply(1:length(obj$eig),function(i) bquote(lambda[.(i)]))
  # draw the plot
  xmax <- eval(as.list(obj$call)$dudi)$eig[1]*1.1
  par(las=1)
  var <- sum.obj[,2]
  moran <- sum.obj[,3]
  plot(x=var,y=moran,type='n',xlab='Inertia',ylab="Spatial autocorrelation (I)",
       xlim=c(0,xmax),ylim=c(Imin*1.1,Imax*1.1),yaxt='n')
  text(x=var,y=moran,do.call(expression,labels))
  ytick <- c(I0,round(seq(Imin,Imax,le=5),1))
  ytlab <- as.character(round(seq(Imin,Imax,le=5),1))
  ytlab <- c(as.character(round(I0,1)),as.character(round(Imin,1)),
             ytlab[2:4],as.character(round(Imax,1)))
  axis(side=2,at=ytick,labels=ytlab)
  rect(0,Imin,xmax,Imax,lty=2)
  segments(0,I0,xmax,I0,lty=2)
  abline(v=0)
  title("Spatial and inertia components of the eigenvalues")
}


# 1. Loading Files ----
#A. Load snps data set to take the snps and individuals ID
snps =  vcfLink("vcf/pilocarpus_adap_filtered.vcf", overwriteID=T)
snps@meta
VCFsummary(snps) ## 264 individuals and 19025 SNPs.


#B. Load candidate SNPs in RDA
rda_snps = vcfLink("GPA/partial_RDA/snps_candidate_gea_RDA_mahalanobis_pilocarpus_subset.vcf", overwriteID=F) # here overwrite need to be F, otherwise we will loss the snps name
VCFsummary(rda_snps) #264 individuals and 133 SNPs.

#C. Load candidate SNPs in LFMM2
lfmm_snps = vcfLink("GPA/LFMM/snps_candidate_lfmm_GPA_pilocarpus.vcf", overwriteID=F)
VCFsummary(lfmm_snps)# 264 individuals and 204 SNPs.

#D. Merging all candidates:
all_candidate = unique(c(rda_snps@site_id,lfmm_snps@site_id))

#E.Subseting a VCF
snps_cand = Subset(snps, sites = all_candidate)
VCFsummary(snps_cand)# 264 individuals and 335 SNPs.

#F. Saving the subsetted VCF:
Save(snps_cand, "./sPCA/GPA/Pilocarpus_All_SNPs_Candidates.vcf")

#G. Create a Genind object based on VCFR
vcf_cand = read.vcfR("./sPCA/GPA/Pilocarpus_All_SNPs_Candidates.vcf", verbose = FALSE)
vcf_cand
genind = vcfR2genind(vcf_cand)
genind


#2. Geographical Information -----
#You can test different connection networks:
#Delaunay triangulation (type=1) = more suited to uniform sampling. Other softwares use this one, for comparasion it is a good one.
#
#Gabriel graph (type=2). The Gabriel graph is a subgraph of the Delaunay triangulation. It can be found in linear time if the Delaunay triangulation is given (Matula & Sokal 1980). The Gabriel graph contains, as subgraphs, the Euclidean minimum spanning tree (type=4), the relative neighborhood graph (type=3), and the nearest neighbor graph (type=6). 
#
#Relative neighbours (type=3)
#
#Minimum spanning tree (type=4)
#
#Neighbourhood by distance (type=5) =  You need to inform minimum (d0) and maximum distance (d1) IN NUMBER OF INDIVIDUALS 0 AND 2 IS A PATTERN. Biologically-defined distance would be better in this case (home-range, dispersion for plants)
#
#K nearests neighbours (type=6) = You need to inform the number of neighbours per point (argument k).Shoul be lower than one-third of the of the number of data points
#
#Inverse distances (type=7) = recommended when clustering of samples is present. This is not a true neighbouring graph: all sites are neighbours, but the spatial weights are directly proportional to the inversed spatial distances. You need to inform the minimum distance between any two distinct points (dmin) and the exponent of the inverse distance matrix (a).


#3. GEOGRAPHIC DATA AND CONNECTION NETWORK ----
#A. Geographical information for all samples
coordinates = snps_cand@meta[,4:5]
head(coordinates)

#B. Create Connection network. Choose Neighbourhood by distance (type=5) for most cases.
#CN = chooseCN(coordinates, type=1, plot=T)
#CN2 = chooseCN(coordinates, type=2, plot=T)
#CN3 = chooseCN(coordinates, type=3, plot=T)
#CN4 = chooseCN(coordinates, type=4, plot=T)
#CN = chooseCN(coordinates, type=5, d1 = 0, d2 = 2, plot=T) #not positive
CN = chooseCN(coordinates, type=6, k = 5, plot=T)
#CN = chooseCN(coordinates, type=7, dmin = 2, a = 2, plot=T)

#C. Then convert your CN to a listw using the nb2listw function.
cn = nb2listw(CN)



# 4. Running sPCA Analysis ----
###4.1. RUNNING sPCA:
#Here we will run sPCA following the tutorial provide by Jombart.
#  http://adegenet.r-forge.r-project.org/files/tutorial-spca.pdf

#A. Number of PCs
n_pcs = nrow(genind@tab)
  
#B. Create a PCA object
input_scaled = scaleGen (genind, center = TRUE, scale = TRUE, NA.method = "mean")
pca_input = dudi.pca(input_scaled, center = TRUE, scannf = FALSE, nf = n_pcs)

#C. Running sPCA:
mySpca = multispati(pca_input, cn, scannf = F)
#verify number of PC in global and local structure
mySpca
plot(mySpca)
#verify lag vector onto the principal axesdata frame with n rows and (nfposi + nfnega) columns
head(mySpca$ls)


###4.2. TEST LOCAL AND GLOBAL STRUCTURE
# Lollipop plot interpretation: significant if lollipop is far from distribution of permutation values, not significant if lollipop falls along the distribution of permutation values.
# GLOBAL TEST: neighbouring individuals are more similar than expected
# LOCAL TEST: neighbouring individuals are more dissimilar than expected

#A. Global test
myGtest = global.rtest(pca_input$tab, cn, nperm=9999)
myGtest
plot(myGtest)

#B. Local test
myLtest = local.rtest(pca_input$tab, cn, nperm=9999)
myLtest
plot(myLtest)


#------------------------------------------------------------------------------
#                         5. Save the Results 
#------------------------------------------------------------------------------
###5.1. EXTRACT RESULTS TO DO THE INTERPOLATON ON QGIS
#verify lag vector onto the principal axesdata frame with n rows and (nfposi + nfnega) columns
head(mySpca$ls)
head(cbind(mySpca$ls,coordinates))
write.csv(cbind(mySpca$ls,coordinates), "./sPCA/GPA/Pilocarpus_sPCA_Axis_QGIS.csv")

###5.2. SAVE RESULTS AS GRAPHS AND TABLES
#A. Tables
results_sPCA = summary(mySpca)
write.csv(results_sPCA,"./sPCA/GPA/Pilocarpus_sPCA_Summary.csv")

#B. Spatial and Inertial Components:
pdf("./sPCA/GPA/Spatial_Inertia.pdf", onefile =F)
plot.new()
fgraph(mySpca)
dev.off()

#C. Barplot with Global and Local Structure in spectral color pattern
pdf("./sPCA/GPA/Global_Local_Structure_Spectral.pdf", onefile =F)
barplot(mySpca$eig, main="A variant of the plot\n of sPCA eigenvalues",
        col=spectral(length(mySpca$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="black")
dev.off()

#D. Barplot with Global and Local Structure in red color for PC significants
#Global Significant PCs:
n_posi = 2
#For Local PCs test diferent values to put color on graph, for no one PC let ir on 0. Usually you need to plus number of positive PC and Negative to put correctly colors:
n_eig = length(mySpca$eig)-0
#verify colors:
barplot(mySpca$eig, main="Eigenvalues of sPCA", col=rep(c("red","grey"),c(n_posi,n_eig)))

pdf("./sPCA/GPA/Global_Local_Structure_Significant.pdf", onefile =F)
barplot(mySpca$eig, main="Eigenvalues of sPCA", col=rep(c("red","grey"),c(n_posi,n_eig)))
dev.off()

#E. Colorplot of mySpca. It represents a cloud of points with colors corresponding to a combination of 1,2 or 3 quantitative variables, assigned to RGB (Red, Green, Blue) channels. For instance, this can be useful to represent up to 3 principal components in space.
#define the number of variables UP TO THREE (Positive + Negatives PCs)
variables = 2
pdf("./sPCA/GPA/Colorplot_Significant_PC.pdf", onefile =F)
colorplot(coordinates, as.data.frame(mySpca$ls[,1:variables]), cex=3, main="Colorplot of sPCA of Neutral SNPs")
dev.off()

#F. Global Test of mySpca
pdf("./sPCA/GPA/Global_Test_Significant.pdf", onefile =F)
plot(myGtest, main="Simulation for \n Moran's Eigenvector Maps (MEMs)",  xlab = "Observed R²", ylab = "Frequency")
dev.off()

#G. Local Test of mySpca
pdf("./sPCA/GPA/Local_Test_Significant.pdf", onefile =F)
plot(myLtest, main="Simulation for \n Moran's Eigenvector Maps (MEMs)",  xlab = "Observed R²", ylab = "Frequency")
dev.off()


## Make Interpolated Maps using QGIS 3.4     
#1. Ins     tall "Processing" plug-in in QGis
#2. Go to "Add a delimited text file" and open the file with spca axes and coordinates
#3. Save the opened file as a shapefile
#4. Interpolate each axis separately using the Interpolation tool in the Processing Toolbox:
    #  Processing Toolbox > I nterpolation > Interpolation IDW.
#5. At the Interpolation window: 
# a) select the shapefile with spca axes; 
# b) add the axis in the second box as a vertorial entry layer;
# c) select the "extension"; hillshade_carajas
# d) Adjust the number of columns and rows; 500 rows  
# e) run the interpolation (the raster will app ear in the main window);
#6. Repeat the previous steps to each axis
#7. Combine the three interpolated axis raster in a rgb raster using Raster > Micellaneous > Mosaic:
# a) Select the three interpolated axis raster;
# b) Select the option "put each file in a separate band";
# c) Run to build your mosaic.

##END