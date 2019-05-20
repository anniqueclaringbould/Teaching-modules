##################################
# GWAS tutorial
# 
# Steps to take:
# 
#  1. Import the data in R and explore them
 
#2. Quality control per SNP
#     i) Remove SNPs with minor allele frequency (MAF) < 0.01
#     i) Remove SNPs that deviate from Hardy Weinberg equilibrium (HWE) using a P threshold < 1e-6
#The HWE is a principle stating that the genetic variation in a population will remain constant from one generation to the next in  the absence of disturbing factors. Departure from this equilibrium can be an indicator of potential genotyping errors, population stratification or even actual association with the phenotype under study
#      i) Remove SNPs with missing genotype rate < 0.99
 
# 3. Quality control per individual
#     i) Identify and remove related individuals
#     i) Perfomr principal component analysis (PCA)and plot the first two components to identify and remove population outliers
 
# 4. Perform genome-wide association testing
#     i) do initial association analysis, linear
#     i) do association analysis, logistic (Bonus)
#     i) prepare a manhattan plot
# 
##################################

#Step 1: Import the data in R and explore them
install.packages("BEDMatrix")

setwd("~/PhD/gwas_tutorial/simulate_phenotypes/initial_dataset")

require(BEDMatrix)

#Import the bed file in R
genotype_matrix <- BEDMatrix("celiac_gwas.bed")

#Explore the data
genotype_matrix[1,] #shows all the SNPs for the first individual
genotype_matrix[,1] #shows all individual genotypes per SNP

dim(genotype_matrix) # shows how many individuals and genetic variants are present

colnames(genotype_matrix) # shows the names of all SNPs
rownames(genotype_matrix) # shows the names of all individuals 

genotype_matrix[1:3, 1:5] # shows the first three individuals and first five SNPs 

#Note that SNP genotypes are encoded as 0 (AA), 1 (AB or BA) and 2 (BB)
##################################################################################################
#Step 2: Quality control per SNP
#Step 2.1: Remove SNPs with low minor allele frequency (MAF)

#First, identify the minor allele frequency per SNP

minor_allele_frequency <- function(row){
	a1_count = 2* sum(row == 2) + sum(row ==1)
	return(a1_count / (2*length(row)))
}

#Identify the allele frequency of the first SNP
minor_allele_frequency(genotype_matrix[,1])

#Now, do it for all SNPs
minor_allele_frequency_per_snp <- apply(genotype_matrix, 2, minor_allele_frequency)

#Plot the minor allele frequency (MAF) for all SNPs
hist(minor_allele_frequency_per_snp)

#Identify SNPs that have a MAF lower than 0.01. These SNPs will be removed in the next steps

maf_filter <- minor_allele_frequency_per_snp < 0.01

#Check how many SNPs show a MAF < 0.01 that will be removed
sum(maf_filter, na.rm=T) # 4098 SNPs

#Now, we can remove the SNPs showing a MAF < 0.01 
genotype_matrix_maf_filtered <- genotype_matrix[,!maf_filter]

#Check how many SNPs remained
dim(genotype_matrix_maf_filtered)

# For simplicity's sake, we will just continue with the 
# unfiltered genotype matrix, and we only filter variants out
# when we known the total number oF SNPs that should be removed

# Step 2.2: Remove SNPs that deviate from HWE using a P threshold < 1e-6

#Hardy Weinberq equilibrium (HWE)
#Check the number of individuals and SNPs before quality control
num_individuals <- nrow(genotype_matrix)
head(num_individuals)
num_variants <- ncol(genotype_matrix)
head(num_variants)

#Find the expected genotype counts per SNP
#Make a matrix per column with the number of expected homozygote (AA), hetererozygote (AB) and homozygote (BB) genotypes per SNP
expected_genotype_count <- cbind( 
	(1-minor_allele_frequency_per_snp)^2 * num_individuals,
	2*(1-minor_allele_frequency_per_snp) * minor_allele_frequency_per_snp * num_individuals,
	minor_allele_frequency_per_snp^2 * num_individuals
	)

#Have a look at the first 10 SNPs
expected_genotype_count[1:10,]
#Now, lets have a look at the observed counts of genotypes:
find_oberserved_count <- function(row){
	counts <- c(
			sum(row == 0, na.rm=T), #number of individuals with genotype AA
			sum(row == 1, na.rm=T), #number of individuals with genotype AB
			sum(row == 2, na.rm=T)  #number of individuals with genotype BB
		)
	return(counts)
}

#Try this function for SNP 42
find_oberserved_count(genotype_matrix[,42])

#Calculate observed genotype counts for all SNPs
observed_genotype_count <- t(matrix(apply(genotype_matrix, 2, find_oberserved_count), nrow=3))

#Now, compare the observed with expected genotype counts and calculate the chi sq statistics
chi_sq_statistics <- apply((observed_genotype_count - expected_genotype_count)^2 / expected_genotype_count, 1, function(x) sum(x, na.rm=T))
#Show the chi sq value for SNP 1
chi_sq_statistics[1]

# Plot an histogram of chi sq statistics
hist(chi_sq_statistics)

#These values are chi squared distributed with one degree of freedom, so we can then identify the P value
hardy_weinberg_p_values <- 1 - pchisq(chi_sq_statistics,1)

#Now, we filter for P values that deviate from HWE (P < 1e-6)
hwe_p_val_filter <- hardy_weinberg_p_values < 1e-6
			   
#Check how many SNPs deviate from HWE and should be removed later
genotype_matrix_hwe_filtered <- genotype_matrix[,!hwe_p_val_filter]
dim(genotype_matrix_hwe_filtered)			   
			   
# Similarly as above when filtering out SNPs with MAF < 0.01, we can filter out SNPs that deviate from HWE
# but to keep the number of datasets limited, we don't do that at this point  

#Step 2.3: Remove SNPs with high missing genotype rate (missingness)

#SNPs with high missingness are encoded as NA in our data.So, we count them and divide them by the number of individuals
missingness <- apply(genotype_matrix, 2, function(x) sum(is.na(x))) / num_individuals

#Plot the distribution of missingness
hist(missingness)

#Now, filter for missingness
missingness_filter <- missingness != 0
##Check how many SNPs show high missingness and should be removed later
genotype_matrix_missingness_filtered <- genotype_matrix[,!missingness_filter]
dim(genotype_matrix_missingness_filtered)
		     
###
#Step 2.4: Filter out all SNPs that follow the below criteria:
#MAF < 0.01
#HWE P < 1e-6
#high missingness (encoded as NA)		     
##		
# We can now filter out all snps, and save the filtered matrix
genotype_matrix_geno_qc <- genotype_matrix[,!missingness_filter & !hwe_p_val_filter & !maf_filter]
minor_allele_frequency_geno_qc <- minor_allele_frequency_per_snp[!missingness_filter & !hwe_p_val_filter & !maf_filter]

#Check how many SNPs and individuals remain after quality control
dim(genotype_matrix_geno_qc)

##################################################################################
##Step 3: Quality control per individual

##Step 3.1: To identify related individuals
#Make a genetic relatedness matrix (complicated in the math)

maf_normalizer <- sqrt(2*(1-minor_allele_frequency_geno_qc)*minor_allele_frequency_geno_qc)

#Following Yang et al 2010:
# grm = W %*% W / N
# where w_{ij} = (x_{ij} - 2*p) / sqrt(2*p*(1-p))
# p is the allele frequency.

w_mat <- t(apply(genotype_matrix_geno_qc, 
				1,
			 	function(x) as.vector((x - 2*(minor_allele_frequency_geno_qc)) / maf_normalizer)
			 	)
			)

# Takes about a minute but 2 GB of memory
grm <- (w_mat %*% t(w_mat)) / num_variants
		     
		     
#Plot an histogram
hist(grm[lower.tri(grm)])
		     
	
# Remove individuals that show a relatedness value more than 0.1 
relatedness_indice <- which(grm > 0.1, arr.ind=TRUE) #find the indices in the matrix.

relatednes_individuals <- unique(relatedness_indice[relatedness_indice[,1] != relatedness_indice[,2],][,1]) #make sure not to remove individuals which are compared as the same.

# For simplicity's sake, we remove all related individuals.
# and apply the relatedness filter.
relatedness_filter <- rep(FALSE, num_individuals)
relatedness_filter[relatednes_individuals] <- TRUE

#Check the number of related individuals
genotype_matrix_relatedness_filtered <- genotype_matrix_geno_qc[!relatedness_filter,]
dim(genotype_matrix_relatedness_filtered)

## Step 3.2: Perform PCA on the individuals and identify population outliers
```{r, message=FALSE, warning=FALSE, eval=TRUE}
principal_components <- prcomp(grm) # relatively easy in R.
```
#Now, we can plot the first two principal components.
plot(principal_components$x[,1], principal_components$x[,2]) # plot the first two principal components
```

#Using some arbitrary threshold, we can now filter the full matrix
#We decide on a threhold for the principal components to filter on
pca_filter <- principal_components$x[,2] < -0.25

#Identify the individuals that showing PCAs < -0.25
rownames(genotype_matrix_geno_qc)[pca_filter]

#Check how many individuals are identified as population outliers
genotype_matrix_pcas_filtered <- genotype_matrix_geno_qc[!pca_filter,]
dim(genotype_matrix_pcas_filtered)

#Remove individuals from the genotype matrix based on relatedness and PCA, and we have concluded our QC steps
genotype_matrix_post_qc <- genotype_matrix_geno_qc[!pca_filter & !relatedness_filter,]

#Also, remove these individuals from the phenotypes.
phenotypes_post_qc <- phenotypes[!pca_filter & !relatedness_filter]

dim(genotype_matrix_post_qc) # leaves 500 individuals! 114404 variants!

```

## Step 4: Testing for genetic associations
#Now we can test whether our genetic variants are associated with celiac disease. Using a linear model, we can associate all the genotypes to the disease.
#Let's examine if the third SNP is associated with celiac disease.
summary(lm(genotype_matrix_post_qc[,1]~phenotypes_post_qc))
```

#Next, we create the association table using an R function, which outputs the slope. Then, we can associate the disease phenotype to any genotype:
do_quantative_association <- function(genotypes, phenotypes){

	sumdat <- summary(lm(genotypes~phenotypes))
	return(as.vector(sumdat$coefficients[2,])) #second row is the column.
}

# It takes a while (about 2 minutes)
associations <- t(matrix(apply(genotype_matrix_post_qc, 2, do_quantative_association, phenotypes=phenotypes_post_qc), nrow=4))

```

#We have calculated the associations and we would like to plot the significance level over the chromosomal positions.
#SNP names are in the format <chr>:<position>_<effect_allelle>. Therefore, we split this format and turn it into a dataframe.
positions_of_snps <- do.call(rbind,strsplit(colnames(genotype_matrix_post_qc), ":|_"))

associations_with_position <- cbind(positions_of_snps, associations)

#Make a dataframe for plotting with ffplot
associations_with_position_df <- data.frame(chr = as.numeric(associations_with_position[,1]), 
                                 position = as.numeric(associations_with_position[,2]),allele = associations_with_position[,3],
                                 beta = as.numeric(associations_with_position[,4]),se = as.numeric(associations_with_position[,5]),
                                 t_stat = as.numeric(associations_with_position[,6]),p_val = as.numeric(associations_with_position[,7]))
```


#Finally, we plot the manhattan plot using ggplot2
require(ggplot2)
ggplot(associations_with_position_df, aes(x=position, y=-log10(p_val), col=as.factor(chr))) + 
  facet_grid(.~chr, scales="free_x") + 
  geom_point()
 
```

#Save the plot in PDF format
pdf("assocplot.pdf", width=6, height=4)
ggplot(associations_with_position_df, aes(x=position, y=-log10(p_val), col=as.factor(chr))) + 
  facet_grid(.~chr, scales="free_x") + 
  geom_point()
  dev.off()
		     
