####################################################
#
# R script for Project. 
# by Dan Jin. 
# last update: Nov. 20, 2012
# 
# To run this script, please copy ProjectScript.r to the same fold
# where containing the input data files (project_genotype.txt, 
# project_phenotype.txt, project_information.txt), 
# then type the following lines in R.
#
# >source("ProjectScript.r")
# >result.ls=list()
#
# For Pre-test (take really very long time to code Xa):
# >result.ls=control.fun(-1,-1,-1,0)
#
# For GWAS:
# >result.ls=control.fun(0.25,0.15,0.05,0.1)
#
# Result of this script:
# 1. Create two .txt files containing filtered genotype and phenotype data.
# 2. return one list, result.ls. See below for more details about the contents in the list.
# 3. create 11 plots in png format. See below for more details about each png file.
#
# Summary of the content in result.ls:
# 1. $pval: p values calculated without covariate
# 2. $pval.covariate: p values calculated with covariate
# 3. $hits.ls: all genotype marker (marker number and their p value)significant than Bonferroni correction
# 4. $hits.info: chromosome and SNP position information for each hit in hits.ls
#
# Structure of result.ls:
# result.ls
# |-$pval
# |-$pval.covariate
# |-$hits.ls
# | |-$position
# | |-$pval
# |-$hits.info
#   |-$Chr
#   |-$SNP
#   |-$Chr.SNP
#
#List of .png image files:
# 1.	01_Phenotype distribution.png
# 2.	02_Percentage of missing data.png
# 3.	03_Percentage of missing individuals.png
# 4.	04_MAF.png
# 5.	05_Post-filter phenotype distribution.png
# 6.    06_PCA.png
# 6.	07.1_QQ Plot.png
# 7.	07.2_QQ Plot with covariate.png
# 8.	08.1_Manhattan Plot with Bonferroni correction.png
# 9.	08.2_Manhattan Plot with covariate with Bonferroni correction.png
# 10.	09_Zoomed in manhattan Plot with Bonferroni correction.png
#
####################################################




control.fun <- function(threshold.geno = -1, threshold.indi = -1, threshold.MAF = -1, threshold.PCA = 0) {

	threshold.geno = threshold.geno
	threshold.indi = threshold.indi
	threshold.MAF = threshold.MAF
	threshold.PCA = threshold.PCA
	
	
	
	
	######################### Import, check and filter data #########################
	
	#:::::::::::::::::Import the phenotype and genotype data:::::::::::::::::
	
	sample.import.fun <- function() {
		import.ls = list()
		import.ls$genotype.df <- read.table("project_genotype.txt", header = F, sep = " ")
		#Note that there are 162 individuals (sample.size; rows) and 19604 columns (6 non-SNP columns and 19598 alleles).
		size = dim(import.ls$genotype.df)
		sample.size = size[1] #=162 individuals
		num.alleles = size[2] - 6 #=19598 alleles
		num.sample = num.alleles/2 #=9799 SNPs

		import.ls$phenotype.df <- read.table("project_phenotype.txt", header = T, sep = " ")
		#Note that there are 137 phenotypes (rows).
		return(import.ls)

	}
	import.ls = list()
	import.ls = sample.import.fun()
	genotype.df = import.ls$genotype.df
	phenotype.df = import.ls$phenotype.df


	# Check data quality.
	# PHENOTYPE: Check whether phenotype is normally distributed.
cat("Check whether phenotype is normally distributed.\n")
	png("01_Phenotype distribution.png", width = 500, height = 500)
	hist(phenotype.df[, 3], main = "Phenotype distribution", ylab = "Frequency", xlab = "Average phenotype", breaks = 20)
	dev.off()
	cat("\"01_Phenotype distribution.png\" created.\n\n")

	# Re-order phenotype list in the order of IID
	ordered.pheno = phenotype.df[order(phenotype.df[, 2]), ]


	#:::::::::::::::::Clean data -- reorder phenotypes and filter phenotypes and genotypes:::::::::::::::::
	
	# Create a list to store filtered genotype.
	filtered.data = list()

	# Filter: GENOTYPE
	filter.genotype.fun <- function(phenotype.df, genotype.df, threshold.geno, threshold.indi, threshold.MAF) {
		cat("Filter data: Genotype\n")
		size = dim(genotype.df)
		sample.size = size[1] #=162 individuals
		num.alleles = size[2] - 6 #=19598 alleles
		num.sample = num.alleles/2 #=9799 SNPs


		# 1. Remove individuals with >20% missing data across all genotypes or without phenotypes
		if (threshold.geno == -1) {
			threshold.geno = (num.alleles - 2)/num.alleles # This threshold will only remove genotypes that have no data at all. Otherwise, it will cause problem when filtering genotypes with MAF!
		}

		cat("# 1. Remove individuals with >", threshold.geno * 100, "% missing data across all genotypes or without phenotypes.\n")
		filtered.geno.df = NULL # store the filtered result

		#Find individuals with >20% missing data across all genotypes
		miss.geno <- function(row) {
			return(sum(row == 0) - 4) #exclude the first 4 zeros in the first 6 columns
		}
		geno.mis.num = apply(genotype.df, 1, miss.geno) # store the number of missing data for each individual
		geno.mis.per = geno.mis.num/num.alleles
		png("02_Percentage of missing data.png", width = 500, height = 500)
		hist(geno.mis.per, main = "Percentage of missing data across all genotypes per person", ylab = "Number of individuals", 
			xlab = "Percentage of missing data")
		dev.off()
		cat("\"02_Percentage of missing data.png\" created.\n")
		geno.mis.bool = (geno.mis.num <= num.alleles * threshold.geno)

		#Find individuals without phenotypes-remove individuals
		miss.pheno <- function(row) {
			return((sum(row == phenotype.df[, 2])) != 0)
		}
		pheno.mis.bool = apply(as.matrix(genotype.df[, 2]), 1, miss.pheno)

		#Create filtered genotype data
		filtered.geno.bool = geno.mis.bool * pheno.mis.bool
		cat(sum(filtered.geno.bool == 1), "individuals(", sum(filtered.geno.bool == 1)/length(filtered.geno.bool) * 
			100, "%) have missing genotype data <=", threshold.geno * 100, "%", "with phenotype data.\n\n")
		filtered.geno.df = genotype.df[which(filtered.geno.bool == 1, arr.ind = TRUE), ]

		# update sample.size
		sample.size = nrow(filtered.geno.df)



		# 2. Remove genotype with >5% missing data across entire individual
		trunc.filtered.geno.df = filtered.geno.df[, 7:dim(filtered.geno.df)[2]] # Remove the first 6 colunms

		if (threshold.indi == -1) {
			threshold.indi = (sample.size - 1)/sample.size # This threshold will only remove genotypes that have no data at all. Otherwise, it will cause problem when filtering genotypes with MAF!
		}

		cat("# 2. Remove genotype with >", threshold.indi * 100, "% missing data across entire individual.\n")
		filtered.geno.df2 = NULL
		indi.mis.num = NULL # store the number of missing individuals for each genotype (threshold=5%)
		miss.indi <- function(col) {
			return(sum(col == 0))
		}
		indi.mis.num = apply(trunc.filtered.geno.df, 2, miss.indi) # store the number of missing individuals for each genotype
		indi.mis.per = indi.mis.num/sample.size
		png("03_Percentage of missing individuals.png", width = 500, height = 500)
		hist(indi.mis.per, main = "Percentage of missing individuals for each genotype", ylab = "Number of genotypes", 
			xlab = "Percentage of missing individuals")
		dev.off()
		cat("\"03_Percentage of missing individuals.png\" created.\n")
		indi.mis.bool = (indi.mis.num <= sample.size * threshold.indi)
		cat(sum(indi.mis.bool == 1), "SNPs (", sum(indi.mis.bool == 1)/length(indi.mis.bool) * 100, "%) have missing individuals <=", 
			threshold.indi * 100, "%.\n\n")
		filtered.geno.df2 = trunc.filtered.geno.df[, which(indi.mis.bool == 1, arr.ind = TRUE)]



		# 3. Remove genotype with MAF<5%
		if (threshold.MAF == -1) {
			threshold.MAF = 1 - (sample.size - 1)/sample.size # This threshold will only remove genotypes that have MAF=0. Otherwise, it will cause problem when calculate pval (system is exactly singular)!
		}

		cat("# 3. Remove genotype with MAF<", threshold.MAF * 100, "%.\n")
		filtered.geno.df3 = NULL
		MAF.fun <- function(col) {
			unique.allel = unique(col) # find the unique allele
			unique.allel.non0 = unique.allel[which(unique(col) != 0)] # remove zero
			if (length(unique.allel.non0) == 1) {
				return(0)
			} else {
				allele.freq = c(sum(col == unique.allel.non0[1]), sum(col == unique.allel.non0[2])) # calculate allele frequence
				return(min(allele.freq)/length(col))
			}
		}
		MAF.vec = apply(filtered.geno.df2, 2, MAF.fun)
		MAF.bool = MAF.vec >= threshold.MAF
		png("04_MAF.png", width = 500, height = 500)
		hist(MAF.vec, main = "MAF", ylab = "Number of genotypes", xlab = "MAF")
		dev.off()
		cat("\"04_MAF.png\" created.\n")
		cat(sum(MAF.bool == 1), "SNPs (", sum(MAF.bool == 1)/length(MAF.bool) * 100, "%) have MAF >=", threshold.MAF * 
			100, "%.\n\n") #=4266
		filtered.geno.df3 = cbind(filtered.geno.df[, 2], filtered.geno.df2[, which(MAF.bool == 1, arr.ind = TRUE)][, 
			1:sum(MAF.bool == 1)])
		#filtered.geno.df3 = cbind(filtered.geno.df[, 2], filtered.geno.df2[, which(MAF.bool == 1, arr.ind = TRUE)])
		
		return(filtered.geno.df3)

	}

	filtered.data$filtered.geno.df = filter.genotype.fun(phenotype.df, genotype.df, threshold.geno, threshold.indi, 
		threshold.MAF)
	write.table(filtered.data$filtered.geno.df, file = "filtered_data_genotype.txt", sep = " ") # write the filtered genotype to a txt file.
	cat("Write the filtered genotype to filtered_data_genotype.txt file.\n\n")



	# Filter: PHENOTYPES
	# Remove the phenotypes with no associated individuals.
cat("Filter data: Phenotype\n")
	cat("Remove the phenotypes with no associated individuals.\n\n")
	filter.genotype.fun <- function(ordered.pheno) {
		# Find individuals without phenotypes
		miss.geno2 <- function(row) {
			return((sum(row == filtered.data$filtered.geno.df[, 1])) != 0)
		}
		geno.mis.bool2 = apply(as.matrix(ordered.pheno[, 2]), 1, miss.geno2)
		filtered.pheno = ordered.pheno[which(geno.mis.bool2 == 1, arr.ind = TRUE), ]

		# Check phenotype distributed again
		cat("Check whether filtered phenotype is normally distributed.\n")
		png("05_Post-filter phenotype distribution.png", width = 500, height = 500)
		hist(filtered.pheno[, 3], breaks = 20, main = "Post-fitler phenotype distribution", ylab = "Frequency", xlab = "Average phenotype")
		dev.off()
		cat("\"05_Post-filter phenotype distribution.png\" created.\n\n")

		return(filtered.pheno)
	}

	filtered.data$filtered.pheno.df = filter.genotype.fun(ordered.pheno)
	write.table(filtered.data$filtered.pheno.df, file = "filtered_data_phenotype.txt", sep = " ") # write the filtered phenotype to a txt file.
	cat("Write the filtered phenotype to filtered_data_phenotype.txt file.\n\n")




	######################### Convert genotype to Xa #########################
	
	#Update sample.size and num.sample after data filtering
	size = dim(filtered.data$filtered.geno.df)
	sample.size = size[1]
	num.alleles = size[2] - 1
	num.sample = num.alleles/2

	Xa = matrix(0, sample.size, num.sample)
	for (col in 1:num.sample) {
		geno.col = filtered.data$filtered.geno.df[, (2 * col):(2 * col + 1)]
		unique.allele = unique(geno.col[1]) # find the unique allele
		unique.allele.non0 = t(unique.allele)[which(unique.allele != 0)] # remove zero
		if (length(unique.allele.non0) == 1) {
			allele.freq = c(sum(geno.col == unique.allele.non0[1]), 0)
		} else {
			allele.freq = c(sum(geno.col == unique.allele.non0[1]), sum(geno.col == unique.allele.non0[2]))
		}


		for (row in 1:sample.size) {
			if (geno.col[row, 1] == 0 | geno.col[row, 2] == 0) {
				Xa[row, col] = 10
			} else if (geno.col[row, 1] == (geno.col[row, 2])) {
				if (geno.col[row, 1] == unique.allele.non0[match(max(allele.freq), allele.freq)]) {
					Xa[row, col] = 1
				} else if (geno.col[row, 1] == unique.allele.non0[match(min(allele.freq), allele.freq)]) {
					Xa[row, col] = -1
				}
			}

		}
	}

	sample.ls = list()
	sample.ls$Xa = Xa
	sample.ls$y = filtered.data$filtered.pheno.df[, 3]





	######################### Perform a PCA to detect population structure #########################
	
	X = as.matrix(t(sample.ls$Xa))
	Y = as.matrix(sample.ls$y)

	PCA.fun <- function(X, Y) {
		#scale the genotypes by their mean and standard deviation
		W <- (X - rowMeans(X))/sqrt(diag((cov(X))))

		#perform a PCA
		geno.pc <- princomp(W)

		#plot the loadings of the first two PCs
		png("06_PCA.png", width = 1000, height = 500)
		par(mfrow = c(1, 2))
		plot(geno.pc$loadings[, c(1, 2)], main = "PCA")
		hist(geno.pc$loadings[, 1], breaks = length(Y), main = "Distribution of Comp.1", xlab = "Comp.1", ylab = "Frequency")
		dev.off()
		cat("\"06_PCA.png\" created.\n\n")
		return(geno.pc)
	}

	geno.pc = PCA.fun(X, Y)

	# Create covariate matrix Xz
	Xz = NULL
	threshold.comp1 <- 0.05 #X-axis value which appears to split group 0 from group 1
	group0 <- which(geno.pc$loadings[, 1] < threshold.comp1)
	Xz <- array(1, sample.size)
	Xz[group0] <- 0





	######################### GWAS #########################
	
	result.ls = list() # store all the results together in one list! Please refer to the top of this file for a summary and the structure of result.ls.

	#:::::::::::::::::F-statistics & p-value:::::::::::::::::
	
	Ftest.fun <- function(sample.ls.y, sample.ls.Xa) {
		#Create a list
		Fstats_pvals.ls = list()

		ones = rep(1, length(sample.ls.y))
		data = cbind(ones, sample.ls.Xa, sample.ls.y)

		#remove individuals (and their phenotypes) with no genotypes.
		geno.bool.fun <- function(row) {
			return(row[2] != 10)
		}
		geno.bool = apply(data, 1, geno.bool.fun)
		filtered.data = data[which(geno.bool == 1, arr.ind = TRUE), ]


		x = filtered.data[, 1:2]
		y = filtered.data[, 3]
		sample.size = length(y)

		#Calculate F statistic
		##Calculate y bar
ybar = sum(y)/sample.size

		##Calculate beta.MLE
		Fstats_pvals.ls$beta.MLE = solve((t(x) %*% x)) %*% t(x) %*% y

		##Calculate y hat
		yhat = x %*% Fstats_pvals.ls$beta.MLE

		## SSM, SSE, MSM, MSE
		SSM = sum((yhat - ybar)^2)
		SSE = sum((y - yhat)^2)
		df.M = 1
		df.E = sample.size - 2
		MSM = SSM/df.M
		MSE = SSE/df.E
		Fstats_pvals.ls$Fstats <- MSM/MSE


		#Calculate p value
		Fstats_pvals.ls$pval <- pf(Fstats_pvals.ls$Fstats, df1 = df.M, df2 = df.E, lower.tail = FALSE, log.p = FALSE)
		Fstats_pvals.ls$log10pval <- -log10(Fstats_pvals.ls$pval)

		#Return
		return(Fstats_pvals.ls)
	}



	LRM.fun <- function(sample.ls.y, sample.ls.Xa) {
		Ftest.ls = list()
		beta.MLE = matrix(0, 2, num.sample)
		for (i in 1:num.sample) {
			Fstats_pvals.ls <- Ftest.fun(sample.ls.y, sample.ls.Xa[, i])
			Ftest.ls$Fstats[i] = Fstats_pvals.ls$Fstats
			Ftest.ls$pval[i] = Fstats_pvals.ls$pval
			beta.MLE[, i] = Fstats_pvals.ls$beta.MLE
		}
		Ftest.ls$beta.MLE = beta.MLE
		return(Ftest.ls)
	}
	Ftest.ls = list()
	Ftest.ls = LRM.fun(sample.ls$y, sample.ls$Xa)

	result.ls$pval = Ftest.ls$pval

	#::::::::::::::::: QQ plot:::::::::::::::::
	
	log10pval.expected = -log10((1:num.sample)/num.sample)
	ranklog10pval.expected = log10pval.expected[order(log10pval.expected)]

	log10pval = -log10(result.ls$pval)
	ranklog10pval = log10pval[order(log10pval)]

	png("07.1_QQ Plot.png", width = 500, height = 500)
	plot(ranklog10pval.expected, ranklog10pval, main = "QQ Plot", xlab = "Expected", ylab = "pval")
	dev.off()
	cat("\"07.1_QQ Plot.png\" created.\n\n")

	#::::::::::::::::: Manhattan plot with Bonferroni correction:::::::::::::::::
	
	alpha = 0.05
	bonferroni = alpha/num.sample
	bonferroni.log10 = -log10(bonferroni) #=4.630021
	position = seq(1, num.sample)
	y.bonferroni = rep(bonferroni.log10, num.sample)
	png("08.1_Manhattan Plot with Bonferroni correction.png", width = 500, height = 500)
	plot(position, -log10(result.ls$pval), main = "Manhattan Plot with Bonferroni Correction", ylab = "-log(p)", 
		xlab = "position in genotype markers")
	lines(position, y.bonferroni, type = "l", col = "red")
	dev.off()
	cat("\"08.1_Manhattan Plot with Bonferroni correction.png\" created.\n\n")





	######################### GWAS with covariate #########################
	
	LRM.covariate.fun <- function(sample.ls.y, sample.ls.Xa, geno.pc.loadings) {
		# Adapted from Jason's code for Lab #9
		# Perform a GWAS with the covariate

		sample.size = nrow(sample.ls.Xa)
		pval.covariate <- NULL

		#calculate the p-values for the genotype for each marker
		for (i in 1:ncol(Xa)) {
			# Xmu
			Xmu <- rep(1, nrow(Xa))

			# Covariate
			xvalue <- 0 #X-axis value which appears to split group 0 from group 1
			group0 <- which(geno.pc.loadings[, 1] < xvalue)
			Xz <- array(1, sample.size)
			Xz[group0] <- 0

			Y = sample.ls.y
			Xa.SNP <- sample.ls.Xa[, i]
			combined.data.df <- data.frame(Y, Xmu, Xa.SNP, Xz)
			names(combined.data.df) <- c("Y", "Xmu", "Xa", "Xz")

			#remove individuals (and their phenotypes) with no genotypes.
			geno.bool.fun <- function(row) {
				return(row[3] != 10)
			}
			geno.bool = apply(combined.data.df, 1, geno.bool.fun)
			testdat.df = combined.data.df[which(geno.bool == 1, arr.ind = TRUE), ]
			Y = testdat.df$Y
			Xmu = testdat.df$Xmu
			Xa.SNP = testdat.df$Xa
			Xz = testdat.df$Xz

			#with covariate 
			L0 <- lm(Y ~ Xz, data = testdat.df)
			L1 <- lm(Y ~ Xa.SNP + Xz, data = testdat.df)
			out.aov <- anova(L0, L1)
			pval.covariate <- c(pval.covariate, out.aov$"Pr(>F)"[2])

		}

		#a quick fix for p-values equal to zero
		pval.covariate[which(pval.covariate == 0)] <- 1e-15

		# QQ plot
		log10pval.expected = -log10((1:ncol(Xa))/ncol(Xa))
		ranklog10pval.expected = log10pval.expected[order(log10pval.expected)]

		log10pval.covariate = -log10(pval.covariate)
		ranklog10pval.covariate = log10pval.covariate[order(log10pval.covariate)]


		png("07.2_QQ Plot with covariate.png", width = 500, height = 500)
		plot(ranklog10pval.expected, ranklog10pval.covariate, main = "QQ Plot with covariate", xlab = "Expected", 
			ylab = "pval.covariate")
		dev.off()
		cat("\"07.2_QQ plot with covariate.png\" created.\n\n")

		return(pval.covariate)

	}

	pval.covariate = LRM.covariate.fun(sample.ls$y, sample.ls$Xa, geno.pc$loadings)
	result.ls$pval.covariate = pval.covariate

	#::::::::::::::::: Manhattan plot with Bonferroni correction :::::::::::::::::
	
	alpha = 0.05
	bonferroni = alpha/num.sample
	bonferroni.log10 = -log10(bonferroni) #=4.630021
	position = seq(1, num.sample)
	y.bonferroni = rep(bonferroni.log10, num.sample)
	png("08.2_Manhattan Plot with covariate with Bonferroni correction.png", width = 500, height = 500)
	plot(position, -log10(result.ls$pval.covariate), main = "Manhattan Plot of pval.covariate with Bonferroni Correction", 
		ylab = "-log(p.covariate)", xlab = "position in genotype markers")
	lines(position, y.bonferroni, type = "l", col = "red")
	dev.off()
	cat("\"08.2_Manhattan Plot with covariate with Bonferroni correction.png\" created.\n\n")






	######################### Potential hits #########################
	
	#::::::::::::::::: finding hits :::::::::::::::::
	
	#List of hits
	alpha = 0.05
	bonferroni = alpha/num.sample
	bonferroni.log10 = -log10(bonferroni)

	hits.ls = list()
	h = 1
	for (p in 1:num.sample) {
		if (-log10(result.ls$pval[p]) > bonferroni.log10) {
			hits.ls$positon[h] = p
			hits.ls$pval[h] = -log10(result.ls$pval[p])
			h = h + 1
		}
	}

	result.ls$hits.ls = hits.ls


	#::::::::::::::::: Zoomed in Manhattan plot :::::::::::::::::
	hits.index = c()
	interval = 25
	log10pval = -log10(result.ls$pval)
	hits.index[1] = which(log10pval == max(log10pval))
	y.bonferroni = rep(bonferroni.log10, interval * 2 + 1)

	png("09_Zoomed in manhattan Plot with Bonferroni correction.png", width = 500 * length(hits.index), height = 500)
	par(mfrow = c(length(hits.index), 1))
	for (i in 1:length(hits.index)) {
		ind = hits.index[i]
		position = seq(ind - interval, ind + interval)
		plot(position, -log10(result.ls$pval)[(ind - interval):(ind + interval)], main = "Zoomed in manhattan Plot with Bonferroni correction", 
			xlab = "Position", ylab = "-log(p)")
		lines(position, y.bonferroni, type = "l", col = "red")
	}
	dev.off()


	#::::::::::::::::: Find hit's information :::::::::::::::::
	map = read.table("project_information.txt")
	cat("Find and store chromosome and SNP information for significant hits into result.ls$hits.info.\n")
	map.vec = c()
	for (i in 1:length(hits.index)) {
		col.name = colnames(filtered.data$filtered.geno.df)[hits.index[i] * 2]
		col.num = as.numeric(substring(col.name, 2))
		map.row = (col.num + 1)/2
		map.vec = rbind(map.vec, map[map.row, -3])
	}
	map.vec = map.vec[, c(1, 3, 2)]
	names(map.vec) <- c("Chr", "SNP", "Chr.SNP")
	result.ls$hits.info = map.vec

	return(result.ls)
}
