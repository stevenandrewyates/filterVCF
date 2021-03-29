
##################################################
######            Details                   ######
##################################################

# author:   Steven Yates
# contact:  steven.yates@usys.ethz.ch
# Year:     2021
# Citation: TBA

##################################################
######            Description               ######
##################################################

# A R script for filtering a vcf file

# this script requires the library "vcfR" to be 
# installed and available!

# https://cran.r-project.org/web/packages/vcfR/index.html

# This script will filter a vcf file based on Minor
# Alelle Frequency and minimum number of individuals.
# When genotyping it filters for read depth in the 
# following way. If a genotype is homozygous it must have a 
# depth > 7. The reason being, that seeing only the 
# same read eight times is > 99% confidence it is homozygous.
# For a heterozygous individual the depth must be 
# greater than two. The reason being, that if both 
# alleles are present then it is safe to assume it 
# is heterozygous. Perhaps you disagree? Maybe you 
# want to filter at a depth greater than ten. If 
# so please see the other script.

# https://github.com/stevenandrewyates/filerVCF/filerVCF10.R

##################################################
######              Usage                   ######
##################################################

# this script takes four inputs:
# 1) an input file (input.vcf)
# 2) an output file (output.vcf)
# 3) the minor allele frequency (as a %)
# 4) the minimum number of genotypes (integar)

# R --vanilla --slave "--args input.vcf output.vcf 5 50" < filterVCF.R

##################################################
######              Script                  ######
##################################################

# load the library
library(vcfR)

# read in the arguments 
args <- commandArgs(trailingOnly = T)

# the INPUT file
INFILE <- args[1]
# OUTPUT file
OUTFILE <- args[2]
# 

# read ijn the samtools output
vcf <- read.vcfR( INFILE, verbose = FALSE )
# check the input
head(vcf@gt)
# extract the depth information
dp <- extract.gt(vcf, element="DPR", as.numeric = FALSE)
dp <- as.data.frame(dp)
# get the positional data
FIX <- vcf@fix
# get the genotypes
gt <- extract.gt(vcf, element="GT", as.numeric = FALSE)

# the minor allele frequency
MAF <- args[3]
#MAF <- 20

# lowest number of genotypes
MG <- args[4]
#MG <- 100


MAF <- as.numeric(as.character(MAF))
MG <- as.numeric(as.character(MG))

OUTPUT <- NULL

for(x in 1:(dim(dp)[1]))
#for(x in 1:100)
{
print(paste("Row...",x))
# make a table of the genotype calls
TAB <- table(gt[x,])
print(sum(TAB))
# if the number of genotype calls is less than the jmber of genotypes then skip
if(sum(TAB) < MG) {next}
# check that the allele is polymorphic
print(length(TAB))
if (length(TAB) < 2) {next}
# now check that the variants are above the MAF
TAB <- TAB[TAB > (100/sum(TAB))*MAF]
if (length(TAB) < 2) {next}
print("found variant")

#  now check that the coverage matches
	LINE <- t(dp[x,])
	COVERAGE <- rep(0,length(LINE))
	GTO <- rep("./.",length(LINE))
	for (y in 1:length(LINE))
		{
#		print(paste(gt[x,y],"...",dp[x,y],"...",sum(as.numeric(strsplit(LINE[y],",")[[1]]))))
		COVERAGE[y] <- sum(as.numeric(strsplit(LINE[y],",")[[1]]))
		if(is.na(gt[x,y])) {next}
		GIN <- strsplit(as.character(gt[x,y]),"/")[[1]]
	# make a call if homozygous
		if(GIN[1] == GIN[2]) 
			{
			if(COVERAGE[y] > 7)
				{
				GTO[y] <- gt[x,y]
				}
			}
	
	# make a call if heterozygous
		if(GIN[1] != GIN[2]) 
			{
			if(COVERAGE[y] > 1)
				{
				GTO[y] <- gt[x,y]
				}
			}
		}
# repeat the first section using the recalled data
TAB <- table(GTO)
# remove the missing value
TAB <- TAB[names(TAB) != "./."]
print(sum(TAB))
# if the number of genotype calls is less than the jmber of genotypes then skip
if(sum(TAB) < MG) {next}
# check that the allele is polymorphic
print(length(TAB))
if (length(TAB) < 2) {next}
# now check that the variants are above the MAF
ODD <- names(TAB)[TAB <= (100/sum(TAB))*MAF]
TAB <- TAB[TAB > (100/sum(TAB))*MAF]
# make sure it is still polymorphic
if (length(TAB) < 2) {next}
# replace the low frequency alleles with unknowns
GTO[GTO == ODD] <- "./."
print(paste(x, "meets criteria"))
# combine the data
DR <- (c(FIX[x,],paste(GTO,LINE,sep=":")))
#save the data to the output
OUTPUT <- rbind(OUTPUT,DR)
}

OUTCOL <- length(colnames(OUTPUT))
OUTDP <- length(colnames(dp))
colnames(OUTPUT)[OUTCOL-OUTDP:OUTCOL]


colnames(OUTPUT)[(OUTCOL-OUTDP+1):OUTCOL] <- colnames(dp)

VCFOUT <- data.frame(OUTPUT[,1:8],FORMAT="GT:DPR",OUTPUT[,9:OUTCOL])


colnames(VCFOUT)[1] <- "##fileformat=VCFv\n#CHROM"
write.table(VCFOUT,file=OUTFILE,sep="\t",col.names=T,quote=F,row.names=F)



