# filterVCF
A R script for filtering a vcf file

# Description

This script requires the library "vcfR" to be installed and available!
https://cran.r-project.org/web/packages/vcfR/index.html

This script will filter a vcf file based on Minor Alelle Frequency and minimum number of individuals. When genotyping it filters for read depth in the following way. If a genotype is homozygous it must have a depth > 7. The reason being, that seeing only the same read eight times is > 99% confidence it is homozygous. For a heterozygous individual the depth must be greater than two. The reason being, that if both alleles are present then it is safe to assume it is heterozygous. Perhaps you disagree? Maybe you want to filter at a depth greater than ten. If so please see the other script.

https://github.com/stevenandrewyates/filerVCF/filerVCF10.R

# Usage

This script takes four inputs:
1) an input file (input.vcf)
2) an output file (output.vcf)
3) the minor allele frequency (as a %)
4) the minimum number of genotypes (integar)

```
R --vanilla --slave "--args input.vcf output.vcf 5 50" < filterVCF.R
```

or for the filterVCF10.R use:

```
R --vanilla --slave "--args input.vcf output.vcf 5 50" < filterVCF10.R
```
