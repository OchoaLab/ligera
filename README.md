# ligera

The goal of LIGERA (LIght GEnetic Robust Association) is to perform extremely fast genetic association studies, while fully modeling relatedness via an estimated population kinship matrix, and being robust to unmodeled environmental effects by employing a reverse regression model (where the genotype, rather than the trait, is the response variable).

## Installation

<!-- You can install the released version of ligera from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->
<!-- install.packages("ligera") -->
<!-- ``` -->

Install the latest development version from GitHub:
```R
install.packages("devtools") # if needed
library(devtools)
install_github("OchoaLab/ligera", build_opts = c())
```

## Example

This pipeline works if the input data are in plink format (BED/BIM/FAM files)

``` r
# this package, to perform the genetic association study
library(ligera)
# the package that estimates the population kinship matrix
library(popkin)
# a package that reads genotype matrices using low memory
library(BEDMatrix)
# a package that reads the BIM (variant annotations) and FAM (individual annotations) files, and reads and writes other genetics files
library(genio)
# for manipulating data.frames/tibbles
library(dplyr)

# write here the base name in common for the plink BED, BIM, and FAM files
name <- 'my-genotype-data'
# create a BEDMatrix object, which allows random access to the genotypes
X <- BEDMatrix( name )
# read the full BIM and FAM tables too
bim <- read_bim( name )
fam <- read_fam( name )

# if the trait of interest is in the fam file, extract it now
trait <- fam$pheno
# otherwise, if the phenotype is in a separate PHEN file, load that
phen <- read_phen( name )
# (reorder phen to match fam, if needed)
stopifnot( phen$id == fam$id )
trait <- phen$pheno

# to estimate the kinship matrix, it's best to group individuals into subpopulations (so the minimum value can be estimated through averaging).
# here we assume that these labels are in the first column of the FAM file
labs <- fam$fam
# estimate the kinship matrix
# (if there are no useful subpopulation labels, `labs` can be omitted)
kinship <- popkin( X, labs )
# NOTE: kinship estimates should always be inspected for correctness
# (visualize using `plot_popkin` from package `popkin`)

# now we have all the parts we need, run LIGERA!
tib <- ligera( X, trait, kinship )
# association p-values
tib$pval
# can merge BIM into table to have a more complete report
tib <- bind_cols( bim, tib )
```

## Command-line script

A script with an interphase just like GCTA's and other related genetic association software is available in the `scripts/` subdirectory of this repository (on GitHub).
The prerequisite R packages are:

- ligera
- popkin
- BEDMatrix
- genio
- dplyr
- readr
- optparse

The options can be seen by calling the script with the help `-h` option:
```bash
Rscript ligera.R -h
```
```
Usage: ligera.R [options]

Options:
	--bfile=CHARACTER
		Base name for input plink files (.BED/BIM/FAM)

	--pheno=CHARACTER
		Base name for phenotype file (.PHEN).  If missing, phenotype in FAM table is used

	--out=CHARACTER
		Base name for report output file (default same as --bfile)

	-h, --help
		Show this help message and exit
```

So if the input files are `name_in.bed`, `name_in.bim`, `name_in.fam`, and `name_pheno.phen`, then the script runs this way:
```bash
Rscript ligera.R --bfile name_in --pheno name_pheno --out name_out
```
The script essentially performs two steps, just as in the earlier example code above:

- First, it estimates the kinship matrix with `popkin`, using the first column of `name_in.fam` as subpopulation labels (will stop if they are all one value).
  For later inspection, the kinship estimate is saved in GCTA's binary GRM format, creating files `name_out.grm.bim` and `name_out.grm.id`.
  These can be loaded back into R using `read_grm` from the package `genio`, and visualized using `plot_popkin` from the package `popkin`.
- Second, it runs `ligera` and writes the statistics table to `name_out.txt`.
