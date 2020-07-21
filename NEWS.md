# 2020-06-26 - ligera 1.0.0.9000

* First GitHub release!  Includes `ligera` function and `scripts/ligera.R` command line interface script, which supports plink BED/BIM/FAM files.

# 2020-07-09 - ligera 1.0.1.9000

* Added `ligera2`, which implements full "BOLT trick" to avoid having to calculate kinship explicitly, and which scales as $O( m n \sqrt{\kappa} )$ instead of $O( mn^2 + n^2 \sqrt{\kappa} )$, which is far better as n increases.

# 2020-07-13 - ligera 1.0.2.9000

* Added `ligera2_bed`, a version of `ligera2` specific to plink BED files sped up with Rcpp code for critical components!
* Added covariate support to all the main functions (`ligera`, `ligera2`, and `ligera2_bed`).

# 2020-07-14 - ligera 1.0.3.9000

* Added multiscan wrapper for all versions of LIGERA (`ligera_multi`, `ligera2_multi`, `ligera2_bed_multi`), which iterates the genetic association scans and adds significant loci to covariates (options for one per iteration or all significant per iteration) to condition upon them and increase power on the weaker signals.

# 2020-07-21 - ligera 1.0.4.9000

* Added support for missing values in covariates matrix
  * In each column of covariates matrix, missing values are replaced with the within-column mean value.
  * All core versions updated (`ligera`, `ligera2`, and `ligera2_bed`).
  * This was necessary for multiscan versions (`ligera_multi`, `ligera2_multi`, `ligera2_bed_multi`) to work when the genotype matrix has missing values, as these may be incorporated as covariates in the iteration.
* Multiscan versions stop iterations if and when the number of covariates exceed the number of individuals, to avoid inverting singular matrices (problem becomes overdetermined).
* Fixed `class` usage now that matrices return a two-element array in R version 4.
