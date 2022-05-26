# ligera 1.0.0.9000 (2020-06-26)

* First GitHub release!  Includes `ligera` function and `scripts/ligera.R` command line interface script, which supports plink BED/BIM/FAM files.

# ligera 1.0.1.9000 (2020-07-09)

* Added `ligera2`, which implements full "BOLT trick" to avoid having to calculate kinship explicitly, and which scales as $O( m n \sqrt{\kappa} )$ instead of $O( mn^2 + n^2 \sqrt{\kappa} )$, which is far better as n increases.

# ligera 1.0.2.9000 (2020-07-13)

* Added `ligera2_bed`, a version of `ligera2` specific to plink BED files sped up with Rcpp code for critical components!
* Added covariate support to all the main functions (`ligera`, `ligera2`, and `ligera2_bed`).

# ligera 1.0.3.9000 (2020-07-14)

* Added multiscan wrapper for all versions of LIGERA (`ligera_multi`, `ligera2_multi`, `ligera2_bed_multi`), which iterates the genetic association scans and adds significant loci to covariates (options for one per iteration or all significant per iteration) to condition upon them and increase power on the weaker signals.

# ligera 1.0.4.9000 (2020-07-21)

* Added support for missing values in covariates matrix
  * In each column of covariates matrix, missing values are replaced with the within-column mean value.
  * All core versions updated (`ligera`, `ligera2`, and `ligera2_bed`).
  * This was necessary for multiscan versions (`ligera_multi`, `ligera2_multi`, `ligera2_bed_multi`) to work when the genotype matrix has missing values, as these may be incorporated as covariates in the iteration.
* Multiscan versions stop iterations if and when the number of covariates exceed the number of individuals, to avoid inverting singular matrices (problem becomes overdetermined).
* Fixed `class` usage now that matrices return a two-element array in R version 4.

# ligera 1.0.5.9000 (2020-10-16)

* All 6 ligera versions (`ligera`, `ligera2`, `ligera2_bed`, and `*_multi` variants)
  - Added option `m_chunk_max` to further reduce memory usage (compared to previous default) without sacrificing speed.
  - Added script wrappers (located in `scripts/`) for main 3 versions of ligera, each with a `--multi` option to also cover the 3 multiscan versions (actually added 2020-09-29)
  - Added sample files to go with scripts (also located in `scripts/`; also added 2020-09-29)

# ligera 1.0.6.9000 (2020-12-08)

* Minor internal improvement (in function `conj_grad_scan_bed_wcpp`) affecting `ligera2_bed` and `ligera2_bed_multi` versions only.

# ligera 1.0.6.9000 (2021-01-08)

- Internal improvements: replaced `sweep` cases with `matrix` version (redid Tiffany Tu devel branch commits e9629dc and 8493282 with minor edits, thanks Tiffany!!!)
- Forgot to update version.

# ligera 1.0.7.9000 (2022-05-26)

- Added functions `ligera_f`, `ligera_f_multi`
  - This version calculates p-values using an F-test, which gives calibrated statistics under both quantitative and binary traits.
  - Compared to the original `ligera`, which uses the faster Wald test (calibrated for quantitative but not binary traits), the F-test version is quite a bit slower, and is optimized for `m >> n`, so it is a work in progress.
- Unit tests: fixed issues arising from apearance of fixed loci when individuals with missing traits are excluded
- Reformatted this `NEWS.md` slightly to improve its automatic parsing.
- Removed "LazyData: true" from DESCRIPTION (to avoid a new "NOTE" on CRAN).

# ligera 1.0.8.9000 (2022-05-26)

- Unit tests: reduced simulation size (`n = 20`, from 100; left `m = 1000`)
  - Reduced its runtime to 3.4s (before it was highly variable but excessive, between 60s-131s).
  - Reduction does not appear to increase testing artifacts (such as singular matrix inversions).
