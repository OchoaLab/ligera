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

# ligera 1.0.9.9000 (2022-05-26)

- Function `ligera` fixed bug that `inbr` wasn't subset for non-NA individuals when there is missingness in the input trait
  - Only occurred if `inbr` was explicitly provided, otherwise lazy evaluation from kinship (which was subset) ensured that it was subset when it was not provided as explicit parameter.
  - Other versions (`ligera2`, `ligera2_bed`, `ligera_f`) do not have the option of providing `inbr` so were not affected by similar bugs.
- Function `ligera_f` 
  - Added option `V` to specify algorithm (there are four versions, numbered 0-3, roughly falling into two scalability classes).  This is a temporary option that will be removed if one version emerges as clearly superior.
  - Default algorithm changed to what is now `V=0` (default was `V=1` before, but internal tests show that it has the lowest numerical stability of all the choices, below machine precision).  Previously all versions were commented out except `V=1`.
- Unit tests
  - Replaced trivial true kinship of unstructured simulations (`I/2`) with popkin estimates in all tests.  This was necessary so WG-biased estimates would differ from their unbiased counterparts, and helped identify the numeric stability issues of the previous `ligera_f` version.
  - Added tests that confirm that `ligera_f` is invariant to WG bias (a type of bias that always results in non-singular kinship matrices).
    - Oddly, `ligera` is not invariant to WG, and it is unclear at the moment why.

# ligera 1.0.10.9000 (2022-09-18)

- Added function `ligera2_f`, a version of `ligera_f` that calculates performs all calculations without an explicit kinship matrix estimate (like `ligera2` was to `ligera`).
  - Retains four `V` versions as before, and (except for `V=1`) all and agree with each other and with `ligera_f` up to machine precision.  Runtime has not been evaluated, might not be good yet for any version (f versions remain experimental).
- Internal function changes:
  - `popkin_prod` added option `want_inbr` (to skip calculation of `inbr`, as it is unused in f versions; default is to calculate and return it as before).
  - `conj_grad_scan` added options `want_inbr` (to skip calculation of `inbr`), `want_b` (to return `b`, default does not as before), and `b` (allows for providing `b` as calculated in a previous run).  Additionally `mean_kinship` is now `NA` by default and it does not have to be provided if `b` is provided.  Lastly, if `want_inbr` and `want_b` are both false, then only `Z` is returned (instead of a list with these objects and `Z`).
- Corrected minor errors in non-f unit tests (the wrong functions or objects were being tested, though the changes were inconsequential).

# ligera 1.0.11.9000 (2022-09-18)

- Functions `ligera_f`, `ligera2_f`: removed previous `V=1` algorithm, which was numerically unstable.  Shifted down previous versions `2:3` to `1:2`.

# ligera 1.0.12.9000 (2022-09-18)

- Added function `ligera2_f_multi`.
- Internal function changes:
  - `conj_grad_scan` and `cgsolve_mat`: handled a case where input `Y` has a column of zeroes, which should give an output `Z` with zeroes in the same columns but which the conjugate gradient algorithm gave NaN's for those same columns of `Z` (because of a 0/0 factor).  This special case is now recognized and handled correctly.

# ligera 1.0.13.9000 (2022-09-19)

- Added function `ligera2_bed_f`.
