# 2020-06-26 - ligera 1.0.0.9000

* First GitHub release!  Includes `ligera` function and `scripts/ligera.R` command line interface script, which supports plink BED/BIM/FAM files.

# 2020-07-09 - ligera 1.0.1.9000

* Added `ligera2`, which implements full "BOLT trick" to avoid having to calculate kinship explicitly, and which scales as $O( m n \sqrt{\kappa} )$ instead of $O( mn^2 + n^2 \sqrt{\kappa} )$, which is far better as n increases.
