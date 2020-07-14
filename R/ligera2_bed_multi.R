#' LIGERA2 BED multi: LIght GEnetic Robust Association multiscan function
#'
#' This function performs multiple genetic association scans, adding one significant locus per iteration to the model (modeled as a covariate) to increase power in the final model.
#' The function returns a tibble containing association statistics and several intermediates.
#' This optimized version requires the genotypes to be in a file in BED format.
#' 
#' Suppose there are `n` individuals and `m` loci.
#'
#' @param file The path to the BED file containing the genotypes, potentially excluding the BED extension.
#' @param m_loci The number of loci in the BED file.
#' @param n_ind The number of individuals in the BED file.
#' @param trait The length-`n` trait vector, which may be real valued and contain missing values.
#' @param mean_kinship An estimate of the mean kinship produced externally, to ensure internal estimates of kinship and inbreeding are unbiased.
#' @param q_cut The q-value threshold to admit new loci into the polygenic model.
#' @param covar An optional `n`-by-`K` matrix of `K` covariates, aligned with the individuals.
#' @param mem_factor Proportion of available memory to use loading and processing genotypes.
#' Ignored if `mem_lim` is not `NA`.
#' @param mem_lim Memory limit in GB, used to break up genotype data into chunks for very large datasets.
#' Note memory usage is somewhat underestimated and is not controlled strictly.
#' Default in Linux and Windows is `mem_factor` times the free system memory, otherwise it is 1GB (OSX and other systems).
#' @param tol Tolerance value passed to conjugate gradient method solver.
#'
#' @return A tibble containing the following association statistics from the last scan for non-selected loci.
#' For selected loci, these are the values from the scan before each was added to the model (as after addition they get `beta ~= 0` and `pval ~= 1`).
#' 
#' - `pval`: The p-value of the last association scan.
#' - `beta`: The estimated effect size coefficient for the trait vector at this locus.
#' - `beta_std_dev`: The estimated coefficient variance of this locus (varies due to dependence on minor allele frequency).
#' - `p_q`: The allele variance estimate (estimate of `p*(1-p)`).  The number of heterozygotes, weighted by inbreeding coefficient, and with pseudocounts included, is used in this estimate (in other words, it does not equal MAF * ( 1 - MAF ), where MAF is the marginal allele frequency.
#' - `t_stat`: The test statistic, equal to `beta / beta_std_dev`.
#' - `qval`: The q-value of the last association scan.
#' - `sel`: the order in which loci were selected, or zero if they were not selected.
#'
#' @examples
#' # MISSING SAMPLE BED FILE
#'
#' @seealso
#' The `popkin` package.
#' 
#' @export
ligera2_bed_multi <- function(
                              file,
                              m_loci,
                              n_ind,
                              trait,
                              mean_kinship,
                              q_cut = 0.05,
                              covar = NULL,
                              mem_factor = 0.7,
                              mem_lim = NA,
                              tol = 1e-15
                              ) {
    # things to initialize for loop
    # indexes of selected loci, to remember
    loci_selected <- c()
    tib_selected <- tibble::tibble() # ok to be blank
    # were there new selections?  To know if we've converged
    # initialize to TRUE so first iteration occurs
    new_selec <- TRUE
    
    # load genotypes using BEDMatrix, for loading selected loci as covariates
    X <- BEDMatrix::BEDMatrix( file, n = n_ind, p = m_loci )
    
    # loop while there are still significant loci
    while( new_selec ) {
        # assume that covar has been grown to include selected loci already
        
        # run ligera
        # NOTE: selected loci become insigificant when retested
        tib <- ligera2_bed(
            file = file,
            m_loci = m_loci,
            n_ind = n_ind,
            trait = trait,
            mean_kinship = mean_kinship,
            covar = covar,
            mem_factor = mem_factor,
            mem_lim = mem_lim,
            tol = tol
        )
        
        # add q-values
        # only to non-selected loci (selected loci all have p=1, this is terrible for q-value estimation!)
        # unfortunately this doesn't work if loci_selected has length zero
        if ( length( loci_selected ) == 0 ) {
            tib$qval <- qvalue::qvalue( tib$pval )$qvalues
        } else {
            # have to initialize qval here before doing this
            tib$qval <- 1 # all selected loci get 1
            tib$qval[ -loci_selected ] <- qvalue::qvalue( tib$pval[ -loci_selected ] )$qvalues
        }

        # find the most significant locus
        # which.min always returns a scalar (*first* minimum)
        index <- which.min( tib$pval ) # p-values may have more resolution, q-values sometimes tie, otherwise the same as q-values

        # if this is significant, add to covariates and reiterate
        if ( tib$qval[ index ] < q_cut ) {

            # add to list of selected loci, for our info
            loci_selected <- c( loci_selected, index )

            # remember their stats from their last round before getting added as covariates
            tib_selected <- dplyr::bind_rows(
                                       tib_selected,
                                       tib[ index, ]
                                   )

            # add to covariates
            # first extract genotype vector
            X_i <- X[ , index ] # BEDMatrix orientation
            # add to covariates matrix now
            covar <-
                if ( is.null( covar ) ) {
                    # turn into a matrix with a single column
                    cbind( X_i )
                } else {
                    cbind( covar, X_i )
                }

            # `new_selec` stays `TRUE`
        } else {
            # we're done discovering new loci :(
            new_selec <- FALSE
        }
    }

    # it'd be nice to return the last tibble, but again the selected loci will all have p=1
    # instead we should record the last p-value and other stats they had prior to being added to the model
    # first we add a new column that marks selected loci, in order
    # all non-selected loci get a zero
    tib$sel <- 0
    # when there were no selections at all, this fails, so have to test for that
    if ( length( loci_selected ) > 0 ) {
        # this sets the right values for the selected loci (they are in the order in which they were selected)
        tib_selected$sel <- 1 : nrow( tib_selected )
        # this overwrites as desired
        tib[ loci_selected, ] <- tib_selected
    }
    
    # done, return tibble!
    return( tib )
}
