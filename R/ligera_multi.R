#' LIGERA_multi: LIght GEnetic Robust Association multiscan function
#'
#' This function performs multiple genetic association scans, adding one significant locus per iteration to the model (modeled as a covariate) to increase power in the final model.
#' The function returns a tibble containing association statistics and several intermediates.
#' 
#' Suppose there are `n` individuals and `m` loci.
#'
#' @param X The `m`-by-`n` genotype matrix, containing dosage values in (0, 1, 2, NA) for the reference allele at each locus.
#' @param trait The length-`n` trait vector, which may be real valued and contain missing values.
#' @param kinship The `n`-by-`n` kinship matrix, estimated by other methods (i.e. the `popkin` package).
#' @param q_cut The q-value threshold to admit new loci into the polygenic model.
#' @param one_per_iter If true, only the most significant locus per iteration is added to model of next iteration.  Otherwise all significant loci per iteration are added to the model of next iteration.
#' @param kinship_inv The optional matrix inverse of the kinship matrix.  Setting this parameter is not recommended, as internally a conjugate gradient method (`\link[cPCG]{cgsolve}`) is used to implicitly invert this matrix, which is much faster.  However, for very large numbers of traits without missingness and the same kinship matrix, inverting once might be faster.
#' @param covar An optional `n`-by-`K` matrix of `K` covariates, aligned with the individuals.
#' @param loci_on_cols If `TRUE`, `X` has loci on columns and individuals on rows; if false (the default), loci are on rows and individuals on columns.
#' If `X` is a BEDMatrix object, `loci_on_cols = TRUE` is set automatically.
#' @param mem_factor Proportion of available memory to use loading and processing genotypes.
#' Ignored if `mem_lim` is not `NA`.
#' @param mem_lim Memory limit in GB, used to break up genotype data into chunks for very large datasets.
#' Note memory usage is somewhat underestimated and is not controlled strictly.
#' Default in Linux and Windows is `mem_factor` times the free system memory, otherwise it is 1GB (OSX and other systems).
#' @param tol Tolerance value passed to `\link[cPCG]{cgsolve}`.
#' @param maxIter Maximum number of iterations passed to `\link[cPCG]{cgsolve}`.
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
#' # Construct random data
#' # number of individuals we want
#' n_ind <- 5
#' # number of loci we want
#' m_loci <- 100
#' # a not so small random genotype matrix
#' X <- matrix(
#'     rbinom( m_loci * n_ind, 2, 0.5 ),
#'     nrow = m_loci
#' )
#' # random trait
#' trait <- rnorm( n_ind )
#' # add a genetic effect from first locus
#' trait <- trait + X[ 1, ]
#' # kinship matrix
#' kinship <- diag( n_ind ) / 2 # unstructured case
#'
#' tib <- ligera_multi( X, trait, kinship )
#' tib
#'
#' @seealso
#' The `popkin` and `cPCG` packages.
#' 
#' @export
ligera_multi <- function(
                         X,
                         trait,
                         kinship,
                         q_cut = 0.05,
                         one_per_iter = FALSE,
                         kinship_inv = NULL,
                         covar = NULL,
                         loci_on_cols = FALSE,
                         mem_factor = 0.7,
                         mem_lim = NA,
                         # cgsolve options
                         tol = 1e-15, # default 1e-6
                         maxIter = 1e6 # default 1e3
                         ) {
    # things to initialize for loop
    # indexes of selected loci, to remember
    loci_selected <- c()
    tib_selected <- tibble::tibble() # ok to be blank
    # were there new selections?  To know if we've converged
    # initialize to TRUE so first iteration occurs
    new_selec <- TRUE
    
    # need to have this here
    if ('BEDMatrix' %in% class(X))
        loci_on_cols <- TRUE

    # loop while there are still significant loci
    while( new_selec ) {
        # assume that covar has been grown to include selected loci already
        
        # run ligera
        # NOTE: selected loci become insigificant when retested
        tib <- ligera(
            X = X,
            trait = trait,
            kinship = kinship,
            kinship_inv = kinship_inv,
            covar = covar,
            loci_on_cols = loci_on_cols,
            mem_factor = mem_factor,
            mem_lim = mem_lim,
            tol = tol,
            maxIter = maxIter
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

        # decide what to add, if anything
        if ( one_per_iter ) {
            # find the most significant locus
            # which.min always returns a scalar (*first* minimum)
            indexes <- which.min( tib$pval ) # p-values may have more resolution, q-values sometimes tie, otherwise the same as q-values
            # tell the next part of the loop if this was significant (overwrite new_select from before, really updates it)
            new_selec <- tib$qval[ indexes ] < q_cut
        } else {
            # directly select subset that was newly significant on this round
            indexes <- which( tib$qval < q_cut )
            # here we move to the next stage as long as indexes was not empty
            new_selec <- length( indexes ) > 0
        }

        # if there were any newly significant loci, add to covariates and reiterate
        if ( new_selec ) {

            # add to list of selected loci
            loci_selected <- c( loci_selected, indexes )

            # remember their stats from their last round before getting added as covariates
            tib_selected <- dplyr::bind_rows(
                                       tib_selected,
                                       tib[ indexes, ]
                                   )

            # add to covariates
            # first extract genotype submatrix
            # NOTE this is matrix even if there was a single index there
            # here transposition is backwards than usual (individuals in X are usually on the column, but as covariates we need then on the rows!)
            X_i <- if ( loci_on_cols ) X[ , indexes, drop = FALSE ] else t( X[ indexes, , drop = FALSE ] )
            
            # add to covariates matrix now
            covar <- if ( is.null( covar ) ) X_i else cbind( covar, X_i )
            
            # `new_selec` stays `TRUE`
        }
        # otherwise `new_selec == FALSE`, so we're done
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
