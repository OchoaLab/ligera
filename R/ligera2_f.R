#' LIGERA2: LIght GEnetic Robust Association main function
#'
#' This function performs the genetic association tests on every locus of a genotype matrix against a quantitative trait, implicitly computing the kinship matrix in a way that scales better than an explicit kinship estimate.
#' The function returns a tibble containing association statistics and several intermediates.
#' This version calculates p-values using an F-test, which gives calibrated statistics under both quantitative and binary traits.
#' Compared to [ligera2()], which uses the faster Wald test (calibrated for quantitative but not binary traits), this F-test version is quite a bit slower, and is optimized for `m >> n`, so it is a work in progress.
#' 
#' Suppose there are `n` individuals and `m` loci.
#'
#' @param X The `m`-by-`n` genotype matrix, containing dosage values in (0, 1, 2, NA) for the reference allele at each locus.
#' @param trait The length-`n` trait vector, which may be real valued and contain missing values.
#' @param mean_kinship An estimate of the mean kinship produced externally, to ensure internal estimates of kinship are unbiased.
#' @param covar An optional `n`-by-`K` matrix of `K` covariates, aligned with the individuals.
#' @param loci_on_cols If `TRUE`, `X` has loci on columns and individuals on rows; if false (the default), loci are on rows and individuals on columns.
#' If `X` is a BEDMatrix object, `loci_on_cols = TRUE` is set automatically.
#' @param mem_factor Proportion of available memory to use loading and processing genotypes.
#' Ignored if `mem_lim` is not `NA`.
#' @param mem_lim Memory limit in GB, used to break up genotype data into chunks for very large datasets.
#' Note memory usage is somewhat underestimated and is not controlled strictly.
#' Default in Linux and Windows is `mem_factor` times the free system memory, otherwise it is 1GB (OSX and other systems).
#' @param m_chunk_max Sets the maximum number of loci to process at the time.
#' Actual number of loci loaded may be lower if memory is limiting.
#' @param V Algorithm version (0, 1, 2).
#' Experimental features, not worth explaining.
#' @param tol Tolerance value passed to conjugate gradient method solver.
#'
#' @return A tibble containing the following association statistics
#' 
#' - `pval`: The p-value of the association test
#' - `beta`: The estimated effect size coefficient for the trait vector at this locus
#' - `f_stat`: The F statistic
#' - `df`: degrees of freedom: number of non-missing individuals minus number of parameters of full model
#'
#' @examples
#' # Construct toy data
#' # genotype matrix
#' X <- matrix(
#'     c(0, 1, 2,
#'       1, 0, 1,
#'       1, 0, 2),
#'     nrow = 3,
#'     byrow = TRUE
#' )
#' trait <- 1 : 3
#' mean_kinship <- mean( diag( 3 ) / 2 ) # unstructured case
#'
#' tib <- ligera2_f( X, trait, mean_kinship )
#' tib
#'
#' @seealso
#' The `popkin` package.
#' 
#' @export
ligera2_f <- function(
                      X,
                      trait,
                      mean_kinship,
                      covar = NULL,
                      loci_on_cols = FALSE,
                      mem_factor = 0.7,
                      mem_lim = NA,
                      m_chunk_max = 1000,
                      V = 0,
                      tol = 1e-15
                      ) {
    # - supports missingness in trait
    # TODO
    # - support true missingness in genotypes (inversion of matrix subsets, etc; right now only approximate)
    
    # informative errors when things are missing
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    if ( missing( trait ) )
        stop( '`trait` is required!' )
    if ( missing( mean_kinship ) )
        stop( '`mean_kinship` is required!' )

    # validate mean_kinship, must be numeric scalar
    if ( !is.numeric( mean_kinship ) )
        stop( '`mean_kinship` must be numeric!' )
    if ( length( mean_kinship ) != 1 )
        stop( '`mean_kinship` must be scalar!  Got length "', length( mean_kinship ), '"' )

    # NOTE: trait and X may have missing values
    
    # override this for BEDMatrix
    if ( 'BEDMatrix' %in% class(X) ) {
        loci_on_cols <- TRUE
    } else if (!is.matrix(X))
        stop('X has unsupported class: ', toString( class( X ) ) )
    
    # need these dimensions
    if (loci_on_cols) {
        m_loci <- ncol(X)
        n_ind <- nrow(X)
    } else {
        m_loci <- nrow(X)
        n_ind <- ncol(X)
    }

    # check dimensions of other items
    if ( length( trait ) != n_ind )
        stop('Number of individuals in `trait` (', length( trait ), ') does not match genotype matrix (', n_ind , ')')
    if ( !is.null( covar ) ) {
        if ( nrow( covar ) != n_ind )
            stop('Number of individuals in `covar` (', nrow( covar ), ') does not match genotype matrix (', n_ind , ')')
    }
    
    # handle missing values in trait
    # the good thing is that this is shared across loci, so comparably it's not so expensive
    # this NULL means there are no filters to apply
    indexes_ind <- NULL
    if ( anyNA( trait ) ) {
        # indexes to keep (need to subset genotypes at load time)
        indexes_ind <- !is.na( trait )
        # subset trait
        trait <- trait[ indexes_ind ]
        # subset covariates, if present
        if ( !is.null( covar ) )
            covar <- covar[ indexes_ind, ]
        # reduce number of individuals, used in some calculations
        n_ind <- length( trait )
        # NOTE: only genotypes are left to filter with indexes_ind
    }
    
    # gather matrix of trait, intercept, and optional covariates
    Y1 <- cbind( trait, 1 )
    # add covariates, if present
    if ( !is.null( covar ) ) {
        # handle NAs now, so final Y has no missingness whatsoever
        covar <- covar_fix_na( covar )
        Y1 <- cbind( Y1, covar )
    }
    # calculate b first time only, remember after that
    obj <- conj_grad_scan(
        X = X,
        Y = Y1,
        mean_kinship = mean_kinship,
        indexes_ind = indexes_ind,
        tol = tol,
        want_inbr = FALSE,
        want_b = TRUE
    )
    Z1 <- obj$Z
    b <- obj$b
    
    # need null model too
    # just remove first column of these two
    Y0 <- Y1[ , -1, drop = FALSE ]
    Z0 <- Z1[ , -1, drop = FALSE ]
    
    # calculate other intermediate parts
    # all have dimensions n x k
    H1 <- Z1 %*% solve( crossprod( Y1, Z1 ) )
    H0 <- Z0 %*% solve( crossprod( Y0, Z0 ) )

    if ( V == 2 ) {
        HZ1 <- tcrossprod( H1, Z1 )
        HZ0 <- tcrossprod( H0, Z0 )
        # O( n^2*k )
    }
    
    ##############################
    ### COEFFICIENT ESTIMATION ###
    ##############################

    # initialize output vectors
    # Do before get_mem_lim_m so free memory is accounted for properly
    beta <- vector('numeric', m_loci)
    f_stat <- vector('numeric', m_loci)
    df <- vector('numeric', m_loci)
    
    # this overcounts since there are logical branches, not all overlap, but meh seriously
    # as usual, this should be conservative
    # recall that ints count as 0.5, doubles as 1
    #
    # vec_m, int
    # # indexes_loci_chunk
    # vec_m, double
    # # drop( Xi %*% proj )
    # # res1, res0
    #
    # vec_n # int
    # # n_ind_no_NA
    #
    # mat_m_n int
    # # Xi
    # # M # unnamed
    # # ( Xi == 1 )
    # mat_m_n double
    # # Res1, Res0
    # add one more double copy of Xi, this happens in matrix operations, are only temporary unnamed matrices
    
    # estimating total memory usage in bytes
    data <- popkin:::solve_m_mem_lim(
                         n = n_ind,
                         m = m_loci,
                         mat_m_n = 5,
                         vec_m = 3.5,
                         vec_n = 0.5,
                         mem = mem_lim,
                         mem_factor = mem_factor
                     )
    m_chunk <- data$m_chunk
    # cap value to a nice performing value (very good speed, minimal memory)
    if ( m_chunk > m_chunk_max )
        m_chunk <- m_chunk_max

    # navigate chunks
    i_chunk <- 1 # start of first chunk (needed for matrix inputs only; as opposed to function inputs)
    while (TRUE) { # start an infinite loop, break inside as needed
        # this means all SNPs have been covered!
        if (i_chunk > m_loci)
            break
        
        # range of SNPs to extract in this chunk
        indexes_loci_chunk <- i_chunk : min(i_chunk + m_chunk - 1, m_loci)

        if ( !is.null( indexes_ind ) ) {
            # individuals get filtered here (indexes_ind; required when there's missingness in trait)
            if (loci_on_cols) {
                Xi <- t( X[ indexes_ind, indexes_loci_chunk, drop = FALSE ] ) # transpose for our usual setup
            } else {
                Xi <- X[ indexes_loci_chunk, indexes_ind, drop = FALSE ]
            }
        } else {
            # keep all individuals (is this faster in that case?)
            if (loci_on_cols) {
                Xi <- t( X[ , indexes_loci_chunk, drop = FALSE ] ) # transpose for our usual setup
            } else {
                Xi <- X[ indexes_loci_chunk, , drop = FALSE ]
            }
        }
        
        # to have good averages, we need the number of non-NA individuals per row
        n_ind_no_NA <- rowSums( !is.na(Xi) )
        # now we can turn all NAs to zeroes (as ints, lower mem)
        Xi[ is.na(Xi) ] <- 0L
        
        # the coefficients are simply the genotypes projected!
        # these are matrices though, a vector for every locus (entry for every covariate)
        # adjust for the NAs? (not sure if this is reasonable or not yet)
        # store trait coefficients
        # projection for trait coefficient H1[,1] only
        beta[ indexes_loci_chunk ] <- drop( Xi %*% H1[ , 1 ] ) * n_ind / n_ind_no_NA

        if ( V == 0 ) {
            # rest are for getting residuals
            # NOTE: here missing values are just not part of sums, so setting them to zero is perfectly fine
            # alt model first
            R <- Xi - tcrossprod( Xi %*% H1, Y1 )
            # O( m*n*k )
            # then sums of residuals weighted by inverse kinship
            KiR <- conj_grad_scan(
                X = X,
                Y = t(R),
                b = b,
                indexes_ind = indexes_ind,
                tol = tol,
                want_inbr = FALSE
            )
            # O( n^2.5*m + n^2*m )
            # TODO: come up with a smart way to do this without explicit transpose
            ssr1 <- rowSums( t( KiR ) * R )

            # repeat for null model now
            R <- Xi - tcrossprod( Xi %*% H0, Y0 )
            KiR <- conj_grad_scan(
                X = X,
                Y = t(R),
                b = b,
                indexes_ind = indexes_ind,
                tol = tol,
                want_inbr = FALSE
            )
            # O( n^2.5*m + n^2*m )
            # TODO: come up with a smart way to do this without explicit transpose
            ssr0 <- rowSums( t( KiR ) * R )
        } else if ( V == 1 ) {
            ssr1 <- rowSums( ( Xi %*% H1 ) * ( Xi %*% Z1 ) )
            ssr0 <- rowSums( ( Xi %*% H0 ) * ( Xi %*% Z0 ) )
            # O( m*n*k + m*k^2 )
        } else if ( V == 2 ) {
            ssr1 <- rowSums( ( Xi %*% HZ1 ) * Xi )
            ssr0 <- rowSums( ( Xi %*% HZ0 ) * Xi )
            # O( m*n^2 )
        }

        if ( V == 1 || V == 2 ) {
            # then sums of residuals weighted by inverse kinship
            KiXi <- conj_grad_scan(
                X = X,
                Y = t(Xi),
                b = b,
                indexes_ind = indexes_ind,
                tol = tol,
                want_inbr = FALSE
            )
            # O( n^2.5*m + n^2*m )
            ssrx <- rowSums( t( KiXi ) * Xi )
        }
        
        # calculate F statistic now that all the parts are in place!
        # here NAs are accounted for in formula (to be normalized in the end!
        if ( V == 0 ) {
            f_stat[ indexes_loci_chunk ] <- (ssr0 - ssr1) / ssr1
        } else if ( V == 1 || V == 2 ) {
            f_stat[ indexes_loci_chunk ] <- (ssr1 - ssr0) / (ssrx - ssr1)
        }
        df[ indexes_loci_chunk ] <- n_ind_no_NA - ncol( Y1 )
        
        # update starting point for next chunk! (overshoots at the end, that's ok)
        i_chunk <- i_chunk + m_chunk
    }
    # all big-Os simplified assuming small k (detailed is unsimplified)
    #
    # V=0: O( m*n^2.5 )
    # detailed: O( m*n*k + n^2.5*m + n^2*m )
    # does appear worse by trading some n's by m's
    #
    # V=1: O( m*n^2.5 )
    # detailed: O( m*n*k + m*k^2 + n^2.5*m + n^2*m )
    # practically the same as V=0, if not a tad worse :(
    #
    # V=2: O( m*n^2.5 )
    # detailed: O( n^2*k + m*n^2 + n^2.5*m )
    # better than V=1 trading one m by n
    
    ################
    ### P-VALUES ###
    ################

    # normalize stats
    f_stat <- f_stat * df
    
    # Get p-values for F test!!! (should be exact)
    pval <- stats::pf( f_stat, 1, df, lower.tail = FALSE )
    
    # done, return quantities of interest (nice table!)
    return(
        tibble::tibble(
            pval = pval,
            beta = beta,
            f_stat = f_stat,
            df = df
        )
    )
}
